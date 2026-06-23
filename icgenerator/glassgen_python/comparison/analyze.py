"""Apples-to-apples glass-quality comparison of a gasoline.gdicgen (C) glass
against a glassgen (Python) glass.

For each glass file we recompute, with the *same* glassgen pipeline
(kNN h + M4 density + ICGEN diagnostics), the SPH density rho, the
partition-of-unity Q1 (-> 1 for a good glass), and the zero-order force
error |E0|. Comparing these on both outputs removes any bias from each
code's internally reported density.

Glass-quality metrics reported / plotted:
  * density error  |1 - rho/rho0|  vs the target rho0 (uniform mean, or a
    density table / callable for structured cases)
  * partition of unity Q1
  * nearest-neighbour distance distribution (glass regularity)

Usage:
  python analyze.py --label-a C --file-a cs_C.02710 \
                    --label-b Python --file-b cs_py_IC \
                    --param pre.param -s 64 --tag currentsheet \
                    [--table densitytable_xdr --direction r_sph]
"""
import argparse
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_PKG = os.path.join(_HERE, '..')
_SETUP = os.path.join(_HERE, '..', '..', '..', 'setupfiles')
sys.path.insert(0, _PKG)
sys.path.insert(0, _SETUP)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import readtipsy as tip
from glassgen.neighbors import NeighborList
from glassgen.sph_core import density_loop, get_force_loop
from glassgen.density import TableDensity, CallableDensity


def read_box(path):
    period = {'dPeriod': 1.0, 'dxPeriod': None, 'dyPeriod': None,
              'dzPeriod': None}
    with open(path) as f:
        for line in f:
            line = line.split('#')[0]
            if '=' not in line:
                continue
            k, v = (s.strip() for s in line.split('=', 1))
            if k in period:
                period[k] = float(v)
    d = period['dPeriod']
    return np.array([period['dxPeriod'] or d, period['dyPeriod'] or d,
                     period['dzPeriod'] or d])


def measure(path, box, nsmooth, density):
    """Recompute rho, h, Q1, |E0| and nearest-neighbour distance for a glass."""
    gas, _, _, _, _ = tip.readtipsy(path)
    pos = gas[:, 1:4].astype(np.float64)
    mass = gas[:, 0].astype(np.float64)
    box = np.asarray(box, float)
    pos -= box * np.rint(pos / box)            # wrap, box-centered on 0
    n = len(pos)

    nb = NeighborList(box, nsmooth, margin=8, build_transpose=True)
    idx, h, rev_indptr, rev_indices = nb.update(pos)
    rho = np.empty(n)
    density_loop(pos, mass, h, idx, box, True, rho)

    icvel = np.empty((n, 3))
    q1 = np.empty(n)
    e0 = np.empty((n, 3))
    pP = np.ones(n)
    fl = get_force_loop('gdforce', one_sided=False)
    fl(pos, mass, h, rho, pP, idx, rev_indptr, rev_indices, box, True,
       icvel, q1, e0, True)
    e0mag = np.sqrt((e0 ** 2).sum(1))

    # nearest-neighbour distance: idx[:,0] is self, idx[:,1] first neighbour
    j = idx[:, 1]
    d = pos - pos[j]
    d -= box * np.rint(d / box)
    nnd = np.sqrt((d ** 2).sum(1))

    rho0 = density.rho0(pos)
    return dict(pos=pos, mass=mass, h=h, rho=rho, rho0=rho0, q1=q1,
                e0=e0mag, nnd=nnd, n=n)


def summarize(label, m):
    err = np.abs(1.0 - m['rho'] / m['rho0'])
    nnd = m['nnd'] / np.median(m['nnd'])     # normalize to median spacing
    return {
        'label': label, 'n': m['n'],
        'rho_cov': m['rho'].std() / m['rho'].mean(),
        'err_mean': err.mean(), 'err_med': np.median(err),
        'err_p95': np.percentile(err, 95), 'err_max': err.max(),
        'q1_mean': m['q1'].mean(), 'q1_std': m['q1'].std(),
        'e0_mean': m['e0'].mean(), 'e0_p95': np.percentile(m['e0'], 95),
        'nnd_cov': nnd.std(),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--file-a', required=True)
    ap.add_argument('--file-b', required=True)
    ap.add_argument('--label-a', default='A')
    ap.add_argument('--label-b', default='B')
    ap.add_argument('--param', required=True)
    ap.add_argument('-s', dest='nsmooth', type=int, default=64)
    ap.add_argument('--table', default=None)
    ap.add_argument('--direction', default='r_sph')
    ap.add_argument('--tag', default='cmp')
    args = ap.parse_args()

    box = read_box(args.param)
    print(f'box = {box.tolist()}', flush=True)

    if args.table:
        density = TableDensity.from_xdr(args.table, mapping=args.direction)
    else:
        # uniform target = total mass / volume (read from file A)
        gas, _, _, _, _ = tip.readtipsy(args.file_a)
        rho_mean = gas[:, 0].astype(np.float64).sum() / np.prod(box)
        density = CallableDensity(lambda p, r=rho_mean: np.full(len(p), r))
        print(f'uniform target rho0 = {rho_mean:.6g}', flush=True)

    ma = measure(args.file_a, box, args.nsmooth, density)
    mb = measure(args.file_b, box, args.nsmooth, density)
    sa = summarize(args.label_a, ma)
    sb = summarize(args.label_b, mb)

    cols = ['n', 'rho_cov', 'err_mean', 'err_med', 'err_p95', 'err_max',
            'q1_mean', 'q1_std', 'e0_mean', 'e0_p95', 'nnd_cov']
    w = max(len(c) for c in cols)
    print(f'\n=== {args.tag}: glass-quality comparison '
          f'(recomputed with glassgen, nsmooth={args.nsmooth}) ===')
    print(f'{"metric":<{w}} {sa["label"]:>14} {sb["label"]:>14}')
    for c in cols:
        fa = f'{sa[c]:.5g}' if c != 'n' else str(sa[c])
        fb = f'{sb[c]:.5g}' if c != 'n' else str(sb[c])
        print(f'{c:<{w}} {fa:>14} {fb:>14}')

    # figures
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    erra = np.abs(1 - ma['rho'] / ma['rho0'])
    errb = np.abs(1 - mb['rho'] / mb['rho0'])
    bins = np.linspace(0, np.percentile(np.r_[erra, errb], 99.5), 60)
    ax[0].hist(erra, bins, histtype='step', density=True, label=args.label_a)
    ax[0].hist(errb, bins, histtype='step', density=True, label=args.label_b)
    ax[0].set_xlabel('|1 - rho/rho0|'); ax[0].set_ylabel('pdf')
    ax[0].set_title('density error'); ax[0].legend()

    qb = np.linspace(min(ma['q1'].min(), mb['q1'].min()),
                     max(ma['q1'].max(), mb['q1'].max()), 60)
    ax[1].hist(ma['q1'], qb, histtype='step', density=True, label=args.label_a)
    ax[1].hist(mb['q1'], qb, histtype='step', density=True, label=args.label_b)
    ax[1].axvline(1.0, color='k', lw=0.6, ls=':')
    ax[1].set_xlabel('Q1 (partition of unity)'); ax[1].set_title('Q1')
    ax[1].legend()

    na = ma['nnd'] / np.median(ma['nnd'])
    nbn = mb['nnd'] / np.median(mb['nnd'])
    nb_ = np.linspace(0, 2.0, 60)
    ax[2].hist(na, nb_, histtype='step', density=True, label=args.label_a)
    ax[2].hist(nbn, nb_, histtype='step', density=True, label=args.label_b)
    ax[2].set_xlabel('nearest-neighbour dist / median')
    ax[2].set_title('glass regularity'); ax[2].legend()

    fig.suptitle(f'{args.tag}: {args.label_a} vs {args.label_b} glass')
    fig.tight_layout()
    out = f'{args.tag}_quality.png'
    fig.savefig(out, dpi=110)
    print(f'\nwrote {out}')


if __name__ == '__main__':
    main()
