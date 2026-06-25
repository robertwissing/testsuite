#!/usr/bin/env python
"""Spatial diagnostics for a relaxed glass: SPH rho recomputed with glassgen,
then (1) a density profile along the case's density direction vs the target
and (2) a thin-slab slice of the particle distribution coloured by log rho.

Works on any tipsy glass + its .param (target/box/direction taken from the
param's dICdens* fields, like the C generator). Optionally overlays a second
glass for side-by-side comparison (e.g. C-gdicgen vs Python).

    python spatial.py mhd_py_IC --param mhd.param --tag mhd_py
    python spatial.py mhd_C.03692 --param mhd.param \
        --compare mhd_py_IC --labels C,Python --tag mhd_cmp
"""
import argparse
import os
import sys

import numpy as np

_PKG = os.path.dirname(os.path.abspath(__file__))
_TESTSUITE = os.path.dirname(os.path.dirname(_PKG))
sys.path.insert(0, _PKG)
sys.path.insert(0, os.path.join(_TESTSUITE, 'setupfiles'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import readtipsy as tip
from glassgen.kernels import W0, wzero_wc2
from glassgen.neighbors import NeighborList
from glassgen.sph_core import density_loop
from benchmark import read_param, box_from_param, target_from_param

_DIRNAME = {1: 'x', 2: 'y', 3: 'z', 4: 'r_sph', 5: 'r_cyl'}


def mapped_coord(pos, d):
    if d < 4:
        return np.abs(pos[:, d - 1]), _DIRNAME[d]
    if d == 4:
        return np.sqrt((pos ** 2).sum(1)), 'r_sph'
    if d == 5:
        return np.sqrt(pos[:, 0] ** 2 + pos[:, 1] ** 2), 'r_cyl'
    return pos[:, d - 6], 'xyz'[d - 6]


def load(path, box, nsmooth):
    gas, _, _, _, _ = tip.readtipsy(path)
    pos = gas[:, 1:4].astype(np.float64)
    mass = gas[:, 0].astype(np.float64)
    pos -= box * np.rint(pos / box)
    nb = NeighborList(box, nsmooth, margin=8, build_transpose=False)
    idx, h, _, _ = nb.update(pos)
    rho = np.empty(len(pos))
    density_loop(pos, mass, h, idx, box, True, rho, W0 * wzero_wc2(nsmooth))
    return pos, rho


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('glass')
    ap.add_argument('--param', required=True)
    ap.add_argument('--compare', default=None)
    ap.add_argument('--labels', default='A,B')
    ap.add_argument('--tag', default='spatial')
    ap.add_argument('-s', '--nsmooth', type=int, default=64)
    ap.add_argument('--slice-axis', default='z', choices='xyz')
    ap.add_argument('--slice-frac', type=float, default=0.12,
                    help='half-thickness of the slab as a fraction of the box')
    args = ap.parse_args()

    param = read_param(args.param)
    box = box_from_param(param)
    target = target_from_param(param)
    d = int(float(param.get('dICdensdir', 1)))
    labels = args.labels.split(',')

    sets = [(args.glass, labels[0])]
    if args.compare:
        sets.append((args.compare, labels[1]))
    data = [(lab,) + load(p, box, args.nsmooth) for p, lab in sets]

    # ---- profile along the density direction + L1 error vs analytic ----
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(7, 6.4), sharex=True,
                                  gridspec_kw=dict(height_ratios=[2, 1]))
    cmin = min(mapped_coord(pos, d)[0].min() for _, pos, _ in data)
    cmax = max(np.percentile(mapped_coord(pos, d)[0], 99.5)
               for _, pos, _ in data)
    edges = np.linspace(cmin, cmax, 70)
    cc = 0.5 * (edges[:-1] + edges[1:])
    print(f'\n=== {args.tag}: L1 density error vs analytic profile ===')
    for (lab, pos, rho), col in zip(data, ('C0', 'C1')):
        c, cname = mapped_coord(pos, d)
        rho0 = target(pos)
        ib = np.digitize(c, edges) - 1
        med = np.full(len(cc), np.nan)
        lo = np.full(len(cc), np.nan)
        hi = np.full(len(cc), np.nan)
        t_b = np.full(len(cc), np.nan)
        for b in range(len(cc)):
            m = ib == b
            if m.any():
                med[b] = np.median(rho[m])
                lo[b] = np.percentile(rho[m], 16)
                hi[b] = np.percentile(rho[m], 84)
                t_b[b] = np.mean(rho0[m])
        ax.plot(cc, med, '-', color=col, lw=1.2, label=f'{lab} (SPH)')
        ax.fill_between(cc, lo, hi, color=col, alpha=0.2)
        # residual of the binned profile vs the analytic target
        resid = med / t_b - 1.0
        axr.plot(cc, resid, '-', color=col, lw=1.1, label=lab)
        # L1 errors: per-particle relative, and binned-profile (abs + rel)
        ok = np.isfinite(med)
        l1_particle = np.mean(np.abs(rho / rho0 - 1.0))
        l1_prof_rel = np.mean(np.abs(resid[ok]))
        l1_prof_abs = np.mean(np.abs(med[ok] - t_b[ok]))
        print(f'  {lab:<16} L1(per-particle, rel)={l1_particle:.4e}  '
              f'L1(profile, rel)={l1_prof_rel:.4e}  '
              f'L1(profile, abs)={l1_prof_abs:.4e}')
    pos0 = data[0][1]
    c0, cname = mapped_coord(pos0, d)
    o = np.argsort(c0)
    ax.plot(c0[o], target(pos0)[o], 'k-', lw=1.3, label='target')
    ax.set_ylabel('rho (SPH, glassgen)')
    if target(pos0).max() / max(target(pos0).min(), 1e-30) > 20:
        ax.set_yscale('log')
    ax.set_title(f'{args.tag}: density profile along {cname}')
    ax.legend(fontsize=8)
    axr.axhline(0.0, color='k', lw=0.6, ls=':')
    axr.set_xlabel(cname)
    axr.set_ylabel('rho_SPH/rho0 - 1')
    axr.set_ylim(-0.3, 0.3)
    axr.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f'{args.tag}_profile.png', dpi=120)
    print(f'wrote {args.tag}_profile.png')

    # ---- particle-distribution slice ----
    ai = 'xyz'.index(args.slice_axis)
    ax2i = [j for j in range(3) if j != ai]
    half = args.slice_frac * 0.5 * box[ai]
    n = len(data)
    fig, axes = plt.subplots(1, n, figsize=(6.2 * n, 6), squeeze=False)
    vmin = min(np.log10(rho).min() for _, _, rho in data)
    vmax = max(np.percentile(np.log10(rho), 99.5) for _, _, rho in data)
    for k, (lab, pos, rho) in enumerate(data):
        sel = np.abs(pos[:, ai]) < half
        a = axes[0, k]
        sc = a.scatter(pos[sel, ax2i[0]], pos[sel, ax2i[1]],
                       c=np.log10(rho[sel]), s=3, cmap='magma',
                       vmin=vmin, vmax=vmax, linewidths=0)
        a.set_aspect('equal')
        a.set_xlabel('xyz'[ax2i[0]]); a.set_ylabel('xyz'[ax2i[1]])
        a.set_title(f'{lab}: |{args.slice_axis}|<{half:.3g} '
                    f'({sel.sum()} part.)')
        fig.colorbar(sc, ax=a, label='log10 rho', shrink=0.8)
    fig.suptitle(f'{args.tag}: particle slice coloured by SPH log density')
    fig.tight_layout()
    fig.savefig(f'{args.tag}_slice.png', dpi=120)
    print(f'wrote {args.tag}_slice.png')


if __name__ == '__main__':
    main()
