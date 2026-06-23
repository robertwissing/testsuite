"""Standalone split / merge CLI for tipsy snapshots (all aux fields).

Reads a tipsy gas snapshot plus its aux sibling files (BField, BFieldx/y/z,
iord, ...), splits or merges every particle field, and writes a new snapshot
with the changed particle count. Unlike the relax round-trip (cli.py), the
particle count CHANGES here, so aux files are re-emitted with the new count
header rather than byte-copied.

Examples:
    python -m glassgen.resample_cli -I pre.00000 -o pre_x8 --op split \
        --factor 8 --param pre.param
    python -m glassgen.resample_cli -I fine.00000 -o coarse --op merge \
        --factor 8 --param fine.param
"""
import argparse
import glob
import os
import sys

import numpy as np

_SETUPFILES = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           '..', '..', '..', 'setupfiles')
sys.path.insert(0, _SETUPFILES)

from . import resample
from .cli import read_box_from_param

# tipsy gas columns -> resample aux keys (mass=0, pos=1:4 handled explicitly)
_GAS_SCALAR = {'rho': 7, 'u': 8, 'h': 9, 'metals': 10, 'phi': 11}
# fields halved on split / summed on merge (extensive): u + every B-field aux
_GAS_HALVE = {'u'}


def _read_aux(infile):
    """Read all aux sibling files of `infile`. Returns a list of dicts with
    keys: label, key (aux-dict key), data (N,) or (N,3), kind in
    {'scalar','vector','int'}."""
    import readtipsy as tip
    out = []
    for path in sorted(glob.glob(infile + '.*')):
        label = path[len(infile) + 1:]
        if any(b in label for b in tip.banauxfiles):
            continue
        if label.endswith('.param') or label == 'param':
            continue
        if label in tip.vecauxfiles:
            data, _ = tip.readtipsyVecaux(infile, label)
            data = np.asarray(data, dtype=np.float64).reshape(-1, 3)
            kind = 'vector'
        else:
            data, _ = tip.readtipsyaux(infile, label)
            data = np.asarray(data)
            kind = 'int' if label in tip.intauxfiles else 'scalar'
        out.append(dict(label=label, key='aux_' + label, data=data,
                        kind=kind))
    return out


def _write_aux(meta, outfile):
    import readtipsy as tip
    for m in meta:
        data = m['data']
        if m['kind'] == 'vector':
            tip.writetipsyVecaux(np.asarray(data, dtype=np.float64), m['label'],
                                 outfile)
        else:
            tip.writetipsyaux(np.asarray(data), m['label'], outfile)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description='split/merge a tipsy snapshot (gas + all aux fields)')
    ap.add_argument('-I', dest='infile', required=True, help='input tipsy')
    ap.add_argument('-o', dest='out', required=True,
                    help='output prefix (writes ${out}.00000 + aux)')
    ap.add_argument('--op', choices=('split', 'merge'), default='split')
    ap.add_argument('--factor', type=int, default=8,
                    help='multiply (split) or divide (merge) the count by this '
                         'power of 2')
    ap.add_argument('-s', dest='nsmooth', type=int, default=64)
    ap.add_argument('--dist', type=float, default=0.4,
                    help='split offset as a fraction of the neighbour distance')
    ap.add_argument('--margin', type=int, default=8)
    grp = ap.add_mutually_exclusive_group(required=True)
    grp.add_argument('--param', help='gasoline .param file to take the box from')
    grp.add_argument('--box', nargs=3, type=float, metavar=('LX', 'LY', 'LZ'))
    grp.add_argument('--open', action='store_true', help='no periodic box')
    args = ap.parse_args(argv)

    import readtipsy as tip

    gas, dark, star, header, time = tip.readtipsy(args.infile)
    if len(dark) or len(star):
        sys.exit('resample_cli supports gas-only snapshots (dark/star aux '
                 'splitting is not implemented)')
    n = len(gas)
    print(f'read {args.infile}: {n} gas particles', flush=True)

    if args.open:
        box = None
    elif args.box:
        box = np.array(args.box)
    else:
        box = read_box_from_param(args.param)

    pos = gas[:, 1:4].astype(np.float64)
    mass = gas[:, 0].astype(np.float64)
    aux = {'vel': gas[:, 4:7].astype(np.float64)}
    for key, col in _GAS_SCALAR.items():
        aux[key] = gas[:, col].astype(np.float64)

    auxmeta = _read_aux(args.infile)
    halve = set(_GAS_HALVE)
    for m in auxmeta:
        aux[m['key']] = m['data'].astype(np.float64)
        # B-field aux are extensive (flux split), matching the SPH split
        if m['label'].startswith('BField') or m['label'].startswith('CurlB'):
            halve.add(m['key'])

    if args.op == 'split':
        pos2, mass2, aux2 = resample.split(
            pos, mass, aux, box=box, factor=args.factor,
            nsmooth=args.nsmooth, margin=args.margin, dist=args.dist,
            halve_fields=halve)
    else:
        pos2, mass2, aux2 = resample.merge(
            pos, mass, aux, box=box, target_factor=args.factor,
            nsmooth=args.nsmooth, margin=args.margin, halve_fields=halve)
    n2 = len(pos2)
    print(f'{args.op}: {n} -> {n2} gas particles', flush=True)

    # rebuild the 12-column gas array
    g = np.zeros((n2, 12), dtype=gas.dtype)
    g[:, 0] = mass2
    g[:, 1:4] = pos2
    g[:, 4:7] = aux2['vel']
    for key, col in _GAS_SCALAR.items():
        g[:, col] = aux2[key]
    hdr = header.copy()
    hdr[0] = n2        # nbodies
    hdr[2] = n2        # ngas
    outsnap = f'{args.out}.00000'
    tip.writetipsy(g, dark, star, outsnap, hdr, time)

    # re-emit aux files with the new count
    for m in auxmeta:
        data = aux2[m['key']]
        if m['kind'] == 'int':
            # particle ids must stay unique: regenerate after a count change
            data = np.arange(n2, dtype=np.int32)
        m['data'] = data
    _write_aux(auxmeta, outsnap)
    print(f'wrote {outsnap}' + (f' + {len(auxmeta)} aux files'
                                if auxmeta else ''))


if __name__ == '__main__':
    main()
