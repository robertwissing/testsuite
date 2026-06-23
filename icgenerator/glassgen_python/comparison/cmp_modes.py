"""Tabulate glass quality (recomputed with glassgen kNN-h) across the C glass
and the three Python h-modes. Usage: python cmp_modes.py <param> <tag>
[--table T --direction D] file1=label1 file2=label2 ..."""
import sys
import numpy as np
from analyze import read_box, measure, summarize
from glassgen.density import TableDensity, CallableDensity

param = sys.argv[1]
tag = sys.argv[2]
rest = sys.argv[3:]
table = direction = None
pairs = []
i = 0
while i < len(rest):
    if rest[i] == '--table':
        table = rest[i + 1]; i += 2
    elif rest[i] == '--direction':
        direction = rest[i + 1]; i += 2
    else:
        f, lab = rest[i].split('=')
        pairs.append((f, lab)); i += 1

box = read_box(param)
if table:
    density = TableDensity.from_xdr(table, mapping=direction or 'r_sph')
else:
    import readtipsy as tip
    gas, *_ = tip.readtipsy(pairs[0][0])
    rho_mean = gas[:, 0].astype(float).sum() / np.prod(box)
    density = CallableDensity(lambda p, r=rho_mean: np.full(len(p), r))

rows = [summarize(lab, measure(f, box, 64, density)) for f, lab in pairs]
cols = ['rho_cov', 'err_mean', 'err_med', 'err_p95', 'err_max',
        'q1_mean', 'q1_std', 'e0_mean', 'nnd_cov']
w = 12
print(f'\n=== {tag}: glass quality (recomputed kNN-h, nsmooth=64) ===')
hdr = 'metric'.ljust(10) + ''.join(r['label'].rjust(w) for r in rows)
print(hdr)
for c in cols:
    print(c.ljust(10) + ''.join(f'{r[c]:.5g}'.rjust(w) for r in rows))
