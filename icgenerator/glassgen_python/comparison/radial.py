"""Radial density profile overlay for the structured (mhdcollapse) case:
recompute SPH rho(r) with glassgen for the C and Python glasses and plot
both against the target tanh profile."""
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, '..'))
sys.path.insert(0, os.path.join(_HERE, '..', '..', '..', 'setupfiles'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from analyze import read_box, measure
from glassgen.density import TableDensity

box = read_box('mhd.param')
density = TableDensity.from_xdr('densitytable_xdr', mapping='r_sph')

inner, outer, R, Rs = 70.73553026306462, 0.19648758406406838, 0.015, 0.0015

def radial_binned(path):
    m = measure(path, box, 64, density)
    r = np.sqrt((m['pos'] ** 2).sum(1))
    edges = np.linspace(0, 0.06, 80)
    idx = np.digitize(r, edges)
    rho_m, rho_lo, rho_hi, rc = [], [], [], []
    for b in range(1, len(edges)):
        sel = idx == b
        if sel.sum() < 5:
            continue
        rr = m['rho'][sel]
        rho_m.append(np.median(rr))
        rho_lo.append(np.percentile(rr, 16))
        rho_hi.append(np.percentile(rr, 84))
        rc.append(0.5 * (edges[b - 1] + edges[b]))
    return np.array(rc), np.array(rho_m), np.array(rho_lo), np.array(rho_hi)

rgrid = np.linspace(0, 0.06, 500)
target = outer + (inner - outer) / (1.0 + np.exp(-(R - rgrid) / Rs))

fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
for a in ax:
    a.plot(rgrid, target, 'k-', lw=1.3, label='target')
for path, lab, col in [('mhd_C.03692', 'C-gdicgen', 'C0'),
                       ('mhd_py_IC', 'Python-glassgen', 'C1')]:
    rc, rm, rlo, rhi = radial_binned(path)
    for a in ax:
        a.plot(rc, rm, '-', color=col, lw=1.1, label=lab)
        a.fill_between(rc, rlo, rhi, color=col, alpha=0.18)
ax[0].axvline(R, color='gray', ls=':', lw=0.7)
ax[0].set_xlabel('r'); ax[0].set_ylabel('rho (SPH, glassgen)')
ax[0].set_title('radial density profile'); ax[0].legend(fontsize=8)
ax[1].set_yscale('log'); ax[1].axvline(R, color='gray', ls=':', lw=0.7)
ax[1].set_xlabel('r'); ax[1].set_title('radial profile (log rho)')
ax[1].legend(fontsize=8)
fig.suptitle('mhdcollapse50: SPH rho(r) vs target (band = 16-84%)')
fig.tight_layout()
fig.savefig('mhdcollapse50_radial.png', dpi=110)
print('wrote mhdcollapse50_radial.png')
