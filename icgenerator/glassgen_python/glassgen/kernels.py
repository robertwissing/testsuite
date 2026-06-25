"""SPH kernel in the gasoline convention.

q = r/h with support q in [0, 2] (gasoline BALL2 = (2h)^2, ar2 = q^2).

The kernel is the **Wendland C2** (3D), used because it has no pairing
instability at high neighbour counts: a particle whose ball the density-based-h
`nsmoothmin` floor expands toward high neighbour numbers must not clump, or the
ball path's rho->h->rho feedback runs away (see relax.py / ball.py). It needs
~64+ neighbours for good density accuracy (more than a cubic spline would).

The kernel SHAPE carries the full normalization so the density/force loops use
the gasoline fNorm = 1/(pi h^3). With support 2h the 3D Wendland C2 is

    W(r, 2h) = (21 / (16 pi h^3)) (1 - q/2)^4 (1 + 2q)

so w_wc2(ar2) = (pi h^3) * W = (21/16) (1 - q/2)^4 (1 + 2q), whose central value
w_wc2(0) = 21/16 = W0 is used by the explicit self contribution. dw_wc2 =
(1/q) dW_shape/dq = -(105/16) (1 - q/2)^3, so the gradient term is
dw_wc2(ar2) * fNorm1 * dx with fNorm1 = fNorm/h^2.
"""
import math

import numpy as np
from numba import njit

# central value of the kernel shape, w_wc2(0); used by the density/Q0 self terms.
W0 = 21.0 / 16.0

# Wendland C2 self-density bias correction Wzero10, ported verbatim from the
# piecewise-linear ladder in gasoline master.c (the Wzero10 branch, nSmooth
# 2..512): breakpoints nSmooth -> factor, linear between. The bare central value
# W0 = 21/16 over-counts a particle's own contribution, so a relaxed UNIFORM
# glass settles at rho slightly ABOVE the per-particle mean; scaling the self
# term by this factor removes the bias (uniform glass -> rho ~ 1.0). The C
# density self term is (21/16)*Wzero10 (smoothfcn.c KERNEL10).
_WZERO10_NS = np.array([2, 16, 24, 32, 48, 64, 96, 128, 192, 256, 512],
                       dtype=np.float64)
_WZERO10_F = np.array([0.50, 0.80507812, 0.90039062, 0.91347656, 0.93349609,
                       0.95693359, 0.97158203, 0.97675781, 0.98554688,
                       0.98847656, 0.99394531], dtype=np.float64)


def wzero_wc2(nsmooth):
    """Self-density bias correction factor for the Wendland C2 self term
    (gasoline smf->Wzero10). Piecewise-linear in nSmooth over the master.c
    table (0.95693359 at nSmooth=64); the corrected self central value is
    W0 * wzero_wc2(nSmooth), so the density self term becomes
    mass_i * W0 * factor instead of mass_i * W0. Endpoints held outside
    [2, 512] (C leaves nSmooth>512 uncorrected; nSmooth defaults to 64).
    Pure-Python (evaluated once per run, not in the njit loop)."""
    return float(np.interp(nsmooth, _WZERO10_NS, _WZERO10_F))


@njit(cache=True, inline='always')
def w_wc2(ar2):
    if ar2 >= 4.0:
        return 0.0
    q = math.sqrt(ar2)
    t = 1.0 - 0.5 * q            # (1 - q/2)
    t2 = t * t
    return (21.0 / 16.0) * t2 * t2 * (1.0 + 2.0 * q)


@njit(cache=True, inline='always')
def dw_wc2(ar2):
    # (1/q) dW_shape/dq ; finite at q->0 (-105/16), and the dx factor kills the
    # self term anyway, but keep the ar2<=0 guard against a divide-by-zero.
    if ar2 >= 4.0 or ar2 <= 0.0:
        return 0.0
    q = math.sqrt(ar2)
    t = 1.0 - 0.5 * q
    return -(105.0 / 16.0) * t * t * t
