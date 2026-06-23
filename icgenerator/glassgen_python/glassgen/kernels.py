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

from numba import njit

# central value of the kernel shape, w_wc2(0); used by the density/Q0 self terms.
W0 = 21.0 / 16.0


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
