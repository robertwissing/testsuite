"""SPH smoothing kernels — the single source of truth for `sph_interp`.

Both interpolation methods (regular SPH deposit and the Petkova exact integral)
use the SAME kernel so results are directly comparable. The default is the
M4 cubic spline with compact support radius 2h (the convention used throughout
this code: `q = r/h`, support `q < 2`).

Normalization (3D): W(r,h) = cnorm3D(h) * w(q),  cnorm3D(h) = 1 / (pi h^3).
"""

import numpy as np
from numba import njit

# Support radius of the cubic spline in units of h, and its square.
RADKERNEL = 2.0
RADKERNEL2 = RADKERNEL * RADKERNEL


@njit(cache=True)
def wkernel(q2):
    """Dimensionless M4 cubic spline w(q) as a function of q2 = (r/h)^2.

    Multiply by cnorm3D(h) = 1/(pi h^3) for the normalized 3D kernel W(r,h).
    Zero for q >= 2.
    """
    q = np.sqrt(q2)
    if q < 1.0:
        return 1.0 - 1.5 * q2 + 0.75 * q * q2
    elif q < 2.0:
        t = 2.0 - q
        return 0.25 * t * t * t
    else:
        return 0.0


@njit(cache=True)
def cnorm3D(h):
    """3D normalization of the cubic spline: 1 / (pi h^3)."""
    return 1.0 / (np.pi * h * h * h)


# --------------------------------------------------------------------------- #
#  Column-integrated kernel (for 2D projections)
# --------------------------------------------------------------------------- #
# The line-of-sight integral of the 3D kernel:
#     Y(R,h) = int_{-inf}^{inf} W(sqrt(R^2+z^2), h) dz = (1/(pi h^2)) * Sigma(R/h)
# where Sigma(s) = int_{-inf}^{inf} w(sqrt(s^2+t^2)) dt is the dimensionless column
# profile of w(q), nonzero for s < 2. Tabulated once (SPLASH-style) and linearly
# interpolated. Normalized so that 2*int_0^2 Sigma(s) s ds = 1 (=> a column deposit
# of mass conserves the total mass as a surface density).

COL_NS = 1024
COL_SMAX = 2.0


def _w_cubic_vec(q):
    q = np.asarray(q, dtype=np.float64)
    out = np.zeros_like(q)
    m1 = q < 1.0
    m2 = (q >= 1.0) & (q < 2.0)
    out[m1] = 1.0 - 1.5 * q[m1] ** 2 + 0.75 * q[m1] ** 3
    out[m2] = 0.25 * (2.0 - q[m2]) ** 3
    return out


def _build_col_table(ns=COL_NS, nt=4096):
    s = np.linspace(0.0, COL_SMAX, ns)
    vals = np.empty(ns)
    for i in range(ns):
        si = s[i]
        tmax = np.sqrt(max(4.0 - si * si, 0.0))
        if tmax <= 0.0:
            vals[i] = 0.0
            continue
        t = np.linspace(0.0, tmax, nt)
        vals[i] = 2.0 * np.trapz(_w_cubic_vec(np.sqrt(si * si + t * t)), t)
    return vals


COL_TABLE = _build_col_table()


@njit(cache=True)
def col_profile(s, table):
    """Dimensionless column profile Sigma(s) via linear interp of `table`."""
    if s >= COL_SMAX or s < 0.0:
        return 0.0
    ns = table.shape[0]
    ds = COL_SMAX / (ns - 1)
    x = s / ds
    i = int(x)
    if i >= ns - 1:
        return table[ns - 1]
    f = x - i
    return table[i] * (1.0 - f) + table[i + 1] * f


# Cumulative column profile P(t) = int_0^t Sigma(s) s ds, on the same s-grid.
# Used by the EXACT (Petkova) 2D column integral: the mass of a particle's column
# enclosed within radius r is M(r) = (1/pi) P(r/h), and P(2) = 0.5 (so the full
# plane integrates to 1). Saturates at 0.5 for t >= 2.
def _build_pcum_table():
    s = np.linspace(0.0, COL_SMAX, COL_NS)
    integ = COL_TABLE * s
    p = np.zeros(COL_NS)
    p[1:] = np.cumsum(0.5 * (integ[1:] + integ[:-1]) * np.diff(s))
    return p


PCUM_TABLE = _build_pcum_table()


@njit(cache=True)
def pcum(t, table):
    """Cumulative column profile P(t) = int_0^t Sigma(s) s ds (saturates at 0.5)."""
    if t <= 0.0:
        return 0.0
    if t >= COL_SMAX:
        return table[table.shape[0] - 1]
    ns = table.shape[0]
    ds = COL_SMAX / (ns - 1)
    x = t / ds
    i = int(x)
    if i >= ns - 1:
        return table[ns - 1]
    f = x - i
    return table[i] * (1.0 - f) + table[i + 1] * f
