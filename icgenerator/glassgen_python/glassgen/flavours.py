"""SPH flavour registry: pseudo-pressure prefactor pairs per flavour.

A flavour is an integer id dispatched inside ig_pair, which the force loop
calls for the (ig_i, ig_j) prefactors of the pair term, mirroring
ICGenerator (smoothfcn.c): GDFORCE uses pP/(rho_i rho_j), standard uses
pP/rho^2. Integer dispatch (rather than inlining per-flavour closures via a
factory) keeps the force loop a top-level @njit(cache=True) function, so it
compiles once ever instead of on every process start.

Future ISPH flavours additionally replace the separation vector dx by the
matrix-corrected -(C @ dx) and weight pairs by the kernel value instead of
its derivative; the force loop in sph_core reserves the cmat slot for that.
"""
from numba import njit

FLAVOURS = {
    'gdforce': 0,
    'rho2': 1,
}


@njit(cache=True, inline='always')
def ig_pair(iflav, pP_i, pP_j, rho_i, rho_j):
    if iflav == 0:  # gdforce
        ig = 1.0 / (rho_i * rho_j)
        return pP_i * ig, pP_j * ig
    # rho2
    return pP_i / (rho_i * rho_i), pP_j / (rho_j * rho_j)
