"""Numba pair loops: SPH density, ICGenerator force, quality diagnostics.

Ports the ICGEN pair physics of icgenerator/gasoline:
- density: gather sum rho_i = 1/(pi h_i^3) sum_j m_j W(q_ij) including self
  (Density/DenDVDX, smoothfcn.c).
- force (ICGenerator, smoothfcn.c:5794): in C each mutual pair is visited
  twice by the symmetric resmooth (once with each owner's h) and the kernel
  carries fNorm1 = 0.5/(pi h^5), which yields the arithmetic-mean kernel
  symmetrization. Here that is reformulated gather-only: particle i
  accumulates its own-ball terms (kernel h_i, neighbors indices[i]) plus the
  h_j-side terms, so each prange thread writes only icvel[i] - no races.
- diagnostics (CHECKQUALITY): Q0_i = sum m_j/rho_j W_ij(h_i) (partition of
  unity, -> 1; the zeroth-order interpolation consistency - see
  diagnostics.py for the full Q0/Q1/E0/E1 hierarchy) and E0_i = zero-order
  force error (the force with pP == 1, scaled by rho_i h_i), accumulated in
  the same loop.

A SINGLE force loop (_force_loop) serves every backend: it walks CSR
neighbour lists (indptr/indices). The kNN backend passes indptr=arange*k,
indices=idx.ravel() (a view, identical iteration order); the ball backend
passes its native CSR. The h_j (transposed) side is taken from the precomputed
transpose graph (symmetric) or evaluated inline from i's own list (one-sided,
SPH-EXA style, no transpose) per the `one_sided` flag.
"""
import math

import numpy as np
from numba import njit, prange

from .flavours import FLAVOURS, ig_pair
from .kernels import w_wc2, dw_wc2, W0

M_1_PI = 1.0 / math.pi


@njit(cache=True, parallel=True)
def density_loop(pos, mass, h, idx, box, periodic, rho, w0self):
    # positions must be wrapped to [-L/2, L/2] when periodic: min-image is
    # then two predictable compares instead of a divide+rint per axis.
    # w0self is the (possibly bias-corrected) self central value W(0): pass W0
    # for the bare kernel, or W0*wzero_wc2(nsmooth) for the Wendland C2
    # self-density correction (uniform glass -> rho ~ 1.0).
    n, k = idx.shape
    bx, by, bz = box[0], box[1], box[2]
    hbx, hby, hbz = 0.5 * bx, 0.5 * by, 0.5 * bz
    for i in prange(n):
        hi = h[i]
        ih2 = 1.0 / (hi * hi)
        s = mass[i] * w0self  # self contribution, W(0)
        for jj in range(k):
            j = idx[i, jj]
            dx = pos[i, 0] - pos[j, 0]
            dy = pos[i, 1] - pos[j, 1]
            dz = pos[i, 2] - pos[j, 2]
            if periodic:
                if dx > hbx:
                    dx -= bx
                elif dx < -hbx:
                    dx += bx
                if dy > hby:
                    dy -= by
                elif dy < -hby:
                    dy += by
                if dz > hbz:
                    dz -= bz
                elif dz < -hbz:
                    dz += bz
            ar2 = (dx * dx + dy * dy + dz * dz) * ih2
            s += mass[j] * w_wc2(ar2)
        rho[i] = s * M_1_PI / (hi * hi * hi)


@njit(cache=True, parallel=True)
def density_loop_csr(pos, mass, h, indptr, indices, box, periodic, rho, w0self):
    """SPH density over a CSR (ball-gather) neighbour list: each particle i
    sums over indices[indptr[i]:indptr[i+1]] (all within ~2*h_i). w0self is the
    (possibly bias-corrected) self central value W(0) - see density_loop."""
    n = pos.shape[0]
    bx, by, bz = box[0], box[1], box[2]
    hbx, hby, hbz = 0.5 * bx, 0.5 * by, 0.5 * bz
    for i in prange(n):
        hi = h[i]
        ih2 = 1.0 / (hi * hi)
        s = mass[i] * w0self
        for e in range(indptr[i], indptr[i + 1]):
            j = indices[e]
            dx = pos[i, 0] - pos[j, 0]
            dy = pos[i, 1] - pos[j, 1]
            dz = pos[i, 2] - pos[j, 2]
            if periodic:
                if dx > hbx:
                    dx -= bx
                elif dx < -hbx:
                    dx += bx
                if dy > hby:
                    dy -= by
                elif dy < -hby:
                    dy += by
                if dz > hbz:
                    dz -= bz
                elif dz < -hbz:
                    dz += bz
            ar2 = (dx * dx + dy * dy + dz * dz) * ih2
            if ar2 < 4.0:
                s += mass[j] * w_wc2(ar2)
        rho[i] = s * M_1_PI / (hi * hi * hi)


@njit(cache=True, parallel=True)
def _force_loop(iflav, pos, mass, h, rho, pP, indptr, indices,
                rev_indptr, rev_indices, one_sided,
                box, periodic, icvel, q0d, e0, diag):
    """ICGenerator pair force, gather-only, over CSR neighbour lists.

    own-ball side: for each j in indices[indptr[i]:indptr[i+1]] evaluate the
    h_i kernel (i is 'p'). The h_j side (i is 'q', owner j's kernel) is then
    supplied either
      - one_sided=False: from the transposed graph rev_indices[i] (the j's
        whose own ball contains i) - the exact arithmetic-mean symmetrization;
      - one_sided=True : inline from i's OWN list (same j's, kernel h_j), no
        transpose - SPH-EXA style; misses only pairs where i sits in j's ball
        but j is beyond i's list (h_j/h_i large; kernel-edge terms).
    Density is independent of this choice (it is computed separately from the
    own list); only the FORCE differs. Each prange thread writes only icvel[i]
    (+ q0d[i]/e0[i] when diag) - no races."""
    n = pos.shape[0]
    bx, by, bz = box[0], box[1], box[2]
    hbx, hby, hbz = 0.5 * bx, 0.5 * by, 0.5 * bz
    for i in prange(n):
        hi = h[i]
        ih2 = 1.0 / (hi * hi)
        fnorm1_i = 0.5 * M_1_PI * ih2 * ih2 / hi  # 0.5/(pi h^5)
        ax = 0.0
        ay = 0.0
        az = 0.0
        q0 = 0.0
        e0x = 0.0
        e0y = 0.0
        e0z = 0.0
        # own-ball side: i is "p", kernel h_i (+ inline h_j side if one_sided)
        for e in range(indptr[i], indptr[i + 1]):
            j = indices[e]
            dx = pos[i, 0] - pos[j, 0]
            dy = pos[i, 1] - pos[j, 1]
            dz = pos[i, 2] - pos[j, 2]
            if periodic:
                if dx > hbx:
                    dx -= bx
                elif dx < -hbx:
                    dx += bx
                if dy > hby:
                    dy -= by
                elif dy < -hby:
                    dy += by
                if dz > hbz:
                    dz -= bz
                elif dz < -hbz:
                    dz += bz
            d2 = dx * dx + dy * dy + dz * dz
            ar2 = d2 * ih2
            if ar2 < 4.0:
                rq = dw_wc2(ar2) * fnorm1_i * mass[j]
                ig_i, ig_j = ig_pair(iflav, pP[i], pP[j], rho[i], rho[j])
                s = (ig_i + ig_j) * rq
                ax -= s * dx
                ay -= s * dy
                az -= s * dz
                if diag:
                    q0 += mass[j] / rho[j] * w_wc2(ar2)
                    g_i, g_j = ig_pair(iflav, 1.0, 1.0, rho[i], rho[j])
                    se = (g_i + g_j) * rq
                    e0x += se * dx
                    e0y += se * dy
                    e0z += se * dz
            if one_sided:
                # neighbor-h side from i's own list (the term the symmetric
                # walk delivers via the transposed graph)
                hj = h[j]
                ihj2 = 1.0 / (hj * hj)
                ar2j = d2 * ihj2
                if ar2j < 4.0:
                    rp = dw_wc2(ar2j) * (0.5 * M_1_PI * ihj2 * ihj2 / hj) \
                        * mass[j]
                    ig_j, ig_i = ig_pair(iflav, pP[j], pP[i], rho[j], rho[i])
                    s = (ig_j + ig_i) * rp
                    ax -= s * dx
                    ay -= s * dy
                    az -= s * dz
                    if diag:
                        g_j, g_i = ig_pair(iflav, 1.0, 1.0, rho[j], rho[i])
                        se = (g_j + g_i) * rp
                        e0x += se * dx
                        e0y += se * dy
                        e0z += se * dz
        # transposed side: i is "q", ball owner j is "p", kernel h_j
        if not one_sided:
            for ptr in range(rev_indptr[i], rev_indptr[i + 1]):
                j = rev_indices[ptr]
                hj = h[j]
                ihj2 = 1.0 / (hj * hj)
                dx = pos[i, 0] - pos[j, 0]
                dy = pos[i, 1] - pos[j, 1]
                dz = pos[i, 2] - pos[j, 2]
                if periodic:
                    if dx > hbx:
                        dx -= bx
                    elif dx < -hbx:
                        dx += bx
                    if dy > hby:
                        dy -= by
                    elif dy < -hby:
                        dy += by
                    if dz > hbz:
                        dz -= bz
                    elif dz < -hbz:
                        dz += bz
                ar2 = (dx * dx + dy * dy + dz * dz) * ihj2
                if ar2 >= 4.0:
                    continue
                rp = dw_wc2(ar2) * (0.5 * M_1_PI * ihj2 * ihj2 / hj) * mass[j]
                ig_j, ig_i = ig_pair(iflav, pP[j], pP[i], rho[j], rho[i])
                s = (ig_j + ig_i) * rp
                # C: ICvel[q] += (...) * rp * dx_pq with dx_pq = r_p - r_q = -dx
                ax -= s * dx
                ay -= s * dy
                az -= s * dz
                if diag:
                    g_j, g_i = ig_pair(iflav, 1.0, 1.0, rho[j], rho[i])
                    se = (g_j + g_i) * rp
                    e0x += se * dx
                    e0y += se * dy
                    e0z += se * dz
        icvel[i, 0] = ax
        icvel[i, 1] = ay
        icvel[i, 2] = az
        if diag:
            hi3 = hi * hi * hi
            q0d[i] = (q0 + mass[i] / rho[i] * W0) * M_1_PI / hi3
            e0[i, 0] = rho[i] * hi * e0x
            e0[i, 1] = rho[i] * hi * e0y
            e0[i, 2] = rho[i] * hi * e0z


def get_force_loop_csr(flavour, one_sided=False):
    """CSR (ball-gather) force loop bound to a flavour id and one_sided flag."""
    if flavour not in FLAVOURS:
        raise ValueError(f"unknown flavour '{flavour}', have {list(FLAVOURS)}")
    iflav = FLAVOURS[flavour]
    osf = one_sided

    def force_loop(pos, mass, h, rho, pP, indptr, indices, rev_indptr,
                   rev_indices, box, periodic, icvel, q0d, e0, diag):
        return _force_loop(iflav, pos, mass, h, rho, pP, indptr, indices,
                           rev_indptr, rev_indices, osf, box, periodic,
                           icvel, q0d, e0, diag)
    return force_loop


def get_force_loop(flavour, one_sided=False):
    """kNN (dense) force loop bound to a flavour id and one_sided flag. The
    dense (N, k) candidate array is presented to the shared CSR _force_loop as
    indptr=arange*k, indices=idx.ravel() (a view; identical iteration order, so
    results are bit-identical to a per-particle range(k) walk). Adding a flavour
    = new id in flavours.FLAVOURS + branch in flavours.ig_pair."""
    if flavour not in FLAVOURS:
        raise ValueError(f"unknown flavour '{flavour}', have {list(FLAVOURS)}")
    iflav = FLAVOURS[flavour]
    osf = one_sided
    state = {}

    def force_loop(pos, mass, h, rho, pP, idx, rev_indptr, rev_indices,
                   box, periodic, icvel, q0d, e0, diag):
        n, k = idx.shape
        if state.get('key') != (n, k):
            state['indptr'] = np.arange(0, (n + 1) * k, k, dtype=np.int64)
            state['key'] = (n, k)
        indices = idx.reshape(-1)
        return _force_loop(iflav, pos, mass, h, rho, pP, state['indptr'],
                           indices, rev_indptr, rev_indices, osf, box,
                           periodic, icvel, q0d, e0, diag)
    return force_loop
