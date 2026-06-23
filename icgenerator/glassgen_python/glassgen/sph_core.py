"""Numba pair loops: SPH density, ICGenerator force, quality diagnostics.

Ports the ICGEN pair physics of icgenerator/gasoline:
- density: gather sum rho_i = 1/(pi h_i^3) sum_j m_j W(q_ij) including self
  (Density/DenDVDX, smoothfcn.c).
- force (ICGenerator, smoothfcn.c:5794): in C each mutual pair is visited
  twice by the symmetric resmooth (once with each owner's h) and the kernel
  carries fNorm1 = 0.5/(pi h^5), which yields the arithmetic-mean kernel
  symmetrization. Here that is reformulated gather-only: particle i
  accumulates its own-ball terms (kernel h_i, neighbors idx[i]) plus the
  transposed-graph terms (kernel h_j, for each j whose ball contains i), so
  each prange thread writes only icvel[i] - no races.
- diagnostics (CHECKQUALITY): Q0_i = sum m_j/rho_j W_ij(h_i) (partition of
  unity, -> 1; the zeroth-order interpolation consistency - see
  diagnostics.py for the full Q0/Q1/E0/E1 hierarchy) and E0_i = zero-order
  force error (the force with pP == 1, scaled by rho_i h_i), accumulated in
  the same loop.
"""
import math

import numpy as np
from numba import njit, prange

from .flavours import FLAVOURS, ig_pair
from .kernels import w_wc2, dw_wc2, W0

M_1_PI = 1.0 / math.pi


@njit(cache=True, parallel=True)
def density_loop(pos, mass, h, idx, box, periodic, rho):
    # positions must be wrapped to [-L/2, L/2] when periodic: min-image is
    # then two predictable compares instead of a divide+rint per axis
    n, k = idx.shape
    bx, by, bz = box[0], box[1], box[2]
    hbx, hby, hbz = 0.5 * bx, 0.5 * by, 0.5 * bz
    for i in prange(n):
        hi = h[i]
        ih2 = 1.0 / (hi * hi)
        s = mass[i] * W0  # self contribution, W(0) = W0
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
def _force_loop(iflav, pos, mass, h, rho, pP, idx, rev_indptr, rev_indices,
                box, periodic, icvel, q0d, e0, diag):
    n, k = idx.shape
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
        # own-ball side: i is "p", kernel h_i
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
            if ar2 >= 4.0:
                continue
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
        # transposed side: i is "q", ball owner j is "p", kernel h_j
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


@njit(cache=True, parallel=True)
def _force_loop_oneside(iflav, pos, mass, h, rho, pP, idx,
                        box, periodic, icvel, q0d, e0, diag):
    """SPH-EXA-style one-sided pair walk: both kernel sides (h_i and h_j)
    are evaluated from particle i's own candidate list, with no transposed
    graph. Differs from the symmetric walk only for pairs beyond i's list
    yet inside j's kernel, which needs h_j/h_i > (k/nsmooth)^(1/3) across
    one kernel radius - absent in near-uniform regions, and the missed
    terms sit at the kernel edge where dW -> 0."""
    n, k = idx.shape
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
            d2 = dx * dx + dy * dy + dz * dz
            # own-h side
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
            # neighbor-h side (the term the symmetric walk delivers via the
            # transposed graph)
            hj = h[j]
            ihj2 = 1.0 / (hj * hj)
            ar2j = d2 * ihj2
            if ar2j < 4.0:
                rp = dw_wc2(ar2j) * (0.5 * M_1_PI * ihj2 * ihj2 / hj) * mass[j]
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
        icvel[i, 0] = ax
        icvel[i, 1] = ay
        icvel[i, 2] = az
        if diag:
            hi3 = hi * hi * hi
            q0d[i] = (q0 + mass[i] / rho[i] * W0) * M_1_PI / hi3
            e0[i, 0] = rho[i] * hi * e0x
            e0[i, 1] = rho[i] * hi * e0y
            e0[i, 2] = rho[i] * hi * e0z


@njit(cache=True, parallel=True)
def density_loop_csr(pos, mass, h, indptr, indices, box, periodic, rho):
    """SPH density over a CSR (ball-gather) neighbour list: each particle i
    sums over indices[indptr[i]:indptr[i+1]] (all within ~2*h_i)."""
    n = pos.shape[0]
    bx, by, bz = box[0], box[1], box[2]
    hbx, hby, hbz = 0.5 * bx, 0.5 * by, 0.5 * bz
    for i in prange(n):
        hi = h[i]
        ih2 = 1.0 / (hi * hi)
        s = mass[i] * W0
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
def _force_loop_csr(iflav, pos, mass, h, rho, pP, indptr, indices,
                    rev_indptr, rev_indices, box, periodic, icvel, q0d, e0,
                    diag):
    """Symmetric ICGEN force on CSR ball lists: own-ball side over
    indices[i] (kernel h_i), transposed side over rev_indices[i] (kernel
    h_j, j's whose ball contains i). Mirrors _force_loop with CSR own side."""
    n = pos.shape[0]
    bx, by, bz = box[0], box[1], box[2]
    hbx, hby, hbz = 0.5 * bx, 0.5 * by, 0.5 * bz
    for i in prange(n):
        hi = h[i]
        ih2 = 1.0 / (hi * hi)
        fnorm1_i = 0.5 * M_1_PI * ih2 * ih2 / hi
        ax = 0.0
        ay = 0.0
        az = 0.0
        q0 = 0.0
        e0x = 0.0
        e0y = 0.0
        e0z = 0.0
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
            if ar2 >= 4.0:
                continue
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


@njit(cache=True, parallel=True)
def _force_loop_csr_oneside(iflav, pos, mass, h, rho, pP, indptr, indices,
                            box, periodic, icvel, q0d, e0, diag):
    """One-sided CSR force: both kernel sides evaluated from i's own ball
    list, no transpose (see _force_loop_oneside)."""
    n = pos.shape[0]
    bx, by, bz = box[0], box[1], box[2]
    hbx, hby, hbz = 0.5 * bx, 0.5 * by, 0.5 * bz
    for i in prange(n):
        hi = h[i]
        ih2 = 1.0 / (hi * hi)
        fnorm1_i = 0.5 * M_1_PI * ih2 * ih2 / hi
        ax = 0.0
        ay = 0.0
        az = 0.0
        q0 = 0.0
        e0x = 0.0
        e0y = 0.0
        e0z = 0.0
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
            hj = h[j]
            ihj2 = 1.0 / (hj * hj)
            ar2j = d2 * ihj2
            if ar2j < 4.0:
                rp = dw_wc2(ar2j) * (0.5 * M_1_PI * ihj2 * ihj2 / hj) * mass[j]
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
    """CSR (ball-gather) force loop bound to a flavour id."""
    if flavour not in FLAVOURS:
        raise ValueError(f"unknown flavour '{flavour}', have {list(FLAVOURS)}")
    iflav = FLAVOURS[flavour]
    if one_sided:
        def force_loop(pos, mass, h, rho, pP, indptr, indices, rev_indptr,
                       rev_indices, box, periodic, icvel, q0d, e0, diag):
            return _force_loop_csr_oneside(iflav, pos, mass, h, rho, pP,
                                           indptr, indices, box, periodic,
                                           icvel, q0d, e0, diag)
    else:
        def force_loop(pos, mass, h, rho, pP, indptr, indices, rev_indptr,
                       rev_indices, box, periodic, icvel, q0d, e0, diag):
            return _force_loop_csr(iflav, pos, mass, h, rho, pP, indptr,
                                   indices, rev_indptr, rev_indices, box,
                                   periodic, icvel, q0d, e0, diag)
    return force_loop


def get_force_loop(flavour, one_sided=False):
    """Bind a flavour id onto the cached force loop. Adding a flavour = new
    id in flavours.FLAVOURS + branch in flavours.ig_pair; future ISPH
    flavours will use the reserved cmat slot. one_sided selects the
    transpose-free walk (see _force_loop_oneside)."""
    if flavour not in FLAVOURS:
        raise ValueError(f"unknown flavour '{flavour}', have {list(FLAVOURS)}")
    iflav = FLAVOURS[flavour]

    if one_sided:
        def force_loop(pos, mass, h, rho, pP, idx, rev_indptr, rev_indices,
                       box, periodic, icvel, q0d, e0, diag):
            return _force_loop_oneside(iflav, pos, mass, h, rho, pP, idx,
                                       box, periodic, icvel, q0d, e0, diag)
    else:
        def force_loop(*args):
            return _force_loop(iflav, *args)

    return force_loop
