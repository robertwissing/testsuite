"""Operator-matched SPH consistency diagnostics (diagnostic-only).

These quantify how consistent the *actual ICGenerator force operator* is on a
glass - they reuse the same Wendland C2 kernel, gasoline normalization (fNorm = 1/pi h^3,
fNorm1 = 1/pi h^5) and the two-sided symmetrized gradient that sph_core uses, so
the numbers describe the operator that relaxes the glass, not an idealized one.

Hierarchy (V_j = m_j/rho_j):

  Q0_i = fNorm_i * sum_j V_j W(q_ij; h_i)                          -> 1
         zeroth-order interpolation consistency (partition of unity).
  Q1_i = |fNorm_i * sum_j V_j (r_j - r_i) W(q_ij; h_i)| / h_i      -> 0
         first-order interpolation consistency (first kernel moment).
  E0_i = |sum_pairs (g_i+g_j) rq dx| * rho_i h_i                   -> 0
         zeroth-order gradient/force error: the residual force with uniform
         pseudo-pressure (pP == 1). With gdforce + pP=1 this reduces to
         h_i |sum_j V_j grad_i W|, the textbook 0th-order gradient error.
  E1_i = ||M_i - I||_F                                             -> 0
         first-order gradient consistency: M_i = -rho_i sum_pairs (g_i+g_j) rq
         (dx (x) dx) is the renormalization matrix sum_j V_j (r_j-r_i)(x)grad_iW,
         which -> I for a linearly-consistent gradient.

Q0/Q1/E1 are gather quantities (kernel h_i over particle i's own neighbour
ball). E0 is two-sided, summing i's own-ball terms (kernel h_i) and the
transposed-graph terms (the j whose ball contains i, kernel h_j) - exactly the
gather-only reformulation of the symmetric resmooth in sph_core._force_loop,
so E0 is the operator's actual force residual. E1 (the renormalization matrix)
is single-kernel by construction: two-siding it double-counts each symmetric
pair and drives M -> 2I instead of I. Input neighbour lists are the ball-gather
CSR (own: indptr/indices, transposed: tindptr/tindices).
"""
import math

import numpy as np
from numba import njit, prange

from .flavours import ig_pair
from .kernels import dw_wc2, w_wc2, W0

M_1_PI = 1.0 / math.pi


@njit(cache=True, parallel=True)
def _consistency(iflav, pos, mass, h, rho, indptr, indices, tindptr, tindices,
                 box, periodic, q0, q1, e0, e1):
    n = pos.shape[0]
    bx, by, bz = box[0], box[1], box[2]
    hbx, hby, hbz = 0.5 * bx, 0.5 * by, 0.5 * bz
    for i in prange(n):
        hi = h[i]
        ih2 = 1.0 / (hi * hi)
        fnorm_i = M_1_PI * ih2 / hi          # 1/(pi h_i^3)
        fnorm1_i = 0.5 * fnorm_i * ih2        # 0.5/(pi h_i^5) (per-side weight)
        # Q0/Q1: single-kernel gather (own ball), self term included
        q0s = mass[i] / rho[i] * W0           # V_i W(0) = V_i W0
        q1x = 0.0
        q1y = 0.0
        q1z = 0.0
        # E0/E1: two-sided operator residual (pP == 1)
        e0x = 0.0
        e0y = 0.0
        e0z = 0.0
        m00 = 0.0; m01 = 0.0; m02 = 0.0
        m10 = 0.0; m11 = 0.0; m12 = 0.0
        m20 = 0.0; m21 = 0.0; m22 = 0.0
        # --- own-ball side: kernel h_i ---
        for ptr in range(indptr[i], indptr[i + 1]):
            j = indices[ptr]
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
            vj = mass[j] / rho[j]
            wij = w_wc2(ar2)
            q0s += vj * wij
            # first moment uses (r_j - r_i) = -dx
            q1x -= vj * dx * wij
            q1y -= vj * dy * wij
            q1z -= vj * dz * wij
            rq = dw_wc2(ar2) * fnorm1_i * mass[j]
            g_i, g_j = ig_pair(iflav, 1.0, 1.0, rho[i], rho[j])
            se = (g_i + g_j) * rq
            e0x += se * dx
            e0y += se * dy
            e0z += se * dz
            # renorm matrix contribution: -rho_i * se * (dx (x) dx); the rho_i
            # and final sign are applied after both sides are summed
            m00 += se * dx * dx; m01 += se * dx * dy; m02 += se * dx * dz
            m10 += se * dy * dx; m11 += se * dy * dy; m12 += se * dy * dz
            m20 += se * dz * dx; m21 += se * dz * dy; m22 += se * dz * dz
        # --- transposed side: kernel h_j (only E0/E1, the force operator) ---
        for ptr in range(tindptr[i], tindptr[i + 1]):
            j = tindices[ptr]
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
            g_j, g_i = ig_pair(iflav, 1.0, 1.0, rho[j], rho[i])
            se = (g_j + g_i) * rp
            e0x += se * dx
            e0y += se * dy
            e0z += se * dz
            # NB: the renorm matrix M (-> E1) is a single-kernel gather
            # quantity (Sum_j V_j (r_j-r_i)(x)grad_i W with h_i, -> I); it is
            # NOT summed on the transposed side - doing so double-counts each
            # symmetric pair and drives M -> 2I. Only E0 (the actual symmetric
            # force residual) is genuinely two-sided.
        q0[i] = q0s * fnorm_i
        q1[i] = fnorm_i * math.sqrt(q1x * q1x + q1y * q1y + q1z * q1z) / hi
        e0[i] = rho[i] * hi * math.sqrt(e0x * e0x + e0y * e0y + e0z * e0z)
        # M = -rho_i * sum se (dx (x) dx); first-order error = ||M - I||_F
        r = rho[i]
        a00 = -r * m00 - 1.0; a11 = -r * m11 - 1.0; a22 = -r * m22 - 1.0
        a01 = -r * m01; a02 = -r * m02; a12 = -r * m12
        a10 = -r * m10; a20 = -r * m20; a21 = -r * m21
        e1[i] = math.sqrt(a00 * a00 + a11 * a11 + a22 * a22
                          + a01 * a01 + a02 * a02 + a12 * a12
                          + a10 * a10 + a20 * a20 + a21 * a21)


def residual_force(pos, mass, h, rho, rho0, box, rhopow=2.0, nsmooth=64,
                   flavour='gdforce', nsmoothmin=32, uniform_pP=False):
    """Per-particle magnitude of the actual residual pseudo-pressure force on
    the glass - the descent force the relaxation moves along and whose RMS
    (frms) drives the force-equilibrium convergence / stop. This is the SAME
    two-sided ICGenerator force operator as sph_core (built on the ball-gather
    graph, exactly like E0 in operator_consistency), but with the REAL
    pP = (rho/rho0)^rhopow rather than the uniform pP that defines E0 - so it
    -> 0 only as the glass reaches the *target* density, not merely operator
    self-consistency. pos must be wrapped to [-L/2, L/2] when periodic; pass
    the force h (density-based) to match what the relaxation used.

    uniform_pP=True replaces pP with ones, returning the raw |icvel| of the
    UNIFORM-pP operator force - the irreducible operator noise floor of frms in
    the same units (the E0 force without the rho*h normalization). frms can
    never fall below the RMS of this, since at the converged glass rho->rho0 so
    pP->1 and the actual force reduces to exactly this floor."""
    from .ball import BallNeighborList
    from .flavours import FLAVOURS
    from .sph_core import get_force_loop_csr
    if nsmoothmin >= nsmooth:        # keep the floor below nsmooth (see relax)
        nsmoothmin = nsmooth - 1
    periodic = box is not None
    boxarr = np.asarray(box if periodic else (1e30, 1e30, 1e30),
                        dtype=np.float64)
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    mass = np.ascontiguousarray(mass, dtype=np.float64)
    h = np.ascontiguousarray(h, dtype=np.float64)
    rho = np.ascontiguousarray(rho, dtype=np.float64)
    bn = BallNeighborList(boxarr, nsmoothmin=nsmoothmin, skin=0.0,
                          build_transpose=True)
    bn._build(pos, h)
    n = len(pos)
    if uniform_pP:
        pP = np.ones(n)
    else:
        pP = (rho / np.asarray(rho0, dtype=np.float64)) ** rhopow
    icvel = np.empty((n, 3))
    dummy = np.empty((1, 1))
    fl = get_force_loop_csr(flavour, one_sided=False)
    fl(pos, mass, h, rho, pP, bn.indptr, bn.indices, bn.tip, bn.tii, boxarr,
       periodic, icvel, dummy[0], dummy, False)
    return np.sqrt((icvel * icvel).sum(1))


def operator_consistency(pos, mass, h, rho, box, nsmooth=64, flavour='gdforce',
                         nsmoothmin=32):
    """Build a ball neighbour graph on the glass and return per-particle
    (Q0, Q1, E0, E1) arrays for the operator-matched consistency hierarchy.
    pos must be wrapped to [-L/2, L/2] when periodic. Pass the *force* h
    (density-based) to match what the relaxation used."""
    from .ball import BallNeighborList
    from .flavours import FLAVOURS
    if nsmoothmin >= nsmooth:        # keep the floor below nsmooth (see relax)
        nsmoothmin = nsmooth - 1
    periodic = box is not None
    boxarr = np.asarray(box if periodic else (1e30, 1e30, 1e30), dtype=np.float64)
    iflav = FLAVOURS[flavour]
    bn = BallNeighborList(boxarr, nsmoothmin=nsmoothmin, skin=0.0,
                          build_transpose=True)
    bn._build(np.ascontiguousarray(pos), np.ascontiguousarray(h))
    n = len(pos)
    q0 = np.empty(n); q1 = np.empty(n); e0 = np.empty(n); e1 = np.empty(n)
    _consistency(iflav, pos, mass, h, rho, bn.indptr, bn.indices, bn.tip,
                 bn.tii, boxarr, periodic, q0, q1, e0, e1)
    return q0, q1, e0, e1


def plot_history(history, path, title=None, stop_reason=None):
    """Plot the relaxation convergence trajectory recorded in
    ``RelaxResult.history`` - the errors over iterations.

    history : list of (it, frms, dens_err, icr0, icrloop0, polish) tuples (what
        relax() records every iteration). For the progressive path this is the
        final (fine) stage. Left log-y axis carries the two error signals -
        ``frms`` (the residual-force RMS the run actually descends, the stop
        basis) and ``dens_err`` (mean |1-rho/rho0|); the right log-y axis
        carries ``icr0`` (the step size). The polish-entry iteration is marked
        (dashed), and the final iteration / stop_reason annotated - so one can
        read where the errors flatten vs where the icr0_stop finish fires.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    h = np.asarray(history, dtype=float)
    if h.size == 0:
        return
    it, frms, denserr, icr0, polish = h[:, 0], h[:, 1], h[:, 2], h[:, 3], h[:, 5]
    polish_start = it[polish > 0][0] if np.any(polish > 0) else None

    fig, axL = plt.subplots(figsize=(9, 5.5))
    axL.semilogy(it, frms, color='C3', lw=1.2, label='frms (residual force)')
    axL.semilogy(it, denserr, color='C0', lw=1.2,
                 label='denserr  <|1-rho/rho0|>')
    axL.set_xlabel('iteration')
    axL.set_ylabel('error (log)')
    axL.grid(True, which='both', alpha=0.25)

    axR = axL.twinx()
    axR.semilogy(it, icr0, color='C2', lw=1.0, alpha=0.7, label='icr0 (step)')
    axR.set_ylabel('icr0', color='C2')
    axR.tick_params(axis='y', labelcolor='C2')

    if polish_start is not None:
        axL.axvline(polish_start, color='0.5', ls='--', lw=1.0)
        axL.annotate(f'polish @ {int(polish_start)}',
                     xy=(polish_start, axL.get_ylim()[1]),
                     xytext=(3, -10), textcoords='offset points',
                     fontsize=8, color='0.4', rotation=90, va='top')
    stop_it = int(it[-1])
    axL.axvline(stop_it, color='k', ls=':', lw=1.0)
    lbl = f'stop @ {stop_it}' + (f'  ({stop_reason})' if stop_reason else '')
    axL.annotate(lbl, xy=(stop_it, axL.get_ylim()[0]),
                 xytext=(-3, 10), textcoords='offset points',
                 fontsize=8, ha='right')

    lines = axL.get_lines() + axR.get_lines()
    axL.legend(lines, [l.get_label() for l in lines], loc='upper right',
               fontsize=8)
    if title:
        axL.set_title(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=110, bbox_inches='tight')
    plt.close(fig)
    return path


def plot_error_history(tail_history, path, stat='p95', title=None,
                       stop_reason=None):
    """Plot a per-iteration error STAT trajectory from RelaxResult.tail_history
    (recorded when relax(diagnostics=True)). ``stat`` in {'mean','p95','max'}
    selects which order statistic of the per-particle density error and
    residual force to plot vs iteration. The p95/max are interface-dominated
    (the 360:1 shell), so they reveal whether the deep polish keeps tightening
    the TAIL after the mean denserr has flattened - the crux of where to stop.

    tail_history rows: (it, denserr_mean, denserr_p95, denserr_max, force_mean,
    force_p95, force_max, polish). Left log-y axis = density error, right
    log-y axis = residual force; polish entry + final/stop iteration marked.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    h = np.asarray(tail_history, dtype=float)
    if h.size == 0:
        return
    col = {'mean': (1, 4), 'p95': (2, 5), 'max': (3, 6)}[stat]
    it = h[:, 0]
    denserr = h[:, col[0]]
    force = h[:, col[1]]
    polish = h[:, 7]
    polish_start = it[polish > 0][0] if np.any(polish > 0) else None

    fig, axL = plt.subplots(figsize=(9, 5.5))
    axL.semilogy(it, denserr, color='C0', lw=1.2,
                 label=f'denserr ({stat})  |1-rho/rho0|')
    axL.set_xlabel('iteration')
    axL.set_ylabel(f'density error ({stat}, log)', color='C0')
    axL.tick_params(axis='y', labelcolor='C0')
    axL.grid(True, which='both', alpha=0.25)

    axR = axL.twinx()
    axR.semilogy(it, force, color='C3', lw=1.2,
                 label=f'force ({stat})  |icvel|')
    axR.set_ylabel(f'residual force ({stat}, log)', color='C3')
    axR.tick_params(axis='y', labelcolor='C3')

    if polish_start is not None:
        axL.axvline(polish_start, color='0.5', ls='--', lw=1.0)
        axL.annotate(f'polish @ {int(polish_start)}',
                     xy=(polish_start, axL.get_ylim()[1]),
                     xytext=(3, -10), textcoords='offset points',
                     fontsize=8, color='0.4', rotation=90, va='top')
    stop_it = int(it[-1])
    axL.axvline(stop_it, color='k', ls=':', lw=1.0)
    lbl = f'stop @ {stop_it}' + (f'  ({stop_reason})' if stop_reason else '')
    axL.annotate(lbl, xy=(stop_it, axL.get_ylim()[0]),
                 xytext=(-3, 10), textcoords='offset points',
                 fontsize=8, ha='right')

    lines = axL.get_lines() + axR.get_lines()
    axL.legend(lines, [l.get_label() for l in lines], loc='upper right',
               fontsize=8)
    axL.set_title(title or f'error tails ({stat})', fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=110, bbox_inches='tight')
    plt.close(fig)
    return path
