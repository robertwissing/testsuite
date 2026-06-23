"""Displacement normalization and the adaptive step-size state machine.

normalize_displacements ports pkdICGetFmax / pkdICNormalizeFmax /
pkdICNormalizeVel (pkd.c): per-particle force magnitudes are log-rescaled
relative to the global maximum (T = 0.1) and the vectors normalized so the
largest component moves speednorm in [1/maxspeed, 1] (maxspeed = 1000);
the displacement applied later is icvel * ICR0 * h.

StepController ports msrICGetFmax (master.c) and the restart logic in
msrCalcEandL: ICR0 shrinks on "bad" moves (rms force grew) and grows
slightly on good ones; once ICR0 < nfac (~1/20 interparticle spacing) it
re-anneals to ICRLOOP0 while the density error keeps improving, halving
ICRLOOP0 each restart, and enters a gentle polish phase (icfmax = -1) when
ICRLOOP0 < 2*nfac. Constants match the C code.
"""
import math
from collections import deque

import numpy as np
from numba import njit, prange

T_LOG = 0.1
MAXSPEED = 1000.0


@njit(cache=True, parallel=True)
def normalize_displacements(icvel):
    """In-place normalization. Returns (frms, fmax): frms = sqrt(sum
    |icvel|^2 / N) measured before normalization, i.e. the *true* per-particle
    RMS residual pseudo-pressure force - same per-particle normalization as the
    mean/p95/max force diagnostics (diagnostics.residual_force), so the logged
    frms is directly comparable to those (RMS sits just above the mean). NB the
    C ICvelAvg is sqrt(sum)/N = this RMS / sqrt(N); we drop that extra 1/sqrt(N)
    here because frms's absolute scale is never used - the step controller only
    reads the SIGN of its change (update_step) and the run stops on icr0, not
    frms - so this is behaviour-preserving and only fixes the displayed value.
    frms drives the step controller. NB this is the *force-equilibrium* descent
    signal, NOT
    a partition/consistency diagnostic: it is distinct from the Q0
    partition-of-unity (q0d, -> 1), the Q1 first-moment error and the E0/E1
    gradient errors in diagnostics.py. frms -> 0 as the glass reaches force
    equilibrium. fmax = global max force component."""
    n = icvel.shape[0]
    fmaxpart = np.empty(n)
    frmstot = 0.0
    fmax = 0.0
    for i in prange(n):
        vx = icvel[i, 0]
        vy = icvel[i, 1]
        vz = icvel[i, 2]
        frmstot += vx * vx + vy * vy + vz * vz
        fm = max(abs(vx), max(abs(vy), abs(vz)))
        fmaxpart[i] = fm
        fmax = max(fmax, fm)
    if fmax <= 0.0:
        # no forces at all: leave displacements zero (C never hits this)
        return 0.0, 0.0
    lognorm = math.log(fmax * T_LOG + 1.0)
    for i in prange(n):
        icfmax = math.log(fmaxpart[i] * T_LOG + 1.0) / lognorm
        speednorm = (1.0 - 1.0 / MAXSPEED) * icfmax + 1.0 / MAXSPEED
        if fmaxpart[i] > 0.0:
            f = speednorm / fmaxpart[i]
            icvel[i, 0] *= f
            icvel[i, 1] *= f
            icvel[i, 2] *= f
        else:
            icvel[i, 0] = 1.0
            icvel[i, 1] = 1.0
            icvel[i, 2] = 1.0
    return math.sqrt(frmstot / n), fmax


@njit(cache=True, parallel=True)
def apply_step(pos, icvel, h, icr0, box, periodic, disp, max_step_frac,
               disp_prev, momentum):
    """Fused move: pos += icvel * icr0 * h, periodic wrap, per-particle
    displacement accumulation (disp[i] += step, for the per-particle
    Verlet-skin staleness test) and the max-step reduction in one pass.

    max_step_frac : 0.0 disables the cap (the coarse / single-relax default):
        early-sweep steps can exceed a box length, hence the while-wrap. When
        > 0 (the fine, confined phase) each component is clamped to
        max_step_frac * h[i] so a particle moves at most that fraction of its
        own smoothing length per step - this both forbids cross-box jumps and
        keeps the relaxation local. With the cap a step is <= h << L/2, so a
        single image correction (if/elif) is always sufficient (no while).

    momentum : 0.0 disables it (byte-identical to the memoryless steepest
        descent move). When > 0 a heavy-ball term gives the descent memory:
        d_i = momentum*d_prev_i + icvel_i ; x_i += d_i*icr0*h ; d_prev_i = d_i.
        icvel is ALREADY the per-particle normalised descent direction (post
        normalize_displacements), so momentum blends normalised directions:
        the cross-valley zig-zag cancels while the along-valley push accumulates,
        improving the slow frms tail (~kappa -> ~sqrt(kappa)). disp_prev is the
        per-particle velocity workspace (input-order, permuted alongside pos on
        reorder; NOT reset on neighbour rebuild - positions are continuous). The
        blended direction is stored UNCLAMPED (bounded by 1/(1-momentum) since
        |icvel| <= 1); the confine cap clamps only the applied position move."""
    n = pos.shape[0]
    bx, by, bz = box[0], box[1], box[2]
    hbx, hby, hbz = 0.5 * bx, 0.5 * by, 0.5 * bz
    capped = max_step_frac > 0.0
    use_mom = momentum > 0.0
    maxstep = 0.0
    for i in prange(n):
        s = icr0 * h[i]
        if use_mom:
            dx = momentum * disp_prev[i, 0] + icvel[i, 0]
            dy = momentum * disp_prev[i, 1] + icvel[i, 1]
            dz = momentum * disp_prev[i, 2] + icvel[i, 2]
            disp_prev[i, 0] = dx
            disp_prev[i, 1] = dy
            disp_prev[i, 2] = dz
            sx = dx * s
            sy = dy * s
            sz = dz * s
        else:
            sx = icvel[i, 0] * s
            sy = icvel[i, 1] * s
            sz = icvel[i, 2] * s
        if capped:
            cap = max_step_frac * h[i]
            sx = min(cap, max(-cap, sx))
            sy = min(cap, max(-cap, sy))
            sz = min(cap, max(-cap, sz))
        m = max(abs(sx), max(abs(sy), abs(sz)))
        disp[i] += m
        maxstep = max(maxstep, m)
        x = pos[i, 0] + sx
        y = pos[i, 1] + sy
        z = pos[i, 2] + sz
        if periodic and capped:
            # capped step < L/2: at most one image boundary crossed
            if x > hbx:
                x -= bx
            elif x < -hbx:
                x += bx
            if y > hby:
                y -= by
            elif y < -hby:
                y += by
            if z > hbz:
                z -= bz
            elif z < -hbz:
                z += bz
        elif periodic:
            while x > hbx:
                x -= bx
            while x < -hbx:
                x += bx
            while y > hby:
                y -= by
            while y < -hby:
                y += by
            while z > hbz:
                z -= bz
            while z < -hbz:
                z += bz
        pos[i, 0] = x
        pos[i, 1] = y
        pos[i, 2] = z
    return maxstep


@njit(cache=True, parallel=True)
def pseudo_pressure(rho, rho0, rhopow, pP):
    """Fused pP_i = (rho_i/rho0_i)^rhopow and mean |1 - rho/rho0| return."""
    n = rho.shape[0]
    err = 0.0
    for i in prange(n):
        r = rho[i] / rho0[i]
        if rhopow == 2.0:
            pP[i] = r * r
        else:
            pP[i] = r ** rhopow
        err += abs(1.0 - r)
    return err / n


class StepController:
    def __init__(self, icr0, icrloop0, icr0rate, nsmooth, box=None,
                 reanneal=True, polish_down=0.1, polish_up=0.002):
        # master.c: if bPeriodic, ICR0 starts at 0.99*maxPeriod*dICR0
        if box is not None:
            icr0 = 0.99 * float(np.max(box)) * icr0
        self.periodic = box is not None
        # reanneal : when False the ICRLOOP0 re-anneal/halve ladder is skipped
        # (the fine progressive phase - large-scale errors are already removed
        # by the coarse stage, so icr0 just decays through sweep into polish).
        # Polish entry + the icr0_stop finish still work; only the ladder that
        # re-inflates icr0 back to ICRLOOP0 is disabled. Default True = the
        # validated single-relax behaviour.
        self.reanneal = reanneal
        # polish-phase icr0 adaptation coefficients (scaled by icr0rate):
        # shrink factor 1/(1+polish_down*rate) when frms rises, grow factor
        # (1+polish_up*rate) when frms falls. Defaults 0.1/0.002 = the C
        # values (~46:1 shrink:grow asymmetry). Exposed to tune how fast frms
        # is driven down in the refine region (a gentler shrink / stronger grow
        # keeps icr0 mobile). The sweep-phase pair stays the C 0.8/0.16.
        self.polish_down = polish_down
        self.polish_up = polish_up
        self.icr0 = icr0
        self.icrloop0 = icrloop0
        self.rate = icr0rate
        self.nfac = 0.025 * (1.0 / nsmooth) ** (1.0 / 3.0)
        self.icfmax = 1e13  # dICFmax default; -1 marks the polish phase
        # prev RMS residual force (the C dICQ1Avg / ICvelAvg force-equilibrium
        # criterion - NOT the Q1 first-moment consistency in diagnostics.py)
        self.frms_prev = 0.0
        self.step_zero = True

    @property
    def polish(self):
        return self.icfmax <= 0.0

    def update_step(self, frms, h_avg):
        """Port of msrICGetFmax: called each iteration with the
        pre-normalization RMS residual force frms and mean smoothing length."""
        if self.polish:
            rate = 1.0 + self.polish_down * self.rate
            rateup = 1.0 + self.polish_up * self.rate
        else:
            rate = 1.0 + 0.8 * self.rate
            rateup = 1.0 + 0.16 * self.rate
        if self.step_zero:
            # periodic: amplify so the first sweeps cross the whole box
            # (displacements wrap). In an open domain that diverges, so
            # icr0 stays in units of h there (deviation from C, which is
            # only ever used with periodic boxes).
            if self.periodic:
                self.icr0 = self.icr0 / h_avg
            self.icrloop0 = self.icr0 * self.icrloop0
            self.step_zero = False
        if self.frms_prev - frms < 0.0:
            self.icr0 /= rate
        else:
            self.icr0 *= rateup
        self.frms_prev = frms

    def update_restart(self, dens_err):
        """Port of the ICGEN block in msrCalcEandL: called once per
        iteration with the mean relative density error <|1 - rho/rho0|>."""
        if not self.reanneal:
            # No re-anneal ladder: icrloop0 is never halved, so the usual
            # icrloop0 < 2*nfac polish trigger can't fire. Instead drop into
            # polish as soon as icr0 has decayed below nfac (fine-scale force
            # equilibrium); the icr0_stop finish then ends the run as normal.
            if self.icr0 < self.nfac:
                self.icfmax = -1.0
            return
        if self.icrloop0 < 2.0 * self.nfac:
            self.icfmax = -1.0
        if self.icr0 < self.nfac and self.icfmax > 0.0:
            if self.icfmax > dens_err:
                self.icr0 = self.icrloop0
                self.icfmax = dens_err
            elif self.icrloop0 < 2.0 * self.nfac:
                self.icr0 = 2.0 * self.nfac
                self.icfmax = -1.0
            else:
                self.icrloop0 /= 2.0
                self.icr0 = self.icrloop0
                self.icfmax = dens_err


class Icr0Stop:
    """Absolute step-size finish - stop when the maximum step has shrunk to a
    negligible fraction of the smoothing length.

    icr0 IS the per-iteration maximum displacement in units of h: the move is
    icvel*icr0*h and icvel is normalised so its largest component is
    speednorm_i <= 1, so the most-active particle moves at most icr0*h - i.e.
    icr0 is the largest step as a FRACTION of h. That makes it resolution-,
    gradient- and particle-count INVARIANT (h scales out), unlike the residual
    force frms (which scales with the per-particle mass/h and the density
    gradient and so can't carry a fixed absolute threshold across cases). So
    'stop once icr0 < threshold' = 'stop once even the most-active particle
    moves less than threshold*h per step' - a portable convergence criterion.

    In the confined refine stage icr0 self-regulates to the critical-step
    SETPOINT (~2.8-3e-3*h, slowly declining as the glass stiffens - see the
    rate-asymmetry analysis in PLAN task #6d), so threshold ~3e-3 stops right as
    icr0 settles onto / drops below that setpoint, i.e. when the glass has
    reached its force-equilibrium stiffness.

    The controller sawtooths icr0 every iteration (grow/shrink on the sign of
    the frms change), so a BARE icr0 < threshold crossing is a noisy, one-sample
    event (it can dip below on one iter and bounce back the next - e.g. momentum
    beta=0.3 keeps icr0 mobile and oscillating around the setpoint). A trailing-
    median smooth over `window` removes the sawtooth so the crossing is
    reproducible. Feed only polish-phase samples (relax does)."""

    def __init__(self, threshold=3e-3, window=10):
        self.threshold = threshold
        self.window = window
        self.buf = deque(maxlen=window)

    def update(self, icr0):
        """Feed one icr0 sample; return True once the smoothed icr0 <= threshold."""
        self.buf.append(icr0)
        if len(self.buf) < self.window:
            return False
        return float(np.median(self.buf)) <= self.threshold
