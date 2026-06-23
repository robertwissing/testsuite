"""Glass relaxation driver: the per-iteration loop of the C icgenerator.

Per iteration: neighbors/h -> SPH density -> target density -> pseudo
pressure pP = (rho/rho0)^rhopow -> pair forces -> displacement
normalization -> adaptive-step move x += icvel * ICR0 * h -> step control.
"""
import time
from dataclasses import dataclass, field

import numpy as np

from .ball import BallNeighborList
from .density import DensityField
from .neighbors import NeighborList, _gather3, cell_order
from .sph_core import (density_loop, density_loop_csr, get_force_loop,
                       get_force_loop_csr)
from .stepper import (Icr0Stop, StepController,
                      apply_step, normalize_displacements, pseudo_pressure)

_DUMMY = np.empty((1, 1))


@dataclass
class IterStats:
    it: int
    frms: float          # true per-particle RMS residual pseudo-pressure force -> 0
    dens_err: float      # <|1 - rho/rho0|> density error
    icr0: float
    icrloop0: float
    polish: bool
    pos: np.ndarray
    h: np.ndarray
    rho: np.ndarray
    rho0: np.ndarray
    q0: float = float('nan')   # mean partition of unity <sum m_j/rho_j W> -> 1 (zeroth-order)
    e0: float = float('nan')   # mean |zeroth-order force error| (diagnostics=True)


@dataclass
class RelaxResult:
    pos: np.ndarray
    h: np.ndarray
    rho: np.ndarray
    niter: int
    converged: bool
    history: list = field(default_factory=list)
    timings: dict = field(default_factory=dict)
    # per-iteration error tails (when diagnostics=True): (it, denserr_mean,
    # denserr_p95, denserr_max, force_mean, force_p95, force_max, polish)
    tail_history: list = field(default_factory=list)
    # which finish criterion fired: 'icr0_abs' (the median-smoothed icr0 < the
    # icr0_stop quality threshold - the default finish), 'polish_iter' (the fixed
    # polish-iteration budget), or 'max_iter' (the hard cap)
    stop_reason: str = 'max_iter'
    # per-stage breakdown for relax_progressive: one dict per level with
    # name/niter/converged/wall/n (empty for a plain single relax)
    levels: list = field(default_factory=list)
    # glass state captured at polish ENTRY (the end of the sweep phase, before
    # the gentle polish ladder runs) - lets the caller measure the pre-polish
    # diagnostic and see how much the polish phase changed it. None if polish
    # was never entered. niter at entry is polish_iter0.
    pos_prepolish: np.ndarray = None
    h_prepolish: np.ndarray = None
    rho_prepolish: np.ndarray = None
    polish_iter0: int = -1


def relax(pos, mass, density, box=None, nsmooth=64,
          flavour='gdforce', rhopow=2.0,
          icr0=1.0, icrloop0=0.5, icr0rate=1.5, max_iter=4800,
          margin=8, reorder=True, one_sided=False,
          hmode='knn', nsmoothmin=32, confine=False,
          reanneal=True, momentum=0.0, icr0_stop=3e-3, stop_window=10,
          polish_iter=0, diagnostics=False,
          callback=None, callback_every=100, verbose=False):
    """Relax particles into a glass matching a target density field.

    pos : (N,3) positions (box centered on 0 when periodic); modified copy
        is returned. mass : (N,). density : DensityField (or callable taken
        as rho(pos) over (N,3)). box : (3,) periodic box lengths or None
        for an open domain. reorder : sort the particle
        arrays into cell order at each neighbor rebuild (cache locality;
        results are returned in the input order regardless). one_sided :
        SPH-EXA-style transpose-free force walk - faster, but pair forces
        are no longer exactly antisymmetric (negligible for relaxation;
        keep False when validating against the C generator). hmode : which
        smoothing length the force kernel (and the move step) use -
        'knn' (default, h = 1/2 dist to the nSmooth-th neighbour, as v1),
        'measured' (h = (m*hfact3/rho_sph)^1/3, the density-based h the C
        DenDVDX sets on fBall2 before the force resmooth, hfact3 =
        3*nSmooth/32pi), or 'target' (h = (m*hfact3/rho0)^1/3 from the
        analytic/table target density - noise-free, the cleanest choice for
        a glass generator). Density and pP are always computed with kNN h
        (matching the C, which reuses the DenDVDX density), so only the
        force/move length changes. nsmoothmin : floor on the density-based
        h (measured/target modes only) - h is never allowed below the kNN
        radius enclosing ~nsmoothmin real neighbours, so a particle that
        sits in a sparse region while rho0 says "dense" (analytic h -> ~0)
        still gathers neighbours instead of an empty ball / zero force.
        confine : when True, cap each move component to one smoothing length
        (max_step_frac=1.0) and forbid cross-box jumps - the fine progressive
        phase, where the large-scale field is already in place and only local
        relaxation is wanted (default False = the unbounded coarse/single-relax
        move). reanneal : when False, skip the ICRLOOP0 re-anneal/halve ladder
        in the step controller and drop straight into polish once icr0 decays
        below nfac (also the fine progressive phase; default True = validated
        behaviour). icr0_stop : PRIMARY finish criterion (default 3e-3) - stop
        once the median-smoothed icr0 (the maximum step as a FRACTION of h, the
        move being icvel*icr0*h with icvel normalised so its max component <= 1)
        falls below this threshold. icr0 is resolution-, gradient- and
        particle-count invariant (h scales out), so a fixed threshold travels
        across cases - unlike the residual force frms, whose absolute scale
        depends on mass/h and the density gradient. It IS the quality dial:
        ~3e-3 stops as icr0 settles onto its critical-step setpoint in the
        confined refine stage (denserr ~+7% above asymptote, the fast/"good"
        glass); tighten to ~1e-3 for a near-asymptote glass. The trailing-median
        smooth over stop_window kills the controller's per-iteration icr0
        sawtooth so the crossing is reproducible. This is the SINGLE finish for
        every regime (confined or box-amplified single relax): icr0 is the step
        as a fraction of h, so the threshold is amplification-independent. Set
        icr0_stop=0.0 to disable (then polish_iter or max_iter ends the run).
        polish_iter :
        when > 0, a FIXED polish-iteration budget - stop after exactly this many
        iterations in the polish phase (the analog of the cascade's inter_iter
        fixed-length heal, applied to the full-resolution run: the deep polish
        past entry barely improves quality, so just do a short fixed number).
        Composes with icr0_stop - whichever fires first wins (an early icr0_stop
        convergence still cuts it short); set icr0_stop=0 to make polish_iter the
        sole finish. 0 = off,
        stop_reason='polish_iter'. momentum : heavy-ball coefficient on the move
        step (0.0 = memoryless steepest descent, today's behaviour); > 0 blends
        the previous normalised descent direction in (see stepper.apply_step),
        accelerating the slow frms tail that pure steepest descent cannot beat
        on the stiff 360:1 gradient. Composes with confine (the cap clamps the
        final blended move) and the icr0 feedback controller (which then sees
        momentum-smoothed moves). Remaining parameters mirror the gdicgen
        ICCONFIG flags.
    """
    if hmode not in ('knn', 'measured', 'target'):
        raise ValueError(f"hmode must be knn/measured/target, got {hmode!r}")
    if nsmoothmin >= nsmooth:
        nsmoothmin = nsmooth - 1
    hfact3 = 3.0 * nsmooth / (32.0 * np.pi)
    pos = np.array(pos, dtype=np.float64)
    mass = np.ascontiguousarray(mass, dtype=np.float64)
    if not isinstance(density, DensityField):
        fn = density
        density = type('_Wrap', (DensityField,), {'rho0': lambda self, p: np.asarray(fn(p), dtype=np.float64)})()
    boxarr = None if box is None else np.asarray(box, dtype=np.float64)
    periodic = boxarr is not None
    if periodic:
        # the njit pair loops use compare-based min-image, which requires
        # wrapped coordinates
        pos -= boxarr * np.rint(pos / boxarr)

    n = len(pos)
    # current slot -> original particle id (identity unless reorder swaps)
    order = np.arange(n)

    def to_input_order(arr):
        out = np.empty_like(arr)
        out[order] = arr
        return out
    rho = np.empty(n)
    icvel = np.empty((n, 3))
    q0d = np.empty(n) if diagnostics else _DUMMY[0]
    e0 = np.empty((n, 3)) if diagnostics else _DUMMY
    tail_history = []

    # gather backend: 'ball' for the density-based h modes (no kNN - h sets a
    # mean nSmooth via the h-rho relation and the force/density are a fixed-
    # radius ball gather, as in the CORRPARTITION gasoline path), 'knn' for
    # the kNN candidate lists. 'auto' = ball iff hmode is measured/target.
    # the density-based-h modes (measured/target) use a fixed-radius ball gather
    # (CORRPARTITION gasoline path, consistent density-based h, no kNN); knn h
    # uses the kNN gather. The choice follows hmode directly.
    use_ball = hmode != 'knn'
    if use_ball and not periodic:
        raise ValueError('the measured/target h modes need a periodic box; '
                         'use hmode="knn" for an open domain')

    nblist = NeighborList(boxarr, nsmooth, margin=margin,
                          build_transpose=not one_sided)
    force_loop = get_force_loop(flavour, one_sided=one_sided)
    ballnl = BallNeighborList(boxarr, nsmoothmin=nsmoothmin,
                              build_transpose=not one_sided) if use_ball \
        else None
    ball_force = get_force_loop_csr(flavour, one_sided=one_sided) \
        if use_ball else None
    # In confine mode the per-step move is capped at one h, so the controller's
    # periodic box-amplification (icr0 -> 0.99*box/h_avg, meant to cross the
    # whole box in the early sweep) is exactly what we forbid; passing box=None
    # to the controller disables both that amplification and the step_zero
    # /h_avg, so capped steps start at ~icr0*h and anneal down. The actual move
    # + wrap in apply_step still use the real periodic box.
    ctrl = StepController(icr0, icrloop0, icr0rate, nsmooth,
                          box=(None if confine else boxarr), reanneal=reanneal)
    dummybox = boxarr if periodic else np.ones(3)
    max_step_frac = 1.0 if confine else 0.0

    # finish criterion (all gated to the polish phase):
    #   icr0_stop (default 3e-3) : stop when the median-smoothed icr0 (the max
    #     step as a FRACTION of h - resolution/gradient/amplification invariant)
    #     falls below icr0_stop. The smooth makes the crossing reproducible. This
    #     is the single finish for every regime (confined or box-amplified).
    #   polish_iter (independent) : a fixed polish-iteration budget - stop after
    #     polish_iter iterations in polish regardless of icr0_stop; composes,
    #     first to fire wins.
    # max_iter is the hard cap (set icr0_stop=0 and polish_iter=0 to run it out).
    if icr0_stop > 0.0:
        pstop = Icr0Stop(threshold=icr0_stop, window=stop_window)
        pstop_reason = 'icr0_abs'
    else:
        pstop = None
        pstop_reason = None
    history = []
    converged = False
    stop_reason = 'max_iter'
    prepolish = None   # (pos, h, rho) snapshot at the first polish iteration
    polish_count = 0   # iterations spent in the polish phase (for polish_iter)
    it = 0
    t_start = time.perf_counter()
    if verbose:
        print(f'relax: n={n} nsmooth={nsmooth} flavour={flavour} '
              f'box={None if boxarr is None else boxarr.tolist()}', flush=True)
    pP = np.empty(n)
    pdisp = np.zeros(n)
    # heavy-ball velocity workspace (input-order, permuted alongside pos on
    # reorder so momentum stays particle-aligned; NOT reset on rebuild).
    disp_prev = np.zeros((n, 3))
    nblist.track_displacements(pdisp)
    if use_ball:
        ballnl.track_displacements(pdisp)
        # bootstrap density-based h (measured uses the previous rho; iter 1
        # seeds from the target, matching the C bStarting step)
        rho[:] = density.rho0(pos)
    # per-phase wall-time accumulators (printed in the verbose summary)
    pc = time.perf_counter
    tprof = {'nbr': 0.0, 'dens': 0.0, 'pP': 0.0, 'frc': 0.0, 'step': 0.0,
             'move': 0.0}

    def _hin(pos, mass, rho, rho0):
        return np.cbrt(mass * hfact3 / (rho0 if hmode == 'target' else rho))

    for it in range(1, max_iter + 1):
        if use_ball:
            # ball gather: h = (m*hfact3/rho)^1/3 sets a *mean* nSmooth via
            # the h-rho relation; density, pP and force all use this same h
            # (consistent, no kNN). target h is known a priori from rho0;
            # measured h from the previous iteration's SPH density. Reorder
            # into cell order at each rebuild for cache-local CSR gathers.
            t = pc()
            rho0 = density.rho0(pos) if hmode == 'target' else None
            h_in = _hin(pos, mass, rho, rho0)
            # decide staleness once, on the pre-reorder state; the reorder
            # permutes pos/h_in but not the ball's cached rcached/pdisp, so a
            # second needs_rebuild() inside update() would test misaligned
            # arrays. Pass the decision through instead.
            rebuild = ballnl.needs_rebuild(h_in)
            if reorder and rebuild:
                perm = cell_order(pos, boxarr)
                if perm is not None:
                    newpos = np.empty_like(pos)
                    _gather3(pos, perm, newpos)
                    pos = newpos
                    if momentum > 0.0:
                        newdp = np.empty_like(disp_prev)
                        _gather3(disp_prev, perm, newdp)
                        disp_prev = newdp
                    mass = np.ascontiguousarray(mass[perm])
                    rho = np.ascontiguousarray(rho[perm])
                    order = order[perm]
                    rho0 = density.rho0(pos) if hmode == 'target' else None
                    h_in = _hin(pos, mass, rho, rho0)
            indptr, indices, hforce, rev_indptr, rev_indices = \
                ballnl.update(pos, h_in, rebuild=rebuild)
            tprof['nbr'] += pc() - t
            t = pc()
            density_loop_csr(pos, mass, hforce, indptr, indices, dummybox,
                             periodic, rho)
            tprof['dens'] += pc() - t
            t = pc()
            if hmode != 'target':
                rho0 = density.rho0(pos)
            dens_err = pseudo_pressure(rho, rho0, rhopow, pP)
            tprof['pP'] += pc() - t
            t = pc()
            ball_force(pos, mass, hforce, rho, pP, indptr, indices,
                       rev_indptr, rev_indices, dummybox, periodic, icvel,
                       q0d, e0, diagnostics)
            tprof['frc'] += pc() - t
        else:
            t = pc()
            if reorder and periodic and nblist.needs_rebuild():
                perm = cell_order(pos, boxarr)
                if perm is not None:
                    newpos = np.empty_like(pos)
                    _gather3(pos, perm, newpos)
                    pos = newpos
                    if momentum > 0.0:
                        newdp = np.empty_like(disp_prev)
                        _gather3(disp_prev, perm, newdp)
                        disp_prev = newdp
                    mass = np.ascontiguousarray(mass[perm])
                    order = order[perm]
            idx, h, rev_indptr, rev_indices = nblist.update(pos)
            tprof['nbr'] += pc() - t
            t = pc()
            density_loop(pos, mass, h, idx, dummybox, periodic, rho)
            tprof['dens'] += pc() - t
            t = pc()
            rho0 = density.rho0(pos)
            dens_err = pseudo_pressure(rho, rho0, rhopow, pP)
            # the kNN-gather path is taken only for hmode='knn' (measured/target
            # use the ball path above), so the force kernel uses the same kNN h
            # as the density.
            hforce = h
            tprof['pP'] += pc() - t
            t = pc()
            force_loop(pos, mass, hforce, rho, pP, idx, rev_indptr,
                       rev_indices, dummybox, periodic, icvel, q0d, e0,
                       diagnostics)
            tprof['frc'] += pc() - t
        if diagnostics:
            # per-particle tails (before normalize clobbers icvel): the RAW
            # residual force magnitude and the relative density error. The mean
            # is in history; the p95/max are interface-dominated, so these show
            # whether deep polish is still tightening the TAIL (the 360:1 shell)
            # after the mean has flattened. percentile cost is why it's gated.
            fmag = np.sqrt((icvel * icvel).sum(1))
            err = np.abs(1.0 - rho / rho0)
            tail_history.append((
                it, dens_err, float(np.percentile(err, 95)), float(err.max()),
                float(fmag.mean()), float(np.percentile(fmag, 95)),
                float(fmag.max()), ctrl.polish))
        # snapshot the glass the moment we ENTER polish (this iteration's
        # incoming pos/h/rho, before any polish move), so the caller can
        # measure how much the polish ladder changes the sweep-phase glass.
        if ctrl.polish and prepolish is None:
            prepolish = (to_input_order(pos).copy(),
                         to_input_order(hforce).copy(),
                         to_input_order(rho).copy(), it)
        t = pc()
        frms, _ = normalize_displacements(icvel)
        ctrl.update_step(frms, float(np.mean(hforce)))
        tprof['step'] += pc() - t

        t = pc()
        max_step = apply_step(pos, icvel, hforce, ctrl.icr0, dummybox,
                              periodic, pdisp, max_step_frac, disp_prev,
                              momentum)
        (ballnl if use_ball else nblist).record_move(max_step)
        tprof['move'] += pc() - t
        if not np.isfinite(dens_err):
            raise RuntimeError(
                f"relaxation diverged at iteration {it} (density is no "
                "longer finite). With box=None there is no confinement: "
                "the pseudo-pressure is positive, so a cloud with a free "
                "surface expands without bound - use a periodic box with "
                "a floor density as the C icgenerator does.")
        ctrl.update_restart(dens_err)

        if verbose and (it % 25 == 0 or it == 1):
            print(f"iter {it:5d}  frms {frms:.3e}  denserr {dens_err:.3e}  "
                  f"icr0 {ctrl.icr0:.3e}  icrloop0 {ctrl.icrloop0:.3e}  "
                  f"polish {ctrl.polish}  "
                  f"rebuilds {(ballnl if use_ball else nblist).nbuilds}  "
                  f"t {time.perf_counter() - t_start:.0f}s",
                  flush=True)
        history.append((it, frms, dens_err, ctrl.icr0, ctrl.icrloop0, ctrl.polish))
        if callback is not None and it % callback_every == 0:
            callback(IterStats(it, frms, dens_err, ctrl.icr0, ctrl.icrloop0,
                               ctrl.polish, to_input_order(pos),
                               to_input_order(hforce), to_input_order(rho),
                               to_input_order(rho0)))
        # icr0_stop finish in polish; pstop.update has state, feed only in polish
        if pstop is not None and ctrl.polish:
            if pstop.update(ctrl.icr0):
                converged = True
                stop_reason = pstop_reason
                break
        if polish_iter > 0 and ctrl.polish:
            # fixed polish-iteration budget: do exactly polish_iter polish iters
            # then stop (composes with icr0_stop - first to fire wins)
            polish_count += 1
            if polish_count >= polish_iter:
                converged = True
                stop_reason = 'polish_iter'
                break

    # final state at the relaxed positions, mapped back to input order
    if use_ball:
        if hmode == 'target':
            h_in = np.cbrt(mass * hfact3 / density.rho0(pos))
        else:
            h_in = np.cbrt(mass * hfact3 / rho)
        indptr, indices, hret, _, _ = ballnl.update(pos, h_in)
        density_loop_csr(pos, mass, hret, indptr, indices, dummybox,
                         periodic, rho)
    else:
        idx, h, rev_indptr, rev_indices = nblist.update(pos)
        density_loop(pos, mass, h, idx, dummybox, periodic, rho)
        hret = h   # kNN-gather path is hmode='knn' only
    if verbose:
        tot = sum(tprof.values()) or 1.0
        nb = (ballnl if use_ball else nblist).nbuilds
        print('phase timings (s, %% of summed loop): ' + '  '.join(
            f'{k} {v:.1f} ({100 * v / tot:.0f}%)' for k, v in tprof.items())
            + f'  | rebuilds {nb}/{it} ({100 * nb / it:.0f}%)'
            + f'  | stop {stop_reason} @ {it}', flush=True)
    return RelaxResult(pos=to_input_order(pos), h=to_input_order(hret),
                       rho=to_input_order(rho), niter=it,
                       converged=converged, history=history,
                       tail_history=tail_history,
                       timings=tprof, stop_reason=stop_reason,
                       pos_prepolish=prepolish[0] if prepolish else None,
                       h_prepolish=prepolish[1] if prepolish else None,
                       rho_prepolish=prepolish[2] if prepolish else None,
                       polish_iter0=prepolish[3] if prepolish else -1)
