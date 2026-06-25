"""Multi-resolution (progressive) glass relaxation.

Relax a COARSE particle set first to capture the large-scale density field
cheaply, then `split` particles up to full resolution and run a FINE, LOCAL
relaxation (move capped at one smoothing length, no cross-box jumps, no
re-anneal loop) that only fixes small-scale structure. Because the coarse glass
already encodes the large-scale field, the fine phase needs neither the big
box-crossing sweep nor the ICRLOOP0 re-anneal ladder.

Schedule is a list of per-level configs (see DEFAULT level builder); a MEDIUM
stage is added simply by inserting another level dict. The coarse stage uses a
FIXED ABSOLUTE nx_coarse: if nx_full <= nx_coarse the coarse stage is skipped
entirely (the driver reduces to a single plain relax of the input).
"""
import time

import numpy as np

from .density import DensityField
from .diagnostics import residual_force
from .relax import relax
from . import resample


def _level_diag(pos, mass, h, rho, density, box, nsmooth, flavour, rhopow,
                nsmoothmin):
    """End-of-stage error on a level's glass: denserr |1-rho/rho0| and the
    residual descent force |icvel|, each as mean/p95/max. Lets us see how well
    the coarse stage already fit the field before the fine stage refines it."""
    pos = np.asarray(pos, float)
    if box is not None:
        box = np.asarray(box, float)
        pos = pos - box * np.rint(pos / box)
    rho = np.asarray(rho, float)
    rho0 = density.rho0(pos)
    err = np.abs(1.0 - rho / rho0)
    force = residual_force(pos, mass, h, rho, rho0, box, rhopow=rhopow,
                           nsmooth=nsmooth, flavour=flavour,
                           nsmoothmin=nsmoothmin)
    return dict(denserr_mean=float(err.mean()),
                denserr_p95=float(np.percentile(err, 95)),
                denserr_max=float(err.max()),
                force_mean=float(force.mean()),
                force_p95=float(np.percentile(force, 95)),
                force_max=float(force.max()))


def _measure(pos, mass, density, box, nsmooth, relax_kw):
    """Measure (pos_wrapped, h, rho) on a STATIC particle set via one zero-step
    relax iteration (icr0=0 -> particles do not move), so the post-split state
    is measured with exactly the same hmode/flavour conventions the fine
    stage uses (the split itself only produces positions + mass; rho/h are not
    defined until a neighbour pass)."""
    mkw = {k: v for k, v in relax_kw.items()
           if k not in ('icr0', 'max_iter', 'icr0_stop',
                        'confine', 'reanneal', 'callback', 'callback_every')}
    r = relax(pos, mass, density, box=box, nsmooth=nsmooth, max_iter=1,
              icr0=0.0, icr0_stop=0.0, **mkw)
    return r.pos, r.h, r.rho


def _nearest_pow2(x):
    """Nearest power of two to x (>=1)."""
    x = max(1.0, float(x))
    lo = 1 << int(np.floor(np.log2(x)))
    hi = lo << 1
    return lo if (x - lo) <= (hi - x) else hi


def _floor_pow2(x):
    """Largest power of two <= x (>=1)."""
    return 1 << max(0, int(np.floor(np.log2(max(1.0, float(x))))))


def auto_coarse_count(density, box, alpha=2.0, grid=128, eps_frac=2e-4):
    """Minimum coarse particle count that still resolves the target density
    field everywhere, for an EQUAL-MASS glass:

        N_coarse_min = max_x  M * |grad rho(x)|^3 / (alpha^3 * rho(x)^4),
        M = integral rho dV.

    Derivation: equal mass m=M/N -> local spacing Delta=(m/rho)^(1/3); requiring
    Delta <= alpha * (rho/|grad rho|) (resolve the gradient length with ~1/alpha
    particles) and solving for N gives the bound above. The max picks whichever
    point is hardest to resolve (interface flank, sphere edge, ...) with NO
    geometry assumed - it is purely a function of the density table + alpha, and
    is INDEPENDENT of the input particle count (that only enters the merge
    factor F = floor_pow2(N_full / N_coarse_min)).

    Computed on a grid^3 lattice over the box; |grad rho| by central differences
    with step eps_frac*box, DECOUPLED from the grid spacing so a thin interface
    is not under-resolved (converged from grid=96 on the mhdcollapse 360:1 case).
    Uniform target (grad rho == 0) -> 0 (the caller then merges to its floor).
    alpha=2.0 is calibrated on mhdcollapse 360:1 -> N_coarse_min ~34k (the
    nx_coarse~25-32 sweet spot; nx12 < that -> F=1, no merge)."""
    rho0 = density.rho0 if isinstance(density, DensityField) else density
    box = np.asarray(box, dtype=np.float64)
    # UNIFORM physical spacing (anisotropic per-axis counts), not an isotropic
    # grid^3 - on an elongated box (e.g. shocktube x=2.0) an isotropic grid
    # over-samples the thin axes and leaves the long axis too coarse to resolve
    # a thin feature (the shock interface), badly under-estimating |grad rho|.
    # Same total budget (grid^3 points), spacing = (vol/budget)^(1/3).
    budget = grid ** 3
    spacing = (float(np.prod(box)) / budget) ** (1.0 / 3.0)
    ns = np.maximum(2, np.round(box / spacing).astype(np.int64))
    ax = [(np.arange(ns[d]) + 0.5) / ns[d] * box[d] - 0.5 * box[d]
          for d in range(3)]
    X, Y, Z = np.meshgrid(*ax, indexing='ij')
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    rho = np.asarray(rho0(pts), dtype=np.float64)
    eps = eps_frac * box   # local-derivative step (decoupled from the spacing)
    g2 = np.zeros(len(pts))
    for d in range(3):
        o = np.zeros(3)
        o[d] = eps[d]
        gd = (np.asarray(rho0(pts + o), dtype=np.float64)
              - np.asarray(rho0(pts - o), dtype=np.float64)) / (2 * eps[d])
        g2 += gd * gd
    gmag = np.sqrt(g2)
    M = float(rho.mean()) * float(np.prod(box))
    S = float(np.max(gmag ** 3 / np.maximum(rho, 1e-30) ** 4))
    return M * S / alpha ** 3


def default_levels(nx_full, nx_coarse, coarse_kw=None, fine_kw=None,
                   momentum=0.0, icr0_stop=3e-3, polish_iter=0,
                   cap=None, start_icr0=1.0, cap_mode='h'):
    """Build the default coarse->fine schedule for a target nx_full given a
    fixed nx_coarse. Returns a list of level dicts. When nx_full <= nx_coarse
    the coarse stage is dropped and a single full-resolution level is returned.

    Each level dict: name, split_factor (applied BEFORE the level's relax),
    plus relax() kwargs (confine, reanneal, icr0, icrloop0, max_iter,
    icr0_stop, ...). coarse_kw / fine_kw override the respective defaults.

    momentum : heavy-ball coefficient applied to the POST-SPLIT fine stage only
        (confined, reanneal=False - where it is validated to accelerate the slow
        frms tail, ~2x fewer iters at beta=0.5). The coarse box-wide reanneal=True
        stage is left at beta=0 (its wide icr0 range + re-anneal ladder are an
        untested regime for momentum, and it already moves fast via box
        amplification). 0.0 = off (memoryless steepest descent everywhere).
    icr0_stop : median-smoothed icr0 stop (max step as a fraction of h), the
        finish for the coarse AND fine stages (default 3e-3; it is
        amplification-independent so it applies to the box-amplified coarse
        stage too). The short intermediate heals are limited by max_iter=
        inter_iter instead. 0.0 = off (then max_iter / polish_iter ends the run).
    polish_iter : fixed polish-iteration budget on the POST-SPLIT fine stage
        only (the full-resolution run) - stop after this many polish iters
        regardless of icr0_stop (which still cuts it short if it converges
        first). 0 = off; the coarse/intermediate stages keep relax()'s off
        default.
    cap : per-step move cap on the POST-SPLIT fine stage, in units of h
        (relax's max_step_frac). None = follow confine = clamp each move
        component to 1*h (the default; only the momentum term can then exceed
        1*h, see relax/stepper). 0.0 = no cap. The fine stage keeps confine's
        box=None controller either way, so start_icr0 stays in units of h.
    start_icr0 : starting icr0 for the POST-SPLIT fine stage in units of h
        (default 1.0). Because the fine stage is confined (box=None controller,
        no box amplification) this IS the iter-0 step as a fraction of h; the
        controller then shrinks it (/~2.2 on iter 1) and self-regulates to its
        critical-step setpoint, so this mostly sets the initial transient."""
    # icr0_stop (step < icr0_stop*h) is the finish for EVERY stage - it is
    # amplification-independent, so the box-amplified coarse stage uses it too.
    coarse_defaults = dict(confine=False, reanneal=True, icr0=1.0,
                           icrloop0=0.5, max_iter=4800, icr0_stop=icr0_stop)
    fine_defaults = dict(confine=True, reanneal=False, icr0=start_icr0,
                         icrloop0=0.5, max_iter=4800,
                         momentum=momentum, icr0_stop=icr0_stop,
                         polish_iter=polish_iter, cap_mode=cap_mode)
    if cap is not None:
        fine_defaults['max_step_frac'] = cap
    if coarse_kw:
        coarse_defaults.update(coarse_kw)
    if fine_kw:
        fine_defaults.update(fine_kw)
    if nx_full <= nx_coarse:
        # no coarse stage: one full-res level (split_factor 1), behaving like a
        # plain relax but keeping the move cap / reanneal choice configurable.
        lvl = dict(name='full', split_factor=1,
                   icr0_stop=icr0_stop, polish_iter=polish_iter)
        lvl.update(coarse_defaults)
        return [lvl]
    factor = _nearest_pow2((nx_full / nx_coarse) ** 3)
    coarse = dict(name='coarse', split_factor=1)
    coarse.update(coarse_defaults)
    fine = dict(name='fine', split_factor=factor)
    fine.update(fine_defaults)
    return [coarse, fine]


def cascade_levels(nx_full=None, nx_coarse=None, *, factor=None, inter_iter=10,
                   coarse_kw=None, fine_kw=None,
                   inter_kw=None, momentum=0.0, icr0_stop=3e-3, polish_iter=0,
                   cap=None, start_icr0=1.0, cap_mode='h'):
    """Coarse -> CASCADE -> fine schedule (a multigrid-style prolongation).

    Instead of the single big x(nx_full/nx_coarse)^3 jump that default_levels
    does, relax the coarse set, then repeatedly binary-split (x2) and run a few
    confined icgen iterations, climbing one resolution level at a time up to
    nx_full, then finish with a full confined relax. One x2 split at a time
    keeps each prolongation a SMALL local perturbation (children land ~0.4 dr
    off the parent) that the handful of intermediate iters heals, so the final
    relax starts from a far healthier glass than after one x8/x64 jump (whose
    post-split denserr is ~0.2 and residual force ~50x the converged value).

    inter_iter : icgen iterations run after each intermediate x2 split (the
        final split is followed by a full relax, not this short smooth). The
        intermediates are fixed short heals: icr0_stop is off, so max_iter=
        inter_iter limits them.

    When nx_full <= nx_coarse a single full-resolution level is returned (as in
    default_levels). coarse_kw/fine_kw/inter_kw override the respective defaults.

    momentum : heavy-ball coefficient applied to the FINAL fine relax only. The
        short intermediate heals (inter_iter iters) and the coarse box-wide
        stage stay beta=0 - momentum needs more iterations than the heals run
        to settle, and on a high-force post-split state it overshoots. 0.0=off.
    """
    # icr0_stop (step < icr0_stop*h, amplification-independent) finishes the
    # coarse and fine stages. The intermediates are fixed short heals limited by
    # max_iter=inter_iter (icr0_stop off - they won't converge in inter_iter).
    coarse_defaults = dict(confine=False, reanneal=True, icr0=1.0,
                           icrloop0=0.5, max_iter=4800, icr0_stop=icr0_stop)
    fine_defaults = dict(confine=True, reanneal=False, icr0=start_icr0,
                         icrloop0=0.5, max_iter=4800,
                         momentum=momentum, icr0_stop=icr0_stop,
                         polish_iter=polish_iter, cap_mode=cap_mode)
    if cap is not None:
        fine_defaults['max_step_frac'] = cap
    # NB momentum is NOT applied to the intermediate heals: they run only
    # inter_iter (~10) iterations from a high-force post-split state, too few
    # for the heavy-ball term to settle, so it overshoots (leaving the next
    # split a worse start). Momentum is reserved for the final fine relax,
    # which has enough iterations to benefit. (Steepest descent here.)
    inter_defaults = dict(confine=True, reanneal=False, icr0=1.0,
                          icrloop0=0.5, max_iter=inter_iter,
                          momentum=0.0, icr0_stop=0.0)
    if coarse_kw:
        coarse_defaults.update(coarse_kw)
    if fine_kw:
        fine_defaults.update(fine_kw)
    if inter_kw:
        inter_defaults.update(inter_kw)
    # total resolution gain as k binary (x2) splits - either an explicit
    # power-of-two `factor` (the merge path, count-based) or the nx ratio.
    if factor is not None:
        total = _floor_pow2(factor)
    else:
        total = _nearest_pow2((nx_full / nx_coarse) ** 3)
    if total <= 1:
        # single full-res level (factor==1, or nx_full <= nx_coarse): honour the
        # caller's icr0_stop / polish_iter (a passed icr0_stop=0 must disable the
        # finish, not inherit relax()'s 3e-3 default) - see default_levels.
        lvl = dict(name='full', split_factor=1,
                   icr0_stop=icr0_stop, polish_iter=polish_iter)
        lvl.update(coarse_defaults)
        return [lvl]
    k = int(round(np.log2(total)))
    coarse = dict(name='coarse', split_factor=1)
    coarse.update(coarse_defaults)
    levels = [coarse]
    # k-1 intermediate (x2 split + short smooth) stages ...
    for i in range(1, k):
        lv = dict(name=f'cascade{i}', split_factor=2)
        lv.update(inter_defaults)
        levels.append(lv)
    # ... then the final x2 split brings us to nx_full and gets the full relax
    fine = dict(name='fine', split_factor=2)
    fine.update(fine_defaults)
    levels.append(fine)
    return levels


def relax_progressive(pos, mass, density, box, *, nx_full=None, nx_coarse=64,
                      levels=None, nsmooth=64, coarse_seed='merge',
                      coarse_target='auto', alpha=2.0, inter_iter=10,
                      coarse_floor=13824, max_iter=4800,
                      split_margin=8, split_dist=0.4,
                      momentum=0.0, icr0_stop=3e-3, polish_iter=0,
                      cap=None, start_icr0=1.0, cap_mode='h',
                      verbose=False, **relax_kw):
    """Progressive coarse->split->fine glass relaxation.

    pos, mass : the FULL-resolution pre-IC (N = nx_full**3). In the default
        two-level path these positions are NOT relaxed directly - the coarse
        stage starts from `coarse_seed` (or a uniform-random fallback at
        nx_coarse) and the fine stage is produced by splitting the coarse glass;
        the input is used only for the total mass and as the seed when the
        coarse stage is skipped (nx_full <= nx_coarse).
    density : DensityField or callable rho0(pos). box : (3,) periodic lengths.
    nx_full : target resolution (elements/side). nx_coarse : fixed absolute
        coarse resolution. levels : explicit schedule (overrides the default).
    coarse_seed : (M,3) array OR callable(n_target, box, total_mass) ->
        (pos, mass) supplying the coarse seed; None -> uniform random. The
        string 'merge' instead COARSENS the full-res input down by the cascade's
        total split factor (resample.coarsen) - so splitting it back restores the
        input count EXACTLY (1:1), with the coarse seed derived from the high-res
        pre-IC itself (no separate low-res file).
    relax_kw : forwarded to every relax() call (e.g. flavour, hmode,
        one_sided, diagnostics). Per-level keys (confine, reanneal, icr0,
        icrloop0, max_iter, icr0_stop, momentum) come from the level dict.
    momentum : heavy-ball coefficient injected into the POST-SPLIT levels only
        (via default_levels) when `levels` is not given explicitly; the coarse
        box-wide stage stays beta=0. When you pass `levels` yourself, bake
        momentum into the level dicts via default_levels/cascade_levels(
        momentum=...) instead - this arg is then ignored.

    Returns the final full-resolution RelaxResult."""
    box = np.asarray(box, dtype=np.float64)
    pos = np.asarray(pos, dtype=np.float64)
    mass = np.asarray(mass, dtype=np.float64)
    total_mass = float(mass.sum())
    if not isinstance(density, DensityField):
        fn = density
        density = type('_Wrap', (DensityField,),
                       {'rho0': lambda self, p: np.asarray(fn(p),
                                                           dtype=np.float64)})()
    if levels is None:
        if isinstance(coarse_seed, str) and coarse_seed == 'merge':
            # DEFAULT: pick the merge/cascade factor from the target density.
            # N_coarse_min = max_x M|grad rho|^3/(alpha^3 rho^4) (resolves the
            # steepest feature; density-table only); F = floor_pow2(N_full /
            # N_coarse_min), capped so the coarse set keeps >= coarse_floor
            # particles. F==1 => no merge => a single full-res relax (e.g. a
            # low-res run, or a near-discontinuous interface). An explicit
            # coarse_target (a particle count) overrides the 'auto' criterion.
            n_full = len(pos)
            if coarse_target == 'auto':
                # coarse count = max(density criterion, floor): the floor
                # (~24^3) is the minimum coarse resolution, so a weak/zero
                # gradient (incl. uniform) coarsens only to the floor, never
                # below; a stronger gradient bumps it up. An input already
                # below the floor (low-res) gives F=1 (no merge).
                target = max(auto_coarse_count(density, box, alpha=alpha),
                             coarse_floor)
            else:
                target = float(coarse_target)   # explicit count bypasses floor
            factor = max(1, _floor_pow2(n_full / max(target, 1.0)))
            levels = cascade_levels(
                factor=factor, inter_iter=inter_iter, momentum=momentum,
                icr0_stop=icr0_stop, polish_iter=polish_iter,
                cap=cap, start_icr0=start_icr0, cap_mode=cap_mode,
                coarse_kw=dict(max_iter=max_iter),
                fine_kw=dict(max_iter=max_iter))
            if verbose:
                print(f"[progressive] auto merge: coarse target={target:.0f} "
                      f"-> factor={factor} (coarse N~{n_full // factor})",
                      flush=True)
        else:
            levels = default_levels(nx_full, nx_coarse, momentum=momentum,
                                    icr0_stop=icr0_stop, polish_iter=polish_iter,
                                    cap=cap, start_icr0=start_icr0,
                                    cap_mode=cap_mode)

    skip_coarse = (len(levels) == 1 and levels[0].get('split_factor', 1) == 1)
    if skip_coarse:
        # full-res single relax of the input directly
        cur_pos, cur_mass, aux = pos.copy(), mass.copy(), {}
    elif isinstance(coarse_seed, str) and coarse_seed == 'merge':
        # coarsen the FULL-RES input down by the cascade's total split factor,
        # so splitting it back up restores the input count exactly (1:1) - no
        # separate low-res file, coarse seed derived from the high-res pre-IC.
        total_factor = 1
        for lv in levels:
            total_factor *= int(lv.get('split_factor', 1))
        cur_pos, cur_mass, aux = resample.coarsen(
            pos, mass, box=box, factor=total_factor, margin=split_margin)
        # coarsen's mutual-nearest merge absorbs a VARYING number of originals
        # per particle (mass std/mean ~0.5), so the merged masses are unequal.
        # The glass must be EQUAL-mass (the coarse positions, denser where the
        # input was denser, already encode the target density), so reset to a
        # uniform coarse mass; the splits then keep it equal down to full res.
        cur_mass = np.full(len(cur_pos), total_mass / len(cur_pos))
    else:
        n_coarse = nx_coarse ** 3
        if coarse_seed is None:
            rng = np.random.default_rng(0)
            cur_pos = (rng.random((n_coarse, 3)) - 0.5) * box
            cur_mass = np.full(n_coarse, total_mass / n_coarse)
        elif callable(coarse_seed):
            cur_pos, cur_mass = coarse_seed(n_coarse, box, total_mass)
            cur_pos = np.asarray(cur_pos, dtype=np.float64)
            cur_mass = np.asarray(cur_mass, dtype=np.float64)
        else:
            cur_pos = np.asarray(coarse_seed, dtype=np.float64)
            cur_mass = np.full(len(cur_pos), total_mass / len(cur_pos))
        aux = {}

    result = None
    level_stats = []
    diag_flavour = relax_kw.get('flavour', 'gdforce')
    diag_rhopow = relax_kw.get('rhopow', 2.0)
    diag_nsmoothmin = relax_kw.get('nsmoothmin', 32)
    for lvl in levels:
        sf = int(lvl.get('split_factor', 1))
        if sf > 1:
            cur_pos, cur_mass, aux = resample.split(
                cur_pos, cur_mass, aux, box=box, factor=sf, nsmooth=nsmooth,
                margin=split_margin, dist=split_dist)
            # post-split diagnostic: the quality the fine stage INHERITS from
            # the split (before any fine relaxation), measured with the same
            # conventions the fine stage will use. Shows how much the split
            # degrades the coarse glass and the starting point the fine stage
            # refines from.
            ps_pos, ps_h, ps_rho = _measure(cur_pos, cur_mass, density, box,
                                            nsmooth, relax_kw)
            ps = dict(name=f"{lvl.get('name')}.split", n=len(cur_pos),
                      niter=0, converged=True, wall=0.0,
                      stop_reason='post-split')
            ps.update(_level_diag(ps_pos, cur_mass, ps_h, ps_rho, density, box,
                                  nsmooth, diag_flavour, diag_rhopow,
                                  diag_nsmoothmin))
            level_stats.append(ps)
            if verbose:
                print(f"[progressive] post-split (N={len(cur_pos)}): "
                      f"denserr(mean/p95/max)={ps['denserr_mean']:.4f}/"
                      f"{ps['denserr_p95']:.4f}/{ps['denserr_max']:.4f}  "
                      f"force(mean/p95/max)={ps['force_mean']:.2e}/"
                      f"{ps['force_p95']:.2e}/{ps['force_max']:.2e}", flush=True)
        lvl_kw = {k: v for k, v in lvl.items()
                  if k not in ('name', 'split_factor')}
        kw = dict(relax_kw)
        kw.update(lvl_kw)
        if verbose:
            print(f"[progressive] level '{lvl.get('name')}' : "
                  f"N={len(cur_pos)} split_factor={sf} "
                  f"confine={kw.get('confine')} reanneal={kw.get('reanneal')}",
                  flush=True)
        t_lvl = time.perf_counter()
        result = relax(cur_pos, cur_mass, density, box=box, nsmooth=nsmooth,
                       verbose=verbose, **kw)
        wall = time.perf_counter() - t_lvl
        stats = dict(name=lvl.get('name'), n=len(cur_pos), niter=result.niter,
                     converged=result.converged, wall=wall,
                     stop_reason=result.stop_reason)
        # end-of-stage error (denserr + residual force) so we can see how well
        # the coarse stage fit the field before the fine stage refined it
        stats.update(_level_diag(result.pos, cur_mass, result.h, result.rho,
                                 density, box, nsmooth, diag_flavour,
                                 diag_rhopow, diag_nsmoothmin))
        level_stats.append(stats)
        if verbose:
            print(f"[progressive] level '{lvl.get('name')}' done: "
                  f"niter={result.niter} converged={result.converged} "
                  f"wall={wall:.1f}s  denserr(mean/p95/max)="
                  f"{stats['denserr_mean']:.4f}/{stats['denserr_p95']:.4f}/"
                  f"{stats['denserr_max']:.4f}  force(mean/p95/max)="
                  f"{stats['force_mean']:.2e}/{stats['force_p95']:.2e}/"
                  f"{stats['force_max']:.2e}", flush=True)
        cur_pos = result.pos
    result.levels = level_stats
    # mass-consistency sanity (legacy nx-based path only; the merge path is
    # exact by construction - split(coarsen(x,F),F) restores the input count).
    if not skip_coarse and nx_full is not None and not (
            isinstance(coarse_seed, str) and coarse_seed == 'merge'):
        expect = total_mass / (nx_full ** 3)
        got = total_mass / len(cur_pos)
        if not np.isclose(got, expect, rtol=0.5):
            import warnings
            warnings.warn(
                f"progressive: final per-particle mass {got:.3e} differs from "
                f"the nx_full target {expect:.3e} (final N={len(cur_pos)} vs "
                f"nx_full**3={nx_full**3}); check nx_coarse / split factor.")
    return result
