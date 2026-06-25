#!/usr/bin/env python
"""glassgen icgenerator benchmark - performance + accuracy, end-to-end.

Relaxes a set of standard test cases with the Python glass generator and
reports, per (case, variant): wall time, iterations, per-iteration cost and
particle throughput (performance) plus the SPH density error vs the target
profile and the operator consistency hierarchy - Q0 (partition of unity), Q1
(first-moment error), E0/E1 (0th/1st-order gradient error) - measured
consistently with each variant's smoothing-length operator (accuracy).
Variants exercise the different options - kNN / measured / target
smoothing-length modes and the symmetric vs one-sided force walk.

Every case is driven through the REAL test-suite pipeline rather than a synth
box or a committed pre-IC file:

  1. front end - the `rand` pre-IC is generated on demand with the runtest.sh
     contract `IC_createsetup.py <test> <nx> <vm> 1 <out> [extra]`, into a work
     dir (default <pkg>/_pipeline). If the pre-IC already exists it is reused
     (no regeneration), mirroring runtest.sh's "Skipping creation".
  2. relax  - glassgen relaxes it; box + target density come straight from the
     generated .param's dICdens*/dPeriod fields (port of pkdUpdateDensity
     profiles 1/2/3) - no analytic profiles hardcoded.
  3. back end (--analyze) - the relaxed glass is written as a tipsy snapshot and
     analyzed IN-PROCESS with the framework's render + profile building blocks
     (render_panels density map + binned_profile of SPH rho vs the TARGET rho0
     along the dICdensdir coordinate + a density-error histogram). This is
     glass-appropriate: a glass is a t=0 IC validated against its target density,
     not the per-test evolved analytic solution, so IC_analysis_<test>.py (which
     assumes t>0) is deliberately NOT used.

Run everything (each case to convergence, all variants):
    python benchmark.py
Quick perf probe (fixed iteration budget, no convergence wait):
    python benchmark.py --max-iter 200 --no-stop
Subset, with a resolution override and end-to-end analysis:
    python benchmark.py --cases sedov64 --variants knn --gen-nx 24 --analyze
    python benchmark.py --cases mhdcollapse12 --variants knn,target
Optionally time the C gasoline.gdicgen on the same case (slow):
    python benchmark.py --cases sedov64 --with-c
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys
import time

import numpy as np

_PKG = os.path.dirname(os.path.abspath(__file__))
_TESTSUITE = os.path.dirname(os.path.dirname(_PKG))
_SETUPFILES = os.path.join(_TESTSUITE, 'setupfiles')
_CREATESETUP = os.path.join(_SETUPFILES, 'IC_createsetup.py')
_DEFAULT_WORKDIR = os.path.join(_PKG, '_pipeline')
sys.path.insert(0, _PKG)
sys.path.insert(0, _SETUPFILES)

from glassgen.relax import relax
from glassgen.progressive import relax_progressive
from glassgen.diagnostics import operator_consistency, residual_force
from glassgen.params import (read_param, box_from_param,
                            target_rho0 as target_from_param)

# Each case names a real test, its IC_createsetup resolution (nx = elements
# along x), and any density-structure `extra` args - the positional create()
# parameters AFTER (nx, distri, vm, entry); see setup_<test>().create. Only the
# density structure matters for the glass, so velocity/B-only extras are omitted
# (they keep their create() defaults).
CASES = {
    # fast low-res structured stress case: the mhdcollapse tanh-sphere (360:1)
    'mhdcollapse12': dict(test='mhdcollapse', nx=12, extra=['10', '360.0']),
    # uniform sedov target (profile 1)
    'sedov64':       dict(test='sedov',       nx=64),
    'mhdcollapse50': dict(test='mhdcollapse', nx=50, extra=['10', '360.0']),
    # Sod smoothed step (choice=1, the create() default -> no extra needed)
    'shocktube256':  dict(test='shocktube',   nx=256),
    # uniform current-sheet box (profile 1)
    'currentsheet64':dict(test='currentsheet',nx=64),
    # accdisk (profile 7, disk, dir=5 r_cyl) is the hardest case but its target
    # is not yet ported - see PLAN.md "Future tasks" (needs a density table
    # generated in setupfiles). Enable here once that lands:
    # 'accdisk64': dict(test='accdisk', nx=64, extra=['0.1']),
}

# name -> relax kwargs exercising the options
VARIANTS = {
    'knn':           dict(hmode='knn', one_sided=False),
    'knn-1side':     dict(hmode='knn', one_sided=True),
    'measured':      dict(hmode='measured', one_sided=False),
    'measured-1side':dict(hmode='measured', one_sided=True),
    'target':        dict(hmode='target', one_sided=False),
    'target-1side':  dict(hmode='target', one_sided=True),
}

# --variants help: each name encodes an h-mode + an optional one-sided suffix.
VARIANTS_HELP = (
    'comma list of force-kernel variants (default: all 6). Each name = an '
    'h-mode + an optional "-1side" suffix. h-mode (the smoothing length the '
    'FORCE kernel + move step use; density & pP always use kNN h): '
    'knn = kNN h, 1/2 the distance to the nSmooth-th neighbour (the v1 '
    'default, kNN gather); '
    'measured = density-based h = (m*hfact3/rho_sph)^(1/3) from the SPH '
    'density (the C DenDVDX choice; fixed-radius ball gather); '
    'target = h = (m*hfact3/rho0)^(1/3) from the analytic/table TARGET '
    'density (noise-free; ball gather). '
    '"-1side" = SPH-EXA-style transpose-free one-sided force walk: faster, '
    'but pair forces are not exactly antisymmetric (negligible for '
    'relaxation; drop it to validate against the C generator). '
    'Names: ' + ', '.join(VARIANTS) + '.'
)


# --------------------------------------------------------------------------- #
#  End-to-end pipeline plumbing: IC_createsetup (front) + IC_analysis (back)
# --------------------------------------------------------------------------- #

def _slug(extra):
    """Compact, filesystem-safe tag for the extra create() args (keys distinct
    density structures, e.g. mhdcollapse rhodiff)."""
    return '_'.join(str(e).replace('/', '-').replace(' ', '') for e in extra)


def _ensure_workdir(workdir):
    """The work dir holds datafiles/ + initruns/ exactly as runtest.sh mkdir's
    them; IC_createsetup writes the pre-IC into datafiles/ relative to cwd."""
    df = os.path.join(workdir, 'datafiles')
    os.makedirs(df, exist_ok=True)
    os.makedirs(os.path.join(workdir, 'initruns'), exist_ok=True)
    return df


def gen_pre_ic(test, nx, vm, extra, workdir, force=False):
    """Generate (or reuse) a `rand` pre-IC via the runtest.sh IC_createsetup
    contract. Returns (pre_path, param_path), both absolute.

    If datafiles/<prefile>.00000 already exists and not `force`, it is reused
    (no regeneration) - same test/nx/extra always maps to the same file."""
    nx = int(nx)
    extra = [str(e) for e in (extra or [])]
    df = _ensure_workdir(workdir)
    prefile = f'{test}{nx}_rand'
    if extra:
        prefile += '_' + _slug(extra)
    pre = os.path.join(df, prefile + '.00000')
    param = os.path.join(df, prefile + '.param')
    if os.path.exists(pre) and not force:
        print(f'  reusing pre-IC {os.path.relpath(pre, workdir)}', flush=True)
        return pre, param
    cmd = [sys.executable, _CREATESETUP, test, str(nx), str(vm), '1',
           prefile, *extra]
    env = dict(os.environ, PYTHONPATH=_SETUPFILES)
    print(f'  IC_createsetup: {" ".join(cmd[2:])}', flush=True)
    r = subprocess.run(cmd, cwd=workdir, env=env)
    if r.returncode != 0:
        raise RuntimeError(f'IC_createsetup failed (rc={r.returncode}) for '
                           f'{test} nx={nx} extra={extra}')
    if not os.path.exists(pre):
        raise RuntimeError(f'IC_createsetup produced no {pre}')
    return pre, param


def write_glass_snapshot(pre_path, out_prefix, res, mass):
    """Write the relaxed glass as a tipsy snapshot <out_prefix>.00000, reusing
    the pre-IC's header/dark/star and aux fields (relax preserves particle
    identity + order, so per-particle aux stays 1:1). Returns the snapshot path.

    Mirrors glassgen/cli.py:write_snapshot - positions/density/hsmooth updated,
    velocities zeroed - then copies the pre-IC's sibling aux files so MHD field
    renders remain meaningful."""
    import readtipsy as tip
    gas, dark, star, header, time_ = tip.readtipsy(pre_path)
    snap = out_prefix + '.00000'
    nfin = len(res.pos)
    if nfin == len(gas):
        g = gas.copy()
        g[:, 1:4] = res.pos
        g[:, 4:7] = 0.0
        g[:, 7] = res.rho
        g[:, 9] = res.h
        tip.writetipsy(g.astype(gas.dtype), dark, star, snap, header, time_)
        # relax preserves count -> per-particle aux stays 1:1, byte-copy it
        for src in glob.glob(pre_path + '.*'):
            ext = src[len(pre_path):]         # e.g. '.BField'
            shutil.copyfile(src, snap + ext)
    else:
        # progressive: count changed, build a fresh equal-mass gas array and
        # skip the aux copy (no longer 1:1 with the pre-IC)
        g = np.zeros((nfin, 12), dtype=gas.dtype)
        g[:, 0] = float(np.asarray(mass).sum()) / nfin if np.ndim(mass) else mass
        g[:, 1:4] = res.pos
        g[:, 7] = res.rho
        g[:, 9] = res.h
        hdr = header.copy()
        hdr[0] = nfin
        hdr[2] = nfin
        tip.writetipsy(g, dark[:0], star[:0], snap, hdr, time_)
    return snap


def mapped_coord(pos, param):
    """The coordinate the target density varies along, per the .param's
    dICdensdir (C direction codes 1-6): signed x/y/z, spherical r, cylindrical r,
    or x/y/z again. Matches the `r`/`x` mapping in target_from_param so the
    profile axis is the one the target is structured on."""
    d = int(float(param.get('dICdensdir', 1)))
    if d < 4:
        return pos[:, d - 1], 'xyz'[d - 1]
    if d == 4:
        return np.sqrt((pos ** 2).sum(1)), 'r'
    if d == 5:
        return np.sqrt(pos[:, 0] ** 2 + pos[:, 1] ** 2), 'r_cyl'
    return pos[:, d - 6], 'xyz'[d - 6]


def analyze_glass(snap_prefix, label, box, param, target, acc, save_dir,
                  res=512, backend='grid'):
    """Glass-appropriate analysis using the framework's render + profile building
    blocks directly (no per-test IC_analysis_<test>.py - that assumes an evolved
    t>0 run). Produces, under save_dir:
      * <name>_rho_render.png   - density field map (render_panels, grid SPH
        deposit) of the glass at the mid-plane;
      * <name>_profile.png      - binned SPH rho vs the analytic target rho0 along
        the dICdensdir coordinate;
      * <name>_denserr.png      - histogram of |1 - rho/rho0|, annotated with the
        mean/p95 error and the Q0/Q1/E0/E1 consistency hierarchy.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from IC_analysis_framework import (
        RenderSpec, render_panels, binned_profile, loaddata)

    os.makedirs(save_dir, exist_ok=True)
    name = os.path.basename(snap_prefix)
    save = os.path.join(save_dir, name)

    # 1. density field map - generic rho quantity (tgdata col 7), grid SPH deposit
    spec = RenderSpec(1, 2, lambda fn, tg, ng: tg[:, 7], r'$\rho$',
                      True, 'slice', None, 'inferno', False, None, 'rho', None)
    try:
        render_panels([(label, snap_prefix)], spec, n=1, res=res, save=save,
                      backend=backend, extent=None, project='slice')
    except Exception as e:                       # render is best-effort
        print(f'  (render failed: {e})', flush=True)

    # 2 + 3. profile + density-error histogram from the stored SPH density
    tgdata, *_ = loaddata(snap_prefix + '.00000')
    pos = tgdata[:, 1:4].astype(float)
    if box is not None:
        pos = pos - np.asarray(box, float) * np.rint(pos / np.asarray(box, float))
    rho = tgdata[:, 7].astype(float)
    rho0 = target(pos)
    coord, clabel = mapped_coord(pos, param)

    xc, prof = binned_profile(coord, {'rho': rho, 'rho0': rho0},
                              nbins=min(len(rho) // 50 + 8, 128), stat='median')
    if xc is not None:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(xc, prof['rho'], 'o-', ms=4, label=r'SPH $\rho$ (glass)')
        ax.plot(xc, prof['rho0'], 'k--', lw=2, label=r'target $\rho_0$')
        ax.set_xlabel(clabel); ax.set_ylabel(r'$\rho$')
        ax.set_title(f'{name}: glass density vs target')
        ax.legend()
        fig.tight_layout(); fig.savefig(save + '_profile.png', bbox_inches='tight')
        plt.close(fig)
        print(f'  saved profile -> {os.path.relpath(save + "_profile.png", save_dir)}',
              flush=True)

    err = np.abs(1.0 - rho / rho0)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(err, bins=60, color='C0', alpha=0.8)
    ax.set_xlabel(r'$|1 - \rho/\rho_0|$'); ax.set_ylabel('count')
    ax.set_title(f'{name}: density error  '
                 f'(mean {acc["denserr_mean"]:.4f}, p95 {acc["denserr_p95"]:.4f}, '
                 f'force_mn {acc["force_mean"]:.2e}, Q1 {acc["q1_mean"]:.2e}, '
                 f'E0 {acc["e0_mean"]:.2e}, E1 {acc["e1_mean"]:.2e})')
    fig.tight_layout(); fig.savefig(save + '_denserr.png', bbox_inches='tight')
    plt.close(fig)
    print(f'  saved denserr -> {os.path.relpath(save + "_denserr.png", save_dir)}',
          flush=True)


def _probe_positions(coord, param, n=400):
    """3D positions sampling the dICdensdir coordinate over coord's range, so
    the target rho0 can be evaluated as a smooth analytic curve (the EXACT
    profile) along that axis. Mirrors mapped_coord's direction codes."""
    xs = np.linspace(float(coord.min()), float(coord.max()), n)
    probe = np.zeros((n, 3))
    d = int(float(param.get('dICdensdir', 1)))
    if d < 4:
        probe[:, d - 1] = xs            # signed x/y/z
    elif d == 4 or d == 5:
        probe[:, 0] = xs                # spherical / cylindrical r -> put on x
    else:
        probe[:, d - 6] = xs
    return xs, probe


def plot_run(res, target, param, box, path, title=None, show=True):
    """Combined per-run figure for --plot:
      (left)  the convergence history (RelaxResult.history): frms + denserr on a
              log-y axis, icr0 on a twin axis, polish entry + stop marked;
      (right) the relaxed glass's binned SPH density profile overlaid on the
              EXACT target rho0 (evaluated analytically along the dICdensdir
              coordinate), so the glass density can be compared to the exact one.
    Saves to `path` and, when show=True, also opens an interactive window
    (plt.show()). On a headless node (no DISPLAY) it falls back to Agg and just
    saves."""
    import matplotlib
    import matplotlib.pyplot as plt
    from IC_analysis_framework import binned_profile

    # Pick an INTERACTIVE backend when a display is available so plt.show() opens
    # a window (the user asked for it); fall back to the always-available Agg
    # (save-only) on a headless node. switch_backend actually imports the
    # backend, so a missing GUI toolkit (e.g. no Qt/Tk bindings) is caught here.
    can_show = show
    if show and os.environ.get('DISPLAY'):
        for be in ('TkAgg', 'QtAgg', 'Qt5Agg', 'GTK3Agg'):
            try:
                plt.switch_backend(be)
                break
            except Exception:
                continue
        else:
            plt.switch_backend('Agg')
            can_show = False
    elif show:
        plt.switch_backend('Agg')
        can_show = False

    h = np.asarray(res.history, dtype=float)
    if h.size == 0:
        return None
    it, frms, denserr = h[:, 0], h[:, 1], h[:, 2]
    icr0, polish = h[:, 3], h[:, 5]
    polish_start = it[polish > 0][0] if np.any(polish > 0) else None

    fig, (axH, axP) = plt.subplots(1, 2, figsize=(14, 5.5))

    # --- left: convergence history (frms/denserr + icr0 twin) ---
    l_frms, = axH.semilogy(it, frms, color='C3', lw=1.2,
                           label='frms (residual force)')
    l_de, = axH.semilogy(it, denserr, color='C0', lw=1.2,
                         label=r'denserr  $\langle|1-\rho/\rho_0|\rangle$')
    axH.set_xlabel('iteration'); axH.set_ylabel('error (log)')
    axH.grid(True, which='both', alpha=0.25)
    axHr = axH.twinx()
    l_icr0, = axHr.semilogy(it, icr0, color='C2', lw=1.0, alpha=0.7,
                            label='icr0 (step)')
    axHr.set_ylabel('icr0', color='C2')
    axHr.tick_params(axis='y', labelcolor='C2')
    if polish_start is not None:
        axH.axvline(polish_start, color='0.5', ls='--', lw=1.0)
    axH.axvline(int(it[-1]), color='k', ls=':', lw=1.0)
    # legend only the labelled data lines (not the axvline markers)
    handles = [l_frms, l_de, l_icr0]
    axH.legend(handles, [h.get_label() for h in handles], loc='upper right',
               fontsize=8)
    axH.set_title(f'convergence (stop {res.stop_reason} @ {res.niter})',
                  fontsize=10)

    # --- right: glass density profile vs the EXACT target ---
    pos = np.asarray(res.pos, float)
    if box is not None:
        b = np.asarray(box, float)
        pos = pos - b * np.rint(pos / b)
    rho = np.asarray(res.rho, float)
    coord, clabel = mapped_coord(pos, param)
    xc, prof = binned_profile(coord, {'rho': rho},
                              nbins=min(len(rho) // 50 + 8, 128), stat='median')
    if xc is not None:
        axP.plot(xc, prof['rho'], 'o', ms=4, color='C0',
                 label=r'SPH $\rho$ (glass, binned median)')
    xs, probe = _probe_positions(coord, param)
    axP.plot(xs, target(probe), 'k-', lw=1.8, label=r'target $\rho_0$ (exact)')
    axP.set_xlabel(clabel); axP.set_ylabel(r'$\rho$')
    axP.set_title('glass density vs exact target', fontsize=10)
    axP.legend(fontsize=8); axP.grid(alpha=0.25)

    if title:
        fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=110, bbox_inches='tight')
    if can_show:
        try:
            plt.show()
        except Exception as e:
            print(f'  (plt.show failed, figure saved only: {e})', flush=True)
    plt.close(fig)
    return path


def load_case(name, spec, workdir, gen_nx=None, gen_vm=0.0, force_gen=False):
    """Generate (or reuse) the case's pre-IC and load it. Returns
    (pos, mass, box, target, pre_path, param_dict)."""
    nx = gen_nx if gen_nx is not None else spec['nx']
    pre, param_path = gen_pre_ic(spec['test'], nx, gen_vm, spec.get('extra'),
                                 workdir, force=force_gen)
    import readtipsy as tip
    gas, _, _, _, _ = tip.readtipsy(pre)
    pos = gas[:, 1:4].astype(np.float64)
    mass = gas[:, 0].astype(np.float64)
    param = read_param(param_path)
    return (pos, mass, box_from_param(param), target_from_param(param), pre,
            param)


def accuracy(res, mass, box, target, nsmooth, flavour='gdforce', nsmoothmin=32,
             rhopow=2.0):
    """Density error + residual force + the operator consistency hierarchy on
    the final glass, measured CONSISTENTLY with the variant.

    Rather than re-deriving a fresh kNN density (which the measured/target
    variants never optimized), it reuses the density (res.rho) and force
    smoothing length (res.h) the relaxation finished with - i.e. the same
    h-mode operator that produced the glass - so the reported numbers reflect
    what each variant actually achieved:

      denserr = |1 - rho/rho0|          SPH density vs the target rho0.
      force   = |icvel|                 residual pseudo-pressure (descent)
                force per particle, real pP=(rho/rho0)^rhopow; its RMS is the
                frms the controller drives to force-equilibrium / stop.
      Q0   = sum_j V_j W(h_i)      -> 1 partition of unity (0th-order interp).
      Q1   = first kernel moment   -> 0 1st-order interpolation consistency.
      E0   = 0th-order force error -> 0 residual force at uniform pP.
      E1   = ||M - I||_F           -> 0 1st-order gradient consistency.

    Q0/Q1/E0/E1 come from glassgen.diagnostics.operator_consistency, which
    matches the actual two-sided ICGenerator force operator (same kernel,
    normalization and symmetrization as sph_core)."""
    pos = np.asarray(res.pos, float)
    if box is not None:
        box = np.asarray(box, float)
        pos = pos - box * np.rint(pos / box)
    rho = np.asarray(res.rho, float)
    h = np.asarray(res.h, float)
    mass = np.asarray(mass, float)
    rho0 = target(pos)
    err = np.abs(1.0 - rho / rho0)
    force = residual_force(pos, mass, h, rho, rho0, box, rhopow=rhopow,
                           nsmooth=nsmooth, flavour=flavour,
                           nsmoothmin=nsmoothmin)
    q0, q1, e0, e1 = operator_consistency(
        pos, mass, h, rho, box, nsmooth=nsmooth, flavour=flavour,
        nsmoothmin=nsmoothmin)
    return dict(denserr_mean=err.mean(), denserr_p95=np.percentile(err, 95),
                denserr_max=err.max(),
                force_mean=force.mean(), force_p95=np.percentile(force, 95),
                force_max=force.max(),
                q0_mean=q0.mean(), q0_std=q0.std(),
                q1_mean=q1.mean(), e0_mean=e0.mean(), e1_mean=e1.mean())


def _fmt_acc(acc):
    """One diagnostic block (denserr + force + the operator consistency line)
    as it appears under each variant - shared by the non-polish and polish
    reports."""
    return (f'denserr(mean/p95/max)={acc["denserr_mean"]:.4f}/'
            f'{acc["denserr_p95"]:.4f}/{acc["denserr_max"]:.4f}  '
            f'force(mean/p95/max)={acc["force_mean"]:.2e}/'
            f'{acc["force_p95"]:.2e}/{acc["force_max"]:.2e}\n'
            f'          Q0={acc["q0_mean"]:.4f}+-{acc["q0_std"]:.4f}  '
            f'Q1={acc["q1_mean"]:.2e}  E0={acc["e0_mean"]:.2e}  '
            f'E1={acc["e1_mean"]:.2e}')


def main():
    ap = argparse.ArgumentParser(
        description='glassgen icgenerator performance + accuracy benchmark')
    ap.add_argument('--cases', default=','.join(CASES),
                    help='comma list: ' + ','.join(CASES))
    ap.add_argument('--variants', default=','.join(VARIANTS),
                    help=VARIANTS_HELP)
    ap.add_argument('-s', '--nsmooth', type=int, default=64)
    ap.add_argument('--max-iter', type=int, default=4800)
    ap.add_argument('--margin', type=int, default=8)
    ap.add_argument('--no-stop', action='store_true',
                    help='disable all early-stop criteria (icr0_stop=0, '
                         'polish_iter=0) so the run goes the full --max-iter '
                         'budget.')
    ap.add_argument('--plot', dest='plot', action='store_true',
                    help='show (plt.show) AND save the per-run figures: a '
                         'convergence+profile figure (frms/denserr/icr0 history '
                         'plus the relaxed glass density profile overlaid on the '
                         'EXACT target rho0 along the dICdensdir coordinate) AND '
                         'the in-process analysis (grid SPH density render + '
                         'density-error histogram). All under '
                         '<work-dir>/analysis/.')
    ap.add_argument('--history-every', dest='history_every', type=int,
                    default=0, metavar='N',
                    help='print the per-iteration convergence trajectory '
                         '(it, denserr, frms/force, icr0, polish) sampled '
                         'every N iters for the last (fine) relax stage - lets '
                         'you read where denserr/force actually plateau vs '
                         'where the stop fires (0 = off)')
    ap.add_argument('--work-dir', default=_DEFAULT_WORKDIR,
                    help='where generated pre-ICs / glasses / figures land '
                         '(default <pkg>/_pipeline)')
    ap.add_argument('--gen-nx', type=int, default=None,
                    help='override every case resolution (nx for IC_createsetup) '
                         '- e.g. sweep a case at several resolutions')
    ap.add_argument('--gen-vm', type=float, default=0.0,
                    help='varied-mass arg passed to IC_createsetup (0 = equal '
                         'mass)')
    ap.add_argument('--force-gen', action='store_true',
                    help='regenerate the pre-IC even if it already exists '
                         '(default: reuse the cached file)')
    ap.add_argument('--no-progressive', dest='progressive',
                    action='store_false',
                    help='disable the DEFAULT multi-resolution path and run a '
                         'single full-resolution relax per variant instead. The '
                         'default (progressive) coarsens the input to the auto '
                         'density-chosen coarse count, relaxes it, then cascades '
                         '(x2 splits + short confined heals) back up to full '
                         'res with a final confined relax.')
    ap.add_argument('--coarse-target', dest='coarse_target', default='auto',
                    help='coarse particle count for the default merge path: '
                         '"auto" picks N_coarse_min from the target density '
                         '(max M|grad rho|^3/(alpha^3 rho^4)); an integer forces '
                         'a specific coarse count. The merge factor F is then '
                         'floor_pow2(N_full/target). (default auto)')
    ap.add_argument('--alpha', type=float, default=4.0,
                    help='resolution tolerance for the auto coarse count '
                         '(coarse spacing <= alpha * gradient length; default '
                         '4.0)')
    ap.add_argument('--icr0-stop', dest='icr0_stop', type=float, default=3e-3,
                    help='absolute icr0 stop (the DEFAULT finish): stop once '
                         'the median-smoothed icr0 (max step as a fraction of '
                         'h, resolution-invariant) falls below this. The '
                         'quality dial - 3e-3 = fast/good glass, ~1e-3 = tight. '
                         '0 = off (then max_iter / polish_iter ends the run).')
    ap.add_argument('--polish-iter', dest='polish_iter', type=int, default=0,
                    help='fixed polish-iteration budget on the full-resolution '
                         '(fine) run: stop after this many polish iters (the '
                         'cascade inter-iter heal applied to the final stage). '
                         'Composes with --icr0-stop (earlier convergence still '
                         'wins). 0 = off.')
    ap.add_argument('--momentum', type=float, default=0.5,
                    help='heavy-ball coefficient on the post-split (confined, '
                         'reanneal=False) progressive stages; coarse stays '
                         'beta=0 (default 0.5, accelerates the fine-stage frms '
                         'tail ~2x; 0 = off)')
    ap.add_argument('--inter-iter', dest='inter_iter', type=int, default=10,
                    help='icgen iterations after each intermediate x2 split in '
                         'the progressive cascade (default 10)')
    ap.add_argument('--no-self-bias', dest='self_bias', action='store_false',
                    help='disable the Wendland C2 self-density bias correction '
                         '(gasoline Wzero10) - density then uses the bare W0 '
                         'self term. Default: correction ON (uniform glass '
                         'relaxes to rho ~ 1.0).')
    args = ap.parse_args()

    cases = [c.strip() for c in args.cases.split(',') if c.strip()]
    variants = [v.strip() for v in args.variants.split(',') if v.strip()]
    if args.no_stop:
        # disable every early-stop criterion -> run the full --max-iter budget
        args.icr0_stop = 0.0
        args.polish_iter = 0
    workdir = os.path.abspath(args.work_dir)
    df = _ensure_workdir(workdir)

    rows = []
    for cname in cases:
        if cname not in CASES:
            print(f'unknown case {cname!r}, skipping'); continue
        spec = CASES[cname]
        nx = args.gen_nx if args.gen_nx is not None else spec['nx']
        pos0, mass, box, target, pre, param = load_case(
            cname, spec, workdir, gen_nx=args.gen_nx, gen_vm=args.gen_vm,
            force_gen=args.force_gen)
        n = len(pos0)
        print(f'\n### {cname}: test={spec["test"]} nx={nx} N={n:,}  '
              f'box={None if box is None else box.tolist()}', flush=True)
        for vname in variants:
            if vname not in VARIANTS:
                print(f'  unknown variant {vname!r}, skipping'); continue
            kw = dict(VARIANTS[vname])
            kw['self_bias'] = args.self_bias   # -> relax (all stages)
            t0 = time.perf_counter()
            if args.progressive:
                # multi-resolution: coarsen the full-res input down to the auto
                # density-chosen coarse count (N_coarse_min, alpha), relax it,
                # then cascade (x2 splits + short confined heals) back up to the
                # input count 1:1 with a final confined relax. coarse_target
                # overrides 'auto' with an explicit coarse particle count.
                ct = args.coarse_target
                ct = ct if ct == 'auto' else int(ct)
                res = relax_progressive(
                    pos0, mass, target, box, nx_full=nx, coarse_seed='merge',
                    coarse_target=ct, alpha=args.alpha,
                    inter_iter=args.inter_iter, nsmooth=args.nsmooth,
                    margin=args.margin, max_iter=args.max_iter,
                    momentum=args.momentum,
                    icr0_stop=args.icr0_stop, polish_iter=args.polish_iter,
                    cap=2.0, start_icr0=2.0,
                    verbose=args.history_every > 0, **kw)
            else:
                res = relax(pos0.copy(), mass, target, box=box,
                            nsmooth=args.nsmooth, margin=args.margin,
                            max_iter=args.max_iter,
                            icr0_stop=args.icr0_stop,
                            polish_iter=args.polish_iter,
                            verbose=False, **kw)
            wall = time.perf_counter() - t0
            nfin = len(res.pos)
            mass_fin = np.full(nfin, mass.sum() / nfin)
            acc = accuracy(res, mass_fin, box, target, args.nsmooth)
            spi = wall / res.niter
            thr = nfin * res.niter / wall / 1e6   # Mparticle-iterations / s
            rows.append((cname, vname, nfin, res.niter, res.converged, wall,
                         spi, thr, acc))
            print(f'  {vname:<15} niter={res.niter:<5} wall={wall:7.1f}s  '
                  f'{spi*1e3:6.1f} ms/it  {thr:6.1f} Mp*it/s  N={nfin}',
                  flush=True)
            # ALWAYS report both the pre-polish (end-of-sweep) and the final
            # (post-polish) diagnostic. A polish phase that DEGRADES the glass
            # (icr0 collapse -> freeze) is then visible directly instead of
            # hidden behind the final number.
            if res.pos_prepolish is not None:
                from types import SimpleNamespace
                prepol = SimpleNamespace(pos=res.pos_prepolish,
                                         h=res.h_prepolish, rho=res.rho_prepolish)
                acc_pre = accuracy(prepol, mass_fin, box, target, args.nsmooth)
                print(f'      non-polish @it{res.polish_iter0:<5} '
                      f'{_fmt_acc(acc_pre)}', flush=True)
                print(f'      polish     @it{res.niter:<5} {_fmt_acc(acc)}'
                      f'{"" if res.converged else "  [NOT CONVERGED]"}',
                      flush=True)
            else:
                print(f'      polish-not-reached   {_fmt_acc(acc)}'
                      f'{"" if res.converged else "  [NOT CONVERGED]"}',
                      flush=True)
            for l in getattr(res, 'levels', None) or []:
                print(f'      stage {l["name"]:<6} N={l["n"]:<8} '
                      f'niter={l["niter"]:<5} {l["wall"]:5.0f}s '
                      f'stop={l.get("stop_reason","?"):<13}'
                      f'{"" if l["converged"] else " *nc"}  '
                      f'denserr(mn/p95/mx)={l["denserr_mean"]:.4f}/'
                      f'{l["denserr_p95"]:.4f}/{l["denserr_max"]:.4f}  '
                      f'force(mn/p95/mx)={l["force_mean"]:.2e}/'
                      f'{l["force_p95"]:.2e}/{l["force_max"]:.2e}', flush=True)
            if args.history_every and getattr(res, 'history', None):
                # res.history entries: (it, frms, dens_err, icr0, icrloop0,
                # polish). For the progressive path this is the FINAL (fine)
                # stage trajectory - the dominant stage we care about.
                hist = res.history
                polish_start = next((h[0] for h in hist if h[5]), None)
                print(f'      trajectory (fine stage, every '
                      f'{args.history_every}): polish@'
                      f'{polish_start if polish_start else "-"}', flush=True)
                print('        {:>6} {:>10} {:>10} {:>10} {:>6}'.format(
                    'it', 'denserr', 'force(frms)', 'icr0', 'polish'),
                    flush=True)
                for h in hist:
                    if h[0] % args.history_every == 0 or h is hist[-1]:
                        print('        {:6d} {:10.4e} {:10.4e} {:10.3e} '
                              '{:>6}'.format(h[0], h[2], h[1], h[3],
                                             'yes' if h[5] else 'no'),
                              flush=True)
            if args.plot:
                pdir = os.path.join(workdir, 'analysis')
                os.makedirs(pdir, exist_ok=True)
                # 1) convergence history + glass density profile vs the EXACT
                #    target (shown + saved).
                if getattr(res, 'history', None):
                    ppath = os.path.join(pdir, f'{cname}_{vname}_plot.png')
                    plot_run(res, target, param, box, ppath,
                             title=f'{cname}/{vname}', show=True)
                    print(f'      plot -> {ppath}', flush=True)
                # 2) in-process analysis: grid SPH density render + density-error
                #    histogram (writes the glass snapshot first).
                out_prefix = os.path.join(
                    df, f'{spec["test"]}{nx}_glass_{vname}')
                write_glass_snapshot(pre, out_prefix, res, mass_fin)
                analyze_glass(out_prefix, f'{cname}/{vname}', box, param,
                              target, acc, pdir, res=512, backend='grid')

    # summary table. denserr = SPH density vs target; force = residual descent
    # force (the convergence/stop basis). Both as mean/p95/max.
    print('\n=== summary ===')
    hdr = (f'{"case":<14}{"variant":<15}{"N":>9}{"niter":>6}{"wall_s":>9}'
           f'{"ms/it":>8}{"denserr_mn":>11}{"denserr_p95":>12}'
           f'{"denserr_mx":>11}{"force_mn":>10}{"force_p95":>11}'
           f'{"force_mx":>10}{"Q1":>10}{"E1":>10}')
    print(hdr)
    for (c, v, n, ni, conv, w, spi, thr, acc) in rows:
        print(f'{c:<14}{v:<15}{n:>9}{ni:>6}{w:>9.1f}{spi*1e3:>8.1f}'
              f'{acc["denserr_mean"]:>11.4f}{acc["denserr_p95"]:>12.4f}'
              f'{acc["denserr_max"]:>11.4f}{acc["force_mean"]:>10.2e}'
              f'{acc["force_p95"]:>11.2e}{acc["force_max"]:>10.2e}'
              f'{acc["q1_mean"]:>10.2e}{acc["e1_mean"]:>10.2e}'
              f'{"" if conv else "  *not-converged"}')


if __name__ == '__main__':
    main()
