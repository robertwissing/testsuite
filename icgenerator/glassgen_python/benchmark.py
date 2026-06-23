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
from glassgen.progressive import (relax_progressive, default_levels,
                                  cascade_levels)
from glassgen.diagnostics import operator_consistency, residual_force

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


def read_param(path):
    d = {}
    with open(path) as f:
        for line in f:
            line = line.split('#')[0]
            if '=' not in line:
                continue
            k, v = (s.strip() for s in line.split('=', 1))
            d[k] = v
    return d


def box_from_param(p):
    g = lambda k, dv=None: float(p[k]) if k in p else dv  # noqa: E731
    if int(float(p.get('bPeriodic', 1))) == 0:
        return None
    dp = g('dPeriod', 1.0)
    return np.array([g('dxPeriod', dp), g('dyPeriod', dp), g('dzPeriod', dp)])


def target_from_param(p):
    """Faithful port of pkdUpdateDensity profiles 1/2/3 from a .param."""
    prof = int(float(p.get('dICdensprofile', 1)))
    d = int(float(p.get('dICdensdir', 1)))
    R = float(p.get('dICdensR', 0.0))
    Rs = float(p.get('dICdensRsmooth', 1.0))
    inner = float(p.get('dICdensinner', 1.0))
    outer = float(p.get('dICdensouter', 1.0))

    def rho0(pos):
        if d < 4:
            r = np.abs(pos[:, d - 1])
        elif d == 4:
            r = np.sqrt((pos ** 2).sum(1))
        elif d == 5:
            r = np.sqrt(pos[:, 0] ** 2 + pos[:, 1] ** 2)
        else:
            r = pos[:, d - 6]
        if prof >= 4:
            raise NotImplementedError(
                f'density profile {prof} (disk/accdisk) not ported - see '
                'PLAN.md "Future tasks": accdisk density-table setup')
        if prof == 1:
            return np.full(len(pos), outer)
        if prof == 2:
            return np.where(r > R, outer, inner)
        # profile 3
        if d > 3:
            return outer + (inner - outer) / (1.0 + np.exp(-(R - r) / Rs))
        x = pos[:, d - 1]
        return outer + (inner - outer) * 0.5 * (
            np.tanh((x + R) / Rs) - np.tanh((x - R) / Rs))

    return rho0


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


def freeze_probe(res, mass, box, target, nsmooth, save=None, label='',
                 flavour='gdforce', rhopow=2.0):
    """Quantify the spatial concentration of the residual force on the final
    glass - the basis for a per-particle FREEZE / active-set optimisation: most
    of the per-particle residual force |icvel| sits in a thin high-gradient
    shell (the density interface), while the bulk is already at the operator
    noise floor and is just doing noise moves. Reports what fraction of
    particles could be frozen below several multiples of the operator floor and
    how concentrated the total force^2 (the descent 'energy') is, and (if save)
    plots the per-particle force histogram + force vs radius.

    The floor estimate is the uniform-pP residual force per particle (the E0
    operator force in raw |icvel| units): pP == 1, so it is the force a PERFECT
    glass still carries from kernel inconsistency - the level below which a
    particle's force is indistinguishable from SPH discretisation noise."""
    pos = np.asarray(res.pos, float)
    if box is not None:
        box = np.asarray(box, float)
        pos = pos - box * np.rint(pos / box)
    rho = np.asarray(res.rho, float)
    h = np.asarray(res.h, float)
    mass = np.asarray(mass, float)
    rho0 = target(pos)
    # actual residual force (real pP) and the uniform-pP operator floor, both
    # per-particle |icvel| in the same units.
    f = residual_force(pos, mass, h, rho, rho0, box, rhopow=rhopow,
                       nsmooth=nsmooth, flavour=flavour)
    f_floor = residual_force(pos, mass, h, rho, rho0, box, rhopow=rhopow,
                             nsmooth=nsmooth, flavour=flavour, uniform_pP=True)
    floor = np.median(f_floor)               # typical operator-noise force
    n = len(f)
    order = np.argsort(f)[::-1]
    f2 = f * f
    cum = np.cumsum(f2[order]) / f2.sum()
    print(f'  freeze-probe [{label}]: N={n}  frms={np.sqrt(f2.mean()):.3e}  '
          f'operator-floor(med |icvel|,pP=1)={floor:.3e}', flush=True)
    for c in (1.0, 2.0, 3.0, 5.0):
        frac = float(np.mean(f < c * floor))
        print(f'      |icvel| < {c:.0f}x floor : {100*frac:5.1f}% of particles '
              f'freezable', flush=True)
    for topf in (0.01, 0.05, 0.10):
        k = max(1, int(topf * n))
        print(f'      top {100*topf:4.1f}% by force carry '
              f'{100*cum[k-1]:5.1f}% of total force^2', flush=True)
    if save:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        r = np.sqrt((pos * pos).sum(1))
        fig, ax = plt.subplots(1, 2, figsize=(13, 5))
        ax[0].hist(np.log10(np.maximum(f, 1e-30)), bins=80, color='C3',
                   alpha=0.8)
        ax[0].axvline(np.log10(floor), color='k', ls='--', lw=1,
                      label=f'operator floor {floor:.2e}')
        ax[0].set_xlabel('log10 |icvel| (per-particle residual force)')
        ax[0].set_ylabel('count'); ax[0].legend(fontsize=8)
        ax[0].set_title(f'{label}: residual-force distribution')
        ax[1].scatter(r, f, s=2, alpha=0.25, color='C0', rasterized=True)
        ax[1].axhline(floor, color='k', ls='--', lw=1)
        ax[1].set_yscale('log')
        ax[1].set_xlabel('radius from box centre')
        ax[1].set_ylabel('|icvel|')
        ax[1].set_title(f'{label}: force vs radius (interface concentration)')
        fig.tight_layout(); fig.savefig(save, dpi=110, bbox_inches='tight')
        plt.close(fig)
        print(f'      freeze-probe plot -> {save}', flush=True)
    return f, f_floor


def run_c(pre, box, nsmooth, max_iter):
    """Optionally time gasoline.gdicgen on the same pre-IC (profile from its
    own .param). Returns (wall, niter) or None if it can't run."""
    icgen = os.path.join(_TESTSUITE, 'icgenerator', 'gasoline',
                         'gasoline.gdicgen')
    if not os.path.exists(icgen):
        print('  (no gasoline.gdicgen binary - skipping C)')
        return None
    out = os.path.join(_PKG, 'comparison', '_bench_C')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    param = os.path.splitext(pre)[0] + '.param'
    nproc = max(1, (os.cpu_count() or 2) // 2)
    cmd = ['mpirun', '-np', str(nproc), icgen, '-o', out, '-I', pre,
           '-s', str(nsmooth), '-ICR0', '1.0', '-n', str(max_iter),
           '-oi', str(max_iter), '-ICRLOOP0', '0.5', '-ICrhopow', '2.0',
           '-ICR0Rate', '1.5', param]
    for f in glob.glob(out + '*'):
        os.remove(f)
    t0 = time.perf_counter()
    r = subprocess.run(cmd, capture_output=True, text=True)
    wall = time.perf_counter() - t0
    snaps = sorted(glob.glob(out + '.[0-9]*'),
                   key=lambda s: int(s.rsplit('.', 1)[1])
                   if s.rsplit('.', 1)[1].isdigit() else -1)
    niter = int(snaps[-1].rsplit('.', 1)[1]) if snaps else -1
    if r.returncode != 0:
        print(f'  (C run failed rc={r.returncode})')
        return None
    return wall, niter


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
    ap.add_argument('--freeze-probe', dest='freeze_probe', action='store_true',
                    help='after each relax, report how spatially concentrated '
                         'the residual force is (freezable fraction below '
                         'multiples of the operator floor + force^2 '
                         'concentration) and save a force-distribution / '
                         'force-vs-radius plot - sizes the per-particle '
                         'freeze/active-set optimisation.')
    ap.add_argument('--tail-plots', dest='tail_plots', action='store_true',
                    help='record per-iteration p95/max of denserr + residual '
                         'force (relax diagnostics=True) and save '
                         '<case>_<variant>_history_{p95,max}.png - shows whether '
                         'deep polish keeps tightening the interface tail after '
                         'the mean has flattened. Adds percentile cost per iter.')
    ap.add_argument('--history-plot', dest='history_plot', action='store_true',
                    help='save a convergence plot (frms + denserr + icr0 vs '
                         'iteration, polish entry + stop marked) per run to '
                         '<work-dir>/analysis/<case>_<variant>_history.png '
                         '(uses RelaxResult.history; the fine stage for '
                         'progressive runs).')
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
    ap.add_argument('--analyze', action='store_true',
                    help='after relaxing, write the glass snapshot and analyze it '
                         'in-process (density render + rho-vs-target profile + '
                         'density-error histogram) under <work-dir>/analysis')
    ap.add_argument('--analyze-res', type=int, default=512,
                    help='render pixel resolution for --analyze (default 512)')
    ap.add_argument('--analyze-backend', default='grid',
                    help='render backend for --analyze (default grid: sph_interp '
                         '2-D SPH deposit, pynbody-free)')
    ap.add_argument('--progressive', action='store_true',
                    help='multi-resolution: relax a coarse set (nx-coarse), '
                         'split up to full res, then a fine confined relax '
                         '(instead of one single-resolution relax per variant)')
    ap.add_argument('--nx-coarse', dest='nx_coarse', type=int, default=32,
                    help='coarse resolution for the LEGACY --legacy-coarse path '
                         '(default 32)')
    ap.add_argument('--coarse-target', dest='coarse_target', default='auto',
                    help='coarse particle count for the default merge path: '
                         '"auto" picks N_coarse_min from the target density '
                         '(max M|grad rho|^3/(alpha^3 rho^4)); an integer forces '
                         'a specific coarse count. The merge factor F is then '
                         'floor_pow2(N_full/target). (default auto)')
    ap.add_argument('--alpha', type=float, default=2.0,
                    help='resolution tolerance for the auto coarse count '
                         '(coarse spacing <= alpha * gradient length; default '
                         '2.0, calibrated on mhdcollapse 360:1 -> nx_coarse~25)')
    ap.add_argument('--legacy-coarse', dest='legacy_coarse',
                    action='store_true',
                    help='use the old IC_createsetup nx_coarse coarse seed '
                         '(separate low-res file) instead of merging the input')
    ap.add_argument('--cascade', action='store_true',
                    help='progressive: climb to full res one x2 split at a '
                         'time, running a few confined iters after each split '
                         '(multigrid-style), instead of one big jump. Implies '
                         '--progressive.')
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
    ap.add_argument('--momentum', type=float, default=0.0,
                    help='heavy-ball coefficient on the post-split (confined, '
                         'reanneal=False) progressive stages; coarse stays '
                         'beta=0 (default 0.0 = off; ~0.3-0.5 accelerates the '
                         'fine-stage frms tail)')
    ap.add_argument('--inter-iter', dest='inter_iter', type=int, default=10,
                    help='icgen iterations after each intermediate x2 split '
                         'in --cascade (default 10)')
    ap.add_argument('--with-c', action='store_true',
                    help='also time gasoline.gdicgen on each case (slow)')
    args = ap.parse_args()
    # --cascade is a flavour of the progressive multi-resolution path; enable it
    # implicitly so `--cascade` alone works (otherwise the single-relax branch
    # runs and the flag is silently ignored).
    if args.cascade:
        args.progressive = True

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
            if args.tail_plots:
                kw['diagnostics'] = True
            t0 = time.perf_counter()
            if args.progressive and args.legacy_coarse:
                # LEGACY: coarse seed = a fresh IC_createsetup pre-IC at
                # nx_coarse (separate low-res file), à la option A in PLAN.md.
                def coarse_seed(_n, _box, total_mass, _test=spec['test'],
                                _extra=spec.get('extra')):
                    cpre, _ = gen_pre_ic(_test, args.nx_coarse, args.gen_vm,
                                         _extra, workdir, force=args.force_gen)
                    import readtipsy as tip
                    cg, _, _, _, _ = tip.readtipsy(cpre)
                    cpos = cg[:, 1:4].astype(np.float64)
                    return cpos, np.full(len(cpos), total_mass / len(cpos))
                if args.cascade:
                    levels = cascade_levels(
                        nx, args.nx_coarse, inter_iter=args.inter_iter,
                        momentum=args.momentum, icr0_stop=args.icr0_stop,
                        polish_iter=args.polish_iter,
                        coarse_kw=dict(max_iter=args.max_iter),
                        fine_kw=dict(max_iter=args.max_iter))
                else:
                    levels = default_levels(
                        nx, args.nx_coarse, momentum=args.momentum,
                        icr0_stop=args.icr0_stop, polish_iter=args.polish_iter,
                        coarse_kw=dict(max_iter=args.max_iter),
                        fine_kw=dict(confine=True, reanneal=False,
                                     max_iter=args.max_iter))
                res = relax_progressive(
                    pos0, mass, target, box, nx_full=nx,
                    nx_coarse=args.nx_coarse, levels=levels,
                    coarse_seed=coarse_seed, nsmooth=args.nsmooth,
                    margin=args.margin,
                    verbose=False, **kw)
            elif args.progressive:
                # DEFAULT progressive: merge the full-res input down to the
                # density-chosen coarse count (auto N_coarse_min, alpha=2),
                # relax, cascade back to the input count 1:1. coarse_target
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
                pre = SimpleNamespace(pos=res.pos_prepolish,
                                      h=res.h_prepolish, rho=res.rho_prepolish)
                acc_pre = accuracy(pre, mass_fin, box, target, args.nsmooth)
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
            if args.freeze_probe:
                pdir = os.path.join(workdir, 'analysis')
                os.makedirs(pdir, exist_ok=True)
                freeze_probe(res, mass_fin, box, target, args.nsmooth,
                             save=os.path.join(
                                 pdir, f'{cname}_{vname}_freeze.png'),
                             label=f'{cname}/{vname}')
            if args.history_plot and getattr(res, 'history', None):
                from glassgen.diagnostics import plot_history
                pdir = os.path.join(workdir, 'analysis')
                os.makedirs(pdir, exist_ok=True)
                ppath = os.path.join(
                    pdir, f'{cname}_{vname}_history.png')
                plot_history(res.history, ppath,
                             title=f'{cname}/{vname}  (stop {res.stop_reason} '
                                   f'@ {res.niter})',
                             stop_reason=res.stop_reason)
                print(f'      history plot -> {ppath}', flush=True)
            if args.tail_plots and getattr(res, 'tail_history', None):
                from glassgen.diagnostics import plot_error_history
                pdir = os.path.join(workdir, 'analysis')
                os.makedirs(pdir, exist_ok=True)
                for stat in ('p95', 'max'):
                    ppath = os.path.join(
                        pdir, f'{cname}_{vname}_history_{stat}.png')
                    plot_error_history(
                        res.tail_history, ppath, stat=stat,
                        title=f'{cname}/{vname} {stat} errors  '
                              f'(stop {res.stop_reason} @ {res.niter})',
                        stop_reason=res.stop_reason)
                    print(f'      {stat} error plot -> {ppath}', flush=True)
            if args.analyze:
                out_prefix = os.path.join(
                    df, f'{spec["test"]}{nx}_glass_{vname}')
                write_glass_snapshot(pre, out_prefix, res, mass_fin)
                analyze_glass(out_prefix, f'{cname}/{vname}', box, param,
                              target, acc, os.path.join(workdir, 'analysis'),
                              res=args.analyze_res, backend=args.analyze_backend)
        if args.with_c:
            print('  running gasoline.gdicgen (C) ...', flush=True)
            c = run_c(pre, box, args.nsmooth, args.max_iter)
            if c:
                cw, cn = c
                print(f'  {"C-gdicgen":<13} niter={cn:<5} wall={cw:7.1f}s  '
                      f'{cw/max(cn,1)*1e3:6.1f} ms/it  '
                      f'{n*max(cn,1)/cw/1e6:6.1f} Mp*it/s', flush=True)

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
