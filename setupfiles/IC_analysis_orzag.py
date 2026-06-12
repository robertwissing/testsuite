#!/usr/bin/env python3
"""
Analysis for the Orszag-Tang vortex (IC_setup_orzag).

The classic 2-D MHD turbulence test: crossed velocity and magnetic vortices
(v = (-v0 sin 2pi y, v0 sin 2pi x, 0.01), B = (-B0 sin 2pi y, B0 sin 4pi x, 0),
B0 = 1/sqrt(4pi), beta = 10/3) wind up into interacting current sheets. There
is no closed-form solution, so the diagnostics are the field structure and the
energy evolution:

  * E_mag(t)/E_mag(0) vs time -- the magnetic energy first GROWS as the field is
    wound up, then decays as reconnection/dissipation set in (a standard code
    comparison curve),
  * field-map renders (--render) in HORIZONTAL SLICES through the centre (x-y
    plane at z=0) at chosen times (default t = 0.25, 0.5, 1.0, 2.5; override with
    --render-times). 2-D: density rho, the signed field components Bx and By, and
    the divergence error h|div B|/|B|. 3-D: density, |B|, divBerr, and the thermal
    pressure P = (gamma-1) rho u (mid-plane slices of the cube).

The TRUE 3-D Orszag-Tang vortex (Tu et al. 2022, arXiv:2202.03761; Rossazza
2025; original 3-D extension Helzel, Rossmanith & Taetz 2011, arXiv:1007.2606)
is selected with `-a 1` (matching `runtest.sh`'s setup_orzag `-a` flag, which
maps to IC_setup_orzag.create's dim3). It lives on the periodic UNIT cube
[0,1]^3 with the SAME normalisation as the 2-D vortex (rho = 25/36pi, p = 5/12pi,
B0 = 1/sqrt(4pi), k = 2pi), plus a z-perturbation eps = 0.2 of the in-plane
velocities and v_z = eps sin(2pi z). Because the box/velocity scale matches the
2-D test, the evolution time is the same (t = 1 is ~one eddy turnover, as in
Tu's figures) -- NOT Helzel's [0,2pi]^3 scaling, which would be ~2pi slower. The
only analysis change vs 2-D is that the field maps are mid-plane SLICES through
the cube instead of the thin-slab z-average. The energy diagnostic is unchanged.

Standalone (reuses the framework loaders + render_panels, not run_cli).

Usage:
    python IC_analysis_orzag.py <run-dir> [--save out]
    python IC_analysis_orzag.py <run-dir> --render --render-times 0.5 1.0 --save out
    python IC_analysis_orzag.py <run-dir> -a 1 --render --save out   # 3-D vortex
"""

import os
import sys

import numpy as np
import readtipsy as tip

from IC_analysis_framework import (
    Param, build_parser, params_from_args, check_residual_args,
    RenderSpec, spec_slug, render_panels, render_movie,
    aux_magnitude, aux_divB_error, aux_component, thermal_pressure,
    horizontal_cut, mhd_field_cut_quants,
    reference_path_for, save_reference, find_reference, residuals,
    render_reference_path_for, reference_dir_for, resolve_reference_dir,
    copy_run_log, snapshot_times, missing_render_times,
    series_from_log_or_snapshots,
)
from IC_analysis_general import (
    read_vector_aux,
    setup_rcparams, finish_figure_with_legend,
)

GAMMA = 5.0 / 3.0               # IC_setup_orzag.gamma
DEFAULT_RENDER_TIMES = [0.25, 0.5, 1.0, 2.5]      # default panel times (2-D)
DEFAULT_RENDER_TIMES_3D = [0.25, 0.5, 1.0, 2.5]   # default panel times (3-D)

# runtest.sh setup_orzag(): -a selects the true 3-D vortex (dim3, Helzel et al.
# 2011), -b/-c are the viscosities. Only dim3 changes the analysis (slice vs
# slab-average renders + render-time defaults); the viscosities are carried for
# CLI parity.
PARAM_SPEC = [
    Param("a", "dim3", 0.0, float, "0: 2-D slab, 1: true 3-D vortex (Helzel et al. 2011)"),
]


# --------------------------------------------------------------------------- #
#  Field-map renders: density, |B|, divBerr, internal energy
# --------------------------------------------------------------------------- #

# Fixed per-quantity color ranges (same for 2-D and 3-D, every plane), so runs
# and planes are directly comparable instead of each auto-scaling to its own
# percentiles. The Orszag-Tang vortex has well-known scales: rho starts at
# 25/36pi~0.221 and shocks compress it (0.1..0.5), |B| starts at B0=1/sqrt(4pi)
# ~0.282 and amplifies in the current sheets (~0.05..1.0 on a log scale), and the
# thermal pressure P=(gamma-1) rho u starts at p=5/12pi~0.13 and is amplified in
# the shocks/current sheets (~0..0.5 linear).
# Each callback returns None for an integrated `column` projection (different
# units) so that case falls back to percentile auto-scaling.
def _fixed_clim(lo, hi):
    def clim(params, proj):
        return None if proj == "column" else (lo, hi)
    return clim


RHO_CLIM = _fixed_clim(0.1, 0.5)        # linear
BMAG_CLIM = _fixed_clim(0.05, 1.0)      # linear
P_CLIM = _fixed_clim(0.0, 0.5)          # linear (thermal pressure)
# divBerr: fixed LOG range 1e-3 .. 1.0 (3 decades).
DIVBERR_CLIM = _fixed_clim(1e-3, 1.0)


# Slice planes. For the 2-D vortex only the x-y plane (z-slab) is rendered; the
# true 3-D vortex additionally renders the y-z (x=0) and x-z (y=0) mid-planes
# with the SAME quantities. Each entry is (slug, axis1, axis2).
PLANE_XY = ("z0", 1, 2)   # x-y plane (z = mid-plane)
PLANES_3D = [PLANE_XY, ("x0", 2, 3), ("y0", 1, 3)]   # +y-z (x=0), x-z (y=0)


def render_specs(project="slice", axis1=1, axis2=2, is3d=False):
    """The four field maps in the (axis1, axis2) plane.

    2-D vortex: density rho, the signed field components Bx and By, and the
    divergence error |div B| h/|B|, in HORIZONTAL SLICES through the centre
    (project="slice" at z=0). 3-D vortex: rho, |B|, divBerr, thermal pressure P
    (mid-plane slices of the cube). Fixed color ranges: rho/|B|/P LINEAR, Bx/By
    symmetric diverging, divBerr LOG (1e-3..1.0)."""
    # NOTE: rho (bwr) and |B| (rainbow) keep their OT-benchmark cmaps rather than
    # the canonical inferno; each spec carries an explicit slug.
    if is3d:
        return [
            RenderSpec(axis1, axis2, 7, r"$\rho$", False, project, RHO_CLIM,
                       "bwr", False, slug="rho"),
            RenderSpec(axis1, axis2, aux_magnitude("BField"), r"$|B|$", False,
                       project, BMAG_CLIM, "rainbow", False, slug="B"),
            RenderSpec(axis1, axis2, aux_divB_error(),
                       r"$|\nabla\!\cdot\!B|\,h/|B|$", True, project,
                       DIVBERR_CLIM, "inferno", False, slug="divBerr"),
            RenderSpec(axis1, axis2, thermal_pressure(GAMMA), r"$P$", False,
                       project, P_CLIM, "viridis", False, slug="P"),
        ]
    return [
        RenderSpec(axis1, axis2, 7, r"$\rho$", False, project, RHO_CLIM,
                   "bwr", False, slug="rho"),
        RenderSpec(axis1, axis2, aux_component("BField", 0), r"$B_x$", False,
                   project, None, "RdBu_r", True, slug="Bx"),
        RenderSpec(axis1, axis2, aux_component("BField", 1), r"$B_y$", False,
                   project, None, "RdBu_r", True, slug="By"),
        RenderSpec(axis1, axis2, aux_divB_error(), r"$|\nabla\!\cdot\!B|\,h/|B|$",
                   True, project, DIVBERR_CLIM, "inferno", False, slug="divBerr"),
    ]


def do_render(inputs, labels, args, save, is3d=False, bless=False,
              all_times=False, reference=None, movie=False):
    """Render the four field maps at the requested times, forwarding --render-*.

    Each quantity goes to its own figure, suffixed by the field name. The 2-D
    vortex renders only the x-y plane; the true 3-D vortex additionally renders
    the y-z (x=0) and x-z (y=0) mid-plane slices with the same quantities, the
    plane appended to the filename (e.g. `_divBerr_x0`). The 3-D vortex defaults
    to a mid-plane `slice` (cube) and the Tu (2022) figure times; both remain
    overridable via --render-project / --render-times.

    `bless`: write the first input's rendered grids as a render baseline (one npz
    per quantity/plane in the run's references/ folder) instead of drawing.
    `reference`: the --reference value (None -> no overlay; REFERENCE_AUTO ->
    the run's own references/ folder; a folder/file -> that reference set). When
    set and there is a single input, the matching grids are overlaid as
    sim | reference | (sim - reference).
    `movie`: instead of a static multi-time figure, animate each quantity over
    ALL snapshots (one movie per quantity/plane); with `reference` it becomes a
    sim|reference|difference movie when the baseline covers every frame time."""
    entries = list(zip(labels, inputs))
    default_times = DEFAULT_RENDER_TIMES_3D if is3d else DEFAULT_RENDER_TIMES
    # Default to a slice through the CENTRE for both 2-D (thin z-slab -> z=0
    # mid-plane) and 3-D (cube mid-plane); --render-project overrides.
    default_project = "slice"
    times = args.render_times if args.render_times else default_times
    # Bless times: the default/--render-times (--save-reference), or EVERY dump
    # (--save-reference-all). The saved grids are fixed-resolution images,
    # decoupled from N, so blessing all dumps stays cheap and makes the reference
    # complete for any later --render-times. Fall back to `times` if no dumps.
    if bless and all_times:
        bless_times = snapshot_times(inputs[0]) or times
    else:
        bless_times = times
    planes = PLANES_3D if is3d else [PLANE_XY]
    # Reference overlay is opt-in (--reference) and only for a single input.
    ref_dir = (resolve_reference_dir(inputs[0], reference)
               if reference is not None and len(inputs) == 1 else None)

    # --movie: one animation per quantity/plane over ALL snapshots of the first
    # input (movies are single-run); sim|reference|diff when the reference covers
    # every frame time (render_movie falls back to sim-only otherwise).
    if movie:
        for plane_slug, axis1, axis2 in planes:
            for sp in render_specs(default_project, axis1, axis2, is3d=is3d):
                slug = spec_slug(sp)
                if save:
                    mv_save = f"{save}_{slug}_movie"
                    if is3d:
                        mv_save = f"{mv_save}_{plane_slug}"
                else:
                    mv_save = None
                ref_path = None
                if ref_dir is not None:
                    cand = os.path.join(ref_dir,
                                        f"render_{slug}_{plane_slug}.npz")
                    if os.path.isfile(cand):
                        ref_path = cand
                    else:
                        print(f"orzag: no render reference for {slug}/"
                              f"{plane_slug} in {ref_dir}", file=sys.stderr)
                render_movie(inputs[0], labels[0], sp, save=mv_save,
                             res=args.render_res, backend=args.render_backend,
                             project=args.render_project,
                             slice_frac=args.render_slice_frac,
                             extent=args.render_extent, aspect=args.render_aspect,
                             render_reference=ref_path, fps=args.movie_fps,
                             fmt=args.movie_format)
        return
    for plane_slug, axis1, axis2 in planes:
        for sp in render_specs(default_project, axis1, axis2):
            slug = spec_slug(sp)
            if bless:
                # Bless the FIRST input's grids at ALL dump times; one npz per
                # quantity/plane in references/<runname>/ (independent of --save).
                ref_path = render_reference_path_for(inputs[0], slug, plane_slug)
                render_panels(entries[:1], sp, times=bless_times,
                              res=args.render_res,
                              backend=args.render_backend,
                              project=args.render_project,
                              slice_frac=args.render_slice_frac,
                              extent=args.render_extent, aspect=args.render_aspect,
                              save_render_reference=ref_path)
                continue
            if save:
                # 2-D: keep the original `<save>_<field>` names (single plane);
                # 3-D: append the plane (`<save>_<field>_<plane>`).
                sp_save = f"{save}_{slug}"
                if is3d:
                    sp_save = f"{sp_save}_{plane_slug}"
            else:
                sp_save = None
            ref_path = None
            if ref_dir is not None:
                cand = os.path.join(ref_dir, f"render_{slug}_{plane_slug}.npz")
                if os.path.isfile(cand):
                    ref_path = cand
                else:
                    print(f"orzag: no render reference for {slug}/{plane_slug} "
                          f"in {ref_dir}", file=sys.stderr)
            render_panels(entries, sp, times=times, n=args.render_n,
                          res=args.render_res, save=sp_save,
                          backend=args.render_backend, project=args.render_project,
                          slice_frac=args.render_slice_frac,
                          extent=args.render_extent, aspect=args.render_aspect,
                          render_reference=ref_path,
                          residual=getattr(args, "residual", False))


# --------------------------------------------------------------------------- #
#  1-D horizontal cut (rho, Bx, By, div B vs x along y=y0)
# --------------------------------------------------------------------------- #

# The classic Orszag-Tang quantitative comparison is a 1-D cut along a horizontal
# line y=y0 at t=0.5 (Dai & Woodward 1998; Ryu et al. 1998: y=0.3125). Default y0
# is the centre, except y0=0.3125 at t=0.5; override with --y0.
def _y0_for_time(t, override):
    if override is not None:
        return float(override)
    return 0.3125 if abs(t - 0.5) < 1e-3 else 0.0


def do_y0_cut(inputs, labels, args, save):
    """1-D horizontal cuts of rho, Bx, By, div B vs x along y=y0, at each render
    time (default y0 = centre, 0.3125 at t=0.5). Delegates to the framework's
    generic `horizontal_cut`."""
    times = args.render_times if args.render_times else DEFAULT_RENDER_TIMES
    horizontal_cut(inputs, labels, mhd_field_cut_quants(), times,
                   lambda t: _y0_for_time(t, args.y0), save=save,
                   title_prefix="Orszag-Tang ")


# --------------------------------------------------------------------------- #
#  Magnetic energy vs time
# --------------------------------------------------------------------------- #

def emag(fn):
    """(time, E_mag) for snapshot fn, E_mag = sum 1/2 |B|^2 V_i (V_i = m/rho).

    NaN if the BField aux is missing."""
    tgdata, _td, _ts, hdr, time = tip.readtipsy(fn)
    ngas = int(hdr[2])
    t = float(np.asarray(time).ravel()[0])
    B = read_vector_aux(fn, "BField", ngas)
    if B is None:
        return t, np.nan
    vol = tgdata[:, 0] / tgdata[:, 7]
    b2 = np.sum(np.asarray(B, float) ** 2, axis=1)
    return t, float(np.sum(0.5 * b2 * vol))


def emag_series(inp):
    """(t, E_mag, source) for run `inp`.

    Prefers the gasoline '.log' Emag column (sampled every timestep -- far denser
    than the snapshot dumps), and falls back to per-snapshot integration of the
    BField aux when no usable '.log' Emag series is present. Returns (None, None,
    None) if neither source yields data."""
    return series_from_log_or_snapshots(inp, ["dTime", "Emag"], emag)


def plot_emag_vs_time(inputs, labels, save=None, ref_table=None):
    """E_mag(t)/E_mag(0) for each run (grows as the field winds up, then decays).

    Uses the dense '.log' Emag series when available (see emag_series). When
    `ref_table` (a blessed {x, Emag} baseline) is given, it is overlaid as a
    dashed reference curve."""
    import matplotlib.pyplot as plt
    setup_rcparams()
    fig, ax = plt.subplots(figsize=(8, 6))
    any_data = False
    if ref_table is not None and "Emag" in ref_table:
        ax.plot(ref_table["x"], ref_table["Emag"], "--", color="0.35", lw=1.5,
                label="reference", zorder=5)
    for inp, lab in zip(inputs, labels):
        ts, em, source = emag_series(inp)
        if ts is None:
            print(f"orzag: no snapshots/.log for '{inp}'", file=sys.stderr)
            continue
        if not np.isfinite(em[0]) or em[0] == 0:
            print(f"orzag[{lab}]: no/zero initial E_mag (missing BField aux?); "
                  f"skipping.", file=sys.stderr)
            continue
        style = "-" if source == "log" else "o-"   # dense log -> line, sparse -> markers
        ax.plot(ts, em / em[0], style, ms=3, label=lab)
        print(f"orzag[{lab}]: E_mag(0)={em[0]:.4g}  "
              f"E_mag/E_mag0 max={np.max(em) / em[0]:.4g}  "
              f"end={em[-1] / em[0]:.4g}  ({len(ts)} pts from {source})")
        any_data = True
    ax.axhline(1.0, color="0.7", ls="--", lw=1.0)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$E_{mag}(t)/E_{mag}(0)$")
    ax.set_title("Orszag-Tang: magnetic energy vs time")
    fig.tight_layout()
    if not any_data:
        print("orzag: no magnetic energy series to plot.", file=sys.stderr)
    finish_figure_with_legend(fig, ax,
                              save=(f"{save}_emag.png" if save else None))


# --------------------------------------------------------------------------- #
#  Regression baseline (bless / compare the E_mag(t)/E_mag(0) series)
# --------------------------------------------------------------------------- #

def reference_table(inp):
    """Regression table {x: t, Emag: E_mag(t)/E_mag(0)} for run `inp`.

    The normalised magnetic-energy history (wind-up to ~t=0.5 then reconnection
    decay) is the standard Orszag-Tang comparison curve and is resolution-
    independent -> the natural CI baseline. None if no usable E_mag series."""
    ts, em, _src = emag_series(inp)
    if ts is None or not np.isfinite(em[0]) or em[0] == 0:
        return None
    return {"x": ts, "Emag": em / em[0]}


def compare_to_baseline(inp, explicit=None, ref_table=None):
    """Print residuals of the E_mag series vs the baseline; return worst RMS
    (None if no baseline / no series). `ref_table` may be passed in to avoid
    re-discovering the same baseline (and re-printing 'reference found')."""
    if ref_table is None:
        ref_table, _p = find_reference(inp, explicit)
    if ref_table is None:
        return None
    table = reference_table(inp)
    if table is None:
        return None
    worst = 0.0
    for key in ("Emag",):
        res = residuals(table, ref_table, key)
        if res is None:
            continue
        mx, rms = res
        worst = max(worst, rms)
        print(f"orzag[regression] {key} vs baseline: "
              f"max|d|={mx:.4g}  rms={rms:.4g}")
    return worst


def main():
    parser = build_parser(PARAM_SPEC,
                          description="Orszag-Tang vortex analysis.")
    parser.add_argument("--reg-tol", type=float, default=None,
                        help="CI gate: exit non-zero if the worst RMS residual "
                             "vs the regression baseline exceeds this tolerance")
    args = parser.parse_args()
    check_residual_args(parser, args)
    params = params_from_args(args, PARAM_SPEC)

    labels = args.labels if args.labels else [
        os.path.basename(os.path.normpath(i)) for i in args.inputs]
    if len(labels) != len(args.inputs):
        parser.error("--labels must match the number of inputs")

    is3d = int(params["dim3"]) != 0

    # --save-reference / --save-reference-all: bless the FIRST input's E_mag
    # series + render grids + a copy of the run .log. --all blesses render grids
    # at EVERY dump time; the plain flag uses the default/--render-times and
    # refuses if one of those times has no snapshot (would otherwise silently
    # bless the nearest dump at the wrong time).
    if args.save_reference or args.save_reference_all:
        all_times = args.save_reference_all
        if not all_times:
            render_times = (args.render_times if args.render_times else
                            (DEFAULT_RENDER_TIMES_3D if is3d
                             else DEFAULT_RENDER_TIMES))
            miss = missing_render_times(args.inputs[0], render_times)
            if miss:
                print(f"orzag: --save-reference aborted -- requested render "
                      f"time(s) {', '.join(f'{t:g}' for t in miss)} have no "
                      f"snapshot in '{args.inputs[0]}' (run too short / not "
                      f"dumped). Nothing saved; use --save-reference-all to "
                      f"bless every dump, or --render-times within the dumped "
                      f"range.", file=sys.stderr)
                sys.exit(1)
        table = reference_table(args.inputs[0])
        if table is None:
            print("orzag: no E_mag series to save as reference "
                  "(missing BField aux?)", file=sys.stderr)
            sys.exit(1)
        save_reference(table, reference_path_for(args.inputs[0]), params,
                       "orzag")
        copy_run_log(args.inputs[0], reference_dir_for(args.inputs[0]))
        # Bless the render grids too (part of the baseline), regardless of
        # whether --render was passed.
        do_render(args.inputs, labels, args, save=args.save, is3d=is3d,
                  bless=True, all_times=all_times)
        return

    print(f"[orzag] {'3-D (Tu 2022 / Helzel et al. 2011)' if is3d else '2-D'} "
          f"Orszag-Tang vortex")

    # Reference comparison is opt-in: only when --reference is given. Discover the
    # metric baseline ONCE (for the Emag overlay + the regression gate).
    ref_table = None
    if args.reference is not None:
        ref_table, _p = find_reference(args.inputs[0], args.reference)

    plot_emag_vs_time(args.inputs, labels, save=args.save, ref_table=ref_table)
    if args.render:
        do_render(args.inputs, labels, args, save=args.save, is3d=is3d,
                  reference=args.reference)
        if not is3d:                    # 1-D horizontal cut (OT benchmark line)
            do_y0_cut(args.inputs, labels, args, save=args.save)
    if args.movie:
        do_render(args.inputs, labels, args, save=args.save, is3d=is3d,
                  reference=args.reference, movie=True)

    worst = None
    if args.reference is not None:
        for inp in args.inputs:
            w = compare_to_baseline(inp, args.reference, ref_table=ref_table)
            if w is not None:
                worst = w if worst is None else max(worst, w)
    if args.reg_tol is not None and worst is not None and worst > args.reg_tol:
        print(f"orzag: REGRESSION FAILED -- worst rms {worst:.4g} > tol "
              f"{args.reg_tol:.4g}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
