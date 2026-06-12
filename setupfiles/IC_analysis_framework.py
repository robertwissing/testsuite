#!/usr/bin/env python3
"""
Generalized analysis driver for the IC test suite.

Every test's analysis script (IC_analysis_<test>.py) is thin: it declares
  * a PARAM spec -- the same -a..-h flags, with the same defaults, as that
    test's setup_<test>() function in runtest.sh,
  * an `analyze(tgdata, time, params, snap)` callback returning a flat dict of
    per-snapshot scalars (one of which is the x-axis key "x"),
  * an `analytic(metric, params)` callback returning (x, y) for the analytic
    overlay (or None where no analytic solution exists),
  * METRICS metadata describing how to plot each scalar,
and then calls run_cli(...) here.

The driver handles, identically for every test:
  * argument parsing (-a..-h plus --save / --save-reference / --reference / ...),
  * looping snapshots and collecting the metric time-series,
  * ALWAYS overlaying the analytic solution,
  * the *reference* (regression-baseline) workflow: --save-reference writes the
    computed series next to the run as "<runname>_reference.json"; a normal run
    auto-discovers a sibling "*_reference.json", overlays it, and prints the
    residuals of the input vs. the reference.

No pass/fail threshold / exit code yet -- residuals are reported, not gated.
That gate is added when CI is wired up.
"""

import argparse
import os
import re
import sys
from collections import namedtuple

import numpy as np

from IC_analysis_general import (
    loaddata, find_files, setup_rcparams, finish_figure,
    try_aux, read_vector_aux, read_bfield, read_log_param, read_tufac,
    read_log_series, code_units, CodeUnits,
    _find_logs,
)
# Reference (regression-baseline) I/O lives in IC_analysis_reference; re-export it
# here so callers keep doing `from IC_analysis_framework import find_reference, ...`.
from IC_analysis_reference import (
    REFERENCE_AUTO, reference_dir_for, reference_path_for, resolve_reference_dir,
    copy_run_log, save_reference, load_reference, find_reference,
    _monotone_xy, residuals,
)
# Field-map rendering lives in IC_analysis_render; re-export it here so
# callers keep doing `from IC_analysis_framework import RenderSpec, ...`.
from IC_analysis_render import (
    _cached_readtipsy, RenderSpec, ContourSpec,
    RenderStyle, RENDER_STYLE, styled, StreamSpec,
    aux_scalar, aux_component, aux_magnitude, aux_divB_error,
    internal_energy, thermal_pressure, entropy_function, magnetic_pressure,
    ke_density, speed, _read_time, read_box_periods,
    _have_pynbody, _field_image, _label_slug, spec_slug,
    _quantity_values, _field_values, _pynbody_image, _smoothing_length,
    _grid_image, _crop_to_view, select_snapshots, missing_render_times,
    snapshot_times, _resolve_render_mode, _render_one_input, _tile_axes,
    save_render_grids, load_render_grids, render_panels, vector_grid,
    streamline_panels, _color_scale_from_pool, _resolve_movie_writer, render_movie,
)


# A single -a..-h parameter: which letter, its long name, default value, the
# type to parse it as, and a help string. `default` must match runtest.sh.
Param = namedtuple("Param", ["letter", "name", "default", "type", "help"])

# A single plotted quantity returned by a test's analyze() callback.
#   norm: optional metric key whose FIRST-snapshot value normalizes this metric
#         (e.g. blob plots M_cloud/M_cloud,0 with norm="mass_cloud"); None = raw.
Metric = namedtuple("Metric",
                    ["key", "ylabel", "xlabel", "yscale", "marker", "norm"])
Metric.__new__.__defaults__ = ("value", "x", "linear", "o", None)  # ylabel..onward


# --------------------------------------------------------------------------- #
#  Argument parsing
# --------------------------------------------------------------------------- #

def build_parser(param_spec, description=""):
    """Build the argparse parser shared by every analysis script.

    param_spec : list[Param]  -- the -a..-h flags for this test.
    """
    # add_help is disabled because runtest.sh uses -h as the 8th test parameter;
    # we expose a long-only --help instead so -h stays free for the param spec.
    p = argparse.ArgumentParser(description=description, add_help=False)
    p.add_argument("--help", action="help", default=argparse.SUPPRESS,
                   help="show this help message and exit")
    p.add_argument(
        "inputs", nargs="+",
        help="one or more run directories / file prefixes / globs to analyze "
             "(files ending with a 5-digit snapshot suffix)",
    )
    for spec in param_spec:
        p.add_argument(
            f"-{spec.letter}", f"--{spec.name}",
            dest=spec.name, type=spec.type, default=spec.default,
            help=f"{spec.help} (runtest -{spec.letter}; default {spec.default})",
        )
    p.add_argument("--labels", nargs="+",
                   help="optional labels for each input (same order)")
    p.add_argument("--save",
                   help="save figure(s) to this filename prefix instead of showing")
    p.add_argument("--save-reference", action="store_true",
                   help="bless the FIRST input as the reference baseline into "
                        "'<parent>/references/<runname>/' (metric reference.json, "
                        "render grids at the default/--render-times, and a copy "
                        "of the run .log)")
    p.add_argument("--save-reference-all", action="store_true",
                   help="like --save-reference but bless render grids at EVERY "
                        "snapshot time (not just the default render times); the "
                        "grids are fixed-resolution images so this stays cheap "
                        "regardless of particle count")
    p.add_argument("--reference", nargs="?", const=REFERENCE_AUTO, default=None,
                   help="compare against a reference (opt-in: with no value, "
                        "auto-discover the run's 'references/<runname>/' folder; "
                        "give a folder or .json to use a specific reference set). "
                        "If omitted entirely, NO reference comparison is done.")
    p.add_argument("--residual", action="store_true",
                   help="draw residual render grids: with ONE input and "
                        "--reference, a times x [sim | reference | difference] "
                        "figure from the blessed render grids; with TWO inputs, "
                        "a times x [run1 | run2 | run1-run2] figure differencing "
                        "the runs directly (no reference needed). Errors with "
                        "more than two inputs. Implies --render.")
    p.add_argument("--no-analytic", action="store_true",
                   help="skip the analytic overlay (debugging)")
    p.add_argument("--movie", action="store_true",
                   help="animate the renders/profiles over ALL snapshots as a "
                        "movie (mp4 via ffmpeg, else gif). Combine with "
                        "--reference for a sim|reference|difference movie when the "
                        "reference covers every frame time (bless with "
                        "--save-reference-all).")
    p.add_argument("--movie-fps", type=float, default=10.0,
                   help="movie frames per second (default 10)")
    p.add_argument("--movie-format", choices=["auto", "mp4", "gif"],
                   default="auto",
                   help="movie container: auto (mp4 if ffmpeg, else gif), mp4, "
                        "or gif")
    p.add_argument("--y0", type=float, default=None,
                   help="height of the 1-D horizontal cut (value vs x along y=y0); "
                        "default is the test's own (usually the centre)")
    p.add_argument("--render", action="store_true",
                   help="also render 2-D field maps at selected snapshots")
    p.add_argument("--render-times", nargs="+", type=float, default=None,
                   help="times to render (nearest snapshot is chosen), e.g. "
                        "--render-times 2 4; default: --render-n evenly spaced")
    p.add_argument("--render-n", type=int, default=4,
                   help="number of evenly-spaced snapshots to render when "
                        "--render-times is not given (default 4)")
    p.add_argument("--render-res", type=int, default=512,
                   help="render pixel resolution (default 512)")
    p.add_argument("--render-backend",
                   choices=["auto", "pynbody", "grid", "grid-petkova",
                            "hist", "particles"],
                   default="grid",
                   help="render backend: grid (default; sph_interp 2-D SPH "
                        "kernel deposit -- pynbody-free SPH-smoothed maps), "
                        "grid-petkova (sph_interp 2-D EXACT Petkova column -- "
                        "mass-conserving per pixel even when h<pixel; slower; "
                        "column/average/rhocolumn only, slice falls back to SPH), "
                        "auto (pynbody SPH if available, else histogram), pynbody, "
                        "hist (gridded), or particles (per-particle scatter "
                        "colored by the field)")
    p.add_argument("--render-extent", nargs=4, type=float, default=None,
                   metavar=("X0", "X1", "Y0", "Y1"),
                   help="explicit render domain (overrides the auto box), e.g. "
                        "the IC_setup domain: --render-extent -0.5 0.5 -1 1")
    p.add_argument("--render-project",
                   choices=["slice", "column", "rhocolumn", "average"],
                   default=None,
                   help="line-of-sight reduction along the 3rd axis: slice "
                        "(thin mid-plane slab), column (field integral), "
                        "rhocolumn (density-weighted column average), or average "
                        "(plain LOS average). Default: the test's RenderSpec.")
    p.add_argument("--render-slice-frac", type=float, default=None,
                   help="for --render-project slice: slab HALF-thickness as a "
                        "fraction of the 3rd-axis extent (the slab kept is "
                        "|z-zmid| < frac*extent). If given it applies to every "
                        "slab backend (pynbody/hist/particles). If omitted, each "
                        "backend uses its natural default: grid/grid-petkova take "
                        "a true mid-plane (3-D kernel at zmid, no slab); "
                        "pynbody/hist use 0.1; particles show each particle "
                        "within half a smoothing length of the plane "
                        "(|z-zmid| < 0.5*h_i).")
    p.add_argument("--render-rows", type=int, default=None,
                   help="instead of one combined grid (dirs=cols, times=rows), "
                        "make ONE figure per time with the directories tiled into "
                        "this many rows (e.g. --render-rows 2 -> 2 rows of folders)")
    p.add_argument("--render-aspect", choices=["equal", "auto"], default="equal",
                   help="panel aspect ratio: 'equal' (default; 1 data-unit in x "
                        "= 1 in y, so circles look circular) or 'auto' (stretch "
                        "each panel to fill its box)")
    p.add_argument("--render-nsmooth", "--ns", type=int, default=None,
                   metavar="N",
                   help="(grid backend only) neighbour number for the smoothing "
                        "length: by default the grid backend uses the run's "
                        "'smoothlength' aux when present, else getsmooth2 with "
                        "N=64. Setting --ns OVERRIDES that and forces h = "
                        "getsmooth2(mass, rho, N) for every particle, e.g. "
                        "--ns 128 for a smoother map than the run used.")
    return p


def params_from_args(args, param_spec):
    """Extract the {name: value} parameter dict from parsed args."""
    return {spec.name: getattr(args, spec.name) for spec in param_spec}


def check_residual_args(parser, args):
    """Validate `--residual` right after parsing (run_cli / standalone_cli):
    it compares exactly two things -- one run vs the blessed reference, or two
    runs against each other -- so more than two inputs is an error, and a
    single input needs --reference. Implies --render (the residual IS a
    render figure)."""
    if not getattr(args, "residual", False):
        return
    if len(args.inputs) > 2:
        parser.error(f"--residual compares at most TWO inputs (one run vs "
                     f"--reference, or two runs against each other); got "
                     f"{len(args.inputs)}")
    if len(args.inputs) == 1 and args.reference is None:
        parser.error("--residual with a single input needs --reference "
                     "(or give a second input to difference against)")
    args.render = True


# --------------------------------------------------------------------------- #
#  Series extraction
# --------------------------------------------------------------------------- #

def run_series(input_arg, analyze, params, seed=None):
    """Run `analyze` over every snapshot matched by `input_arg`.

    Returns an ordered dict {key: np.ndarray} (sorted by the "x" key) or None
    if nothing could be processed. `analyze` returns a flat dict per snapshot,
    including an "x" key, or None to skip that snapshot.

    `seed` : optional flat dict (same keys as analyze rows) prepended before
    sorting -- e.g. an analytic initial-condition anchor point at x=0.
    """
    try:
        files = find_files(input_arg)
    except FileNotFoundError as e:
        print(f"Input '{input_arg}': {e}", file=sys.stderr)
        return None

    print(f"Input '{input_arg}': {len(files)} snapshot(s); reading in order...")
    rows = [dict(seed)] if seed is not None else []
    for fn in files:
        try:
            tgdata, _td, _ts, _hdr, time, _N, _ng, _nd, _ns, _h = loaddata(fn)
        except Exception as e:
            print(f"  failed to load {fn}: {e}", file=sys.stderr)
            continue
        try:
            row = analyze(tgdata, time, params, fn)
        except Exception as e:
            print(f"  analyze failed for {fn}: {e}", file=sys.stderr)
            continue
        if row is not None:
            rows.append(row)

    if not rows:
        print(f"Input '{input_arg}': no snapshots successfully processed",
              file=sys.stderr)
        return None

    keys = rows[0].keys()
    table = {k: np.array([r[k] for r in rows], dtype=float) for k in keys}
    order = np.argsort(table["x"])
    return {k: v[order] for k, v in table.items()}


def apply_normalization(table, metrics):
    """In-place: divide each metric with a `norm` by that key's first value.

    Denominators are captured from the original table first, so a metric can be
    normalized by another that is itself being normalized (e.g. both M_cloud and
    M_mix divided by M_cloud,0).
    """
    denoms = {}
    for m in metrics:
        if m.norm and m.norm in table and len(table[m.norm]):
            denoms[m.key] = float(table[m.norm][0])
    for m in metrics:
        d = denoms.get(m.key)
        if d:  # skip None and 0
            table[m.key] = table[m.key] / d


def parse_resolution(input_arg):
    """Linear resolution n_x from a run-folder name: the digits immediately after
    the leading test-name (runtest builds runname = '<test><$2>_N<$1>...', where
    $2 is the setup nx passed to IC_createsetup), e.g. 'alfven64_N32...' -> 64,
    'shocktube256_N64_..._Sod_shock' -> 256. None if absent.

    This is the leading '<test><nx>' number ($2 = the resolution that varies in a
    convergence series), NOT the later '_N<...>' ($1) token."""
    base = os.path.basename(os.path.normpath(input_arg))
    m = re.match(r"[A-Za-z]+(\d+)", base)
    return int(m.group(1)) if m else None


def binned_profile(x, qs, nbins=64, lo=None, hi=None, percentile=None,
                   edges=None, stat="mean", drop_empty=True, min_bins=2,
                   clip_outside=True, geom_centres=False):
    """Binned profile(s) of per-particle data on a fixed grid -- THE shared
    binning helper (replaces the per-script _binned/_binned_profile/
    _binned_multi/binned copies).

    x   : per-particle coordinate.
    qs  : ONE array, a list/tuple of arrays, or a dict {key: array}; all share
          the same bin assignment. The returned `values` mirrors the container.
    Grid priority: explicit `edges` > [lo, hi] > `percentile=(p_lo, p_hi)` of x
          > the full [min(x), max(x)] range.
    stat : "mean" or "median" per bin.
    clip_outside : True lumps particles beyond the grid into the edge bins
          (np.clip on the bin index); False EXCLUDES them (the digitize-without-
          clip semantics of the fixed log-grid regression tables).
    drop_empty : True drops unpopulated bins (NaN-free output; returns
          (None, None) if fewer than `min_bins` bins are populated); False keeps
          the full fixed-length grid with NaN in empty bins (element-wise
          regression comparison on a common grid).
    geom_centres : geometric bin centres sqrt(e_i e_i+1) (log-spaced grids).

    Returns (centres, values), or (None, None) on a degenerate grid."""
    x = np.asarray(x, dtype=float)
    if edges is not None:
        edges = np.asarray(edges, dtype=float)
        nbins = len(edges) - 1
    else:
        if lo is None or hi is None:
            if percentile is not None:
                glo, ghi = np.percentile(x, list(percentile))
            else:
                glo, ghi = float(np.min(x)), float(np.max(x))
            lo = glo if lo is None else lo
            hi = ghi if hi is None else hi
        if not (hi > lo):
            return None, None
        edges = np.linspace(lo, hi, nbins + 1)
    centres = (np.sqrt(edges[:-1] * edges[1:]) if geom_centres
               else 0.5 * (edges[:-1] + edges[1:]))

    if isinstance(qs, dict):
        keys, arrays, container = list(qs.keys()), list(qs.values()), "dict"
    elif isinstance(qs, (list, tuple)):
        keys, arrays, container = None, list(qs), "list"
    else:
        keys, arrays, container = None, [qs], "single"
    arrays = [np.asarray(a, dtype=float) for a in arrays]

    idx = np.digitize(x, edges) - 1
    if clip_outside:
        idx = np.clip(idx, 0, nbins - 1)
        valid = np.ones(len(x), dtype=bool)
    else:
        valid = (idx >= 0) & (idx < nbins)
    fstat = np.median if stat == "median" else np.mean
    outs = [np.full(nbins, np.nan) for _ in arrays]
    pop = np.zeros(nbins, dtype=bool)
    for b in range(nbins):
        mb = valid & (idx == b)
        if np.any(mb):
            pop[b] = True
            for o, a in zip(outs, arrays):
                o[b] = fstat(a[mb])
    if drop_empty:
        if pop.sum() < min_bins:
            return None, None
        centres = centres[pop]
        outs = [o[pop] for o in outs]
    if container == "dict":
        return centres, dict(zip(keys, outs))
    if container == "list":
        return centres, outs
    return centres, outs[0]


def reg_nbins(input_arg, fallback=64, lo=16, hi=512):
    """Resolution-aware bin count for a binned profile: ~n_x (the run's setup
    resolution parsed from its folder name) so the reference/overlay curve tracks
    the native resolution -- one bin per resolution element rather than a fixed
    grid that over-smooths hi-res runs and under-fills lo-res ones. Clipped to
    [lo, hi]; `fallback` (64) when n_x can't be read from the name."""
    nx = parse_resolution(input_arg)
    if nx is None:
        return fallback
    return int(np.clip(nx, lo, hi))


def series_from_log_or_snapshots(inp, log_cols, per_snapshot, valid=None):
    """A time series for run `inp`, PREFERRING the dense gasoline '.log' and
    falling back to per-snapshot computation -- THE shared helper for the
    read-dense-log-else-recompute pattern (mhdloop/orzag/mhdrotor/gresho/sedov/
    evrard/... each hand-rolled this).

    log_cols     : LOG_COLUMNS names to read; log_cols[0] is the time column
                   ("dTime"). Read via `read_log_series`.
    valid(log)   : whether the '.log' series is usable. Default: the FIRST
                   non-time column has any value > 0. (Pass e.g.
                   `lambda lg: np.any(lg["Etot"] != 0)` for a signed total.)
    per_snapshot(fn) -> (time, v1, v2, ...) : the fallback computed per snapshot
                   (one row, same column order as log_cols). Snapshots come from
                   `find_files(inp)`; rows are sorted by time.

    Returns `(col0, col1, ..., source)` -- one array per `log_cols` entry plus
    `source` in {"log", "snapshots"} -- or `(None,)*(len(log_cols)+1)` if neither
    source yields data. So a 2-column call returns (t, v, source); 3-column
    returns (t, v1, v2, source); etc."""
    log = read_log_series(inp, list(log_cols))
    ok = (valid(log) if valid is not None
          else np.any(np.asarray(log[log_cols[1]]) > 0)) if log is not None else False
    if ok:
        return (*(np.asarray(log[c]) for c in log_cols), "log")
    try:
        files = find_files(inp)
    except FileNotFoundError:
        files = []
    rows = [tuple(per_snapshot(fn)) for fn in files]
    if not rows:
        return (*(None for _ in log_cols), None)
    a = np.array(sorted(rows, key=lambda r: r[0]), dtype=float)
    return (*(a[:, j] for j in range(a.shape[1])), "snapshots")


# --------------------------------------------------------------------------- #
#  Plotting
# --------------------------------------------------------------------------- #

def _legend_figure(handles, labels):
    """Build a standalone figure holding just the legend (no axes)."""
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(6, 0.5 + 0.35 * len(labels)))
    fig.legend(handles, labels, loc="center", ncol=min(3, len(labels)),
               frameon=False, handlelength=2.5, handletextpad=0.8,
               columnspacing=1.5)
    fig.tight_layout(pad=0.1)
    return fig


def plot_compare(results, ref_table, analytic, params, metrics, save=None):
    """One figure per metric (+ a separate legend figure) for each metric.

    results : list of (label, table, ...) for each input. The legend is placed in
    its own figure rather than on the axes (publication-friendly).
    """
    import matplotlib.pyplot as plt
    setup_rcparams()

    for met in metrics:
        fig, ax = plt.subplots(figsize=(9, 6))
        plotfn = ax.semilogy if met.yscale == "log" else ax.plot

        for label, table, *_ in results:
            if met.key in table:
                plotfn(table["x"], table[met.key], marker=met.marker,
                       linestyle="-", label=label)

        if ref_table is not None and met.key in ref_table:
            plotfn(ref_table["x"], ref_table[met.key], marker="x",
                   linestyle="--", color="0.4", label="reference")

        if analytic is not None:
            curve = analytic(met.key, params)
            if curve is not None:
                xa, ya = curve
                plotfn(xa, ya, "k:", linewidth=2.0, label="analytic")

        ax.set_xlabel(met.xlabel)
        ax.set_ylabel(met.ylabel)
        fig.tight_layout()

        # Legend goes in its own figure, not on the axes.
        handles, labels = ax.get_legend_handles_labels()
        legfig = _legend_figure(handles, labels) if handles else None

        if save:
            out = f"{save}_{met.key}.png"
            fig.savefig(out, bbox_inches="tight")
            print(f"saved figure -> {out}")
            if legfig is not None:
                legout = f"{save}_{met.key}_legend.png"
                legfig.savefig(legout, bbox_inches="tight")
                print(f"saved legend -> {legout}")
                plt.close(legfig)
            plt.close(fig)
        else:
            plt.show()


def select_at_time(input_arg, t, test="analysis", warn_outside=True):
    """The single snapshot of `input_arg` nearest time `t`, as (filename, time).

    THE shared single-snapshot selector (replaces the ~12 per-script copies).
    `t is None` -> the LAST snapshot. Otherwise the nearest dump to `t`; when
    `warn_outside` and `t` falls outside the run's snapshot span, a note is
    printed to stderr (the comparison happens at the nearest available time, not
    the requested one). `test` only labels the stderr messages. Returns None
    (after a message) when the run has no snapshots."""
    files = find_files(input_arg)
    if not files:
        print(f"{test}: no snapshots for '{input_arg}'", file=sys.stderr)
        return None
    if t is None:
        return select_snapshots(files, n=len(files))[-1]
    pairs = select_snapshots(files, n=len(files))
    tmin, tmax = min(p[1] for p in pairs), max(p[1] for p in pairs)
    fn, tt = min(pairs, key=lambda p: abs(p[1] - t))
    if warn_outside and (t < tmin - 1e-9 or t > tmax + 1e-9):
        print(f"{test}: requested t={t:g} -> using nearest snapshot t={tt:g} "
              f"(run only spans [{tmin:g},{tmax:g}]).", file=sys.stderr)
    return fn, tt


def render_reference_path_for(input_arg, slug, plane="z0"):
    """Render-reference path inside the run's reference folder:
    '<references>/<runname>/render_<slug>_<plane>.npz'.

    A single run produces SEVERAL render figures (one per quantity, and per plane
    in 3-D), so each gets its own npz keyed by `slug`/`plane` within the shared
    'references/<runname>/' folder. Independent of the user's `--save` output
    prefix, so the baseline is discoverable from the run alone (auto) or from an
    explicit reference folder."""
    return os.path.join(reference_dir_for(input_arg),
                        f"render_{slug}_{plane}.npz")



# --------------------------------------------------------------------------- #
#  1-D horizontal cuts (field vs x along a line y=y0)
# --------------------------------------------------------------------------- #

def mhd_field_cut_quants():
    """Standard (slug, label, qfn) set for a 1-D MHD field cut: density rho, the
    signed field components Bx, By, and the signed divergence div(B). Each qfn is
    `qfn(fn, tgdata, ngas) -> per-particle array`, reusing the aux factories."""
    return [
        ("rho", r"$\rho$", lambda fn, tg, ng: tg[:, 7]),
        ("Bx", r"$B_x$", aux_component("BField", 0)),
        ("By", r"$B_y$", aux_component("BField", 1)),
        ("divB", r"$\nabla\!\cdot\!B$", aux_scalar("DivB")),
    ]


def horizontal_cut(inputs, labels, quants, times, y0_of, save=None,
                   time_label="t", title_prefix="", axis_x=1, axis_y=2,
                   strip_cells=1.0):
    """1-D cuts of fields vs x along the horizontal line y=y0, one figure per time
    (one stacked panel per quantity).

    quants  : [(slug, ylabel, qfn(fn, tgdata, ngas))] -- the fields to cut.
    times   : the times to cut at (nearest snapshot chosen).
    y0_of   : callable t -> y0 (the cut height); annotated in the title.
    Particles within |y - y0| < strip_cells / n_x (n_x from the run name) are
    scattered vs x and overlaid with a binned-mean line. axis_x/axis_y are the
    tgdata column indices of the in-plane coordinates (default x=1, y=2)."""
    import matplotlib.pyplot as plt
    setup_rcparams()
    xname = "xyz"[axis_x - 1]
    for t in times:
        y0 = y0_of(t)
        fig, axes = plt.subplots(len(quants), 1,
                                 figsize=(8, 2.2 * len(quants) + 1.0),
                                 sharex=True, squeeze=False)
        axes = axes[:, 0]
        any_data = False
        for ci, (inp, lab) in enumerate(zip(inputs, labels)):
            try:
                files = find_files(inp)
            except FileNotFoundError:
                continue
            chosen = select_snapshots(files, times=[t], n=1)
            if not chosen:
                continue
            fn, _tt = chosen[0]
            tgdata, _td, _ts, hdr, _tm = _cached_readtipsy(fn)
            ngas = int(hdr[2])
            nx = parse_resolution(inp) or 64
            strip = np.abs(tgdata[:, axis_y] - y0) < strip_cells / nx
            x = np.asarray(tgdata[strip, axis_x], dtype=float)
            order = np.argsort(x)
            x = x[order]
            for qi, (slug, ylabel, qfn) in enumerate(quants):
                q = np.asarray(qfn(fn, tgdata, ngas), dtype=float)[strip][order]
                axes[qi].scatter(x, q, s=4, alpha=0.4, color=f"C{ci}",
                                 rasterized=True,
                                 label=lab if qi == 0 else None)
                if x.size > 8:                       # binned-mean line
                    xc, qm = binned_profile(x, q, nbins=min(nx, 128), min_bins=1)
                    if xc is not None:
                        axes[qi].plot(xc, qm, "-", color=f"C{ci}", lw=1.5)
                axes[qi].set_ylabel(ylabel)
            any_data = True
        for ax in axes:
            ax.axhline(0.0, color="0.85", lw=0.8, zorder=0)
        axes[-1].set_xlabel(xname)
        axes[0].set_title(f"{title_prefix}horizontal cut at "
                          f"{'xyz'[axis_y - 1]}0 = {y0:g}  ({time_label} = {t:g})")
        if any_data and len(inputs) > 1:
            axes[0].legend(fontsize=8, loc="best")
        fig.tight_layout()
        if not any_data:
            print(f"horizontal_cut: no data at {time_label}={t:g}",
                  file=sys.stderr)
        if save:
            out = f"{save}_cut_t{t:g}.png"
            fig.savefig(out, bbox_inches="tight")
            print(f"saved cut -> {out}")
            plt.close(fig)
        else:
            plt.show()


def bless_render_specs(inp, specs, args, test, plane="z0", extent=None,
                       times=None):
    """Bless run `inp`'s render-grid baselines -- one npz per RenderSpec slug, in
    the named `plane` -- into the run's references/<runname>/ folder, so a test's
    reference set carries the visual render baseline alongside the metric one
    (the shared helper behind run_cli's and the standalone render tests' bless;
    mirrors orzag).

    `specs` are RenderSpecs (each disambiguated by spec_slug). `times` are the
    times to bless at -- pass the SAME times the test renders at so the overlay
    (matched by time) has a grid to compare against; None -> --render-times. With
    --save-reference-all the grids are blessed at EVERY dump time instead (a
    time-agnostic reference, needed for --movie --reference). The plain flag warns
    if a requested render time has no snapshot. `extent` overrides --render-extent
    for the bless (e.g. a test's fixed IC box); None falls through to
    --render-extent. Scatter-only backends can't be blessed (no grid) and are
    skipped by render_panels with a warning."""
    if not specs:
        return
    want = times if times is not None else args.render_times
    if args.save_reference_all:
        bless_times = snapshot_times(inp) or want
    else:
        bless_times = want
        if bless_times:
            miss = missing_render_times(inp, bless_times)
            if miss:
                print(f"[{test}] --save-reference: render time(s) "
                      f"{', '.join(f'{t:g}' for t in miss)} have no snapshot in "
                      f"'{inp}'; rendering the nearest dumps for the baseline.",
                      file=sys.stderr)
    label = os.path.basename(os.path.normpath(inp))
    ext = extent if extent is not None else args.render_extent
    for sp in specs:
        ref_path = render_reference_path_for(inp, spec_slug(sp), plane)
        render_panels([(label, inp)], sp,
                      times=bless_times, n=args.render_n, res=args.render_res,
                      backend=args.render_backend, project=args.render_project,
                      slice_frac=args.render_slice_frac,
                      extent=ext, aspect=args.render_aspect,
                      nsmooth=args.render_nsmooth,
                      save_render_reference=ref_path)


def render_ref_overlay_path(inp, reference, slug, plane="z0"):
    """Path of a blessed render-grid npz to OVERLAY on a fresh render, or None.

    Opt-in: `reference` None -> None (no overlay). Otherwise resolve the
    reference folder (REFERENCE_AUTO -> the run's own; a dir/file -> that set) and
    return '<dir>/render_<slug>_<plane>.npz' if it exists, else None (a warning is
    left to the caller, which knows the test name)."""
    if reference is None:
        return None
    ref_dir = resolve_reference_dir(inp, reference)
    if ref_dir is None:
        return None
    cand = os.path.join(ref_dir, f"render_{slug}_{plane}.npz")
    return cand if os.path.isfile(cand) else None


def movie_render_specs(inp, label, specs, args, save, test, plane="z0",
                       extent=None, time_to_x=None, time_label="t", params=None):
    """Animate each RenderSpec over ALL snapshots of run `inp` -- one movie per
    spec slug -- the shared helper behind run_cli's and the standalone render
    tests' `--movie` (mirrors orzag). Movies are single-run.

    With `--reference` each frame is sim | reference | (sim - reference) when the
    blessed render grids cover EVERY frame time (bless with --save-reference-all);
    render_movie falls back to a sim-only movie with a warning otherwise. `extent`
    overrides --render-extent; None falls through. Image backends only (a
    `particles`/scatter backend has no grid to animate -- render_movie warns)."""
    if not specs:
        print(f"[{test}] --movie given but this test defines no RenderSpec; "
              f"skipping.", file=sys.stderr)
        return
    ext = extent if extent is not None else args.render_extent
    overlay = args.reference is not None
    for sp in specs:
        slug = spec_slug(sp)
        mv_save = f"{save}_{slug}_movie" if save else None
        ref_path = (render_ref_overlay_path(inp, args.reference, slug, plane)
                    if overlay else None)
        render_movie(inp, label, sp, save=mv_save,
                     res=args.render_res, backend=args.render_backend,
                     project=args.render_project,
                     slice_frac=args.render_slice_frac, extent=ext,
                     aspect=args.render_aspect, time_to_x=time_to_x,
                     time_label=time_label, params=params,
                     nsmooth=args.render_nsmooth, render_reference=ref_path,
                     fps=args.movie_fps, fmt=args.movie_format)


def run_cli(test, param_spec, analyze, analytic, metrics, description="",
            seed=None, render_spec=None, stream_spec=None, time_to_x=None,
            time_label="t"):
    """Entry point a test's main() delegates to.

    test        : str, e.g. "kh" (used for reference metadata / messages).
    param_spec  : list[Param].
    analyze     : callable(tgdata, time, params, snap) -> dict | None.
    analytic    : callable(metric_key, params) -> (x, y) | None.
    metrics     : list[Metric].
    seed        : optional initial-condition anchor row (see run_series).
    render_spec : enables `--render` field maps. One RenderSpec, a list of them
                  (one figure each, output files suffixed by the spec label), or a
                  callable render_spec(params) -> RenderSpec | list | None so a
                  test can pick what to render from the run params (e.g. kh adds
                  magnetic-field maps only when B0>0).
    stream_spec : streamline figures drawn with `--render` alongside the field
                  maps. One StreamSpec, a list, or a callable
                  stream_spec(params) -> StreamSpec | list | None (e.g. kh adds
                  B streamlines only when B0>0). Output files are suffixed
                  '_stream_<labelslug>'.
    time_to_x   : optional callable(raw_time, params) -> x, so `--render-times`
                  is expressed in the analysis x-axis (e.g. t/t_crush for blob).
    time_label  : axis label used in render panel titles (e.g. "t/t_crush").
    """
    parser = build_parser(param_spec, description=description)
    # CI gate (added here, not in build_parser, so it doesn't collide with the
    # standalone analyses that declare their own --reg-tol).
    parser.add_argument("--reg-tol", type=float, default=None,
                        help="CI gate: exit non-zero if the worst residual L2 "
                             "vs the regression baseline exceeds this tolerance")
    args = parser.parse_args()
    check_residual_args(parser, args)
    params = params_from_args(args, param_spec)
    if args.no_analytic:
        analytic = None

    print(f"[{test}] parameters: " +
          "  ".join(f"{k}={v}" for k, v in params.items()))

    results = []  # (label, table, input_arg)
    for i, inp in enumerate(args.inputs):
        table = run_series(inp, analyze, params, seed=seed)
        if table is None:
            continue
        label = (args.labels[i] if args.labels and i < len(args.labels)
                 else os.path.basename(os.path.normpath(inp)))
        apply_normalization(table, metrics)
        results.append((label, table, inp))

    if not results:
        print("No data to plot; exiting.", file=sys.stderr)
        sys.exit(1)

    # Resolve the render field-map specs once (a single spec, a list, or a
    # callable spec(params) -> spec | list | None). Reused for blessing render
    # grids (--save-reference) and for drawing/overlaying them (--render).
    def _resolved_specs():
        sp = render_spec(params) if callable(render_spec) else render_spec
        if sp is None:
            return []
        if isinstance(sp, RenderSpec):
            return [sp]
        return list(sp)

    # --save-reference[-all]: bless the FIRST input as the baseline and stop.
    # Besides the metric series this also blesses the render-grid baselines (one
    # npz per RenderSpec) + a copy of the run .log, mirroring the standalone
    # render tests, so a run_cli test's reference set is complete for both the
    # numeric gate and the visual render-reference overlay.
    if args.save_reference or args.save_reference_all:
        ref_path = reference_path_for(args.inputs[0])
        save_reference(results[0][1], ref_path, params, test)
        copy_run_log(args.inputs[0], reference_dir_for(args.inputs[0]))
        bless_render_specs(args.inputs[0], _resolved_specs(), args, test)
        return

    # Normal run: discover (opt-in via --reference) the reference for the first
    # input, report residuals, and overlay everything.
    ref_table, _ref_path = find_reference(args.inputs[0], explicit=args.reference)
    worst_l2 = None
    if ref_table is not None:
        label0, table0, _ = results[0]
        print(f"residuals ({label0} vs reference):")
        for met in metrics:
            stat = residuals(table0, ref_table, met.key)
            if stat is not None:
                mx, l2 = stat
                worst_l2 = l2 if worst_l2 is None else max(worst_l2, l2)
                print(f"  {met.key}: max|res|={mx:.4g}  L2={l2:.4g}")

    plot_compare(results, ref_table, analytic, params, metrics, save=args.save)

    # Optional: render field maps. All inputs go into ONE figure (directories as
    # columns, times as rows when there are several).
    if args.render:
        # render_spec/stream_spec may each be: a single spec, a list of them, or
        # a callable spec(params) -> spec | list | None (so a test can choose
        # what to render from the run params, e.g. kh adds B maps/streamlines
        # only if B0>0).
        specs = _resolved_specs()
        sspecs = stream_spec(params) if callable(stream_spec) else stream_spec
        if sspecs is None:
            sspecs = []
        elif isinstance(sspecs, StreamSpec):
            sspecs = [sspecs]
        entries = [(label, inp) for (label, _table, inp) in results]
        if not specs and not sspecs:
            print("--render given but this test defines no RenderSpec; skipping.",
                  file=sys.stderr)
        multi = len(specs) > 1
        # Render-reference overlay is opt-in (--reference) and only for a single
        # input (the sim|reference|difference grid compares one run to a baseline).
        overlay = args.reference is not None and len(results) == 1
        for sp in specs:
            # Disambiguate output files per quantity when several are rendered.
            sp_save = (f"{args.save}_{spec_slug(sp)}"
                       if args.save and multi else args.save)
            ref_path = (render_ref_overlay_path(args.inputs[0], args.reference,
                                                spec_slug(sp))
                        if overlay else None)
            if overlay and ref_path is None:
                print(f"[{test}] no render reference for {spec_slug(sp)}",
                      file=sys.stderr)
            render_panels(entries, sp,
                          times=args.render_times, n=args.render_n,
                          res=args.render_res, save=sp_save,
                          backend=args.render_backend,
                          extent=args.render_extent,
                          project=args.render_project,
                          slice_frac=args.render_slice_frac,
                          time_to_x=time_to_x, time_label=time_label,
                          params=params, render_rows=args.render_rows,
                          aspect=args.render_aspect,
                          nsmooth=args.render_nsmooth,
                          render_reference=ref_path, residual=args.residual)
        for sp in sspecs:
            # Streamline grids read best at moderate resolution; cap at 128.
            sp_save = (f"{args.save}_stream_{_label_slug(sp.label)}"
                       if args.save else None)
            streamline_panels(entries, sp,
                              times=args.render_times, n=args.render_n,
                              res=min(args.render_res or 64, 128),
                              save=sp_save, extent=args.render_extent,
                              time_label=time_label)

    # Optional: animate each field map over ALL snapshots (one movie per
    # quantity). Movies are single-run -> the FIRST input; with --reference each
    # frame is sim|reference|difference when the baseline covers every frame time.
    if args.movie:
        movie_render_specs(args.inputs[0], results[0][0], _resolved_specs(), args,
                           args.save, test, time_to_x=time_to_x,
                           time_label=time_label, params=params)

    # CI gate: all plots/renders are produced first, then the exit code signals
    # pass/fail (worst residual L2 of the first input vs its blessed baseline).
    if args.reg_tol is not None and worst_l2 is not None and worst_l2 > args.reg_tol:
        print(f"[{test}] REGRESSION FAILED -- worst L2 {worst_l2:.4g} > tol "
              f"{args.reg_tol:.4g}", file=sys.stderr)
        sys.exit(2)


# --------------------------------------------------------------------------- #
#  Standalone CLI scaffold (the analogue of run_cli for the standalone scripts
#  that build their own figures rather than the metric-vs-time plot_compare)
# --------------------------------------------------------------------------- #

def compare_to_baseline(test, inp, ref_table_fn, reg_keys=None, explicit=None,
                        ref_table=None):
    """Generic input-vs-baseline residual report -- replaces the ~19 near-identical
    per-script `compare_to_baseline` copies.

    `ref_table_fn()` is a zero-arg closure returning this run's regression table
    ({"x": grid, key: series, ...}); the script binds its own `inp`/time/params.
    `reg_keys` selects which series to compare (default: every key except "x").
    Discovers the baseline with `find_reference(inp, explicit)`, prints
    `<test>[regression] <key> vs baseline: max|d|=… rms=…` per key, and returns
    the worst RMS (None if no baseline / no table / no comparable key).

    `ref_table` may be passed pre-discovered (e.g. standalone_cli already found
    the first input's baseline for the overlay) to avoid re-reading it and
    re-printing 'reference found'."""
    if ref_table is None:
        ref_table, _p = find_reference(inp, explicit)
    if ref_table is None:
        return None
    table = ref_table_fn()
    if table is None:
        return None
    keys = reg_keys if reg_keys is not None else [k for k in table if k != "x"]
    worst, found = 0.0, False
    for key in keys:
        res = residuals(table, ref_table, key)
        if res is None:
            continue
        mx, rms = res
        worst = max(worst, rms)
        found = True
        print(f"{test}[regression] {key} vs baseline: "
              f"max|d|={mx:.4g}  rms={rms:.4g}")
    return worst if found else None


def standalone_cli(test, param_spec, plot, reference_table=None, *,
                   description="", reg_keys=None, add_time=True,
                   time_default=None, times_default=None, time_help=None,
                   add_arguments=None, ref_params=None, compare=None,
                   save_extra=None):
    """The standalone analogue of `run_cli`: absorbs the ~40-line main()
    boilerplate (parser + labels + `--save-reference` bless-and-exit +
    compare/`--reg-tol` gate) so each script's `main` is one call plus its
    test-specific `plot`/`reference_table` callbacks.

    test            : name for messages and the saved reference.
    param_spec       : the -a..-h Param list.
    plot(inputs, labels, params, args)            -- draw every figure (the only
                       test-specific code; branch on params/args as needed).
    reference_table(inp, params, args) -> table | None
                       -- this run's regression table; None disables the
                       bless/compare/gate entirely (a plot-only analysis).
    reg_keys         : which table keys to gate on (default: all but "x").
    add_time         : add the standard comparison-time flag `--times` (nargs +,
                       one profile panel each where supported -- a single
                       comparison time is just `--times <t>`). It is resolved so
                       the callbacks ALWAYS see both `args.times` (list) and
                       `args.time` (its first entry, the scalar). False for tests
                       with no comparison time.
    time_default     : the single-time default (-> `--times` default [t]).
    times_default    : the multi-time default list (e.g. [1.0, 3.0] for gresho);
                       takes precedence over time_default for the `--times`
                       default. `args.time` then defaults to its first entry.
    time_help        : help text for `--times`.
    add_arguments(parser)  : add test-specific CLI flags (--convergence, --y0,
                       --refdir, ...).
    ref_params(params, args) -> dict : the params recorded in the blessed JSON
                       (default: the params dict).
    compare(inp, params, args) -> worst_rms | None : a custom per-input
                       comparison, used INSTEAD of the generic residual loop (for
                       tests with a bespoke residual, e.g. evrard/polytrope's
                       element-wise NaN-grid `reg_residual`). reference_table is
                       still used for `--save-reference`.
    save_extra(args, params)  : called on `--save-reference[-all]` AFTER the
                       metric bless (for render tests: bless the render-grid
                       baselines via bless_render_specs + copy the run .log) so a
                       standalone render test's reference set is complete."""
    parser = build_parser(param_spec, description=description)
    if add_time:
        dt = (times_default if times_default is not None else
              ([time_default] if time_default is not None else None))
        parser.add_argument("--times", type=float, nargs="+", default=dt,
                            help=(time_help or
                                  "comparison time(s); one profile panel each "
                                  "where the test supports it -- pass a single "
                                  "value for one time")
                                 + f" (default {dt})")
    parser.add_argument("--reg-tol", type=float, default=None,
                        help="CI gate: exit 2 if the worst RMS residual vs the "
                             "regression baseline exceeds this tolerance")
    if add_arguments is not None:
        add_arguments(parser)
    args = parser.parse_args()
    check_residual_args(parser, args)
    params = params_from_args(args, param_spec)

    # `--times` is the canonical comparison-time list; expose its first entry as
    # the scalar `args.time` for the single-time callbacks (so both are always
    # available). A single comparison time is just `--times <t>`.
    if add_time:
        args.time = args.times[0] if args.times else None

    labels = args.labels if args.labels else [
        os.path.basename(os.path.normpath(i)) for i in args.inputs]
    if len(labels) != len(args.inputs):
        parser.error("--labels must match the number of inputs")

    # --save-reference[-all]: bless the FIRST input and stop.
    if args.save_reference or args.save_reference_all:
        if reference_table is None:
            parser.error("this analysis has no regression baseline to save")
        table = reference_table(args.inputs[0], params, args)
        if table is None:
            print(f"{test}: nothing to save as a reference for "
                  f"'{args.inputs[0]}'", file=sys.stderr)
            sys.exit(1)
        rp = ref_params(params, args) if ref_params else params
        save_reference(table, reference_path_for(args.inputs[0]), rp, test)
        # Render tests additionally bless their render-grid baselines + the .log.
        if save_extra is not None:
            copy_run_log(args.inputs[0], reference_dir_for(args.inputs[0]))
            save_extra(args, params)
        return args

    # Reference comparison is OPT-IN (only when --reference is given). Discover
    # the FIRST input's baseline ONCE here -- it is stashed on `args.ref_table`
    # so the test's `plot` can OVERLAY it on every non-render figure, and reused
    # for the first input's regression gate below (so the baseline is read and
    # 'reference found' printed only once). The resolved comparison time is passed
    # so a time mismatch vs the blessed reference warns. args.ref_table is None
    # when --reference is absent or the test has no regression baseline.
    args.ref_table = None
    if args.reference is not None and reference_table is not None:
        tcmp = getattr(args, "time", None)
        if add_time:
            sel0 = select_at_time(args.inputs[0], tcmp, test=test)
            if sel0 is not None:
                tcmp = sel0[1]
        args.ref_table, _rp = find_reference(args.inputs[0],
                                             explicit=args.reference, time=tcmp)

    plot(args.inputs, labels, params, args)

    # Compare each input to its baseline, then the CI gate. The first input reuses
    # the already-discovered args.ref_table; the others discover their own sibling.
    worst = None
    if compare is not None:
        for inp in args.inputs:
            w = compare(inp, params, args)
            if w is not None:
                worst = w if worst is None else max(worst, w)
    elif reference_table is not None:
        for i, inp in enumerate(args.inputs):
            w = compare_to_baseline(
                test, inp,
                lambda inp=inp: reference_table(inp, params, args),
                reg_keys=reg_keys, explicit=args.reference,
                ref_table=(args.ref_table if i == 0 else None))
            if w is not None:
                worst = w if worst is None else max(worst, w)
    if args.reg_tol is not None and worst is not None and worst > args.reg_tol:
        print(f"{test}: REGRESSION FAILED -- worst rms {worst:.4g} > tol "
              f"{args.reg_tol:.4g}", file=sys.stderr)
        sys.exit(2)
    return args
