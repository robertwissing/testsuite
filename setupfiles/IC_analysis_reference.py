#!/usr/bin/env python3
"""
Reference (regression-baseline) I/O for the analysis framework.

The metric-reference subsystem extracted from IC_analysis_framework so the
framework module is smaller and the reference policy lives in one place. This
module is self-contained (depends only on IC_analysis_general + numpy/stdlib);
IC_analysis_framework re-exports every public name here, so callers keep doing
`from IC_analysis_framework import save_reference, find_reference, ...`.

Layout of a run's reference set (all artifacts together, per configuration):
  <parent>/references/<runname>/
      reference.json              -- the blessed metric series (this module)
      render_<slug>_<plane>.npz   -- render-grid baselines (IC_analysis_render)
      <runname>.log               -- a copy of the run log (copy_run_log)
"""

import json
import os
import shutil
import sys

import numpy as np

from IC_analysis_general import _find_logs

# Sentinel for `--reference` given WITHOUT a value: "compare, auto-discovering the
# run's references/<runname>/ folder". `--reference <path>` carries the path;
# `--reference` absent is None -> NO reference comparison (opt-in, not automatic).
REFERENCE_AUTO = object()


# --------------------------------------------------------------------------- #
#  Reference-set paths
# --------------------------------------------------------------------------- #

def reference_dir_for(input_arg):
    """Folder holding the full reference set for a run directory / prefix.

    All reference artifacts for one run live together under
    '<parent>/references/<runname>/' -- the metric 'reference.json', the render
    grid baselines 'render_<slug>_<plane>.npz', and a copy of the run's '.log'
    (so the blessed parameters can be inspected later). The run-folder basename
    already encodes test/resolution/params/code (that's how runtest.sh names it),
    so the per-run subfolder is unique per configuration. Trailing slashes are
    stripped so a directory and 'directory/' map to the same reference.
    """
    base = os.path.normpath(input_arg)
    parent = os.path.dirname(base)
    name = os.path.basename(base)
    return os.path.join(parent, "references", name)


def reference_path_for(input_arg):
    """Metric-reference JSON path for a run: '<references>/<runname>/reference.json'
    (the per-run folder disambiguates, so the filename itself is generic)."""
    return os.path.join(reference_dir_for(input_arg), "reference.json")


def resolve_reference_dir(input_arg, reference):
    """Reference FOLDER to compare `input_arg` against, given a --reference value.

    OPT-IN: None -> None (no reference); REFERENCE_AUTO -> the run's own
    '<references>/<runname>/'; a directory -> itself; a file (e.g. a
    reference.json) -> its parent directory. This lets '--reference' point at
    another run's blessed reference set (folder) for a cross-run comparison."""
    if reference is None:
        return None
    if reference is REFERENCE_AUTO:
        return reference_dir_for(input_arg)
    if os.path.isdir(reference):
        return reference
    return os.path.dirname(reference)


def copy_run_log(input_arg, dest_dir):
    """Copy the run's gasoline '.log' into the reference folder so the blessed
    parameters (box, EOS, viscosity, dTuFac, ...) can be inspected later."""
    logs = _find_logs(input_arg)
    if not logs:
        print(f"reference: no '.log' found for '{input_arg}' to copy",
              file=sys.stderr)
        return
    os.makedirs(dest_dir, exist_ok=True)
    for lg in logs:
        dst = os.path.join(dest_dir, os.path.basename(lg))
        shutil.copy2(lg, dst)
        print(f"copied run log -> {dst}")


# --------------------------------------------------------------------------- #
#  Metric reference: save / load / discover
# --------------------------------------------------------------------------- #

def save_reference(table, path, params, test, time=None):
    """Write a computed series table to `path` as JSON.

    `time` (when given) records the snapshot/comparison time the reference was
    blessed at, so a later compare can warn if it is run at a different time (a
    profile-vs-time test blessed at t=0.1 is meaningless against a run at t=0.2).
    """
    payload = {
        "test": test,
        "params": params,
        "series": {k: v.tolist() for k, v in table.items()},
    }
    if time is not None:
        payload["time"] = float(time)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2)
    msg = f"saved reference -> {path}"
    if time is not None:
        msg += f"  (blessed at t={float(time):g})"
    print(msg)


def load_reference(path):
    """Load a reference JSON, returning (table, payload) or (None, None)."""
    if not path or not os.path.isfile(path):
        return None, None
    with open(path) as fh:
        payload = json.load(fh)
    table = {k: np.asarray(v, dtype=float) for k, v in payload["series"].items()}
    return table, payload


def find_reference(input_arg, explicit=None, time=None):
    """Resolve the metric reference table for an input run, or (None, None).

    OPT-IN policy (uniform across every test): `explicit` may be
      * None            -> NO comparison; return (None, None) silently. This is
                           the default when `--reference` is absent, so a plain
                           run never auto-compares.
      * REFERENCE_AUTO  -> (bare `--reference`) auto-discover the run's own
                           '<references>/<runname>/reference.json'.
      * a directory     -> '<dir>/reference.json' (e.g. another run's set);
      * a file          -> that JSON.
    When `time` is given and the reference recorded its own blessed time, warn if
    they differ (a comparison between profiles at different evolution stages).
    """
    if explicit is None:
        return None, reference_path_for(input_arg)   # opt-in: no auto-compare
    if explicit is REFERENCE_AUTO:
        path = reference_path_for(input_arg)
        was_explicit = False
    elif os.path.isdir(explicit):
        path = os.path.join(explicit, "reference.json")
        was_explicit = True
    else:
        path = explicit
        was_explicit = True
    table, payload = load_reference(path)
    if table is not None:
        print(f"reference found: {path}")
        tref = payload.get("time") if payload else None
        if (time is not None and tref is not None
                and not np.isclose(time, tref, rtol=1e-3, atol=1e-6)):
            print(f"WARNING: comparison time t={float(time):g} differs from the "
                  f"reference's blessed time t={float(tref):g}; the profiles are "
                  f"at different evolution stages, so the overlay/residuals may "
                  f"be misleading. Re-bless at this time or pass the matching "
                  f"--time.", file=sys.stderr)
    elif was_explicit:
        print(f"reference '{explicit}' not found", file=sys.stderr)
    else:
        print(f"no reference found (looked for {path})")
    return table, path


# --------------------------------------------------------------------------- #
#  Residuals (input vs reference)
# --------------------------------------------------------------------------- #

def _monotone_xy(x, y):
    """Sort (x, y) by x and collapse duplicate x (keeping the LAST value), so the
    x-grid is strictly increasing -- the precondition np.interp needs.

    Duplicate x arises e.g. when a run dumps at t=0 and a test prepends a seed
    row at x=0 (kh); keeping the last value picks the real snapshot over the seed
    anchor. For the common already-unique series this is a no-op.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    order = np.argsort(x, kind="stable")
    x, y = x[order], y[order]
    keep = np.ones(len(x), dtype=bool)
    keep[:-1] = x[1:] != x[:-1]            # at equal x, drop all but the last
    return x[keep], y[keep]


def residuals(input_table, ref_table, key):
    """Residual stats of input vs reference for one metric key.

    Interpolates the input metric onto the reference x-grid and returns
    (max_abs, l2) where l2 is the RMS deviation. Returns None if either side
    lacks the metric. Both grids are made strictly-increasing first
    (`_monotone_xy`) so a duplicate x (e.g. a t=0 dump coinciding with a seed
    anchor) can't produce a spurious non-zero residual against identical data.
    """
    if key not in input_table or ref_table is None or key not in ref_table:
        return None
    xr, yr = _monotone_xy(ref_table["x"], ref_table[key])
    xi, yi = _monotone_xy(input_table["x"], input_table[key])
    yi_on_r = np.interp(xr, xi, yi)
    res = yi_on_r - yr
    return float(np.max(np.abs(res))), float(np.sqrt(np.mean(res ** 2)))
