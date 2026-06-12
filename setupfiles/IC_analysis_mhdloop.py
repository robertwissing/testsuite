#!/usr/bin/env python3
"""
Analysis for the MHD field-loop advection test (IC_setup_mhdloop).

The Gardiner & Stone (2005) test: a weak magnetic field loop (radius rloop=0.3,
|B| = Azero = 1e-3 inside, 0 outside) is advected diagonally across a periodic
box at v = sqrt(5)*(cos t, sin t, 0.1). The box (x in [-1,1], y in [-0.5,0.5])
and velocity are chosen so the loop returns to its start after each unit time;
with nsteps=200, deltastep=0.1 the run ends at t=20 (20 crossings). An ideal
scheme advects the loop unchanged, so:

  * the magnetic energy E_mag = sum 1/2 |B|^2 V_i is conserved -- its decay is a
    direct measure of numerical dissipation,
  * |B| keeps its top-hat-in-radius loop shape; smearing / a drop below |B|_0 is
    the dissipation in real space,
  * div(B) stays ~0; |div B| h / |B| is the divergence-cleanliness diagnostic.

This script produces:
  * E_mag(t) / E_mag(0) vs time, against the ideal conserved line 1 (always),
  * --render: |B|/|B|_0 at --time (default 20) and the normalised divergence
    error |div B| h / |B|, both in the x-y advection plane.

Usage:
    python IC_analysis_mhdloop.py <run-dir> [-a 1.0] [--time 20] [--save out]
    python IC_analysis_mhdloop.py <run-dir> --render --save out
"""

import os
import sys

import numpy as np
import readtipsy as tip

from IC_analysis_framework import (
    Param, RenderSpec, render_panels, aux_divB_error, standalone_cli,
    series_from_log_or_snapshots, spec_slug,
    bless_render_specs, render_ref_overlay_path, movie_render_specs,
)
from IC_analysis_general import (
    read_vector_aux,
    setup_rcparams, finish_figure_with_legend,
)

GAMMA = 5.0 / 3.0              # IC_setup_mhdloop.gamma
AZERO = 1.0e-3                 # initial |B| inside the loop (setup hardcodes this)
RLOOP = 0.3                    # loop radius
DEFAULT_TIME = 20.0           # nsteps*deltastep = 200*0.1 (loop back at origin)

# runtest.sh setup_mhdloop() parameter order. Only rhodiff (the loop/background
# density contrast) is physical here.
PARAM_SPEC = [
    Param("a", "rhodiff", 1.0, float,
          "loop/background density contrast rhodisk/rhozero"),
]


# --------------------------------------------------------------------------- #
#  Render quantities
# --------------------------------------------------------------------------- #

def q_bmag_norm(fn, tgdata, ngas):
    """|B| / |B|_0 per particle (|B|_0 = Azero), from the BField aux.

    1 inside the freshly-advected loop, 0 outside; a peak below 1 (and a smeared
    edge) is numerical dissipation. Zeros where the aux is absent."""
    B = read_vector_aux(fn, "BField", ngas)
    if B is None:
        return np.zeros(len(tgdata))
    bmag = np.sqrt(np.sum(np.asarray(B, dtype=float) ** 2, axis=1))
    return bmag / AZERO


def _clim_bnorm(params, proj):
    """Fixed [0,1] color range for the field-valued |B|/|B|_0 projections so the
    loop's decay is read off directly; None for the integrated `column`."""
    return None if proj == "column" else (0.0, 1.0)


def render_specs():
    """The two render fields: |B|/|B|_0 (linear, [0,1]) and the normalised
    divergence error |div B| h / |B| (log)."""
    return [
        RenderSpec(1, 2, q_bmag_norm, r"$|B|/|B|_0$", False, "average",
                   _clim_bnorm, "inferno", False, slug="Bnorm"),
        RenderSpec(1, 2, aux_divB_error(), r"$|\nabla\!\cdot\!B|\,h/|B|$", True,
                   "average", None, "inferno", False, slug="divBerr"),
    ]


def do_render(inputs, labels, t, args, save, reference=None):
    """Render |B|/|B|_0 and divBerr at time t, forwarding the --render-* flags.

    Each spec goes to its own figure (file suffixed by the field). The CLI
    overrides win: --render-backend, --render-project (None -> the spec's
    'average' default), --render-slice-frac, --render-res, --render-extent
    (default: the box auto-read from the run .log). When `reference` is given and
    there is a single input, the blessed render grids are overlaid (sim | ref |
    difference)."""
    entries = list(zip(labels, inputs))
    overlay = reference is not None and len(inputs) == 1
    for sp in render_specs():
        slug = spec_slug(sp)
        sp_save = f"{save}_{slug}" if save else None
        ref_path = (render_ref_overlay_path(inputs[0], reference, slug)
                    if overlay else None)
        render_panels(entries, sp, times=[t], n=1, res=args.render_res,
                      save=sp_save, backend=args.render_backend,
                      project=args.render_project,
                      slice_frac=args.render_slice_frac,
                      extent=args.render_extent, aspect=args.render_aspect,
                      render_reference=ref_path,
                      residual=getattr(args, "residual", False))


def _save_extra(args, params):
    """Bless the render-grid baselines (at the render time) for --save-reference."""
    bless_render_specs(args.inputs[0], render_specs(), args, "mhdloop",
                       times=[args.time])


# --------------------------------------------------------------------------- #
#  Magnetic energy vs time
# --------------------------------------------------------------------------- #

def emag(fn):
    """(time, E_mag) for snapshot fn, E_mag = sum 1/2 |B|^2 V_i (V_i = m_i/rho_i).

    NaN if the BField aux is missing."""
    tgdata, _td, _ts, hdr, time = tip.readtipsy(fn)
    ngas = int(hdr[2])
    t = float(np.asarray(time).ravel()[0])
    B = read_vector_aux(fn, "BField", ngas)
    if B is None:
        return t, np.nan
    vol = tgdata[:, 0] / tgdata[:, 7]            # m / rho
    b2 = np.sum(np.asarray(B, dtype=float) ** 2, axis=1)
    return t, float(np.sum(0.5 * b2 * vol))


def emag_series(inp):
    """(t, E_mag, source) for run `inp`: the dense '.log' Emag column, else
    per-snapshot integration of the BField aux. See series_from_log_or_snapshots."""
    return series_from_log_or_snapshots(inp, ["dTime", "Emag"], emag)


def plot_emag_vs_time(inputs, labels, save=None, ref_table=None):
    """E_mag(t)/E_mag(0) for each run, against the ideal conserved line (=1).

    Uses the dense '.log' Emag series when available (see emag_series). When a
    reference is supplied its blessed E_mag(t)/E_mag(0) is overlaid (dashed)."""
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
            print(f"mhdloop: no snapshots/.log for '{inp}'", file=sys.stderr)
            continue
        if not np.isfinite(em[0]) or em[0] == 0:
            print(f"mhdloop[{lab}]: no/zero initial E_mag (missing BField aux?); "
                  f"skipping.", file=sys.stderr)
            continue
        style = "-" if source == "log" else "o-"   # dense log -> line
        ax.plot(ts, em / em[0], style, ms=3, label=lab)
        print(f"mhdloop[{lab}]: E_mag(0)={em[0]:.4g}  "
              f"E_mag({ts[-1]:.3g})/E_mag(0)={em[-1] / em[0]:.4g}  "
              f"({len(ts)} pts from {source})")
        any_data = True
    ax.axhline(1.0, color="k", ls="--", lw=1.3, label="ideal (conserved)")
    ax.set_xlabel("t")
    ax.set_ylabel(r"$E_{mag}(t)/E_{mag}(0)$")
    ax.set_title("MHD loop advection: magnetic energy vs time")
    fig.tight_layout()
    if not any_data:
        print("mhdloop: no magnetic energy series to plot.", file=sys.stderr)
    finish_figure_with_legend(fig, ax,
                              save=(f"{save}_emag.png" if save else None))


# --------------------------------------------------------------------------- #
#  Regression baseline (bless / compare the E_mag(t)/E_mag(0) series)
# --------------------------------------------------------------------------- #

def reference_table(inp):
    """Regression table {x: t, Emag: E_mag(t)/E_mag(0)} for run `inp`.

    The magnetic-energy decay normalised to its initial value is the test's
    headline diagnostic and is resolution-independent, so it makes the natural
    CI baseline. None if there is no usable (finite, non-zero) E_mag series."""
    ts, em, _src = emag_series(inp)
    if ts is None or not np.isfinite(em[0]) or em[0] == 0:
        return None
    return {"x": ts, "Emag": em / em[0]}


def _plot(inputs, labels, params, args):
    print(f"[mhdloop] rhodiff={params['rhodiff']}  |B|_0={AZERO:g}  "
          f"rloop={RLOOP:g}")
    plot_emag_vs_time(inputs, labels, save=args.save,
                      ref_table=getattr(args, "ref_table", None))
    if args.render:
        do_render(inputs, labels, args.time, args, save=args.save,
                  reference=args.reference)
    if args.movie:
        movie_render_specs(inputs[0], labels[0], render_specs(), args, args.save,
                           "mhdloop")


def main():
    standalone_cli(
        "mhdloop", PARAM_SPEC, _plot,
        reference_table=lambda inp, params, args: reference_table(inp),
        description="MHD field-loop advection analysis.",
        reg_keys=("Emag",), time_default=DEFAULT_TIME,
        time_help=f"render time (default {DEFAULT_TIME:g}; the loop is back at "
                  f"its starting position)",
        save_extra=_save_extra)


if __name__ == "__main__":
    main()
