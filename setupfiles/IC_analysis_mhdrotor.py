#!/usr/bin/env python3
"""
Analysis for the MHD rotor test (IC_setup_mhdrotor; Balsara & Spicer 1999, Toth
2000).

A dense rotating disk (rho_disk = rhodisk, default 10; radius 0.1; solid-body
rotation v = 2) sits in ambient gas (rho = 1, P = 1) threaded by a uniform field
Bx = 5/sqrt(4pi). The rotor winds the field into a tightly wrapped spiral,
launching torsional Alfven waves that magnetically brake the disk. There is no
closed-form solution, so the diagnostics are the field structure and the energy
exchange (kinetic -> magnetic), exactly as for the Orszag-Tang vortex:

  * E_mag(t)/E_mag0 and E_kin(t)/E_kin0 vs time -- the field winds up (E_mag
    rises) while the rotor spins down (E_kin falls): the magnetic-braking curve.
    Normalised by the ANALYTIC IC energies (E_mag0 = 1/2 B0^2 V, E_kin0 =
    (pi/4) rhodisk v0^2 r0^2 dz; box from the '.log'), NOT each run's first
    output -- so a run whose first dump is after t=0 (e.g. SPH-EXA) shares the
    same baseline as one that dumps the IC,
  * field-map renders (--render) of density rho, magnetic pressure
    Pmag = 1/2 |B|^2 (with 30 Pmag iso-contours over the fixed range [0, 2.642],
    shared across the panels), the signed field components Bx and By, and the
    divergence error h|div B|/|B| (fixed log range 1e-3..1.0), in
    HORIZONTAL SLICES through the centre (x-y plane at z=0) at chosen times
    (default t = 0.05, 0.15, 0.25, 0.5; the classic comparison time is 0.15),
  * magnetic-field STREAMLINES (Bx, By) in the x-y plane at the same times
    (with --render, like currentsheet's `_bstream`) -- the initially uniform
    horizontal field wound into the spiral around the rotor (initial disk
    radius marked).

The classic Balsara gas-pressure / Mach panels are omitted because they need the
tipsy temperature->u factor (dTuFac); every field rendered here is dTuFac-free.

Standalone (reuses the framework loaders + render_panels, not run_cli).

Usage:
    python IC_analysis_mhdrotor.py <run-dir> [-a rhodisk] [--save out]
    python IC_analysis_mhdrotor.py <run-dir> --render --render-times 0.15 --save out
"""

import os
import sys

import numpy as np
import readtipsy as tip

from IC_analysis_framework import (
    Param, standalone_cli,
    styled, spec_slug, render_panels, ContourSpec,
    StreamSpec, streamline_panels,
    aux_divB_error, aux_component, magnetic_pressure,
    horizontal_cut, mhd_field_cut_quants, series_from_log_or_snapshots,
    bless_render_specs, render_ref_overlay_path, movie_render_specs,
)
from IC_analysis_general import (
    read_vector_aux, read_log_param,
    setup_rcparams, finish_figure_with_legend,
)

# Constants fixed by IC_setup_mhdrotor.create() (not -a..-c parameters).
GAMMA = 1.4
RDISK = 0.1                    # disk radius
VZERO = 2.0                    # solid-body rotation speed (|v| max in the disk)
RHOZERO = 1.0                  # ambient density
BZERO = 5.0 / np.sqrt(4.0 * np.pi)   # uniform initial Bx
DEFAULT_RENDER_TIMES = [0.05, 0.15, 0.25, 0.5]   # 0.15 = the standard rotor time

# runtest.sh setup_mhdrotor(): -a rhodisk (disk density), -b/-c viscosities.
PARAM_SPEC = [
    Param("a", "rhodisk", 10.0, float, "rotating-disk density (ambient = 1)"),
]


# --------------------------------------------------------------------------- #
#  Field-map renders
# --------------------------------------------------------------------------- #

def _fixed_clim(lo, hi):
    """clim callback: a fixed range for field-valued projections, None (auto) for
    an integrated `column` projection (different units)."""
    def clim(params, proj):
        return None if proj == "column" else (lo, hi)
    return clim


# divBerr: fixed LOG range 1e-3 .. 1.0 (3 decades).
DIVBERR_CLIM = _fixed_clim(1e-3, 1.0)


def render_specs(params):
    """The five x-y field maps -- density rho, magnetic pressure Pmag, the signed
    field components Bx, By, and the divergence error |div B| h/|B| -- in
    HORIZONTAL SLICES through the centre (project="slice" at z=0; the box is a
    thin z-slab so this is the mid-plane). Density spans ambient..disk; Bx/By
    use a symmetric diverging scale; divBerr is fixed log 1e-3..1.0. The Pmag
    map carries 30 fixed iso-contours of Pmag itself over [0, 2.642] (the
    Toth 2000 rotor figure convention, shared by all panels) -- the wound-up
    spiral and the torsional Alfven front."""
    rhodisk = float(params.get("rhodisk", 10.0))

    def rho_clim(p, proj):
        return (RHOZERO, rhodisk) if proj in ("slice", "average", "rhocolumn") else None

    return [
        styled("rho", 1, 2, 7, label=r"$\rho$", clim=rho_clim),
        styled("Pmag", 1, 2, magnetic_pressure(), cmap="rainbow", log=False,
               clim=_fixed_clim(0.0, 2.642),
               label=r"$P_{\rm mag}=\frac{1}{2}|B|^2$",
               contour=ContourSpec(magnetic_pressure(), label=r"$P_{\rm mag}$",
                                   levels=np.linspace(0.0, 2.642, 30))),
        styled("Bx", 1, 2, aux_component("BField", 0), label=r"$B_x$"),
        styled("By", 1, 2, aux_component("BField", 1), label=r"$B_y$"),
        styled("divBerr", 1, 2, aux_divB_error(),
               label=r"$|\nabla\!\cdot\!B|\,h/|B|$", clim=DIVBERR_CLIM),
    ]


def _render_times(args):
    return args.render_times if args.render_times else DEFAULT_RENDER_TIMES


def do_render(inputs, labels, params, args, save, reference=None):
    """Render every field map at the requested times, forwarding all --render-*.

    With --residual the panels become a [A | B | A-B] comparison grid: B is the
    blessed render baseline (single input + --reference) or the second input
    (two inputs, differenced directly)."""
    entries = list(zip(labels, inputs))
    times = _render_times(args)
    overlay = reference is not None and len(inputs) == 1
    for sp in render_specs(params):
        slug = spec_slug(sp)
        sp_save = f"{save}_{slug}" if save else None
        ref_path = (render_ref_overlay_path(inputs[0], reference, slug)
                    if overlay else None)
        render_panels(entries, sp, times=times, n=args.render_n,
                      res=args.render_res, save=sp_save,
                      backend=args.render_backend, project=args.render_project,
                      slice_frac=args.render_slice_frac,
                      extent=args.render_extent, aspect=args.render_aspect,
                      render_reference=ref_path,
                      residual=getattr(args, "residual", False))


def plot_bfield_streamlines(inputs, labels, args, save):
    """Magnetic-field streamlines (Bx, By) in the x-y plane at the render times
    -- directly shows the initially uniform horizontal field wound into the
    spiral around the rotor (drawn with --render, like currentsheet's).

    Delegates to the framework `streamline_panels`; the per-test bits are the
    StreamSpec (Bx, By from the BField aux) and the initial disk-radius marker."""
    spec = StreamSpec(1, 2, aux_component("BField", 0),
                      aux_component("BField", 1), "B")

    def disk(ax):                                     # initial rotor radius
        from matplotlib.patches import Circle
        ax.add_patch(Circle((0.0, 0.0), RDISK, fill=False, color="0.6",
                            ls=":", lw=0.8, zorder=1))

    streamline_panels(list(zip(labels, inputs)), spec,
                      times=_render_times(args), n=args.render_n,
                      res=min(args.render_res or 64, 128),
                      save=(f"{save}_bstream" if save else None),
                      extent=args.render_extent, decorate=disk,
                      title="MHD rotor: magnetic-field streamlines")


def _save_extra(args, params):
    """Bless the render-grid baselines (at the render times) for --save-reference."""
    bless_render_specs(args.inputs[0], render_specs(params), args, "mhdrotor",
                       times=_render_times(args))


def do_cut(inputs, labels, args, save):
    """1-D horizontal cuts of rho, Bx, By, div B vs x along y=y0 (default centre;
    --y0 overrides), at each render time. Delegates to the framework's generic
    `horizontal_cut`."""
    times = args.render_times if args.render_times else DEFAULT_RENDER_TIMES
    horizontal_cut(inputs, labels, mhd_field_cut_quants(), times,
                   lambda t: (args.y0 if args.y0 is not None else 0.0), save=save,
                   title_prefix="MHD rotor ")


# --------------------------------------------------------------------------- #
#  Energy exchange (magnetic braking) vs time
# --------------------------------------------------------------------------- #

def _energy_snapshot(fn):
    """(time, E_mag, E_kin) for snapshot fn. E_mag NaN if the BField aux is absent."""
    tgdata, _td, _ts, hdr, time = tip.readtipsy(fn)
    ngas = int(hdr[2])
    t = float(np.asarray(time).ravel()[0])
    m = tgdata[:, 0]
    vol = m / tgdata[:, 7]
    e_kin = float(np.sum(0.5 * m * np.sum(tgdata[:, 4:7] ** 2, axis=1)))
    B = read_vector_aux(fn, "BField", ngas)
    if B is None:
        return t, np.nan, e_kin
    b2 = np.sum(np.asarray(B, float) ** 2, axis=1)
    return t, float(np.sum(0.5 * b2 * vol)), e_kin


def energy_series(inp):
    """(t, E_mag, E_kin, source) for run `inp`.

    Prefers the dense gasoline '.log' Emag/Ekin columns; falls back to
    per-snapshot integration. Returns (None, ...) if neither yields data."""
    return series_from_log_or_snapshots(inp, ["dTime", "Emag", "Ekin"],
                                        _energy_snapshot)


def analytic_E0(inp, rhodisk):
    """Exact IC energies (E_mag0, E_kin0) of the sharp-edged rotor.

    E_mag0 = 1/2 B0^2 V (uniform Bx over the box) and E_kin0 = (pi/4) rhodisk
    v0^2 r0^2 dz (solid-body disk, ambient at rest), with the box dimensions
    read from the run '.log'. Used to anchor the E(t)/E0 normalisation at the
    TRUE t=0 even when the run has no t=0 output (e.g. SPH-EXA's first dump at
    t=0.05, which would otherwise become the baseline and hide the early
    growth/decay). Returns (None, None) when the '.log' has no box."""
    dz = read_log_param(inp, "dzPeriod")
    if dz is None or dz <= 0:
        return None, None
    dx = read_log_param(inp, "dxPeriod", 1.0)
    dy = read_log_param(inp, "dyPeriod", 1.0)
    e_mag0 = 0.5 * BZERO ** 2 * (dx * dy * dz)
    e_kin0 = 0.25 * np.pi * rhodisk * VZERO ** 2 * RDISK ** 2 * dz
    return e_mag0, e_kin0


def _normalised_energy_series(inp, lab, rhodisk):
    """(t, E_mag/E_mag0, E_kin/E_kin0, source) with the ANALYTIC IC energies as
    the normalisation baseline (falling back to each series' first point when
    the '.log' carries no box, with a stderr note). None keys are skipped by
    the callers via NaN entries."""
    ts, em, ek, source = energy_series(inp)
    if ts is None:
        return None, None, None, None
    em0, ek0 = analytic_E0(inp, rhodisk)
    if em0 is None:
        print(f"mhdrotor[{lab}]: no box in '.log' -> normalising by the first "
              f"point (t={ts[0]:.4g}), not the analytic E0", file=sys.stderr)
        em0 = em[0] if np.isfinite(em[0]) and em[0] > 0 else None
        ek0 = ek[0] if np.isfinite(ek[0]) and ek[0] > 0 else None
    em = em / em0 if em0 else np.full_like(np.asarray(em, float), np.nan)
    ek = ek / ek0 if ek0 else np.full_like(np.asarray(ek, float), np.nan)
    return ts, em, ek, source


def plot_energy_vs_time(inputs, labels, params, save=None, ref_table=None):
    """E_mag(t)/E_mag0 and E_kin(t)/E_kin0 for each run (magnetic braking:
    E_mag rises as the field winds up, E_kin falls as the rotor spins down),
    normalised by the ANALYTIC IC energies (analytic_E0) so runs whose output
    starts after t=0 share the same baseline as runs that dump the IC.

    When a reference is supplied its blessed E_mag/E_kin braking curves are
    overlaid as dashed grey lines."""
    import matplotlib.pyplot as plt
    setup_rcparams()
    rhodisk = float(params.get("rhodisk", 10.0))
    fig, ax = plt.subplots(figsize=(8, 6))
    any_data = False
    if ref_table is not None:
        for key in ("Emag", "Ekin"):
            if key in ref_table:
                ax.plot(ref_table["x"], ref_table[key], "--", color="0.35",
                        lw=1.4, zorder=5,
                        label=("reference" if key == "Emag" else None))
    for inp, lab in zip(inputs, labels):
        ts, em, ek, source = _normalised_energy_series(inp, lab, rhodisk)
        if ts is None:
            print(f"mhdrotor: no snapshots/.log for '{inp}'", file=sys.stderr)
            continue
        style = "-" if source == "log" else "o-"
        if np.any(np.isfinite(em)):
            ax.plot(ts, em, style, ms=3, label=f"{lab}  $E_{{mag}}$")
            print(f"mhdrotor[{lab}]: E_mag/E_mag0 max={np.nanmax(em):.4g}  "
                  f"end={em[-1]:.4g}  ({len(ts)} pts from {source})")
        else:
            print(f"mhdrotor[{lab}]: no/zero E_mag (missing BField aux?).",
                  file=sys.stderr)
        if np.any(np.isfinite(ek)):
            ax.plot(ts, ek, style.replace("-", "--", 1), ms=3,
                    label=f"{lab}  $E_{{kin}}$")
            print(f"mhdrotor[{lab}]: E_kin/E_kin0 end={ek[-1]:.4g}")
        any_data = True
    ax.axhline(1.0, color="0.7", ls=":", lw=1.0)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$E(t)/E_0^{\rm analytic}$")
    ax.set_title("MHD rotor: magnetic braking (solid $E_{mag}$, dashed $E_{kin}$)")
    fig.tight_layout()
    if not any_data:
        print("mhdrotor: no energy series to plot.", file=sys.stderr)
    finish_figure_with_legend(fig, ax,
                              save=(f"{save}_energy.png" if save else None))


# --------------------------------------------------------------------------- #
#  Regression baseline (bless / compare the normalised energy series)
# --------------------------------------------------------------------------- #

def reference_table(inp, params):
    """Regression table {x: t, Emag: E_mag/E_mag0, Ekin: E_kin/E_kin0} for `inp`,
    normalised by the SAME analytic IC energies as the plot.

    The normalised magnetic-braking curves (E_mag rises, E_kin falls) are the
    rotor's headline diagnostic and resolution-independent -> the CI baseline.
    Each key is included only when its series has finite values; None if
    neither is usable."""
    rhodisk = float(params.get("rhodisk", 10.0))
    ts, em, ek, _src = _normalised_energy_series(inp, str(inp), rhodisk)
    if ts is None:
        return None
    table = {"x": ts}
    if np.any(np.isfinite(em)):
        table["Emag"] = em
    if np.any(np.isfinite(ek)):
        table["Ekin"] = ek
    return table if len(table) > 1 else None


def _plot(inputs, labels, params, args):
    print("[mhdrotor] Balsara-Spicer MHD rotor")
    plot_energy_vs_time(inputs, labels, params, save=args.save,
                        ref_table=getattr(args, "ref_table", None))
    if args.render:
        do_render(inputs, labels, params, args, save=args.save,
                  reference=args.reference)
        plot_bfield_streamlines(inputs, labels, args, save=args.save)
        do_cut(inputs, labels, args, save=args.save)
    if args.movie:
        movie_render_specs(inputs[0], labels[0], render_specs(params), args,
                           args.save, "mhdrotor")


def main():
    standalone_cli(
        "mhdrotor", PARAM_SPEC, _plot,
        reference_table=lambda inp, params, args: reference_table(inp, params),
        description="MHD rotor analysis.", reg_keys=("Emag", "Ekin"),
        add_time=False, save_extra=_save_extra)


if __name__ == "__main__":
    main()
