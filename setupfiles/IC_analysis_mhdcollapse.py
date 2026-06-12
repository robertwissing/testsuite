#!/usr/bin/env python3
"""
Analysis for the magnetized rotating cloud collapse (IC_setup_mhdcollapse).

This IS the Wissing & Shen (2020, A&A 638, A140; SPMHD in Gasoline2/ChaNGa)
magnetized-cloud-collapse test: a 1 Msol cloud, radius R_c = 0.015 pc, cloud/
medium density contrast 360, solid-body rotation (period 4.7e5 yr, E_K/E_P =
0.045), a uniform field B || z (the rotation axis) set by the mass-to-flux ratio
mu, a barotropic EOS (isothermal at c_s = 0.2 km/s, stiffening above
rho_0 = 1e-14 g/cm^3), in a 0.15 pc box. Under self-gravity the cloud collapses,
flattens into a pseudo-disk, and -- in the weak-field regime -- winds the field
into a "magnetic tower" that launches a bipolar jet.

Following Wissing & Shen's collapse analysis, the diagnostics are:

  * field-map renders (--render) zoomed to the central ~2000 AU (their window),
    as mid-plane SLICES (the cloud is a 3-D sphere):
      - density face-on (x-y) and edge-on (x-z),
      - poloidal |B_pol| = sqrt(B_R^2 + B_z^2) edge-on  (the magnetic tower),
      - toroidal B_phi = (x B_y - y B_x)/R edge-on       (the wound-up field),
      - spherical radial velocity v_r edge-on            (infall / the jet),
      - divergence error h|div B|/|B| edge-on,
  * the B-rho field-amplification relation (binned median |B| vs rho at the final
    snapshot) against the flux-freezing references rho^(1/2) and rho^(2/3),
  * peak density rho_peak(t)/rho_peak(0) vs time (run to first-core densities),
  * the energy budget (E_kin, E_mag, |E_pot|) vs time (barotropic EOS -> no E_th).

Standalone (reuses the framework loaders + render_panels, not run_cli). All
diagnostics are dimensionless ratios / code units, avoiding the IC's units; the
2000 AU render window assumes the setup's length unit is 1 pc (dkpcunit=0.001).

Usage:
    python IC_analysis_mhdcollapse.py <run-dir> [-a mu -b rhodiff -c Erat] [--save out]
    python IC_analysis_mhdcollapse.py <run-dir> --render --render-n 4 --save out
"""

import os
import sys

import numpy as np
import readtipsy as tip

from IC_analysis_framework import (
    Param, standalone_cli,
    RenderSpec, spec_slug, render_panels, aux_divB_error, code_units,
    bless_render_specs, render_ref_overlay_path, movie_render_specs,
)
from IC_analysis_general import (
    find_files, read_vector_aux, read_log_series,
    setup_rcparams, finish_figure_with_legend,
)

# runtest.sh setup_mhdcollapse(): the three IC physics parameters. Carried for
# labelling/context; the (reference) diagnostics here are parameter-agnostic.
PARAM_SPEC = [
    Param("a", "mu",      10.0,  float, "mass-to-flux ratio (sets Bz; mu=critical units)"),
    Param("b", "rhodiff", 360.0, float, "cloud/medium density contrast"),
    Param("c", "Erat",    0.045, float, "rotational energy ratio (sets the spin)"),
]

# Central render window: Wissing & Shen use L = 2000 AU (face-on and edge-on).
# The setup's length unit is 1 pc (dkpcunit = 0.001 kpc), so 1000 AU (the half
# window) = 1000 / 206264.806 pc in code units. Overridable via --render-extent.
_AU_PER_PC = 206264.806
RENDER_HALF = 1000.0 / _AU_PER_PC          # half of the 2000 AU window, in pc/code units


# --------------------------------------------------------------------------- #
#  Cylindrical / spherical field decompositions (callable render quantities)
# --------------------------------------------------------------------------- #

def _Rcyl(tgdata):
    return np.hypot(tgdata[:, 1], tgdata[:, 2])


def _bpol(fn, tgdata, ngas):
    """Poloidal field magnitude |B_pol| = sqrt(B_R^2 + B_z^2) (cylindrical R about
    the z rotation axis) -- the magnetic-tower diagnostic. Zeros if no BField aux."""
    B = read_vector_aux(fn, "BField", ngas)
    if B is None:
        return np.zeros(len(tgdata))
    B = np.asarray(B, float)
    R = _Rcyl(tgdata)
    with np.errstate(invalid="ignore", divide="ignore"):
        B_R = np.where(R > 0, (tgdata[:, 1] * B[:, 0] + tgdata[:, 2] * B[:, 1]) / R, 0.0)
    return np.sqrt(B_R ** 2 + B[:, 2] ** 2)


def _btor(fn, tgdata, ngas):
    """Toroidal (azimuthal) field B_phi = (x B_y - y B_x)/R -- the wound-up field
    (signed). Zeros if no BField aux."""
    B = read_vector_aux(fn, "BField", ngas)
    if B is None:
        return np.zeros(len(tgdata))
    B = np.asarray(B, float)
    R = _Rcyl(tgdata)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(R > 0, (tgdata[:, 1] * B[:, 1] - tgdata[:, 2] * B[:, 0]) / R, 0.0)


def _vr(fn, tgdata, ngas):
    """Spherical radial velocity v_r = (r . v)/|r| (signed: <0 infall, >0 outflow/jet)."""
    r = np.sqrt(np.sum(tgdata[:, 1:4] ** 2, axis=1))
    vr = np.sum(tgdata[:, 1:4] * tgdata[:, 4:7], axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(r > 0, vr / r, 0.0)


def _bmag(fn, tgdata, ngas):
    """Total |B| per particle (for the B-rho relation). NaN array if no aux."""
    B = read_vector_aux(fn, "BField", ngas)
    if B is None:
        return np.full(len(tgdata), np.nan)
    return np.sqrt(np.sum(np.asarray(B, float) ** 2, axis=1))


# --------------------------------------------------------------------------- #
#  Peak (central) density vs time -- the collapse history
# --------------------------------------------------------------------------- #

def density_series(inp):
    """(t, rho_peak) arrays for run `inp`, rho_peak = 99.9th percentile of rho."""
    try:
        files = find_files(inp)
    except FileNotFoundError:
        return None, None
    if not files:
        return None, None
    ts, peaks = [], []
    for fn in files:
        tgdata, _td, _ts, _hdr, time = tip.readtipsy(fn)
        ts.append(float(np.asarray(time).ravel()[0]))
        peaks.append(float(np.percentile(tgdata[:, 7], 99.9)))
    o = np.argsort(ts)
    return np.array(ts)[o], np.array(peaks)[o]


def plot_density_vs_time(inputs, labels, save=None, ref_table=None):
    """rho_peak(t)/rho_peak(0) for each run (the collapse amplification).

    When a reference is supplied its blessed rho_peak(t)/rho_peak(0) is overlaid
    (dashed grey)."""
    import matplotlib.pyplot as plt
    setup_rcparams()
    fig, ax = plt.subplots(figsize=(8, 6))
    any_data = False
    if ref_table is not None and "rho_peak" in ref_table:
        ax.semilogy(ref_table["x"], ref_table["rho_peak"], "--", color="0.35",
                    lw=1.5, zorder=5, label="reference")
    for inp, lab in zip(inputs, labels):
        ts, peaks = density_series(inp)
        if ts is None:
            print(f"mhdcollapse: no snapshots for '{inp}'", file=sys.stderr)
            continue
        if peaks[0] <= 0:
            print(f"mhdcollapse[{lab}]: zero initial peak density; skipping.",
                  file=sys.stderr)
            continue
        ax.semilogy(ts, peaks / peaks[0], "o-", ms=3, label=lab)
        print(f"mhdcollapse[{lab}]: rho_peak amplification end="
              f"{peaks[-1]/peaks[0]:.4g} over {len(ts)} snapshots")
        any_data = True
    ax.set_xlabel("t")
    ax.set_ylabel(r"$\rho_{peak}(t)/\rho_{peak}(0)$")
    ax.set_title("MHD collapse: central density amplification vs time")
    fig.tight_layout()
    if not any_data:
        print("mhdcollapse: no density series to plot.", file=sys.stderr)
    finish_figure_with_legend(fig, ax,
                              save=(f"{save}_density.png" if save else None))


# --------------------------------------------------------------------------- #
#  B-rho field-amplification relation (Wissing & Shen / Price & Bate)
# --------------------------------------------------------------------------- #

def plot_brho(inputs, labels, save=None):
    """Binned median |B| vs rho at the final snapshot, against flux-freezing
    references |B| ~ rho^(1/2) (1-D / strong field) and rho^(2/3) (spherical).

    Whether the simulated amplification tracks rho^(2/3) (ideal spherical flux
    freezing) or falls below it (numerical/physical resistivity, or disk geometry)
    is the headline field-amplification diagnostic for the collapse."""
    import matplotlib.pyplot as plt
    setup_rcparams()
    fig, ax = plt.subplots(figsize=(8, 6))
    any_data = False
    anchor = None
    for inp, lab in zip(inputs, labels):
        try:
            files = find_files(inp)
        except FileNotFoundError:
            files = []
        if not files:
            print(f"mhdcollapse: no snapshots for '{inp}'", file=sys.stderr)
            continue
        fn = files[-1]
        tgdata, _td, _ts, hdr, _time = tip.readtipsy(fn)
        rho = tgdata[:, 7]
        Bmag = _bmag(fn, tgdata, int(hdr[2]))
        good = np.isfinite(Bmag) & (Bmag > 0) & (rho > 0)
        if not np.any(good):
            print(f"mhdcollapse[{lab}]: no BField aux at the final snapshot; "
                  f"skipping B-rho.", file=sys.stderr)
            continue
        rho, Bmag = rho[good], Bmag[good]
        # Median |B| in log-rho bins.
        edges = np.logspace(np.log10(rho.min()), np.log10(rho.max()), 25)
        idx = np.digitize(rho, edges)
        rc, bc = [], []
        for b in range(1, len(edges)):
            m = idx == b
            if np.count_nonzero(m) >= 5:
                rc.append(np.sqrt(edges[b - 1] * edges[b]))
                bc.append(np.median(Bmag[m]))
        if not rc:
            continue
        rc, bc = np.array(rc), np.array(bc)
        ax.loglog(rc, bc, "o-", ms=3, label=lab)
        anchor = anchor or (rc[0], bc[0])
        print(f"mhdcollapse[{lab}]: B-rho at final snapshot, "
              f"rho in [{rc[0]:.3g},{rc[-1]:.3g}], |B| in [{bc.min():.3g},{bc.max():.3g}]")
        any_data = True
    if anchor is not None:
        r0, b0 = anchor
        rg = np.logspace(np.log10(r0), np.log10(r0 * 1e6), 50)
        ax.loglog(rg, b0 * (rg / r0) ** (2.0 / 3.0), "k--", lw=1.2,
                  label=r"$\propto\rho^{2/3}$ (spherical)")
        ax.loglog(rg, b0 * (rg / r0) ** (1.0 / 2.0), "k:", lw=1.2,
                  label=r"$\propto\rho^{1/2}$")
    ax.set_xlabel(r"$\rho$")
    ax.set_ylabel(r"$|B|$")
    ax.set_title("MHD collapse: field amplification (final snapshot)")
    fig.tight_layout()
    if not any_data:
        print("mhdcollapse: no B-rho data to plot.", file=sys.stderr)
    finish_figure_with_legend(fig, ax,
                              save=(f"{save}_brho.png" if save else None))


# --------------------------------------------------------------------------- #
#  Energy budget vs time
# --------------------------------------------------------------------------- #

def _energy_snapshot(fn):
    """(t, E_kin, E_mag) for snapshot fn (per-snapshot fallback when no .log)."""
    tgdata, _td, _ts, hdr, time = tip.readtipsy(fn)
    ngas = int(hdr[2])
    t = float(np.asarray(time).ravel()[0])
    m = tgdata[:, 0]
    e_kin = float(np.sum(0.5 * m * np.sum(tgdata[:, 4:7] ** 2, axis=1)))
    B = read_vector_aux(fn, "BField", ngas)
    if B is None:
        return t, e_kin, np.nan
    vol = m / tgdata[:, 7]
    b2 = np.sum(np.asarray(B, float) ** 2, axis=1)
    return t, e_kin, float(np.sum(0.5 * b2 * vol))


def energy_series_norm(inp):
    """(t, {|E_pot|, E_kin, E_mag}) each normalised by |E_pot(0)| from the dense
    '.log' (Ekin/Epot/Emag), or (None, None). The blessed energy-overlay source."""
    log = read_log_series(inp, ["dTime", "Ekin", "Epot", "Emag"])
    if log is None or not len(log["dTime"]):
        return None, None
    escale = abs(log["Epot"][0]) or abs(log["Ekin"][0]) or 1.0
    return log["dTime"], {"|E_pot|": np.abs(log["Epot"]) / escale,
                          "E_kin": np.abs(log["Ekin"]) / escale,
                          "E_mag": np.abs(log["Emag"]) / escale}


def plot_energy_vs_time(inputs, labels, save=None, ref_table=None):
    """Energy components vs time, normalized by |E_pot(0)|. Prefers the dense .log
    (Ekin/Epot/Eth/Emag); falls back to per-snapshot E_kin/E_mag. When a reference
    is supplied its blessed |E_pot|/E_kin/E_mag are overlaid (grey)."""
    import matplotlib.pyplot as plt
    setup_rcparams()
    fig, ax = plt.subplots(figsize=(8, 6))
    # Barotropic EOS -> no independent thermal energy; E_th is omitted.
    comp_style = {"|E_pot|": "-", "E_kin": "--", "E_mag": "-."}
    if ref_table is not None and "t_energy" in ref_table:
        te = ref_table["t_energy"]
        for key, rk in (("|E_pot|", "Epot_n"), ("E_kin", "Ekin_n"),
                        ("E_mag", "Emag_n")):
            if rk in ref_table:
                ax.semilogy(te, ref_table[rk], comp_style[key], color="0.5",
                            lw=1.2, zorder=5,
                            label=("reference" if rk == "Epot_n" else None))
    any_data = False
    for inp, lab in zip(inputs, labels):
        log = read_log_series(inp, ["dTime", "Ekin", "Epot", "Emag"])
        if log is not None and len(log["dTime"]):
            t = log["dTime"]
            escale = abs(log["Epot"][0]) or abs(log["Ekin"][0]) or 1.0
            series = {"|E_pot|": np.abs(log["Epot"]), "E_kin": log["Ekin"],
                      "E_mag": log["Emag"]}
            src = "log"
        else:
            ts, ek, em = [], [], []
            try:
                files = find_files(inp)
            except FileNotFoundError:
                files = []
            for fn in files:
                tt, e_kin, e_mag = _energy_snapshot(fn)
                ts.append(tt); ek.append(e_kin); em.append(e_mag)
            if not ts:
                print(f"mhdcollapse: no snapshots/.log for '{inp}'", file=sys.stderr)
                continue
            o = np.argsort(ts)
            t = np.array(ts)[o]
            escale = abs(np.array(ek)[o][0]) or 1.0
            series = {"E_kin": np.array(ek)[o], "E_mag": np.array(em)[o]}
            src = "snapshots"
        for name, vals in series.items():
            if np.any(np.isfinite(vals)) and np.any(vals != 0):
                ax.semilogy(t, np.abs(vals) / escale, comp_style.get(name, "-"),
                            label=f"{lab} {name}")
        print(f"mhdcollapse[{lab}]: energy budget from {src} ({len(t)} pts)")
        any_data = True
    ax.set_xlabel("t")
    ax.set_ylabel(r"$|E|/|E_{pot}(0)|$")
    ax.set_title("MHD collapse: energy budget vs time")
    fig.tight_layout()
    if not any_data:
        print("mhdcollapse: no energy series to plot.", file=sys.stderr)
    finish_figure_with_legend(fig, ax,
                              save=(f"{save}_energy.png" if save else None))


# --------------------------------------------------------------------------- #
#  Field-map renders (Wissing & Shen collapse figures)
# --------------------------------------------------------------------------- #

def _fixed_clim(lo, hi):
    def clim(params, proj):
        return None if proj == "column" else (lo, hi)
    return clim


def render_specs(units=None):
    """Density face-on (x-y) + edge-on (x-z), poloidal |B| (tower) and toroidal
    |B_phi| (wound field) edge-on, the radial-velocity magnitude |v_r| (infall/
    jet) edge-on, and the divergence error edge-on, all as mid-plane SLICES of the
    3-D cloud. Fields are shown in PHYSICAL units on LOG color scales matched to
    the collapse: rho 1e-17..1e-11 g/cm^3, |B_pol|/|B_tor| 1e2..1e6 uG, |v_r|
    0.1..10 km/s, divBerr 1e-2..1.0. The code->cgs factors come from the generic
    framework `code_units` (None -> render in code units). (B_tor and v_r are shown
    as magnitudes so they share the log scale; the signed/diverging versions can be
    restored if the winding/infall sign is wanted.)"""
    have = units is not None
    rf = units.density_gcc if have else 1.0
    vf = units.velocity_kms if have else 1.0
    bf = units.b_uG if have else 1.0
    rlab = r"$\rho$ [g cm$^{-3}$]" if have else r"$\rho$ [code]"
    vlab = r"$|v_r|$ [km s$^{-1}$]" if have else r"$|v_r|$ [code]"
    blab_pol = r"$|B_{pol}|$ [$\mu$G]" if have else r"$|B_{pol}|$ [code]"
    blab_tor = r"$|B_{tor}|$ [$\mu$G]" if have else r"$|B_{tor}|$ [code]"

    def rho_q(fn, tg, ng):
        return tg[:, 7] * rf

    def bpol_q(fn, tg, ng):
        return _bpol(fn, tg, ng) * bf

    def btor_q(fn, tg, ng):
        return np.abs(_btor(fn, tg, ng)) * bf

    def vr_q(fn, tg, ng):
        return np.abs(_vr(fn, tg, ng)) * vf

    rho_clim = _fixed_clim(1e-17, 1e-11) if have else None
    v_clim = _fixed_clim(0.1, 10.0) if have else None
    b_clim = _fixed_clim(1e2, 1e6) if have else None
    db_clim = _fixed_clim(1e-2, 1.0)
    # NOTE: B_tor and v_r are rendered as MAGNITUDES on a shared log scale
    # (deliberate -- see the docstring), so they keep inferno/viridis here rather
    # than the canonical signed RdBu_r style. Each spec carries an explicit slug.
    return [
        RenderSpec(1, 2, rho_q, rlab + " [face-on]", True, "slice", rho_clim,
                   "inferno", False, slug="rho_faceon"),
        RenderSpec(1, 3, rho_q, rlab + " [edge-on]", True, "slice", rho_clim,
                   "inferno", False, slug="rho_edgeon"),
        RenderSpec(1, 3, bpol_q, blab_pol + " [edge-on]", True, "slice",
                   b_clim, "inferno", False, slug="bpol"),
        RenderSpec(1, 3, btor_q, blab_tor + " [edge-on]", True, "slice",
                   b_clim, "inferno", False, slug="btor"),
        RenderSpec(1, 3, vr_q, vlab + " [edge-on]", True, "slice", v_clim,
                   "viridis", False, slug="vr"),
        RenderSpec(1, 3, aux_divB_error(), r"$|\nabla\!\cdot\!B|\,h/|B|$ [edge-on]",
                   True, "slice", db_clim, "inferno", False, slug="divBerr"),
    ]


def _render_extent(args):
    """Render domain: --render-extent, else the central ~2000 AU window."""
    if args.render_extent is not None:
        return args.render_extent
    return [-RENDER_HALF, RENDER_HALF, -RENDER_HALF, RENDER_HALF]


def do_render(inputs, labels, args, save, reference=None):
    """Render every field map at the requested times (default: render-n evenly
    spaced), forwarding all --render-* flags. The view defaults to the central
    ~2000 AU window (Wissing & Shen); --render-extent overrides it. When
    `reference` is given and there is a single input, the blessed render grids are
    overlaid (sim | reference | difference)."""
    entries = list(zip(labels, inputs))
    times = args.render_times if args.render_times else None
    extent = _render_extent(args)
    if args.render_extent is None:
        print(f"mhdcollapse: default render window = central 2000 AU "
              f"(+-{RENDER_HALF:.4g} code units); override with --render-extent.")
    units = code_units(inputs[0])       # generic code->cgs factors from the .log
    if units is not None:
        print(f"mhdcollapse: physical units rho*{units.density_gcc:.4g} g/cm^3, "
              f"v*{units.velocity_kms:.4g} km/s, B*{units.b_uG:.4g} uG")
    else:
        print(f"mhdcollapse: no dMsolUnit/dKpcUnit in '{inputs[0]}' .log; "
              f"rendering in code units.", file=sys.stderr)
    overlay = reference is not None and len(inputs) == 1
    for sp in render_specs(units):
        slug = spec_slug(sp)
        sp_save = f"{save}_{slug}" if save else None
        ref_path = (render_ref_overlay_path(inputs[0], reference, slug)
                    if overlay else None)
        render_panels(entries, sp, times=times, n=args.render_n,
                      res=args.render_res, save=sp_save,
                      backend=args.render_backend, project=args.render_project,
                      slice_frac=args.render_slice_frac,
                      extent=extent, aspect=args.render_aspect,
                      render_reference=ref_path,
                      residual=getattr(args, "residual", False))


def _save_extra(args, params):
    """Bless the render-grid baselines for --save-reference (same physical units
    and zoom window the renders use, so the overlay aligns)."""
    units = code_units(args.inputs[0])
    bless_render_specs(args.inputs[0], render_specs(units), args, "mhdcollapse",
                       extent=_render_extent(args))


# --------------------------------------------------------------------------- #
#  Regression baseline (bless / compare the peak-density amplification series)
# --------------------------------------------------------------------------- #

def reference_table(inp):
    """Regression table {x: t, rho_peak: rho_peak(t)/rho_peak(0)} for run `inp`.

    The normalised peak-density amplification is the headline collapse
    diagnostic and dimensionless (avoids the IC's intricate unit system) -> the
    natural CI baseline. None if no usable density series."""
    ts, peaks = density_series(inp)
    if ts is None or not np.isfinite(peaks[0]) or peaks[0] == 0:
        return None
    table = {"x": ts, "rho_peak": peaks / peaks[0]}
    # Energy block (own time grid) for the energy-budget plot overlay.
    te, comps = energy_series_norm(inp)
    if te is not None:
        table["t_energy"] = te
        table["Epot_n"] = comps["|E_pot|"]
        table["Ekin_n"] = comps["E_kin"]
        table["Emag_n"] = comps["E_mag"]
    return table


def _plot(inputs, labels, params, args):
    print("[mhdcollapse] magnetized rotating cloud collapse (Wissing & Shen 2020)")
    ref = getattr(args, "ref_table", None)
    plot_density_vs_time(inputs, labels, save=args.save, ref_table=ref)
    plot_brho(inputs, labels, save=args.save)
    plot_energy_vs_time(inputs, labels, save=args.save, ref_table=ref)
    if args.render:
        do_render(inputs, labels, args, save=args.save,
                  reference=args.reference)
    if args.movie:
        units = code_units(inputs[0])
        movie_render_specs(inputs[0], labels[0], render_specs(units), args,
                           args.save, "mhdcollapse", extent=_render_extent(args))


def main():
    standalone_cli(
        "mhdcollapse", PARAM_SPEC, _plot,
        reference_table=lambda inp, params, args: reference_table(inp),
        description="MHD cloud-collapse analysis (Wissing & Shen 2020).",
        reg_keys=("rho_peak",), add_time=False, save_extra=_save_extra)


if __name__ == "__main__":
    main()
