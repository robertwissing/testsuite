#!/usr/bin/env python3
"""
Analysis for the divergence-cleaning test (IC_setup_mhddivtest), from Tricco &
Price (2012, JCP 231, 7214). A localised magnetic-field bump

    Bx = Bz0 (1 - r^2/r0^2)^2   (r^2/r0^2 < 1),   Bz = Bz0,  By = 0,  v = 0

(Bz0 = 0.01/sqrt(4 pi), r0 = 1/sqrt(8)) seeds a peak of non-zero div(B) on an
otherwise static background (with an optional x-density contrast `rhodiff`).
Constrained hyperbolic/parabolic cleaning should propagate the divergence away
and dissipate it, so the divergence error must DECAY in time. This script
produces:

  * div(B) field-map renders at selected times (--render),
  * the divergence error h|div B|/|B| vs time -- volume-weighted mean AND maximum
    (the headline Tricco & Price diagnostic; both should fall),
  * the energy budget vs time, normalised by the initial magnetic energy:
    magnetic E_mag = sum 1/2|B|^2 V, cleaning E_psi = sum 1/2 psi^2 V / c_h^2
    (c_h = fast magnetosonic speed), kinetic E_kin, and the thermal CHANGE
    dE_therm -- the divergent field energy is cleaned into psi and dissipated as
    heat, while E_kin stays ~0.

Standalone (reuses the framework loaders + render_panels, not run_cli).

Usage:
    python IC_analysis_mhddivtest.py <run-dir> [-a 1] [--save out]
    python IC_analysis_mhddivtest.py <run-dir> --render --render-n 4 --save out
"""

import os
import sys

import numpy as np
import readtipsy as tip

from IC_analysis_framework import (
    Param, standalone_cli,
    RenderSpec, render_panels, aux_divB_error, spec_slug,
    bless_render_specs, render_ref_overlay_path, movie_render_specs,
)
from IC_analysis_general import (
    find_files, try_aux, read_vector_aux, read_log_series,
    setup_rcparams, finish_figure_with_legend,
)

GAMMA = 5.0 / 3.0                       # IC_setup_mhddivtest.gamma
BZ0 = 0.01 / np.sqrt(4.0 * np.pi)       # background / peak field amplitude
R0 = 1.0 / np.sqrt(8.0)                # bump radius
PRZERO = 6.0                           # initial (uniform) pressure

# runtest.sh setup_mhddivtest() order. Only rhodiff/square matter physically;
# smooth/cosmo are carried for CLI parity.
PARAM_SPEC = [
    Param("a", "rhodiff", 1.0, float, "x-density contrast rhodens/rhozero"),
    Param("b", "smooth", 1, int, "smooth (1) vs sharp (0) density step"),
    Param("c", "square", 1, int, "2-D square (1) vs 3-D sphere (0) bump"),
    Param("d", "cosmo", 0, int, "cosmological run flag (unused in analysis)"),
]


# --------------------------------------------------------------------------- #
#  div(B) render
# --------------------------------------------------------------------------- #

def render_spec():
    """Signed div(B) in the x-y plane (diverging cmap, zero-centred)."""
    def q_divB(fn, tgdata, ngas):
        d = try_aux(fn, "DivB", ngas)
        return np.zeros(len(tgdata)) if d is None else np.asarray(d, float)
    return RenderSpec(1, 2, q_divB, r"$\nabla\!\cdot\!B$", False, "average",
                      None, "RdBu_r", True, slug="divB")


def do_render(inputs, labels, args, save, reference=None):
    """Render div(B) at the requested times, forwarding the --render-* flags.

    When `reference` is given and there is a single input, the blessed div(B)
    render grids are overlaid (sim | reference | difference)."""
    entries = list(zip(labels, inputs))
    sp = render_spec()
    overlay = reference is not None and len(inputs) == 1
    ref_path = (render_ref_overlay_path(inputs[0], reference, spec_slug(sp))
                if overlay else None)
    render_panels(entries, sp,
                  times=args.render_times, n=args.render_n,
                  res=args.render_res,
                  save=(f"{save}_{spec_slug(sp)}" if save else None),
                  backend=args.render_backend, project=args.render_project,
                  slice_frac=args.render_slice_frac, extent=args.render_extent,
                  aspect=args.render_aspect, render_reference=ref_path,
                  residual=getattr(args, "residual", False))


def _save_extra(args, params):
    """Bless the div(B) render-grid baseline for --save-reference."""
    bless_render_specs(args.inputs[0], [render_spec()], args, "mhddivtest")


# --------------------------------------------------------------------------- #
#  Per-snapshot diagnostics
# --------------------------------------------------------------------------- #

def _snapshot(fn):
    """(time, tgdata, ngas) for a snapshot."""
    tgdata, _td, _ts, hdr, time = tip.readtipsy(fn)
    return float(np.asarray(time).ravel()[0]), tgdata, int(hdr[2])


def divberr_series(fn):
    """(time, mean, max) of the divergence error h|div B|/|B| for snapshot fn.

    Mean is volume-weighted (V_i = m_i/rho_i); both are NaN if the DivB/BField
    aux is missing. Uses the framework's standard divBerr callable."""
    t, tgdata, ngas = _snapshot(fn)
    err = aux_divB_error()(fn, tgdata, ngas)        # h|div B|/|B|, 0 where |B|=0
    err = np.asarray(err, float)
    if not np.any(err > 0):
        return t, np.nan, np.nan
    vol = tgdata[:, 0] / tgdata[:, 7]
    mean = float(np.sum(vol * err) / np.sum(vol))
    return t, mean, float(np.max(err))


def energies(fn):
    """(time, E_mag, E_psi, E_kin, E_therm) for snapshot fn.

    E_mag = sum 1/2 |B|^2 V, E_psi = sum 1/2 psi^2 V / c_h^2 (c_h^2 = gamma(gamma
    -1)u + |B|^2/rho, the fast magnetosonic speed), E_kin = sum 1/2 m|v|^2,
    E_therm = sum m u. E_mag/E_psi are NaN if the relevant aux is missing."""
    t, tgdata, ngas = _snapshot(fn)
    m, rho = tgdata[:, 0], tgdata[:, 7]
    vol = m / rho
    u = tgdata[:, 8]                                # specific internal energy
    e_kin = float(np.sum(0.5 * m * np.sum(tgdata[:, 4:7] ** 2, axis=1)))
    e_therm = float(np.sum(m * u))
    B = read_vector_aux(fn, "BField", ngas)
    if B is None:
        return t, np.nan, np.nan, e_kin, e_therm
    b2 = np.sum(np.asarray(B, float) ** 2, axis=1)
    e_mag = float(np.sum(0.5 * b2 * vol))
    psi = try_aux(fn, "BClean", ngas)
    if psi is None:
        return t, e_mag, np.nan, e_kin, e_therm
    ch2 = GAMMA * (GAMMA - 1.0) * u + b2 / rho       # fast magnetosonic speed^2
    with np.errstate(divide="ignore", invalid="ignore"):
        e_psi = float(np.sum(0.5 * np.asarray(psi, float) ** 2 * vol
                             / np.where(ch2 > 0, ch2, np.nan)))
    return t, e_mag, e_psi, e_kin, e_therm


def _snapshot_rows(inp):
    """Sorted-by-time per-snapshot array (t, divb_mean, divb_max, Emag, Epsi,
    Ekin, Etherm) for run `inp`, or None. The fallback / E_psi source."""
    try:
        files = find_files(inp)
    except FileNotFoundError:
        return None
    rows = []
    for fn in files:
        try:
            t, dmean, dmax = divberr_series(fn)
            _t, em, ep, ek, et = energies(fn)
        except Exception as e:
            print(f"  failed on {fn}: {e}", file=sys.stderr)
            continue
        rows.append((t, dmean, dmax, em, ep, ek, et))
    if not rows:
        return None
    a = np.array(rows, float)
    return a[np.argsort(a[:, 0])]


def _series(inp):
    """Time-series diagnostics for one run as a dict, or None.

    Prefers the dense gasoline '.log' (logged every timestep) for the divergence
    error (divBerrAvg/divBerrMax) and the energies (Emag/Ekin/Eth); the cleaning
    energy E_psi is not logged, so it is always taken from the snapshots (BClean
    aux). Keys:
      'divberr' : (t, mean, max)             '.log' preferred
      'energy'  : (t, Emag, Ekin, Etherm)    '.log' preferred
      'psi'     : (t, Epsi) or None          snapshots only
      'src'     : {'divberr': ..., 'energy': ...}  source per series ('log'/'snapshots')
    """
    log = read_log_series(
        inp, ["dTime", "divBerrAvg", "divBerrMax", "Emag", "Ekin", "Eth"])
    snap = _snapshot_rows(inp)
    if log is None and snap is None:
        print(f"mhddivtest: no snapshots/.log for '{inp}'", file=sys.stderr)
        return None

    out = {"src": {}}
    if log is not None and np.any(log["divBerrAvg"] > 0):
        out["divberr"] = (log["dTime"], log["divBerrAvg"], log["divBerrMax"])
        out["src"]["divberr"] = "log"
    elif snap is not None:
        out["divberr"] = (snap[:, 0], snap[:, 1], snap[:, 2])
        out["src"]["divberr"] = "snapshots"

    if log is not None and np.any(log["Emag"] > 0):
        out["energy"] = (log["dTime"], log["Emag"], log["Ekin"], log["Eth"])
        out["src"]["energy"] = "log"
    elif snap is not None:
        out["energy"] = (snap[:, 0], snap[:, 3], snap[:, 5], snap[:, 6])
        out["src"]["energy"] = "snapshots"

    if snap is not None and np.any(np.isfinite(snap[:, 4])):
        out["psi"] = (snap[:, 0], snap[:, 4])      # E_psi from the BClean aux
    else:
        out["psi"] = None
    return out


# --------------------------------------------------------------------------- #
#  Plots
# --------------------------------------------------------------------------- #

def _plot_divberr_one(series, labels, which, save=None, ref_table=None):
    """One divergence-error curve (`which` = 'mean' or 'max') of h|div B|/|B|
    vs time, log-y. The data summary is printed once, on the 'mean' pass. When a
    reference is supplied its blessed curve (divBerrAvg/divBerrMax) is overlaid."""
    import matplotlib.pyplot as plt
    setup_rcparams()
    fig, ax = plt.subplots(figsize=(8, 6))
    refkey = "divBerrAvg" if which == "mean" else "divBerrMax"
    if ref_table is not None and refkey in ref_table:
        ax.semilogy(ref_table["x"], ref_table[refkey], "--", color="0.35",
                    lw=1.5, zorder=5, label="reference")
    for s, lab in zip(series, labels):
        if s is None or "divberr" not in s:
            continue
        t, mean, mx = s["divberr"]
        if np.all(np.isnan(mean)):
            if which == "mean":
                print(f"mhddivtest[{lab}]: no DivB/BField aux; skipping.",
                      file=sys.stderr)
            continue
        src = s["src"].get("divberr", "snapshots")
        style = "-" if src == "log" else "o-"
        y = mean if which == "mean" else mx
        ax.semilogy(t, y, style, ms=3, label=f"{lab} {which}")
        if which == "mean":
            print(f"mhddivtest[{lab}]: divBerr mean {mean[0]:.3g}->{mean[-1]:.3g}  "
                  f"max {mx[0]:.3g}->{mx[-1]:.3g}  ({len(t)} pts from {src})")
    ax.set_xlabel("t")
    ax.set_ylabel(r"$h\,|\nabla\!\cdot\!B|/|B|$")
    ax.set_title(f"Divergence error ({which}) vs time (cleaning should reduce it)")
    fig.tight_layout()
    finish_figure_with_legend(
        fig, ax, save=(f"{save}_divberr_{which}.png" if save else None))


def plot_divberr_vs_time(series, labels, save=None, ref_table=None):
    """Divergence error h|div B|/|B| vs time, as two separate plots: mean-only
    (_divberr_mean.png) and max-only (_divberr_max.png); both should decay.

    Uses the dense '.log' divBerrAvg/divBerrMax columns when available."""
    _plot_divberr_one(series, labels, "mean", save=save, ref_table=ref_table)
    _plot_divberr_one(series, labels, "max", save=save, ref_table=ref_table)


def plot_energy_vs_time(series, labels, save=None, ref_table=None):
    """Magnetic energy E_mag/E_mag(0) vs time -- its decay is the cleaning of
    divergent field energy (only a few %). The cleaning/kinetic/thermal budget
    (E_psi, E_kin, dE_therm) is still reported per run on the console. When a
    reference is supplied its blessed E_mag/E_mag(0) is overlaid (dashed)."""
    import matplotlib.pyplot as plt
    setup_rcparams()
    fig, axM = plt.subplots(figsize=(8, 6))
    if ref_table is not None and "Emag" in ref_table:
        axM.plot(ref_table["t_energy"], ref_table["Emag"], "--", color="0.35",
                 lw=1.5, zorder=5, label="reference")
    for s, lab in zip(series, labels):
        if s is None or "energy" not in s:
            continue
        t, emag, ekin, eth = s["energy"]
        e0 = emag[0]
        if not np.isfinite(e0) or e0 == 0:
            print(f"mhddivtest[{lab}]: no/zero initial E_mag; skipping energy.",
                  file=sys.stderr)
            continue
        # E_mag follows the (dense) '.log'; E_psi is snapshot-only.
        estyle = "-" if s["src"].get("energy") == "log" else "o-"
        pre = f"{lab} " if len(series) > 1 else ""
        axM.plot(t, emag / e0, estyle, ms=3, label=f"{pre}$E_{{mag}}$")
        epm = float("nan")
        if s.get("psi") is not None:
            _tpsi, epsi = s["psi"]
            if np.any(np.isfinite(epsi)):
                epm = np.nanmax(epsi) / e0
        print(f"mhddivtest[{lab}]: E_mag/E_mag0 {emag[0]/e0:.3g}->{emag[-1]/e0:.3g}"
              f"  E_psi/E_mag0 max {epm:.3g}"
              f"  E_kin/E_mag0 max {np.nanmax(ekin)/e0:.3g}"
              f"  ({len(t)} pts from {s['src'].get('energy')})")
    axM.set_ylabel(r"$E_{mag}/E_{mag}(0)$")
    axM.set_xlabel("t")
    axM.set_title("Magnetic energy vs time (cleaning of divergent field energy)")
    fig.tight_layout()
    finish_figure_with_legend(fig, axM,
                              save=(f"{save}_energy.png" if save else None))


# --------------------------------------------------------------------------- #
#  Regression baseline (bless / compare the divergence-error decay series)
# --------------------------------------------------------------------------- #

def reference_table(inp):
    """Regression table {x: t, divBerrAvg, divBerrMax} for run `inp`.

    The (dimensionless) divergence-error decay h|div B|/|B| is the headline
    Tricco & Price diagnostic and shares one x-grid (mean and max together), so
    it makes the natural CI baseline. None if no divBerr series."""
    s = _series(inp)
    if s is None or "divberr" not in s:
        return None
    t, mean, mx = s["divberr"]
    table = {"x": np.asarray(t, dtype=float),
             "divBerrAvg": np.asarray(mean, dtype=float),
             "divBerrMax": np.asarray(mx, dtype=float)}
    # Energy block (its own time grid) for the E_mag/E_mag(0) plot overlay.
    if "energy" in s:
        te, emag, _ek, _et = s["energy"]
        e0 = emag[0]
        if np.isfinite(e0) and e0 != 0:
            table["t_energy"] = np.asarray(te, dtype=float)
            table["Emag"] = np.asarray(emag, dtype=float) / e0
    return table


def _plot(inputs, labels, params, args):
    print(f"[mhddivtest] rhodiff={params['rhodiff']}  square={params['square']}  "
          f"Bz0={BZ0:.4g}  r0={R0:.4g}")
    series = [_series(inp) for inp in inputs]
    if all(s is None for s in series):
        print("mhddivtest: no data to plot; exiting.", file=sys.stderr)
        sys.exit(1)
    ref = getattr(args, "ref_table", None)
    plot_divberr_vs_time(series, labels, save=args.save, ref_table=ref)
    plot_energy_vs_time(series, labels, save=args.save, ref_table=ref)
    if args.render:
        do_render(inputs, labels, args, save=args.save,
                  reference=args.reference)
    if args.movie:
        movie_render_specs(inputs[0], labels[0], [render_spec()], args, args.save,
                           "mhddivtest")


def main():
    standalone_cli(
        "mhddivtest", PARAM_SPEC, _plot,
        reference_table=lambda inp, params, args: reference_table(inp),
        description="Divergence-cleaning test analysis (Tricco & Price 2012).",
        reg_keys=("divBerrAvg", "divBerrMax"), add_time=False,
        save_extra=_save_extra)


if __name__ == "__main__":
    main()
