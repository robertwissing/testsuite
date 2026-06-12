#!/usr/bin/env python3
"""
Analysis for the 2D MHD current-sheet test (IC_setup_currentsheet), from
Gardiner & Stone (2005, JCP 205, 509) and Hopkins & Raives (2016, MNRAS 455,
51). A periodic box (x,y in [-0.5,0.5]) holds two antiparallel field sheets

    By = -B0 for |x|<0.25,  +B0 otherwise,   Bx=Bz=0,  rho=1,  p = beta B0^2/2

with a sinusoidal shear vx = v0 sin(2 pi (y+0.5)) that drives magnetic
reconnection at the two field reversals. There is no analytic solution; this is
a robustness / divergence-handling test. Magnetic islands form and merge, and
the magnetic energy decays as the reversed field reconnects and dissipates.

This script produces:
  * magnetic energy vs time, E_mag/E_mag(0) (with E_kin/E_mag(0)) -- the headline
    Hopkins & Raives diagnostic: reconnection should dissipate E_mag smoothly,
    not violently (a robust scheme avoids spurious blow-ups);
  * the divergence error h|div B|/|B| vs time (mean and max) -- the
    divergence-handling robustness check;
  * field-map renders (--render): the signed By (islands) and the magnetic
    pressure |B|^2/2.

It also supports the regression-baseline workflow (--save-reference / --reg-tol)
so the time series can be locked in for CI/CD.

Usage:
    python IC_analysis_currentsheet.py <run-dir> [-a beta] [-b v0] [--save out]
    python IC_analysis_currentsheet.py <run-dir> --render --render-n 4 --save out
    python IC_analysis_currentsheet.py <run-dir> --save-reference        # bless
    python IC_analysis_currentsheet.py <run-dir> --reg-tol 1e-6          # CI gate
"""

import os
import sys

import numpy as np
import readtipsy as tip

from IC_analysis_framework import (
    Param, standalone_cli,
    RenderSpec, render_panels, aux_divB_error,
    StreamSpec, streamline_panels, aux_component, magnetic_pressure,
)
from IC_analysis_general import (
    find_files, read_vector_aux, read_tufac, read_log_series,
    setup_rcparams, finish_figure_with_legend,
)

GAMMA = 5.0 / 3.0               # IC_setup_currentsheet.gamma
B0 = 1.0                       # field amplitude
REG_NPTS = 101                 # points on the regression time grid


# --------------------------------------------------------------------------- #
#  runtest.sh setup_currentsheet() parameter order (beta, v0).
#  Only carried for CLI/label parity; the analysis normalises by E_mag(0).
# --------------------------------------------------------------------------- #
PARAM_SPEC = [
    Param("a", "beta", 0.1, float, "plasma beta (gas pressure = beta*B0^2/2)"),
    Param("b", "v0", 0.1, float, "shear velocity perturbation amplitude"),
]


# --------------------------------------------------------------------------- #
#  Per-snapshot diagnostics (snapshot fallback for the '.log' series)
# --------------------------------------------------------------------------- #

def _snapshot(fn):
    """(time, tgdata, ngas) for a snapshot."""
    tgdata, _td, _ts, hdr, time = tip.readtipsy(fn)
    return float(np.asarray(time).ravel()[0]), tgdata, int(hdr[2])


def divberr_of(fn):
    """(time, mean, max) of the divergence error h|div B|/|B| for snapshot fn.

    Mean is volume-weighted (V_i = m_i/rho_i); NaN if the DivB/BField aux is
    absent. Uses the framework's standard divBerr callable."""
    t, tgdata, ngas = _snapshot(fn)
    err = np.asarray(aux_divB_error()(fn, tgdata, ngas), float)
    if not np.any(err > 0):
        return t, np.nan, np.nan
    vol = tgdata[:, 0] / tgdata[:, 7]
    return t, float(np.sum(vol*err) / np.sum(vol)), float(np.max(err))


def energies_of(fn, tufac):
    """(time, E_mag, E_kin, E_therm) for snapshot fn.

    E_mag = sum 1/2 |B|^2 V (NaN if BField aux absent); E_kin = sum 1/2 m|v|^2;
    E_therm = sum m u with u = tufac*T (col 8 is temperature)."""
    t, tgdata, ngas = _snapshot(fn)
    m, rho = tgdata[:, 0], tgdata[:, 7]
    e_kin = float(np.sum(0.5 * m * np.sum(tgdata[:, 4:7]**2, axis=1)))
    e_therm = float(tufac * np.sum(m * tgdata[:, 8]))
    B = read_vector_aux(fn, "BField", ngas)
    if B is None:
        return t, np.nan, e_kin, e_therm
    b2 = np.sum(np.asarray(B, float)**2, axis=1)
    e_mag = float(np.sum(0.5 * b2 * (m / rho)))
    return t, e_mag, e_kin, e_therm


def _snapshot_rows(inp):
    """Sorted (t, divb_mean, divb_max, Emag, Ekin, Etherm) per snapshot, or None."""
    try:
        files = find_files(inp)
    except FileNotFoundError:
        return None
    tufac = read_tufac(inp)
    rows = []
    for fn in files:
        try:
            t, dmean, dmax = divberr_of(fn)
            _t, em, ek, et = energies_of(fn, tufac)
        except Exception as e:
            print(f"  failed on {fn}: {e}", file=sys.stderr)
            continue
        rows.append((t, dmean, dmax, em, ek, et))
    if not rows:
        return None
    a = np.array(rows, float)
    return a[np.argsort(a[:, 0])]


def series(inp):
    """Time-series diagnostics for one run as a dict, or None.

    Prefers the dense gasoline '.log' (every timestep) for the divergence error
    (divBerrAvg/Max) and energies (Emag/Ekin/Eth), falling back to the snapshots.
    Keys: 'divberr'=(t,mean,max), 'energy'=(t,Emag,Ekin,Etherm), 'src'={...}."""
    log = read_log_series(
        inp, ["dTime", "divBerrAvg", "divBerrMax", "Emag", "Ekin", "Eth"])
    snap = _snapshot_rows(inp)
    if log is None and snap is None:
        print(f"currentsheet: no snapshots/.log for '{inp}'", file=sys.stderr)
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
        out["energy"] = (snap[:, 0], snap[:, 3], snap[:, 4], snap[:, 5])
        out["src"]["energy"] = "snapshots"
    return out


# --------------------------------------------------------------------------- #
#  Plots
# --------------------------------------------------------------------------- #

def plot_energy_vs_time(runs, labels, save=None, ref_table=None):
    """Energy budget vs time, normalised by E_mag(0). Reconnection dissipates the
    reversed field, so E_mag falls and the released energy is thermalised: the
    thermal GAIN dE_th = E_th(t)-E_th(0) mirrors the E_mag drop (budget closure).
    E_kin (shear + reconnection outflows) is small and shown for completeness.

    When a reference is supplied its blessed E_mag / dE_th / E_kin are overlaid
    (grey: dashed E_mag, dotted dE_th, faint dashed E_kin)."""
    import matplotlib.pyplot as plt
    setup_rcparams()
    fig, ax = plt.subplots(figsize=(8, 6))
    if ref_table is not None and "Emag" in ref_table:
        xr = ref_table["x"]
        ax.plot(xr, ref_table["Emag"], "--", color="0.5", lw=1.4, zorder=5,
                label="reference")
        if "Eth" in ref_table:
            ax.plot(xr, ref_table["Eth"] - ref_table["Eth"][0], ":", color="0.5",
                    lw=1.4, zorder=5)
        if "Ekin" in ref_table:
            ax.plot(xr, ref_table["Ekin"], "--", color="0.5", lw=0.9, alpha=0.7,
                    zorder=5)
    for s, lab in zip(runs, labels):
        if s is None or "energy" not in s:
            continue
        t, emag, ekin, eth = s["energy"]
        e0 = emag[0]
        if not np.isfinite(e0) or e0 == 0:
            print(f"currentsheet[{lab}]: no/zero initial E_mag; skipping.",
                  file=sys.stderr)
            continue
        style = "-" if s["src"].get("energy") == "log" else "o-"
        pre = f"{lab} " if len(runs) > 1 else ""
        line, = ax.plot(t, emag/e0, style, ms=3, label=f"{pre}$E_{{mag}}$")
        c = line.get_color()
        ax.plot(t, (eth - eth[0])/e0, ":", lw=1.8, color=c,
                label=f"{pre}$\\Delta E_{{th}}$")
        ax.plot(t, ekin/e0, "--", lw=1.0, color=c, alpha=0.7,
                label=f"{pre}$E_{{kin}}$")
        print(f"currentsheet[{lab}]: E_mag/E_mag0 {emag[0]/e0:.4g}->{emag[-1]/e0:.4g}"
              f"  dE_th/E_mag0 {(eth[-1]-eth[0])/e0:+.4g}"
              f"  E_kin/E_mag0 max {np.nanmax(ekin)/e0:.4g}"
              f"  ({len(t)} pts from {s['src'].get('energy')})")
    ax.axhline(0.0, color="0.7", lw=0.8)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$E / E_{mag}(0)$")
    ax.set_title("Current sheet: energy budget vs time "
                 r"($E_{mag}$ decay $\to$ heat)")
    fig.tight_layout()
    finish_figure_with_legend(fig, ax, save=(f"{save}_energy.png" if save else None))


def plot_divberr_vs_time(runs, labels, save=None, ref_table=None):
    """Mean and max of h|div B|/|B| vs time (divergence-handling robustness).

    When a reference is supplied its blessed mean/max curves are overlaid (grey
    dashed mean, dotted max)."""
    import matplotlib.pyplot as plt
    setup_rcparams()
    fig, ax = plt.subplots(figsize=(8, 6))
    if ref_table is not None and "divBerrMean" in ref_table:
        ax.semilogy(ref_table["x"], ref_table["divBerrMean"], "--", color="0.5",
                    lw=1.4, zorder=5, label="reference")
        if "divBerrMax" in ref_table:
            ax.semilogy(ref_table["x"], ref_table["divBerrMax"], ":", color="0.5",
                        lw=1.4, zorder=5)
    for s, lab in zip(runs, labels):
        if s is None or "divberr" not in s:
            continue
        t, mean, mx = s["divberr"]
        if np.all(np.isnan(mean)):
            print(f"currentsheet[{lab}]: no DivB/BField aux; skipping.",
                  file=sys.stderr)
            continue
        src = s["src"].get("divberr", "snapshots")
        mstyle, xstyle = ("-", "--") if src == "log" else ("o-", "s--")
        line, = ax.semilogy(t, mean, mstyle, ms=3, label=f"{lab} mean")
        ax.semilogy(t, mx, xstyle, ms=3, color=line.get_color(), label=f"{lab} max")
        print(f"currentsheet[{lab}]: divBerr mean {mean[0]:.3g}->{mean[-1]:.3g}  "
              f"max {mx[0]:.3g}->{mx[-1]:.3g}  ({len(t)} pts from {src})")
    ax.set_xlabel("t")
    ax.set_ylabel(r"$h\,|\nabla\!\cdot\!B|/|B|$")
    ax.set_title("Current sheet: divergence error vs time")
    fig.tight_layout()
    finish_figure_with_legend(fig, ax,
                              save=(f"{save}_divberr.png" if save else None))


# --------------------------------------------------------------------------- #
#  Field-map renders
# --------------------------------------------------------------------------- #

def render_specs():
    return [
        RenderSpec(1, 2, aux_component("BField", 1), r"$B_y$", False, "average",
                   None, "RdBu_r", True),
        RenderSpec(1, 2, magnetic_pressure(), r"$\frac{1}{2}|B|^2$", False,
                   "average", None, "inferno", False),
    ]


def do_render(inputs, labels, args, save):
    """Render By and magnetic pressure at the requested times."""
    entries = list(zip(labels, inputs))
    extent = args.render_extent if args.render_extent else [-0.5, 0.5, -0.5, 0.5]
    for sp, slug in zip(render_specs(), ("By", "magpressure")):
        render_panels(entries, sp, times=args.render_times, n=args.render_n,
                      res=args.render_res,
                      save=(f"{save}_{slug}" if save else None),
                      backend=args.render_backend, project=args.render_project,
                      slice_frac=args.render_slice_frac, extent=extent,
                      aspect=args.render_aspect)


# --------------------------------------------------------------------------- #
#  Magnetic-field streamlines (vector field)
# --------------------------------------------------------------------------- #

def plot_bfield_streamlines(inputs, labels, args, save):
    """Magnetic-field streamlines (Bx, By) in the x-y plane on a white
    background, at the selected snapshot times -- directly shows the reconnection
    islands (X- and O-points) at the two sheets x = +/-0.25.

    Delegates to the framework `streamline_panels` (which was generalized FROM
    this plot); the per-test bits here are the StreamSpec (Bx, By from the
    BField aux) and the initial sheet-location markers."""
    spec = StreamSpec(1, 2, aux_component("BField", 0),
                      aux_component("BField", 1), "B")

    def sheets(ax):                                   # initial sheet locations
        for xs in (-0.25, 0.25):
            ax.axvline(xs, color="0.6", ls=":", lw=0.8, zorder=1)

    extent = args.render_extent if args.render_extent else [-0.5, 0.5, -0.5, 0.5]
    streamline_panels(list(zip(labels, inputs)), spec,
                      times=args.render_times, n=args.render_n,
                      res=min(args.render_res or 64, 128),
                      save=(f"{save}_bstream" if save else None),
                      extent=extent, decorate=sheets,
                      title="Current sheet: magnetic-field streamlines")


# --------------------------------------------------------------------------- #
#  Regression baseline (CI/CD)
# --------------------------------------------------------------------------- #

def reference_table(inp):
    """Regression table on a fixed time grid: {x: t, Emag, Ekin, divBerrMean,
    divBerrMax} with the energies normalised by E_mag(0). Resampled so a saved
    baseline can be compared across runs regardless of log/snapshot sampling.
    Returns None if no usable series."""
    s = series(inp)
    if s is None or "energy" not in s:
        return None
    t, emag, ekin, eth = s["energy"]
    e0 = emag[0] if (np.isfinite(emag[0]) and emag[0] != 0) else 1.0
    tg = np.linspace(float(t[0]), float(t[-1]), REG_NPTS)
    tbl = {"x": tg,
           "Emag": np.interp(tg, t, emag / e0),
           "Ekin": np.interp(tg, t, ekin / e0),
           "Eth": np.interp(tg, t, eth / e0)}
    if "divberr" in s:
        td, dm, dx = s["divberr"]
        tbl["divBerrMean"] = np.interp(tg, td, dm)
        tbl["divBerrMax"] = np.interp(tg, td, dx)
    return tbl


def _plot(inputs, labels, params, args):
    print(f"[currentsheet] beta={params['beta']}  v0={params['v0']}  B0={B0:g}")
    runs = [series(inp) for inp in inputs]
    if all(s is None for s in runs):
        print("currentsheet: no data to plot; exiting.", file=sys.stderr)
        sys.exit(1)
    ref = getattr(args, "ref_table", None)
    plot_energy_vs_time(runs, labels, save=args.save, ref_table=ref)
    plot_divberr_vs_time(runs, labels, save=args.save, ref_table=ref)
    if args.render:
        do_render(inputs, labels, args, save=args.save)
        plot_bfield_streamlines(inputs, labels, args, save=args.save)


def main():
    standalone_cli(
        "currentsheet", PARAM_SPEC, _plot,
        reference_table=lambda inp, params, args: reference_table(inp),
        description="MHD current-sheet analysis "
                    "(Gardiner & Stone 2005 / Hopkins & Raives 2016).",
        reg_keys=("Emag", "Ekin", "Eth", "divBerrMean", "divBerrMax"),
        add_time=False)


if __name__ == "__main__":
    main()
