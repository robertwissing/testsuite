#!/usr/bin/env python3
"""
Analysis for the Evrard (1988) adiabatic spherical-collapse test
(IC_setup_evrard). A cold (u0 = 0.05 G M / R) gas sphere of mass M=1, radius
R=1, with the 1/r density profile

    rho(r) = M / (2 pi R^2 r)        (softened by eps=0.02 in the IC),  v = 0

collapses under self-gravity (G=1, gamma=5/3). An accretion shock forms and
moves outward; by t ~= 0.8 the standard radial profiles (density, radial
velocity, entropy) have developed and are compared here against a high-
resolution 1D PPM reference (8192 grid points, Mandel et al. 2023), stored as
CSVs in reference_solution_highres/evrard/.

This script produces:
  * radial profiles at the comparison time (default t=0.8) vs the PPM reference:
      density   rho(r)                          (log-log)
      radial vel  v_r(r) = (x vx + y vy + z vz)/r
      entropy   A(r) = p/rho^gamma = (gamma-1) u rho^(1-gamma)   (log-log)
    with a volume-weighted L1 error vs the reference over the overlap radii;
  * the energy budget vs time -- kinetic, thermal, potential and total -- the
    classic Evrard conservation diagnostic (potential released by collapse is
    converted to kinetic then thermalised at the shock; the total should be
    ~conserved). Read from the dense gasoline '.log' (Ekin/Epot/Eth/Etot).
  * optional density field-map renders (--render).

Standalone (reuses the framework loaders + render_panels).

Usage:
    python IC_analysis_evrard.py <run-dir> [--time 0.8] [--save out]
    python IC_analysis_evrard.py <run-dir> --refdir <path> --save out
    python IC_analysis_evrard.py <run-dir> --render --render-n 4 --save out
"""

import os
import sys

import numpy as np
import readtipsy as tip

from IC_analysis_framework import (
    Param, standalone_cli, find_reference,
    RenderSpec, render_panels,
    select_at_time, binned_profile, series_from_log_or_snapshots,
)
from IC_analysis_general import (
    loaddata, read_tufac,
    setup_rcparams, finish_figure_with_legend,
)

GAMMA = 5.0 / 3.0               # IC_setup_evrard.gamma
M, R, G = 1.0, 1.0, 1.0        # mass, radius, gravitational constant (code units)
EPS = 0.02                     # density softening in the IC (IC_setup_evrard.eps)
DEFAULT_TIME = 0.8             # standard Evrard comparison time (Mandel et al. 2023)

# Default reference directory: <testsuite>/reference_solution_highres/evrard,
# i.e. one level up from setupfiles/. Overridable with --refdir.
DEFAULT_REFDIR = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "..", "reference_solution_highres", "evrard"))

# Reference CSV (x=radius, y=quantity) for each profile panel.
REF_FILES = {"rho": "evrardrho.csv", "vel": "evrardvel.csv",
             "entrop": "evrardentrop.csv"}

# The setup carries no tunable parameters (runtest setup_evrard sets EXTRACONDIC
# empty); the spec is kept empty for CLI parity with the other analyses.
PARAM_SPEC = []


# --------------------------------------------------------------------------- #
#  Reference (high-res 1D PPM) loading
# --------------------------------------------------------------------------- #

def load_reference(refdir):
    """{name: (x, y)} for the reference CSVs found in `refdir`, sorted by x.

    Missing files are skipped (with a note); returns {} if none are found."""
    ref = {}
    for name, fname in REF_FILES.items():
        path = os.path.join(refdir, fname)
        if not os.path.exists(path):
            print(f"evrard: reference '{fname}' not found in {refdir}; skipping.",
                  file=sys.stderr)
            continue
        data = np.genfromtxt(path, delimiter=",", skip_header=1)
        data = data[np.argsort(data[:, 0])]      # ensure monotone x for interp
        ref[name] = (data[:, 0], data[:, 1])
    return ref


# --------------------------------------------------------------------------- #
#  Per-snapshot fields
# --------------------------------------------------------------------------- #

def profile(fn, tufac):
    """Per-particle (r, rho, v_r, u, entropy, volume) and the snapshot time.

    r is spherical radius from the origin; v_r the outward radial velocity;
    u = tufac*T the specific internal energy (col 8 is temperature); entropy is
    the adiabat A = p/rho^gamma = (gamma-1) u rho^(1-gamma)."""
    tgdata, _td, _ts, _hdr, time, _N, _ng, _nd, _ns, _h = loaddata(fn)
    pos = tgdata[:, 1:4]
    vel = tgdata[:, 4:7]
    rho = tgdata[:, 7]
    u = tufac * tgdata[:, 8]
    r = np.sqrt(np.sum(pos**2, axis=1))
    with np.errstate(invalid="ignore", divide="ignore"):
        v_r = np.where(r > 0, np.sum(pos*vel, axis=1) / r, 0.0)
    entropy = (GAMMA - 1.0) * u * rho ** (1.0 - GAMMA)
    vol = tgdata[:, 0] / rho                      # m/rho = particle volume
    return (r, rho, v_r, u, entropy, vol,
            float(np.asarray(time).ravel()[0]))


# --------------------------------------------------------------------------- #
#  Regression baseline (CI/CD): binned radial profile at the comparison time
# --------------------------------------------------------------------------- #

# Fixed log-spaced radial grid for the regression table -- resolution- and
# code-independent so a saved baseline can be compared across runs.
REG_NBINS = 40
REG_RMIN, REG_RMAX = 0.01, 1.0


def reference_table(inp, t):
    """Binned radial profile at time t as a regression table:
    {x: bin centres, rho, v_r, entropy} (median per log-r bin, NaN where empty).

    This is the run's OWN profile on a fixed grid -- the regression baseline for
    CI (does a code change move the result?) -- distinct from the physical PPM
    reference overlaid in the profile plot. Returns None if no snapshot."""
    sel = select_at_time(inp, t)
    if sel is None:
        return None
    fn, _tt = sel
    r, rho, v_r, _u, entropy, _vol, _time = profile(fn, read_tufac(inp))
    edges = np.geomspace(REG_RMIN, REG_RMAX, REG_NBINS + 1)
    centres, vals = binned_profile(
        r, {"rho": rho, "v_r": v_r, "entropy": entropy},
        edges=edges, stat="median", clip_outside=False, drop_empty=False,
        geom_centres=True)
    table = {"x": centres, "rho": vals["rho"], "v_r": vals["v_r"],
             "entropy": vals["entropy"]}
    # Energy block (own time grid) for the energy-budget plot overlay.
    ts, ek, eth, ep, et, _src = energy_series(inp, read_tufac(inp))
    if ts is not None:
        table["t_energy"] = ts
        table["Etot"] = et
        table["Ekin"] = ek
        table["Eth"] = eth
        table["Epot"] = ep
    return table


def reg_residual(table, ref, key):
    """(max_abs, rms) of table vs reference for `key`, over bins finite in both.

    Both tables use the same fixed grid (REG_*), so compare element-wise."""
    if table is None or ref is None or key not in table or key not in ref:
        return None
    a, b = np.asarray(table[key]), np.asarray(ref[key])
    if a.shape != b.shape:
        return None
    m = np.isfinite(a) & np.isfinite(b)
    if not np.any(m):
        return None
    d = a[m] - b[m]
    return float(np.max(np.abs(d))), float(np.sqrt(np.mean(d ** 2)))


def compare_to_baseline(inp, t, explicit=None):
    """Print per-metric residuals of run `inp` vs its sibling regression baseline.

    Returns the max RMS residual across metrics (None if no baseline found)."""
    ref_table, _path = find_reference(inp, explicit)
    if ref_table is None:
        return None
    table = reference_table(inp, t)
    worst = 0.0
    for key in ("rho", "v_r", "entropy"):
        res = reg_residual(table, ref_table, key)
        if res is None:
            continue
        mx, rms = res
        worst = max(worst, rms)
        print(f"evrard[regression] {key:<8s} vs baseline: max|d|={mx:.4g}  rms={rms:.4g}")
    return worst


# --------------------------------------------------------------------------- #
#  L1 error vs the reference
# --------------------------------------------------------------------------- #

def ref_at(ref_xy, r, logx=False, logy=False):
    """Reference y interpolated onto radii `r`, plus a mask of the radii that
    fall within the reference's x-range (outside -> NaN, excluded from L1)."""
    rx, ry = ref_xy
    inside = (r >= rx[0]) & (r <= rx[-1])
    xi = np.log10(np.clip(r, 1e-30, None)) if logx else r
    xp = np.log10(rx) if logx else rx
    fp = np.log10(np.clip(ry, 1e-30, None)) if logy else ry
    yi = np.interp(xi, xp, fp)
    yi = 10.0**yi if logy else yi
    return np.where(inside, yi, np.nan), inside


def l1_vs_ref(r, field, ref_xy, vol, logx=False, logspace=False):
    """Volume-weighted L1 of `field` vs the reference over the overlap radii.

    `logspace`: compare log10(field) vs log10(ref) (for quantities spanning
    decades, e.g. density / entropy)."""
    ya, inside = ref_at(ref_xy, r, logx=logx, logy=logspace)
    sel = inside & np.isfinite(ya) & np.isfinite(field)
    if logspace:
        sel &= field > 0
    if not np.any(sel):
        return float("nan")
    a = np.log10(field[sel]) if logspace else field[sel]
    b = np.log10(ya[sel]) if logspace else ya[sel]
    return float(np.sum(vol[sel]*np.abs(a - b)) / np.sum(vol[sel]))


# --------------------------------------------------------------------------- #
#  Profiles: rho, v_r, entropy vs r, with reference overlay + L1
# --------------------------------------------------------------------------- #

def plot_profiles(inputs, labels, t, ref, save=None, reg_ref=None):
    import matplotlib.pyplot as plt
    setup_rcparams()
    fig, (axRho, axVel, axEnt) = plt.subplots(3, 1, figsize=(8, 13), sharex=True)
    axRho.set_xscale("log")
    axRho.set_yscale("log")                  # density spans many decades
    axEnt.set_yscale("log")                  # entropy spans many decades

    # Blessed regression-baseline overlay (the run's OWN binned profile on the
    # fixed grid; grey dashed -- distinct from the black physical PPM reference).
    if reg_ref is not None and "rho" in reg_ref:
        for ax, key in ((axRho, "rho"), (axVel, "v_r"), (axEnt, "entropy")):
            if key in reg_ref:
                ax.plot(reg_ref["x"], reg_ref[key], "--", color="0.5", lw=1.4,
                        zorder=4,
                        label=("blessed baseline" if ax is axRho else None))

    l1_lines = []
    for inp, lab in zip(inputs, labels):
        sel = select_at_time(inp, t)
        if sel is None:
            continue
        fn, tt = sel
        tufac = read_tufac(inp)
        r, rho, v_r, u, entropy, vol, _time = profile(fn, tufac)

        axRho.scatter(r, rho, s=3, alpha=0.4, rasterized=True, label=lab)
        axVel.scatter(r, v_r, s=3, alpha=0.4, rasterized=True, label=lab)
        axEnt.scatter(r, entropy, s=3, alpha=0.4, rasterized=True, label=lab)

        bits = [f"{lab}: t={tt:.4g}"]
        if "rho" in ref:
            l1 = l1_vs_ref(r, rho, ref["rho"], vol, logx=True, logspace=True)
            bits.append(f"L1(log rho)={l1:.3g}")
        if "vel" in ref:
            l1 = l1_vs_ref(r, v_r, ref["vel"], vol, logx=True, logspace=False)
            bits.append(f"L1(v_r)={l1:.3g}")
        if "entrop" in ref:
            l1 = l1_vs_ref(r, entropy, ref["entrop"], vol, logx=True, logspace=True)
            bits.append(f"L1(log A)={l1:.3g}")
        print(f"evrard[{lab}]: " + "  ".join(bits[1:]) + f"  (t={tt:.4g})")
        l1_lines.append("  ".join(bits))

    # Reference overlay (high-res 1D PPM), as a black line on each panel.
    if "rho" in ref:
        axRho.plot(ref["rho"][0], ref["rho"][1], "k-", lw=1.8, label="PPM 8192 (Mandel et al. 2023)")
    if "vel" in ref:
        axVel.plot(ref["vel"][0], ref["vel"][1], "k-", lw=1.8, label="PPM 8192 (Mandel et al. 2023)")
    if "entrop" in ref:
        axEnt.plot(ref["entrop"][0], ref["entrop"][1], "k-", lw=1.8,
                   label="PPM 8192 (Mandel et al. 2023)")

    if l1_lines:
        abox = dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85)
        axRho.text(0.02, 0.04, "\n".join(l1_lines), transform=axRho.transAxes,
                   fontsize=8, va="bottom", ha="left", bbox=abox)

    axRho.set_ylabel(r"$\rho$")
    axVel.set_ylabel(r"$v_r$")
    axEnt.set_ylabel(r"$A=p/\rho^{\gamma}$")
    axEnt.set_xlabel(r"$r=\sqrt{x^2+y^2+z^2}$")
    axEnt.set_xlim(5e-3, 1.2)
    fig.suptitle(f"Evrard collapse profiles vs high-res PPM (t={t:g})")
    fig.tight_layout()
    finish_figure_with_legend(fig, axRho,
                              save=(f"{save}_profile.png" if save else None))


# --------------------------------------------------------------------------- #
#  Energy budget vs time (conservation)
# --------------------------------------------------------------------------- #

def energy_series(inp, tufac):
    """(t, Ekin, Eth, Epot, Etot, source) for run `inp`.

    Prefers the dense gasoline '.log' (Ekin/Eth/Epot/Etot logged every step);
    falls back to per-snapshot summation, with the potential energy from the
    snapshot phi column (1/2 sum m phi). Returns Nones if no data."""
    def _snap(fn):
        tgdata, _td, _ts, _hdr, time = tip.readtipsy(fn)
        tt = float(np.asarray(time).ravel()[0])
        m = tgdata[:, 0]
        ek = float(np.sum(0.5 * m * np.sum(tgdata[:, 4:7]**2, axis=1)))
        eth = float(tufac * np.sum(m * tgdata[:, 8]))
        # phi (potential) is the last tipsy gas column when present.
        ep = float(0.5 * np.sum(m * tgdata[:, 11])) if tgdata.shape[1] > 11 else np.nan
        return tt, ek, eth, ep, ek + eth + ep

    return series_from_log_or_snapshots(
        inp, ["dTime", "Ekin", "Eth", "Epot", "Etot"], _snap,
        valid=lambda lg: np.any(lg["Etot"] != 0))


def plot_energy_vs_time(inputs, labels, save=None, ref_table=None):
    """Kinetic / thermal / potential / total energy vs time (Evrard
    conservation diagnostic). When a reference is supplied its blessed
    total/kinetic/thermal/potential histories are overlaid (grey)."""
    import matplotlib.pyplot as plt
    setup_rcparams()
    fig, ax = plt.subplots(figsize=(8, 6))
    if ref_table is not None and "Etot" in ref_table:
        te = ref_table["t_energy"]
        ax.plot(te, ref_table["Etot"], "--", color="0.5", lw=1.3, zorder=5,
                label="reference")
        for key, st in (("Ekin", "--"), ("Eth", ":"), ("Epot", "-.")):
            if key in ref_table:
                ax.plot(te, ref_table[key], st, color="0.5", lw=1.0, alpha=0.8,
                        zorder=5)
    for inp, lab in zip(inputs, labels):
        tufac = read_tufac(inp)
        ts, ek, eth, ep, et, source = energy_series(inp, tufac)
        if ts is None:
            print(f"evrard: no snapshots/.log for '{inp}'", file=sys.stderr)
            continue
        style = "-" if source == "log" else "o-"
        pre = f"{lab} " if len(inputs) > 1 else ""
        line, = ax.plot(ts, et, style, ms=3, label=f"{pre}total")
        c = line.get_color()
        ax.plot(ts, ek, "--", lw=1.2, color=c, alpha=0.8, label=f"{pre}kinetic")
        ax.plot(ts, eth, ":", lw=1.6, color=c, alpha=0.8, label=f"{pre}thermal")
        ax.plot(ts, ep, "-.", lw=1.2, color=c, alpha=0.8, label=f"{pre}potential")
        drift = (et[-1] - et[0]) / abs(et[0]) if et[0] else float("nan")
        print(f"evrard[{lab}]: E_tot {et[0]:.4g}->{et[-1]:.4g} "
              f"(drift {100*drift:+.2f}%)  from {source}")
    ax.axhline(0.0, color="0.7", lw=0.8)
    ax.set_xlabel("t")
    ax.set_ylabel("energy")
    ax.set_title("Evrard collapse energy budget vs time")
    fig.tight_layout()
    finish_figure_with_legend(fig, ax, save=(f"{save}_energy.png" if save else None))


# --------------------------------------------------------------------------- #
#  Density field-map render (optional)
# --------------------------------------------------------------------------- #

def do_render(inputs, labels, t, args, save):
    """Render the density field at the comparison time."""
    entries = list(zip(labels, inputs))
    spec = RenderSpec(1, 2, 7, r"$\rho$", True, "slice", None, "inferno", False)
    times = None if t is None else [t]
    extent = args.render_extent if args.render_extent else [-1.0, 1.0, -1.0, 1.0]
    render_panels(entries, spec, times=times, n=args.render_n,
                  res=args.render_res, save=(f"{save}_rho" if save else None),
                  backend=args.render_backend, project=args.render_project,
                  slice_frac=args.render_slice_frac, extent=extent,
                  aspect=args.render_aspect)


def _add_args(parser):
    parser.add_argument("--refdir", type=str, default=DEFAULT_REFDIR,
                        help="directory with the reference CSVs "
                             "(evrardrho/evrardvel/evrardentrop.csv)")


def _plot(inputs, labels, params, args):
    ref = load_reference(args.refdir)
    reg_ref = getattr(args, "ref_table", None)
    print(f"[evrard] gamma={GAMMA:.4g}  M={M:g} R={R:g} G={G:g}  "
          f"comparison t={args.time:g}  reference: {sorted(ref) or 'none'}")
    plot_profiles(inputs, labels, args.time, ref, save=args.save, reg_ref=reg_ref)
    plot_energy_vs_time(inputs, labels, save=args.save, ref_table=reg_ref)
    if args.render:
        do_render(inputs, labels, args.time, args, save=args.save)


def main():
    standalone_cli(
        "evrard", PARAM_SPEC, _plot,
        reference_table=lambda inp, params, args: reference_table(inp, args.time),
        compare=lambda inp, params, args: compare_to_baseline(inp, args.time,
                                                              args.reference),
        description="Evrard collapse analysis (vs high-res 1D PPM, "
                    "Mandel et al. 2023).",
        time_default=DEFAULT_TIME,
        time_help=f"comparison time (default {DEFAULT_TIME:g})",
        add_arguments=_add_args)


if __name__ == "__main__":
    main()
