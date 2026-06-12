#!/usr/bin/env python3
"""
Analysis for the Gresho-Chan vortex test (IC_setup_gresho).

The Gresho vortex is a STEADY-STATE rotating flow: the centrifugal force is
balanced by a radial pressure gradient, so the azimuthal velocity profile is
time-independent. With rcyl = sqrt(x^2 + y^2) the analytic azimuthal velocity is

    v_phi(r) = 5 r            (r <= 0.2)
             = 2 - 5 r        (0.2 < r <= 0.4)
             = 0              (r > 0.4)

(IC_setup_gresho.getveli). The simulated azimuthal velocity is recovered from the
Cartesian (vx, vy) by v_phi = (x vy - y vx) / rcyl (the inverse of the cyl->cart
transform in IC_setup_gresho.transform_coordinate). Any deviation from v_phi(r)
is numerical diffusion of the vortex.

Diagnostics (mirroring IC_analysis_alfven):
  default       compare one (or more) runs at --times (default 1 and 3, one
                panel per time): v_phi vs rcyl, overlaid with the analytic
                profile; prints the (volume-weighted) L1 error per time.
  --convergence treat the positional inputs as the SAME test at different
                resolutions (n_x = the '<test><nx>' number in the run-folder
                name, e.g. 'gresho32...' -> 32), and plot L1(v_phi) vs n_x on
                log-log axes against the n_x^-2 line (ideal 2nd-order
                convergence), at the first of --times.
  energy        the kinetic energy E_kin(t) is always plotted vs time against its
                analytic (constant) value -- it should be conserved by the steady
                vortex; numerical dissipation shows up as a decay.
  --render      velocity streamlines (framework streamline_panels, colored by
                |v|), a density field map, and a particle slice (scatter colored
                by density) at the --render-times.

Usage:
    python IC_analysis_gresho.py <run-dir> [--times 1 3] [--save out] [--render]
    python IC_analysis_gresho.py <dir_N32> <dir_N64> <dir_N128> --convergence
"""

import os
import sys

import numpy as np
import readtipsy as tip

from IC_analysis_framework import (
    Param, standalone_cli,
    RenderSpec, render_panels, StreamSpec, streamline_panels,
    parse_resolution, select_at_time, series_from_log_or_snapshots,
    binned_profile, reg_nbins,
)
from IC_analysis_general import (
    find_files, loaddata,
    setup_rcparams, finish_figure_with_legend,
)

# Vortex constants, fixed by IC_setup_gresho (not -a..-h parameters).
RIN = 0.2                       # inner radius: v_phi peaks (=1) here
ROUT = 0.4                      # outer radius: v_phi returns to 0
VPEAK = 1.0                     # peak azimuthal velocity (5*RIN)

DEFAULT_TIMES = (1.0, 3.0)      # comparison times (~1 and ~3 turnovers at RIN)

# -a matches runtest.sh setup_gresho() (runtest.sh:331): the flow Mach number.
# It sets the background pressure but NOT the analytic v_phi(r), so it is carried
# only for labelling / bookkeeping; the velocity comparison is mach-independent.
PARAM_SPEC = [
    Param("a", "mach", float(np.sqrt(3.0 / 25.0)), float,
          "flow Mach number (sets pressure; v_phi profile is mach-independent)"),
]


def analytic_vphi(r):
    """Analytic azimuthal velocity v_phi(r) of the Gresho vortex (time-independent)."""
    r = np.asarray(r, dtype=float)
    return np.where(r <= RIN, 5.0 * r,
                    np.where(r <= ROUT, 2.0 - 5.0 * r, 0.0))


def profile(fn):
    """Per-particle (rcyl, v_phi, volume) and the snapshot time for run file fn.

    v_phi = (x vy - y vx) / rcyl  inverts the IC's cyl->cart velocity transform;
    the per-particle volume V_i = m_i / rho_i weights the spatial L1 norm."""
    tgdata, _td, _ts, _hdr, time, _N, _ngas, _nd, _ns, _h = loaddata(fn)
    x, y = tgdata[:, 1], tgdata[:, 2]
    vx, vy = tgdata[:, 4], tgdata[:, 5]
    rcyl = np.sqrt(x * x + y * y)
    with np.errstate(invalid="ignore", divide="ignore"):
        v_phi = np.where(rcyl > 0, (x * vy - y * vx) / rcyl, 0.0)
    vol = tgdata[:, 0] / tgdata[:, 7]               # m / rho (== particle volume)
    return rcyl, v_phi, vol, float(np.asarray(time).ravel()[0])


def l1_error(rcyl, v_phi, vol):
    """Volume-weighted L1 error of v_phi: sum V_i |v_phi - v_phi_analytic| / sum V_i.

    Volume-weighting (V_i = m_i/rho_i) makes this a proper spatial average,
    robust to non-uniform particle masses (the varying-mass IC option)."""
    err = np.abs(v_phi - analytic_vphi(rcyl))
    return float(np.sum(vol * err) / np.sum(vol))


# --------------------------------------------------------------------------- #
#  Kinetic-energy diagnostic vs time
# --------------------------------------------------------------------------- #

def kinetic_energy(fn):
    """(time, E_kin) for snapshot fn, with E_kin = sum 1/2 m_i |v_i|^2.

    Analytically constant for the steady Gresho vortex; a decay measures the
    numerical dissipation of the rotation."""
    tgdata, _td, _ts, _hdr, time = tip.readtipsy(fn)
    t = float(np.asarray(time).ravel()[0])
    m = tgdata[:, 0]
    v2 = np.sum(tgdata[:, 4:7] ** 2, axis=1)
    return t, float(np.sum(0.5 * m * v2))


def analytic_kinetic_energy(fn):
    """Analytic E_kin on a run's particle set: sum 1/2 m_i v_phi_analytic(r_i)^2.

    Evaluated on the actual particles (positions/masses) so it matches the run's
    thin-slab geometry and resolution exactly -- the constant the measured E_kin
    is conserved against."""
    tgdata, _td, _ts, _hdr, _time = tip.readtipsy(fn)
    rcyl = np.sqrt(tgdata[:, 1] ** 2 + tgdata[:, 2] ** 2)
    return float(np.sum(0.5 * tgdata[:, 0] * analytic_vphi(rcyl) ** 2))


def kinetic_energy_series(inp):
    """(t, E_kin, source) for run `inp`.

    Prefers the dense gasoline '.log' Ekin column (logged every timestep),
    falling back to per-snapshot computation. Returns (None, None, None) if
    neither source yields data."""
    return series_from_log_or_snapshots(inp, ["dTime", "Ekin"], kinetic_energy)


def plot_energy_vs_time(inputs, labels, save=None, ref_table=None):
    """E_kin(t) for each run against its analytic (constant) value.

    Uses the dense '.log' Ekin series when available (see kinetic_energy_series).
    The analytic constant is always evaluated on the run's particle set. When a
    reference is supplied its blessed (raw) E_kin(t) is overlaid as a dashed line."""
    import matplotlib.pyplot as plt
    setup_rcparams()
    fig, ax = plt.subplots(figsize=(8, 5))
    if ref_table is not None and "Ekin_raw" in ref_table:
        ax.plot(ref_table["te_raw"], ref_table["Ekin_raw"], "--", color="0.35",
                lw=1.5, zorder=5, label="reference")
    for inp, lab in zip(inputs, labels):
        ts, ek, source = kinetic_energy_series(inp)
        if ts is None:
            print(f"gresho: no snapshots/.log for '{inp}'", file=sys.stderr)
            continue
        style = "-" if source == "log" else "o-"   # dense log -> line
        line, = ax.plot(ts, ek, style, ms=3, label=lab)
        # Per-run analytic constant (depends on the thin-slab volume / resolution);
        # needs a snapshot's particle set, so skip it if only the .log is present.
        try:
            files = find_files(inp)
        except FileNotFoundError:
            files = []
        if files:
            e_an = analytic_kinetic_energy(files[0])
            ax.axhline(e_an, color=line.get_color(), ls="--", lw=1.0, alpha=0.7)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$E_{kin}=\sum\frac{1}{2} m_i|v_i|^2$")
    ax.set_title("Gresho vortex kinetic energy vs time "
                 "(dashed = analytic constant)")
    fig.tight_layout()
    finish_figure_with_legend(fig, ax,
                              save=(f"{save}_energy.png" if save else None))


# --------------------------------------------------------------------------- #
#  Default mode: azimuthal-velocity profile at --time vs analytic
# --------------------------------------------------------------------------- #

def plot_profiles(inputs, labels, times, save=None, ref_table=None):
    """v_phi vs r_cyl scatter + analytic, ONE PANEL PER TIME (default t=1, 3).

    The per-run L1 is annotated at the panel's top-LEFT corner (the legend --
    drawn on the axes in interactive mode -- prefers the data-free top-right of
    this profile, which used to cover the annotation there). When a reference is
    supplied its blessed binned v_phi(r) curve is overlaid (the steady vortex is
    time-independent, so one blessed profile serves every panel)."""
    import matplotlib.pyplot as plt
    setup_rcparams()
    times = list(times)
    fig, axes = plt.subplots(1, len(times), figsize=(7 * len(times), 6),
                             sharey=True, squeeze=False)
    axes = axes[0]

    abox = dict(boxstyle="round", fc="white", ec="0.7", alpha=0.8)
    has_ref = ref_table is not None and "v_phi" in ref_table
    for pi, (ax, t) in enumerate(zip(axes, times)):
        rmax = ROUT
        l1_txt = []
        for inp, lab in zip(inputs, labels):
            sel = select_at_time(inp, t)
            if sel is None:
                continue
            fn, tt = sel
            rcyl, v_phi, vol, _time = profile(fn)
            rmax = max(rmax, float(np.percentile(rcyl, 99.5)))
            ax.scatter(rcyl, v_phi, s=3, alpha=0.4, rasterized=True,
                       label=(lab if pi == 0 else None))
            l1 = l1_error(rcyl, v_phi, vol)
            print(f"gresho[{lab}]: t={tt:.4g}  L1(v_phi)={l1:.6g}")
            l1_txt.append(f"{lab}: {l1:.3g}")

        # Reference overlay: the blessed binned v_phi(r) (dashed grey).
        if has_ref:
            ax.plot(ref_table["r_prof"], ref_table["v_phi"], "--", color="0.35",
                    lw=1.6, zorder=5, label=("reference" if pi == 0 else None))
            if "L1_vphi" in ref_table:
                l1_txt.append(f"ref: {float(ref_table['L1_vphi'][0]):.3g}")

        # Analytic overlay (time-independent).
        rg = np.linspace(0.0, rmax, 400)
        ax.plot(rg, analytic_vphi(rg), "k-", lw=1.8, zorder=6,
                label=("analytic" if pi == 0 else None))

        if l1_txt:
            ax.text(0.02, 0.97, "L1(v_phi)\n" + "\n".join(l1_txt),
                    transform=ax.transAxes, fontsize=9, va="top", ha="left",
                    bbox=abox)
        ax.set_xlabel(r"$r_{cyl}=\sqrt{x^2+y^2}$")
        ax.set_title(f"t = {t:g}")
    axes[0].set_ylabel(r"$v_\phi$")
    fig.suptitle("Gresho vortex azimuthal velocity")
    fig.tight_layout()
    finish_figure_with_legend(fig, list(axes),
                              save=(f"{save}_profile.png" if save else None))


# --------------------------------------------------------------------------- #
#  Renders (--render): velocity streamlines, density map, particle slice
# --------------------------------------------------------------------------- #

def do_render(inputs, labels, args, save):
    """Three --render figures (the vortex is a thin z slab, x-y plane):
      * velocity streamlines (framework streamline_panels, colored by |v| --
        the vortex ring and any numerical break-up of it);
      * density field map ("average" LOS projection; rho=1 everywhere
        analytically, so structure is discretisation noise);
      * particle slice (per-particle scatter colored by density through the
        z mid-plane -- shows the particle distribution itself).
    All forward the generic --render-* flags."""
    entries = list(zip(labels, inputs))
    streamline_panels(entries, StreamSpec(1, 2, 4, 5, "velocity"),
                      times=args.render_times, n=args.render_n,
                      res=min(args.render_res or 64, 128),
                      save=(f"{save}_stream_velocity" if save else None),
                      extent=args.render_extent,
                      title="Gresho vortex velocity streamlines")
    rho_spec = RenderSpec(1, 2, 7, "density", log=False, project="average",
                          cmap="viridis")
    render_panels(entries, rho_spec, times=args.render_times,
                  n=args.render_n, res=args.render_res,
                  save=(f"{save}_rho" if save else None),
                  backend=args.render_backend, project=args.render_project,
                  slice_frac=args.render_slice_frac,
                  extent=args.render_extent, aspect=args.render_aspect,
                  nsmooth=args.render_nsmooth)
    render_panels(entries, rho_spec, times=args.render_times,
                  n=args.render_n, res=args.render_res,
                  save=(f"{save}_particles" if save else None),
                  backend="particles", project="slice",
                  slice_frac=args.render_slice_frac,
                  extent=args.render_extent, aspect=args.render_aspect,
                  nsmooth=args.render_nsmooth)


# --------------------------------------------------------------------------- #
#  Convergence mode: L1(v_phi) vs resolution, against n_x^-2
# --------------------------------------------------------------------------- #

def plot_convergence(inputs, t, save=None, ref_table=None):
    import matplotlib.pyplot as plt
    nxs, l1s = [], []
    for inp in inputs:
        nx = parse_resolution(inp)
        if nx is None:
            print(f"gresho: cannot read '<test><nx>' resolution from '{inp}'; "
                  f"skipping (needed for --convergence).", file=sys.stderr)
            continue
        sel = select_at_time(inp, t)
        if sel is None:
            continue
        fn, tt = sel
        rcyl, v_phi, vol, _time = profile(fn)
        l1 = l1_error(rcyl, v_phi, vol)
        nxs.append(nx)
        l1s.append(l1)
        print(f"gresho: n_x={nx:<5d} t={tt:.4g}  L1(v_phi)={l1:.6g}")

    if len(nxs) < 2:
        print("gresho: need >=2 resolutions with a valid '<test><nx>' number "
              "for a convergence plot.", file=sys.stderr)
        return

    order = np.argsort(nxs)
    nxs = np.array(nxs)[order]
    l1s = np.array(l1s)[order]

    setup_rcparams()
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.loglog(nxs, l1s, "o-", label=r"$L_1(v_\phi)$")
    # Ideal 2nd-order reference, anchored at the coarsest resolution.
    ref = l1s[0] * (nxs / nxs[0]) ** (-2.0)
    ax.loglog(nxs, ref, "k--", label=r"$\propto n_x^{-2}$")
    # Blessed-baseline L1 point at its own resolution (if available).
    if (ref_table is not None and "n_x" in ref_table
            and "L1_vphi" in ref_table):
        ax.loglog([float(ref_table["n_x"][0])], [float(ref_table["L1_vphi"][0])],
                  "D", color="C3", ms=9, mfc="none", label="reference")
    ax.set_xlabel(r"$n_x$")
    ax.set_ylabel(r"$L_1(v_\phi)$")
    ax.set_title(f"Gresho vortex convergence (t={t:g})")
    fig.tight_layout()
    finish_figure_with_legend(fig, ax,
                              save=(f"{save}_convergence.png" if save else None))


# --------------------------------------------------------------------------- #
#  Regression baseline (bless / compare the normalised kinetic-energy series)
# --------------------------------------------------------------------------- #

def reference_table(inp, t=None):
    """Multi-block regression table feeding ALL of gresho's non-render plots:

      gate / energy   x        : time grid
                      Ekin     : E_kin(t)/E_kin(0)  (normalised, the gate key --
                                 resolution-independent dissipation probe)
                      te_raw   : time grid (the raw-axis energy plot's x)
                      Ekin_raw : E_kin(t)           (raw, for the energy overlay)
      profile         r_prof   : r_cyl bin centres
                      v_phi    : binned v_phi(r) at the comparison time
                      L1_vphi  : volume-weighted L1 vs analytic (1-elt array)
      convergence     n_x      : the run's '<test><nx>' resolution

    E_kin is analytically constant for the steady Gresho vortex, so the
    normalised history is a flat-at-1 dissipation probe (the resolution-
    independent CI gate). Blocks are included only when computable. None if no
    usable E_kin series AND no profile snapshot."""
    table = {}
    ts, ek, _src = kinetic_energy_series(inp)
    if ts is not None and np.isfinite(ek[0]) and ek[0] != 0:
        table.update({"x": ts, "Ekin": ek / ek[0], "te_raw": ts, "Ekin_raw": ek})
    # Profile block (steady vortex -> blessed at the comparison time t).
    sel = select_at_time(inp, t)
    if sel is not None:
        fn, _tt = sel
        rcyl, v_phi, vol, _time = profile(fn)
        cx, vb = binned_profile(rcyl, v_phi, nbins=reg_nbins(inp))
        if cx is not None:
            table["r_prof"] = cx
            table["v_phi"] = vb
            table["L1_vphi"] = np.array([l1_error(rcyl, v_phi, vol)])
        nx = parse_resolution(inp)
        if nx is not None:
            table["n_x"] = np.array([float(nx)])
    return table or None


def _add_args(parser):
    parser.add_argument("--convergence", action="store_true",
                        help="treat inputs as different resolutions and plot "
                             "L1(v_phi) vs n_x against the n_x^-2 line")


def _plot(inputs, labels, params, args):
    # args.times is the resolved comparison-time list (one profile panel each);
    # args.time is its first entry (used for the single-time convergence plot).
    # args.ref_table is the blessed baseline (None unless --reference), overlaid
    # on every figure.
    ref = getattr(args, "ref_table", None)
    if args.convergence:
        plot_convergence(inputs, args.time, save=args.save, ref_table=ref)
    else:
        plot_profiles(inputs, labels, args.times, save=args.save, ref_table=ref)
    plot_energy_vs_time(inputs, labels, save=args.save, ref_table=ref)
    if args.render:
        do_render(inputs, labels, args, args.save)


def main():
    standalone_cli(
        "gresho", PARAM_SPEC, _plot,
        reference_table=lambda inp, params, args: reference_table(inp, args.time),
        description="Gresho-Chan vortex analysis.",
        reg_keys=("Ekin",), times_default=list(DEFAULT_TIMES),
        add_arguments=_add_args)


if __name__ == "__main__":
    main()
