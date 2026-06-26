#!/usr/bin/env python3
"""
Analysis for the Sedov-Taylor blast wave test (IC_setup_sedov).

A point/region of energy E is deposited at the origin of a uniform medium
(rho0, gamma=5/3); the resulting spherical blast follows the self-similar
Sedov-Taylor solution. With spherical radius r = sqrt(x^2+y^2+z^2) this script
compares three radial profiles at a chosen time against that analytic solution
and reports the (volume-weighted) L1 error of each:

    density     rho(r)
    radial vel  v_r(r) = (x vx + y vy + z vz) / r
    entropy     A(r)   = p / rho^gamma   (p = (gamma-1) rho u)

It also renders 2-D maps (mid-plane slice) of density, thermal pressure,
magnetic pressure (|B|^2/2) and kinetic-energy density (rho|v|^2/2).

The analytic Sedov-Taylor solution is computed by integrating the self-similar
ODEs (see _sedov_tables); it is the STRONG-EXPLOSION solution and so is overlaid
only for the non-magnetic, low-ambient-pressure cases. IC_setup_sedov has several
configurations (runtest.sh -a select, -b betain):
  * betain != 0  -> a magnetic field is present; there is NO analytic solution,
                    so only the profiles/renders are shown (no overlay, no L1).
  * select 0     -> finite ambient pressure Pout; the strong-shock analytic is
                    approximate there (a caveat is printed).

The blast energy E and the ambient state are read from the FIRST snapshot (E is
conserved), so no hand-tuned constants are needed.

Usage:
    python IC_analysis_sedov.py <run-dir> [--time T] [--save out] [--render]
    python IC_analysis_sedov.py <run-dir> -b 0 --save out        # hydro + analytic
"""

import os
import sys

import numpy as np
import readtipsy as tip
from scipy.integrate import solve_ivp

from IC_analysis_framework import (
    Param, RenderSpec, render_panels, aux_divB_error, spec_slug,
    bless_render_specs, render_ref_overlay_path, movie_render_specs,
    select_at_time, binned_profile, standalone_cli,
    thermal_pressure, magnetic_pressure, ke_density,
    series_from_log_or_snapshots, emit_warning,
)
from IC_analysis_general import (
    find_files, loaddata, read_vector_aux, read_tufac,
    setup_rcparams, finish_figure, finish_figure_with_legend,
)

GAMMA = 5.0 / 3.0               # IC_setup_sedov.gamma
NU = 3                          # spherical blast
DEFAULT_TIME = 0.02             # default comparison time (before the shock reaches
                                # the periodic box edge ~0.5)
ENT_RMIN = 0.1                  # inner radius for the entropy comparison: A=p/rho^g
                                # diverges in the near-vacuum core, so the log-A
                                # overlay and its L1 are taken over r >= ENT_RMIN

# runtest.sh setup_sedov() parameter order (runtest.sh:383-411). Only `select`
# and `betain` change the analysis (magnetic vs hydro); the rest are carried for
# CLI/label parity with runtest.
PARAM_SPEC = [
    Param("a", "select", 1, int, "energy-injection mode (0 Pin/Rin, 1 point, 2 kernel)"),
    Param("b", "betain", 0.0, float, "blast plasma-beta (0 -> no B, analytic valid)"),
    Param("c", "ublast", 10.0, float, "blast energy (select 1/2)"),
    Param("d", "Rin", 0.125, float, "inner blast radius (select 0)"),
    Param("e", "Pin", 100.0, float, "inner pressure (select 0)"),
    Param("f", "Pout", 1.0, float, "ambient pressure (select 0)"),
]


# --------------------------------------------------------------------------- #
#  Analytic Sedov-Taylor self-similar solution
# --------------------------------------------------------------------------- #
# Reduced variables in xi = r/R(t) (R the shock radius, D = dR/dt = (2/5)R/t):
#   v   = D V(xi),   rho = rho0 G(xi),   c^2 = gamma p/rho = D^2 Z(xi).
# Integrating the Euler equations in self-similar form gives the ODEs in
# _sedov_tables; validated against the strong-shock jumps (G=4, V=3/4, Z=5/16
# for gamma=5/3) and the dimensionless constant xi0=1.1517 (gamma=5/3 spherical).

_SEDOV_CACHE = {}


def _sedov_tables(gamma=GAMMA, xi_min=1e-6, n=4000):
    """Integrate the self-similar Sedov ODEs from the shock (xi=1) inward.

    Returns (xi, V, G, Z, xi0, I_E) sorted in increasing xi, where
    R = xi0 (E t^2 / rho0)^(1/5) is the shock radius and I_E the energy integral.
    Cached per gamma."""
    key = round(gamma, 10)
    if key in _SEDOV_CACHE:
        return _SEDOV_CACHE[key]
    g = gamma
    V1, G1, Z1 = 2.0/(g+1.0), (g+1.0)/(g-1.0), 2.0*g*(g-1.0)/(g+1.0)**2

    def rhs(xi, y):
        V, G, Z = y
        s = V - xi
        gg = (3.0*V + 4.0*V*s/xi - 6.0*Z/(g*s)) / (2.0*Z - 2.0*s*s)  # G'/G
        return [-2.0*V/xi - s*gg,                # V'
                G*gg,                            # G'
                Z*(3.0/s - (1.0-g)*gg)]          # Z'

    xi_eval = np.linspace(1.0, xi_min, n)
    sol = solve_ivp(rhs, (1.0, xi_min), [V1, G1, Z1], t_eval=xi_eval,
                    rtol=1e-9, atol=1e-12, max_step=0.01)
    xi, V, G, Z = sol.t[::-1], sol.y[0][::-1], sol.y[1][::-1], sol.y[2][::-1]
    # Energy integral E = 4 pi R^3 D^2 rho0 I_E, with D=(2/5)R/t:
    #   I_E = int_0^1 [1/2 G V^2 + G Z/(gamma(gamma-1))] xi^2 dxi
    integ = (0.5*G*V**2 + G*Z/(g*(g-1.0))) * xi**2
    I_E = float(np.trapz(integ, xi))
    xi0 = (25.0/(16.0*np.pi*I_E))**0.2
    out = (xi, V, G, Z, xi0, I_E)
    _SEDOV_CACHE[key] = out
    return out


def sedov_profiles(r, t, E, rho0, p_amb, gamma=GAMMA):
    """Analytic (rho, v_r, p) at spherical radii `r` for blast energy E at time t.

    Inside the shock (r<R) interpolates the self-similar tables; outside it is the
    undisturbed ambient (rho0, 0, p_amb). Returns (rho, v, p, R)."""
    xi, V, G, Z, xi0, _I = _sedov_tables(gamma)
    R = xi0 * (E * t**2 / rho0) ** 0.2
    D = 0.4 * R / t
    lam = np.clip(np.asarray(r, dtype=float) / R, 0.0, 1.0)
    Vi = np.interp(lam, xi, V)
    Gi = np.interp(lam, xi, G)
    Zi = np.interp(lam, xi, Z)
    rho = rho0 * Gi
    v = D * Vi
    p = rho0 * Gi * D**2 * Zi / gamma
    outside = np.asarray(r, dtype=float) > R
    rho = np.where(outside, rho0, rho)
    v = np.where(outside, 0.0, v)
    p = np.where(outside, p_amb, p)
    return rho, v, p, R


def energy_fractions(gamma=GAMMA):
    """(f_kin, f_therm): the kinetic / thermal fractions of the Sedov energy.

    The blast is self-similar, so these fractions are CONSTANT in time (depend
    only on gamma): E_kin = f_kin E, E_therm = f_therm E, f_kin+f_therm=1.
      I_kin   = int_0^1 (1/2) G V^2 xi^2 dxi
      I_therm = int_0^1 G Z/(gamma(gamma-1)) xi^2 dxi   (E_int = p/(gamma-1))
    and f = I/I_E. For gamma=5/3 spherical: ~0.283 kinetic, ~0.717 thermal."""
    xi, V, G, Z, _xi0, I_E = _sedov_tables(gamma)
    I_kin = float(np.trapz(0.5 * G * V**2 * xi**2, xi))
    I_therm = float(np.trapz(G * Z / (gamma * (gamma - 1.0)) * xi**2, xi))
    return I_kin / I_E, I_therm / I_E


# --------------------------------------------------------------------------- #
#  Per-snapshot fields and blast energy
# --------------------------------------------------------------------------- #

def profile(fn, tufac):
    """Per-particle (r, rho, v_r, u, entropy, volume) and the snapshot time.

    r is spherical radius from the origin (the blast is centred there by
    IC_setup_sedov.centermybox); v_r is the outward radial velocity; u = tufac*T
    is the specific internal energy (col 8 is temperature); entropy is the adiabat
    A = p/rho^gamma = (gamma-1) u rho^(1-gamma)."""
    tgdata, _td, _ts, _hdr, time, _N, _ng, _nd, _ns, _h = loaddata(fn)
    pos = tgdata[:, 1:4]
    vel = tgdata[:, 4:7]
    rho = tgdata[:, 7]
    u = tufac * tgdata[:, 8]                      # T (col 8) -> internal energy
    r = np.sqrt(np.sum(pos**2, axis=1))
    with np.errstate(invalid="ignore", divide="ignore"):
        v_r = np.where(r > 0, np.sum(pos*vel, axis=1) / r, 0.0)
    entropy = (GAMMA - 1.0) * u * rho ** (1.0 - GAMMA)
    vol = tgdata[:, 0] / rho                      # m/rho = particle volume
    return (r, rho, v_r, u, entropy, vol,
            float(np.asarray(time).ravel()[0]))


def ambient_state(fn, tufac):
    """(rho0, p_amb): the undisturbed background density and pressure (medians)."""
    tgdata, _td, _ts, _hdr, _time = tip.readtipsy(fn)
    rho0 = float(np.median(tgdata[:, 7]))
    u_amb = tufac * float(np.median(tgdata[:, 8]))
    return rho0, (GAMMA - 1.0) * rho0 * u_amb


def snapshot_energies(fn, tufac):
    """(time, E_kin, E_therm, E_mag) for snapshot fn.

    E_kin = sum 1/2 m v^2,  E_therm = sum m u = tufac sum m T,
    E_mag = sum 1/2 |B|^2 V_i (V_i = m/rho; 0 if no BField aux)."""
    tgdata, _td, _ts, hdr, time = tip.readtipsy(fn)
    t = float(np.asarray(time).ravel()[0])
    m = tgdata[:, 0]
    e_kin = float(np.sum(0.5 * m * np.sum(tgdata[:, 4:7]**2, axis=1)))
    e_therm = float(tufac * np.sum(m * tgdata[:, 8]))
    B = read_vector_aux(fn, "BField", int(hdr[2]))
    e_mag = 0.0
    if B is not None:
        b2 = np.sum(np.asarray(B, dtype=float)**2, axis=1)
        e_mag = float(np.sum(0.5 * b2 * m / tgdata[:, 7]))
    return t, e_kin, e_therm, e_mag


def energy_series(inp, tufac):
    """(t, E_kin, E_therm, E_mag, source) for run `inp`.

    Prefers the dense gasoline '.log' Ekin/Eth/Emag columns (logged every
    timestep -- far denser than the snapshot dumps), falling back to
    per-snapshot summation (E_mag from the BField aux there). Returns
    (None,)*5 if neither source has data."""
    return series_from_log_or_snapshots(
        inp, ["dTime", "Ekin", "Eth", "Emag"],
        lambda fn: snapshot_energies(fn, tufac),
        valid=lambda lg: np.any((lg["Ekin"] + lg["Eth"]) != 0))


def analytic_E(params, fn, tufac):
    """(E, rho0, p_amb) for the analytic solution.

    For point/kernel injection (select 1/2) the Sedov energy is the input
    `ublast` (-c): it drives the shock, whereas the snapshot's summed energy can
    differ from numerical losses. For select 0 it is the snapshot total in excess
    of ambient."""
    rho0, p_amb = ambient_state(fn, tufac)
    if int(params["select"]) in (1, 2):
        return float(params["ublast"]), rho0, p_amb
    tgdata, *_ = tip.readtipsy(fn)
    _t, e_kin, e_therm, _e = snapshot_energies(fn, tufac)
    u_amb = tufac * float(np.median(tgdata[:, 8]))
    return e_kin + e_therm - u_amb * float(np.sum(tgdata[:, 0])), rho0, p_amb


def l1_error(r, field, analytic, vol, rmax):
    """Volume-weighted L1 error of `field` vs `analytic`, over particles r<rmax."""
    sel = r < rmax
    if not np.any(sel):
        return float("nan")
    err = np.abs(field[sel] - analytic[sel])
    return float(np.sum(vol[sel]*err) / np.sum(vol[sel]))




# --------------------------------------------------------------------------- #
#  Profiles: rho, v_r, entropy vs r, with analytic overlay + L1
# --------------------------------------------------------------------------- #

def plot_profiles(inputs, labels, t, params, save=None, ref_table=None):
    import matplotlib.pyplot as plt
    setup_rcparams()
    fig, axarr = plt.subplots(3, 2, figsize=(14, 13), sharex=True)
    # 2 columns x 3 rows.
    axRho, axRhoLog = axarr[0]
    axVel, axP = axarr[1]
    axEnt, axU = axarr[2]
    axes = (axRho, axRhoLog, axVel, axP, axEnt, axU)
    axRhoLog.set_yscale("log")
    axEnt.set_yscale("log")                  # entropy spans many decades -> log A
    axU.set_yscale("log")                    # u = p/((g-1)rho) diverges at centre
    have_analytic = (params["betain"] == 0)  # magnetic case has no analytic solution
    rmax_box = 0.5                           # inscribed sphere of the [-.5,.5]^3 box

    # Reference overlay: the blessed binned profile of each panel (dashed grey).
    if ref_table is not None and "rho" in ref_table:
        xr = ref_table["x"]
        for ax, key in ((axRho, "rho"), (axRhoLog, "rho"), (axVel, "v_r"),
                        (axP, "P"), (axEnt, "entropy"), (axU, "u")):
            if key in ref_table:
                ax.plot(xr, ref_table[key], "--", color="red", lw=1.5, zorder=5,
                        label=("reference" if ax is axRho else None))

    l1_lines = []
    for inp, lab in zip(inputs, labels):
        sel = select_at_time(inp, t)
        if sel is None:
            continue
        fn, tt = sel
        tufac = read_tufac(inp)
        r, rho, v_r, u, entropy, vol, _time = profile(fn, tufac)
        P = (GAMMA - 1.0) * rho * u

        axRho.scatter(r, rho, s=3, alpha=0.4, rasterized=True, label=lab)
        axRhoLog.scatter(r, rho, s=3, alpha=0.4, rasterized=True, label=lab)
        axVel.scatter(r, v_r, s=3, alpha=0.4, rasterized=True, label=lab)
        axP.scatter(r, P, s=3, alpha=0.4, rasterized=True, label=lab)
        axEnt.scatter(r, entropy, s=3, alpha=0.4, rasterized=True, label=lab)
        axU.scatter(r, u, s=3, alpha=0.4, rasterized=True, label=lab)

        if have_analytic:
            E, rho0, p_amb = analytic_E(params, find_files(inp)[0], tufac)
            rho_a, v_a, p_a, R = sedov_profiles(r, tt, E, rho0, p_amb)
            ent_a = p_a / rho_a ** GAMMA
            u_a = p_a / ((GAMMA - 1.0) * rho_a)
            l1r = l1_error(r, rho, rho_a, vol, rmax_box)
            l1v = l1_error(r, v_r, v_a, vol, rmax_box)
            l1p = l1_error(r, P, p_a, vol, rmax_box)
            # Entropy A and internal energy u both diverge in the near-vacuum core
            # (rho->0); compare log10 over r >= ENT_RMIN.
            shell = r >= ENT_RMIN
            la_s = np.log10(np.clip(entropy, 1e-30, None))
            la_a = np.log10(np.clip(ent_a, 1e-30, None))
            l1e = l1_error(r[shell], la_s[shell], la_a[shell], vol[shell], rmax_box)
            lu_s = np.log10(np.clip(u, 1e-30, None))
            lu_a = np.log10(np.clip(u_a, 1e-30, None))
            l1u = l1_error(r[shell], lu_s[shell], lu_a[shell], vol[shell], rmax_box)
            print(f"sedov[{lab}]: t={tt:.4g}  E={E:.4g}  R_shock={R:.4g}  "
                  f"L1(rho)={l1r:.4g}  L1(v_r)={l1v:.4g}  L1(P)={l1p:.4g}  "
                  f"L1(logA)={l1e:.4g}  L1(logu)={l1u:.4g}")
            l1_lines.append(f"{lab}: rho {l1r:.3g}, v {l1v:.3g}, P {l1p:.3g}, "
                            f"logA {l1e:.3g}, logu {l1u:.3g}")
        else:
            print(f"sedov[{lab}]: t={tt:.4g}  (magnetic / no analytic: profiles only)")

    # Single analytic overlay (smooth curve) using the first input's E/ambient.
    if have_analytic and inputs:
        tufac = read_tufac(inputs[0])
        E, rho0, p_amb = analytic_E(params, find_files(inputs[0])[0], tufac)
        tt = select_at_time(inputs[0], t)[1]
        rg = np.linspace(1e-4, rmax_box, 800)
        rho_a, v_a, p_a, R = sedov_profiles(rg, tt, E, rho0, p_amb)
        ent_a = p_a / rho_a ** GAMMA
        u_ag = p_a / ((GAMMA - 1.0) * rho_a)
        entg = rg >= ENT_RMIN                       # A and u shown over r >= ENT_RMIN
        axRho.plot(rg, rho_a, "k-", lw=1.8, label="analytic")
        axRhoLog.plot(rg, rho_a, "k-", lw=1.8, label="analytic")
        axVel.plot(rg, v_a, "k-", lw=1.8, label="analytic")
        axP.plot(rg, p_a, "k-", lw=1.8, label="analytic")
        axEnt.plot(rg[entg], ent_a[entg], "k-", lw=1.8, label="analytic")
        axU.plot(rg[entg], u_ag[entg], "k-", lw=1.8, label="analytic")
        for ax in axes:
            ax.axvline(R, color="0.5", ls=":", lw=1.2)
        es = ent_a[entg]
        if es.size:
            axEnt.set_ylim(max(es.min() * 0.3, 1e-5), es.max() * 3.0)
        us = u_ag[entg]
        if us.size:
            axU.set_ylim(max(us.min() * 0.3, 1e-5), us.max() * 3.0)

    if l1_lines:
        abox = dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85)
        axRho.text(0.98, 0.96, "\n".join(l1_lines), transform=axRho.transAxes,
                   fontsize=8, va="top", ha="right", bbox=abox)

    axRho.set_ylabel(r"$\rho$")
    axRhoLog.set_ylabel(r"$\rho$ (log)")
    axRhoLog.set_ylim(1e-4, 10)
    axVel.set_ylabel(r"$v_r$")
    axP.set_ylabel(r"$P$")
    axEnt.set_ylabel(r"$A=p/\rho^{\gamma}$ (log)")
    axU.set_ylabel(r"$u=p/((\gamma-1)\rho)$ (log)")
    rlabel = r"$r=\sqrt{x^2+y^2+z^2}$"
    axEnt.set_xlabel(rlabel)            # bottom-left cell
    axU.set_xlabel(rlabel)             # bottom-right cell
    axEnt.set_xlim(0.0, rmax_box)            # shared x: limit r to the box inscribed sphere
    suffix = "" if have_analytic else "  (magnetic: no analytic)"
    fig.suptitle(f"Sedov blast profiles (t={t:g}){suffix}")
    fig.tight_layout()
    finish_figure_with_legend(fig, axRho,
                              save=(f"{save}_profile.png" if save else None))


# --------------------------------------------------------------------------- #
#  Total energy vs time (conservation; compare against SPLASH / input ublast)
# --------------------------------------------------------------------------- #

def plot_energy_vs_time(inputs, labels, t, params, save=None, ref_table=None,
                        warn_tol=None):
    """Energy vs time for each run.

    For the hydro point/kernel cases (select 1/2, betain=0) each energy is plotted
    NORMALIZED by its analytic value -- E/E_analytic, E_kin/E_kin,analytic,
    E_therm/E_therm,analytic (analytic total=ublast, kinetic=f_kin ublast,
    thermal=f_therm ublast) -- so a perfectly conserved blast with the Sedov
    partition sits on 1. The per-component relative difference at the comparison
    time `t` is also printed. Magnetic runs (no analytic) plot the raw energies."""
    import matplotlib.pyplot as plt
    setup_rcparams()
    fig, ax = plt.subplots(figsize=(8, 5))
    hydro_blast = (int(params["select"]) in (1, 2) and params["betain"] == 0)
    E0 = float(params["ublast"])
    f_kin, f_therm = energy_fractions() if hydro_blast else (None, None)
    magnetic = params["betain"] != 0
    # Reference overlay: the blessed (normalised) total/kinetic/thermal energy.
    if ref_table is not None and "Etot_n" in ref_table:
        te = ref_table["t_energy"]
        ax.plot(te, ref_table["Etot_n"], "--", color="red", lw=1.3, zorder=5,
                label="reference")
        if "Ekin_n" in ref_table:
            ax.plot(te, ref_table["Ekin_n"], "--", color="red", lw=0.9,
                    alpha=0.7, zorder=5)
        if "Eth_n" in ref_table:
            ax.plot(te, ref_table["Eth_n"], ":", color="red", lw=1.1, zorder=5)
    for inp, lab in zip(inputs, labels):
        tufac = read_tufac(inp)
        ts, ek, eth, em, source = energy_series(inp, tufac)
        if ts is None:
            print(f"sedov: no snapshots/.log for '{inp}'", file=sys.stderr)
            continue
        et = ek + eth + (em if magnetic else 0.0)   # total incl. E_mag when MHD
        if hydro_blast:                         # normalize each by its analytic
            yt, yk, yth = et / E0, ek / (f_kin * E0), eth / (f_therm * E0)
            yem = None
        elif magnetic:
            # No analytic -> normalize EVERY component by the initial TOTAL energy
            # so they share a common, physically-meaningful scale: the total
            # stays ~1 (conservation check) and each component reads as its
            # fraction of the budget (E_kin starts at 0, E_mag small, E_th ~1).
            et0 = et[0] if (np.isfinite(et[0]) and et[0] != 0) else 1.0
            yt, yk, yth, yem = et / et0, ek / et0, eth / et0, em / et0
        else:
            yt, yk, yth, yem = et, ek, eth, None
        style = "-" if source == "log" else "o-"   # dense log -> line
        line, = ax.plot(ts, yt, style, ms=3, label=f"{lab} total")
        c = line.get_color()
        ax.plot(ts, yk, "--", lw=1.2, color=c, alpha=0.8, label="kinetic")
        ax.plot(ts, yth, ":", lw=1.4, color=c, alpha=0.8, label="thermal")
        if magnetic:                            # E_mag / E_tot(0)
            ax.plot(ts, yem, "-.", lw=1.4, color=c, alpha=0.8, label="magnetic")
            print(f"sedov[{lab}] energy (/ E_tot(0), t={ts[0]:.4g}): "
                  f"total {yt[0]:.4g}->{yt[-1]:.4g}  "
                  f"kinetic {yk[0]:.4g}->{yk[-1]:.4g}  "
                  f"thermal {yth[0]:.4g}->{yth[-1]:.4g}  "
                  f"magnetic {yem[0]:.4g}->{yem[-1]:.4g}")

        # Relative difference of each energy vs analytic at the comparison time.
        if hydro_blast:
            j = int(np.argmin(np.abs(ts - t)))
            refs = (("total", et[j], E0),
                    ("kinetic", ek[j], f_kin * E0),
                    ("thermal", eth[j], f_therm * E0))
            print(f"sedov[{lab}] energy vs analytic at t={ts[j]:.4g}:")
            for name, sim, ana in refs:
                rel = 100.0 * (sim - ana) / ana if ana else float("nan")
                print(f"   {name:<7s} sim={sim:8.4g}  analytic={ana:8.4g}  "
                      f"rel.diff={rel:+6.2f}%")
                # Soft energy-conservation warning on the TOTAL energy only --
                # informational, never changes the exit code (the --reg-tol gate
                # on the density profile is the hard pass/fail).
                if (name == "total" and warn_tol is not None and ana
                        and abs((sim - ana) / ana) > warn_tol):
                    emit_warning("Sedov energy conservation",
                                 f"{lab}: total energy {rel:+.2f}% vs analytic at "
                                 f"t={ts[j]:.4g} exceeds the {100 * warn_tol:g}% "
                                 f"soft tolerance")

    if hydro_blast:
        ax.axhline(1.0, color="k", ls="--", lw=1.2, label="analytic (=1)")
        ax.set_ylabel(r"$E / E_{\rm analytic}$")
        ax.set_title("Sedov energy / analytic vs time (total / kinetic / thermal)")
    elif magnetic:
        ax.axhline(1.0, color="k", ls="--", lw=1.2, label="initial total (=1)")
        ax.set_ylabel(r"$E / E_{\rm tot}(0)$")
        ax.set_title("Sedov energy / initial total vs time "
                     "(total / kinetic / thermal / magnetic)")
    else:
        if int(params["select"]) in (1, 2):
            ax.axhline(E0, color="k", ls="--", lw=1.3,
                       label=f"ublast = {E0:g} (total)")
        ax.set_ylabel("energy")
        ax.set_title("Sedov energy vs time (total / kinetic / thermal)")
    ax.set_xlabel("t")
    fig.tight_layout()
    finish_figure_with_legend(fig, ax, save=(f"{save}_energy.png" if save else None))


# --------------------------------------------------------------------------- #
#  Field-map renders: density, pressure, magnetic pressure, KE density
# --------------------------------------------------------------------------- #

# Fixed color limits for `-a 2` (uniform-field MHD blast): the same limits as
# Stone et al. (2008) and Price et al. (2018), for a DIRECT visual comparison
# with their figures (linear scales there, so log=False for these specs).
STONE_CLIMS = {
    "rho":     (0.19, 2.98),
    "Ekin":    (0.0, 33.1),
    "P":       (1.0, 42.4),
    "Pmag":    (25.2, 65.9),
}


def _fixed_clim(lo, hi):
    """clim callback returning a fixed range for field-valued projections and
    None (auto) for integrated 'column' renders."""
    def clim(params, proj):
        return None if proj == "column" else (lo, hi)
    return clim


def render_specs(betain, select):
    """The unified render set, identical across the -a cases: rho, KE density,
    P, and -- only when betain>0 (no field otherwise) -- Pmag and divBerr (LOG,
    fixed 1e-3..1.0). For -a 2 the rho/Ekin/P/Pmag maps carry the Stone et al. (2008)
    / Price et al. (2018) fixed limits (linear scale) for direct comparison."""
    stone = int(select) == 2

    def spec(q, label, slug, cmap):
        if stone and slug in STONE_CLIMS:
            return RenderSpec(1, 2, q, label, False, "slice",
                              _fixed_clim(*STONE_CLIMS[slug]), cmap, False,
                              slug=slug)
        return RenderSpec(1, 2, q, label, True, "slice", None, cmap, False,
                          slug=slug)

    specs = [
        (spec(7, r"$\rho$", "rho", "inferno"), "rho"),
        (spec(ke_density(), r"$\frac{1}{2}\rho v^2$", "Ekin", "viridis"), "Ekin"),
        (spec(thermal_pressure(GAMMA), r"$P$", "P", "inferno"), "P"),
    ]
    if betain != 0:
        specs.append((spec(magnetic_pressure(), r"$\frac{1}{2}|B|^2$", "Pmag",
                           "magma"), "Pmag"))
        # divBerr is only meaningful when a field is present (betain != 0).
        specs.append((RenderSpec(1, 2, aux_divB_error("DivB", "BField"),
                                 "divBerr", True, "slice",
                                 _fixed_clim(1e-3, 1.0), "inferno", False,
                                 slug="divBerr"),
                      "divBerr"))
    return specs


IC_BOX = [-0.5, 0.5, -0.5, 0.5]    # default render domain (the [-0.5,0.5]^2 box)


def do_renders(inputs, labels, t, params, args, save, reference=None):
    """Render the field maps, forwarding the --render-* CLI flags.

    The CLI overrides win: --render-backend, --render-project (None -> each
    RenderSpec's "slice" default), --render-slice-frac, --render-res, --render-
    aspect, and --render-extent (defaults to the IC box [-0.5,0.5]^2 when not
    given). For -a 1 the density map is ALSO drawn with the particles backend
    (per-particle scatter -- the particle plot). When `reference` is given and
    there is a single input, the blessed grid renders are overlaid (sim |
    reference | difference); the scatter particle plot is never overlaid."""
    entries = list(zip(labels, inputs))
    times = None if t is None else [t]
    extent = args.render_extent if args.render_extent else IC_BOX
    betain, select = params["betain"], int(params["select"])
    overlay = reference is not None and len(inputs) == 1
    common = dict(times=times, n=4, res=args.render_res,
                  project=args.render_project,
                  slice_frac=args.render_slice_frac, extent=extent,
                  aspect=args.render_aspect, nsmooth=args.render_nsmooth)
    for sp, slug in render_specs(betain, select):
        ref_path = (render_ref_overlay_path(inputs[0], reference, slug)
                    if overlay else None)
        render_panels(entries, sp, save=(f"{save}_{slug}" if save else None),
                      backend=args.render_backend, render_reference=ref_path,
                      residual=getattr(args, "residual", False),
                      **common)
        if select == 1 and slug == "rho":      # -a 1: particle plot of density
            render_panels(entries, sp,
                          save=(f"{save}_{slug}_particles" if save else None),
                          backend="particles", **common)


def _save_extra(args, params):
    """Bless the grid render baselines (at the comparison time) for
    --save-reference, on the IC box so the overlay aligns."""
    betain, select = params["betain"], int(params["select"])
    extent = args.render_extent if args.render_extent else IC_BOX
    specs = [sp for sp, _slug in render_specs(betain, select)]
    bless_render_specs(args.inputs[0], specs, args, "sedov",
                       times=(None if args.time is None else [args.time]),
                       extent=extent)


# --------------------------------------------------------------------------- #
#  MHD radial profiles (betain != 0): |B| vs r and divBerr vs r
# --------------------------------------------------------------------------- #

def plot_mhd_profiles(inputs, labels, t, save=None):
    """Two-panel radial profile at the comparison time for magnetic runs:
    field strength |B|(r) and the normalised divergence error
    divBerr = h |div B| / |B| (r) -- the field-amplification and divergence-
    cleanliness diagnostics (no analytic; scatter + binned-median line)."""
    import matplotlib.pyplot as plt
    setup_rcparams()
    fig, (axb, axd) = plt.subplots(1, 2, figsize=(14, 6))
    divb_of = aux_divB_error("DivB", "BField")
    any_data = False
    for ci, (inp, lab) in enumerate(zip(inputs, labels)):
        sel = select_at_time(inp, t)
        if sel is None:
            continue
        fn, tt = sel
        tgdata, _td, _ts, hdr, _tm = tip.readtipsy(fn)
        ngas = int(hdr[2])
        B = read_vector_aux(fn, "BField", ngas)
        if B is None:
            print(f"sedov[{lab}]: no BField aux; skipping MHD profiles",
                  file=sys.stderr)
            continue
        r = np.sqrt(np.sum(tgdata[:, 1:4]**2, axis=1))
        bmag = np.sqrt(np.sum(np.asarray(B, dtype=float)**2, axis=1))
        derr = divb_of(fn, tgdata, ngas)
        for ax, q in ((axb, bmag), (axd, derr)):
            ax.scatter(r, q, s=3, alpha=0.25, color=f"C{ci}", rasterized=True,
                       label=(f"{lab} (t={tt:.3g})" if ax is axb else None))
            cx, qb = binned_profile(r, q, percentile=(0.5, 99.5))
            if cx is not None:
                ax.plot(cx, qb, "-", lw=1.8, color=f"C{ci}")
        any_data = True
    if not any_data:
        plt.close(fig)
        return
    axb.set_xlabel("r")
    axb.set_ylabel(r"$|B|$")
    axd.set_xlabel("r")
    axd.set_ylabel(r"$h\,|\nabla\cdot B|/|B|$")
    axd.set_yscale("log")
    fig.suptitle(f"Sedov MHD radial profiles (t={t:g})")
    fig.tight_layout()
    finish_figure_with_legend(fig, [axb, axd],
                              save=(f"{save}_mhdprofile.png" if save else None))


# --------------------------------------------------------------------------- #
#  Regression baseline (bless / compare the binned density profile at --time)
# --------------------------------------------------------------------------- #

REG_NBINS = 64


def reference_table(inp, t, params):
    """Multi-block regression table at the comparison time t:
      profile  x        : r-grid
               rho/v_r/P/entropy/u : binned radial profiles (one per panel)
      energy   t_energy : time grid
               Etot_n/Ekin_n/Eth_n : energies / analytic (hydro_blast only)

    The radial density profile of the blast is the headline Sedov diagnostic and
    the gate key. None if no snapshot / too few bins."""
    sel = select_at_time(inp, t)
    if sel is None:
        return None
    fn, _tt = sel
    tufac = read_tufac(inp)
    r, rho, v_r, u, entropy, _vol, _time = profile(fn, tufac)
    P = (GAMMA - 1.0) * rho * u
    cx, qb = binned_profile(r, {"rho": rho, "v_r": v_r, "P": P,
                                "entropy": entropy, "u": u},
                            nbins=REG_NBINS, percentile=(0.5, 99.5))
    if cx is None:
        return None
    table = {"x": cx, **qb}
    # Energy block (own time grid), normalised by the analytic partition for the
    # hydro blast (matching the energy plot's y-axis); raw otherwise.
    hydro_blast = (int(params["select"]) in (1, 2) and params["betain"] == 0)
    if hydro_blast:
        ts, ek, eth, em, _src = energy_series(inp, tufac)
        if ts is not None:
            E0 = float(params["ublast"])
            f_kin, f_therm = energy_fractions()
            table["t_energy"] = ts
            table["Etot_n"] = (ek + eth) / E0
            table["Ekin_n"] = ek / (f_kin * E0)
            table["Eth_n"] = eth / (f_therm * E0)
    return table


def _plot(inputs, labels, params, args):
    betain = params["betain"]
    if betain == 0:
        if params["select"] == 0 and params["Pout"] != 0:
            print("sedov: select 0 has finite ambient pressure; the strong-shock "
                  "Sedov analytic is approximate.", file=sys.stderr)
    else:
        print("sedov: betain != 0 (magnetic): no analytic solution; plotting "
              "profiles and renders only.", file=sys.stderr)
    ref = getattr(args, "ref_table", None)
    plot_profiles(inputs, labels, args.time, params, save=args.save, ref_table=ref)
    plot_energy_vs_time(inputs, labels, args.time, params, save=args.save,
                        ref_table=ref, warn_tol=getattr(args, "energy_warn_tol",
                                                         None))
    if betain != 0:                 # field amplification + divergence cleanliness
        plot_mhd_profiles(inputs, labels, args.time, save=args.save)
    if args.render:
        do_renders(inputs, labels, args.time, params, args, save=args.save,
                   reference=args.reference)
    if args.movie:
        betain, select = params["betain"], int(params["select"])
        extent = args.render_extent if args.render_extent else IC_BOX
        specs = [sp for sp, _slug in render_specs(betain, select)]
        movie_render_specs(inputs[0], labels[0], specs, args, args.save, "sedov",
                           extent=extent)


def _add_args(parser):
    parser.add_argument(
        "--energy-warn-tol", type=float, default=0.05,
        help="soft gate: emit a WARNING (and, under GitHub Actions, a yellow "
             "::warning:: annotation) when |total energy / analytic - 1| at the "
             "comparison time exceeds this fraction (default 0.05 = 5%%). This is "
             "informational only -- it NEVER changes pass/fail (the --reg-tol "
             "density gate does). Set to a large value to silence it.")


def main():
    standalone_cli(
        "sedov", PARAM_SPEC, _plot,
        reference_table=lambda inp, params, args: reference_table(inp, args.time,
                                                                  params),
        description="Sedov-Taylor blast analysis.",
        add_arguments=_add_args,
        reg_keys=("rho",), time_default=DEFAULT_TIME,
        time_help=f"comparison time (default {DEFAULT_TIME:g}; before the shock "
                  f"reaches the box edge)",
        save_extra=_save_extra)


if __name__ == "__main__":
    main()
