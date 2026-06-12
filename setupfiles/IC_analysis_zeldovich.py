#!/usr/bin/env python3
"""
Analysis for the Zeldovich pancake (IC_setup_zeldovich).

A single sinusoidal density perturbation in an expanding (Einstein-de Sitter)
universe collapses into a "pancake" caustic at redshift z_c = 2. The Zeldovich
approximation is the EXACT nonlinear solution until shell crossing: with the
Lagrangian coordinate q and growing-mode amplitude b(a) = (1+z_c) a (= (1+z_c)/
(1+z); a = snapshot scale factor),

    x(q)        = q - b sin(k q) / k                 (Eulerian position)
    rho(q)/rho0 = 1 / (1 - b cos(k q))               (mass conservation dq/dx)
    v(q)        = -H0 (1+z_c) sqrt(a) sin(k q) / k    (IC velocity convention)

with k = 2 pi / L, L = 1 (the comoving box [-0.5, 0.5] in code units). The map is
single-valued (invertible) while b < 1, i.e. a < a_c = 1/(1+z_c) = 1/3; at the
caustic the density diverges. The simulation snapshot time IS the scale factor a
(comoving run), so the analytic is evaluated at each snapshot's a.

Diagnostics (standalone, sedov-style: reuses the framework parser/loaders, not
run_cli):
  default       density rho/rho0 (and velocity vx) vs x at --time (the scale
                factor a; default 0.3, a well-developed pre-caustic pancake),
                overlaid with the parametric Zeldovich solution; prints the
                volume-weighted L1 density error.
  peak          rho_peak/rho0 vs a is always plotted against the analytic
                1/(1 - (1+z_c) a) -> infinity at the caustic a_c = 1/3.

The DENSITY comparison is exact; the velocity overlay follows the IC's own
(comoving, v ~ sqrt(a)) convention and is convention-dependent -- treat it as a
guide, not a rigorous check.

Usage:
    python IC_analysis_zeldovich.py <run-dir> [--time 0.3] [--save out]
"""

import os
import sys

import numpy as np

from IC_analysis_framework import (
    Param, standalone_cli,
    select_at_time, binned_profile,
)
from IC_analysis_general import (
    find_files, loaddata, setup_rcparams, finish_figure_with_legend,
)

# Constants fixed by IC_setup_zeldovich (not -a/-b parameters).
ZC = 2.0                       # collapse redshift; caustic at a_c = 1/(1+ZC) = 1/3
L = 1.0                        # comoving box length (code units, [-0.5, 0.5])
K = 2.0 * np.pi / L            # perturbation wavenumber
H0 = 2.894405                  # Hubble constant in code units (IC velocity)
RHO0 = 1.0                     # mean comoving density (rhozero)
A_CAUSTIC = 1.0 / (1.0 + ZC)   # = 1/3

DEFAULT_TIME = 0.3             # scale factor a of a well-developed, pre-caustic pancake

# -a/-b match runtest.sh setup_zeldovich(); neither enters the code-unit density
# analytic (k and z_c are fixed). Lboxinkpc is the length unit; zi sets the
# initial scale factor a_i = 1/(1+zi) (carried for context / the velocity).
PARAM_SPEC = [
    Param("a", "Lboxinkpc", 1000.0, float, "box size in kpc (length unit; density analytic is in code units)"),
    Param("b", "zi",        100.0,  float, "initial redshift (a_i = 1/(1+zi); must exceed z_c=2)"),
]


def growth_b(a):
    """Zeldovich growing-mode amplitude b(a) = (1+z_c) a (EdS, D ~ a)."""
    return (1.0 + ZC) * a


def analytic_profile(a, nq=4000):
    """Parametric Zeldovich solution at scale factor a, pre-caustic.

    Returns (x(q), rho/rho0, v) on a Lagrangian q-grid. x(q) is monotonic for
    b<1, so it can be inverted (np.interp) to evaluate the analytic at particle
    positions. Returns (None, None, None) if a >= a_c (caustic; map multivalued)."""
    b = growth_b(a)
    if b >= 1.0:
        return None, None, None
    q = np.linspace(-L / 2.0, L / 2.0, nq)
    x = q - b * np.sin(K * q) / K
    rho = 1.0 / (1.0 - b * np.cos(K * q))
    v = -H0 * (1.0 + ZC) * np.sqrt(a) * np.sin(K * q) / K
    return x, rho, v


def profile(fn):
    """Per-particle (x, rho/rho0, vx, volume) and the scale factor a for run fn.

    rho is normalized by the volume-weighted mean rho0 = sum(m)/sum(m/rho) (=
    total mass / total volume), so the comparison is unit-free."""
    tgdata, _td, _ts, _hdr, time, _N, _ngas, _nd, _ns, _h = loaddata(fn)
    a = float(np.asarray(time).ravel()[0])
    x = tgdata[:, 1]
    vx = tgdata[:, 4]
    rho = tgdata[:, 7]
    vol = tgdata[:, 0] / rho                       # m / rho == particle volume
    rho_bar = float(np.sum(tgdata[:, 0]) / np.sum(vol))
    return x, rho / rho_bar, vx, vol, a


def l1_density(x, rho_norm, vol, a):
    """Volume-weighted L1 of rho/rho0 vs the Zeldovich analytic (np.interp on the
    monotonic x(q) map). None if a >= a_c (no single-valued analytic)."""
    xa, rhoa, _va = analytic_profile(a)
    if xa is None:
        return None
    rho_pred = np.interp(x, xa, rhoa)              # xa monotonic increasing (b<1)
    err = np.abs(rho_norm - rho_pred)
    return float(np.sum(vol * err) / np.sum(vol))


# --------------------------------------------------------------------------- #
#  Default mode: density (+velocity) profile at --time vs analytic
# --------------------------------------------------------------------------- #

def plot_profiles(inputs, labels, t, save=None, ref_table=None):
    import matplotlib.pyplot as plt
    setup_rcparams()
    fig, (axr, axv) = plt.subplots(1, 2, figsize=(13, 6))

    # Reference overlay: the blessed binned rho/rho0(x) and v_x(x) (dashed grey).
    if ref_table is not None and "rho" in ref_table:
        axr.plot(ref_table["x"], ref_table["rho"], "--", color="0.35", lw=1.6,
                 zorder=5, label="reference")
        if "vx" in ref_table:
            axv.plot(ref_table["x"], ref_table["vx"], "--", color="0.35", lw=1.6,
                     zorder=5, label="reference")

    l1_txt = []
    for inp, lab in zip(inputs, labels):
        sel = select_at_time(inp, t)
        if sel is None:
            continue
        fn, a = sel
        x, rho_norm, vx, vol, _a = profile(fn)
        axr.scatter(x, rho_norm, s=3, alpha=0.4, rasterized=True,
                    label=f"{lab} (a={a:.3g})")
        axv.scatter(x, vx, s=3, alpha=0.4, rasterized=True,
                    label=f"{lab} (a={a:.3g})")
        l1 = l1_density(x, rho_norm, vol, a)
        if l1 is not None:
            print(f"zeldovich[{lab}]: a={a:.4g} (z={1/a-1:.4g})  L1(rho/rho0)={l1:.6g}")
            l1_txt.append(f"{lab}: {l1:.3g}")
        else:
            print(f"zeldovich[{lab}]: a={a:.4g} >= a_caustic={A_CAUSTIC:.4g}; "
                  f"past shell crossing -- no single-valued analytic.",
                  file=sys.stderr)

    # Analytic overlay at the requested time (pre-caustic only).
    xa, rhoa, va = analytic_profile(t)
    if xa is not None:
        axr.plot(xa, rhoa, "k-", lw=1.8, label="Zeldovich")
        axv.plot(xa, va, "k-", lw=1.8, label="Zeldovich (IC conv.)")
    else:
        print(f"zeldovich: requested a={t:.4g} is past the caustic a_c={A_CAUSTIC:.4g};"
              f" analytic overlay omitted.", file=sys.stderr)

    if l1_txt:
        abox = dict(boxstyle="round", fc="white", ec="0.7", alpha=0.8)
        axr.text(0.98, 0.96, "L1(rho/rho0)  " + ",  ".join(l1_txt),
                 transform=axr.transAxes, fontsize=9, va="top", ha="right", bbox=abox)

    axr.set_xlabel(r"$x$"); axr.set_ylabel(r"$\rho/\rho_0$"); axr.set_yscale("log")
    axv.set_xlabel(r"$x$"); axv.set_ylabel(r"$v_x$")
    fig.suptitle(f"Zeldovich pancake (a={t:g}, z={1/t-1:.3g}; "
                 rf"$a_c$={A_CAUSTIC:.3g})")
    fig.tight_layout()
    finish_figure_with_legend(fig, axr,
                              save=(f"{save}_profile.png" if save else None))


# --------------------------------------------------------------------------- #
#  Always: peak density vs scale factor (caustic formation)
# --------------------------------------------------------------------------- #

def plot_peak_vs_time(inputs, labels, save=None, ref_table=None):
    """rho_peak/rho0 (99.5th percentile) vs a against 1/(1 - (1+z_c) a)."""
    import matplotlib.pyplot as plt
    setup_rcparams()
    fig, ax = plt.subplots(figsize=(8, 6))
    if ref_table is not None and "rho_peak" in ref_table:
        ax.plot(ref_table["a_peak"], ref_table["rho_peak"], "--", color="0.35",
                lw=1.5, zorder=5, label="reference")
    a_all = []
    for inp, lab in zip(inputs, labels):
        files = find_files(inp)
        if not files:
            print(f"zeldovich: no snapshots for '{inp}'", file=sys.stderr)
            continue
        a_s, peak_s = [], []
        for fn in files:
            x, rho_norm, vx, vol, a = profile(fn)
            a_s.append(a)
            peak_s.append(float(np.percentile(rho_norm, 99.5)))
        o = np.argsort(a_s)
        a_s = np.array(a_s)[o]; peak_s = np.array(peak_s)[o]
        a_all.extend(a_s.tolist())
        ax.plot(a_s, peak_s, "o-", ms=3, label=lab)
        print(f"zeldovich[{lab}]: rho_peak/rho0 at a_max={a_s[-1]:.4g} = {peak_s[-1]:.4g}")

    # Analytic peak 1/(1-b), b=(1+z_c)a, over the pre-caustic range present.
    if a_all:
        a_grid = np.linspace(min(a_all), min(max(a_all), 0.995 * A_CAUSTIC), 300)
        ax.plot(a_grid, 1.0 / (1.0 - growth_b(a_grid)), "k-", lw=1.8,
                label=r"$1/(1-(1+z_c)a)$")
    ax.axvline(A_CAUSTIC, color="0.6", ls=":", lw=1.0,
               label=rf"$a_c={A_CAUSTIC:.3g}$")
    ax.set_xlabel(r"$a$ (scale factor)")
    ax.set_ylabel(r"$\rho_{peak}/\rho_0$")
    ax.set_yscale("log")
    ax.set_title("Zeldovich pancake: peak density vs scale factor")
    fig.tight_layout()
    finish_figure_with_legend(fig, ax,
                              save=(f"{save}_peak.png" if save else None))


# --------------------------------------------------------------------------- #
#  Regression baseline (bless / compare the binned density profile at --time)
# --------------------------------------------------------------------------- #

REG_NBINS = 64


def reference_table(inp, t):
    """Regression table at scale factor t feeding both non-render plots:
      profile  x        : x-grid
               rho      : binned rho/rho0(x)   (headline diagnostic + gate key)
               vx       : binned v_x(x)        (the v_x profile panel)
      peak     a_peak   : scale-factor grid
               rho_peak : rho_peak/rho0 vs a   (the caustic-formation plot)
    None if no snapshot / too few bins."""
    sel = select_at_time(inp, t)
    if sel is None:
        return None
    fn, _a = sel
    x, rho_n, vx, _vol, _aa = profile(fn)
    cx, qb = binned_profile(x, {"rho": rho_n, "vx": vx},
                            nbins=REG_NBINS, percentile=(0.5, 99.5))
    if cx is None:
        return None
    table = {"x": cx, "rho": qb["rho"], "vx": qb["vx"]}
    # Peak-density-vs-a block (its own scale-factor grid).
    files = find_files(inp)
    if files:
        a_s, peak_s = [], []
        for f in files:
            xx, rn, _v, _vol, a = profile(f)
            a_s.append(a)
            peak_s.append(float(np.percentile(rn, 99.5)))
        o = np.argsort(a_s)
        table["a_peak"] = np.array(a_s)[o]
        table["rho_peak"] = np.array(peak_s)[o]
    return table


def _plot(inputs, labels, params, args):
    ref = getattr(args, "ref_table", None)
    plot_profiles(inputs, labels, args.time, save=args.save, ref_table=ref)
    plot_peak_vs_time(inputs, labels, save=args.save, ref_table=ref)


def main():
    standalone_cli(
        "zeldovich", PARAM_SPEC, _plot,
        reference_table=lambda inp, params, args: reference_table(inp, args.time),
        description="Zeldovich pancake analysis.",
        reg_keys=("rho",), time_default=DEFAULT_TIME,
        time_help=f"scale factor a to compare profiles at (default "
                  f"{DEFAULT_TIME:g}; must be < a_c={A_CAUSTIC:.3g})")


if __name__ == "__main__":
    main()
