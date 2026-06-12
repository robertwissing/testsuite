#!/usr/bin/env python3
"""
Analysis for the cosmowave test (Berlok 2022, MNRAS 515 3492 + this work's
divB-cleaning extension): linear (M)HD waves in an expanding Einstein-de Sitter
universe. The "time" axis is the SCALE FACTOR a (the Tipsy snapshot time in a
comoving run); wave amplitudes evolve as power laws of a, so the diagnostic is
the perturbation amplitude vs ln(a/a_i) overlaid with linear theory -- NOT a
metric-vs-time / convergence test.

Per snapshot the amplitude of each perturbed quantity (dv, dB, drho) is the
Fourier projection onto sin/cos(k(x-x0)) (`amp_mag`); these are compared to the
Berlok 2022 closed-form envelopes scaled by the analytic a-powers
((a/a_i)^{1/4} for dB_c/B and drho/rho, (a/a_i)^{-3/4} for du/V), under which an
exact linear mode is flat.

Cases (-a, must match IC_setup_cosmowave.create's `case`):
    0 - standing Alfven                (paper Sec. 5.1.1, eqs. 75-76)
    1 - traveling Alfven, n=2          (paper Sec. 5.1.2, eqs. 79-80)
    2 - standing compressible          (paper Sec. 5.2.1/2/3, eqs. 81-82)
    3 - traveling compressible g=4/3   (paper Sec. 5.2.1, eqs. 83-84; n=3)
    4 - divergence-cleaning test       (this work)

Standalone like alfven/gresho/sedov/noh/wave: reuses the framework's
Param/build_parser/params_from_args/find_files/loaddata, not run_cli (the x-axis
is the scale factor, not the per-snapshot time, and there is no single metric).
Several inputs overlay as separate marker series (one shared analytic curve).

Usage:
    python IC_analysis_cosmowave.py <dir> [-a 3] [-b 1] [-c 0] [--save out]
    python IC_analysis_cosmowave.py <dir32> <dir64> -a 4 --ch 1.0 --save out
"""

import os
import sys

import numpy as np

from IC_analysis_framework import (
    Param, build_parser, params_from_args,
    reference_path_for, save_reference, find_reference, residuals,
)
from IC_analysis_general import (
    loaddata, find_files, try_aux, read_bfield,
    amp_mag, setup_rcparams, finish_figure_with_legend,
)

# Fixed IC constants (match IC_setup_cosmowave.create): amplitude, Hubble, the
# wavenumber k=2pi (wavelength = 2*dx = 1), phase origin xdot0 = -dx, omega_S.
A_U = 1e-1
H0 = 2.894405
KWAVE = 2.0 * np.pi
X0 = -0.5
OMEGA_S = np.pi

# runtest.sh setup_cosmowave() parameter order (runtest.sh:setup_cosmowave).
PARAM_SPEC = [
    Param("a", "case", 3, int,
          "0 standing Alfven, 1 traveling Alfven, 2 standing compr., "
          "3 traveling compr. g=4/3, 4 div-clean EdS"),
    Param("b", "factor_A", 1.0, float,
          "omega_A = factor_A*pi (0 = MHD off, 1 = paper default)"),
    Param("c", "factor_G", 0.0, float,
          "omega_G = factor_G*pi/2 (0 = self-gravity off, 1 = paper default)"),
]


# --------------------------------------------------------------------------- #
#  Per-case linear-theory constants (mirrors IC_setup_cosmowave.create)
# --------------------------------------------------------------------------- #

def case_params(case, factor_A=1.0, factor_G=0.0):
    """The per-case setup constants computed exactly as IC_setup_cosmowave."""
    p = dict(A_u=A_U, H0=H0, k=KWAVE, x0=X0, omega_S=OMEGA_S)
    p['case'] = int(case)
    p['omega_A'] = factor_A * np.pi
    p['omega_G'] = factor_G * np.pi / 2.0
    p['gamma'] = 4. / 3.
    p['alfven'] = int(case) in (0, 1)
    sigma2 = p['omega_A']**2 + p['omega_S']**2 - p['omega_G']**2
    p['kappa'] = np.sqrt(max(sigma2 - 1./16., 1e-30))
    p['kappa_A'] = np.sqrt(max(p['omega_A']**2 - 1./16., 1e-30))
    if int(case) in (0, 2, 4):
        p['a_i'] = 1. / 128.
    elif int(case) == 1:
        p['a_i'] = np.exp(-2 * 2 * np.pi / p['kappa_A'])   # n=2
    else:                                                  # case 3
        p['a_i'] = np.exp(-2 * 3 * np.pi / p['kappa'])     # n=3
    p['V_A'] = (p['H0'] / p['k']) * p['omega_A']
    p['V_S'] = (p['H0'] / p['k']) * p['omega_S']
    p['V_G'] = (p['H0'] / p['k']) * p['omega_G']
    p['rhozero'] = (((p['k'] * p['V_G'])**2) / (4 * np.pi)
                    if p['omega_G'] > 0 else 1.0)
    p['Bzero_i'] = p['V_A'] * np.sqrt(p['rhozero']) * p['a_i']**(-0.5)
    return p


# --------------------------------------------------------------------------- #
#  Linear-theory amplitude predictions (Berlok 2022), all comoving
# --------------------------------------------------------------------------- #

def _alfven_envelope(a, p):
    """Standing-Alfven modal envelope (paper eq. 75-76), psi=kappa_A ln(a/a_i)."""
    A_u, ai, kA = p['A_u'], p['a_i'], p['kappa_A']
    psi = kA * np.log(a / ai)
    dB_over_B = A_u * (a / ai)**0.25 * np.abs(np.sin(psi)) / kA
    du_over_VA = A_u * (a / ai)**(-0.75) * np.abs(np.cos(psi) + np.sin(psi) / (4.0 * kA))
    return dB_over_B, du_over_VA


def _alfven_traveling(a, p):
    """Traveling-Alfven eigenmode (paper eq. 79-80): flat under the a-scalings."""
    A_u, ai = p['A_u'], p['a_i']
    return A_u * (a / ai)**0.25, A_u * (a / ai)**(-0.75)


def _compressible_envelope(a, p, traveling):
    """Compressible-mode envelope/eigenmode (paper eq. 81-84)."""
    A_u, ai, kp = p['A_u'], p['a_i'], p['kappa']
    psi = kp * np.log(a / ai)
    if traveling:
        env_rho = env_v = env_B = A_u
    else:
        env_rho = A_u * np.abs(np.sin(psi)) / kp
        env_v = A_u * np.abs(np.cos(psi) + np.sin(psi) / (4.0 * kp))
        env_B = env_rho
    return (env_rho * (a / ai)**0.25,
            env_v * (a / ai)**(-0.75),
            env_B * (a / ai)**0.25)


def _divclean_hyperbolic(a, p, c_h=None):
    """Hyperbolic-only closed form of the divB-cleaning test (this work):
        max|div B_c|(a) = (a/a_i)^2 * A_u B0_c k * cos[2 Omega_h(sqrt a - sqrt a_i)]
    with Omega_h = k c_h / H_0. Without c_h only the (a/a_i)^2 envelope is known."""
    ai = p['a_i']
    envelope = (a / ai)**2
    if c_h is None:
        return envelope, None
    Omega_h = p['k'] * c_h / p['H0']
    return envelope, np.cos(2.0 * Omega_h * (np.sqrt(a) - np.sqrt(ai)))


# --------------------------------------------------------------------------- #
#  Per-run amplitude scan vs scale factor
# --------------------------------------------------------------------------- #

def scan_input(input_arg, p):
    """Amplitude table for one run, sorted by scale factor a.

    Columns: a, dv_amp, dB_amp, drho_amp, dBx, dBy, dBz, divB_l2, divB_max,
    psi_max. dv/dB are the projections onto the mode each case excites."""
    files = find_files(input_arg)
    print(f"[scan] {len(files)} snapshot(s) in {input_arg}")
    rows = []
    for fn in files:
        try:
            tgdata, _td, _ts, _hdr, time, _N, ngas, _nd, _ns, _h = loaddata(fn)
        except Exception as e:
            print(f"  skip {fn}: {e}", file=sys.stderr)
            continue
        a = float(np.asarray(time).ravel()[0])
        if a <= 0:
            a = p['a_i']
        m, x = tgdata[:, 0], tgdata[:, 1]
        vx, vy = tgdata[:, 4], tgdata[:, 5]
        rho = tgdata[:, 7]
        rho_mean = float(np.average(rho, weights=m))
        v_pert = vy if p['alfven'] else vx        # Alfven dv_y, compressible dv_x
        dv_amp = amp_mag(m, x, v_pert, p['k'], p['x0'])
        drho_amp = amp_mag(m, x, rho - rho_mean, p['k'], p['x0'])

        B = read_bfield(fn, ngas)
        if B is not None:
            dBx = amp_mag(m, x, B[:, 0] - np.average(B[:, 0], weights=m), p['k'], p['x0'])
            dBy = amp_mag(m, x, B[:, 1] - np.average(B[:, 1], weights=m), p['k'], p['x0'])
            dBz = amp_mag(m, x, B[:, 2] - np.average(B[:, 2], weights=m), p['k'], p['x0'])
            dB_amp = {0: dBy, 1: dBy, 2: dBz, 3: dBz, 4: dBx}[p['case']]
        else:
            dBx = dBy = dBz = dB_amp = np.nan

        divB = try_aux(fn, "DivB", ngas)           # gasoline aux is 'DivB'
        divB_l2 = float(np.sqrt(np.mean(divB**2))) if divB is not None else np.nan
        divB_max = float(np.max(np.abs(divB))) if divB is not None else np.nan
        psi = try_aux(fn, "BClean", ngas)
        psi_max = float(np.max(np.abs(psi))) if psi is not None else np.nan

        rows.append((a, dv_amp, dB_amp, drho_amp, dBx, dBy, dBz,
                     divB_l2, divB_max, psi_max))
        print(f"  a={a:.5g}  dv={dv_amp:.3e}  dB={dB_amp:.3e}  drho={drho_amp:.3e}"
              + (f"  divB_l2={divB_l2:.3e}" if not np.isnan(divB_l2) else ""))

    if not rows:
        print(f"cosmowave: no usable snapshots in {input_arg}", file=sys.stderr)
        return None
    out = np.array(rows, dtype=float)
    return out[np.argsort(out[:, 0])]


# --------------------------------------------------------------------------- #
#  Plots (amplitude vs ln(a/a_i)), legend in its own figure
# --------------------------------------------------------------------------- #

def _a_curve(results, p):
    a_all = np.concatenate([d[:, 0] for _, d in results])
    return np.geomspace(p['a_i'], a_all.max(), 400)


def plot_alfven(results, p, save=None, ref_table=None):
    import matplotlib.pyplot as plt
    setup_rcparams()
    ac = _a_curve(results, p)
    lac = np.log(ac / p['a_i'])
    if p['case'] == 0:
        dBth, duth = _alfven_envelope(ac, p)
        title = "Standing Alfven (Berlok 2022 Fig. 1)"
    else:
        dBth, duth = _alfven_traveling(ac, p)
        title = "Traveling Alfven (Berlok 2022 Fig. 2)"

    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    # Reference overlay: the blessed amplitudes rescaled to each panel's axes.
    if ref_table is not None and "dB" in ref_table:
        ar = ref_table["x"]; lnar = np.log(ar / p['a_i'])
        ax[0].plot(lnar, (ref_table["dB"] / p['Bzero_i']) * (ar / p['a_i'])**0.25,
                   "--", color="0.5", lw=1.4, zorder=5, label="reference")
        if "dv" in ref_table:
            ax[1].plot(lnar, (ref_table["dv"] / p['V_A']) * (ar / p['a_i'])**0.75,
                       "--", color="0.5", lw=1.4, zorder=5, label="reference")
    for lab, d in results:
        a, dv, dB = d[:, 0], d[:, 1], d[:, 2]
        lna = np.log(a / p['a_i'])
        ax[0].plot(lna, (dB / p['Bzero_i']) * (a / p['a_i'])**0.25,
                   marker='o', ls='none', label=lab)
        ax[1].plot(lna, (dv / p['V_A']) * (a / p['a_i'])**0.75,
                   marker='o', ls='none', label=lab)
    ax[0].plot(lac, dBth * (ac / p['a_i'])**0.25, 'k-', label='linear theory')
    ax[1].plot(lac, duth * (ac / p['a_i'])**0.75, 'k-', label='linear theory')
    ax[0].set_xlabel(r"$\ln(a/a_i)$")
    ax[0].set_ylabel(r"$(\delta B_c / B_{c,0}) \, (a/a_i)^{1/4}$")
    ax[1].set_xlabel(r"$\ln(a/a_i)$")
    ax[1].set_ylabel(r"$(\delta u / V_A) \, (a/a_i)^{3/4}$")
    fig.suptitle(title)
    fig.tight_layout()
    finish_figure_with_legend(fig, (ax[0], ax[1]),
                              save=(f"{save}_alfven.png" if save else None))


def plot_compressible(results, p, save=None, ref_table=None):
    import matplotlib.pyplot as plt
    setup_rcparams()
    ac = _a_curve(results, p)
    lac = np.log(ac / p['a_i'])
    traveling = (p['case'] == 3)
    drth, duth, dBth = _compressible_envelope(ac, p, traveling=traveling)
    title = ("Traveling compressible (Berlok 2022 Fig. 4)"
             if traveling else "Standing compressible (Berlok 2022 Fig. 3)")

    fig, ax = plt.subplots(1, 3, figsize=(16, 5))
    # Reference overlay: the blessed amplitudes rescaled to each panel's axes.
    if ref_table is not None and "drho" in ref_table:
        ar = ref_table["x"]; lnar = np.log(ar / p['a_i'])
        ax[0].plot(lnar, (ref_table["drho"] / p['rhozero']) * (ar / p['a_i'])**0.25,
                   "--", color="0.5", lw=1.4, zorder=5, label="reference")
        if "dv" in ref_table:
            ax[1].plot(lnar, (ref_table["dv"] / p['V_S']) * (ar / p['a_i'])**0.75,
                       "--", color="0.5", lw=1.4, zorder=5, label="reference")
        if "dB" in ref_table:
            ax[2].plot(lnar, (ref_table["dB"] / p['Bzero_i']) * (ar / p['a_i'])**0.25,
                       "--", color="0.5", lw=1.4, zorder=5, label="reference")
    for lab, d in results:
        a, dv, dB, drho = d[:, 0], d[:, 1], d[:, 2], d[:, 3]
        lna = np.log(a / p['a_i'])
        ax[0].plot(lna, (drho / p['rhozero']) * (a / p['a_i'])**0.25,
                   marker='o', ls='none', label=lab)
        ax[1].plot(lna, (dv / p['V_S']) * (a / p['a_i'])**0.75,
                   marker='o', ls='none', label=lab)
        ax[2].plot(lna, (dB / p['Bzero_i']) * (a / p['a_i'])**0.25,
                   marker='o', ls='none', label=lab)
    ax[0].plot(lac, drth * (ac / p['a_i'])**0.25, 'k-', label='linear theory')
    ax[1].plot(lac, duth * (ac / p['a_i'])**0.75, 'k-', label='linear theory')
    ax[2].plot(lac, dBth * (ac / p['a_i'])**0.25, 'k-', label='linear theory')
    ax[0].set_xlabel(r"$\ln(a/a_i)$")
    ax[0].set_ylabel(r"$(\delta\rho/\rho_0) \, (a/a_i)^{1/4}$")
    ax[1].set_xlabel(r"$\ln(a/a_i)$")
    ax[1].set_ylabel(r"$(\delta u / V_S) \, (a/a_i)^{3/4}$")
    ax[2].set_xlabel(r"$\ln(a/a_i)$")
    ax[2].set_ylabel(r"$(\delta B_c / B_{c,0}) \, (a/a_i)^{1/4}$")
    fig.suptitle(title)
    fig.tight_layout()
    finish_figure_with_legend(fig, (ax[0], ax[1], ax[2]),
                              save=(f"{save}_compressible.png" if save else None))


def plot_divclean(results, p, save=None, c_h=None):
    import matplotlib.pyplot as plt
    setup_rcparams()
    ac = _a_curve(results, p)
    lac = np.log(ac / p['a_i'])
    envelope, modulation = _divclean_hyperbolic(ac, p, c_h=c_h)
    analytic_dBx = p['A_u'] * p['Bzero_i'] * envelope          # max|dB_x|_c
    analytic_divBc = analytic_dBx * p['k']                     # max|div B_c|

    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    for lab, d in results:
        a = d[:, 0]
        lna = np.log(a / p['a_i'])
        ax[0].plot(lna, d[:, 4], marker='o', ls='none', label=f"{lab} sim")
        if not np.all(np.isnan(d[:, 7])):
            ax[1].plot(lna, d[:, 7], marker='o', ls='-',
                       label=r"$\|\nabla\!\cdot\!B_c\|_2$ " + lab)
        if not np.all(np.isnan(d[:, 8])):
            ax[1].plot(lna, d[:, 8], marker='s', ls='--',
                       label=r"$\max|\nabla\!\cdot\!B_c|$ " + lab)
        if not np.all(np.isnan(d[:, 9])):
            ax[1].plot(lna, d[:, 9], marker='^', ls=':',
                       label=r"$\max|\psi|$ " + lab)

    if modulation is None:
        ax[0].plot(lac, analytic_dBx, 'k-', label=r'envelope $A_u B_{c,0}(a/a_i)^2$')
        ax[1].plot(lac, analytic_divBc, 'k-', label=r'envelope $A_u B_{c,0} k (a/a_i)^2$')
    else:
        ax[0].plot(lac, np.abs(analytic_dBx * modulation), 'k-', label='analytic')
        ax[0].plot(lac, analytic_dBx, 'k--', alpha=0.4, label='envelope')
        ax[1].plot(lac, np.abs(analytic_divBc * modulation), 'k-', label='analytic')
        ax[1].plot(lac, analytic_divBc, 'k--', alpha=0.4, label='envelope')

    ax[0].set_yscale('log'); ax[1].set_yscale('log')
    ax[0].set_xlabel(r"$\ln(a/a_i)$"); ax[0].set_ylabel(r"$|\delta B_x|_c$")
    ax[1].set_xlabel(r"$\ln(a/a_i)$"); ax[1].set_ylabel("divB cleaning diagnostics")
    fig.suptitle("divB cleaning test (this work)"
                 + (f"  c_h={c_h}" if c_h is not None else "  (no c_h supplied)"))
    fig.tight_layout()
    finish_figure_with_legend(fig, (ax[0], ax[1]),
                              save=(f"{save}_divclean.png" if save else None))


# --------------------------------------------------------------------------- #
#  CLI
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
#  Regression baseline (bless / compare the amplitude-vs-a series)
# --------------------------------------------------------------------------- #

def reference_table(inp, p):
    """Regression table {x: a, dv, drho, [dB]} for run `inp`.

    The mode amplitudes (dv, drho, and dB when magnetic) vs the scale factor a
    are the headline Berlok diagnostics, so they make the CI baseline. dB is
    included only when finite for every snapshot (a magnetic run). None if no
    usable snapshots."""
    d = scan_input(inp, p)
    if d is None:
        return None
    table = {"x": d[:, 0], "dv": d[:, 1], "drho": d[:, 3]}
    if np.all(np.isfinite(d[:, 2])):
        table["dB"] = d[:, 2]
    return table


def compare_to_baseline(inp, p, explicit=None, ref_table=None):
    """Print residuals of the amplitude series vs the sibling baseline; return
    worst RMS (None if no baseline / no series). `ref_table` may be passed
    pre-discovered (the overlay already found the first input's baseline)."""
    if ref_table is None:
        ref_table, _path = find_reference(inp, explicit)
    if ref_table is None:
        return None
    table = reference_table(inp, p)
    if table is None:
        return None
    worst = 0.0
    for key in ("dv", "drho", "dB"):
        res = residuals(table, ref_table, key)
        if res is None:
            continue
        mx, rms = res
        worst = max(worst, rms)
        print(f"cosmowave[regression] {key:<4s} vs baseline: "
              f"max|d|={mx:.4g}  rms={rms:.4g}")
    return worst


def main():
    parser = build_parser(PARAM_SPEC,
                          description="Cosmological (M)HD wave analysis "
                                      "(Berlok 2022).")
    parser.add_argument("--ch", type=float, default=None,
                        help="run-time cleaning speed c_h (case 4 only); enables "
                             "the cos[2 Omega_h(sqrt a - sqrt a_i)] modulation, "
                             "Omega_h = k*c_h/H_0")
    parser.add_argument("--reg-tol", type=float, default=None,
                        help="CI gate: exit non-zero if the worst RMS residual "
                             "vs the regression baseline exceeds this tolerance")
    args = parser.parse_args()
    params = params_from_args(args, PARAM_SPEC)
    case = int(params["case"])

    labels = args.labels if args.labels else [
        os.path.basename(os.path.normpath(i)) for i in args.inputs]
    if len(labels) != len(args.inputs):
        parser.error("--labels must match the number of inputs")

    p = case_params(case, factor_A=params["factor_A"], factor_G=params["factor_G"])
    print(f"[cosmowave case {case}] kappa_A={p['kappa_A']:.4g} kappa={p['kappa']:.4g} "
          f"a_i={p['a_i']:.4g} V_A={p['V_A']:.4g} V_S={p['V_S']:.4g} "
          f"Bzero_i={p['Bzero_i']:.4g} rhozero={p['rhozero']:.4g}")

    # --save-reference: bless the FIRST input's amplitude series as the baseline.
    if args.save_reference:
        table = reference_table(args.inputs[0], p)
        if table is None:
            print("cosmowave: no amplitude series to save as reference",
                  file=sys.stderr)
            sys.exit(1)
        save_reference(table, reference_path_for(args.inputs[0]), params,
                       "cosmowave")
        return

    results = []
    for inp, lab in zip(args.inputs, labels):
        d = scan_input(inp, p)
        if d is not None:
            results.append((lab, d))
    if not results:
        print("cosmowave: no data to plot; exiting.", file=sys.stderr)
        sys.exit(1)

    # Reference comparison is opt-in (only when --reference is given). Discover
    # the FIRST input's baseline once -- overlaid on the plots and reused for the
    # first input's regression gate.
    ref_table = None
    if args.reference is not None:
        ref_table, _rp = find_reference(args.inputs[0], args.reference)

    if case in (0, 1):
        plot_alfven(results, p, save=args.save, ref_table=ref_table)
    elif case in (2, 3):
        plot_compressible(results, p, save=args.save, ref_table=ref_table)
    else:
        plot_divclean(results, p, save=args.save, c_h=args.ch)

    worst = None
    for i, inp in enumerate(args.inputs):
        w = compare_to_baseline(inp, p, args.reference,
                                ref_table=(ref_table if i == 0 else None))
        if w is not None:
            worst = w if worst is None else max(worst, w)
    if args.reg_tol is not None and worst is not None and worst > args.reg_tol:
        print(f"cosmowave: REGRESSION FAILED -- worst rms {worst:.4g} > tol "
              f"{args.reg_tol:.4g}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
