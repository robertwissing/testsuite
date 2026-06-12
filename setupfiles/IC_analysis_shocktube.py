#!/usr/bin/env python3
"""
Analysis for the shock-tube / Riemann-problem suite (IC_setup_shocktube).

A 1-D Riemann problem: a left state and a right state meet at x=0 and evolve
self-similarly. The runs use PERIODIC boundaries, so the IC actually has TWO
discontinuities -- the L|R jump at x=0 AND its periodic image (an R|L jump) at the
wrapped box edge x=+-L/2 -- and a periodic run develops a second Riemann fan from
the boundary. The hydro exact overlay (`exact_riemann_periodic`) draws BOTH fans
(valid until they collide). `choice` (-a) selects one of the presets ported
verbatim from IC_setup_shocktube.choose_shock
(state = [dens, press, vx, vy, vz, Bx, By, Bz]):

    1  Sod shock              (hydro)        7  C-shock        (isothermal MHD)
    2  Ryu 1a                 (MHD)          8  Steady shock   (isothermal MHD)
    3  Ryu 1b                 (MHD)          9  Toth (=Ryu 1a) (MHD)
    4  Ryu 2a                 (MHD)         10  MHD discontinuity
    5  Ryu 2b                 (MHD)         11  strong blast    (hydro, g=1.4)
    6  Brio-Wu (Ryu 5a)       (MHD)

The analysis plots the profiles of rho, v and P (plus the transverse v and B
components for the MHD presets) vs x at the final snapshot, scattered over all
particles. For the PURE-HYDRO presets (B=0, gamma>1: Sod and strong-blast) it
overlays the EXACT adiabatic Riemann solution (Toro 2009, ch. 4) and prints the
volume-weighted L1 error of rho, v and P. For the standard MHD presets
(Brio-Wu #6, Ryu 1a #2/Toth #9, Ryu 2a #4) it overlays the tabulated
piecewise-constant reference solution extracted from SPLASH's exact_mhdshock.f90
(Ryu & Jones 1995 / Brio & Wu 1988; the MHD Riemann problem has no simple closed
form) and L1s against it. The remaining MHD / isothermal presets get profiles
only -- compare by eye, or pass a high-resolution run as a second input.

Standalone like sedov/noh: reuses the framework Param/build_parser/loaders, not
run_cli. Several inputs overlay as separate marker series.

Usage:
    python IC_analysis_shocktube.py <run-dir> [-a 1] [--time T] [--save out]
    python IC_analysis_shocktube.py <lowres> <hires> -a 6 --labels low high
"""

import os
import sys

import numpy as np
import readtipsy as tip

from IC_analysis_framework import (
    Param, build_parser, params_from_args,
    reference_path_for, save_reference, find_reference, residuals, reg_nbins,
    select_at_time, binned_profile,
)
from IC_analysis_general import (
    find_files, loaddata, read_bfield, read_tufac,
    setup_rcparams, finish_figure_with_legend,
)

CONST = np.sqrt(4.0 * np.pi)     # IC_setup_shocktube.choose_shock normalisation

# Per-choice presets, ported verbatim from IC_setup_shocktube.choose_shock.
# state = [dens, press, vx, vy, vz, Bx, By, Bz]; xshock = 0 for every preset
# (xshock = (xleft+xright)/2, xright=-xleft). `tmax` is informational only --
# the exact solution is sampled at the actual snapshot time.
SHOCKS = {
    1:  dict(name="Sod shock", gamma=5./3., tmax=0.2,
             L=[1.0, 1.0, 0., 0., 0., 0., 0., 0.],
             R=[0.125, 0.1, 0., 0., 0., 0., 0., 0.]),
    2:  dict(name="Ryu 1a", gamma=5./3., tmax=0.08,
             L=[1., 20., 10., 0., 0., 5./CONST, 5./CONST, 0.],
             R=[1., 1., -10., 0., 0., 5./CONST, 5./CONST, 0.]),
    3:  dict(name="Ryu 1b", gamma=5./3., tmax=0.06,
             L=[1.0, 1., 0., 0., 0., 5./CONST, 5./CONST, 0.],
             R=[0.1, 10., 0., 0., 0., 5./CONST, 2./CONST, 0.]),
    4:  dict(name="Ryu 2a", gamma=5./3., tmax=0.2,
             L=[1.08, 0.95, 1.2, 0.01, 0.5, 2./CONST, 3.6/CONST, 2./CONST],
             R=[1., 1., 0., 0., 0., 2./CONST, 4.0/CONST, 2./CONST]),
    5:  dict(name="Ryu 2b", gamma=5./3., tmax=0.035,
             L=[1.0, 1., 0., 0., 0., 3./CONST, 6./CONST, 0.],
             R=[0.1, 10., 0., 2., 1., 3./CONST, 1./CONST, 0.]),
    6:  dict(name="Brio-Wu", gamma=2.0, tmax=0.1,
             L=[1.0, 1.0, 0., 0., 0., 0.75, 1., 0.],
             R=[0.125, 0.1, 0., 0., 0., 0.75, -1., 0.]),
    7:  dict(name="C-shock", gamma=1.0, cs2=0.01, tmax=0.8,
             L=[1., 0.006, 4.45, 0., 0., 1./np.sqrt(2.), 1./np.sqrt(2.), 0.],
             R=[1., 0.006, -4.45, 0., 0., 1./np.sqrt(2.), 1./np.sqrt(2.), 0.]),
    8:  dict(name="Steady shock", gamma=1.0, cs2=0.01, tmax=0.8,
             L=[1.7942, 0.017942, -0.9759, -0.6561, 0., 1., 1.74885, 0.],
             R=[1., 0.01, -1.7510, 0., 0., 1., 0.6, 0.]),
    9:  dict(name="Toth", gamma=5./3., tmax=0.08,
             L=[1., 20., 10., 0., 0., 5./CONST, 5./CONST, 0.],
             R=[1., 1., -10., 0., 0., 5./CONST, 5./CONST, 0.]),
    10: dict(name="MHD discontinuity", gamma=5./3., tmax=0.2,
             L=[1., 1., 0., 0., 0., 20./CONST, 20./CONST, 0.],
             R=[1., 1., 0., 0., 0., 5./CONST, 5./CONST, 0.]),
    11: dict(name="strong blast", gamma=1.4, tmax=0.012,
             L=[1., 1000., 0., 0., 0., 0., 0., 0.],
             R=[1., 0.1, 0., 0., 0., 0., 0., 0.]),
    12: dict(name="fast slow shock", gamma=5./3., tmax=0.15,
             L=[1.0, 1.0, 0., 0., 0., 1.0, 1.0, 0.],
             R=[0.2, 0.1, 0., 0., 0., 1.0, 0.0, 0.]),
    13: dict(name="isothermal B98", gamma=1.0, cs2=1.0, tmax=0.2,
             L=[1.08, 1.08, 1.2, 0.01, 0.5, 2./CONST, 3.6/CONST, 2./CONST],
             R=[1.0, 1.0, 0., 0., 0., 2./CONST, 4.0/CONST, 2./CONST]),
    14: dict(name="Mach 25 DW94", gamma=5./3., tmax=0.03,
             L=[1.0, 1.0, 36.87, -0.1546, -0.03864, 4./CONST, 4./CONST, 1./CONST],
             R=[0.1, 1.0, -36.87, 0., 0., 4./CONST, 4./CONST, 1./CONST]),
    15: dict(name="rarefaction", gamma=5./3., tmax=0.1,
             L=[1.0, 1.0, -1.0, 0., 0., 0., 1.0, 0.],
             R=[1.0, 1.0, 1.0, 0., 0., 0., 1.0, 0.]),
}

XSHOCK = 0.0
DEFAULT_TIME = None             # None -> the last (final) snapshot

PARAM_SPEC = [
    Param("a", "choice", 1, int, "shock preset (1..11; see module docstring)"),
]


def shock_params(choice):
    c = int(choice)
    if c not in SHOCKS:
        print(f"shocktube: unknown choice {c}; valid: {sorted(SHOCKS)}",
              file=sys.stderr)
        sys.exit(1)
    s = SHOCKS[c]
    L, R = np.asarray(s["L"], float), np.asarray(s["R"], float)
    mhd = bool(np.any(L[5:] != 0) or np.any(R[5:] != 0))
    hydro_exact = (not mhd) and s["gamma"] > 1.0 + 1e-6   # exact solver applies
    ref = CHOICE_TO_REF.get(c)                            # tabulated MHD ref, or None
    return dict(s, choice=c, L=L, R=R, mhd=mhd, hydro_exact=hydro_exact, ref=ref)


# --------------------------------------------------------------------------- #
#  Exact adiabatic Riemann solver (Toro 2009, ch. 4) -- hydro presets only
# --------------------------------------------------------------------------- #

def _fK(p, rhoK, pK, aK, g):
    """Pressure function and derivative for one side (shock if p>pK else fan)."""
    if p > pK:                                   # shock
        AK = 2.0 / ((g + 1.0) * rhoK)
        BK = (g - 1.0) / (g + 1.0) * pK
        f = (p - pK) * np.sqrt(AK / (p + BK))
        df = np.sqrt(AK / (BK + p)) * (1.0 - 0.5 * (p - pK) / (BK + p))
    else:                                        # rarefaction
        f = 2.0 * aK / (g - 1.0) * ((p / pK) ** ((g - 1.0) / (2.0 * g)) - 1.0)
        df = (1.0 / (rhoK * aK)) * (p / pK) ** (-(g + 1.0) / (2.0 * g))
    return f, df


def _star_state(WL, WR, g):
    """Star-region (p*, u*) by Newton iteration on the pressure function."""
    rhoL, uL, pL = WL
    rhoR, uR, pR = WR
    aL, aR = np.sqrt(g * pL / rhoL), np.sqrt(g * pR / rhoR)
    p = max(1e-8, 0.5 * (pL + pR)
            - 0.125 * (uR - uL) * (rhoL + rhoR) * (aL + aR))   # PVRS guess
    for _ in range(100):
        fL, dfL = _fK(p, rhoL, pL, aL, g)
        fR, dfR = _fK(p, rhoR, pR, aR, g)
        pnew = p - (fL + fR + (uR - uL)) / (dfL + dfR)
        if pnew < 0:
            pnew = 1e-8
        if abs(pnew - p) < 1e-12 * 0.5 * (pnew + p):
            p = pnew
            break
        p = pnew
    fL, _ = _fK(p, rhoL, pL, aL, g)
    fR, _ = _fK(p, rhoR, pR, aR, g)
    ustar = 0.5 * (uL + uR) + 0.5 * (fR - fL)
    return p, ustar


def _sample(WL, WR, g, pstar, ustar, xi):
    """(rho, u, p) of the exact solution at self-similar speed xi=(x-xshock)/t."""
    rhoL, uL, pL = WL
    rhoR, uR, pR = WR
    aL, aR = np.sqrt(g * pL / rhoL), np.sqrt(g * pR / rhoR)
    gm, gp = g - 1.0, g + 1.0
    if xi <= ustar:                              # left of the contact
        if pstar > pL:                           # left shock
            SL = uL - aL * np.sqrt(gp / (2 * g) * pstar / pL + gm / (2 * g))
            if xi <= SL:
                return rhoL, uL, pL
            rho = rhoL * (pstar / pL + gm / gp) / (gm / gp * pstar / pL + 1.0)
            return rho, ustar, pstar
        aLs = aL * (pstar / pL) ** (gm / (2 * g))   # left rarefaction
        SHL, STL = uL - aL, ustar - aLs
        if xi <= SHL:
            return rhoL, uL, pL
        if xi >= STL:
            return rhoL * (pstar / pL) ** (1.0 / g), ustar, pstar
        u = 2.0 / gp * (aL + 0.5 * gm * uL + xi)        # inside fan
        a = 2.0 / gp * (aL + 0.5 * gm * (uL - xi))
        return rhoL * (a / aL) ** (2.0 / gm), u, pL * (a / aL) ** (2.0 * g / gm)
    else:                                        # right of the contact
        if pstar > pR:                           # right shock
            SR = uR + aR * np.sqrt(gp / (2 * g) * pstar / pR + gm / (2 * g))
            if xi >= SR:
                return rhoR, uR, pR
            rho = rhoR * (pstar / pR + gm / gp) / (gm / gp * pstar / pR + 1.0)
            return rho, ustar, pstar
        aRs = aR * (pstar / pR) ** (gm / (2 * g))   # right rarefaction
        SHR, STR = uR + aR, ustar + aRs
        if xi >= SHR:
            return rhoR, uR, pR
        if xi <= STR:
            return rhoR * (pstar / pR) ** (1.0 / g), ustar, pstar
        u = 2.0 / gp * (-aR + 0.5 * gm * uR + xi)
        a = 2.0 / gp * (aR - 0.5 * gm * (uR - xi))
        return rhoR * (a / aR) ** (2.0 / gm), u, pR * (a / aR) ** (2.0 * g / gm)


def exact_riemann(p, t, x):
    """Exact (rho, vx, P) of the adiabatic Riemann problem at time t, positions x.

    Uses the preset's left/right primitive state (density, pressure, vx) and
    gamma. Only meaningful for the hydro presets (B=0, gamma>1). This is the
    single-discontinuity (non-periodic) solution centred at XSHOCK=0."""
    g = p["gamma"]
    L, R = p["L"], p["R"]
    WL = (L[0], L[2], L[1])                      # (rho, vx, P)
    WR = (R[0], R[2], R[1])
    pstar, ustar = _star_state(WL, WR, g)
    xi = (np.asarray(x, float) - XSHOCK) / t
    rho = np.empty_like(xi)
    vx = np.empty_like(xi)
    P = np.empty_like(xi)
    for i, xii in enumerate(xi):
        rho[i], vx[i], P[i] = _sample(WL, WR, g, pstar, ustar, xii)
    return rho, vx, P


def exact_riemann_periodic(p, t, x, xmin, xmax):
    """Exact (rho, vx, P) on a PERIODIC domain [xmin, xmax] (xmin = -xmax).

    The IC has the L|R jump at XSHOCK=0 AND its periodic image -- an R|L jump --
    at the wrapped boundary x = xmin == xmax. So a periodic run develops a SECOND
    Riemann fan emanating from the box edge (states swapped). Before the two fans
    collide each point is governed by its nearer discontinuity, so we sample the
    central fan (L|R, centred at 0) near the middle and the boundary fan (R|L,
    centred at the edge) near the edges -- the two share the same constant
    plateaus, so the switch at the midpoint |x| = xmax/2 is seamless. Valid until
    the fans interact (~ t < (xmax/2)/max wave speed); past that the run itself is
    no longer a clean pair of Riemann problems either.
    """
    g = p["gamma"]
    L, R = p["L"], p["R"]
    WL = (L[0], L[2], L[1])                      # (rho, vx, P)  central left
    WR = (R[0], R[2], R[1])                      #               central right
    ps_c, us_c = _star_state(WL, WR, g)          # central fan (L|R) star state
    ps_b, us_b = _star_state(WR, WL, g)          # boundary fan (R|L), states swapped
    x = np.asarray(x, float)
    rho = np.empty_like(x)
    vx = np.empty_like(x)
    P = np.empty_like(x)
    for i, xx in enumerate(x):
        xb = xmax if xx >= XSHOCK else xmin      # nearest periodic-image boundary
        if abs(xx - XSHOCK) <= abs(xx - xb):     # nearer the centre -> central fan
            rho[i], vx[i], P[i] = _sample(WL, WR, g, ps_c, us_c,
                                          (xx - XSHOCK) / t)
        else:                                    # nearer the edge -> boundary fan
            rho[i], vx[i], P[i] = _sample(WR, WL, g, ps_b, us_b, (xx - xb) / t)
    return rho, vx, P


# --------------------------------------------------------------------------- #
#  Tabulated MHD reference solutions (from SPLASH exact_mhdshock.f90)
# --------------------------------------------------------------------------- #
# The MHD Riemann problem has no simple closed form, so -- exactly as SPLASH
# does -- we ship piecewise-constant reference states from the literature
# (Ryu & Jones 1995; Brio & Wu 1988). Each solution is a polyline: wave
# positions scale linearly in time (x = coeff * t/tref), constant state between
# them; the two ends run to the plot domain. B is in the same 1/sqrt(4pi)-
# normalised units as the gasoline IC (e.g. Ryu By,L = 5/sqrt(4pi) = 1.4105).
# `xc` holds the interior wave-position coefficients (at t=tref); the endpoints
# (None) extend to xmin/xmax. Transcribed verbatim from exact_mhdshock.f90.

_cc = 1.0 / CONST   # = 1/sqrt(4*pi)

MHD_REF = {
    "briowu": dict(   # SPLASH ishk=1 (Brio & Wu, gamma=2), tref=0.1
        tref=0.1, Bx0=0.75,
        xc=[None, -0.18, -0.08, -0.03, -0.03, -0.03, -0.005, 0.06, 0.06,
            0.147, 0.147, 0.33, 0.36, None],
        rho=[1, 1, 0.67623, 0.67623, 0.827, 0.775, 0.6962, 0.6962, 0.2352,
             0.2352, 0.117, 0.117, 0.125, 0.125],
        P=[1, 1, 0.447, 0.447, 0.727219, 0.6, 0.516, 0.516, 0.516, 0.516,
           0.0876, 0.0876, 0.1, 0.1],
        vx=[0, 0, 0.63721, 0.63721, 0.48, 0.52, 0.6, 0.6, 0.6, 0.6,
            -0.24, -0.24, 0, 0],
        vy=[0, 0, -0.23345, -0.23345, -1.3, -1.4, -1.584, -1.584, -1.584,
            -1.584, -0.166, -0.166, 0, 0],
        vz=[0] * 14,
        By=[1, 1, 2.1 * _cc, 2.1 * _cc, -1.2 * _cc, -1.3 * _cc, -1.9 * _cc,
            -1.9 * _cc, -1.9 * _cc, -1.9 * _cc, -3.25 * _cc, -3.25 * _cc, -1, -1],
        Bz=[0] * 14),
    "ryu2a": dict(   # SPLASH ishk=3 (RJ95 7-discontinuity), tref=0.2
        tref=0.2, Bx0=2 * _cc,
        xc=[None, -0.19, -0.19, 0.03, 0.03, 0.051, 0.051, 0.12, 0.12,
            0.18, 0.18, 0.205, 0.205, 0.45, 0.45, None],
        rho=[1.08, 1.08, 1.4903, 1.4903, 1.6343, 1.6343, 1.6343, 1.6343,
             1.4735, 1.4735, 1.3090, 1.3090, 1.3090, 1.3090, 1.0, 1.0],
        P=[0.95, 0.95, 1.6558, 1.6558, 1.9317, 1.9317, 1.9317, 1.9317,
           1.9317, 1.9317, 1.5844, 1.5844, 1.5844, 1.5844, 1.0, 1.0],
        vx=[1.2, 1.2, 0.60588, 0.60588, 0.60588, 0.57538, 0.57538, 0.57538,
            0.57538, 0.57538, 0.57538, 0.53432, 0.53432, 0.53432, 0, 0],
        vy=[0.01, 0.01, 0.11235, 0.11235, 0.22157, 0.22157, 0.047602, 0.047602,
            0.047601, 0.047601, -0.18411, -0.18411, -0.094572, -0.094572, 0, 0],
        vz=[0.5, 0.5, 0.55686, 0.55686, 0.30125, 0.30125, 0.24734, 0.24734,
            0.24734, 0.24734, 0.17554, 0.17554, -0.047286, -0.047286, 0, 0],
        By=[1.0155, 1.0155, 1.4383, 1.4383, 1.5716, 1.5716, 1.4126, 1.4126,
            1.4126, 1.4126, 1.6103, 1.6103, 1.5078, 1.5078, 1.1284, 1.1284],
        Bz=[0.56419, 0.56419, 0.79907, 0.79907, 0.48702, 0.48702, 0.43772,
            0.43772, 0.43772, 0.43772, 0.49899, 0.49899, 0.75392, 0.75392,
            0.56419, 0.56419]),
    "ryu1a": dict(   # SPLASH ishk=7 (RJ95 problem 1A), tref=0.08
        tref=0.08, Bx0=5 * _cc,
        xc=[None, -0.386, -0.386, -0.01, -0.01, 0.0505, 0.0505, 0.12, 0.12,
            0.37, 0.37, None],
        rho=[1, 1, 2.6797, 2.6797, 2.6713, 2.6713, 3.8508, 3.8508, 3.7481,
             3.7481, 1, 1],
        P=[20, 20, 150.98, 150.98, 150.19, 150.19, 150.19, 150.19, 143.57,
           143.57, 1, 1],
        vx=[10, 10, 0.72113, 0.72113, 0.72376, 0.72376, 0.72376, 0.72376,
            0.70505, 0.70505, -10, -10],
        vy=[0, 0, 0.23139, 0.23139, 0.35684, 0.35684, 0.35684, 0.35684,
            -0.38804, -0.38804, 0, 0],
        vz=[0] * 12,
        By=[1.4105, 1.4105, 3.8389, 3.8389, 4.0380, 4.0380, 4.0380, 4.0380,
            5.4272, 5.4272, 1.4105, 1.4105],
        Bz=[0] * 12),
    "fastslow": dict(   # SPLASH ishk=2 (RJ95 fast/slow); B in RAW units, tref=0.15
        tref=0.15, Bx0=1.0,
        xc=[None, -0.27, -0.09, -0.03, -0.01, 0.135, 0.135, 0.25, 0.25,
            0.35, 0.35, None],
        rho=[1, 1, 0.5955, 0.5955, 0.55151, 0.55151, 0.41272, 0.41272,
             0.2337, 0.2337, 0.2, 0.2],
        P=[1, 1, 0.42629, 0.42629, 0.37090, 0.37090, 0.37090, 0.37090,
           0.12402, 0.12402, 0.1, 0.1],
        vx=[0, 0, 0.81237, 0.81237, 0.89416, 0.89416, 0.89416, 0.89416,
            0.24722, 0.24722, 0, 0],
        vy=[0, 0, -0.59961, -0.59961, -0.5447, -0.5447, -0.5447, -0.5447,
            -0.91164, -0.91164, 0, 0],
        vz=[0] * 12,
        By=[1, 1, 0.28431, 0.28431, 0.31528, 0.31528, 0.31528, 0.31528,
            0.43086, 0.43086, 0, 0],
        Bz=[0] * 12),
    "isob98": dict(   # SPLASH ishk=4 (Balsara 1998 isothermal); P=rho, tref=0.2
        tref=0.2, Bx0=2 * _cc,
        xc=[None, -0.15, -0.15, 0.035, 0.035, 0.07, 0.07, 0.17, 0.17,
            0.2, 0.2, 0.41, 0.41, None],
        rho=[1.08, 1.08, 1.515, 1.515, 1.515, 1.515, 1.745, 1.745,
             1.36, 1.36, 1.36, 1.36, 1.0, 1.0],
        vx=[1.2, 1.2, 0.65, 0.65, 0.65, 0.65, 0.62, 0.62,
            0.54, 0.54, 0.54, 0.54, 0, 0],
        vy=[0.01, 0.01, 0.13, 0.13, 0.24, 0.24, 0.071, 0.071,
            -0.215, -0.215, -0.125, -0.125, 0, 0],
        vz=[0.5, 0.5, 0.57, 0.57, 0.31, 0.31, 0.255, 0.255,
            0.165, 0.165, -0.06, -0.06, 0, 0],
        By=[3.6 * _cc, 3.6 * _cc, 5.2 * _cc, 5.2 * _cc, 5.7 * _cc, 5.7 * _cc,
            5.22 * _cc, 5.22 * _cc, 5.96 * _cc, 5.96 * _cc, 5.58 * _cc,
            5.58 * _cc, 4.0 * _cc, 4.0 * _cc],
        Bz=[2 * _cc, 2 * _cc, 2.885 * _cc, 2.885 * _cc, 1.76 * _cc, 1.76 * _cc,
            1.62 * _cc, 1.62 * _cc, 1.85 * _cc, 1.85 * _cc, 2.79 * _cc,
            2.79 * _cc, 2 * _cc, 2 * _cc],
        # P = rho (isothermal, cs^2=1) -- filled in below
        P=None),
    "mach25": dict(   # SPLASH ishk=6 (Dai & Woodward 1994 Mach 25), tref=0.03
        tref=0.03, Bx0=4 * _cc,
        xc=[None, -0.35, -0.35, 0.35, 0.35, None],
        rho=[1.0, 1.0, 3.982, 3.982, 0.1, 0.1],
        P=[1.0, 1.0, 1806.0, 1806.0, 1.0, 1.0],
        vx=[36.87, 36.87, 0.0, 0.0, -36.87, -36.87],
        vy=[-0.1546, -0.1546, -0.07727, -0.07727, 0.0, 0.0],
        vz=[-0.03864, -0.03864, -0.01932, -0.01932, 0.0, 0.0],
        By=[4 * _cc, 4 * _cc, 15.95 * _cc, 15.95 * _cc, 4 * _cc, 4 * _cc],
        Bz=[1 * _cc, 1 * _cc, 3.988 * _cc, 3.988 * _cc, 1 * _cc, 1 * _cc]),
    "rarefaction": dict(   # SPLASH ishk=5 (RJ95 two-rarefaction); Bx=0, tref=0.1
        tref=0.1, Bx0=0.0,
        xc=[None, -0.27, -0.12, 0.12, 0.27, None],
        rho=[1, 1, 0.49653, 0.49653, 1, 1],
        P=[1, 1, 0.31134, 0.31134, 1, 1],
        vx=[-1, -1, 0, 0, 1, 1],          # middle vx~0 (approximate, RJ95)
        vy=[0] * 6,
        vz=[0] * 6,
        By=[1, 1, 0.49638, 0.49638, 1, 1],
        Bz=[0] * 6),
}
MHD_REF["isob98"]["P"] = list(MHD_REF["isob98"]["rho"])   # isothermal: P = rho

# my preset choice -> SPLASH tabulated reference (same initial state)
CHOICE_TO_REF = {2: "ryu1a", 9: "ryu1a", 4: "ryu2a", 6: "briowu",
                 12: "fastslow", 13: "isob98", 14: "mach25", 15: "rarefaction"}


# --------------------------------------------------------------------------- #
#  Exact oblique C-shock (ported from SPLASH exact_Cshock.f90)
# --------------------------------------------------------------------------- #
# Isothermal MHD shock with ambipolar diffusion, Mac-Low, Norman, Konigl &
# Wardle (1995), ApJ 442, 726. The neutral-density ratio D(x) satisfies a 1st-
# order ODE integrated from the upstream side; the field and velocities follow
# from conservation. Parametrised by the sonic & Alfven Mach numbers (machs,
# macha) and the field angle theta; `shockl` is the ambipolar shock thickness
# (it only rescales x -- the profile SHAPE is set by machs/macha/theta). Fixed
# constants match SPLASH: cs=0.1, rho_n0=1, B0=1 (-> Bx0=By0=1/sqrt2 for
# theta=pi/4), upstream vx0=-5, inflow vxin=-4.45.
CSHOCK_CS = 0.1
CSHOCK_RHON0 = 1.0
CSHOCK_B0 = 1.0
CSHOCK_VX0 = -5.0
CSHOCK_VXIN = -4.45


def _cshock_getb(b0, macha, machs, D):
    arg = b0 ** 2 + 2. * macha ** 2 * (D - 1.) * (1. / D - 1. / machs ** 2)
    return np.sqrt(np.maximum(arg, 0.0))           # clamp at the sonic point


def _cshock_D(machs, macha, theta, shockl, xshock, xgrid, n=4000):
    """Integrate the neutral-density ratio D(x) (=1 upstream, x>=xshock).

    Midpoint-rule integration from the shock downstream, mirroring SPLASH's
    `integrate`. D rises monotonically from 1 to the post-shock value D2 (where
    the source term vanishes) over a width ~shockl; once equilibrium is reached
    (or the sub-magnetosonic D=machs point is approached) D is held constant, so
    the profile stays finite for any Mach numbers."""
    b0, sinth, cos2 = np.sin(theta), np.sin(theta), np.cos(theta) ** 2
    Dmax = 0.999 * machs ** 2                       # stay below the sonic point

    def f(D):
        b = _cshock_getb(b0, macha, machs, D)
        num = b / macha * (b - D * ((b - b0) / macha ** 2 * cos2 + sinth))
        den = (b ** 2 + cos2) * (1. / D ** 2 - 1. / machs ** 2)
        return num / den if den > 1e-30 else 0.0

    xmin = min(xgrid.min(), xshock)
    xs = np.linspace(xshock, xmin, n)              # descending from the shock
    dx = abs(xs[1] - xs[0])
    D = np.empty(n)
    D[0] = 1.0 + 1e-6
    for i in range(1, n):
        if D[i - 1] >= Dmax or f(D[i - 1]) <= 0.0:  # reached post-shock D2
            D[i] = D[i - 1]
            continue
        Dh = min(D[i - 1] + 0.5 * dx / shockl * f(D[i - 1]), Dmax)
        D[i] = min(max(D[i - 1] + dx / shockl * f(Dh), 1.0), Dmax)
    Dg = np.interp(xgrid, xs[::-1], D[::-1])
    return np.where(xgrid >= xshock, 1.0, Dg)


def cshock_reference(machs, macha, theta, shockl, t, x):
    """Polyline-free C-shock reference at time t on positions x: dict of
    (rho, vx, vy, By, Bx). The shock sits at xshock = (6/8) v_A t."""
    x = np.asarray(x, float)
    va = CSHOCK_B0 / np.sqrt(CSHOCK_RHON0)
    xshock = 0.75 * va * t
    b0 = np.sin(theta)
    Bx = CSHOCK_B0 * np.cos(theta)
    D = _cshock_D(machs, macha, theta, shockl, xshock, x)
    # upstream invariants (far upstream, D=1)
    By0 = CSHOCK_B0 * _cshock_getb(b0, macha, machs, 1.0)
    rhon_pre = CSHOCK_RHON0
    Pr0 = rhon_pre * CSHOCK_CS ** 2
    K1 = Pr0 + 0.5 * By0 ** 2 + rhon_pre * CSHOCK_VX0 ** 2
    K2 = rhon_pre * CSHOCK_VX0 * 0.0 - Bx * By0
    dvx = CSHOCK_VXIN - CSHOCK_VX0
    rhon = D * CSHOCK_RHON0
    By = CSHOCK_B0 * _cshock_getb(b0, macha, machs, D)
    vx2 = (K1 - 0.5 * By ** 2 - rhon * CSHOCK_CS ** 2) / rhon
    with np.errstate(invalid="ignore", divide="ignore"):
        vx = np.where(vx2 > 0, -np.sqrt(np.abs(vx2)), 0.0)
        vy = np.where(vx2 > 0, (K2 + Bx * By) / (rhon * vx), 0.0)
    vx = vx + dvx
    return dict(rho=rhon, vx=vx, vy=vy, By=By, Bx=np.full_like(x, Bx))


def steady_step(p, x):
    """Stationary-step exact solution for the steady shock (#8, Falle 2003).

    The left/right states satisfy the steady MHD jump conditions (mass,
    momentum and induction fluxes are equal across -- verified), so the exact
    solution is the initial discontinuity held FIXED at xshock=0 for all time.
    P is the isothermal cs^2*rho. Returns a dict of state arrays vs x."""
    x = np.asarray(x, float)
    left = x < XSHOCK
    L, R = p["L"], p["R"]
    cs2 = p.get("cs2", 1.0)
    rho = np.where(left, L[0], R[0])
    return dict(rho=rho,
                vx=np.where(left, L[2], R[2]),
                vy=np.where(left, L[3], R[3]),
                vz=np.where(left, L[4], R[4]),
                By=np.where(left, L[6], R[6]),
                Bz=np.where(left, L[7], R[7]),
                P=cs2 * rho)


def mhd_reference(refname, t, xmin, xmax):
    """Polyline (x, {rho,vx,vy,vz,P,By,Bz}) of the tabulated MHD reference at t.

    Wave positions scale as x = coeff * t/tref; the two endpoints extend to
    [xmin, xmax] at the left/right states. Returns x and a dict of state arrays
    (same length), suitable both for plotting and for np.interp-based L1."""
    r = MHD_REF[refname]
    tfac = t / r["tref"]
    x = np.array([(xmin if c is None and i == 0 else
                   xmax if c is None else c * tfac)
                  for i, c in enumerate(r["xc"])], float)
    state = {k: np.asarray(r[k], float)
             for k in ("rho", "vx", "vy", "vz", "P", "By", "Bz")}
    return x, state


# --------------------------------------------------------------------------- #
#  Profile extraction
# --------------------------------------------------------------------------- #

def profile(fn, params):
    """Per-particle profile dict (x, rho, vx, vy, vz, P, u, By, Bz) + time.

    P = (gamma-1) rho u with u = dTuFac * T (col 8); B from the BField aux."""
    g = params["gamma"]
    tufac = read_tufac(os.path.dirname(fn) or ".")
    tgdata, _td, _ts, _hdr, time, _N, ngas, _nd, _ns, _h = loaddata(fn)
    u = tufac * tgdata[:, 8]
    rho = tgdata[:, 7]
    # P = (gamma-1) rho u, except isothermal (gamma=1) where P = cs^2 rho.
    P = (params.get("cs2", 1.0) * rho if g <= 1.0 + 1e-6
         else (g - 1.0) * rho * u)
    d = dict(x=tgdata[:, 1], rho=rho,
             vx=tgdata[:, 4], vy=tgdata[:, 5], vz=tgdata[:, 6],
             u=u, P=P, vol=tgdata[:, 0] / rho)
    B = read_bfield(fn, ngas)
    if B is not None:
        B = np.asarray(B, float)
        d["By"], d["Bz"] = B[:, 1], B[:, 2]
    else:
        d["By"] = d["Bz"] = np.zeros(len(d["x"]))
    return d, float(np.asarray(time).ravel()[0])


def l1_error(x, field, exact, vol):
    """Volume-weighted L1 error mean_V |field - exact|."""
    return float(np.sum(vol * np.abs(field - exact)) / np.sum(vol))


# --------------------------------------------------------------------------- #
#  Plot
# --------------------------------------------------------------------------- #

def plot_profiles(inputs, labels, params, t, save=None, cshock_opts=None,
                  ref_table=None):
    import matplotlib.pyplot as plt
    setup_rcparams()
    p = params
    # Panel set: hydro -> rho,vx,P,u ; MHD -> rho,vx,vy,P,By,Bz.
    if p["mhd"]:
        keys = [("rho", r"$\rho$"), ("vx", r"$v_x$"), ("vy", r"$v_y$"),
                ("P", r"$P$"), ("By", r"$B_y$"), ("Bz", r"$B_z$")]
        nrow, ncol = 2, 3
    else:
        keys = [("rho", r"$\rho$"), ("vx", r"$v_x$"),
                ("P", r"$P$"), ("u", r"$u$")]
        nrow, ncol = 2, 2
    fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 4 * nrow),
                             sharex=True)
    axd = {k: axes.ravel()[i] for i, (k, _lab) in enumerate(keys)}

    xs_all = []
    for i, (inp, lab) in enumerate(zip(inputs, labels)):
        sel = select_at_time(inp, t)
        if sel is None:
            continue
        fn, tt = sel
        d, _time = profile(fn, p)
        xs_all.append((d["x"], tt))
        color = f"C{i % 10}"
        # Binned-mean curve of THIS run (same binning used for the reference), so
        # every panel shows a clean line through its own scatter. The scatter is
        # kept faint (low alpha, behind) so the line is highlighted on top.
        # Bin count tracks the run's own resolution (reg_nbins ~ n_x).
        cx, binned = _binned_multi(d["x"], d, [k for k, _l in keys],
                                   nbins=reg_nbins(inp))
        for k, _lab in keys:
            axd[k].scatter(d["x"], d[k], s=4, alpha=0.12, color=color,
                           rasterized=True, zorder=1,
                           label=f"{lab} (t={tt:.3g})")
            if cx is not None and k in binned:
                axd[k].plot(cx, binned[k], "-", color=color, lw=1.8, zorder=4)

    # Reference overlay: the blessed binned profile on EVERY panel (dashed line).
    if ref_table is not None and "x" in ref_table:
        for k, _lab in keys:
            if k in ref_table:
                axd[k].plot(ref_table["x"], ref_table[k], "--", color="0.35",
                            lw=1.6, zorder=5, label="reference")

    # Exact overlay + L1 (hydro presets only). PERIODIC domain: include the
    # boundary (wrapped) Riemann fan as well as the central one.
    if p["hydro_exact"] and xs_all:
        xx = np.concatenate([x for x, _ in xs_all])
        tt = xs_all[0][1]
        xmin, xmax = float(xx.min()), float(xx.max())
        xg = np.linspace(xmin, xmax, 2000)       # denser: two fans + sharp edges
        rho_a, vx_a, P_a = exact_riemann_periodic(p, tt, xg, xmin, xmax)
        ana = {"rho": rho_a, "vx": vx_a, "P": P_a,
               "u": P_a / ((p["gamma"] - 1.0) * rho_a)}
        for k, _lab in keys:
            if k in ana:
                axd[k].plot(xg, ana[k], "k-", lw=1.6, zorder=6,
                            label="exact (periodic)")
        # L1 against the exact (periodic) solution at each particle, first input.
        sel = select_at_time(inputs[0], t)
        if sel is not None:
            d, tt0 = profile(sel[0], p)
            xmin0, xmax0 = float(d["x"].min()), float(d["x"].max())
            r_a, v_a, P_a = exact_riemann_periodic(p, tt0, d["x"], xmin0, xmax0)
            l1r = l1_error(d["x"], d["rho"], r_a, d["vol"])
            l1v = l1_error(d["x"], d["vx"], v_a, d["vol"])
            l1p = l1_error(d["x"], d["P"], P_a, d["vol"])
            print(f"shocktube[{p['name']}]: t={tt0:.4g}  L1(rho)={l1r:.4g}  "
                  f"L1(vx)={l1v:.4g}  L1(P)={l1p:.4g}")
            abox = dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85)
            axd["rho"].text(0.03, 0.95,
                            f"L1: rho {l1r:.2g}\n     vx {l1v:.2g}\n      P {l1p:.2g}",
                            transform=axd["rho"].transAxes, fontsize=8,
                            va="top", ha="left", bbox=abox)
    elif p["ref"] and xs_all:
        # MHD preset with a tabulated reference (SPLASH exact_mhdshock): overlay
        # the piecewise-constant solution and L1 the first input against it.
        xx = np.concatenate([x for x, _ in xs_all])
        tt = xs_all[0][1]
        rx, rstate = mhd_reference(p["ref"], tt, xx.min(), xx.max())
        for k, _lab in keys:
            if k in rstate:
                axd[k].plot(rx, rstate[k], "k-", lw=1.6, zorder=6,
                            label="reference (RJ95/BW)")
        sel = select_at_time(inputs[0], t)
        if sel is not None:
            d, tt0 = profile(sel[0], p)
            rx0, rstate0 = mhd_reference(p["ref"], tt0, d["x"].min(), d["x"].max())
            l1s = {}
            for k in ("rho", "vx", "P", "By"):
                ref_at = np.interp(d["x"], rx0, rstate0[k])
                l1s[k] = l1_error(d["x"], d[k], ref_at, d["vol"])
            print(f"shocktube[{p['name']}]: t={tt0:.4g}  "
                  + "  ".join(f"L1({k})={v:.4g}" for k, v in l1s.items())
                  + "  (vs SPLASH/RJ95 reference)")
            abox = dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85)
            axd["rho"].text(0.03, 0.95,
                            "L1 vs ref\n" + "\n".join(f"{k} {v:.2g}"
                                                      for k, v in l1s.items()),
                            transform=axd["rho"].transAxes, fontsize=8,
                            va="top", ha="left", bbox=abox)
    elif p["choice"] == 7 and xs_all:
        # C-shock (#7): overlay the Mac-Low et al. (1995) ambipolar C-shock
        # (SPLASH exact_Cshock) using the supplied Mach numbers / shock length.
        co = cshock_opts or {}
        xx = np.concatenate([x for x, _ in xs_all])
        tt = xs_all[0][1]
        xg = np.linspace(xx.min(), xx.max(), 1500)
        cs = cshock_reference(co.get("machs", 50.0), co.get("macha", 5.0),
                              co.get("theta", np.pi / 4.0),
                              co.get("shockl", 0.3), tt, xg)
        for k, _lab in keys:
            if k in cs:
                axd[k].plot(xg, cs[k], "k-", lw=1.6, zorder=6,
                            label="C-shock (Mac-Low 95)")
        print(f"shocktube[{p['name']}]: C-shock overlay "
              f"(machs={co.get('machs', 50.0)}, macha={co.get('macha', 5.0)}, "
              f"shockl={co.get('shockl', 0.3)}). NOTE: matches the Mac-Low "
              f"oblique C-shock only if the run uses that configuration/frame.")
    elif p["choice"] == 8 and xs_all:
        # Steady shock (#8, Falle 2003): the IC satisfies the steady MHD jump
        # conditions, so the exact solution is the stationary step at xshock=0.
        xx = np.concatenate([x for x, _ in xs_all])
        xg = np.linspace(xx.min(), xx.max(), 1000)
        st = steady_step(p, xg)
        for k, _lab in keys:
            if k in st:
                axd[k].plot(xg, st[k], "k-", lw=1.6, zorder=6,
                            label="steady (stationary)")
        sel = select_at_time(inputs[0], t)
        if sel is not None:
            d, tt0 = profile(sel[0], p)
            st0 = steady_step(p, d["x"])
            l1s = {k: l1_error(d["x"], d[k], st0[k], d["vol"])
                   for k in ("rho", "vx", "vy", "By")}
            print(f"shocktube[{p['name']}]: t={tt0:.4g}  "
                  + "  ".join(f"L1({k})={v:.4g}" for k, v in l1s.items())
                  + "  (vs the stationary steady-shock step)")
            abox = dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85)
            axd["rho"].text(0.03, 0.95, "L1 vs steady\n"
                            + "\n".join(f"{k} {v:.2g}" for k, v in l1s.items()),
                            transform=axd["rho"].transAxes, fontsize=8,
                            va="top", ha="left", bbox=abox)
    elif xs_all:
        print(f"shocktube[{p['name']}]: "
              + ("isothermal" if p["gamma"] <= 1.0 + 1e-6 else "MHD")
              + " preset -- no closed-form or tabulated solution; profiles only.")

    for k, lab in keys:
        axd[k].set_ylabel(lab)
    for ax in axes.ravel()[ncol * (nrow - 1):]:
        ax.set_xlabel("x")
    fig.suptitle(f"Shock tube: {p['name']} (choice {p['choice']}, "
                 f"gamma={p['gamma']:.3g})")
    fig.tight_layout()
    finish_figure_with_legend(fig, list(axes.ravel()),
                              save=(f"{save}_profile.png" if save else None))


# --------------------------------------------------------------------------- #
#  Regression baseline (bless / compare the binned density profile at --time)
# --------------------------------------------------------------------------- #

REG_NBINS = 64


def panel_keys(p):
    """The quantities plotted (and now blessed/overlaid) for this preset:
    hydro -> rho,vx,P,u ; MHD -> rho,vx,vy,P,By,Bz (mirrors plot_profiles)."""
    if p["mhd"]:
        return ["rho", "vx", "vy", "P", "By", "Bz"]
    return ["rho", "vx", "P", "u"]


def _binned_multi(x, d, keys, nbins=REG_NBINS):
    """Per-x-bin mean of each `keys` field over the FULL [min,max] x range, via the
    framework `binned_profile` -> (centres, {k: binned}) or (None, None). The
    full-extent grid (not a trimmed percentile window) reaches the periodic
    boundary shock that lives at the box edges, so the blessed reference covers
    them."""
    return binned_profile(x, {k: d[k] for k in keys}, nbins=nbins)


def reference_table(inp, params, t):
    """Regression table {x: x-grid, <q>: binned q(x) for every plotted panel} at
    the comparison time t -- rho,vx,P,u (hydro) or rho,vx,vy,P,By,Bz (MHD).

    Binning every plotted quantity (not just density) lets the saved reference be
    overlaid on ALL profile panels and L1'd per panel. None if no snapshot / too
    few bins."""
    sel = select_at_time(inp, t)
    if sel is None:
        return None
    fn, _tt = sel
    d, _time = profile(fn, params)
    cx, binned = _binned_multi(d["x"], d, panel_keys(params),
                               nbins=reg_nbins(inp))
    if cx is None:
        return None
    table = {"x": cx}
    table.update(binned)
    return table


def compare_to_baseline(inp, params, t, explicit=None, ref_table=None):
    """Print per-panel residuals (every blessed profile field) vs the sibling
    baseline; return the worst RMS (None if no baseline / no profile). `ref_table`
    may be passed pre-loaded to avoid re-discovering it."""
    if ref_table is None:
        ref_table, _p = find_reference(inp, explicit)
    if ref_table is None:
        return None
    table = reference_table(inp, params, t)
    if table is None:
        return None
    worst = 0.0
    for key in panel_keys(params):
        res = residuals(table, ref_table, key)
        if res is None:
            continue
        mx, rms = res
        worst = max(worst, rms)
        print(f"shocktube[regression] {key:<3s} vs baseline: "
              f"max|d|={mx:.4g}  rms={rms:.4g}")
    return worst


def main():
    parser = build_parser(PARAM_SPEC, description="Shock-tube / Riemann analysis.")
    parser.add_argument("--times", type=float, nargs="+", default=DEFAULT_TIME,
                        help="comparison time(s) (default: the final snapshot). "
                             "A single comparison time is just --times <t>.")
    parser.add_argument("--machs", type=float, default=50.0,
                        help="C-shock (#7) sonic Mach number (default 50, "
                             "consistent with the SPLASH vx0=-5, cs=0.1)")
    parser.add_argument("--macha", type=float, default=5.0,
                        help="C-shock (#7) Alfven Mach number (default 5)")
    parser.add_argument("--cshock-shockl", type=float, default=0.3,
                        help="C-shock (#7) ambipolar shock thickness (x-scale; "
                             "default 0.3)")
    parser.add_argument("--cshock-theta", type=float, default=np.pi / 4.0,
                        help="C-shock (#7) field angle in radians (default pi/4)")
    parser.add_argument("--reg-tol", type=float, default=None,
                        help="CI gate: exit non-zero if the worst RMS residual "
                             "vs the regression baseline exceeds this tolerance")
    args = parser.parse_args()
    # --times is the comparison-time list; shocktube compares at a single time
    # (its first entry; None -> the final snapshot).
    args.time = args.times[0] if args.times else None
    raw_params = params_from_args(args, PARAM_SPEC)
    params = shock_params(raw_params["choice"])

    labels = args.labels if args.labels else [
        os.path.basename(os.path.normpath(i)) for i in args.inputs]
    if len(labels) != len(args.inputs):
        parser.error("--labels must match the number of inputs")

    # --save-reference: bless the FIRST input's binned profiles (all panels).
    if args.save_reference:
        table = reference_table(args.inputs[0], params, args.time)
        if table is None:
            print("shocktube: no profile to save as reference",
                  file=sys.stderr)
            sys.exit(1)
        sel0 = select_at_time(args.inputs[0], args.time)
        tref = sel0[1] if sel0 else args.time
        save_reference(table, reference_path_for(args.inputs[0]), raw_params,
                       "shocktube", time=tref)
        return

    overlay = ("exact Riemann" if params["hydro_exact"]
               else f"tabulated ref ({params['ref']})" if params["ref"]
               else "C-shock (Mac-Low 95)" if params["choice"] == 7
               else "steady step (Falle 2003)" if params["choice"] == 8
               else "none (profiles only)")
    print(f"[shocktube] choice {params['choice']} = {params['name']}  "
          f"gamma={params['gamma']:.4g}  mhd={params['mhd']}  overlay={overlay}")

    # Auto-discover (or use --reference) the baseline once; overlay it on the
    # profile panels and reuse it for the first input's residual gate. Pass the
    # resolved comparison time so a time-mismatch vs the blessed reference warns.
    sel0 = select_at_time(args.inputs[0], args.time)
    tcmp = sel0[1] if sel0 else args.time
    ref_table, _refpath = find_reference(args.inputs[0], explicit=args.reference,
                                         time=tcmp)

    cshock_opts = dict(machs=args.machs, macha=args.macha,
                       shockl=args.cshock_shockl, theta=args.cshock_theta)
    plot_profiles(args.inputs, labels, params, args.time, save=args.save,
                  cshock_opts=cshock_opts, ref_table=ref_table)

    worst = None
    for i, inp in enumerate(args.inputs):
        rt = ref_table if i == 0 else None
        w = compare_to_baseline(inp, params, args.time, args.reference,
                                ref_table=rt)
        if w is not None:
            worst = w if worst is None else max(worst, w)
    if args.reg_tol is not None and worst is not None and worst > args.reg_tol:
        print(f"shocktube: REGRESSION FAILED -- worst rms {worst:.4g} > tol "
              f"{args.reg_tol:.4g}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
