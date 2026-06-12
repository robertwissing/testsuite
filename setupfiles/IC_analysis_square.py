#!/usr/bin/env python3
"""
Analysis for the hydrostatic-square test (square; Saitoh & Makino 2013, Hopkins
2013/2015, Hu et al. 2014; 3D cube = Sandnes 2025).

A high-density square (rho_in = rhodiff, default 4) sits in a low-density medium
(rho_out = rhozero = 1) at UNIFORM pressure P = przero with zero velocity. With a
constant pressure and no velocity there are no forces, so the exact solution is
STATIC -- the square stays a sharp square forever. Standard SPH instead develops
spurious "surface tension" at the contact discontinuity, which rounds the corners
and drives spurious flows; this test measures that departure from the static IC.

Two complementary diagnostics:
  * `--render` density in the image plane (the canonical visual: does the square
    stay square or round off / blob up?).
  * spurious-velocity Mach number vs time -- the RMS and maximum |v| over the
    medium sound speed c_s = sqrt(gamma P / rho_out). The exact solution is
    static, so the analytic overlay is zero; any rise is pure numerical error.

Parameters mirror runtest.sh setup_square() flag-for-flag with the same defaults,
so you declare the run the same way you launched it, e.g.

    python IC_analysis_square.py test_cases/square/..._GASOLINE -a 4.0 -b 2.5

Reference (regression-baseline) workflow + the density --render are handled by
IC_analysis_framework (run_cli), exactly as for kh.

Usage:
    python IC_analysis_square.py <dir-or-prefix> [-a..-f ...] [--render]
"""

import numpy as np

from IC_analysis_framework import (
    Param, Metric, RenderSpec, run_cli,
)


# -a..-f match runtest.sh setup_square() flag-for-flag with the same defaults.
PARAM_SPEC = [
    Param("a", "rhodiff", 4.0,           float, "density contrast chi = rho_in/rho_out"),
    Param("b", "przero",  2.5,           float, "uniform pressure P (sets the sound speed)"),
    Param("c", "gamma",   5.0 / 3.0,     float, "adiabatic index"),
    Param("d", "is3d",    0,             int,   "1 = 3D inner cube, 0 = 2D thin-slab square"),
]

# Spurious-velocity diagnostics: ideal (static) = 0, so the analytic is a flat
# zero line and any growth is the spurious surface-tension flow. Linear y-axis.
METRICS = [
    Metric("mach_rms", r"RMS spurious Mach $\langle|v|^2\rangle^{1/2}/c_s$", "t", "linear", "o"),
    Metric("mach_max", r"max spurious Mach $|v|_{max}/c_s$",                 "t", "linear", "s"),
]

# Constants fixed by IC_setup_square (not -a..-f): the outer/background density,
# and the nominal run end = nsteps*deltastep = 200*0.02 (span of the analytic line).
RHOZERO = 1.0
T_END = 4.0


def sound_speed(params):
    """Medium sound speed c_s = sqrt(gamma P / rho_out) (rho_out = rhozero = 1).

    Used as the Mach-number reference; the ambient medium fills most of the box
    and is the faster of the two (uniform P, lower rho)."""
    gamma = float(params.get("gamma", 5.0 / 3.0))
    przero = float(params.get("przero", 2.5))
    return np.sqrt(gamma * przero / RHOZERO)


def analyze(tgdata, time, params, snap):
    """Framework callback: RMS and maximum spurious Mach number for a snapshot."""
    m = tgdata[:, 0]
    v2 = np.sum(tgdata[:, 4:7] ** 2, axis=1)
    cs = sound_speed(params)
    vrms = np.sqrt(np.sum(m * v2) / np.sum(m))
    vmax = np.sqrt(np.max(v2))
    return {"x": float(np.asarray(time).ravel()[0]),
            "mach_rms": float(vrms / cs), "mach_max": float(vmax / cs)}


def analytic(metric_key, params):
    """Analytic overlay: the exact solution is static, so both spurious-Mach
    metrics are identically zero."""
    t = np.linspace(0.0, T_END, 2)
    return t, np.zeros_like(t)


# Field map for --render: density in the image plane. For the 2D thin-slab square
# the structure is z-invariant -> plain LOS average along z ("average"); for the
# 3D cube average would blur the faces, so take the mid-plane slice. Density spans
# rho_out=rhozero=1 .. rho_in=rhodiff, the natural color range for field-valued
# projections (None for an integrated "column").
def render_specs(params):
    proj = "slice" if int(params.get("is3d", 0)) else "average"

    def clim(p, projection):
        if projection in ("slice", "average", "rhocolumn"):
            return (RHOZERO, float(p.get("rhodiff", 4.0)))
        return None

    return RenderSpec(axis1=1, axis2=2, quantity=7, label="density",
                      log=False, project=proj, clim=clim)


def main():
    run_cli(
        test="square",
        param_spec=PARAM_SPEC,
        analyze=analyze,
        analytic=analytic,
        metrics=METRICS,
        description="Hydrostatic-square (Saitoh-Makino) spurious-velocity analysis.",
        render_spec=render_specs,
    )


if __name__ == "__main__":
    main()
