#!/usr/bin/env python3
"""
Analysis for the Taylor-Green vortex (taylorgreen) testcase.

The 2-D Taylor-Green vortex (IC_setup_taylorgreen) on the unit periodic box is
    vx =  v0 sin(2 pi x) cos(2 pi y)
    vy = -v0 sin(2 pi y) cos(2 pi x) ,   vz = 0,   rho = rhozero = 1
which is a Laplacian eigenmode of wavenumber k = 2 pi. Under Navier-Stokes shear
viscosity nu = eta/rho the field keeps its spatial form and decays exactly:
    v(t)   = v(0) exp(-2 nu k^2 t)
    E_k(t) = E_k(0) exp(-4 nu k^2 t)            (energy ~ v^2)
This run is launched WITH a physical shear viscosity (runtest -a etavisc ->
-PhysViscEta), so exp(-4 nu k^2 t) is the exact analytic; the artificial
viscosity adds EXCESS dissipation, so the measured decay running faster than the
analytic is the numerical-dissipation diagnostic. The TG flow is divergence-free,
so the physical bulk viscosity xi (-b) does not act on it.

Extracts, per snapshot, the kinetic energy and the amplitude of the TG velocity
eigenmode (projection of vx onto sin(2 pi x) cos(2 pi y)), each NORMALIZED by its
t=0 value (so the comparison is resolution-independent), versus time. The analytic
exp(-4 nu k^2 t) / exp(-2 nu k^2 t) decay is overlaid on a log y-axis.

Parameters mirror runtest.sh setup_taylorgreen() flag-for-flag with the same
defaults; you declare the run the same way you launched it, e.g.

    python IC_analysis_taylorgreen.py test_cases/taylorgreen/..._GASOLINE -a 0.1

Reference (regression-baseline) workflow + the |v| --render are handled by
IC_analysis_framework (run_cli), exactly as for kh.

Usage:
    python IC_analysis_taylorgreen.py <dir-or-prefix> [-a..-d ...] [--render]
"""

import numpy as np

from IC_analysis_framework import (
    Param, Metric, RenderSpec, run_cli,
)


# -a..-d match runtest.sh setup_taylorgreen() flag-for-flag with the same defaults.
PARAM_SPEC = [
    Param("a", "etavisc", 0.1,  float, "physical shear viscosity eta (-PhysViscEta); nu = eta/rho"),
    Param("b", "xivisc",  0.0,  float, "physical bulk viscosity xi (no effect on the div-free TG flow)"),
]

# Both metrics are normalized by their own t=0 value -> the analytic decay starts
# at 1 and the comparison is resolution-independent (energy scales with the thin
# z-slab volume otherwise).
METRICS = [
    Metric("Ekin",  r"$E_{kin}/E_{kin,0}$", "t", "log", "o", norm="Ekin"),
    Metric("TGamp", r"$A_v/A_{v,0}$",       "t", "log", "s", norm="TGamp"),
]

# Constants fixed by IC_setup_taylorgreen (not -a..-d parameters):
KTG = 2.0 * np.pi      # TG eigenmode wavenumber (velocity ~ sin(2 pi x) cos(2 pi y))
RHOZERO = 1.0          # uniform background density (nu = eta / RHOZERO)
V0 = 0.1               # IC velocity amplitude (create() default; for the render clim)
T_END = 2.0            # nominal run end = nsteps*deltastep = 200*0.01 (analytic span)


def _nu(params):
    """Kinematic shear viscosity nu = eta / rho."""
    return float(params.get("etavisc", 0.1)) / RHOZERO


def analyze_tg(tgdata):
    """Kinetic energy and TG eigenmode amplitude from one snapshot.

    TGamp projects vx onto the eigenmode S = sin(2 pi x) cos(2 pi y):
    amp = sum(m vx S) / sum(m S^2) = v0 for the exact field (a mass-weighted,
    resolution-independent recovery of the eigenmode amplitude)."""
    m = tgdata[:, 0]
    x = tgdata[:, 1]
    y = tgdata[:, 2]
    vx = tgdata[:, 4]
    v2 = np.sum(tgdata[:, 4:7] ** 2, axis=1)
    ekin = float(np.sum(0.5 * m * v2))

    S = np.sin(KTG * x) * np.cos(KTG * y)
    denom = float(np.sum(m * S * S))
    tgamp = float(np.sum(m * vx * S) / denom) if denom > 0 else 0.0
    return ekin, abs(tgamp)


def analyze(tgdata, time, params, snap):
    """Framework callback: one row of metrics for a snapshot."""
    ekin, tgamp = analyze_tg(tgdata)
    return {"x": float(np.asarray(time).ravel()[0]), "Ekin": ekin, "TGamp": tgamp}


def analytic(metric_key, params):
    """Analytic viscous decay (normalized to t=0): exp(-4 nu k^2 t) for the
    energy, exp(-2 nu k^2 t) for the velocity amplitude (energy ~ v^2)."""
    nu = _nu(params)
    rate_v = 2.0 * nu * KTG ** 2          # velocity-amplitude decay rate
    t = np.linspace(0.0, T_END, 200)
    if metric_key == "Ekin":
        return t, np.exp(-2.0 * rate_v * t)   # energy decays at twice the v rate
    if metric_key == "TGamp":
        return t, np.exp(-rate_v * t)
    return None


# Field map for --render: speed |v| = sqrt(vx^2 + vy^2) in the x-y plane (the TG
# four-vortex pattern). The slab is thin and the flow is z-invariant, so the plain
# line-of-sight average along z ("average") is the natural projection. |v| spans
# 0 .. v0 (declared as the natural color range for field-valued projections).
def _speed(fn, tgdata, ngas):
    return np.sqrt(tgdata[:, 4] ** 2 + tgdata[:, 5] ** 2)


def render_specs(params):
    def clim(p, proj):
        if proj in ("slice", "average", "rhocolumn"):
            return (0.0, V0)
        return None
    return RenderSpec(axis1=1, axis2=2, quantity=_speed, label="|v|",
                      log=False, project="average", clim=clim, cmap="viridis")


def main():
    run_cli(
        test="taylorgreen",
        param_spec=PARAM_SPEC,
        analyze=analyze,
        analytic=analytic,
        metrics=METRICS,
        description="Taylor-Green vortex viscous-decay analysis vs time.",
        render_spec=render_specs,
    )


if __name__ == "__main__":
    main()
