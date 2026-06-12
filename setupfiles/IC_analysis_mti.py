#!/usr/bin/env python3
"""
Analysis for the magnetothermal instability (MTI) test (IC_setup_mti).

This is the MTI test of Owens et al. (2026, "Fast, stable, and physical:
hyperbolic, magnetic field-aligned diffusion in SPH"; the gasoline-lineage
anisotropic-conduction paper) -- itself the Parrish & Stone (2005)/Balbus (2000)
MTI. A stratified atmosphere has temperature decreasing in the direction of
gravity over an unstable layer (here u = u0(1 - (z-b1)/H), rho ~ (1-(z-b1)/H)^2,
H = 3, isothermal boundary layers outside), threaded by a weak HORIZONTAL field
(Bx = B0, beta ~ 2e8). Anisotropic conduction along the field converts the
buoyantly-unstable temperature gradient into convective motions: in the linear
phase the kinetic energy grows as exp(2 t / t_buoy) (Balbus 2000), saturating
around 10 t_buoy, while the convection bends the initially horizontal field
toward vertical.

The buoyancy time is t_buoy = |g d(ln T)/dz|^(-1/2); with g = 1 and
|d ln T/dz| = 1/H = 1/3 (the unstable-layer profile), t_buoy = sqrt(3). The seed
is delta vz = v0 sin(4 pi x / L), v0 = machno c_s.

Following Owens et al., the diagnostics (vs t/t_buoy) are:
  * the kinetic energy E_kin(t)/E_kin(0) -- exponential growth exp(2 t/t_buoy)
    in the linear phase, then saturation (semilogy, with the analytic overlay),
  * the fraction of magnetic energy in the vertical direction Bz^2/B^2 -- it
    starts at 0 (purely horizontal field) and rises toward ~0.5 (x-z isotropy)
    as the convection reorients the field,
  * field-map renders (--render): temperature, vertical velocity vz, and the
    composition tracer in the x-z plane (the convective plumes / layer mixing).

Built on run_cli (metric-vs-time + analytic + reference + render), like kh/rt.

Usage:
    python IC_analysis_mti.py <run-dir> [--save out] [--render]
"""

import numpy as np

from IC_analysis_framework import (
    Param, Metric, RenderSpec, run_cli,
)
from IC_analysis_general import read_vector_aux

# Constants fixed by IC_setup_mti.create() (not parameters):
#   t_buoy = |g d(ln T)/dz|^(-1/2) = sqrt(H/g); g = 1, H = 3 -> sqrt(3).
GRAV = 1.0
H_SCALE = 3.0
T_BUOY = np.sqrt(H_SCALE / GRAV)

# runtest.sh setup_mti() passes NO parameters (machno uses the create() default).
# Carried for documentation only; the analysis self-normalizes so the seed
# amplitude (machno) does not affect the growth-rate comparison.
PARAM_SPEC = [
    Param("a", "machno", 0.01, float, "seed Mach number v0=machno*c_s (runtest passes none; default 0.01)"),
]

# E_kin normalized to its t=0 value (resolution/machno-independent); Bz^2/B^2 is
# already a fraction (no normalization). Time axis is t/t_buoy.
METRICS = [
    Metric("KE",     r"$E_{kin}/E_{kin,0}$", r"$t/t_{buoy}$", "log",    "o", norm="KE"),
    Metric("Bzfrac", r"$B_z^2/B^2$",         r"$t/t_{buoy}$", "linear", "s"),
]


def analyze(tgdata, time, params, snap):
    """Framework callback: kinetic energy and vertical magnetic-energy fraction.

    The time is reported in buoyancy times (x = t/t_buoy). Bz^2/B^2 is the
    volume-weighted fraction of magnetic energy in z (read from the BField aux);
    NaN if the run is non-magnetic / the aux is absent."""
    m = tgdata[:, 0]
    v2 = np.sum(tgdata[:, 4:7] ** 2, axis=1)
    ke = float(np.sum(0.5 * m * v2))

    ngas = len(tgdata)
    B = read_vector_aux(snap, "BField", ngas)
    if B is None:
        bzfrac = np.nan
    else:
        B = np.asarray(B, float)
        vol = m / tgdata[:, 7]
        b2 = np.sum(B ** 2, axis=1)
        denom = float(np.sum(b2 * vol))
        bzfrac = float(np.sum(B[:, 2] ** 2 * vol) / denom) if denom > 0 else np.nan

    return {"x": float(np.asarray(time).ravel()[0]) / T_BUOY,
            "KE": ke, "Bzfrac": bzfrac}


def analytic(metric_key, params):
    """Linear-MTI overlays (x = t/t_buoy): exp(2 t/t_buoy) growth for the kinetic
    energy (Balbus 2000), and the x-z isotropic saturation level 0.5 for Bz^2/B^2."""
    if metric_key == "KE":
        x = np.linspace(0.0, 5.0, 200)        # linear phase (saturation ~10 t_buoy)
        return x, np.exp(2.0 * x)
    if metric_key == "Bzfrac":
        return np.array([0.0, 12.0]), np.array([0.5, 0.5])
    return None


# Field maps for --render: the x-z plane (axis1=1 x, axis2=3 z; stratification is
# along z). The slab is thin in y and the problem is y-invariant -> "average".
def render_specs(params):
    return [
        RenderSpec(1, 3, 8, "T", False, "average", None, "inferno", False),
        RenderSpec(1, 3, 6, "vz", False, "average", None, "RdBu_r", True),
        RenderSpec(1, 3, 10, "tracer", False, "average", None, "viridis", False),
    ]


def main():
    run_cli(
        test="mti",
        param_spec=PARAM_SPEC,
        analyze=analyze,
        analytic=analytic,
        metrics=METRICS,
        description="Magnetothermal instability (Owens 2026 / Parrish-Stone) analysis.",
        render_spec=render_specs,
    )


if __name__ == "__main__":
    main()
