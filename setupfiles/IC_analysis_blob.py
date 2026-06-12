#!/usr/bin/env python3
"""
Analysis for the cloud-crushing (blob) testcase.

Per snapshot it measures the surviving dense-cloud mass and the mass in the
mixed-temperature phase, versus time in cloud-crushing times t_crush. The two
masses are plotted normalised by the initial cloud mass (M/M_cloud,0), the
standard cloud-crushing diagnostic.

Parameters mirror runtest.sh setup_blob() flag-for-flag with the same defaults,
so the physical scales (cloud pressure, wind speed, t_crush, mixing temperature)
are DERIVED from the run parameters -- not hardcoded. In particular the cloud
pressure follows IC_setup_blob.setphysunits():
    przero = clouddens * Tcloud * dkpcunit*KPCCM*KBOLTZ
             / (molweight * dmsolunit * MSOLG*MHYDR*GCGS)
with dmsolunit = rhophyswind*(dkpcunit*KPCCM)^3/MSOLG and rhophyswind=10^-winddens.

blob has no analytic solution, so validation is against a saved reference
baseline (see IC_analysis_framework / runtest.sh -r).

Usage:
    python IC_analysis_blob.py <run-dir> [-a..-h ...] [--save out] [--render]
"""

import numpy as np

from IC_analysis_framework import Param, Metric, RenderSpec, run_cli

# Geometric / physical constants fixed by IC_setup_blob (not -a..-h params).
GAMMA = 5.0 / 3.0
R_IN = 0.1            # cloud radius (IC_setup_blob.Rin)
RHO_OUT = 1.0         # ambient (wind) density in code units
MOLWEIGHT = 0.59259   # IC_setup_blob.molweight
DKPCUNIT = 1.0        # IC_setup_blob.dkpcunit

# -a..-h match runtest.sh setup_blob() (runtest.sh:474-494) with the same defaults.
PARAM_SPEC = [
    Param("a", "rhodiff",  10.0,  float, "density contrast (cloud/wind)"),
    Param("b", "beta",     0.0,   float, "plasma beta (unused in hydro blob)"),
    Param("c", "mach",     2.7,   float, "wind mach number"),
    Param("d", "inflow",   0.0,   float, "inflow conditions"),
    Param("e", "cool",     0.0,   float, "with cooling"),
    Param("f", "windref",  0.0,   float, "reference frame of wind"),
    Param("g", "winddens", 26.0,  float, "wind density 10^-g g/cm^3"),
    Param("h", "tcloud",   4000.0, float, "cloud temperature [K]"),
]

# Both masses are normalised by the initial cloud mass (M_cloud at t=0).
METRICS = [
    Metric("mass_cloud", r"$M_{cloud}/M_{cloud,0}$", r"$t/t_{crush}$",
           "linear", "o", norm="mass_cloud"),
    Metric("mass_mix",   r"$M_{mix}/M_{cloud,0}$",   r"$t/t_{crush}$",
           "linear", "s", norm="mass_cloud"),
]

def render_clim(params, proj):
    """Natural density color range for the blob render.

    The density field runs from the ambient medium (rho=RHO_OUT) up to the cloud
    (rho=rhodiff*RHO_OUT), so for a field-valued projection those are the natural
    limits. A `column` projection integrates density along the line of sight, so
    its units differ -> return None and let the framework auto-scale instead.
    """
    if proj == "column":
        return None
    return (RHO_OUT, params["rhodiff"] * RHO_OUT)


# Field maps: density in the x-z plane (wind blows along z); slice through y.
# clim gives the field's natural color range (medium..cloud density).
RENDER = RenderSpec(axis1=1, axis2=3, quantity=7, label="density",
                    log=True, project="slice", clim=render_clim)


def blob_pressure(clouddens, tcloud, winddens):
    """Cloud pressure, exactly as IC_setup_blob.setphysunits() computes przero."""
    MSOLG = 1.99e33
    KBOLTZ = 1.38e-16
    MHYDR = 1.67e-24
    KPCCM = 3.085678e21
    GCGS = 6.67e-8
    rhophyswind = 10.0 ** (-winddens)
    dmsolunit = rhophyswind * (DKPCUNIT * KPCCM) ** 3 / MSOLG
    return (clouddens * tcloud * DKPCUNIT * KPCCM * KBOLTZ
            / (MOLWEIGHT * dmsolunit * MSOLG * MHYDR * GCGS))


def t_crush(params):
    """Cloud-crushing time from the run params (constant over the run)."""
    clouddens = params["rhodiff"]
    przero = blob_pressure(clouddens, params["tcloud"], params["winddens"])
    v_wind = params["mach"] * np.sqrt(GAMMA * przero / RHO_OUT)
    return R_IN * np.sqrt(clouddens / RHO_OUT) / v_wind


def time_to_x(raw_time, params):
    """Map a snapshot's raw time to t/t_crush (used for --render-times)."""
    return float(np.asarray(raw_time).ravel()[0]) / t_crush(params)


def analyze(tgdata, time, params, snap):
    """Framework callback: surviving cloud mass and mixed-phase mass."""
    clouddens = params["rhodiff"]
    tcloud = params["tcloud"]
    tcrush = t_crush(params)

    rho = tgdata[:, 7]
    mass_cloud = float(np.sum(tgdata[rho > clouddens / 3.0, 0]))

    # Mixed phase: geometric-mean temperature of cloud (T~Tcloud) and pressure-
    # matched wind (T~Tcloud*clouddens) => T_mix = Tcloud*sqrt(clouddens), within
    # a +-1/4 log10(clouddens) band.
    T_mix = tcloud * np.sqrt(clouddens)
    logT = np.log10(tgdata[:, 8])
    lo = np.log10(T_mix) - 0.25 * np.log10(clouddens)
    hi = np.log10(T_mix) + 0.25 * np.log10(clouddens)
    mass_mix = float(np.sum(tgdata[(logT > lo) & (logT < hi), 0]))

    t = float(np.asarray(time).ravel()[0])
    return {"x": t / tcrush, "mass_cloud": mass_cloud, "mass_mix": mass_mix}


def main():
    run_cli(
        test="blob",
        param_spec=PARAM_SPEC,
        analyze=analyze,
        analytic=None,            # no analytic solution; validate vs reference
        metrics=METRICS,
        description="Cloud-crushing (blob) mass-survival analysis vs t_crush.",
        render_spec=RENDER,
        time_to_x=time_to_x,
        time_label="t/t_crush",
    )


if __name__ == "__main__":
    main()
