#!/usr/bin/env python3
"""
Analysis for the Rayleigh-Taylor (rt) testcase.

A heavy fluid sits above a light fluid in a gravitational field; the interface is
RT-unstable and a seeded velocity perturbation grows exponentially in the linear
regime before rolling up into bubbles/spikes. Extracts, per snapshot, the linear
RT mode amplitude M and the mean vertical kinetic-energy density of the most
energetic 5% of particles (the growing fingers), versus time. The classic
inviscid linear-growth line exp(n t) -- n = sqrt(A g k) -- is overlaid on the
mode-amplitude panel (and exp(2 n t) on the kinetic-energy panel, since
KE ~ vz^2).

WHAT THE MODE-GROWTH METRIC M SHOWS: the IC seeds a vertical-velocity
perturbation vz ~ dvz (1 + cos(8 pi x)) localized about the interface(s). M is
the amplitude of the cos(8 pi x) Fourier component of vz (sin/cos projection,
so it is phase-agnostic), mass-weighted and LOCALIZED to the interfaces by the
RT eigenmode's own vertical envelope exp(-k |z' - z_iface|) (z' = |z| for the
symmetric case 1, which folds the two interfaces at +-z_iface coherently --
their seeded eigenmodes are mirror images, so both add in phase). In the linear
regime this is precisely the amplitude of the seeded RT eigenmode and grows as
exp(n t); it says nothing about the constant (k=0) part of the seed or about
other modes until nonlinear coupling feeds them.

Geometry/physics fixed by IC_setup_rt.create():
  * horizontal seed mode  vz ~ (1 + cos(2 C x)),  C = 4 pi  ->  k = 8 pi
  * interface at |z| = dICdensR * nlowdens = 0.25 * nlowdens  (case 1: two
    symmetric interfaces at +-z_iface, light |z|<z_iface / heavy outside, gravity
    toward the midplane; case 2: a single interface at z = z_iface, heavy on top)
  * Atwood number  A = (rho_heavy - rho_light)/(rho_heavy + rho_light)
                     = (rhodiff - 1)/(rhodiff + 1)        (rho_light = rhozero = 1)
  * gravity magnitude g = 1.0, FIXED: gasoline applies the constant body force
    a_z = -dBodyForceConst sign(z) with dBodyForceConst = 1.0 (written by
    IC_createparamfile; these are specific test cases, and changing g is just
    a time rescale t -> t/sqrt(g), so it is deliberately not a parameter). The
    IC hydrostatic pressure profile is integrated with the SAME g = 1.0.
    (Before 2026-06-09 the profile hardcoded g=0.5 against the param file's
    1.0 -- only half the run's gravity was supported, which made the whole
    column free-fall/compress.)
  * MHD (beta > 0, -e): uniform B with strength from the plasma beta at the
    interface, B0 = sqrt(2 p(z_iface)/beta); geometry -f bdir: 1 = x horizontal
    (default), 2 = y, 3 = z vertical.

Parameters mirror runtest.sh setup_rt() flag-for-flag with the same defaults, so
you declare the run the same way you launched it, e.g.

    python IC_analysis_rt.py test_cases/rt/rt128_..._GASOLINE -a 2.0 -b 1 -c 2

Reference (regression-baseline) workflow + the density --render are handled by
IC_analysis_framework (run_cli), exactly as for kh.

Usage:
    python IC_analysis_rt.py <dir-or-prefix> [-a..-d ...] [--save out] [--render]
"""

import numpy as np

from IC_analysis_framework import (
    Param, Metric, RenderSpec, StreamSpec, run_cli,
    aux_component, aux_magnitude, aux_divB_error,
)
from IC_analysis_general import amp_mag


# -a..-f match runtest.sh setup_rt() flag-for-flag with the same defaults.
PARAM_SPEC = [
    Param("a", "rhodiff",  2.0, float, "density contrast (rho_heavy/rho_light)"),
    Param("b", "case",     1,   int,   "1 = symmetric (two interfaces), 2 = single interface"),
    Param("c", "nlowdens", 2,   float, "number of low-density segments (sets z_iface = 0.25*nlowdens)"),
    Param("d", "smooth",   0,   float, "smooth (tanh-ramped) interface density"),
    Param("e", "beta",     0.0, float, "plasma beta at the interface (0 = hydro)"),
    Param("f", "bdir",     1,   int,   "B direction: 1=x horizontal (default), 2=y, 3=z vertical"),
]

# Plotted quantities returned by analyze().
METRICS = [
    Metric("M",   "Linear RT mode amplitude",            "t",    "log", "o"),
    Metric("Evz", "Max vertical kinetic-energy density", "time", "log", "o"),
]

# Geometric/physical constants of the RT IC (fixed by setup_rt, not -a..-f):
#   KX   horizontal wavenumber of the seeded mode, cos(2*C*x) with C = 4 pi.
#   GRAV gravity magnitude: the run's bBodyForce dBodyForceConst = 1.0 (fixed
#        in IC_createparamfile), which the IC hydrostatic profile matches.
#   RDENS dICdensR = 0.25 (sets the interface position z_iface = RDENS * nlowdens).
KX = 8.0 * np.pi
GRAV = 1.0
RDENS = 0.25

# Seed (t=0) anchors, set from the IC perturbation amplitude (measured in memory
# from setup_rt; see the module test), so the measured growth curve and the
# analytic exp(n t) / exp(2 n t) lines start where the IC does -- mirrors kh.
SEED_M = 0.043
SEED_EVZ = 7.0e-3
SEED = {"x": 0.0, "M": SEED_M, "Evz": SEED_EVZ}


def _z_iface(params):
    """Interface position z_iface = dICdensR * nlowdens (IC_setup_rt)."""
    return RDENS * float(params.get("nlowdens", 2.0))


def growth_rate(params):
    """Classic inviscid RT linear growth rate n = sqrt(A g k).

    A = (rhodiff-1)/(rhodiff+1) is the Atwood number (rho_light = rhozero = 1)."""
    rhodiff = float(params.get("rhodiff", 2.0))
    atwood = (rhodiff - 1.0) / (rhodiff + 1.0)
    return np.sqrt(max(atwood, 0.0) * GRAV * KX)


def analyze_rt(tgdata, params):
    """Mode amplitude M and top-5% mean vertical KE density from one snapshot.

    M is the modal amplitude of vz onto the seeded horizontal mode cos(KX x),
    mass-weighted and localized about the interface by exp(-KX |zterm - z_iface|)
    (zterm = |z| for the symmetric case 1, folding both interfaces coherently;
    z for case 2). This is the RT analogue of analyze_kh's shear-layer projection.
    """
    case = int(params.get("case", 1))
    z_iface = _z_iface(params)

    mass = tgdata[:, 0]
    x = tgdata[:, 1]
    z = tgdata[:, 3]
    vz = tgdata[:, 6]
    rho = tgdata[:, 7]

    zterm = np.abs(z) if case == 1 else z
    weight = mass * np.exp(-KX * np.abs(zterm - z_iface))
    M = amp_mag(weight, x, vz, KX, 0.0)

    # Vertical kinetic-energy density (mirrors analyze_kh's top-5% Ekin metric,
    # using vz -- the RT-growing component); grows as exp(2 n t) in the linear regime.
    Evz = 0.5 * rho * vz ** 2
    threshold = np.percentile(Evz, 95)
    return M, float(Evz[Evz >= threshold].mean())


def analyze(tgdata, time, params, snap):
    """Framework callback: one row of metrics for a snapshot."""
    M, Evz = analyze_rt(tgdata, params)
    return {"x": float(np.asarray(time).ravel()[0]), "M": M, "Evz": Evz}


def analytic(metric_key, params):
    """Analytic overlay: linear RT growth. exp(n t) for the mode amplitude M,
    exp(2 n t) for the vertical KE density (KE ~ vz^2). Anchored at the IC seed
    amplitude over the early (linear) window, before bubbles/spikes go nonlinear."""
    n = growth_rate(params)
    if n <= 0.0:
        return None                       # no instability (rhodiff <= 1)
    # Linear window: from t=0 up to ~3 e-foldings, where exponential growth holds
    # before the nonlinear roll-up.
    t = np.linspace(0.0, 3.0 / n, 200)
    if metric_key == "M":
        return t, SEED_M * np.exp(n * t)
    if metric_key == "Evz":
        return t, SEED_EVZ * np.exp(2.0 * n * t)
    return None


# Field map for --render: density in the x-z plane (axis1=1 x, axis2=3 z; gravity
# is along z). The slab is thin in y and the RT structure is y-invariant, so the
# natural projection is the plain line-of-sight average along y ("average").
# Density spans rho_light=rhozero=1 .. rho_heavy=rhodiff, declared as the natural
# color range for field-valued projections (None for an integrated "column").
# The render DOMAIN is the particle region, declared as the spec's default
# extent: the run's dzPeriod is 2x the particle extent (setup_rt doubles dzbound
# to pad the periodic box with vacuum), so the .log box would draw the slab into
# half-empty panels (the "renders look tiny" issue).
def _domain(params):
    """Particle-region extent [x0,x1,z0,z1] from the setup_rt geometry."""
    rhodiff = float(params.get("rhodiff", 2.0))
    nlow = float(params.get("nlowdens", 2.0))
    if int(params.get("case", 1)) == 1:
        dz = 0.5 + RDENS * (nlow - 1.0)
        fb = 0.0                       # forceboundary = 0 identically for case 1
    else:
        dz = 0.5 + RDENS * 0.5 * (rhodiff - 1.0) + RDENS * (nlow - 1.0)
        fb = dz - RDENS - RDENS * nlow
    return (-0.25, 0.25, -dz - fb, dz - fb)


def render_specs(params):
    def clim(p, proj):
        if proj in ("slice", "average", "rhocolumn"):
            return (1.0, float(p.get("rhodiff", 2.0)))
        return None
    ext = _domain(params)
    specs = [RenderSpec(axis1=1, axis2=3, quantity=7, label="density",
                        log=False, project="average", clim=clim, extent=ext)]
    if params.get("beta", 0.0) > 0.0:
        bdiv = dict(log=False, project="rhocolumn", cmap="RdBu_r",
                    symmetric=True, extent=ext)
        specs += [
            RenderSpec(1, 3, aux_magnitude("BField"), "|B|",
                       log=True, project="rhocolumn", extent=ext),
            RenderSpec(1, 3, aux_component("BField", 0), "Bx", **bdiv),
            RenderSpec(1, 3, aux_component("BField", 2), "Bz", **bdiv),
            RenderSpec(1, 3, aux_divB_error("DivB", "BField"), "divBerr",
                       log=True, project="rhocolumn", extent=ext,
                       clim=lambda params, proj: (0.001, 0.1)),
        ]
    return specs


# Streamlines drawn with --render (framework streamline_panels): the in-plane
# (x, z) velocity (vx col 4, vz col 6 -- bubbles/spikes circulation) always, and
# the in-plane B (Bx, Bz from the BField aux) for magnetised (beta>0) runs.
def stream_specs(params):
    ext = _domain(params)
    specs = [StreamSpec(1, 3, 4, 6, "velocity", extent=ext)]
    if params.get("beta", 0.0) > 0.0:
        specs.append(StreamSpec(1, 3, aux_component("BField", 0),
                                aux_component("BField", 2), "B", extent=ext))
    return specs


def main():
    run_cli(
        test="rt",
        param_spec=PARAM_SPEC,
        analyze=analyze,
        analytic=analytic,
        metrics=METRICS,
        description="Rayleigh-Taylor mode-amplitude / growth-rate analysis vs time.",
        seed=SEED,
        render_spec=render_specs,
        stream_spec=stream_specs,
    )


if __name__ == "__main__":
    main()
