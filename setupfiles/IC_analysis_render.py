#!/usr/bin/env python3
"""
Field-map rendering for the analysis framework: 2-D field maps, streamlines,
movies, the sph_interp/pynbody/hist/particles backends, render-grid baselines,
the RenderSpec/StreamSpec/RENDER_STYLE style table and the aux/quantity
factories. Extracted from IC_analysis_framework so the framework module is
smaller; IC_analysis_framework re-exports every public name here, so callers
keep doing `from IC_analysis_framework import RenderSpec, render_panels, ...`.

Self-contained: depends only on IC_analysis_general (+ readtipsy, numpy, and --
lazily, inside the functions that need them -- matplotlib and sph_interp). It
does NOT import the reference or framework modules (no import cycle).
"""

import os
import re
import sys
from collections import namedtuple
from functools import lru_cache

import numpy as np

import readtipsy as tip
from IC_analysis_general import (
    loaddata, find_files, setup_rcparams, finish_figure,
    try_aux, read_vector_aux, read_bfield, read_log_param, read_tufac,
    _find_logs,
)

@lru_cache(maxsize=128)
def _cached_readtipsy(fn):
    """Memoized tip.readtipsy: the render path reads each snapshot once per
    (plane, quantity) panel otherwise -- e.g. the orzag 3-D render touches the
    same snapshot 12x (3 planes x 4 fields). Snapshot files are immutable during
    an analysis, so cache by path. The returned arrays are marked read-only so an
    accidental in-place write to the shared cache raises instead of corrupting
    other panels (the render path only reads them)."""
    out = tip.readtipsy(fn)
    for a in out:
        if isinstance(a, np.ndarray):
            a.flags.writeable = False
    return out


# How to render a 2-D field map for a test: which tgdata columns form the image
# plane (axis1, axis2) and which column is the field (quantity), its colorbar
# label, and whether to use a log color scale. tgdata columns:
#   0 mass 1 x 2 y 3 z 4 vx 5 vy 6 vz 7 rho 8 T 9 eps 10 metals 11 potential
#   project: default line-of-sight reduction for this field --
#     "slice"     value in a thin slab through the mid-plane,
#     "column"    integral of the field along the 3rd axis (surface density-like),
#     "rhocolumn" density-weighted column average (mass-weighted mean along LOS),
#     "average"   plain (volume-weighted) average along LOS.
# `quantity` is either a tgdata COLUMN INDEX (int) or a CALLABLE
# quantity(fn, tgdata, ngas) -> per-particle scalar array, for fields not stored
# as a tgdata column (e.g. magnetic-field components/magnitude read from the aux
# file). `clim` is an optional per-spec callback clim(params, proj) ->
# (vmin, vmax) or None, giving this quantity's natural color range (e.g. blob
# density: medium..cloud); None -> percentile scale. `cmap` is the colormap and
# `symmetric` forces a zero-centered color range (for signed fields like B
# components, with a diverging cmap). Each rendered field carries its own
# label/log/cmap/range, so a test can emit several (density, |B|, Bx, ...).
# `extent` is an optional per-test DEFAULT domain [x0, x1, y0, y1] used when no
# --render-extent is given, with priority manual flag > spec.extent > .log box >
# percentiles -- for tests whose particles fill only part of the periodic box
# (e.g. rt's vacuum-padded dzPeriod) or that want a zoomed view.
# `contour` is an optional ContourSpec: iso-contours of a SECOND quantity
# overlaid on every panel (e.g. |B| contours over a density map).
RenderSpec = namedtuple(
    "RenderSpec",
    ["axis1", "axis2", "quantity", "label", "log", "project", "clim",
     "cmap", "symmetric", "extent", "slug", "contour"])
# defaults cover the trailing fields (label..contour); `slug=None` -> the save
# filename token is derived from the label via _label_slug (back-compatible).
RenderSpec.__new__.__defaults__ = ("field", True, "column", None,
                                   "inferno", False, None, None, None)


# Iso-contours of a second quantity drawn OVER every panel of a RenderSpec
# (attach via RenderSpec(..., contour=ContourSpec(...)) or
# styled(..., contour=...)). `quantity` follows the RenderSpec convention
# (tgdata column index or callable quantity(fn, tgdata, ngas), e.g.
# aux_magnitude('BField')); it is gridded with the SAME backend / projection /
# domain as the underlying map, so the contours line up with the image pixels.
# `levels`: an int n -> n automatic levels SHARED BY ALL PANELS (interior
# values between the pooled 1st/99th percentiles of the contour field; log-
# spaced when log=True), or an explicit sequence of level values. Static
# figures only (the movie path ignores it); contour grids are not blessed.
ContourSpec = namedtuple(
    "ContourSpec",
    ["quantity", "label", "levels", "colors", "log", "linewidths", "alpha"])
ContourSpec.__new__.__defaults__ = ("contour", 8, "w", False, 0.8, 0.9)


# --------------------------------------------------------------------------- #
#  Canonical render style: one {quantity -> (cmap, log, symmetric, slug)} map so
#  the SAME physical quantity is drawn with the same colormap/scale (and gets a
#  stable, collision-free filename token) across every test -- figures become
#  directly comparable, and the per-script label->slug string-matching (which had
#  documented "|B|" collisions) is retired. A test builds a spec with
#  `styled(key, axis1, axis2, quantity, ...)`; per-call cmap/log/symmetric/slug
#  overrides are still allowed for a deliberately different rendering.
RenderStyle = namedtuple("RenderStyle", ["cmap", "log", "symmetric", "slug"])

RENDER_STYLE = {
    # scalars
    "rho":      RenderStyle("inferno", True,  False, "rho"),
    "P":        RenderStyle("viridis", False, False, "P"),       # thermal pressure
    "Pmag":     RenderStyle("magma",   True,  False, "Pmag"),    # 1/2 |B|^2
    "Ekin":     RenderStyle("viridis", False, False, "Ekin"),    # 1/2 rho v^2
    "entropy":  RenderStyle("magma",   False, False, "entropy"),
    "speed":    RenderStyle("viridis", False, False, "speed"),   # |v|
    "T":        RenderStyle("inferno", False, False, "T"),       # temperature
    "tracer":   RenderStyle("viridis", False, False, "tracer"),
    "u":        RenderStyle("viridis", True,  False, "u"),       # internal energy
    # magnetic field
    "Bmag":     RenderStyle("inferno", True,  False, "B"),       # |B|
    "Bx":       RenderStyle("RdBu_r",  False, True,  "Bx"),
    "By":       RenderStyle("RdBu_r",  False, True,  "By"),
    "Bz":       RenderStyle("RdBu_r",  False, True,  "Bz"),
    "Bpol":     RenderStyle("inferno", True,  False, "Bpol"),    # |B_pol| (magnitude)
    "Btor":     RenderStyle("RdBu_r",  False, True,  "Btor"),    # signed toroidal
    "divB":     RenderStyle("RdBu_r",  False, True,  "divB"),    # SIGNED div B
    "divBerr":  RenderStyle("inferno", True,  False, "divBerr"),  # |div B| h/|B| (log)
    # velocity components (signed)
    "vx":       RenderStyle("RdBu_r",  False, True,  "vx"),
    "vy":       RenderStyle("RdBu_r",  False, True,  "vy"),
    "vz":       RenderStyle("RdBu_r",  False, True,  "vz"),
    "vr":       RenderStyle("RdBu_r",  False, True,  "vr"),      # radial velocity
}


def styled(key, axis1, axis2, quantity, *, label, project="slice", clim=None,
           extent=None, cmap=None, log=None, symmetric=None, slug=None,
           contour=None):
    """Build a RenderSpec whose cmap/log/symmetric/slug come from the canonical
    RENDER_STYLE[`key`] (so the quantity is drawn consistently across tests),
    with the per-test `quantity`/`label`/`project`/`clim`/`extent`/`contour`.
    Pass cmap/log/symmetric/slug to deliberately override the canonical style."""
    st = RENDER_STYLE[key]
    return RenderSpec(axis1, axis2, quantity, label,
                      st.log if log is None else log,
                      project, clim,
                      st.cmap if cmap is None else cmap,
                      st.symmetric if symmetric is None else symmetric,
                      extent,
                      slug or st.slug,
                      contour)

# How to draw STREAMLINES of a 2-D vector field for a test (e.g. velocity or
# magnetic field), drawn by `streamline_panels`. (axis1, axis2) is the image
# plane (tgdata position columns, like RenderSpec); `qx`/`qy` are the in-plane
# vector components, each either a tgdata COLUMN INDEX (e.g. 4, 5 for vx, vy) or
# a CALLABLE quantity(fn, tgdata, ngas) -> per-particle array (e.g.
# aux_component('BField', 0)). Lines are colored by the in-plane magnitude with
# `cmap`; `density` is matplotlib streamplot's line density.
StreamSpec = namedtuple(
    "StreamSpec",
    ["axis1", "axis2", "qx", "qy", "label", "cmap", "density", "extent"])
StreamSpec.__new__.__defaults__ = ("flow", "viridis", 1.4, None)


# --------------------------------------------------------------------------- #
#  Aux-field render quantities (generic across tests)
# --------------------------------------------------------------------------- #
# A tipsy snapshot can carry extra per-particle fields in sibling aux files named
# "<snapshot>.<Label>", e.g. "kh64_..._GASOLINE.00060.DivB" / ".HeI" / ".BField".
# These factories turn any such aux field into a RenderSpec `quantity` callable
# (and the same callables work for analyze()), so every test can render/analyze
# arbitrary aux output -- not just the magnetic field. The aux is truncated to
# the gas particles (ngas) so it aligns row-for-row with tgdata. A missing aux
# yields zeros (a run with no magnetic field renders a uniform zero panel rather
# than NaN dots / a blank); pass missing=np.nan for the old blank behavior.

def aux_scalar(label, missing=0.0):
    """quantity callable for a per-particle SCALAR aux field, e.g. 'DivB', 'HeI'."""
    def q(fn, tgdata, ngas):
        a = try_aux(fn, label, ngas)
        return np.full(len(tgdata), missing) if a is None \
            else np.asarray(a, dtype=float)
    return q


def aux_component(label, i, missing=0.0):
    """quantity callable for component i (0=x,1=y,2=z) of a per-particle VECTOR
    aux, e.g. aux_component('BField', 0) -> Bx."""
    def q(fn, tgdata, ngas):
        v = read_vector_aux(fn, label, ngas)
        return np.full(len(tgdata), missing) if v is None \
            else np.asarray(v[:, i], dtype=float)
    return q


def aux_magnitude(label, missing=0.0):
    """quantity callable for the magnitude of a per-particle VECTOR aux,
    e.g. aux_magnitude('BField') -> |B|."""
    def q(fn, tgdata, ngas):
        v = read_vector_aux(fn, label, ngas)
        return np.full(len(tgdata), missing) if v is None \
            else np.sqrt(np.sum(np.asarray(v, dtype=float) ** 2, axis=1))
    return q


def aux_divB_error(div_label="DivB", b_label="BField", h_label="smoothlength",
                   nsmooth=64, missing=0.0):
    """quantity callable for the normalised divergence error |div B| * h / |B|.

    The standard divergence-cleanliness diagnostic for MHD runs: dimensionless,
    so it is comparable across tests. Combines the scalar `div_label` aux
    (e.g. 'DivB'), the magnitude of the vector `b_label` aux (|B|), and the SPH
    smoothing length h: the actual run h from the `h_label` aux ('smoothlength')
    when present, else reconstructed from mass/rho with getsmooth2 (nsmooth, the
    framework's standard h). Where |B|=0 (no field), the error is set to
    `missing` (0) rather than inf, so it renders as zero.
    """
    import IC_smoothlength as smth
    def q(fn, tgdata, ngas):
        d = try_aux(fn, div_label, ngas)
        v = read_vector_aux(fn, b_label, ngas)
        if d is None or v is None:
            return np.full(len(tgdata), missing)
        h = try_aux(fn, h_label, ngas)               # actual SPH h, if written
        if h is None:
            h = smth.getsmooth2(tgdata[:, 0], tgdata[:, 7], nsmooth)
        bmag = np.sqrt(np.sum(np.asarray(v, dtype=float) ** 2, axis=1))
        num = np.abs(np.asarray(d, dtype=float)) * np.asarray(h, dtype=float)
        out = np.full(len(tgdata), float(missing))
        nz = bmag > 0
        out[nz] = num[nz] / bmag[nz]                  # |B|=0 -> missing (no inf)
        return out
    return q


def internal_energy():
    """quantity callable for the specific internal energy u.

    gasoline stores u as a TEMPERATURE in tgdata col 8, NOT as energy; the
    conversion is u = dTuFac * T with dTuFac read from the run '.log' per input
    (read_tufac, which defaults to 1.0 -> u==T with a warning if the log lacks
    dTuFac). Use this anywhere u is rendered/analysed instead of reading col 8
    directly (col 8 is off by the dTuFac factor)."""
    def q(fn, tgdata, ngas):
        tufac = read_tufac(os.path.dirname(fn) or ".")
        return tufac * np.asarray(tgdata[:, 8], dtype=float)
    return q


def thermal_pressure(gamma):
    """quantity callable for the thermal pressure P = (gamma-1) rho u.

    rho is tgdata col 7; u = dTuFac * T (see internal_energy -- col 8 is a
    temperature, so this is NOT (gamma-1) rho * col8). Pass the run's adiabatic
    index gamma."""
    u_of = internal_energy()
    def q(fn, tgdata, ngas):
        rho = np.asarray(tgdata[:, 7], dtype=float)
        return (gamma - 1.0) * rho * u_of(fn, tgdata, ngas)
    return q


def entropy_function(gamma):
    """quantity callable for the entropic function A = P / rho^gamma
    = (gamma-1) u rho^(1-gamma), with u = dTuFac * T (see internal_energy).

    A is conserved along adiabatic flow, so it cleanly traces mixing and shock
    heating (e.g. the kh shear-layer mixing)."""
    u_of = internal_energy()
    def q(fn, tgdata, ngas):
        rho = np.asarray(tgdata[:, 7], dtype=float)
        return (gamma - 1.0) * u_of(fn, tgdata, ngas) * rho ** (1.0 - gamma)
    return q


def magnetic_pressure(b_label="BField", missing=0.0):
    """quantity callable for the magnetic pressure 1/2 |B|^2 from a per-particle
    VECTOR aux (zeros when the aux is absent, i.e. a non-magnetic run)."""
    def q(fn, tgdata, ngas):
        v = read_vector_aux(fn, b_label, ngas)
        return np.full(len(tgdata), missing) if v is None \
            else 0.5 * np.sum(np.asarray(v, dtype=float) ** 2, axis=1)
    return q


def ke_density():
    """quantity callable for the kinetic-energy density 1/2 rho |v|^2."""
    def q(fn, tgdata, ngas):
        v2 = np.sum(np.asarray(tgdata[:, 4:7], dtype=float) ** 2, axis=1)
        return 0.5 * np.asarray(tgdata[:, 7], dtype=float) * v2
    return q


def speed():
    """quantity callable for the speed |v|."""
    def q(fn, tgdata, ngas):
        return np.sqrt(np.sum(np.asarray(tgdata[:, 4:7], dtype=float) ** 2,
                              axis=1))
    return q


# --------------------------------------------------------------------------- #
#  Field rendering at selected time points
# --------------------------------------------------------------------------- #

def _read_time(fn):
    """Cheap-ish read of just a snapshot's time (skips smoothing-length calc)."""
    *_rest, time = tip.readtipsy(fn)
    return float(np.asarray(time).ravel()[0])


def read_box_periods(run_arg):
    """Authoritative domain from the run's gasoline `.log`, or None.

    The `.log` lives in the run folder (`<runname>/<runname>.log`) and records
        # BOX PARAMETERS: ... dxPeriod: <Lx> dyPeriod: <Ly> dzPeriod: <Lz>
    which is the periodic box (centered on the origin) written from the IC setup.
    Returns {1: Lx, 2: Ly, 3: Lz} (keyed by tgdata axis index) or None.

    Reads each period through the shared `read_log_param` (so all `.log` parsing
    goes through one place); the keys are unique in the log.
    """
    lx = read_log_param(run_arg, "dxPeriod")
    ly = read_log_param(run_arg, "dyPeriod")
    lz = read_log_param(run_arg, "dzPeriod")
    if None in (lx, ly, lz):
        return None
    return {1: lx, 2: ly, 3: lz}


# pynbody field names for the common tgdata columns, used by the SPH backend.
PYNBODY_QTY = {0: "mass", 4: "vx", 5: "vy", 6: "vz", 7: "rho", 8: "temp"}

_HAVE_PYNBODY = None  # cached probe result


def _have_pynbody():
    global _HAVE_PYNBODY
    if _HAVE_PYNBODY is None:
        try:
            import pynbody  # noqa: F401
            _HAVE_PYNBODY = True
        except Exception:
            _HAVE_PYNBODY = False
    return _HAVE_PYNBODY


def _field_image(a1, a2, q, rho, z, x0, x1, y0, y1, res,
                 proj="column", zmid=0.0, slice_half=None):
    """Histogram fallback field map, as (res,res), origin lower.

    Line-of-sight reduction along the 3rd axis (`z`), to match pynbody:
      slice     -> mean of q over particles in the mid-plane slab |z-zmid|<slice_half,
      column    -> sum of q per pixel (field integral),
      rhocolumn -> sum(q*rho)/sum(rho) per pixel (density-weighted average),
      average   -> mean of q per pixel (count-weighted average).
    Used when pynbody is unavailable; pynbody's kernel-smoothed image is preferred.
    """
    rng = [[x0, x1], [y0, y1]]
    if proj == "slice" and slice_half is not None:
        keep = np.abs(z - zmid) < slice_half
        a1, a2, q = a1[keep], a2[keep], q[keep]
        rho = rho[keep] if rho is not None else None
        proj = "average"

    def H(w=None):
        out, _, _ = np.histogram2d(a1, a2, bins=res, range=rng, weights=w)
        return out

    with np.errstate(invalid="ignore", divide="ignore"):
        if proj == "column":
            img = H(q)
        elif proj == "rhocolumn" and rho is not None:
            num, den = H(q * rho), H(rho)
            img = np.where(den > 0, num / den, np.nan)
        else:  # average (and slice, reduced above)
            num, den = H(q), H()
            img = np.where(den > 0, num / den, np.nan)
    return img.T


def _label_slug(label):
    """Filesystem-safe token from a render label (e.g. '|B|' -> 'B', 'B_x'->'Bx')."""
    s = re.sub(r"[^0-9A-Za-z]+", "", str(label))
    return s or "field"


def spec_slug(spec):
    """The save-filename token for a RenderSpec: its explicit `slug` (canonical,
    collision-free) when set, else derived from the label (`_label_slug`). Use
    this instead of per-script label string-matching."""
    return spec.slug or _label_slug(spec.label)


def _quantity_values(q, fn, tgdata, ngas):
    """Per-particle field array for a quantity: tgdata column index or callable
    quantity(fn, tgdata, ngas) -> 1-D array (for fields derived from the aux
    file, e.g. magnetic-field components). Always a float array of len(tgdata)."""
    if callable(q):
        return np.asarray(q(fn, tgdata, ngas), dtype=float)
    return np.asarray(tgdata[:, q], dtype=float)


def _field_values(spec, fn, tgdata, ngas):
    """Per-particle field array for a RenderSpec (see _quantity_values)."""
    return _quantity_values(spec.quantity, fn, tgdata, ngas)


def _pynbody_image(fn, spec, rh, xmid, ymid, zmid, res,
                   av_z=False, slice_half=None, q=None):
    """SPH kernel-smoothed field map via pynbody, as (res,res), origin lower.

    Position columns are permuted so the requested (axis1, axis2) plane maps to
    pynbody's x-y, recentred on the box (3rd axis -> z, recentred on zmid). The
    line-of-sight reduction is controlled by `av_z` (False=column integral,
    True=LOS average, or a field name like 'rho'=weighted average). For a slice,
    the snapshot is first cut to the mid-plane slab |z|<slice_half and averaged.

    For a named tgdata column (density, ...) pynbody renders its own array. For a
    derived/aux field (callable quantity, e.g. magnetic field), the precomputed
    per-particle values `q` are injected into the snapshot as a temporary array
    and kernel-smoothed -- mirroring the `s.g['Bx']=...; sph.image(qty='Bx')`
    idiom. Returns None if it cannot render (no values, or length mismatch).
    """
    import warnings
    import pynbody
    from pynbody.plot import sph as psph

    qty = PYNBODY_QTY.get(spec.quantity)
    a1i, a2i = spec.axis1 - 1, spec.axis2 - 1
    third = ({0, 1, 2} - {a1i, a2i}).pop()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # silence the "no param file" notice
        s = pynbody.load(fn)
        g = s.gas
        if qty is None:                            # derived/aux field: inject q
            if q is None or len(q) != len(g):
                return None
            g["_renderq"] = pynbody.array.SimArray(np.asarray(q, dtype=float))
            qty = "_renderq"
        pos = np.asarray(g["pos"], dtype=float)
        newpos = np.empty_like(pos)
        newpos[:, 0] = pos[:, a1i] - xmid
        newpos[:, 1] = pos[:, a2i] - ymid
        newpos[:, 2] = pos[:, third] - zmid
        g["pos"] = newpos
        if slice_half is not None:                 # thin mid-plane slab
            g = g[np.abs(g["z"]) < slice_half]
            av_z = True                            # average the field in the slab
        im = psph.image(g, qty=qty, width=2.0 * rh, resolution=res,
                        noplot=True, av_z=av_z, log=False)
    return np.asarray(im, dtype=float)


def _smoothing_length(fn, mass, rho, ngas, nsmooth=None):
    """Per-particle smoothing length for the render backends, as a gas-length
    (ngas) array. With `nsmooth=None` the run's `smoothlength` aux is used when
    present, else h = getsmooth2(mass, rho, 64). When `nsmooth` is given the aux
    is ignored and h = getsmooth2(mass, rho, nsmooth) (the --ns knob). Shared by
    the grid deposit and the particles-backend slice so both see the same h.
    """
    import IC_smoothlength as smth
    if nsmooth is not None:
        h = smth.getsmooth2(mass, rho, int(nsmooth))
    else:
        h = try_aux(fn, "smoothlength", ngas)
        if h is None:
            h = smth.getsmooth2(mass, rho, 64)
    return np.asarray(h, dtype=float).ravel()[:ngas]


def _grid_image(fn, tgdata, ngas, spec, q, x0, x1, y0, y1, res,
                proj, zmid, periods, nsmooth=None, method="sph"):
    """Field map via the sph_interp 2-D deposit, as (res,res), origin lower.
    The pynbody-free SPH backend (`method='sph'`) or its exact mass-conserving
    Petkova column variant (`method='petkova'`).

    Builds a code-agnostic `Particles` container from the already-loaded gas
    rows (pos/mass/rho + the precomputed per-particle field `q`) and a smoothing
    length, then deposits onto a `Projection2D` plane in the requested
    (axis1, axis2) world plane. Smoothing length `h`: with `nsmooth=None`
    (default) the run's `smoothlength` aux is used when present, else
    h = getsmooth2(mass, rho, 64) (the divBerr render convention). When `nsmooth`
    is given (e.g. --ns 128) the aux is IGNORED and h = getsmooth2(mass, rho,
    nsmooth) for every particle -- a knob to smooth more (or less) than the run.
    The projection mode maps 1:1 to the
    framework's `--render-project`: column (LOS field integral / surface density
    for rho), average (kernel-weighted LOS mean), rhocolumn (mass-weighted mean),
    slice (kernel value on the mid-plane). This is a SPLASH-style deposit that
    scales as res^2 * N (no 3-D cube) and reads each snapshot once.

    `method='petkova'` integrates each particle's projected kernel EXACTLY over
    every pixel (mass-conserving per pixel even when h<pixel), at higher cost; it
    has no 'slice' (not a column integral) -- slice silently falls back to 'sph'.
    Returns None (caller falls back to the histogram) if sph_interp is
    unimportable.
    """
    try:
        from sph_interp import Particles, Projection2D, interpolate
    except Exception as e:
        print(f"render: sph_interp unavailable ({e}); using histogram.",
              file=sys.stderr)
        return None

    g = slice(0, ngas)
    pos = np.asarray(tgdata[g, 1:4], dtype=float)
    mass = np.asarray(tgdata[g, 0], dtype=float)
    rho = np.asarray(tgdata[g, 7], dtype=float)
    qg = np.asarray(q[:ngas], dtype=float)
    h = _smoothing_length(fn, mass, rho, ngas, nsmooth)

    keep = np.isfinite(qg) & np.isfinite(rho) & (rho > 0.0) & (h > 0.0)
    pos, mass, rho, h, qg = pos[keep], mass[keep], rho[keep], h[keep], qg[keep]
    if pos.shape[0] == 0:
        return None

    if periods is not None:                      # centered periodic box
        box = np.array([periods[1], periods[2], periods[3]], dtype=float)
        periodic = True
    else:                                        # unknown period -> non-periodic
        box = pos.max(axis=0) - pos.min(axis=0)
        periodic = False

    # Petkova-2D is a column integral with no 'slice' mode -> fall back to SPH.
    if method == "petkova" and proj == "slice":
        print("render: grid-petkova has no 'slice' (not a column integral); "
              "using SPH for this slice.", file=sys.stderr)
        method = "sph"

    p = Particles(pos=pos, mass=mass, rho=rho, h=h, values={"q": qg},
                  box=box, periodic=periodic)
    target = Projection2D(axis1=spec.axis1 - 1, axis2=spec.axis2 - 1,
                          lo=(x0, y0), hi=(x1, y1), npx=int(res),
                          mode=proj, zslice=zmid)
    out = interpolate(p, target, method=method, values=["q"])
    # data["q"] is (nx,ny) indexed [axis1, axis2]; transpose to (row=y, col=x)
    # origin-lower to match the histogram/pynbody panels.
    return np.asarray(out.data["q"], dtype=float).T


def _crop_to_view(im, im_extent, x0, x1, y0, y1):
    """Sub-array of an (origin-lower) image covering the view box [x0,x1]x[y0,y1].

    Used to set the color scale from only the in-domain pixels (pynbody renders a
    square that can include out-of-domain padding we crop away in the view).
    """
    ex0, ex1, ey0, ey1 = im_extent
    ny, nx = im.shape
    xs = ex0 + (np.arange(nx) + 0.5) * (ex1 - ex0) / nx
    ys = ey0 + (np.arange(ny) + 0.5) * (ey1 - ey0) / ny
    cm = (xs >= x0) & (xs <= x1)
    rm = (ys >= y0) & (ys <= y1)
    if not cm.any() or not rm.any():
        return im
    return im[np.ix_(rm, cm)]


def select_snapshots(files, times=None, n=4, time_of=None):
    """Pick snapshot files to render: nearest to `times`, or `n` evenly spaced.

    `time_of(fn) -> float` maps a snapshot to the axis the request is expressed
    in (default: the raw snapshot time). For blob this is t/t_crush, so
    `--render-times 10` selects the t_crush=10 snapshot. Returns (file, value).
    """
    time_of = time_of or _read_time
    all_t = [time_of(fn) for fn in files]
    if times is not None:
        chosen = []
        for tt in times:
            j = int(np.argmin([abs(a - tt) for a in all_t]))
            chosen.append((files[j], all_t[j]))
        return chosen
    k = min(n, len(files))
    idx = sorted(set(np.linspace(0, len(files) - 1, k).round().astype(int)))
    return [(files[i], all_t[i]) for i in idx]


def missing_render_times(input_arg, times, time_of=None, tol=None):
    """Requested render `times` that have NO matching snapshot in the run.

    `select_snapshots` always returns the NEAREST snapshot, so a request past the
    run's end (e.g. t=2.5 for a run that stops at 2.0) silently snaps to the last
    dump. This reports such times so a bless can refuse to save a reference at the
    wrong time. A time matches if the nearest snapshot is within `tol` (default:
    half the median dump spacing -- an in-range request always matches, a time
    beyond the dumped range does not). Returns the unmatched times (empty if all
    present)."""
    try:
        files = find_files(input_arg)
    except FileNotFoundError:
        return list(times)
    tof = time_of or _read_time
    all_t = sorted(tof(fn) for fn in files)
    if not all_t:
        return list(times)
    if tol is None:
        if len(all_t) > 1:
            d = np.diff(all_t)
            d = d[d > 0]
            dt = float(np.median(d)) if d.size else 0.0
        else:
            dt = 0.0
        tol = 0.5 * dt + 1e-6 * (1.0 + abs(all_t[-1]))
    return [tt for tt in times
            if min(abs(a - tt) for a in all_t) > tol]


def snapshot_times(input_arg, time_of=None):
    """All snapshot times of a run, sorted ascending (empty if none).

    Used to bless EVERY dump as a render reference: the saved grids are a
    fixed-resolution image (set by --render-res, NOT the particle count), so the
    storage is independent of simulation resolution and saving all times stays
    cheap even for very high-N runs."""
    try:
        files = find_files(input_arg)
    except FileNotFoundError:
        return []
    tof = time_of or _read_time
    return sorted(tof(fn) for fn in files)


def _resolve_render_mode(backend, spec):
    """Map the requested backend to one of grid / pynbody / hist / particles.

    The DEFAULT backend is 'grid' (sph_interp 2-D SPH kernel deposit) -- a
    pynbody-free, quantity-agnostic, deterministic SPH-smoothed map; 'grid-petkova'
    is its exact-column variant. 'auto' is the legacy pynbody-preferring path:
    pynbody kernel-smooths both named tgdata columns (density, ...) and derived/
    aux fields (callable quantity, e.g. magnetic field -- the computed values are
    injected into the snapshot and smoothed), so under 'auto' pynbody is used
    whenever importable; when it is not, a derived field falls back to the
    per-particle 'particles' backend rather than the coarse histogram (which looks
    bad for a field map), and named columns fall back to the histogram. An explicit
    --render-backend always wins.
    """
    if backend == "particles":
        return "particles"
    if backend in ("grid", "grid-petkova"):
        return backend
    if backend == "hist":
        return "hist"
    if backend == "pynbody":
        if _have_pynbody():
            return "pynbody"
        print("render: pynbody requested but not importable; using fallback.",
              file=sys.stderr)
    elif _have_pynbody():
        return "pynbody"
    # pynbody unavailable: particles for derived fields, histogram for columns.
    return "particles" if callable(spec.quantity) else "hist"


def _render_one_input(input_arg, label, spec, times, n, res, mode, extent,
                      project, slice_frac, time_of, nsmooth=None):
    """Build the rendered panels (one per chosen snapshot) for a single run.

    Returns dict(label, panels, times, box, pt) or None if nothing to render.
    `panels` items are ("image", array, extent) or ("scatter", a1, a2, q).
    """
    try:
        files = find_files(input_arg)
    except FileNotFoundError as e:
        print(f"render: input '{input_arg}': {e}", file=sys.stderr)
        return None
    chosen = select_snapshots(files, times=times, n=n, time_of=time_of)
    if not chosen:
        print(f"render: no snapshots for '{input_arg}'", file=sys.stderr)
        return None

    # Domain box: 1) explicit extent (--render-extent, else the spec's own
    # default), 2) .log BOX PARAMETERS (centered periodic box), 3) inferred
    # per-axis from robust 0.5-99.5 percentiles.
    tg0, *_ = _cached_readtipsy(chosen[0][0])
    if extent is None and getattr(spec, "extent", None) is not None:
        extent = list(spec.extent)
    periods = None if extent is not None else read_box_periods(input_arg)
    if extent is not None:
        x0, x1, y0, y1 = extent
        domain_src = "render extent"
    elif periods is not None:
        x0, x1 = -0.5 * periods[spec.axis1], 0.5 * periods[spec.axis1]
        y0, y1 = -0.5 * periods[spec.axis2], 0.5 * periods[spec.axis2]
        domain_src = "log BOX PARAMETERS"
    else:
        def _pad(lo, hi):
            mid, half = 0.5 * (lo + hi), 0.51 * (hi - lo)
            return mid - half, mid + half
        x0, x1 = _pad(*np.nanpercentile(tg0[:, spec.axis1], [0.5, 99.5]))
        y0, y1 = _pad(*np.nanpercentile(tg0[:, spec.axis2], [0.5, 99.5]))
        domain_src = "inferred (percentile)"
    print(f"render[{label}]: domain {domain_src}: "
          f"{'xyz'[spec.axis1-1]}=[{x0:.4g},{x1:.4g}] "
          f"{'xyz'[spec.axis2-1]}=[{y0:.4g},{y1:.4g}]")
    xmid, ymid = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
    rh = max(0.5 * (x1 - x0), 0.5 * (y1 - y0))   # pynbody square width; we crop

    proj = project or spec.project
    third = ({1, 2, 3} - {spec.axis1, spec.axis2}).pop()
    if periods is not None:
        zmid, zext = 0.0, periods[third]
    else:
        zlo, zhi = np.nanpercentile(tg0[:, third], [0.5, 99.5])
        zmid, zext = 0.5 * (zlo + zhi), (zhi - zlo)
    # Slab half-thickness for the SLAB backends (pynbody/hist, and the particles
    # backend when --render-slice-frac is given). None -> each backend's default:
    # 0.1 for pynbody/hist; the particles branch instead uses per-particle h.
    if proj != "slice":
        slice_half = None
    elif slice_frac is not None:
        slice_half = slice_frac * zext
    else:
        slice_half = 0.1 * zext
    av_z = {"column": False, "average": True, "rhocolumn": "rho"}.get(proj, False)
    if mode in ("pynbody", "grid", "grid-petkova"):
        res_eff = int(min(res, 1024))            # kernel-smoothed: full res is fine
    else:
        res_eff = int(min(res, max(32, round(np.sqrt(len(tg0))))))
    # Histogram deposits get noisy past ~sqrt(N) pixels (used for the contour
    # fallback even when the main map is kernel-smoothed at full res).
    res_hist = int(min(res, max(32, round(np.sqrt(len(tg0))))))

    contour = getattr(spec, "contour", None)

    def _contour_grid(fn, tgdata, ngas):
        """(grid, extent) of the contour quantity for one snapshot, via the SAME
        backend/projection/domain as the main map (histogram fallback)."""
        cq = _quantity_values(contour.quantity, fn, tgdata, ngas)
        cspec = spec._replace(quantity=contour.quantity)
        if mode == "pynbody":
            try:
                cim = _pynbody_image(fn, cspec, rh, xmid, ymid, zmid, res_eff,
                                     av_z=av_z, slice_half=slice_half, q=cq)
                if cim is not None:
                    return cim, (xmid - rh, xmid + rh, ymid - rh, ymid + rh)
            except Exception as e:
                print(f"render: pynbody contour failed on {fn} ({e}); using "
                      f"histogram.", file=sys.stderr)
        elif mode in ("grid", "grid-petkova"):
            try:
                cim = _grid_image(fn, tgdata, ngas, cspec, cq, x0, x1, y0, y1,
                                  res_eff, proj, zmid, periods, nsmooth,
                                  method=("petkova" if mode == "grid-petkova"
                                          else "sph"))
                if cim is not None:
                    return cim, (x0, x1, y0, y1)
            except Exception as e:
                print(f"render: {mode} contour deposit failed on {fn} ({e}); "
                      f"using histogram.", file=sys.stderr)
        cim = _field_image(tgdata[:, spec.axis1], tgdata[:, spec.axis2],
                           cq, tgdata[:, 7], tgdata[:, third],
                           x0, x1, y0, y1, res_hist,
                           proj=proj, zmid=zmid, slice_half=slice_half)
        return cim, (x0, x1, y0, y1)

    panels, contours, ts = [], [], []
    for fn, t in chosen:
        tgdata, _td, _ts, hdr, _tm = _cached_readtipsy(fn)
        ngas = int(hdr[2])
        q = _field_values(spec, fn, tgdata, ngas)
        contours.append(_contour_grid(fn, tgdata, ngas)
                        if contour is not None else None)
        if mode == "particles":
            a1 = np.asarray(tgdata[:, spec.axis1], dtype=float)
            a2 = np.asarray(tgdata[:, spec.axis2], dtype=float)
            keep = (a1 >= x0) & (a1 <= x1) & (a2 >= y0) & (a2 <= y1)
            keep &= np.isfinite(q)     # don't scatter NaN/inf as stray dots
            if proj == "slice":
                zrel = np.abs(np.asarray(tgdata[:, third], dtype=float) - zmid)
                if slice_frac is None:
                    # default: keep each gas particle whose smoothing kernel
                    # reaches the plane within half a smoothing length
                    # (|z-zmid| < 0.5*h_i) -- self-adapting to the local
                    # resolution instead of a fixed slab.
                    h = _smoothing_length(fn, tgdata[:ngas, 0], tgdata[:ngas, 7],
                                          ngas, nsmooth)
                    ksl = np.zeros(len(tgdata), dtype=bool)
                    ksl[:ngas] = zrel[:ngas] < 0.5 * h
                    keep &= ksl
                else:                  # explicit --render-slice-frac: fixed slab
                    keep &= zrel < slice_half
            panels.append(("scatter", a1[keep], a2[keep], q[keep]))
        else:
            im = None
            if mode == "pynbody":
                try:
                    im = _pynbody_image(fn, spec, rh, xmid, ymid, zmid, res_eff,
                                        av_z=av_z, slice_half=slice_half, q=q)
                except Exception as e:
                    print(f"render: pynbody failed on {fn} ({e}); using histogram.",
                          file=sys.stderr)
                    im = None
            elif mode in ("grid", "grid-petkova"):
                gmethod = "petkova" if mode == "grid-petkova" else "sph"
                try:
                    im = _grid_image(fn, tgdata, ngas, spec, q, x0, x1, y0, y1,
                                     res_eff, proj, zmid, periods, nsmooth,
                                     method=gmethod)
                except Exception as e:
                    print(f"render: {mode} deposit failed on {fn} ({e}); using "
                          f"histogram.", file=sys.stderr)
                    im = None
                if im is not None:               # grid fills the view extent exactly
                    panels.append(("image", im, (x0, x1, y0, y1)))
                    ts.append(t)
                    continue
            if im is None:
                im = _field_image(tgdata[:, spec.axis1], tgdata[:, spec.axis2],
                                  q, tgdata[:, 7], tgdata[:, third],
                                  x0, x1, y0, y1, res_eff,
                                  proj=proj, zmid=zmid, slice_half=slice_half)
                panels.append(("image", im, (x0, x1, y0, y1)))
            else:
                panels.append(("image", im, (xmid - rh, xmid + rh,
                                             ymid - rh, ymid + rh)))
        ts.append(t)

    pt = float(np.clip(1500.0 / np.sqrt(len(tg0) + 1.0), 0.4, 5.0))
    return dict(label=label, panels=panels, contours=contours, times=ts,
                box=(x0, x1, y0, y1), pt=pt)


def _tile_axes(axes, same_x, same_y):
    """Drop duplicated inner tick labels so shared-limit panels read as a multiplot.

    Robust to ragged grids: keeps y on the left-most populated cell of each row
    and x on the bottom-most populated cell of each column.
    """
    nr, nc = axes.shape
    if same_y:
        for r in range(nr):
            first = min((c for c in range(nc) if axes[r][c].has_data()), default=0)
            for c in range(nc):
                if c != first:
                    axes[r][c].tick_params(labelleft=False)
                    axes[r][c].set_ylabel("")
    if same_x:
        for c in range(nc):
            last = max((r for r in range(nr) if axes[r][c].has_data()),
                       default=nr - 1)
            for r in range(nr):
                if r != last:
                    axes[r][c].tick_params(labelbottom=False)
                    axes[r][c].set_xlabel("")


def save_render_grids(path, panels, times, spec, project):
    """Persist a render's IMAGE-panel grids as a baseline npz.

    Stores the kernel-smoothed field per snapshot (pre-colormap, as FLOAT32 --
    half the size, ample precision for a visual reference + difference; numeric
    regression uses the JSON metric, not these grids) + its extent + time, plus
    the spec metadata, so a later compare can draw the reference with the CURRENT
    run's colormap/limits and difference the two. The grid resolution is fixed by
    --render-res, independent of the particle count, so file size does not grow
    with simulation resolution. Scatter panels (the `particles` backend) can't be
    blessed and are skipped; if none remain nothing is written."""
    ims, exts, ts = [], [], []
    for p, t in zip(panels, times):
        if p[0] != "image":
            continue
        ims.append(np.asarray(p[1], dtype=np.float32))
        exts.append(np.asarray(p[2], dtype=float))
        ts.append(float(t))
    if not ims:
        print(f"render-reference: no image panels to bless for {path} "
              f"(scatter backend?); skipped.", file=sys.stderr)
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    payload = {
        "n": np.array(len(ims)),
        "times": np.array(ts, dtype=float),
        "label": np.array(spec.label),
        "log": np.array(bool(spec.log)),
        "symmetric": np.array(bool(spec.symmetric)),
        "cmap": np.array(spec.cmap),
        "project": np.array(project),
        "axis1": np.array(int(spec.axis1)),
        "axis2": np.array(int(spec.axis2)),
    }
    for i, (im, ext) in enumerate(zip(ims, exts)):
        payload[f"im_{i}"] = im
        payload[f"ext_{i}"] = ext
    np.savez_compressed(path, **payload)
    print(f"saved render reference -> {path}  ({len(ims)} panels)")


def load_render_grids(path):
    """Load a render-reference npz written by `save_render_grids`, or None.

    Returns dict(grids=[(image, extent_tuple, time), ...], meta={...})."""
    if not path or not os.path.isfile(path):
        return None
    d = np.load(path, allow_pickle=False)
    n = int(d["n"])
    grids = [(d[f"im_{i}"], tuple(float(v) for v in d[f"ext_{i}"]),
              float(d["times"][i])) for i in range(n)]
    meta = dict(label=str(d["label"]), log=bool(d["log"]),
                symmetric=bool(d["symmetric"]), cmap=str(d["cmap"]),
                project=str(d["project"]), axis1=int(d["axis1"]),
                axis2=int(d["axis2"]))
    return dict(grids=grids, meta=meta)


def render_panels(entries, spec, times=None, n=4, res=512, save=None,
                  backend="auto", extent=None, project=None, slice_frac=None,
                  time_to_x=None, time_label="t", params=None, render_rows=None,
                  aspect="equal", nsmooth=None,
                  save_render_reference=None, render_reference=None,
                  residual=False):
    """Render `spec.quantity` field maps for one or more runs in a single figure.

    entries : list of (label, input_arg). A single run -> one row of panels (one
    per snapshot). Several runs -> a grid with DIRECTORIES as columns and TIMES
    as rows (so 6 directories give 6 columns). All panels share one color scale.

    save_render_reference : if a path is given, write the FIRST run's image-panel
    grids there as a baseline npz (`save_render_grids`) and return without drawing
    -- the render analogue of `--save-reference` for the metric series.
    residual : opt-in (--residual) `times x [A | B | A-B]` comparison grid: A and
    B share the field color scale, the difference gets its own diverging scale.
    With TWO entries B is the second run, differenced directly (saved
    `<save>_residual.png`); with ONE entry B is the blessed baseline at
    `render_reference` (saved `<save>_render_ref.png`). More than two entries is
    an error. Without `residual`, a passed `render_reference` is ignored (with a
    hint).

    See module docstring for backend/projection/extent/time semantics.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    mode = _resolve_render_mode(backend, spec)
    time_of = ((lambda fn: time_to_x(_read_time(fn), params)) if time_to_x
               else _read_time)
    proj = project or spec.project
    if proj != "slice":
        print(f"render: projection '{proj}'")
    elif mode in ("grid", "grid-petkova"):
        print("render: projection 'slice' (grid: mid-plane at z=zmid, 3-D kernel)")
    elif mode == "particles" and slice_frac is None:
        print("render: projection 'slice' (particles: |z-zmid| < 0.5 h per particle)")
    else:
        print(f"render: projection 'slice' "
              f"(slab frac {slice_frac if slice_frac is not None else 0.1} "
              f"of 3rd axis)")

    rows = []
    for lab, inp in entries:
        r = _render_one_input(inp, lab, spec, times, n, res, mode, extent,
                              project, slice_frac, time_of, nsmooth=nsmooth)
        if r is not None:
            rows.append(r)
    if not rows:
        print("render: nothing to render", file=sys.stderr)
        return
    ntimes = max(len(r["panels"]) for r in rows)

    # --- Bless: write the first run's image grids as the render baseline, stop.
    if save_render_reference is not None:
        save_render_grids(save_render_reference, rows[0]["panels"],
                          rows[0]["times"], spec, proj)
        return

    # --- Residual comparison (--residual): build the "B" side of the
    # [A | B | A-B] grid -- either the SECOND run's panels (two inputs) or the
    # blessed render baseline (one input + --reference). Opt-in: without
    # `residual` a passed render_reference is ignored, with a hint.
    ref, ref_name = None, None
    if residual:
        if len(rows) > 2:
            raise SystemExit(f"render --residual: at most two inputs (one run "
                             f"vs the reference, or two runs); got {len(rows)}")
        if len(rows) == 2:
            if any(p[0] != "image" for r in rows for p in r["panels"]):
                raise SystemExit("render --residual: the 'particles' backend "
                                 "has no grid to difference; use a gridded "
                                 "backend (grid/pynbody/hist)")
            if render_reference is not None:
                print("render --residual: two inputs given; differencing them "
                      "and ignoring the blessed render reference.",
                      file=sys.stderr)
            ref = dict(grids=[(p[1], p[2], t) for p, t in
                              zip(rows[1]["panels"], rows[1]["times"])])
            ref_name = rows[1]["label"]
            # rows[1] re-enters the color pool below via the ref grids.
            rows = rows[:1]
        elif render_reference is not None:
            ref = load_render_grids(render_reference)
            ref_name = "reference"
            if ref is not None:
                print(f"render reference: {render_reference} "
                      f"({len(ref['grids'])} panels)")
            else:
                print(f"render --residual: no blessed render grids for "
                      f"'{spec.label}'; rendering the simulation only.",
                      file=sys.stderr)
        else:
            print(f"render --residual: no reference resolved for "
                  f"'{spec.label}'; rendering the simulation only.",
                  file=sys.stderr)
    elif render_reference is not None and len(rows) == 1:
        print("render: blessed render grids exist for this quantity; pass "
              "--residual to draw the sim | reference | difference figure.",
              file=sys.stderr)

    # Shared color scale from in-domain rendered values (adapts to projection).
    pool = []
    for r in rows:
        bx = r["box"]
        for p in r["panels"]:
            v = p[3] if p[0] == "scatter" else _crop_to_view(p[1], p[2], *bx)
            v = np.asarray(v, dtype=float).ravel()
            v = v[np.isfinite(v)]
            if spec.log:
                v = v[v > 0]
            if v.size:
                pool.append(v)
    # Include the reference grids so the sim and reference panels share one color
    # scale (essential for a like-for-like visual comparison).
    if ref is not None:
        for im, _ext, _t in ref["grids"]:
            v = np.asarray(im, dtype=float).ravel()
            v = v[np.isfinite(v)]
            if spec.log:
                v = v[v > 0]
            if v.size:
                pool.append(v)
    # Degenerate field (e.g. no BField aux -> all zeros, or no positive values on
    # a log scale): still draw a UNIFORM ZERO panel (linear, bottom of the scale)
    # rather than dots/blank, so the figure shows the field is absent/zero.
    use_log = spec.log
    if not pool:
        cl = spec.clim(params, proj) if spec.clim is not None else None
        vmin, vmax = cl if cl is not None else (0.0, 1.0)
        use_log = False
        print(f"render: '{spec.label}' has no {'positive ' if spec.log else ''}"
              f"values (missing aux / zero field); plotting as zero "
              f"[{vmin:.4g}, {vmax:.4g}].", file=sys.stderr)
    else:
        vals = np.concatenate(pool)
        # Color limits: this quantity's natural range if the RenderSpec declares
        # one (spec.clim(params, proj), e.g. blob density medium..cloud),
        # otherwise robust 1st/99th percentiles of the rendered values.
        vmin, vmax = np.percentile(vals, 1), np.percentile(vals, 99)
        if spec.clim is not None:
            cl = spec.clim(params, proj)
            if cl is not None:
                vmin, vmax = cl
        if spec.symmetric:             # zero-centered range for signed fields
            a = max(abs(vmin), abs(vmax))
            vmin, vmax = -a, a
        if vmax <= vmin:               # uniform field (e.g. all-zero component)
            eps = max(abs(vmin), 1.0)
            vmin, vmax = vmin - eps, vmax + eps   # keeps 0 centered if symmetric
        print(f"render: color limits [{vmin:.4g}, {vmax:.4g}]"
              + (" (log)" if spec.log else "")
              + (" (symmetric)" if spec.symmetric else ""))
    norm = LogNorm(vmin=max(vmin, 1e-30), vmax=vmax) if use_log else None
    lo = None if use_log else vmin
    hi = None if use_log else vmax

    # Contour overlay: one set of levels SHARED by every panel (like the color
    # scale), so the same line means the same field value across times/runs.
    # An explicit sequence in the ContourSpec is used as-is; an int n gives n
    # interior levels between the pooled 1st/99th percentiles (log-spaced when
    # contour.log). A degenerate contour field (all-NaN/zero) skips the overlay.
    contour = getattr(spec, "contour", None)
    clevels = None
    if contour is not None:
        if not isinstance(contour.levels, int):
            clevels = np.asarray(contour.levels, dtype=float)
        else:
            cpool = []
            for r in rows:
                for cp in r.get("contours", []):
                    if cp is None:
                        continue
                    v = np.asarray(cp[0], dtype=float).ravel()
                    v = v[np.isfinite(v)]
                    if contour.log:
                        v = v[v > 0]
                    if v.size:
                        cpool.append(v)
            if cpool:
                cv = np.concatenate(cpool)
                clo, chi = np.percentile(cv, 1), np.percentile(cv, 99)
                if chi > clo:
                    sp_fn = np.geomspace if contour.log else np.linspace
                    clevels = sp_fn(clo, chi, contour.levels + 2)[1:-1]
        if clevels is None:
            print(f"render: contour field '{contour.label}' is degenerate "
                  f"(missing aux / uniform); skipping the contour overlay.",
                  file=sys.stderr)
        else:
            print(f"render: contour levels ({contour.label}): "
                  + ", ".join(f"{v:.4g}" for v in clevels))

    # Empty pixels (no particles -> NaN, or zero on a log scale) are masked by
    # imshow and would otherwise show the figure background as white salt-and-
    # pepper speckle (e.g. the vacuum corners of the Noh slice once the gas has
    # collapsed inward). Render them at the BOTTOM of the color scale instead, so
    # vacuum reads as "zero field" rather than rendering artifacts.
    cmap = plt.get_cmap(spec.cmap).copy()
    cmap.set_bad(cmap(0.0))

    def draw(ax, panel, box, pt, cpanel=None):
        if panel[0] == "scatter":
            _, a1, a2, q = panel
            m = ax.scatter(a1, a2, c=q, cmap=cmap, norm=norm, vmin=lo,
                           vmax=hi, s=pt, marker="o", linewidths=0, rasterized=True)
        else:
            _, im, ext = panel
            m = ax.imshow(im, extent=ext, origin="lower", cmap=cmap,
                          norm=norm, vmin=lo, vmax=hi, aspect=aspect)
        if cpanel is not None and clevels is not None:
            cim, cext = cpanel
            ny_, nx_ = cim.shape
            # pixel CENTRES of the (origin-lower, edge-anchored) contour grid
            xs = cext[0] + (np.arange(nx_) + 0.5) * (cext[1] - cext[0]) / nx_
            ys = cext[2] + (np.arange(ny_) + 0.5) * (cext[3] - cext[2]) / ny_
            ax.contour(xs, ys, np.ma.masked_invalid(cim), levels=clevels,
                       colors=contour.colors, linewidths=contour.linewidths,
                       alpha=contour.alpha)
        ax.set_xlim(box[0], box[1])
        ax.set_ylim(box[2], box[3])
        # 'equal' (default) keeps circles circular; the imshow aspect already
        # does this for image panels, but set it explicitly so scatter panels
        # match too.
        ax.set_aspect(aspect)
        return m

    setup_rcparams()
    a1name, a2name = "xyz"[spec.axis1 - 1], "xyz"[spec.axis2 - 1]
    backend_name = {"pynbody": "pynbody SPH", "grid": "sph_interp grid",
                    "grid-petkova": "sph_interp grid (Petkova exact)",
                    "hist": "histogram", "particles": "particles"}[mode]
    csuffix = (f", contours: {contour.label}"
               if contour is not None and clevels is not None else "")
    m = None

    # --- Reference comparison layout: rows = times, columns = sim | reference |
    # (sim - reference). A single input only (guaranteed: ref is None unless
    # len(rows)==1). Sim and reference share the field color scale/cmap; the
    # difference uses its own zero-centered diverging scale.
    if ref is not None:
        r = rows[0]
        npan = len(r["panels"])           # one row per RENDERED SIM time
        # Match each sim panel to the reference grid at the SAME time (not by
        # index): rendering only t=0.5 must compare against the reference's t=0.5
        # panel, not its first (t=0.25). A match is the nearest reference time
        # within half the closest reference-time spacing (so a sim time with no
        # blessed reference -- e.g. t=1.75 -- is left unmatched, sim-only).
        ref_times = np.array([g[2] for g in ref["grids"]], dtype=float)
        if ref_times.size > 1:
            sp = np.diff(np.sort(ref_times))
            sp = sp[sp > 0]
            ttol = 0.5 * float(sp.min()) + 1e-6 if sp.size else 1e-6
        else:
            ttol = max(1e-6, 1e-3 * (1.0 + abs(ref_times[0])))
        matched = []                      # ref-grid index (or None) per sim panel
        for st in r["times"]:
            j = int(np.argmin(np.abs(ref_times - st)))
            matched.append(j if abs(ref_times[j] - st) <= ttol else None)
        unm = [f"{st:g}" for st, m in zip(r["times"], matched) if m is None]
        if unm:
            print(f"render residual: no {ref_name} panel at time(s) "
                  f"{', '.join(unm)} ({ref_name} has t="
                  f"{', '.join(f'{t:g}' for t in ref_times)}); those rows show "
                  f"the first run only.", file=sys.stderr)
        # Per-time difference and a shared symmetric scale for the diff column.
        diffs, dvals = [], []
        for i in range(npan):
            sp_panel, j = r["panels"][i], matched[i]
            d = None
            if j is not None and sp_panel[0] == "image":
                rim, rext, _rt = ref["grids"][j]
                sim, rr = np.asarray(sp_panel[1], float), np.asarray(rim, float)
                if sim.shape == rr.shape and np.allclose(sp_panel[2], rext,
                                                         rtol=1e-3, atol=1e-9):
                    d = sim - rr
                else:
                    print(f"render reference: panel {i} grid mismatch "
                          f"(sim {sim.shape} vs ref {rr.shape}); no difference.",
                          file=sys.stderr)
            diffs.append(d)
            if d is not None:
                dv = d[np.isfinite(d)]
                if dv.size:
                    dvals.append(np.abs(dv))
        damax = float(np.concatenate(dvals).max()) if dvals else 1.0
        if damax <= 0:
            damax = 1.0
        diffcmap = plt.get_cmap("RdBu_r").copy()
        diffcmap.set_bad(diffcmap(0.5))

        # Layout: [sim | ref | field-cbar | <spacer> | diff | diff-cbar] via
        # gridspec with two dedicated thin colorbar columns, so the colorbars sit
        # at their own column edges instead of stealing space between the panels.
        # The field colorbar is in the MIDDLE, so its right-side label/ticks need a
        # blank spacer column before the diff panel (else the rho label overlaps
        # the diff figure; the diff colorbar at the far edge needs none). Panels
        # tile seamlessly: inner tick labels dropped, x only bottom row, y only left.
        fig = plt.figure(figsize=(13, 4.0 * npan))
        gs = fig.add_gridspec(npan, 6, width_ratios=[1, 1, 0.06, 0.16, 1, 0.06],
                              wspace=0.07, hspace=0.07)
        axes = np.empty((npan, 3), dtype=object)
        for i in range(npan):
            axes[i][0] = fig.add_subplot(gs[i, 0])
            axes[i][1] = fig.add_subplot(gs[i, 1])
            axes[i][2] = fig.add_subplot(gs[i, 4])
        cax_field = fig.add_subplot(gs[:, 2])
        cax_diff = fig.add_subplot(gs[:, 5])
        col_titles = (("simulation", "reference", "sim − reference")
                      if ref_name == "reference" else
                      (r["label"], ref_name, f"{r['label']} − {ref_name}"))
        mfield, mdiff = None, None
        for i in range(npan):
            box = r["box"]
            j = matched[i]
            mfield = draw(axes[i][0], r["panels"][i], box, r["pt"],
                          cpanel=r["contours"][i])
            if j is not None:
                draw(axes[i][1], ("image",) + ref["grids"][j][:2], box, r["pt"])
            else:
                axes[i][1].text(0.5, 0.5, "no reference\nat this time",
                                ha="center", va="center", color="0.5",
                                fontsize=9, transform=axes[i][1].transAxes)
                axes[i][1].set_xlim(box[0], box[1])
                axes[i][1].set_ylim(box[2], box[3])
                axes[i][1].set_aspect(aspect)
            axd = axes[i][2]
            if diffs[i] is not None:
                mdiff = axd.imshow(diffs[i], extent=ref["grids"][j][1],
                                   origin="lower", cmap=diffcmap,
                                   vmin=-damax, vmax=damax, aspect=aspect)
                axd.set_xlim(box[0], box[1])
                axd.set_ylim(box[2], box[3])
                axd.set_aspect(aspect)
            for c in range(3):
                ax = axes[i][c]
                if i == 0:
                    ax.set_title(col_titles[c], fontsize=11)
                if c == 0:
                    ax.set_ylabel(f"{time_label} = {r['times'][i]:.3g}")
                else:
                    ax.tick_params(labelleft=False)
                if i == npan - 1:
                    ax.set_xlabel(a1name)
                else:
                    ax.tick_params(labelbottom=False)
        if mfield is not None:
            fig.colorbar(mfield, cax=cax_field, label=spec.label)
        else:
            cax_field.axis("off")
        if mdiff is not None:
            fig.colorbar(mdiff, cax=cax_diff, label=f"Δ {spec.label}")
        else:
            cax_diff.axis("off")
        fig.suptitle(f"{r['label']}  vs {ref_name}  —  {spec.label}  "
                     f"({backend_name}, {proj})")
        if save:
            out = (f"{save}_render_ref.png" if ref_name == "reference"
                   else f"{save}_residual.png")
            fig.savefig(out, bbox_inches="tight")
            print(f"saved render -> {out}")
        else:
            plt.show()
        return

    # When all cells share an x (or y) range, tile them seamlessly: tighten the
    # spacing on that axis and drop the duplicated inner tick labels.
    boxes = [r["box"] for r in rows]
    same_x = all(np.isclose(b[0], boxes[0][0]) and np.isclose(b[1], boxes[0][1])
                 for b in boxes)
    same_y = all(np.isclose(b[2], boxes[0][2]) and np.isclose(b[3], boxes[0][3])
                 for b in boxes)
    gkw = dict(wspace=0.04 if same_y else 0.25,
               hspace=0.06 if same_x else 0.30)

    # Per-panel figure size: with the (default) equal aspect, a non-square
    # domain drawn into a square panel collapses to a thin strip (e.g. the tall
    # rt slab) -- size each panel by the domain aspect ratio instead (clamped to
    # 3:1 so degenerate boxes stay usable).
    x0_, x1_, y0_, y1_ = boxes[0]
    if aspect == "equal" and (x1_ > x0_) and (y1_ > y0_):
        rat = float(np.clip((y1_ - y0_) / (x1_ - x0_), 1.0 / 3.0, 3.0))
    else:
        rat = 1.0
    pw = 4.0 if rat >= 1.0 else 4.0 / rat       # panel width  (in)
    ph = 4.0 * rat if rat >= 1.0 else 4.0       # panel height (in)

    # --render-rows: one figure PER TIME, directories tiled into `render_rows`
    # rows (rest as columns), instead of the single dirs-x-times grid.
    if render_rows:
        import math
        for ti in range(ntimes):
            cells = [(r["label"], r["panels"][ti], r["box"], r["pt"],
                      (times[ti] if times is not None else r["times"][ti]),
                      r["contours"][ti])
                     for r in rows if ti < len(r["panels"])]
            if not cells:
                continue
            nd = len(cells)
            nr = max(1, min(render_rows, nd))
            nc = math.ceil(nd / nr)
            fig, axes = plt.subplots(nr, nc, squeeze=False,
                                     figsize=(pw * nc, (ph + 0.2) * nr),
                                     gridspec_kw=gkw)
            mm = None
            for idx, (lab, panel, box, pt, _tv, cpan) in enumerate(cells):
                ax = axes[idx // nc][idx % nc]
                mm = draw(ax, panel, box, pt, cpanel=cpan)
                ax.set_title(lab, fontsize=10)
                ax.set_xlabel(a1name)
                ax.set_ylabel(a2name)
            for idx in range(nd, nr * nc):          # blank unused cells
                axes[idx // nc][idx % nc].axis("off")
            if mm is not None:
                fig.colorbar(mm, ax=axes.ravel().tolist(), label=spec.label,
                             fraction=0.046, pad=0.04)
            _tile_axes(axes, same_x, same_y)
            tv = cells[0][4]
            fig.suptitle(f"{time_label} = {tv:.3g}  —  {spec.label}  "
                         f"({backend_name}, {proj}{csuffix})")
            if save:
                out = f"{save}_render_t{tv:g}.png"
                fig.savefig(out, bbox_inches="tight")
                print(f"saved render -> {out}")
            else:
                plt.show()
        return

    if len(rows) == 1:
        # Single run: one row, snapshots as columns (titles = time).
        r = rows[0]
        k = len(r["panels"])
        fig, axes = plt.subplots(1, k, figsize=(pw * k, ph + 0.6),
                                 squeeze=False, gridspec_kw=gkw)
        for c, (ax, panel) in enumerate(zip(axes[0], r["panels"])):
            m = draw(ax, panel, r["box"], r["pt"], cpanel=r["contours"][c])
            ax.set_title(f"{time_label} = {r['times'][c]:.3g}")
            ax.set_xlabel(a1name)
            ax.set_ylabel(a2name)
        fig.suptitle(f"{r['label']}  —  {spec.label}  "
                     f"({backend_name}, {proj}{csuffix})")
    else:
        # Several runs: DIRECTORIES as columns, TIMES as rows.
        ndirs = len(rows)
        fig, axes = plt.subplots(ntimes, ndirs, squeeze=False,
                                 figsize=(pw * ndirs, (ph + 0.2) * ntimes),
                                 gridspec_kw=gkw)
        for c, r in enumerate(rows):
            for rr in range(ntimes):
                ax = axes[rr][c]
                if rr >= len(r["panels"]):
                    ax.axis("off")
                    continue
                m = draw(ax, r["panels"][rr], r["box"], r["pt"],
                         cpanel=r["contours"][rr])
                if rr == 0:
                    ax.set_title(r["label"], fontsize=10)
                if c == 0:
                    rowt = (times[rr] if times is not None
                            else rows[0]["times"][rr])
                    ax.set_ylabel(f"{time_label} = {rowt:.3g}")
                if rr == ntimes - 1:
                    ax.set_xlabel(a1name)
        fig.suptitle(f"{spec.label}  ({backend_name}, {proj}, "
                     f"{a1name}{a2name}-plane{csuffix})")

    if m is not None:
        fig.colorbar(m, ax=axes.ravel().tolist(), label=spec.label,
                     fraction=0.046, pad=0.04)

    # Tile: drop inner tick labels on shared axes so panels read as one multiplot.
    _tile_axes(axes, same_x, same_y)

    if save:
        out = f"{save}_render.png"
        fig.savefig(out, bbox_inches="tight")
        print(f"saved render -> {out}")
    else:
        plt.show()

# --------------------------------------------------------------------------- #
#  Streamlines (2-D vector fields: velocity, magnetic field, ...)
# --------------------------------------------------------------------------- #

def vector_grid(fn, qx, qy, axis1, axis2, res, extent):
    """Interpolate a per-particle 2-D vector field onto a regular res x res grid
    over `extent` for matplotlib's streamplot.

    qx/qy follow the StreamSpec convention (tgdata column index or callable
    quantity(fn, tgdata, ngas)). Linear scattered interpolation, with the
    convex-hull edges filled by nearest-neighbour (streamplot needs a NaN-free
    regular grid). Returns (gx, gy, QX, QY) with QX/QY shaped [ny, nx].
    Generalized from IC_analysis_currentsheet's _bfield_grid."""
    from scipy.interpolate import griddata
    tgdata, _td, _ts, hdr, _tm = _cached_readtipsy(fn)
    ngas = int(hdr[2])
    vx = _quantity_values(qx, fn, tgdata, ngas)
    vy = _quantity_values(qy, fn, tgdata, ngas)
    pts = np.column_stack([np.asarray(tgdata[:, axis1], dtype=float),
                           np.asarray(tgdata[:, axis2], dtype=float)])
    x0, x1, y0, y1 = extent
    gx = np.linspace(x0, x1, res)
    gy = np.linspace(y0, y1, res)
    GX, GY = np.meshgrid(gx, gy)                      # shape (ny, nx)

    def interp(q):
        gi = griddata(pts, q, (GX, GY), method="linear")
        nan = np.isnan(gi)
        if np.any(nan):                               # fill hull edges with nearest
            gi[nan] = griddata(pts, q, (GX, GY), method="nearest")[nan]
        return gi

    return gx, gy, interp(vx), interp(vy)


def streamline_panels(entries, spec, times=None, n=4, res=64, save=None,
                      extent=None, time_of=None, time_label="t",
                      title=None, decorate=None):
    """Streamline figures of a 2-D vector field (StreamSpec): ONE figure per
    input, snapshot times as columns, lines colored by the in-plane magnitude
    on a white background.

    entries  : [(label, input_arg)], like render_panels.
    times/n  : snapshot selection (nearest to `times`, else n evenly spaced).
    res      : streamplot grid resolution (moderate, ~64, reads best).
    extent   : explicit [x0 x1 y0 y1]; default = the run's `.log` BOX
               PARAMETERS, else a robust percentile fallback (the same domain
               priority as the field-map render path).
    title    : figure heading (default '<label> streamlines').
    decorate : optional callable decorate(ax) drawn on every panel, for
               per-test annotations (e.g. currentsheet's initial sheet lines).
    Saves '<save>_<lab>.png' per input when there are several inputs, else
    '<save>.png' (the caller picks the stem, e.g. '<out>_stream_velocity');
    with no `save` the figures are shown interactively."""
    import matplotlib.pyplot as plt
    setup_rcparams()
    for lab, inp in entries:
        try:
            files = find_files(inp)
        except FileNotFoundError as e:
            print(f"stream[{lab}]: {e}", file=sys.stderr)
            continue
        chosen = select_snapshots(files, times=times, n=n, time_of=time_of)
        if not chosen:
            print(f"stream[{lab}]: no snapshots for '{inp}'", file=sys.stderr)
            continue

        # Domain box: explicit extent (--render-extent, else the spec's own
        # default) > .log BOX PARAMETERS > percentiles.
        tg0, *_ = _cached_readtipsy(chosen[0][0])
        ext = extent if extent is not None else spec.extent
        periods = None if ext is not None else read_box_periods(inp)
        if ext is not None:
            x0, x1, y0, y1 = ext
            src = "render extent"
        elif periods is not None:
            x0, x1 = -0.5 * periods[spec.axis1], 0.5 * periods[spec.axis1]
            y0, y1 = -0.5 * periods[spec.axis2], 0.5 * periods[spec.axis2]
            src = "log BOX PARAMETERS"
        else:
            def _pad(lo, hi):
                mid, half = 0.5 * (lo + hi), 0.51 * (hi - lo)
                return mid - half, mid + half
            x0, x1 = _pad(*np.nanpercentile(tg0[:, spec.axis1], [0.5, 99.5]))
            y0, y1 = _pad(*np.nanpercentile(tg0[:, spec.axis2], [0.5, 99.5]))
            src = "inferred (percentile)"
        print(f"stream[{lab}]: domain {src}: "
              f"{'xyz'[spec.axis1-1]}=[{x0:.4g},{x1:.4g}] "
              f"{'xyz'[spec.axis2-1]}=[{y0:.4g},{y1:.4g}]")

        nt = len(chosen)
        fig, axes = plt.subplots(1, nt, figsize=(4.2 * nt, 4.6), squeeze=False)
        fig.patch.set_facecolor("white")
        axes = axes[0]
        strm = None
        for ax, (fn, tt) in zip(axes, chosen):
            ax.set_facecolor("white")
            ax.set_aspect("equal")
            ax.set_xlim(x0, x1)
            ax.set_ylim(y0, y1)
            ax.set_title(f"{time_label} = {tt:.3g}")
            ax.set_xlabel("xyz"[spec.axis1 - 1])
            if decorate is not None:
                decorate(ax)
            gx, gy, QX, QY = vector_grid(fn, spec.qx, spec.qy, spec.axis1,
                                         spec.axis2, int(res), (x0, x1, y0, y1))
            mag = np.hypot(QX, QY)
            if not np.any(mag > 0):           # e.g. missing aux -> zero field
                ax.text(0.5, 0.5, f"zero {spec.label}", ha="center",
                        va="center", transform=ax.transAxes)
                continue
            strm = ax.streamplot(gx, gy, QX, QY, color=mag, cmap=spec.cmap,
                                 density=spec.density, linewidth=0.8,
                                 arrowsize=0.7, zorder=2)
        axes[0].set_ylabel("xyz"[spec.axis2 - 1])
        if strm is not None:
            fig.colorbar(strm.lines, ax=list(axes),
                         label=f"|{spec.label}|", fraction=0.046, pad=0.02)
        fig.suptitle(f"{title or f'{spec.label} streamlines'}  ({lab})")
        outname = None
        if save:
            outname = (f"{save}_{lab}.png" if len(entries) > 1
                       else f"{save}.png")
        finish_figure(fig, save=outname)


# --------------------------------------------------------------------------- #
#  Movies (animate a render over all snapshots)
# --------------------------------------------------------------------------- #

def _color_scale_from_pool(pool, spec, params, proj):
    """(norm, lo, hi, use_log, cmap) from pooled in-domain field values -- the
    same color-scale rule render_panels uses (spec.clim natural range > robust
    1/99 percentiles; symmetric for signed fields; degenerate -> linear [0,1])."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    use_log = spec.log
    if not pool:
        cl = spec.clim(params, proj) if spec.clim is not None else None
        vmin, vmax = cl if cl is not None else (0.0, 1.0)
        use_log = False
    else:
        vals = np.concatenate(pool)
        vmin, vmax = np.percentile(vals, 1), np.percentile(vals, 99)
        if spec.clim is not None:
            cl = spec.clim(params, proj)
            if cl is not None:
                vmin, vmax = cl
        if spec.symmetric:
            a = max(abs(vmin), abs(vmax))
            vmin, vmax = -a, a
        if vmax <= vmin:
            eps = max(abs(vmin), 1.0)
            vmin, vmax = vmin - eps, vmax + eps
    norm = LogNorm(vmin=max(vmin, 1e-30), vmax=vmax) if use_log else None
    lo = None if use_log else vmin
    hi = None if use_log else vmax
    cmap = plt.get_cmap(spec.cmap).copy()
    cmap.set_bad(cmap(0.0))
    return norm, lo, hi, use_log, cmap


def _resolve_movie_writer(fps, fmt="auto"):
    """Pick an animation writer: ffmpeg (mp4) preferred, Pillow (gif) fallback.
    Returns (writer, extension) or (None, None) if neither is available."""
    import matplotlib.animation as anim
    if fmt in ("auto", "mp4") and anim.FFMpegWriter.isAvailable():
        return anim.FFMpegWriter(fps=fps), "mp4"
    if fmt in ("auto", "gif") and anim.PillowWriter.isAvailable():
        return anim.PillowWriter(fps=fps), "gif"
    if fmt == "mp4" and anim.PillowWriter.isAvailable():
        print("movie: ffmpeg unavailable; falling back to gif.", file=sys.stderr)
        return anim.PillowWriter(fps=fps), "gif"
    return None, None


def render_movie(input_arg, label, spec, save=None, res=512, backend="grid",
                 extent=None, project=None, slice_frac=None, time_to_x=None,
                 time_label="t", params=None, aspect="equal", nsmooth=None,
                 render_reference=None, fps=10, fmt="auto"):
    """Animate `spec.quantity` over ALL snapshots of one run as a movie.

    Each frame is the field at one snapshot time (single panel). With
    `render_reference` AND every frame time present in that baseline, each frame
    becomes sim | reference | (sim - reference) sharing the field scale (the diff
    on its own symmetric scale); if any frame time is missing from the reference
    it falls back to a simulation-only movie with a warning. mp4 via ffmpeg, else
    gif via Pillow."""
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    mode = _resolve_render_mode(backend, spec)
    time_of = ((lambda fn: time_to_x(_read_time(fn), params)) if time_to_x
               else _read_time)
    proj = project or spec.project
    times = snapshot_times(input_arg, time_of=time_of)
    if not times:
        print(f"movie: no snapshots for '{input_arg}'", file=sys.stderr)
        return
    r = _render_one_input(input_arg, label, spec, times, len(times), res, mode,
                          extent, project, slice_frac, time_of, nsmooth=nsmooth)
    if r is None:
        return
    panels, box, ptimes = r["panels"], r["box"], r["times"]
    if any(p[0] != "image" for p in panels):
        print("movie: needs image panels (grid/pynbody/hist backend, not "
              "particles); skipping.", file=sys.stderr)
        return
    nframes = len(panels)

    # Optional reference: require EVERY frame time to have a matching grid.
    ref, ref_match = None, None
    if render_reference is not None:
        ref = load_render_grids(render_reference)
        if ref is not None:
            rt = np.array([g[2] for g in ref["grids"]], dtype=float)
            if rt.size > 1:
                sp = np.diff(np.sort(rt))
                sp = sp[sp > 0]
                ttol = 0.5 * float(sp.min()) + 1e-6 if sp.size else 1e-6
            else:
                ttol = max(1e-6, 1e-3 * (1.0 + abs(rt[0])))
            match = [int(np.argmin(np.abs(rt - t))) for t in ptimes]
            ok = [abs(rt[j] - t) <= ttol for j, t in zip(match, ptimes)]
            if all(ok):
                ref_match = match
            else:
                miss = [f"{t:g}" for t, o in zip(ptimes, ok) if not o]
                print(f"movie: reference is missing time(s) {', '.join(miss)} -- "
                      f"making a simulation-only movie (bless with "
                      f"--save-reference-all to cover every dump).",
                      file=sys.stderr)
                ref = None

    pool = []
    for p in panels:
        v = _crop_to_view(p[1], p[2], *box)
        v = np.asarray(v, dtype=float).ravel()
        v = v[np.isfinite(v)]
        if spec.log:
            v = v[v > 0]
        if v.size:
            pool.append(v)
    if ref is not None:
        for im, _e, _t in ref["grids"]:
            v = np.asarray(im, dtype=float).ravel()
            v = v[np.isfinite(v)]
            if spec.log:
                v = v[v > 0]
            if v.size:
                pool.append(v)
    norm, lo, hi, _ul, cmap = _color_scale_from_pool(pool, spec, params, proj)

    writer, ext = _resolve_movie_writer(fps, fmt)
    if writer is None:
        print("movie: no animation writer (ffmpeg/Pillow) available.",
              file=sys.stderr)
        return

    setup_rcparams()
    a1name, a2name = "xyz"[spec.axis1 - 1], "xyz"[spec.axis2 - 1]

    if ref_match is not None:
        diffs, dvals = [], []
        for i in range(nframes):
            sim = np.asarray(panels[i][1], dtype=float)
            rr = np.asarray(ref["grids"][ref_match[i]][0], dtype=float)
            d = sim - rr if sim.shape == rr.shape else None
            diffs.append(d)
            if d is not None:
                dd = d[np.isfinite(d)]
                if dd.size:
                    dvals.append(np.abs(dd))
        damax = float(np.concatenate(dvals).max()) if dvals else 1.0
        if damax <= 0:
            damax = 1.0
        diffcmap = plt.get_cmap("RdBu_r").copy()
        diffcmap.set_bad(diffcmap(0.5))
        z0 = np.zeros_like(panels[0][1])

        fig = plt.figure(figsize=(13, 4.8))
        # Spacer column (gs[*,3]) after the middle field colorbar so its right-side
        # label doesn't overlap the diff panel (see render_panels residual layout).
        gs = fig.add_gridspec(1, 6, width_ratios=[1, 1, 0.06, 0.16, 1, 0.06],
                              wspace=0.07)
        ax0, ax1, ax2 = (fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]),
                         fig.add_subplot(gs[0, 4]))
        caxf, caxd = fig.add_subplot(gs[0, 2]), fig.add_subplot(gs[0, 5])
        im0 = ax0.imshow(panels[0][1], extent=panels[0][2], origin="lower",
                         cmap=cmap, norm=norm, vmin=lo, vmax=hi, aspect=aspect)
        im1 = ax1.imshow(ref["grids"][ref_match[0]][0],
                         extent=ref["grids"][ref_match[0]][1], origin="lower",
                         cmap=cmap, norm=norm, vmin=lo, vmax=hi, aspect=aspect)
        im2 = ax2.imshow(diffs[0] if diffs[0] is not None else z0,
                         extent=panels[0][2], origin="lower", cmap=diffcmap,
                         vmin=-damax, vmax=damax, aspect=aspect)
        for ax, ttl in ((ax0, "simulation"), (ax1, "reference"),
                        (ax2, "sim − reference")):
            ax.set_xlim(box[0], box[1])
            ax.set_ylim(box[2], box[3])
            ax.set_aspect(aspect)
            ax.set_title(ttl, fontsize=11)
            ax.set_xlabel(a1name)
        ax0.set_ylabel(a2name)
        ax1.tick_params(labelleft=False)
        ax2.tick_params(labelleft=False)
        fig.colorbar(im0, cax=caxf, label=spec.label)
        fig.colorbar(im2, cax=caxd, label=f"Δ {spec.label}")
        sttl = fig.suptitle("")

        def update(i):
            im0.set_data(panels[i][1])
            im1.set_data(ref["grids"][ref_match[i]][0])
            im2.set_data(diffs[i] if diffs[i] is not None else z0)
            sttl.set_text(f"{label}  vs reference  —  {spec.label}  "
                          f"({time_label} = {ptimes[i]:.3g})")
            return im0, im1, im2, sttl
    else:
        fig, ax = plt.subplots(figsize=(6.2, 5.6))
        im0 = ax.imshow(panels[0][1], extent=panels[0][2], origin="lower",
                        cmap=cmap, norm=norm, vmin=lo, vmax=hi, aspect=aspect)
        ax.set_xlim(box[0], box[1])
        ax.set_ylim(box[2], box[3])
        ax.set_aspect(aspect)
        ax.set_xlabel(a1name)
        ax.set_ylabel(a2name)
        fig.colorbar(im0, ax=ax, label=spec.label, fraction=0.046, pad=0.04)
        sttl = fig.suptitle("")

        def update(i):
            im0.set_data(panels[i][1])
            sttl.set_text(f"{label}  —  {spec.label}  "
                          f"({time_label} = {ptimes[i]:.3g})")
            return im0, sttl

    ani = animation.FuncAnimation(fig, update, frames=nframes, blit=False)
    if save:
        out = f"{save}.{ext}"
        ani.save(out, writer=writer, dpi=120)
        print(f"saved movie -> {out}  ({nframes} frames @ {fps} fps)")
        plt.close(fig)
    else:
        plt.show()


# --------------------------------------------------------------------------- #
#  Top-level CLI
# --------------------------------------------------------------------------- #

