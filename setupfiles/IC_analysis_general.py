import argparse
import glob
import os
import re
import sys
from collections import namedtuple
from functools import lru_cache
from typing import List, Tuple
import numpy as np
import readtipsy as tip
import IC_smoothlength as smth


# ----------------------- aux-file readers ------------------------------

@lru_cache(maxsize=256)
def _cached_readtipsyaux(snap, label):
    """Memoized tip.readtipsyaux (raw disk read). The render path reads the same
    aux several times for one snapshot -- e.g. orzag reads BField once for |B|
    and again for divBerr, x3 planes. Aux files are immutable during an analysis,
    so cache by (path, label); try_aux hands out a fresh writeable copy so callers
    keep their old (mutable, independent) contract. Raises are NOT cached."""
    arr, _ = tip.readtipsyaux(snap, label)
    arr = np.asarray(arr, dtype=float)
    arr.flags.writeable = False        # protect the shared cache; try_aux copies
    return arr


def try_aux(snap, label, n=None):
    """Read a tipsy aux file, returning a 1-D float array or None if absent.

    `snap`  : path to the snapshot (the .label suffix is added internally).
    `label` : aux label, e.g. "BField", "divB", "BClean".
    `n`     : if given, truncate to the first n entries (typically ngas).
    """
    try:
        arr = _cached_readtipsyaux(snap, label)
    except (FileNotFoundError, OSError):
        return None
    if n is not None:
        arr = arr[:n]
    return np.array(arr)               # fresh writeable copy (cache stays clean)


def read_vector_aux(snap, label, ngas):
    """Read a per-particle vector aux field, returning shape (ngas, 3) or None.

    Tries the packed form (label, length 3*ngas) first, then falls back to
    component files (label+"x", label+"y", label+"z").
    """
    arr = try_aux(snap, label)
    if arr is not None and arr.size >= 3 * ngas:
        return arr[:3 * ngas].reshape(ngas, 3)
    cx = try_aux(snap, label + "x", ngas)
    cy = try_aux(snap, label + "y", ngas)
    cz = try_aux(snap, label + "z", ngas)
    if cx is None or cy is None or cz is None:
        return None
    return np.stack([cx, cy, cz], axis=1)


def read_bfield(snap, ngas):
    """Convenience: read the magnetic field aux. See read_vector_aux."""
    return read_vector_aux(snap, "BField", ngas)


# ----------------------- run .log parameters ---------------------------

_LOG_PARAM_CACHE = {}


def _find_logs(run_arg):
    """The gasoline '<runname>.log' file(s) for a run directory or file-prefix."""
    if os.path.isdir(run_arg):
        return sorted(glob.glob(os.path.join(run_arg, "*.log")))
    cands = sorted(glob.glob(run_arg + "*.log"))
    if os.path.isfile(run_arg + ".log"):
        cands.append(run_arg + ".log")
    return cands


def read_log_param(run_arg, key, default=None):
    """Read a numeric '<key>: <value>' parameter from the run's gasoline '.log'.

    gasoline records its run parameters in '<runname>.log' as '… key: value …'.
    `run_arg` is the run directory or file-prefix. Returns the float value, or
    `default` if no '.log' / key is found. Cached per (run, key)."""
    ck = (os.path.normpath(run_arg), key)
    if ck in _LOG_PARAM_CACHE:
        return _LOG_PARAM_CACHE[ck]
    val = default
    for log in _find_logs(run_arg):
        try:
            with open(log) as fh:
                txt = fh.read()
        except OSError:
            continue
        m = re.search(rf"\b{re.escape(key)}:\s*([-\d.eE+]+)", txt)
        if m:
            val = float(m.group(1))
            break
    _LOG_PARAM_CACHE[ck] = val
    return val


def read_tufac(run_arg, default=1.0):
    """dTuFac from the run '.log': u = dTuFac * T converts the tipsy 'temperature'
    column (tgdata col 8) to specific internal energy.

    gasoline stores specific internal energy as a temperature; dTuFac is the
    conversion (dGasConst = (gamma-1) dTuFac in the same log). Returns `default`
    (1.0, i.e. u=T) with a warning if it can't be determined.

    Codes that don't log `dTuFac` directly (e.g. ChaNGa) still log the ideal-gas
    constituents, so we reconstruct it from the same relation:
        dTuFac = dGasConst / (dConstGamma - 1) / dMeanMolWeight."""
    val = read_log_param(run_arg, "dTuFac", None)
    if val is not None:
        return val
    # Backup: derive from the ideal-gas params (dGasConst = (gamma-1) mu dTuFac).
    gc = read_log_param(run_arg, "dGasConst", None)
    gamma = read_log_param(run_arg, "dConstGamma", None)
    mu = read_log_param(run_arg, "dMeanMolWeight", 1.0)
    if gc is not None and gamma is not None and gamma != 1.0 and mu:
        return gc / (gamma - 1.0) / mu
    print(f"read_tufac: dTuFac not found (and dGasConst/dConstGamma absent) in a "
          f".log for '{run_arg}'; assuming u=T (factor {default}).", file=sys.stderr)
    return default


# Code->physical (cgs) unit conversion for a gasoline run (G=1 system).
_G_CGS = 6.67430e-8                      # gravitational constant, cgs
_MSOL_G = 1.989e33                       # solar mass, g
_KPC_CM = 3.0856775814913673e21         # kpc, cm
_MYR_S = 3.1557e13                       # Myr, s

CodeUnits = namedtuple("CodeUnits", [
    "msol", "kpc",              # the raw .log dMsolUnit / dKpcUnit
    "mass_g", "length_cm", "time_s", "time_myr",
    "velocity_kms", "density_gcc", "erg_per_g",
    "b_gauss", "b_uG",          # magnetic field unit (ideal-MHD code->Gauss)
])


def code_units(run_arg):
    """Code->cgs unit factors for a gasoline run (G=1 system), or None.

    Derived from the run '.log' dMsolUnit (M, in Msol) and dKpcUnit (L, in kpc),
    which fix every other unit when G=1: V = sqrt(G M / L), T = L/V, rho = M/L^3,
    specific energy = V^2. The MAGNETIC unit follows the ideal-MHD (Gaussian)
    Alfven relation v_A = B/sqrt(4 pi rho) -> the code B unit corresponds to V, so
    1 code B = V_cgs * sqrt(4 pi * rho_cgs) Gauss. (The snapshot BField is in this
    internal code unit, NOT Gauss, even when the IC was specified in Gauss.)

    Returns a `CodeUnits` namedtuple; multiply a code-unit quantity by the matching
    field (density_gcc, velocity_kms, b_uG, ...) to get cgs. None if dMsolUnit /
    dKpcUnit are absent from the .log."""
    msol = read_log_param(run_arg, "dMsolUnit")
    kpc = read_log_param(run_arg, "dKpcUnit")
    if msol is None or kpc is None:
        return None
    M, L = msol * _MSOL_G, kpc * _KPC_CM
    v_cms = np.sqrt(_G_CGS * M / L)
    rho_gcc = M / L ** 3
    b_gauss = v_cms * np.sqrt(4.0 * np.pi * rho_gcc)
    return CodeUnits(
        msol=msol, kpc=kpc, mass_g=M, length_cm=L,
        time_s=L / v_cms, time_myr=(L / v_cms) / _MYR_S,
        velocity_kms=v_cms / 1.0e5, density_gcc=rho_gcc, erg_per_g=v_cms ** 2,
        b_gauss=b_gauss, b_uG=b_gauss * 1.0e6,
    )


# Canonical per-step column names -- the vocabulary the analysis scripts request
# from `read_log_series`. Both gasoline and ChaNGa emit a NAMED column-header line
# in their '.log' (see `_parse_log_header`), so the series is parsed BY NAME; this
# survives column additions/reordering across code versions and lets a different
# code (ChaNGa) expose a different subset. LOG_COLUMNS below is only the LEGACY
# positional fallback, used when a log has no recognisable header line (e.g. the
# smoke-test synthetic logs, or very old runs).
#
# gasoline's header is '# [dTime] [z] [Etot] [Ekin] [Epot] [Eth] [Emag] [EClean]
# [totentrop] [totenstro] [Lx] [Ly] [Lz] ...' (the positional map below predates
# the [EClean] insertion, so it is misaligned from totentrop onward for current
# builds -- another reason to parse by name).
LOG_COLUMNS = {
    "dTime": 0, "z": 1, "Etot": 2, "Ekin": 3, "Epot": 4, "Eth": 5, "Emag": 6,
    "totentrop": 7, "totenstro": 8, "Lx": 9, "Ly": 10, "Lz": 11,
    "Llinx": 12, "Liny": 13, "Llinz": 14, "cmx": 15, "cmy": 16, "cmz": 17,
    "MWxy": 18, "MWyz": 19, "MWzx": 20, "RSxy": 21, "RSyz": 22, "RSzx": 23,
    "divBAvg": 24, "divBMax": 25, "divBerrAvg": 26, "divBerrMax": 27,
    "alphaAvg": 28, "alphaMax": 29, "betaMax": 30, "betaAvg": 31, "betaMin": 32,
    "etaresAvg": 33, "kinviscAvg": 34, "rmsmach": 35, "vrms": 36,
    "rhogasAvg": 37, "rhogasMax": 38, "Q1Avg": 39, "Q1Max": 40, "Q2Avg": 41,
    "Q2Max": 42, "E0Avg": 43, "E0Max": 44, "Q4Avg": 45, "Q4Max": 46,
    "WallTime": 47, "dWMax": 48, "dImax": 49, "dEMax": 50, "dMultiEff": 51,
}

# Foreign (per-code) column-name -> canonical name. gasoline already writes the
# canonical names (its header tokens are just '[name]'), so only codes with a
# different vocabulary need entries here. Add a code by adding its column names.
LOG_NAME_ALIASES = {
    # --- ChaNGa (Main::writeOutputLog header: "time redshift TotalEVir TotalE
    #     Kinetic Virial Potential TotalECosmo Ethermal Lx Ly Lz Wallclock") ---
    "time": "dTime", "redshift": "z", "TotalE": "Etot", "Kinetic": "Ekin",
    "Potential": "Epot", "Ethermal": "Eth", "Wallclock": "WallTime",
    # ChaNGa-only columns with no gasoline analogue (TotalEVir, Virial,
    # TotalECosmo) pass through unmapped and are simply not requestable.
}


def _canon_log_name(tok):
    """Canonical column name for a raw header token: strips gasoline's '[ ]'
    ('[Etot]' -> 'Etot') and applies LOG_NAME_ALIASES (ChaNGa 'Kinetic' -> 'Ekin').
    Unknown tokens pass through unchanged."""
    tok = tok.strip().strip("[]")
    return LOG_NAME_ALIASES.get(tok, tok)


def _parse_log_header(line):
    """If `line` is a per-step column header, return {canonical_name: index}; else
    None. A header is a '#' comment line whose first column is the time -- this
    recognises both the gasoline '# [dTime] [z] [Etot] ...' and the ChaNGa
    '# time redshift TotalEVir TotalE ...' forms while rejecting prose/param lines."""
    s = line.strip()
    if not s.startswith("#"):
        return None
    toks = s.lstrip("#").split()
    if not toks:
        return None
    names = [_canon_log_name(t) for t in toks]
    if names[0] != "dTime":
        return None
    return {name: i for i, name in enumerate(names)}


def read_log_series(run_arg, columns):
    """Per-timestep series from the run '.log' step rows, code-agnostically.

    A code (gasoline, ChaNGa, ...) appends one row per timestep to '<runname>.log',
    so the energy/divergence diagnostics there are sampled far more densely than
    the periodic snapshot dumps -- the preferred source for any quantity-vs-time
    plot. The columns are mapped BY NAME from each log's header line (see
    `_parse_log_header`); a log with no recognisable header falls back to the
    legacy positional LOG_COLUMNS.

    `columns` is a sequence of canonical names (from LOG_COLUMNS; include "dTime"
    for the time axis, which is always columns[0]). Returns a dict name ->
    np.ndarray, rows sorted by time, or None if no '.log' yields all requested
    columns (e.g. a ChaNGa log has no Emag/divB) -- letting callers fall back to
    per-snapshot computation via `series_from_log_or_snapshots`. Comment / blank /
    RESTART ('#') lines and short rows are skipped; a restart-appended header
    re-maps the columns for the rows that follow it."""
    rows = []
    for log in _find_logs(run_arg):
        try:
            fh = open(log)
        except OSError:
            continue
        with fh:
            colmap = None      # name->index from this log's header (None = none yet)
            idx = None         # resolved indices for `columns`; False = unsatisfiable
            for line in fh:
                hdr = _parse_log_header(line)
                if hdr is not None:
                    colmap, idx = hdr, None   # re-resolve against the new header
                    continue
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                if idx is None:
                    cmap = colmap if colmap is not None else LOG_COLUMNS
                    try:
                        idx = [cmap[c] for c in columns]
                    except KeyError:
                        idx = False           # a requested column is absent here
                if idx is False:
                    break                     # this log can't satisfy the request
                parts = s.split()
                if len(parts) <= max(idx):
                    continue
                try:
                    vals = [float(parts[i]) for i in idx]
                except ValueError:
                    continue
                rows.append((vals[0], vals))  # columns[0] is the time axis
    if not rows:
        return None
    rows.sort(key=lambda r: r[0])
    arr = np.array([r[1] for r in rows], dtype=float)
    return {name: arr[:, k] for k, name in enumerate(columns)}


# ----------------------- Fourier-mode amplitude ------------------------

def sin_cos_amp(mass, x, f, k, x0):
    """Project f onto sin(k(x-x0)) and cos(k(x-x0)), mass-weighted.

    For a perturbation f = As sin(k(x-x0)) + Ac cos(k(x-x0)) sampled on a
    mass-weighted particle distribution,
        As = 2 <m f sin(.)> / <m>,  Ac = 2 <m f cos(.)> / <m>
    recovers the modal coefficients. Works for any 1-D wave decomposition;
    pass a different coordinate (y, r, ...) for waves along other directions.
    """
    phase = k * (np.asarray(x) - x0)
    m = np.asarray(mass)
    W = float(np.sum(m))
    As = 2.0 * float(np.sum(m * f * np.sin(phase))) / W
    Ac = 2.0 * float(np.sum(m * f * np.cos(phase))) / W
    return As, Ac


def amp_mag(mass, x, f, k, x0):
    """Magnitude of the (sin, cos) Fourier amplitude pair."""
    As, Ac = sin_cos_amp(mass, x, f, k, x0)
    return float(np.hypot(As, Ac))


# ----------------------- plotting helpers ------------------------------

def setup_rcparams():
    """Uniform matplotlib styling for the analysis scripts."""
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.size": 14, "axes.labelsize": 16, "axes.titlesize": 16,
        "legend.fontsize": 12, "xtick.labelsize": 12, "ytick.labelsize": 12,
        "lines.linewidth": 2.0, "lines.markersize": 5,
    })


def finish_figure(fig, save=None, dpi=150):
    """Save the figure to `save` if given, otherwise show interactively."""
    import matplotlib.pyplot as plt
    if save:
        fig.savefig(save, bbox_inches='tight', dpi=dpi)
        print(f"saved {save}")
    else:
        plt.show()


def finish_figure_with_legend(fig, ax, save=None, dpi=150):
    """Like `finish_figure`, but the legend (built from `ax`'s labelled artists)
    is written to a SEPARATE figure rather than drawn on the axes -- the project
    convention so the plot stays clean (publication-friendly).

    `save` is the figure path (e.g. 'out_profile.png'); the legend is saved
    alongside as '<stem>_legend.png'. When `ax` is a sequence of axes, the labels
    are pooled across them (de-duplicated, first occurrence wins). With no `save`,
    the legend is shown as its OWN figure window too -- drawing it on the axes
    (the old behavior) covered in-axes annotations like the L1-error boxes
    (loc="best" avoids data artists but not text)."""
    import matplotlib.pyplot as plt
    axes = ax if isinstance(ax, (list, tuple, np.ndarray)) else [ax]
    handles, labels = [], []
    seen = set()
    for a in axes:
        for h, l in zip(*a.get_legend_handles_labels()):
            if l and l not in seen:
                seen.add(l)
                handles.append(h)
                labels.append(l)

    if not save:
        if handles:
            legfig = plt.figure(figsize=(6, 0.5 + 0.35 * len(labels)))
            legfig.legend(handles, labels, loc="center",
                          ncol=min(3, len(labels)), frameon=False,
                          handlelength=2.5, handletextpad=0.8,
                          columnspacing=1.5)
            legfig.tight_layout(pad=0.1)
        plt.show()
        return

    legfig = None
    if handles:
        legfig = plt.figure(figsize=(6, 0.5 + 0.35 * len(labels)))
        legfig.legend(handles, labels, loc="center", ncol=min(3, len(labels)),
                      frameon=False, handlelength=2.5, handletextpad=0.8,
                      columnspacing=1.5)
        legfig.tight_layout(pad=0.1)
    fig.savefig(save, bbox_inches="tight", dpi=dpi)
    print(f"saved {save}")
    if legfig is not None:
        stem, _, ext = save.rpartition(".")
        legpath = f"{stem}_legend.{ext}" if stem else f"{save}_legend.png"
        legfig.savefig(legpath, bbox_inches="tight", dpi=dpi)
        print(f"saved {legpath}")
        plt.close(legfig)
    plt.close(fig)

def loaddata(entry: str):
    # tgdata structure mass,x,y,z,vx,vy,vz,rho,T,eps,metals,potential
    tgdata, tddata, tsdata, data_header, time = tip.readtipsy(entry)
    N = data_header[0]
    ngas = data_header[2]
    ndark = data_header[3]
    nstar = data_header[4]
    # compute smoothing length for gas (if you need it)
    h = smth.getsmooth2(tgdata[:, 0], tgdata[:, 7], 64)  # mass, dens, nsmooth
    return tgdata, tddata, tsdata, data_header, time, N, ngas, ndark, nstar, h


# ----------------------- helpers ---------------------------------------

from typing import List
import glob, os, re

def find_files(arg: str) -> List[str]:
    """Return a sorted list of filenames matching the given prefix/dir/pattern.

    Matches only files whose BASENAME ends exactly with a 5-digit suffix,
    e.g. 'prefix.00200' (no extra characters after the digits).
    """
    # treat an explicit glob pattern only if it contains wildcard characters
    if any(ch in arg for ch in "*?["):
        candidates = glob.glob(arg)
    else:
        # if argument is a directory, search that directory
        if os.path.isdir(arg):
            dirn = arg
            entries = os.listdir(dirn)
            candidates = [os.path.join(dirn, f) for f in entries]
        else:
            # treat arg as prefix: search in the arg's directory for names that
            # start with the basename of arg
            dirn = os.path.dirname(arg) or "."
            prefix = os.path.basename(arg)
            try:
                entries = os.listdir(dirn)
            except FileNotFoundError:
                candidates = []
            else:
                candidates = [os.path.join(dirn, f) for f in entries if f.startswith(prefix)]

    # Strictly match basenames that end with an optional dot + 5 digits and nothing after
    # Examples matched: prefix.00200, prefix00200 (the optional dot covers both cases)
    rx = re.compile(r"^(.+?)(?:\.)?(\d{5})$")

    files_with_suffix = []
    for fn in candidates:
        base = os.path.basename(fn)
        m = rx.match(base)
        if m:
            suffix = int(m.group(2))
            files_with_suffix.append((suffix, fn))

    if not files_with_suffix:
        raise FileNotFoundError(
            f"No files found matching pattern/prefix '{arg}' with a trailing 5-digit suffix."
        )

    files_with_suffix.sort(key=lambda x: x[0])
    return [fn for _, fn in files_with_suffix]


def find_files2(arg: str) -> List[str]:
    """Return a sorted list of filenames matching the given prefix/dir/pattern.

    The function accepts either:
      - a glob pattern (if it contains * or ?), or
      - a filename prefix or directory. In that case it will search the
        directory for files whose basename contains a *trailing* five-digit
        suffix (e.g. "00001") and return them sorted by that numeric suffix.
    """
    # treat an explicit glob pattern
    if any(ch in arg for ch in "*[0-9]"):
        candidates = glob.glob(arg)
    else:
        # if argument is a directory, search all files inside
        if os.path.isdir(arg):
            candidates = glob.glob(os.path.join(arg, "*[0-9]"))
        else:
            # treat arg as prefix: match anything that starts with it
            candidates = glob.glob(arg + "*[0-9]")

    # filter for names that contain a trailing 5-digit group
    files_with_suffix = []
    rx = re.compile(r"(\d{5})\D*$")
    for fn in candidates:
        base = os.path.basename(fn)
        m = rx.search(base)
        if m:
            suffix = int(m.group(1))
            files_with_suffix.append((suffix, fn))

    if not files_with_suffix:
        raise FileNotFoundError(
            f"No files found matching pattern/prefix '{arg}' with a trailing 5-digit suffix."
        )



import os, glob
import numpy as np
import imageio.v2 as imageio
from matplotlib.colors import LogNorm
from matplotlib import pyplot as plt

# --- helper to build a weighted 2D histogram (no plotting) ---
def _hist2d_weighted(axis1, axis2, render, box_size=1.0, res=512):
    x = np.asarray(axis1, float).ravel()
    y = np.asarray(axis2, float).ravel()
    w = np.asarray(render, float).ravel()
    n = min(len(x), len(y), len(w))
    x, y, w = x[:n], y[:n], w[:n]
    half = box_size/2.0
    m = (
        np.isfinite(x) & np.isfinite(y) & np.isfinite(w) &
        (x >= -half) & (x <= half) &
        (y >= -half) & (y <= half)
    )
    x, y, w = x[m], y[m], w[m]
    H, xedges, yedges = np.histogram2d(
        x, y,
        bins=res,
        range=[[-half, half], [-half, half]],
        weights=w
    )
    return H.T  # origin lower

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def render_particles(axis1, axis2, render, box_size=None, res=512, outname=None, log=True):
    """
    Render a weighted 2D histogram (projection) of `render` onto (axis1, axis2).

    If box_size is None, it is inferred from the data range of axis1 and axis2.
    """
    x = np.asarray(axis1, float).ravel()
    y = np.asarray(axis2, float).ravel()
    w = np.asarray(render, float).ravel()
    n = min(len(x), len(y), len(w))
    x, y, w = x[:n], y[:n], w[:n]

    # --- Infer box size automatically if not provided ---
    if box_size is None:
        # use the largest side length among x and y ranges, add 5% padding
        xmin, xmax = np.nanmin(x), np.nanmax(x)
        ymin, ymax = np.nanmin(y), np.nanmax(y)
        size_x = xmax - xmin
        size_y = ymax - ymin
        box_size = 1.05 * max(size_x, size_y)
        print(f"[info] auto box_size = {box_size:.4g} from data range")

    half = box_size / 2.0

    # center box on median of coordinates
    xmid = 0.5 * (np.nanmin(x) + np.nanmax(x))
    ymid = 0.5 * (np.nanmin(y) + np.nanmax(y))
    x -= xmid
    y -= ymid

    # mask to box limits
    m = (
        np.isfinite(x) & np.isfinite(y) & np.isfinite(w) &
        (np.abs(x) <= half) & (np.abs(y) <= half)
    )
    x, y, w = x[m], y[m], w[m]

    # weighted 2D histogram
    H, xedges, yedges = np.histogram2d(
        x, y,
        bins=res,
        range=[[-half, half], [-half, half]],
        weights=w
    )
    H = H.T  # origin at bottom-left for imshow

    # color scale
    if log:
        pos = H[H > 0]
        eps = np.nanmin(pos) if pos.size else 1e-12
        norm = LogNorm(vmin=eps, vmax=np.nanmax(pos) if pos.size else 1.0)
    else:
        norm = None

    plt.figure(figsize=(6, 6))
    im = plt.imshow(
        H,
        extent=(-half, half, -half, half),
        origin="lower",
        cmap="inferno",
        norm=norm
    )
    plt.xlabel("axis1")
    plt.ylabel("axis2")
    plt.colorbar(im, label="sum(render) per pixel")
    plt.tight_layout()

    if outname:
        plt.savefig(outname, dpi=150, bbox_inches="tight")
        plt.close()
    else:
        plt.show()

def infer_box_size_from_particles(tgdata):
    xmin, xmax = np.nanmin(tgdata[:,1]), np.nanmax(tgdata[:,1])
    ymin, ymax = np.nanmin(tgdata[:,2]), np.nanmax(tgdata[:,2])
    zmin, zmax = np.nanmin(tgdata[:,3]), np.nanmax(tgdata[:,3])
    size = 1.05 * max(xmax - xmin, ymax - ymin, zmax - zmin)
    return size


# --- MOVIE MAKER ---
def render_movie_over_files(files, axis1_idx, axis2_idx, render_idx,
                            box_size=1.0, res=512, fps=15, out="movie.mp4",
                            log=True, cmap="inferno", use_global_norm=True):
    """
    files: list of snapshot paths (sorted)
    axis*_idx, render_idx: integer column indices in tgdata (e.g. 1,2,7)
    """
    if not files:
        raise ValueError("No files provided")
    tgdata, *_ = loaddata(files[0])
    box_size = infer_box_size_from_particles(tgdata)
    # PASS 1: compute global vmin/vmax on the histograms (for stable colors)
    vmin = None; vmax = None
    if use_global_norm:
        vmins = []; vmaxs = []
        for fn in files:
            tgdata, *_ = loaddata(fn)
            H = _hist2d_weighted(tgdata[:,axis1_idx], tgdata[:,axis2_idx], tgdata[:,render_idx],
                                 box_size, res)
            if log:
                pos = H[H > 0]
                if pos.size:
                    vmins.append(np.percentile(pos, 1))   # robust low
                    vmaxs.append(np.percentile(pos, 99))  # robust high
            else:
                vmins.append(np.percentile(H, 1))
                vmaxs.append(np.percentile(H, 99))
        if vmins and vmaxs:
            vmin = float(np.min(vmins))
            vmax = float(np.max(vmaxs))
            if log and vmin <= 0: vmin = 1e-12

    # PASS 2: stream frames directly into a video
    with imageio.get_writer(out, fps=fps, codec="libx264", quality=8) as writer:
        for i, fn in enumerate(files):
            print(f"[{i+1}/{len(files)}] {fn}")
            try:
                tgdata, *_ = loaddata(fn)
            except Exception as e:
                print(f"  skip: {e}", file=sys.stderr)
                continue

            # make a frame to a temporary buffer (no disk)
            H = _hist2d_weighted(tgdata[:,axis1_idx], tgdata[:,axis2_idx], tgdata[:,render_idx],
                                 box_size, res)

            # draw with the same normalization
            plt.figure(figsize=(6,6))
            if log:
                pos = H[H > 0]
                vmin_eff = vmin if vmin is not None else (np.percentile(pos,1) if pos.size else 1e-12)
                vmax_eff = vmax if vmax is not None else (np.percentile(pos,99) if pos.size else 1.0)
                norm = LogNorm(vmin=vmin_eff, vmax=vmax_eff)
            else:
                norm = None

            im = plt.imshow(H,
                            extent=(-box_size/2, box_size/2, -box_size/2, box_size/2),
                            origin="lower", cmap=cmap, norm=norm)
            plt.axis('off'); plt.tight_layout(pad=0)

            # grab the canvas as an RGB array and append to video
            plt.draw()
            frame = np.frombuffer(plt.gcf().canvas.tostring_rgb(), dtype=np.uint8)
            frame = frame.reshape(plt.gcf().canvas.get_width_height()[::-1] + (3,))
            writer.append_data(frame)
            plt.close()

    print(f"Saved movie to {out}")

