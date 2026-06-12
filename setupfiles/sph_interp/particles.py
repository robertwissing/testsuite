"""Code-agnostic particle container — the canonical input for `sph_interp`.

Every simulation code converts *into* this minimal, physics-driven container
rather than the library speaking any one code's on-disk layout. It carries
exactly what interpolation needs and nothing else:

    pos   (N,3)  particle positions
    mass  (N,)   particle masses
    rho   (N,)   densities
    h     (N,)   smoothing lengths
    values dict[str, (N,)]  the per-particle field(s) to interpolate (by name)
    box   (3,)   periodic box size (period) for each axis
    periodic bool | (3,) bool  periodicity per axis

Per-code adapters build it; `from_tipsy` is the first (gasoline/tipsy). Others
(`from_gadget`, ...) slot in the same way.
"""

from dataclasses import dataclass, field
from typing import Dict, Union, Tuple

import numpy as np


@dataclass
class Particles:
    pos: np.ndarray                      # (N, 3)
    mass: np.ndarray                     # (N,)
    rho: np.ndarray                      # (N,)
    h: np.ndarray                        # (N,)
    values: Dict[str, np.ndarray] = field(default_factory=dict)
    box: np.ndarray = None               # (3,) period; None = non-periodic / unknown
    periodic: Union[bool, Tuple[bool, bool, bool]] = False

    def __post_init__(self):
        self.pos = np.asarray(self.pos, dtype=np.float64)
        if self.pos.ndim != 2 or self.pos.shape[1] != 3:
            raise ValueError("pos must have shape (N, 3)")
        n = self.pos.shape[0]
        self.mass = self._as1d(self.mass, n, "mass")
        self.rho = self._as1d(self.rho, n, "rho")
        self.h = self._as1d(self.h, n, "h")
        self.values = {k: self._as1d(v, n, f"values[{k!r}]")
                       for k, v in self.values.items()}
        if self.box is not None:
            self.box = np.asarray(self.box, dtype=np.float64).reshape(3)
        self.periodic = self._as_periodic(self.periodic)

    @staticmethod
    def _as1d(a, n, name):
        a = np.asarray(a, dtype=np.float64).ravel()
        if a.shape[0] != n:
            raise ValueError(f"{name} has length {a.shape[0]}, expected {n}")
        return a

    @staticmethod
    def _as_periodic(p):
        if np.isscalar(p):
            return (bool(p), bool(p), bool(p))
        p = tuple(bool(x) for x in p)
        if len(p) != 3:
            raise ValueError("periodic must be a bool or a length-3 sequence")
        return p

    @property
    def n(self):
        return self.pos.shape[0]

    def value_array(self, names=None):
        """Stack the requested value fields into a (N, K) array + the name list.

        `names=None` uses all fields in insertion order.
        """
        if names is None:
            names = list(self.values.keys())
        if not names:
            raise ValueError("no value fields to interpolate")
        cols = [self.values[k] for k in names]
        return np.ascontiguousarray(np.stack(cols, axis=1)), list(names)


def from_tipsy(snap, fields=("rho",), nsmooth=64, periodic=True, box=None):
    """Build `Particles` from a gasoline/tipsy snapshot.

    `snap`     path to the tipsy snapshot (aux files are siblings `<snap>.<label>`).
    `fields`   names of per-particle scalars to interpolate. Recognised:
                 - "rho", "T"/"temp", "u" (col 8 is TEMPERATURE; "u" needs dTuFac),
                   "vx"/"vy"/"vz", "metals", "phi" (tgdata columns);
                 - any vector-aux magnitude, e.g. "Bmag" (from BField aux);
                 - any scalar aux label, e.g. "DivB", "HeI", "HI", "BClean".
    `nsmooth`  neighbour number for the h estimate when no smoothlength aux exists.
    `periodic` periodicity (default True for these test boxes).
    `box`      period (3,). If None, read from the run `.log` BOX PARAMETERS when
               available, else inferred from the particle extent.

    Smoothing length: uses a `smoothlength` aux if present, else reconstructs
    h = getsmooth2(mass, rho, nsmooth) (matches the divBerr render convention).
    """
    import readtipsy as tip
    import IC_smoothlength as smth
    from IC_analysis_general import try_aux, read_vector_aux

    tgdata, _td, _ts, hdr, _time = tip.readtipsy(snap)
    ngas = int(hdr[2])
    tg = np.asarray(tgdata, dtype=np.float64)[:ngas]

    mass = tg[:, 0]
    pos = tg[:, 1:4]
    rho = tg[:, 7]

    h = try_aux(snap, "smoothlength", ngas)
    if h is None:
        h = smth.getsmooth2(mass, rho, nsmooth)
    h = np.asarray(h, dtype=np.float64).ravel()[:ngas]

    if box is None:
        box = _read_box(snap)
    if box is None:
        lo = pos.min(axis=0)
        hi = pos.max(axis=0)
        box = hi - lo

    col = {"rho": 7, "t": 8, "temp": 8, "vx": 4, "vy": 5, "vz": 6,
           "metals": 10, "phi": 11}
    bcomp = {"bx": 0, "by": 1, "bz": 2}
    values = {}
    for name in fields:
        key = name.lower()
        if key in col:
            values[name] = tg[:, col[key]]
        elif key in ("bmag", "|b|"):
            B = read_vector_aux(snap, "BField", ngas)
            values[name] = (np.zeros(ngas) if B is None
                            else np.sqrt(np.sum(np.asarray(B, float) ** 2, axis=1)))
        elif key in bcomp:                          # magnetic-field components
            B = read_vector_aux(snap, "BField", ngas)
            values[name] = (np.zeros(ngas) if B is None
                            else np.asarray(B, float)[:, bcomp[key]])
        elif key.endswith(("_x", "_y", "_z")):      # any vector aux component
            v = read_vector_aux(snap, name[:-2], ngas)
            c = {"x": 0, "y": 1, "z": 2}[key[-1]]
            values[name] = (np.zeros(ngas) if v is None
                            else np.asarray(v, float)[:, c])
        else:  # treat as a scalar aux label (DivB, HeI, HI, BClean, ...)
            a = try_aux(snap, name, ngas)
            values[name] = np.zeros(ngas) if a is None else np.asarray(a, float)

    return Particles(pos=pos, mass=mass, rho=rho, h=h, values=values,
                     box=box, periodic=periodic)


def _read_box(snap):
    """Periodic box period from the run `.log` BOX PARAMETERS, or None."""
    try:
        from IC_analysis_framework import read_box_periods
    except Exception:
        return None
    import os
    run = os.path.dirname(snap) or "."
    periods = read_box_periods(run)          # {1: Lx, 2: Ly, 3: Lz} or None
    if periods is None:
        return None
    return np.array([periods[1], periods[2], periods[3]], dtype=np.float64)
