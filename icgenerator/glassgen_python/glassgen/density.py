"""Target density fields: table lookup or arbitrary Python callable.

Coordinate mappings follow the C direction codes in pkdUpdateDensity
(pkd.c): 1=x, 2=y, 3=z, 4=spherical r, 5=cylindrical r, 6=signed x.
A 1D table/callable is evaluated at the mapped coordinate; a 2D table at a
pair of mapped coordinates; a 3D table at (x, y, z) directly.

TableDensity reads/writes the densitytable_xdr format used by the C
icgenerator and setupfiles/IC_denstable.py: a 3-line ascii header
(dim, n, R) followed by XDR float arrays (each prefixed by a u32 length):
one length-n array per axis, then the flattened (C-order) density grid.
The reader is implemented with numpy big-endian dtypes rather than xdrlib
(removed in Python 3.13).
"""
import numpy as np
from numba import njit, prange

MAPPINGS = ('x', 'y', 'z', 'r_sph', 'r_cyl', 'signed_x')


def map_coord(pos, mapping):
    if mapping == 'x':
        return np.abs(pos[:, 0])
    if mapping == 'y':
        return np.abs(pos[:, 1])
    if mapping == 'z':
        return np.abs(pos[:, 2])
    if mapping == 'r_sph':
        return np.sqrt(pos[:, 0]**2 + pos[:, 1]**2 + pos[:, 2]**2)
    if mapping == 'r_cyl':
        return np.sqrt(pos[:, 0]**2 + pos[:, 1]**2)
    if mapping == 'signed_x':
        return pos[:, 0].copy()
    raise ValueError(f"unknown mapping '{mapping}', expected one of {MAPPINGS}")


@njit(cache=True, parallel=True)
def _interp1d(axis, values, x, out):
    n = axis.shape[0]
    x0 = axis[0]
    dx = axis[1] - axis[0]
    for i in prange(x.shape[0]):
        t = (x[i] - x0) / dx
        if t <= 0.0:
            out[i] = values[0]
        elif t >= n - 1:
            out[i] = values[n - 1]
        else:
            j = int(t)
            f = t - j
            out[i] = values[j] * (1.0 - f) + values[j + 1] * f


@njit(cache=True, parallel=True)
def _interp2d(ax0, ax1, values, x, y, out):
    n0 = ax0.shape[0]
    n1 = ax1.shape[0]
    for i in prange(x.shape[0]):
        t0 = (x[i] - ax0[0]) / (ax0[1] - ax0[0])
        t1 = (y[i] - ax1[0]) / (ax1[1] - ax1[0])
        if t0 < 0.0:
            t0 = 0.0
        elif t0 > n0 - 1:
            t0 = n0 - 1.0
        if t1 < 0.0:
            t1 = 0.0
        elif t1 > n1 - 1:
            t1 = n1 - 1.0
        j0 = min(int(t0), n0 - 2)
        j1 = min(int(t1), n1 - 2)
        f0 = t0 - j0
        f1 = t1 - j1
        out[i] = (values[j0, j1] * (1 - f0) * (1 - f1)
                  + values[j0, j1 + 1] * (1 - f0) * f1
                  + values[j0 + 1, j1] * f0 * (1 - f1)
                  + values[j0 + 1, j1 + 1] * f0 * f1)


@njit(cache=True, parallel=True)
def _interp3d(ax0, ax1, ax2, values, pos, out):
    n0, n1, n2 = ax0.shape[0], ax1.shape[0], ax2.shape[0]
    for i in prange(pos.shape[0]):
        t0 = (pos[i, 0] - ax0[0]) / (ax0[1] - ax0[0])
        t1 = (pos[i, 1] - ax1[0]) / (ax1[1] - ax1[0])
        t2 = (pos[i, 2] - ax2[0]) / (ax2[1] - ax2[0])
        if t0 < 0.0:
            t0 = 0.0
        elif t0 > n0 - 1:
            t0 = n0 - 1.0
        if t1 < 0.0:
            t1 = 0.0
        elif t1 > n1 - 1:
            t1 = n1 - 1.0
        if t2 < 0.0:
            t2 = 0.0
        elif t2 > n2 - 1:
            t2 = n2 - 1.0
        j0 = min(int(t0), n0 - 2)
        j1 = min(int(t1), n1 - 2)
        j2 = min(int(t2), n2 - 2)
        f0, f1, f2 = t0 - j0, t1 - j1, t2 - j2
        c00 = values[j0, j1, j2] * (1 - f2) + values[j0, j1, j2 + 1] * f2
        c01 = values[j0, j1 + 1, j2] * (1 - f2) + values[j0, j1 + 1, j2 + 1] * f2
        c10 = values[j0 + 1, j1, j2] * (1 - f2) + values[j0 + 1, j1, j2 + 1] * f2
        c11 = values[j0 + 1, j1 + 1, j2] * (1 - f2) + values[j0 + 1, j1 + 1, j2 + 1] * f2
        c0 = c00 * (1 - f1) + c01 * f1
        c1 = c10 * (1 - f1) + c11 * f1
        out[i] = c0 * (1 - f0) + c1 * f0


class DensityField:
    """Interface: rho0(pos) maps (N,3) positions to (N,) target densities."""

    def rho0(self, pos):
        raise NotImplementedError


class TableDensity(DensityField):
    def __init__(self, axes, values, mapping='r_sph'):
        values = np.ascontiguousarray(values, dtype=np.float64)
        self.axes = tuple(np.ascontiguousarray(a, dtype=np.float64)
                          for a in (axes if isinstance(axes, (tuple, list)) else (axes,)))
        self.values = values
        self.dim = len(self.axes)
        if values.ndim != self.dim:
            raise ValueError(f"{self.dim} axes but values has ndim {values.ndim}")
        if self.dim == 1:
            self.mapping = mapping
        elif self.dim == 2:
            if not (isinstance(mapping, (tuple, list)) and len(mapping) == 2):
                raise ValueError("2D table needs a pair of mappings, e.g. ('r_cyl','z')")
            self.mapping = tuple(mapping)
        else:
            self.mapping = None  # 3D tables are indexed by (x, y, z)

    @classmethod
    def from_xdr(cls, filename, mapping='r_sph'):
        with open(filename, 'rb') as f:
            dim = int(f.readline())
            n = int(f.readline())
            float(f.readline())  # R, domain size; axes carry the actual grid
            data = f.read()

        be_u4 = np.dtype('>u4')
        be_f4 = np.dtype('>f4')
        off = 0

        def unpack_array():
            nonlocal off
            count = int(np.frombuffer(data, be_u4, count=1, offset=off)[0])
            off += 4
            arr = np.frombuffer(data, be_f4, count=count, offset=off)
            off += 4 * count
            return arr.astype(np.float64)

        axes = tuple(unpack_array() for _ in range(dim))
        for ax in axes:
            if len(ax) != n:
                raise ValueError(f"axis length {len(ax)} != header n {n}")
        values = unpack_array().reshape((n,) * dim)
        return cls(axes, values, mapping=mapping)

    def rho0(self, pos):
        out = np.empty(len(pos))
        if self.dim == 1:
            _interp1d(self.axes[0], self.values, map_coord(pos, self.mapping), out)
        elif self.dim == 2:
            _interp2d(self.axes[0], self.axes[1], self.values,
                      map_coord(pos, self.mapping[0]),
                      map_coord(pos, self.mapping[1]), out)
        else:
            _interp3d(self.axes[0], self.axes[1], self.axes[2], self.values,
                      np.ascontiguousarray(pos), out)
        return out


class CallableDensity(DensityField):
    """Wrap a callable. With a mapping, fn receives the (N,) mapped
    coordinate; with mapping=None, fn receives the (N,3) positions."""

    def __init__(self, fn, mapping=None):
        self.fn = fn
        self.mapping = mapping

    def rho0(self, pos):
        if self.mapping is None:
            rho = self.fn(pos)
        else:
            rho = self.fn(map_coord(pos, self.mapping))
        return np.asarray(rho, dtype=np.float64)
