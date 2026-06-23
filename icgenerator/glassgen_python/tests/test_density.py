import struct

import numpy as np
import pytest

from glassgen.density import CallableDensity, TableDensity, map_coord


def write_xdr_table(path, dim, n, R, axes, values):
    """Reference writer matching setupfiles/IC_denstable.py (xdrlib format:
    each array prefixed by a big-endian u32 length, big-endian f32 data)."""
    with open(path, 'wb') as f:
        f.write(f"{dim}\n{n}\n{R:.15g}\n".encode('ascii'))
        for ax in axes:
            f.write(struct.pack('>I', len(ax)))
            f.write(np.asarray(ax, dtype='>f4').tobytes())
        flat = np.asarray(values, dtype='>f4').ravel()
        f.write(struct.pack('>I', len(flat)))
        f.write(flat.tobytes())


def test_xdr_roundtrip_1d(tmp_path):
    n = 100
    x = np.linspace(0, 2.0, n)
    rho = 1.0 + np.exp(-x)
    path = tmp_path / 'densitytable_xdr'
    write_xdr_table(path, 1, n, 2.0, (x,), rho)
    tab = TableDensity.from_xdr(str(path), mapping='r_sph')
    assert tab.dim == 1
    assert np.allclose(tab.axes[0], x, atol=1e-6)
    assert np.allclose(tab.values, rho, atol=1e-6)


def test_xdr_roundtrip_2d(tmp_path):
    n = 16
    x = np.linspace(-1, 1, n)
    vals = np.add.outer(x ** 2, x)  # (n, n), C-order
    path = tmp_path / 'tab2d'
    write_xdr_table(path, 2, n, 1.0, (x, x), vals)
    tab = TableDensity.from_xdr(str(path), mapping=('x', 'y'))
    assert tab.dim == 2
    assert np.allclose(tab.values, vals, atol=1e-5)


def test_linear_table_interp_is_exact():
    x = np.linspace(0, 1, 11)
    tab = TableDensity((x,), 2.0 + 3.0 * x, mapping='x')
    pos = np.zeros((50, 3))
    pos[:, 0] = np.linspace(0.01, 0.99, 50)
    assert np.allclose(tab.rho0(pos), 2.0 + 3.0 * pos[:, 0])


def test_interp_clamps_at_edges():
    x = np.linspace(0, 1, 11)
    tab = TableDensity((x,), 2.0 + 3.0 * x, mapping='x')
    pos = np.array([[5.0, 0, 0], [-5.0, 0, 0]])
    assert np.allclose(tab.rho0(pos), [5.0, 5.0])  # |x| mapping + clamp


def test_mappings():
    pos = np.array([[3.0, -4.0, 12.0]])
    assert map_coord(pos, 'x')[0] == 3.0
    assert map_coord(pos, 'y')[0] == 4.0
    assert map_coord(pos, 'z')[0] == 12.0
    assert np.isclose(map_coord(pos, 'r_sph')[0], 13.0)
    assert np.isclose(map_coord(pos, 'r_cyl')[0], 5.0)
    pos2 = np.array([[-3.0, 0.0, 0.0]])
    assert map_coord(pos2, 'signed_x')[0] == -3.0
    with pytest.raises(ValueError):
        map_coord(pos, 'bogus')


def test_callable_vs_table_agree():
    fn = lambda r: 1.0 + np.exp(-r * r)
    r = np.linspace(0, 3, 3000)
    tab = TableDensity((r,), fn(r), mapping='r_sph')
    cal = CallableDensity(fn, mapping='r_sph')
    rng = np.random.default_rng(5)
    pos = rng.uniform(-1, 1, (200, 3))
    assert np.allclose(tab.rho0(pos), cal.rho0(pos), rtol=1e-5)


def test_3d_table():
    x = np.linspace(-1, 1, 21)
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    vals = 1.0 + X + 2 * Y + 3 * Z  # trilinear -> interp exact
    tab = TableDensity((x, x, x), vals)
    rng = np.random.default_rng(6)
    pos = rng.uniform(-0.9, 0.9, (100, 3))
    expect = 1.0 + pos[:, 0] + 2 * pos[:, 1] + 3 * pos[:, 2]
    assert np.allclose(tab.rho0(pos), expect)
