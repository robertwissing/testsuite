"""Read a gasoline .param and build the box + target density it encodes.

The C gdicgen reads the density profile it relaxes toward straight from the
.param (dICdensprofile / dICdensdir / dICdensR / ...). This module ports the
same reading so the Python backend (cli.py) and the benchmark get an identical
target without a separate density table - profiles 1/2/3 (uniform / step /
tanh). Profile >= 4 (disk / accdisk) is not ported and needs a density table
(TableDensity.from_xdr); see PLAN.md "Future tasks".
"""
import numpy as np

from .density import DensityField


def read_param(path):
    """Parse a gasoline .param into a {key: value-string} dict."""
    d = {}
    with open(path) as f:
        for line in f:
            line = line.split('#')[0]
            if '=' not in line:
                continue
            k, v = (s.strip() for s in line.split('=', 1))
            d[k] = v
    return d


def box_from_param(p):
    """Periodic box (3,) from the .param, or None for an open domain."""
    g = lambda k, dv=None: float(p[k]) if k in p else dv  # noqa: E731
    if int(float(p.get('bPeriodic', 1))) == 0:
        return None
    dp = g('dPeriod', 1.0)
    return np.array([g('dxPeriod', dp), g('dyPeriod', dp), g('dzPeriod', dp)])


def target_rho0(p):
    """Callable rho0(pos) -> (N,) : a faithful port of pkdUpdateDensity
    profiles 1/2/3 (uniform / radial step / tanh) from a .param. For uniform
    (profile 1) the absolute scale is irrelevant (only the ratio rho/rho0
    drives the relaxation); for 2/3 the inner/outer set the contrast, matching
    the C generator."""
    prof = int(float(p.get('dICdensprofile', 1)))
    d = int(float(p.get('dICdensdir', 1)))
    R = float(p.get('dICdensR', 0.0))
    Rs = float(p.get('dICdensRsmooth', 1.0))
    inner = float(p.get('dICdensinner', 1.0))
    outer = float(p.get('dICdensouter', 1.0))

    def rho0(pos):
        if d < 4:
            r = np.abs(pos[:, d - 1])
        elif d == 4:
            r = np.sqrt((pos ** 2).sum(1))
        elif d == 5:
            r = np.sqrt(pos[:, 0] ** 2 + pos[:, 1] ** 2)
        else:
            r = pos[:, d - 6]
        if prof >= 4:
            raise NotImplementedError(
                f'density profile {prof} (disk/accdisk) is not ported - pass a '
                'density table (--table densitytable_xdr) or use the C backend; '
                'see PLAN.md "Future tasks".')
        if prof == 1:
            return np.full(len(pos), outer)
        if prof == 2:
            return np.where(r > R, outer, inner)
        # profile 3 (tanh): radial for r-directions, double-tanh slab for x/y/z
        if d > 3:
            return outer + (inner - outer) / (1.0 + np.exp(-(R - r) / Rs))
        x = pos[:, d - 1]
        return outer + (inner - outer) * 0.5 * (
            np.tanh((x + R) / Rs) - np.tanh((x - R) / Rs))

    return rho0


class ParamDensity(DensityField):
    """DensityField backed by a .param's dICdens* profile (profiles 1/2/3)."""

    def __init__(self, p):
        self._fn = target_rho0(p)

    def rho0(self, pos):
        return self._fn(pos)


def density_from_param(p):
    """DensityField for the .param's density profile - the Python backend's
    automatic target, matching what the C gdicgen reads from the same file."""
    return ParamDensity(p)
