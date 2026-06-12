"""Top-level dispatch: interpolate particles onto a target with a chosen method.

    result = interpolate(particles, UniformGrid.centered(box, 128), method="sph")
    rho_grid = result.data["rho"]

The (method, target) pair selects the engine. Implemented now:
(sph, UniformGrid) and (petkova, UniformGrid). Other targets raise until ported.
"""

from dataclasses import dataclass
from typing import Dict
import numpy as np

from .targets import UniformGrid, Projection2D, AMRGrid, VoronoiGrid
from .sph_method import (interpolate_uniform_sph, interpolate_projection2d,
                         interpolate_amr_sph, interpolate_voronoi_sph)
from .petkova_method import (interpolate_uniform_petkova,
                             interpolate_projection2d_petkova,
                             interpolate_amr_petkova, interpolate_voronoi_petkova)


@dataclass
class GridResult:
    """Interpolated grids plus the geometry needed to place/plot them.

    data    : dict[name -> ndarray] interpolated field grids (target-shaped).
    norm    : summed kernel weight per cell (SPH method) or None (Petkova).
    target  : the target geometry used.
    method  : "sph" or "petkova".
    """
    data: Dict[str, np.ndarray]
    norm: np.ndarray
    target: object
    method: str

    def edges(self):
        return self.target.edges()

    def centres(self):
        return self.target.centres()


def interpolate(particles, target, method="sph", values=None, **kwargs):
    """Interpolate `particles` onto `target` using `method`.

    method : "sph" (kernel deposit, fast/approximate) or
             "petkova" (exact kernel-volume integral, mass-conserving).
    values : list of field names to interpolate (default: all particles.values).
    kwargs : forwarded to the engine (e.g. hmin, normalise for sph).
    """
    method = method.lower()
    if isinstance(target, UniformGrid):
        if method == "sph":
            data, norm = interpolate_uniform_sph(particles, target,
                                                 values=values, **kwargs)
            return GridResult(data=data, norm=norm, target=target, method="sph")
        elif method == "petkova":
            data = interpolate_uniform_petkova(particles, target,
                                               values=values, **kwargs)
            return GridResult(data=data, norm=None, target=target,
                              method="petkova")
        raise ValueError(f"unknown method {method!r} (use 'sph' or 'petkova')")
    if isinstance(target, Projection2D):
        if method == "sph":
            data = interpolate_projection2d(particles, target, values=values,
                                            **kwargs)
            return GridResult(data=data, norm=None, target=target, method="sph")
        elif method == "petkova":
            data = interpolate_projection2d_petkova(particles, target,
                                                    values=values, **kwargs)
            return GridResult(data=data, norm=None, target=target,
                              method="petkova")
        raise ValueError(f"unknown method {method!r} (use 'sph' or 'petkova')")
    if isinstance(target, AMRGrid):
        if method == "sph":
            data, norm = interpolate_amr_sph(particles, target, values=values,
                                             **kwargs)
            return GridResult(data=data, norm=norm, target=target, method="sph")
        elif method == "petkova":
            data = interpolate_amr_petkova(particles, target, values=values,
                                           **kwargs)
            return GridResult(data=data, norm=None, target=target,
                              method="petkova")
        raise ValueError(f"unknown method {method!r} (use 'sph' or 'petkova')")
    if isinstance(target, VoronoiGrid):
        if method == "sph":
            data, norm = interpolate_voronoi_sph(particles, target,
                                                 values=values, **kwargs)
            return GridResult(data=data, norm=norm, target=target, method="sph")
        elif method == "petkova":
            data = interpolate_voronoi_petkova(particles, target,
                                               values=values, **kwargs)
            return GridResult(data=data, norm=None, target=target,
                              method="petkova")
        raise ValueError(f"unknown method {method!r} (use 'sph' or 'petkova')")
    raise NotImplementedError(
        f"target {type(target).__name__} not supported yet "
        f"(UniformGrid, Projection2D, AMRGrid and VoronoiGrid are implemented)")
