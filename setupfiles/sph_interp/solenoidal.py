"""Make an interpolated cell-centred B divergence-free for constrained transport.

`project_solenoidal` performs the MAC (staggered Marker-and-Cell) Hodge
projection on a UniformGrid result: it averages the cell-centred B onto the cell
faces, removes the divergent part by solving a discrete Poisson equation, and
returns the six RAMSES-style FACE-centred components (Bx_l/Bx_r/By_l/By_r/
Bz_l/Bz_r). These satisfy the *discrete* constrained-transport divergence

    (Bx_r - Bx_l)/dx + (By_r - By_l)/dy + (Bz_r - Bz_l)/dz = 0

to machine precision per cell — guaranteed algebraically (not by interpolation
quality) because the same staggered grad/div stencils compose into the Poisson
operator used in the solve: div(B - grad phi) = div B - L phi = 0 identically.

Periodic only (FFT). The projection changes the field by exactly its divergent
(monopole) component — typically ~1% for a reasonably div-cleaned SPMHD snapshot,
i.e. it removes the unphysical part and leaves the physical field essentially
intact. See the constrained-transport note in export.py.
"""

import numpy as np

from .targets import UniformGrid


def project_solenoidal(result, fields=("Bx", "By", "Bz"), periodic=True,
                       add_to_result=True):
    """Project a cell-centred B (on a UniformGrid result) to a div-free face B.

    result : a `GridResult` whose target is a `UniformGrid`, carrying the three
             cell-centred magnetic components named by `fields`.
    fields : the (Bx, By, Bz) field names in `result.data`.
    periodic : FFT projection assumes a periodic box (only mode implemented).
    add_to_result : also insert the outputs into `result.data` (so a subsequent
             `to_ramses(..., hydro_fields=RAMSES_MHD_VARS)` finds them).

    Returns dict with the six face components `Bx_l,Bx_r,By_l,By_r,Bz_l,Bz_r`
    (each (nx,ny,nz)); the cleaned cell-centred field `Bx_c,By_c,Bz_c` (face
    averages, for plotting/analysis); and `divB` (residual cell divergence of the
    cleaned face field — ~1e-14).
    """
    grid = getattr(result, "target", None)
    if not isinstance(grid, UniformGrid):
        raise TypeError("project_solenoidal needs a UniformGrid GridResult "
                        f"(got {type(grid).__name__})")
    if not periodic:
        raise NotImplementedError(
            "only the periodic FFT projection is implemented; a non-periodic box "
            "needs Neumann BCs (DCT/DST) or a multigrid solve.")

    bx = np.ascontiguousarray(result.data[fields[0]], dtype=np.float64)
    by = np.ascontiguousarray(result.data[fields[1]], dtype=np.float64)
    bz = np.ascontiguousarray(result.data[fields[2]], dtype=np.float64)
    nx, ny, nz = grid.npx
    dx, dy, dz = grid.pixwidth

    # cell-centred -> face-centred (MAC): face i+1/2 = mean of cells i, i+1
    Bxf = 0.5 * (bx + np.roll(bx, -1, axis=0))
    Byf = 0.5 * (by + np.roll(by, -1, axis=1))
    Bzf = 0.5 * (bz + np.roll(bz, -1, axis=2))

    def celldiv(ax, ay, az):
        return ((ax - np.roll(ax, 1, axis=0)) / dx
                + (ay - np.roll(ay, 1, axis=1)) / dy
                + (az - np.roll(az, 1, axis=2)) / dz)

    D = celldiv(Bxf, Byf, Bzf)

    # discrete Poisson  L phi = D  via FFT, with the 7-point-stencil symbol
    # lambda = -sum (2/dp)^2 sin^2(pi m / n)  (NOT -k^2: must match the stencil)
    sx = (2.0 / dx) * np.sin(np.pi * (np.fft.fftfreq(nx) * nx) / nx)
    sy = (2.0 / dy) * np.sin(np.pi * (np.fft.fftfreq(ny) * ny) / ny)
    sz = (2.0 / dz) * np.sin(np.pi * (np.fft.fftfreq(nz) * nz) / nz)
    lam = -(sx[:, None, None]**2 + sy[None, :, None]**2 + sz[None, None, :]**2)
    lam[0, 0, 0] = 1.0                       # zero-mode of phi := 0
    phi = np.fft.ifftn(np.fft.fftn(D) / lam).real

    # correct the faces: subtract the discrete face-gradient of phi
    Bxf = Bxf - (np.roll(phi, -1, axis=0) - phi) / dx
    Byf = Byf - (np.roll(phi, -1, axis=1) - phi) / dy
    Bzf = Bzf - (np.roll(phi, -1, axis=2) - phi) / dz

    out = {
        # RAMSES face components: right face = Bxf[i]; left face = Bxf[i-1]
        "Bx_l": np.roll(Bxf, 1, axis=0), "Bx_r": Bxf,
        "By_l": np.roll(Byf, 1, axis=1), "By_r": Byf,
        "Bz_l": np.roll(Bzf, 1, axis=2), "Bz_r": Bzf,
        # cleaned cell-centred field (for analysis/plots)
        "Bx_c": 0.5 * (Bxf + np.roll(Bxf, 1, axis=0)),
        "By_c": 0.5 * (Byf + np.roll(Byf, 1, axis=1)),
        "Bz_c": 0.5 * (Bzf + np.roll(Bzf, 1, axis=2)),
        "divB": celldiv(Bxf, Byf, Bzf),
    }
    if add_to_result:
        result.data.update(out)
    return out
