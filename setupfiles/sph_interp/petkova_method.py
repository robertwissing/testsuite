"""Petkova (2018) exact interpolation: integrate the kernel over each cell.

Mass-conserving alternative to the regular SPH deposit: instead of evaluating the
kernel at the cell centre, integrate the (cubic-spline) kernel over the cell
volume analytically (per particle, per overlapping cell), so the deposited mass
equals the particle mass to high precision.

`_vertex_integral` / `_compute_logs` are cleaned ports of the math core from
`interpolate_petkova.py` (guards instead of try/except, no origin hardwiring, no
`/5.0`). `_cube_integral` applies the per-face/edge construction of
`exact_interpolate` to an axis-aligned cube with the PARTICLE at the origin,
returning the fraction of the (unit-normalized) kernel mass contained in the cube.

Currently implements the UniformGrid target. NOTE: the exact integral is
expensive (per particle × overlapping cell × 6 faces × 4 edges); intended for
accuracy/validation. Far cells (entirely beyond 2h) are skipped.

VERIFIED mass-conserving: a single particle deposits its full mass for any npx;
on a real orzag snapshot the summed grid mass matches the particle mass to 0.00%
(vs ~0.5% SPH diffusion). Two things are essential and easy to get wrong:
  * the cube vertex numbering + face windings MUST match pyvoro's convention
    (see `_CUBE_FACES`) — `vertex_integral`'s sign logic depends on it;
  * the per-particle vertex scratch arrays MUST be allocated inside the prange
    loop (thread-private), else parallel threads corrupt each other.
"""

import math
import numpy as np
import numba
from numba import njit, prange

from .kernels import pcum, PCUM_TABLE, gcum3D
from .targets import (UniformGrid, Projection2D, AMRGrid, VoronoiGrid,
                      voronoi_candidates)

# Gauss-Legendre nodes/weights on [0,1] for the per-edge angular quadrature of
# the exact 2D column integral (the integrand is smooth and bounded, so a modest
# order is plenty / effectively exact).
_GLX, _GLW = np.polynomial.legendre.leggauss(32)
_GLX = 0.5 * (_GLX + 1.0)          # map [-1,1] -> [0,1]
_GLW = 0.5 * _GLW


@njit(cache=True)
def _compute_logs(u, tol=1e-12):
    if u >= 1.0 - tol:
        return 1000.0
    elif u <= -1.0 + tol:
        return -1000.0
    else:
        return 2.0 * math.atanh(u)


@njit(cache=True)
def _clamp_acos(x):
    if x < -1.0:
        x = -1.0
    elif x > 1.0:
        x = 1.0
    return math.acos(x)


@njit(cache=True)
def _iblock(phi, a, a2, r0h, r0h2, r0h3, r0h_3):
    """The I-integrals I0..I_5, I1 at angle phi (shared by all regions)."""
    cosp = math.cos(phi)
    cosp2 = cosp * cosp
    mu = cosp / a / math.sqrt(1.0 + cosp2 / a2)
    tanp = math.tan(phi)
    I0 = phi
    I_2 = phi + a2 * tanp
    I_4 = phi + 2.0 * a2 * tanp + (a2 * a2 / 3.0) * tanp * (2.0 + 1.0 / cosp2)
    arg = 1.0 - mu * mu
    u = math.sin(phi) * math.sqrt(arg) if arg >= 0.0 else 0.0
    logs = _compute_logs(u)
    I1 = math.atan(u / a)
    I_1 = (a / 2.0) * logs + I1
    om = 1.0 - u * u
    if om != 0.0:
        I_3 = I_1 + a * (1.0 + a2) / 4.0 * (2.0 * u / om + logs)
        I_5 = I_3 + a * (1.0 + a2) ** 2 / 16.0 * ((10.0 * u - 6.0 * u ** 3) / (om * om) + 3.0 * logs)
    else:
        I_3 = I_1 + a * (1.0 + a2) / 4.0 * logs
        I_5 = I_3 + a * (1.0 + a2) ** 2 / 16.0 * 3.0 * logs
    return I0, I1, I_1, I_2, I_3, I_4, I_5


@njit(cache=True)
def _vertex_integral(phi, r0, R_0, h):
    """Exact kernel wall integral (Petkova 2018), cleaned/guarded. r0,R_0,h>0."""
    if r0 == 0.0 or R_0 == 0.0 or phi == 0.0:
        return 0.0

    r03 = r0 ** 3
    r0h = r0 / h
    r0h2 = r0h ** 2
    r0h3 = r0h ** 3
    r0h_2 = (h / r0) ** 2
    r0h_3 = (h / r0) ** 3

    B1, B2, B3 = 0.0, 0.0, 0.0
    if r0 >= 2.0 * h:
        B3 = (h ** 3) / 4.0
    elif r0 > h:
        B3 = (r03 / 4.0) * (-4.0/3.0 + r0h - 0.3*r0h2 + (1.0/30.0)*r0h3 - (1.0/15.0)*r0h_3 + (8.0/5.0)*r0h_2)
        B2 = (r03 / 4.0) * (-4.0/3.0 + r0h - 0.3*r0h2 + (1.0/30.0)*r0h3 - (1.0/15.0)*r0h_3)
    else:
        B3 = (r03 / 4.0) * (-2.0/3.0 + 0.3*r0h2 - 0.1*r0h3 + (7.0/5.0)*r0h_2)
        B2 = (r03 / 4.0) * (-2.0/3.0 + 0.3*r0h2 - 0.1*r0h3 - (1.0/5.0)*r0h_2)
        B1 = (r03 / 4.0) * (-2.0/3.0 + 0.3*r0h2 - 0.1*r0h3)

    a = R_0 / r0
    a2 = a * a
    linedist = math.sqrt(r0 * r0 + R_0 * R_0)
    R = R_0 / math.cos(phi)
    r = math.sqrt(r0 * r0 + R * R)

    D2 = 0.0
    D3 = 0.0

    if linedist < 1.0 * h:
        arg = h * h - r0 * r0
        phi1 = _clamp_acos(R_0 / math.sqrt(arg)) if arg > 0.0 else 0.0
        I0, I1, I_1, I_2, I_3, I_4, I_5 = _iblock(phi1, a, a2, r0h, r0h2, r0h3, r0h_3)
        D2 = (-1.0/6 * I_2 + 0.25 * r0h * I_3 - 0.15 * r0h2 * I_4 +
              (1.0/30) * r0h3 * I_5 - (1.0/60) * r0h_3 * I1 +
              (B1 - B2) / r03 * I0)
        arg = 4.0 * h * h - r0 * r0
        phi2 = _clamp_acos(R_0 / math.sqrt(arg)) if arg > 0.0 else 0.0
        I0, I1, I_1, I_2, I_3, I_4, I_5 = _iblock(phi2, a, a2, r0h, r0h2, r0h3, r0h_3)
        D3 = (1.0/3 * I_2 - 0.25 * r0h * I_3 + 0.075 * r0h2 * I_4 -
              (1.0/120)*r0h3 * I_5 + (4.0/15)*r0h_3 * I1 +
              (B2 - B3)/r03 * I0 + D2)
    elif linedist < 2.0 * h:
        arg = 4.0 * h * h - r0 * r0
        phi2 = _clamp_acos(R_0 / math.sqrt(arg)) if arg > 0.0 else 0.0
        I0, I1, I_1, I_2, I_3, I_4, I_5 = _iblock(phi2, a, a2, r0h, r0h2, r0h3, r0h_3)
        D3 = (1.0/3 * I_2 - 0.25 * r0h * I_3 + 0.075 * r0h2 * I_4 -
              (1.0/120)*r0h3 * I_5 + (4.0/15)*r0h_3 * I1 +
              (B2 - B3)/r03 * I0 + D2)

    I0, I1, I_1, I_2, I_3, I_4, I_5 = _iblock(phi, a, a2, r0h, r0h2, r0h3, r0h_3)
    if r < 1.0 * h:
        term1 = (1.0/6 * I_2 - 0.075 * r0h2 * I_4 + (1.0/40) * r0h3 * I_5 + (B1 / r03) * I0)
        return (r0h3 / math.pi) * term1
    elif r < 2.0 * h:
        term1 = 0.25 * (4.0/3 * I_2 - r0h * I_3 + 0.3 * r0h2 * I_4 -
                        (1.0/30)*r0h3 * I_5 + (1.0/15)*r0h_3 * I1) + (B2 / r03) * I0 + D2
        return (r0h3 / math.pi) * term1
    else:
        term1 = -0.25 * r0h_3 * I1 + (B3 / r03) * I0 + D3
        return (r0h3 / math.pi) * term1


# Axis-aligned cube topology, MATCHING pyvoro's convention exactly (the working
# exact_interpolate depends on this vertex numbering + face winding for its sign
# logic). Vertex index i is the binary triple (bit0=x, bit1=y, bit2=z):
#   0(-,-,-) 1(+,-,-) 2(-,+,-) 3(+,+,-) 4(-,-,+) 5(+,-,+) 6(-,+,+) 7(+,+,+)
# Face windings are pyvoro's (verified by feeding a single-generator box cell to
# exact_interpolate and getting Mtot=1).
_CUBE_FACES = np.array([
    [1, 3, 2, 0],   # z-lo
    [1, 5, 7, 3],   # x-hi
    [1, 0, 4, 5],   # y-lo
    [2, 6, 4, 0],   # x-lo
    [2, 3, 7, 6],   # y-hi
    [4, 6, 7, 5],   # z-hi
], dtype=np.int64)


@njit(cache=True)
def _cone_tri(ax, ay, az, bx, by, bz, cx, cy, cz, h):
    """Kernel integral over the cone from the PARTICLE (at the origin) through the
    triangle A,B,C (coords RELATIVE to the particle). Always >= 0.

    Building block of the signed-tetrahedralization cell integral:
     int_cell W = sum over faces (front/back sign) * sum over fan-triangles
    _cone_tri. The cone integral is int_{Omega in tri} g(R(Omega)) dOmega with
    g(R) = int_0^R W r^2 dr the radial mass profile and R(Omega) = r0/cos(theta)
    the distance to the triangle's plane (r0 = perp distance particle->plane).

    Computed ANALYTICALLY (exact at any face size): a triangle entirely beyond 2h
    (r0 >= 2h) is SATURATED -> closed-form solid angle (Van Oosterom & Strijbos);
    otherwise the triangle is split into right-triangle wedges about the in-plane
    foot of the particle and each wedge is the exact `_vertex_integral` (the
    support-aware Petkova M-integral, verified == the cone integral to ~1e-16),
    combined per edge by the foot-inside/outside-segment rule and the edge's
    in-plane orientation. This is the per-TRIANGLE wedge decomposition (always
    correct); the OLD bug was the per-POLYGON assembly (cond/sign_M across a whole
    face with the particle outside) — avoided entirely by the signed-tet outer
    loop. Unlike a fixed quadrature it does not under-resolve a huge cube face
    whose kernel-active patch is tiny (the coarse-UniformGrid case).
    """
    # plane normal n = (B-A) x (C-A); area2 = |n| = 2 * triangle area
    e1x = bx - ax; e1y = by - ay; e1z = bz - az
    e2x = cx - ax; e2y = cy - ay; e2z = cz - az
    nx = e1y*e2z - e1z*e2y
    ny = e1z*e2x - e1x*e2z
    nz = e1x*e2y - e1y*e2x
    area2 = math.sqrt(nx*nx + ny*ny + nz*nz)
    if area2 < 1e-300:
        return 0.0
    invn = 1.0 / area2
    nhx = nx*invn; nhy = ny*invn; nhz = nz*invn
    cpl = ax*nhx + ay*nhy + az*nhz       # signed plane offset (particle at origin)
    r0 = abs(cpl)
    if r0 < 1e-300:
        return 0.0
    if r0 >= 2.0 * h:
        # saturated: cone = (1/(4 pi)) * solid angle (Van Oosterom & Strijbos)
        la = math.sqrt(ax*ax + ay*ay + az*az)
        lb = math.sqrt(bx*bx + by*by + bz*bz)
        lc = math.sqrt(cx*cx + cy*cy + cz*cz)
        triple = (ax*(by*cz - bz*cy) + ay*(bz*cx - bx*cz) + az*(bx*cy - by*cx))
        dab = ax*bx + ay*by + az*bz
        dac = ax*cx + ay*cy + az*cz
        dbc = bx*cx + by*cy + bz*cz
        denom = la*lb*lc + dab*lc + dac*lb + dbc*la
        omega = 2.0 * math.atan2(abs(triple), denom)
        return (0.25 / math.pi) * omega
    # foot of the particle on the triangle plane
    fx = cpl*nhx; fy = cpl*nhy; fz = cpl*nhz
    total = 0.0
    for e in range(3):
        if e == 0:
            p1x = ax; p1y = ay; p1z = az; p2x = bx; p2y = by; p2z = bz
        elif e == 1:
            p1x = bx; p1y = by; p1z = bz; p2x = cx; p2y = cy; p2z = cz
        else:
            p1x = cx; p1y = cy; p1z = cz; p2x = ax; p2y = ay; p2z = az
        d1x = p1x - fx; d1y = p1y - fy; d1z = p1z - fz
        d2x = p2x - fx; d2y = p2y - fy; d2z = p2z - fz
        lx = p2x - p1x; ly = p2y - p1y; lz = p2z - p1z
        r23 = math.sqrt(lx*lx + ly*ly + lz*lz)
        if r23 < 1e-300:
            continue
        # R_0 = perpendicular distance from the foot to the edge line = |d1 x L|/|L|
        cx1 = d1y*lz - d1z*ly; cy1 = d1z*lx - d1x*lz; cz1 = d1x*ly - d1y*lx
        R_0 = math.sqrt(cx1*cx1 + cy1*cy1 + cz1*cz1) / r23
        r12 = math.sqrt(d1x*d1x + d1y*d1y + d1z*d1z)
        r13 = math.sqrt(d2x*d2x + d2y*d2y + d2z*d2z)
        phi1 = math.acos(R_0/r12) if (r12 > 0.0 and R_0 < r12) else 0.0
        phi2 = math.acos(R_0/r13) if (r13 > 0.0 and R_0 < r13) else 0.0
        v1 = _vertex_integral(phi1, r0, R_0, h)
        v2 = _vertex_integral(phi2, r0, R_0, h)
        # foot-on-edge-line inside the segment -> wedges add; outside -> subtract
        s1 = math.sqrt(r12*r12 - R_0*R_0) if r12 > R_0 else 0.0
        s2 = math.sqrt(r13*r13 - R_0*R_0) if r13 > R_0 else 0.0
        if s1 >= r23 or s2 >= r23:
            wedge = v1 - v2 if v1 >= v2 else v2 - v1
        else:
            wedge = v1 + v2
        # signed by the edge's in-plane orientation about the foot
        sx = d1y*d2z - d1z*d2y; sy = d1z*d2x - d1x*d2z; sz = d1x*d2y - d1y*d2x
        sgn = 1.0 if (sx*nhx + sy*nhy + sz*nhz) > 0.0 else -1.0
        total += sgn * wedge
    if total < 0.0:
        total = -total
    return total


@njit(cache=True)
def _cube_integral(vx, vy, vz, h, faces):
    """Integral of the unit-normalized kernel over a cube (particle at origin).

    vx,vy,vz : (8,) cube vertex coords RELATIVE to the particle.
    Returns the contained kernel-mass fraction (>= 0). Signed tetrahedralization
    from the particle: one front/back sign per face, fan-triangulate, sum
    `_cone_tri` (see its docstring). Mass-conserving for any cube position/size.
    """
    # interior point = centroid of the 8 vertices (particle-relative), used to
    # orient each face's outward normal.
    cgx = 0.0; cgy = 0.0; cgz = 0.0
    for k in range(8):
        cgx += vx[k]; cgy += vy[k]; cgz += vz[k]
    cgx *= 0.125; cgy *= 0.125; cgz *= 0.125
    total = 0.0
    nfv = faces.shape[1]
    for f in range(faces.shape[0]):
        i0 = faces[f, 0]; i1 = faces[f, 1]; i2 = faces[f, 2]
        v0x = vx[i0]; v0y = vy[i0]; v0z = vz[i0]
        ax = vx[i1]-v0x; ay = vy[i1]-v0y; az = vz[i1]-v0z
        bx = vx[i2]-v0x; by = vy[i2]-v0y; bz = vz[i2]-v0z
        nx = ay*bz - az*by; ny = az*bx - ax*bz; nz = ax*by - ay*bx
        # orient outward (point away from the centroid)
        if nx*(cgx-v0x) + ny*(cgy-v0y) + nz*(cgz-v0z) > 0.0:
            nx = -nx; ny = -ny; nz = -nz
        sgn = 1.0 if (v0x*nx + v0y*ny + v0z*nz) > 0.0 else -1.0
        # fan-triangulate the face polygon from vertex 0
        for k in range(1, nfv - 1):
            ka = faces[f, k]; kb = faces[f, k + 1]
            total += sgn * _cone_tri(
                v0x, v0y, v0z,
                vx[ka], vy[ka], vz[ka],
                vx[kb], vy[kb], vz[kb], h)
    if total < 0.0:
        total = 0.0
    return total


@njit(cache=True)
def _nearest_dist2(xi, yi, zi, x0, x1, y0, y1, z0, z1):
    """Squared distance from a point to an axis-aligned box [x0,x1]x..."""
    dx = 0.0
    if xi < x0:
        dx = x0 - xi
    elif xi > x1:
        dx = xi - x1
    dy = 0.0
    if yi < y0:
        dy = y0 - yi
    elif yi > y1:
        dy = yi - y1
    dz = 0.0
    if zi < z0:
        dz = z0 - zi
    elif zi > z1:
        dz = zi - z1
    return dx*dx + dy*dy + dz*dz


@njit(parallel=True, cache=True)
def _deposit_petkova_uniform(x, y, z, h, mass, volw, vals,
                             nx, ny, nz, lox, loy, loz, pwx, pwy, pwz,
                             perx, pery, perz, faces,
                             mass_local, den_local, num_local):
    npart = x.shape[0]
    K = vals.shape[1]
    for i in prange(npart):
        mi = mass[i]
        hi = h[i]
        if mi <= 0.0 or hi <= 0.0:
            continue
        t = numba.get_thread_id()
        # per-iteration (thread-private) vertex scratch — must NOT be shared
        # across prange iterations or threads race and corrupt each other.
        vx = np.empty(8); vy = np.empty(8); vz = np.empty(8)
        wi = volw[i]
        xi = x[i]; yi = y[i]; zi = z[i]
        kr = 2.0 * hi
        kr2 = kr * kr
        ixmin = int(np.floor((xi - kr - lox) / pwx))
        ixmax = int(np.floor((xi + kr - lox) / pwx)) + 1
        jymin = int(np.floor((yi - kr - loy) / pwy))
        jymax = int(np.floor((yi + kr - loy) / pwy)) + 1
        kzmin = int(np.floor((zi - kr - loz) / pwz))
        kzmax = int(np.floor((zi + kr - loz) / pwz)) + 1
        if not perx:
            ixmin = max(ixmin, 0); ixmax = min(ixmax, nx)
        if not pery:
            jymin = max(jymin, 0); jymax = min(jymax, ny)
        if not perz:
            kzmin = max(kzmin, 0); kzmax = min(kzmax, nz)

        for ix in range(ixmin, ixmax):
            iw = ix % nx if perx else ix
            xe0 = lox + ix * pwx
            xe1 = xe0 + pwx
            for jy in range(jymin, jymax):
                jw = jy % ny if pery else jy
                ye0 = loy + jy * pwy
                ye1 = ye0 + pwy
                for kz in range(kzmin, kzmax):
                    ze0 = loz + kz * pwz
                    ze1 = ze0 + pwz
                    # skip cells entirely outside the kernel sphere
                    if _nearest_dist2(xi, yi, zi, xe0, xe1, ye0, ye1, ze0, ze1) >= kr2:
                        continue
                    kw = kz % nz if perz else kz
                    # cube vertices relative to particle (min-image via raw index)
                    ax0 = xe0 - xi; ax1 = xe1 - xi
                    ay0 = ye0 - yi; ay1 = ye1 - yi
                    az0 = ze0 - zi; az1 = ze1 - zi
                    vx[0]=ax0; vy[0]=ay0; vz[0]=az0
                    vx[1]=ax1; vy[1]=ay0; vz[1]=az0
                    vx[2]=ax0; vy[2]=ay1; vz[2]=az0
                    vx[3]=ax1; vy[3]=ay1; vz[3]=az0
                    vx[4]=ax0; vy[4]=ay0; vz[4]=az1
                    vx[5]=ax1; vy[5]=ay0; vz[5]=az1
                    vx[6]=ax0; vy[6]=ay1; vz[6]=az1
                    vx[7]=ax1; vy[7]=ay1; vz[7]=az1
                    frac = _cube_integral(vx, vy, vz, hi, faces)
                    if frac <= 0.0:
                        continue
                    mass_local[t, iw, jw, kw] += mi * frac
                    den_local[t, iw, jw, kw] += wi * frac
                    for p in range(K):
                        num_local[t, iw, jw, kw, p] += wi * vals[i, p] * frac


def interpolate_uniform_petkova(particles, grid, values=None):
    """Petkova exact interpolation of `particles` onto a `UniformGrid`.

    Returns data: dict[name -> (nx,ny,nz)]. Density-like fields ('rho'/'density')
    are the MASS-CONSERVING estimate sum_j m_j (integral_cell W_j) / V_cell; other
    fields are the kernel-volume-weighted cell average
    sum_j (m_j/rho_j) A_j I / sum_j (m_j/rho_j) I.
    """
    if not isinstance(grid, UniformGrid):
        raise TypeError("interpolate_uniform_petkova requires a UniformGrid target")

    vals, names = particles.value_array(values)
    nx, ny, nz = grid.npx
    lo = grid.lo
    pw = grid.pixwidth
    perx, pery, perz = particles.periodic

    nthreads = numba.get_num_threads()
    K = vals.shape[1]
    mass_local = np.zeros((nthreads, nx, ny, nz), dtype=np.float64)
    den_local = np.zeros((nthreads, nx, ny, nz), dtype=np.float64)
    num_local = np.zeros((nthreads, nx, ny, nz, K), dtype=np.float64)

    pos = particles.pos
    _deposit_petkova_uniform(
        np.ascontiguousarray(pos[:, 0]),
        np.ascontiguousarray(pos[:, 1]),
        np.ascontiguousarray(pos[:, 2]),
        particles.h, particles.mass, particles.mass / particles.rho, vals,
        nx, ny, nz, lo[0], lo[1], lo[2], pw[0], pw[1], pw[2],
        perx, pery, perz, _CUBE_FACES,
        mass_local, den_local, num_local)

    mass = mass_local.sum(axis=0)
    den = den_local.sum(axis=0)
    num = num_local.sum(axis=0)
    cv = grid.cell_volume
    mask = den > 0.0

    data = {}
    for p, name in enumerate(names):
        if name.lower() in ("rho", "density"):
            data[name] = mass / cv
        else:
            g = np.zeros((nx, ny, nz), dtype=np.float64)
            g[mask] = num[..., p][mask] / den[mask]
            data[name] = g
    # always expose the conserved mass grid for convenience
    data.setdefault("_mass", mass)
    return data


# --------------------------------------------------------------------------- #
#  AMR target: exact kernel-volume integral over each octree leaf cube. Same
#  base-cell -> leaves traversal as the SPH-AMR deposit, but each leaf gets the
#  exact contained kernel-mass fraction via `_cube_integral` (mass-conserving on
#  the adaptive cells, like the UniformGrid Petkova path).
# --------------------------------------------------------------------------- #

@njit(parallel=True, cache=True)
def _deposit_petkova_amr(x, y, z, h, mass, volw, vals, cell_lo, cell_size,
                         base_start, base_count, blo0, blo1, blo2,
                         bpw0, bpw1, bpw2, nbx, nby, nbz,
                         perx, pery, perz, boxx, boxy, boxz, faces,
                         mass_local, den_local, num_local):
    npart = x.shape[0]
    K = vals.shape[1]
    for i in prange(npart):
        mi = mass[i]
        hi = h[i]
        if mi <= 0.0 or hi <= 0.0:
            continue
        t = numba.get_thread_id()
        vx = np.empty(8); vy = np.empty(8); vz = np.empty(8)
        wi = volw[i]
        xi = x[i]; yi = y[i]; zi = z[i]
        kr = 2.0 * hi
        kr2 = kr * kr
        ibx0 = int(np.floor((xi - kr - blo0) / bpw0))
        ibx1 = int(np.floor((xi + kr - blo0) / bpw0)) + 1
        iby0 = int(np.floor((yi - kr - blo1) / bpw1))
        iby1 = int(np.floor((yi + kr - blo1) / bpw1)) + 1
        ibz0 = int(np.floor((zi - kr - blo2) / bpw2))
        ibz1 = int(np.floor((zi + kr - blo2) / bpw2)) + 1
        if not perx:
            if ibx0 < 0: ibx0 = 0
            if ibx1 > nbx: ibx1 = nbx
        if not pery:
            if iby0 < 0: iby0 = 0
            if iby1 > nby: iby1 = nby
        if not perz:
            if ibz0 < 0: ibz0 = 0
            if ibz1 > nbz: ibz1 = nbz
        for ibx in range(ibx0, ibx1):
            if perx:
                kx = ibx // nbx; iwx = ibx - kx * nbx; sx = kx * boxx
            else:
                iwx = ibx; sx = 0.0
            for iby in range(iby0, iby1):
                if pery:
                    ky = iby // nby; iwy = iby - ky * nby; sy = ky * boxy
                else:
                    iwy = iby; sy = 0.0
                for ibz in range(ibz0, ibz1):
                    if perz:
                        kz = ibz // nbz; iwz = ibz - kz * nbz; sz = kz * boxz
                    else:
                        iwz = ibz; sz = 0.0
                    bflat = (iwx * nby + iwy) * nbz + iwz
                    s0 = base_start[bflat]
                    n0 = base_count[bflat]
                    for li in range(s0, s0 + n0):
                        xe0 = cell_lo[li, 0] + sx
                        xe1 = xe0 + cell_size[li, 0]
                        ye0 = cell_lo[li, 1] + sy
                        ye1 = ye0 + cell_size[li, 1]
                        ze0 = cell_lo[li, 2] + sz
                        ze1 = ze0 + cell_size[li, 2]
                        if _nearest_dist2(xi, yi, zi, xe0, xe1, ye0, ye1,
                                          ze0, ze1) >= kr2:
                            continue
                        ax0 = xe0 - xi; ax1 = xe1 - xi
                        ay0 = ye0 - yi; ay1 = ye1 - yi
                        az0 = ze0 - zi; az1 = ze1 - zi
                        vx[0]=ax0; vy[0]=ay0; vz[0]=az0
                        vx[1]=ax1; vy[1]=ay0; vz[1]=az0
                        vx[2]=ax0; vy[2]=ay1; vz[2]=az0
                        vx[3]=ax1; vy[3]=ay1; vz[3]=az0
                        vx[4]=ax0; vy[4]=ay0; vz[4]=az1
                        vx[5]=ax1; vy[5]=ay0; vz[5]=az1
                        vx[6]=ax0; vy[6]=ay1; vz[6]=az1
                        vx[7]=ax1; vy[7]=ay1; vz[7]=az1
                        frac = _cube_integral(vx, vy, vz, hi, faces)
                        if frac <= 0.0:
                            continue
                        mass_local[t, li] += mi * frac
                        den_local[t, li] += wi * frac
                        for p in range(K):
                            num_local[t, li, p] += wi * vals[i, p] * frac


def interpolate_amr_petkova(particles, grid, values=None):
    """Petkova exact interpolation of `particles` onto an `AMRGrid`.

    Returns data: dict[name -> (Ncell,)]. Density-like fields ('rho'/'density')
    are the MASS-CONSERVING estimate sum_j m_j (integral_cell W_j) / V_cell; other
    fields are the kernel-volume-weighted leaf average. The summed leaf mass
    equals the particle mass to high precision (the headline Petkova property),
    now on the adaptive cells. `_mass` (conserved per-leaf mass) is also exposed.
    """
    if not isinstance(grid, AMRGrid):
        raise TypeError("interpolate_amr_petkova requires an AMRGrid target")

    vals, names = particles.value_array(values)
    Ncell = grid.ncell
    nbx, nby, nbz = grid.base_npx
    blo = grid.lo
    bpw = grid.base_pw
    box = grid.box
    perx, pery, perz = grid.periodic

    nthreads = numba.get_num_threads()
    K = vals.shape[1]
    mass_local = np.zeros((nthreads, Ncell), dtype=np.float64)
    den_local = np.zeros((nthreads, Ncell), dtype=np.float64)
    num_local = np.zeros((nthreads, Ncell, K), dtype=np.float64)

    pos = particles.pos
    _deposit_petkova_amr(
        np.ascontiguousarray(pos[:, 0]),
        np.ascontiguousarray(pos[:, 1]),
        np.ascontiguousarray(pos[:, 2]),
        particles.h, particles.mass, particles.mass / particles.rho, vals,
        grid.cell_lo, grid.cell_size, grid.base_start, grid.base_count,
        blo[0], blo[1], blo[2], bpw[0], bpw[1], bpw[2],
        nbx, nby, nbz, perx, pery, perz, box[0], box[1], box[2], _CUBE_FACES,
        mass_local, den_local, num_local)

    mass = mass_local.sum(axis=0)
    den = den_local.sum(axis=0)
    num = num_local.sum(axis=0)
    vol = grid.volumes()
    mask = den > 0.0

    data = {}
    for p, name in enumerate(names):
        if name.lower() in ("rho", "density"):
            data[name] = mass / vol
        else:
            g = np.zeros(Ncell, dtype=np.float64)
            g[mask] = num[mask, p] / den[mask]
            data[name] = g
    data.setdefault("_mass", mass)
    return data


# --------------------------------------------------------------------------- #
#  Voronoi target: the NATIVE Petkova case — integrate the kernel exactly over an
#  arbitrary convex polyhedron. `_voro_cell_integral` is `_cube_integral`
#  generalized from the fixed 8-vertex cube to the CSR face/vertex layout of a
#  VoronoiGrid cell (vertices shifted to the particle frame); the per-face/edge
#  divergence-theorem construction and pyvoro's face winding are identical.
# --------------------------------------------------------------------------- #

@njit(cache=True)
def _voro_cell_integral(px, py, pz, h, gx, gy, gz, c,
                        cf_start, fv_start, fvx, fvy, fvz):
    """Integral of the unit-normalized kernel over Voronoi cell `c`, with the
    particle at (px,py,pz) and the cell's generator at (gx,gy,gz). Returns the
    contained kernel-mass fraction (>= 0).

    Signed tetrahedralization from the particle: one front/back sign per face
    (the generator gives the interior side to orient the outward normal),
    fan-triangulate each polygonal face, sum `_cone_tri`. This replaces the old
    foot-point wedge assembly (`_vertex_integral` + cond1/cond2 + sign_M), which
    failed to telescope for cells straddling the kernel break (h/2h) when the
    particle is outside the cell — the ~1-1.5% Voronoi conservation error on
    high-dynamic-range (collapsed-core) meshes. Now machine-precision conserving.
    """
    total = 0.0
    f0 = cf_start[c]; f1 = cf_start[c + 1]
    for f in range(f0, f1):
        i0 = fv_start[f]; i1 = fv_start[f + 1]
        nfv = i1 - i0
        if nfv < 3:
            continue
        # first three face vertices (particle-relative) -> face normal
        v0x = fvx[i0] - px;   v0y = fvy[i0] - py;   v0z = fvz[i0] - pz
        ax = fvx[i0+1] - px - v0x; ay = fvy[i0+1] - py - v0y; az = fvz[i0+1] - pz - v0z
        bx = fvx[i0+2] - px - v0x; by = fvy[i0+2] - py - v0y; bz = fvz[i0+2] - pz - v0z
        nx = ay*bz - az*by; ny = az*bx - ax*bz; nz = ax*by - ay*bx
        # generator relative to the particle, to orient the outward normal
        grx = gx - px; gry = gy - py; grz = gz - pz
        if nx*(grx - v0x) + ny*(gry - v0y) + nz*(grz - v0z) > 0.0:
            nx = -nx; ny = -ny; nz = -nz
        sgn = 1.0 if (v0x*nx + v0y*ny + v0z*nz) > 0.0 else -1.0
        # fan-triangulate the face polygon from vertex 0
        for k in range(1, nfv - 1):
            ka = i0 + k; kb = i0 + k + 1
            total += sgn * _cone_tri(
                v0x, v0y, v0z,
                fvx[ka] - px, fvy[ka] - py, fvz[ka] - pz,
                fvx[kb] - px, fvy[kb] - py, fvz[kb] - pz, h)
    if total < 0.0:
        total = 0.0
    return total


@njit(parallel=True, cache=True)
def _deposit_petkova_voronoi(x, y, z, h, mass, volw, vals,
                             genx, geny, genz, gen_rad,
                             cand_start, cand_count, cand_cell,
                             cf_start, fv_start, fvx, fvy, fvz,
                             mass_local, den_local, num_local):
    npart = x.shape[0]
    K = vals.shape[1]
    for i in prange(npart):
        mi = mass[i]
        hi = h[i]
        if mi <= 0.0 or hi <= 0.0:
            continue
        t = numba.get_thread_id()
        wi = volw[i]
        xi = x[i]; yi = y[i]; zi = z[i]
        kr = 2.0 * hi
        s0 = cand_start[i]; n0 = cand_count[i]
        for jj in range(s0, s0 + n0):
            c = cand_cell[jj]
            # prune: the kernel reaches cell c only if |particle-gen| < 2h + rad
            dx = genx[c] - xi; dy = geny[c] - yi; dz = genz[c] - zi
            reach = kr + gen_rad[c]
            if dx*dx + dy*dy + dz*dz >= reach * reach:
                continue
            frac = _voro_cell_integral(xi, yi, zi, hi,
                                       genx[c], geny[c], genz[c], c,
                                       cf_start, fv_start, fvx, fvy, fvz)
            if frac <= 0.0:
                continue
            mass_local[t, c] += mi * frac
            den_local[t, c] += wi * frac
            for p in range(K):
                num_local[t, c, p] += wi * vals[i, p] * frac


def interpolate_voronoi_petkova(particles, grid, values=None):
    """Petkova exact interpolation of `particles` onto a `VoronoiGrid`.

    The native Petkova case: each particle's kernel is integrated EXACTLY over
    every overlapping Voronoi polyhedron, so the summed cell mass equals the
    particle mass to high precision (the headline mass-conserving property), now
    on an unstructured mesh. Density-like fields ('rho'/'density') are
    sum_j m_j (integral_cell W_j)/V_cell; other fields are the kernel-volume-
    weighted cell average. `_mass` (conserved per-cell mass) is also exposed.

    Returns data: dict[name -> (Ncell,)].
    """
    if not isinstance(grid, VoronoiGrid):
        raise TypeError("interpolate_voronoi_petkova requires a VoronoiGrid target")

    vals, names = particles.value_array(values)
    Ntot = grid.ntotal                # deposit over ALL cells (incl. ghosts)
    # candidate cells: generator within 2h + max circumradius (necessary overlap
    # condition |particle-gen| < 2h + cell_rad; the integral itself -> 0 for any
    # marginal cell the kernel does not actually reach).
    Rmax = float(grid.cell_rad.max()) if Ntot else 0.0
    radii = 2.0 * np.asarray(particles.h, dtype=np.float64) + Rmax
    cand_start, cand_count, cand_cell = voronoi_candidates(
        grid.gen, particles.pos, radii)

    nthreads = numba.get_num_threads()
    K = vals.shape[1]
    mass_local = np.zeros((nthreads, Ntot), dtype=np.float64)
    den_local = np.zeros((nthreads, Ntot), dtype=np.float64)
    num_local = np.zeros((nthreads, Ntot, K), dtype=np.float64)

    pos = particles.pos
    _deposit_petkova_voronoi(
        np.ascontiguousarray(pos[:, 0]),
        np.ascontiguousarray(pos[:, 1]),
        np.ascontiguousarray(pos[:, 2]),
        particles.h, particles.mass, particles.mass / particles.rho, vals,
        np.ascontiguousarray(grid.gen[:, 0]),
        np.ascontiguousarray(grid.gen[:, 1]),
        np.ascontiguousarray(grid.gen[:, 2]), grid.cell_rad,
        cand_start, cand_count, cand_cell,
        grid.cf_start, grid.fv_start, grid.fvx, grid.fvy, grid.fvz,
        mass_local, den_local, num_local)

    # fold periodic ghost-cell deposits back onto the real (central) cells
    mass = grid.fold_to_central(mass_local.sum(axis=0))
    den = grid.fold_to_central(den_local.sum(axis=0))
    num = grid.fold_to_central(num_local.sum(axis=0))
    Ncell = grid.ncell
    vol = grid.volumes()
    mask = den > 0.0

    data = {}
    for p, name in enumerate(names):
        if name.lower() in ("rho", "density"):
            data[name] = mass / vol
        else:
            g = np.zeros(Ncell, dtype=np.float64)
            g[mask] = num[mask, p] / den[mask]
            data[name] = g
    data.setdefault("_mass", mass)
    return data


# --------------------------------------------------------------------------- #
#  Exact 2D column (Petkova projection): integrate the projected kernel Y(R,h)
#  exactly over each pixel square, instead of point-sampling it (SPH 2D).
# --------------------------------------------------------------------------- #

@njit(cache=True)
def _edge_angular(doh, psi, ptable, glx, glw):
    """E(d,psi) = int_0^psi M(d/cos t) dt with M(r) = (1/pi) P(r/h).

    doh = d/h (perpendicular distance / smoothing length); integrand bounded
    because the argument d/(h cos t) <= r/h <= 2 over the subtended range."""
    if psi <= 0.0:
        return 0.0
    acc = 0.0
    for k in range(glx.shape[0]):
        t = psi * glx[k]
        arg = doh / math.cos(t)
        acc += glw[k] * pcum(arg, ptable)
    return psi * acc / math.pi


@njit(cache=True)
def _pixel_column_integral(cx, cy, h, ptable, glx, glw):
    """Exact fraction of a particle's column mass falling in a pixel.

    cx,cy : (4,) pixel-corner coords RELATIVE to the particle's projected
            position (particle at origin), ordered around the pixel.
    Returns int_pixel Y(R,h) dA (in [0,1]); sum over a covering grid -> 1.
    """
    eps = 1e-12
    total = 0.0
    for e in range(4):
        ax = cx[e]; ay = cy[e]
        bx = cx[(e + 1) % 4]; by = cy[(e + 1) % 4]
        lx = bx - ax; ly = by - ay
        ll = math.sqrt(lx * lx + ly * ly)
        if ll < eps:
            continue
        d = abs(ax * ly - ay * lx) / ll          # perp distance O->line(A,B)
        r12 = math.sqrt(ax * ax + ay * ay)
        r13 = math.sqrt(bx * bx + by * by)
        psiA = math.acos(min(d / r12, 1.0)) if r12 > eps else 0.0
        psiB = math.acos(min(d / r13, 1.0)) if r13 > eps else 0.0
        EA = _edge_angular(d / h, psiA, ptable, glx, glw)
        EB = _edge_angular(d / h, psiB, ptable, glx, glw)
        # perpendicular foot inside the segment -> the two wedges add; outside ->
        # they subtract (same construction as the 3D vertex_integral combine).
        outside = (r12 * math.sin(psiA) + eps >= ll) or (r13 * math.sin(psiB) + eps >= ll)
        contrib = abs(EA - EB) if outside else (EA + EB)
        sgn = ax * by - ay * bx                   # signed triangle (O,A,B) area*2
        if sgn > 0.0:
            total += contrib
        elif sgn < 0.0:
            total -= contrib
    if total < 0.0:
        total = 0.0
    return total


@njit(parallel=True, cache=True)
def _deposit_petkova_proj(a1, a2, h, mass, volw, vals, nx, ny, lo0, lo1,
                          pw0, pw1, per0, per1, ptable, glx, glw,
                          a_vol, w_vol, a_mass, w_mass):
    """Exact 2D column deposit. Per-pixel weight = (int_pixel Y dA)/pixarea, i.e.
    the exact pixel-averaged projected kernel (cf. the point-sampled SPH version).
    Accumulators per-thread (nthreads, nx, ny[, K])."""
    npart = a1.shape[0]
    K = vals.shape[1]
    pixarea = pw0 * pw1
    for i in prange(npart):
        mi = mass[i]
        hi = h[i]
        if mi <= 0.0 or hi <= 0.0:
            continue
        t = numba.get_thread_id()
        vi = volw[i]
        p1 = a1[i]; p2 = a2[i]
        kr = 2.0 * hi
        imin = int(np.floor((p1 - kr - lo0) / pw0))
        imax = int(np.floor((p1 + kr - lo0) / pw0)) + 1
        jmin = int(np.floor((p2 - kr - lo1) / pw1))
        jmax = int(np.floor((p2 + kr - lo1) / pw1)) + 1
        if not per0:
            imin = max(imin, 0); imax = min(imax, nx)
        if not per1:
            jmin = max(jmin, 0); jmax = min(jmax, ny)
        cx = np.empty(4); cy = np.empty(4)
        for ix in range(imin, imax):
            iw = ix % nx if per0 else ix
            xe0 = lo0 + ix * pw0 - p1
            xe1 = xe0 + pw0
            for iy in range(jmin, jmax):
                ye0 = lo1 + iy * pw1 - p2
                ye1 = ye0 + pw1
                # CCW corners relative to particle
                cx[0] = xe0; cy[0] = ye0
                cx[1] = xe1; cy[1] = ye0
                cx[2] = xe1; cy[2] = ye1
                cx[3] = xe0; cy[3] = ye1
                ipix = _pixel_column_integral(cx, cy, hi, ptable, glx, glw)
                if ipix <= 0.0:
                    continue
                wgt = ipix / pixarea            # exact pixel-averaged column kernel
                jw = iy % ny if per1 else iy
                mw = mi * wgt
                vw = vi * wgt
                w_mass[t, iw, jw] += mw
                w_vol[t, iw, jw] += vw
                for p in range(K):
                    a = vals[i, p]
                    a_vol[t, iw, jw, p] += vw * a
                    a_mass[t, iw, jw, p] += mw * a


def interpolate_projection2d_petkova(particles, target, values=None):
    """Exact (Petkova) 2D column interpolation onto a Projection2D target.

    Supports the column-type modes ('column', 'average', 'rhocolumn'); 'slice' is
    not a column integral (use the SPH slice, or the 3D Petkova on a thin slab).
    Returns data: dict[name -> (nx,ny)] with the same conventions as the SPH 2D
    projection, but each particle's column is integrated EXACTLY over every pixel
    (mass-conserving per pixel, even when h < pixel).
    """
    if not isinstance(target, Projection2D):
        raise TypeError("requires a Projection2D target")
    if target.mode == "slice":
        raise NotImplementedError(
            "Petkova-2D is a column integral; use method='sph' for slice (or the "
            "3D Petkova on a thin slab).")

    vals, names = particles.value_array(values)
    nx, ny = target.npx
    lo = target.lo
    pw = target.pixwidth
    a1c = np.ascontiguousarray(particles.pos[:, target.axis1])
    a2c = np.ascontiguousarray(particles.pos[:, target.axis2])
    per0 = particles.periodic[target.axis1]
    per1 = particles.periodic[target.axis2]

    nthreads = numba.get_num_threads()
    K = vals.shape[1]
    a_vol = np.zeros((nthreads, nx, ny, K))
    w_vol = np.zeros((nthreads, nx, ny))
    a_mass = np.zeros((nthreads, nx, ny, K))
    w_mass = np.zeros((nthreads, nx, ny))

    _deposit_petkova_proj(a1c, a2c, particles.h, particles.mass,
                          particles.mass / particles.rho, vals, nx, ny,
                          lo[0], lo[1], pw[0], pw[1], per0, per1,
                          PCUM_TABLE, _GLX, _GLW, a_vol, w_vol, a_mass, w_mass)

    a_vol = a_vol.sum(axis=0); w_vol = w_vol.sum(axis=0)
    a_mass = a_mass.sum(axis=0); w_mass = w_mass.sum(axis=0)
    mv = w_vol > 0.0; mm = w_mass > 0.0

    data = {}
    for p, name in enumerate(names):
        is_dens = name.lower() in ("rho", "density")
        if target.mode == "column":
            data[name] = w_mass.copy() if is_dens else a_vol[..., p]
        elif target.mode == "average":
            g = np.zeros((nx, ny)); g[mv] = a_vol[..., p][mv] / w_vol[mv]
            data[name] = g
        else:  # rhocolumn
            g = np.zeros((nx, ny)); g[mm] = a_mass[..., p][mm] / w_mass[mm]
            data[name] = g
    return data
