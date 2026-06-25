#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 15:46:10 2017

@author: robertwi
"""

import numpy as np
from numba import jit

"""
Routines containing functions to provide cubic/closed/random IC distributions for
particles. D is data array and ranges are the domain ranges and delta is the given
spacing.
"""

def setcubicdist_old(D,xmin,xmax,ymin,ymax,zmin,zmax,delta):
        dxbound = xmax - xmin
        dybound = ymax - ymin
        dzbound = zmax - zmin
        nx = int((1-np.finfo(float).eps)*dxbound/delta)
        ny = int((1-np.finfo(float).eps)*dybound/delta)
        nz = int((1-np.finfo(float).eps)*dzbound/delta)
        deltax = dxbound/nx
        deltay = dybound/ny
        deltaz = dzbound/nz
        print ('x',xmin,xmax,'y',ymin,ymax,'z',zmin,zmax)
        print ('dx',deltax,'dy',deltay,'dz',deltaz)
        xcentre=(xmax+xmin)/2.0
        ycentre=(ymax+ymin)/2.0
        zcentre=(zmax+zmin)/2.0
        xstart = 0.00001*deltax
        ystart = 0.00001*deltay
        zstart = 0.00001*deltaz
        ipart = len(D.x)
        print('dxbound', dxbound, 'dybound', dybound, 'dzbound', dzbound)
        print('xmin', xmin, 'ymin', ymin, 'zmin', zmin)
        print('xmax', xmax, 'ymax', ymax, 'zmax', zmax)
        print('dx', deltax, 'dy', deltay, 'dz', deltaz)
        print('dy/dx:', deltay / deltax, 'dz/dx:', deltaz / deltax)
        print('nx', nx, 'ny', ny, 'nz', nz)

        for k in range (nz):
           zi = zmin + k*deltaz + zstart
           for j in range (ny):
              yi = ymin + j*deltay + ystart
              for i in range(nx):
                 xi = xmin + i*deltax + xstart

                 if (checkboundary(D,xi,yi,zi,xcentre,ycentre,zcentre)):
                     D.x.append(xi)
                     D.y.append(yi)
                     D.z.append(zi)
                     ipart=ipart+1
                     D.rho.append(D.getrhoi(ipart-1))
                 else:
                     pass
        return(D)

def setcubicdist(D, xmin, xmax, ymin, ymax, zmin, zmax, delta, *, force_odd=False):
    def fit_axis(a_min, a_max, delta):
        L = float(a_max - a_min)
        # choose integer number of cells closest to request
        n = max(1, int(round(L / delta)))
        if force_odd and n % 2 == 0:  # ensure a center point if desired
            n += 1
        d = L / n                    # adjusted spacing that fits exactly
        # centers: a_min + (i+0.5)*d, i=0..n-1
        return n, d, a_min + (np.arange(n) + 0.5) * d

    nx, dx, xs = fit_axis(xmin, xmax, delta)
    ny, dy, ys = fit_axis(ymin, ymax, delta)
    nz, dz, zs = fit_axis(zmin, zmax, delta)

    xcentre = 0.5*(xmax + xmin)
    ycentre = 0.5*(ymax + ymin)
    zcentre = 0.5*(zmax + zmin)

    print('x', xmin, xmax, 'y', ymin, ymax, 'z', zmin, zmax)
    print('dx', dx, 'dy', dy, 'dz', dz)
    print('dy/dx:', dy/dx, 'dz/dx:', dz/dx)
    print('nx', nx, 'ny', ny, 'nz', nz)

    ipart0 = len(D.x)
    for zk in zs:
        for yj in ys:
            for xi in xs:
                if checkboundary(D, xi, yj, zk, xcentre, ycentre, zcentre):
                    D.x.append(xi); D.y.append(yj); D.z.append(zk)
                    # compute density for this new particle
                    D.rho.append(D.getrhoi(len(D.x)-1))
                # else: skip

    return D



def setcloseddist2(D, xmin, xmax, ymin, ymax, zmin, zmax, delta):
    if xmax**2 > D.rmax2:
        xmax = np.sqrt(D.rmax2 + 0.01 * D.rmax2)
        xmin = -xmax
    if ymax**2 > D.rmax2:
        ymax = np.sqrt(D.rmax2 + 0.01 * D.rmax2)
        ymin = -ymax
    if zmax**2 > D.rmax2:
        zmax = np.sqrt(D.rmax2 + 0.01 * D.rmax2)
        zmin = -zmax

    deltax = delta
    deltay = delta * np.sqrt(3. / 4.)
    deltaz = delta * np.sqrt(6.) / 3.

    dxbound = xmax - xmin
    dybound = ymax - ymin
    dzbound = zmax - zmin
    nx = int(dxbound / deltax)
    ny = 2 * (int(dybound / deltay)) // 2
    nz = 3 * (int(dzbound / deltaz)) // 3

    # Recalculate deltax to fit dxbound exactly
    deltax = dxbound / nx
    deltay = dybound / ny
    deltaz = dzbound / nz
    delx = 0.5 * deltax
    dely = deltay / 3.  # psep * sqrt(3.) / 6.

    npnew = nx * ny * nz
    print('dxbound', dxbound, 'dybound', dybound, 'dzbound', dzbound)
    print('xmin', xmin, 'ymin', ymin, 'zmin', zmin)
    print('xmax', xmax, 'ymax', ymax, 'zmax', zmax)
    print('dx', deltax, 'dy', deltay, 'dz', deltaz)
    print('dy/dx:', deltay / deltax, 'dz/dx:', deltaz / deltax)
    print('nx', nx, 'ny', ny, 'nz', nz)

    ipart = len(D.x)
    xstart = xmin + 0.5 * delx
    ystart = ymin + 0.5 * dely
    zstart = zmin + 0.5 * deltaz
    xcentre=(xmax+xmin)/2.0
    ycentre=(ymax+ymin)/2.0
    zcentre=(zmax+zmin)/2.0
    for m in range(1, nz + 1):
        jz = m % 3
        for l in range(1, ny + 1):
            jy = l % 2
            xadjust = xstart
            yadjust = ystart

            # Adjust based on layer (jz) and row (jy)
            if jz == 0:  # 3rd layer
                yadjust += 2. * dely
                if jy == 0:
                    xadjust += delx
            elif jz == 2:  # 2nd layer
                yadjust += dely
                if jy == 1:
                    xadjust += delx
            elif jz == 1:  # first layer
                 if jy == 0:  # even rows get offset
                    xadjust += delx

            for k in range(1, nx+1):
                xi = xadjust + (k - 1) * deltax
                yi = yadjust + (l - 1) * deltay
                zi = zstart + (m - 1) * deltaz
               
                if checkboundary(D, xi, yi, zi,xcentre,ycentre,zcentre):
                    D.x.append(xi)
                    D.y.append(yi)
                    D.z.append(zi)
                    ipart += 1
                    D.rho.append(D.getrhoi(ipart - 1))

    print("MAX x", np.max(D.x), xmax, "MIN x", np.min(D.x), xmin)
    print("MAX y", np.max(D.y), ymax, "MIN y", np.min(D.y), ymin)
    print("MAX z", np.max(D.z), zmax, "MIN z", np.min(D.z), zmin)
    print('npart =', len(D.x))
    return D

def setcloseddist(D, xmin, xmax, ymin, ymax, zmin, zmax, delta, *, clip_to_rmax=True, periodic=True):
    # Optional spherical clipping like your version
    if clip_to_rmax and xmax**2 > D.rmax2:
        xmax = np.sqrt(1.01 * D.rmax2); xmin = -xmax
    if clip_to_rmax and ymax**2 > D.rmax2:
        ymax = np.sqrt(1.01 * D.rmax2); ymin = -ymax
    if clip_to_rmax and zmax**2 > D.rmax2:
        zmax = np.sqrt(1.01 * D.rmax2); zmin = -zmax

    # Nearest-neighbour spacing along x
    dx_req = float(delta)
    dy_req = dx_req * np.sqrt(3.0) / 2.0         # hex row spacing
    dz_req = dx_req * np.sqrt(6.0) / 3.0         # layer spacing for close packing

    Lx = xmax - xmin; Ly = ymax - ymin; Lz = zmax - zmin
    if Lx <= 0 or Ly <= 0 or Lz <= 0:
        raise ValueError("Non-positive box extent")

    # Choose integer counts; keep tiling constraints (ny even, nz multiple of 3)
    nx = max(1, int(round(Lx / dx_req)))
    ny = max(2, 2 * int(round((Ly / dy_req) / 2.0)))     # enforce even
    nz = max(3, 3 * int(round((Lz / dz_req) / 3.0)))     # enforce ABC

    # Recompute exact spacings so centers land at min+0.5*Δ and max-0.5*Δ
    dx = Lx / nx
    dy = Ly / ny
    dz = Lz / nz

    # Center start points
    x0 = xmin + 0.5 * dx
    y0 = ymin + 0.5 * dy
    z0 = zmin + 0.5 * dz

    half_dx = 0.5 * dx

    print('dxbound', Lx, 'dybound', Ly, 'dzbound', Lz)
    print('dx', dx, 'dy', dy, 'dz', dz)
    print('dy/dx:', dy/dx, 'dz/dx:', dz/dx)
    print('nx', nx, 'ny', ny, 'nz', nz)

    xcentre = 0.5*(xmax + xmin)
    ycentre = 0.5*(ymax + ymin)
    zcentre = 0.5*(zmax + zmin)

    start_idx = len(D.x)
    added = 0
    for m in range(nz):
        # ABC layer offsets (base offsets are independent of row parity)
        jz = m % 3
        if jz == 0:      # A
            x_layer_off = 0.0
            y_layer_off = 0.0
        elif jz == 1:    # B
            x_layer_off = half_dx
            y_layer_off = dy / 3.0
        else:            # C
            x_layer_off = 0.0
            y_layer_off = 2.0 * dy / 3.0

        z = z0 + m * dz

        for l in range(ny):
            # Standard hex row shift: even rows get +dx/2 (always)
            x_row_off = half_dx if (l % 2 == 0) else 0.0
            y = y0 + l * dy + y_layer_off

            # Base x for this row and layer
            x_base = x0 + x_layer_off + x_row_off

            for k in range(nx):
                x = x_base + k * dx
                if checkboundary(D, x, y, z, xcentre, ycentre, zcentre):
                    D.x.append(x); D.y.append(y); D.z.append(z)
                    D.rho.append(D.getrhoi(len(D.x) - 1))
                    added += 1

    if periodic:
        # The close-packed hex/ABC offsets push the last columns/rows slightly
        # past xmax/ymax (up to half a cell). For a periodic box those points are
        # the periodic images of cells inside the box, so wrap them back into
        # [min,max). This keeps the particle count and the lattice spacing, and
        # avoids out-of-box particles that the stretch mapping cannot place
        # (target enclosed mass > total mass -> "not converged").
        Lx = xmax - xmin; Ly = ymax - ymin; Lz = zmax - zmin
        for idx in range(start_idx, len(D.x)):
            D.x[idx] = xmin + (D.x[idx] - xmin) % Lx
            D.y[idx] = ymin + (D.y[idx] - ymin) % Ly
            D.z[idx] = zmin + (D.z[idx] - zmin) % Lz

    print('npart =', len(D.x), '(added', added, ')')
    print("MAX x", np.max(D.x), xmax, "MIN x", np.min(D.x), xmin)
    print("MAX y", np.max(D.y), ymax, "MIN y", np.min(D.y), ymin)
    print("MAX z", np.max(D.z), zmax, "MIN z", np.min(D.z), zmin)
    return D

def setrandomdist(D,xmin,xmax,ymin,ymax,zmin,zmax,delta):
        # Vectorised random fill. The acceptance region (checkboundary) is often
        # MUCH smaller than the full [xmin,xmax]^3 domain - e.g. a sphere of
        # radius sqrt(rmax2) << the box - so drawing uniformly over the whole
        # domain and rejecting wastes ~(domain/region) of the draws, which is
        # catastrophic when the region is small AND delta is fine (the inner,
        # dense component). Instead draw inside the TIGHT bounding box of the
        # region (the domain intersected with [-Rmax, Rmax]^3, Rmax=sqrt(rmax2),
        # the radius cut checkboundary imposes) at the SAME number density, so
        # almost nothing is rejected, and apply the boundary test vectorised.
        dxbound = xmax - xmin
        dybound = ymax - ymin
        dzbound = zmax - zmin
        nx = int(dxbound/delta); ny = int(dybound/delta); nz = int(dzbound/delta)
        npnew = nx*ny*nz
        if npnew <= 0:
            return D
        ndens = npnew / (dxbound*dybound*dzbound)        # ~ 1/delta^3
        # tight bounding box of the acceptance region (origin-centred radius cut
        # |coord| < Rmax); Rmax unset/0 -> falls back to the full domain.
        rmax2 = float(getattr(D, 'rmax2', 0.0))   # may be a huge sentinel int
        Rmax = np.sqrt(rmax2) if rmax2 > 0 else np.inf
        gx0 = max(xmin, -Rmax); gx1 = min(xmax, Rmax)
        gy0 = max(ymin, -Rmax); gy1 = min(ymax, Rmax)
        gz0 = max(zmin, -Rmax); gz1 = min(zmax, Rmax)
        gdx = gx1 - gx0; gdy = gy1 - gy0; gdz = gz1 - gz0
        if gdx <= 0 or gdy <= 0 or gdz <= 0:
            return D
        npts = max(1, int(round(ndens*gdx*gdy*gdz)))
        xi = gx0 + np.random.rand(npts)*gdx
        yi = gy0 + np.random.rand(npts)*gdy
        zi = gz0 + np.random.rand(npts)*gdz
        # vectorised checkboundary (mirrors the scalar version below)
        x2 = xi*xi; y2 = yi*yi; z2 = zi*zi
        if D.shape == 1:        # cube / box-shell
            mask = ((x2 < D.rmax2) & (y2 < D.rmax2) & (z2 < D.rmax2) &
                    ((D.rmin2 <= x2) | (D.rmin2 <= y2) | (D.rmin2 <= z2)))
        else:                   # sphere / cylinder
            rcyl2 = x2 + y2; rr2 = rcyl2 + z2
            mask = ((D.rmin2 <= rr2) & (rr2 < D.rmax2) &
                    (D.rcylmin2 <= rcyl2) & (rcyl2 < D.rcylmax2))
        xs = xi[mask]; ys = yi[mask]; zs = zi[mask]
        start = len(D.x)
        D.x.extend(xs.tolist())
        D.y.extend(ys.tolist())
        D.z.extend(zs.tolist())
        # per-particle rho (getrhoi reads D.x[i]/y/z, now populated); only over
        # the ACCEPTED particles, as before - no longer over the rejected draws.
        for i in range(start, len(D.x)):
            D.rho.append(D.getrhoi(i))
        return(D)

def setbound(D,xmin,xmax,ymin,ymax,zmin,zmax):
           D.dxbound = xmax - xmin
           D.dybound = ymax - ymin
           D.dzbound = zmax - zmin
           D.totvol = D.dxbound*D.dybound*D.dzbound
           return(D)

def get_ny_nz_closepacked(delta,ymin,ymax,zmin,zmax):
    deltay = (delta)*np.sqrt(3./4.)
    deltaz = (delta)*np.sqrt(6.)/3.
    ny = 2 * (int((ymax-ymin) / deltay)) // 2
    nz = 3 * (int((zmax-zmin) / deltaz)) // 3
    #ny = 2*int((int((1-np.finfo(float).eps)*(ymax-ymin)/deltay)+1)/2)
    #nz = 3*int((int((1-np.finfo(float).eps)*(zmax-zmin)/deltaz)+1)/3)
    return ny,nz

def checkboundary(D,xi,yi,zi,xcentre,ycentre,zcentre):
    #boundarybool = (2*(xi-xcentre))**2 < D.dxbound**2 and (2*(yi-ycentre))**2 < D.dybound**2 and (2*(zi-zcentre))**2 < D.dzbound**2
    boundarybool = True
    if(D.shape == 1): ## FOR CUBE
        if ((xi**2 < D.rmax2 and yi**2 < D.rmax2 and zi**2 < D.rmax2) and (D.rmin2 <= xi**2 or D.rmin2 <= yi**2 or D.rmin2 <= zi**2) and boundarybool):
            return True
        else:
            return False
    else:
        rcyl2 = xi**2 + yi**2
        rr2   = rcyl2 + zi**2
        if (D.rmin2 <= rr2 < D.rmax2 and D.rcylmin2 <= rcyl2 < D.rcylmax2 and boundarybool):
            return True
        else:
            return False
