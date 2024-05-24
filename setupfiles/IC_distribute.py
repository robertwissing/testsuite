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

def setcubicdist(D,xmin,xmax,ymin,ymax,zmin,zmax,delta):

        nx = int(D.dxbound/delta)
        ny = int(D.dybound/delta)
        nz = int(D.dzbound/delta)
        deltaz = D.dzbound/nz
        deltay = D.dybound/ny
        deltax = D.dxbound/nx
        print ('x',xmin,xmax,'y',ymin,ymax,'z',zmin,zmax)
        print ('dx',deltax,'dy',deltay,'dz',deltaz)

        xstart = 0.
        ystart = 0.
        zstart = 0.
        ipart = len(D.x)

        for k in range (nz):
           zi = zmin + (k+0.5)*deltaz + zstart
           for j in range (ny):
              yi = ymin + (j+0.5)*deltay + ystart
              for i in range(nx):
                 xi = xmin + (i+0.5)*deltax + xstart

                 if (checkboundary(D,xi,yi,zi)):
                     D.x.append(xi)
                     D.y.append(yi)
                     D.z.append(zi)
                     ipart=ipart+1
                     D.rho.append(D.getrhoi(ipart-1))
                 else:
                     pass
        return(D)

def setcloseddist(D,xmin,xmax,ymin,ymax,zmin,zmax,delta):

        if(xmax**2>D.rmax2):
            xmax=np.sqrt(D.rmax2+0.01*D.rmax2)
            xmin=-xmax
        if(ymax**2>D.rmax2):
            ymax=np.sqrt(D.rmax2+0.01*D.rmax2)
            ymin=-ymax
        if(zmax**2>D.rmax2):
            zmax=np.sqrt(D.rmax2+0.01*D.rmax2)
            zmin=-zmax
        deltax = delta
        deltay = delta*np.sqrt(3./4.)
        deltaz = delta*np.sqrt(6.)/3.
        delx = 0.5*delta
        dxbound=xmax-xmin
        dybound=ymax-ymin
        dzbound=zmax-zmin
        nx = int(int((1-np.finfo(float).eps)*dxbound/deltax)+1)
        ny = 2*int((int((1-np.finfo(float).eps)*dybound/deltay)+1)/2)
        nz = 3*int((int((1-np.finfo(float).eps)*dzbound/deltaz)+1)/3)
        deltay = dybound/ny
        deltaz = dzbound/nz
        dely = 1./3.*deltay  #psep*sqrt(3.)/6.
        npnew = nx*ny*nz
        print ('dxbound',dxbound,'dybound',dybound,'dzbound',dzbound)
        print('xmin',xmin,'ymin',ymin,'zmin',zmin)
        print ('dx',deltax,'dy',deltay,'dz',deltaz)
        print ('dy/dx: ',deltay/deltax,' dz/dx: ',deltaz/deltax)
        print ('nx', nx, 'ny', ny, 'nz', nz)

        k = 0
        l = 1
        m = 1
        ipart = len(D.x)
        for i in range(0,npnew):
           k = k + 1
           if (k > nx):
              k = 1
              l = l + 1
              if (l > ny):
                 l = 1
                 m = m + 1
                 if (m > nz):
                    k = 1
                    l = 1
                    nz = nz + 1

           xstart = xmin + 0.5*delx
           ystart = ymin + 0.5*dely
           zstart = zmin + 0.5*deltaz

           jy = np.mod(l, 2)
           jz = np.mod(m, 3)

           if (jz==0):  # 3rd layer
              ystart = ystart + 2.*dely
              if (jy==0):
                  xstart = xstart + delx
           elif (jz==2):# 2nd layer
              ystart = ystart + dely
              if (jy==1):
                  xstart = xstart + delx
           elif (jy==0):  # first layer, jz=1
              xstart = xstart + delx

           xi = xstart + float(k - 1)*deltax
           yi = ystart + float(l - 1)*deltay
           zi = zstart + float(m - 1)*deltaz
           
           #print*,'part spacing with the edges of the box ','x',(xpartmin-xmin)/deltax,(xmax-xpartmax)/deltax, &
           #    'y',(ypartmin-ymin)/deltay,(ymax-ypartmax)/deltay, &
           #    'z',(zpartmin-zmin)/deltaz,(zmax-zpartmax)/deltaz

       #--trim to fit radius. do not allow particles to have *exactly* rmax
       #  (this stops round-off error from giving non-zero centre of mass)

           if (checkboundary(D,xi,yi,zi)):
              D.x.append(xi)
              D.y.append(yi)
              D.z.append(zi)
              ipart=ipart+1
              D.rho.append(D.getrhoi(ipart-1))
                   
        print((-np.max(D.x)+xmax)/delx+(np.min(D.x)-xmin)/delx)
        print((-np.max(D.y)+ymax)/dely+(np.min(D.y)-ymin)/dely)
        print((-np.max(D.z)+zmax)/deltaz+(np.min(D.z)-zmin)/deltaz)
        print('npart = ',len(D.x))
        return(D)

def setrandomdist(D,xmin,xmax,ymin,ymax,zmin,zmax,delta):
        #
        #--initialise random number generator
        #
        iseed = -43587
        dxbound=xmax-xmin
        dybound=ymax-ymin
        dzbound=zmax-zmin
        nx = int(dxbound/delta)
        ny = int(dybound/delta)
        nz = int(dzbound/delta)
        deltax=delta
        npnew = nx*ny*nz
        ipart = len(D.x)
        for i in range (0,npnew):

           xi = xmin + np.random.rand()*dxbound
           yi = ymin + np.random.rand()*dybound
           zi = zmin + np.random.rand()*dzbound
        
           if (checkboundary(D,xi,yi,zi)):
              D.x.append(xi)
              D.y.append(yi)
              D.z.append(zi)
              ipart=ipart+1
              D.rho.append(D.getrhoi(ipart-1))
        return(D)

def setbound(D,xmin,xmax,ymin,ymax,zmin,zmax):
           D.dxbound = xmax - xmin
           D.dybound = ymax - ymin
           D.dzbound = zmax - zmin
           D.totvol = D.dxbound*D.dybound*D.dzbound
           return(D)

def get_ny_nz_closepacked(delta,ymin,ymax,zmin,zmax):
    deltay = delta*np.sqrt(3./4.)
    deltaz = delta*np.sqrt(6.)/3.

    ny = 2*int((int((1-np.finfo(float).eps)*(ymax-ymin)/deltay)+1)/2)
    nz = 3*int((int((1-np.finfo(float).eps)*(zmax-zmin)/deltaz)+1)/3)
    return ny,nz

def checkboundary(D,xi,yi,zi):
    if(D.shape == 1): ## FOR CUBE
        if ((xi**2 < D.rmax2 and yi**2 < D.rmax2 and zi**2 < D.rmax2) and (D.rmin2 <= xi**2 or D.rmin2 <= yi**2 or D.rmin2 <= zi**2)):
            return True
        else:
            return False
    else:
        rcyl2 = xi**2 + yi**2
        rr2   = rcyl2 + zi**2
        if (D.rmin2 <= rr2 < D.rmax2 and D.rcylmin2 <= rcyl2 < D.rcylmax2):
            return True
        else:
            return False
