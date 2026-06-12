#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 27 14:57:12 2018

@author: robertwi
"""

## Setup for ORZAG-TANG
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
from numba import njit
import readtipsy as tip
import IC_denstable as denstable

class setup_mti(object):
    dICdensRsmooth = 0.025
    dICdensprofile = 8
    dICdensdir = 3
    dICdensR = 0.25
    dICdensinner = 1.0 # calculated depending on mass
    dICdensouter = 2.0
    rhoit=1
    molweight=1.0
    cosmo = 0
    periodic = 1
    ns=64
    deltastep=0.02
    adi = 1
    grav = 1
    nsteps=500
    freqout=10
    dmsolunit=1.0
    dkpcunit=1.0
    npart=0
    ngas=0
    ndark=0
    nstar=0
    ndim=3
    time=0
    shape=2
    inflow=0
    gamma=1.66666667
    mass=[]
    x=[]
    y=[]
    z=[]
    vx=[]
    vy=[]
    vz=[]
    rho=[]
    u=[]
    metals=[]
    h=[]
    P=[]
    Bx=[]
    By=[]
    Bz=[]
    dxbound=[]
    dybound=[]
    dzbound=[]
    totvol=[]
    gam1 = 0.66667
    _p_z_grid=None
    _rho_z_grid=None
    rcylmin2= 0
    rmin2= 0
    rcylmax2= 10**22
    rmax2= 10**22
    R = 0.2
    distri = 1

    def __init__(self):
        pass
    def create(self,nx,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000',machno=0.01):
        self.gam1 = self.gamma-1
        self.distri = distri
        grav=1.0
        csound_bound = np.sqrt(1.5*self.gamma*self.gam1) #sound at bottom of unstable region
        v0 = machno*csound_bound
        beta0=2.0*10**8
        B0=np.sqrt(0.5*self.gam1*1.1*1.5 / beta0)
        t_buoy = 1.0/np.sqrt(grav*1.0/3.0)
        self.deltastep = 0.01*t_buoy
        self.nsteps = 5000
        print('Setup for MTI...')
        self.dxbound=0.1
        deltax = self.dxbound/nx
        dx=self.dxbound/2
        dy=2.0*np.sqrt(6)*deltax
        dz=0.2
        
        pzero = 0.886
        zgrid = np.linspace(-dz, dz, 20001)
        self.compute_pressure_profile_u(grav,-dz,dz,pzero,zgrid)
        self._rho_z_grid =(zgrid, self.getrho_vec(zgrid))

        self.R=dz #domain size table
        distribute.setbound(self,-dx,dx,-dy,dy,-dz,dz)
        if self.distri==0:
            distribute.setcloseddist(self,-dx,dx,-dy,dy,-dz,dz,deltax)
            #fixdens.set_density_profile(self,-dz,dz,1,icoord=2)
            fixdens.set_density_profile(self,-dz,dz,rhotab=self.getrho_vec(zgrid),xtab=zgrid,icoord=2)
        elif self.distri==1:
            distribute.setrandomdist(self,-dx,dx,-dy,dy,-dz,dz,deltax)
        else:
            tgdata,tddata,tsdata,data_header,time=tip.readtipsy(entry);
            self.x=tgdata[:,1]
            self.y=tgdata[:,2]
            self.z=tgdata[:,3]
            self.rho=tgdata[:,7]
            self.mass= tgdata[:,0]
        self.npart=len(self.x)
        self.ngas=self.npart

        totmass = np.trapz(self.getrho_vec(zgrid),zgrid)*self.dxbound*self.dybound;
        print("Numerically integrated mass " , totmass)
        ntab=10000
        dimtab=1
        self.Rtab = np.sqrt(self.dxbound*self.dxbound+self.dybound*self.dybound+self.dzbound*self.dzbound)
        gridtab, rhotab = denstable.generate_density_table(self,n=ntab,dim=dimtab)
        denstable.write_density_table_xdr(self,"densitytable_xdr",gridtab,rhotab,n=ntab,dim=dimtab)
        self.mass = [totmass/self.npart]*self.npart
        
        if vm==1:
            print("Varying masses -> ONLY FOR RAND + RELAXING TO GLASS");
            for i in range(0,self.npart,2):
                dmrat=0.25*self.mass[i]*np.random.rand()
                self.mass[i]=self.mass[i]+dmrat
                self.mass[i+1]=self.mass[i+1]-dmrat

        #self.compute_pressure_profile(grav,-dz,dz,pzero)
        self.h=smth.getsmooth(self,200)
        print('npart = ',self.npart,' particle mass = ',self.mass[0])
        for i in range(self.npart):
            #self.vx.append(v0*np.sin(3*np.pi*self.z[i]/self.dxbound))
            self.vx.append(0.0)
            self.vy.append(0.0)
            self.vz.append(v0*np.sin(4*np.pi*self.x[i]/self.dxbound))
            #self.vz.append(v0)
            #self.vz.append(v0*np.sin(4*np.pi*self.x[i]/self.dxbound))
            self.Bx.append(B0)
            self.By.append(0.0)
            self.Bz.append(0.0)
            u_i = self.getui(i)
            bound_i = self.getbound(self.z[i])
            self.metals.append(bound_i)
            self.u.append(u_i)
        print('rhomean',np.mean(self.rho))

    def getbound(self, z):
        zterm=np.abs(z)
        bound1 = 0.1*0.5
        bound2 = 0.3*0.5
        if(zterm < bound1):
            r = 1.0
        elif (zterm > bound2):
            r = 1.0
        else:
            r = 0.0

        return r

    def getrho(self, z):
        zterm=np.abs(z)
        bound1 = 0.1*0.5
        bound2 = 0.3*0.5
        if(zterm < bound1):
            r = 1.1*(np.exp(-zterm))
        elif (zterm > bound2):
            r = 0.977758*np.exp(-30./29.*(zterm-bound2))
        else:
            r = 1.1*np.exp(-bound1)*( 1 - (zterm-bound1)/3.0 )**2
        
        if self.distri != 55:
            r = self.p_at_z(z)/( self.gam1 * self.getu(z) )

        return r

    def getu(self, z):
        zterm=np.abs(z)
        bound1 = 0.1*0.5
        bound2 = 0.3*0.5
        if(zterm < bound1):
            u = 1.5
        elif (zterm > bound2):
            u = 1.45
        else:
            u = 1.5*( 1 - (zterm-bound1)/3.0 )
        return u

    def getui(self, i):
        return self.getu(self.z[i])

    def getu_vec(self, z):
        zterm = np.abs(z)
        bound1 = 0.1*0.5
        bound2 = 0.3*0.5

        u = np.empty_like(zterm, dtype=float)
        m1 = zterm < bound1
        m2 = zterm > bound2
        mmid = ~(m1 | m2)

        u[m1] = 1.5
        u[m2] = 1.45
        u[mmid] = 1.5*( 1 - (zterm[mmid]-bound1)/3.0 )

        return u

    def getrhoi(self, i):
        if self.distri != 55:
            return self.p_at_z(self.z[i])/( self.gam1 * self.getu(self.z[i]) )
        else:
            return self.getrho(self.z[i])

    def getrho_vec(self, z):
        zterm = np.abs(z)
        bound1 = 0.1*0.5
        bound2 = 0.3*0.5

        r = np.empty_like(zterm, dtype=float)
        m1 = zterm < bound1
        m2 = zterm > bound2
        mmid = ~(m1 | m2)

        r[m1] = 1.1 * np.exp(-zterm[m1])
        r[m2] = 0.977758 * np.exp(-30./29.*(zterm[m2] - bound2))
        r[mmid] = 1.1 * np.exp(-bound1) * (1 - (zterm[mmid]-bound1)/3.0)**2

        if self.distri != 55:
            r = self.p_at_z(z) / ( self.getu_vec(z) * self.gam1)
        return r

    def compute_pressure_profile(self, grav, zmin, zmax, pzero, nz=10001):
        print("ZMIN",zmin,"ZMAX",zmax,"DZBOUND",self.dzbound)

        zgrid = np.linspace(zmin, zmax, nz)
        a_grid = -grav * np.sign(zgrid)
        rho_grid = np.array([self.getrho(z) for z in zgrid])
        # p(z) = p0 + integral_0^z rho(z') * a(z') dz'
        integrand = rho_grid * a_grid;
        dp = np.cumsum(0.5*(integrand[1:] + integrand[:-1]) * np.diff(zgrid))
        p0 = pzero
        dp = np.insert(dp, 0, 0.0)
        p_grid = p0 + dp
        if p_grid[-1] < 0.0:
            p_grid -= p_grid[-1]
        print("p(zmin), p(zmax), jump:", p_grid[0], p_grid[-1], p_grid[-1]-p_grid[0])
        # store interpolator
        self._p_z_grid = (zgrid, p_grid)

    def compute_pressure_profile_u(self, grav, zmin, zmax, p0, zgrid):
        """
        Compute and store hydrostatic-equilibrium p, rho, u on a z-grid.
        Also stores a tuple self._p_z_grid for compatibility with p_at_z.
        """
        nz = len(zgrid)
        a_grid = -grav * np.sign(zgrid)
        u_grid = self.getu_vec(zgrid)

        p_grid = np.empty_like(zgrid, dtype=float)
        p_grid[0] = p0
        for i in range(1, nz):
            dz = zgrid[i] - zgrid[i-1]
            coeff = 0.5 * (
                a_grid[i] / ((self.gamma - 1.0) * u_grid[i]) +
                a_grid[i-1] / ((self.gamma - 1.0) * u_grid[i-1])
            )
            p_grid[i] = p_grid[i-1] * np.exp(coeff * dz)

        rho_grid = p_grid / ((self.gamma - 1.0) * u_grid)

        # interpolator
        self._p_z_grid = (zgrid, p_grid)


    def p_at_z(self, z):
            zgrid, p_grid = self._p_z_grid
            return np.interp(z, zgrid, p_grid)
