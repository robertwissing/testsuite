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

class setup_rt(object):
    dICdensRsmooth = 0.025
    dICdensprofile = 3
    dICdensdir = 3
    dICdensR = 0.25
    dICdensinner = 1.0 # calculated depending on mass
    dICdensouter = 2.0
    rhoit=1
    molweight=1.0
    inflow=0
    cosmo = 0
    periodic = 1
    ns=64
    deltastep=0.02
    adi = 1
    grav = 1
    nsteps=200
    freqout=20
    dmsolunit=1.0
    dkpcunit=1.0
    npart=0
    ngas=0
    ndark=0
    nstar=0
    ndim=3
    time=0
    gamma=1.4
    mass=[]
    x=[]
    y=[]
    z=[]
    vx=[]
    vy=[]
    vz=[]
    rho=[]
    u=[]
    h=[]
    P=[]
    Bx=[]
    By=[]
    Bz=[]
    dxbound=[]
    dybound=[]
    dzbound=[]
    totvol=[]
    rcylmin2= 0
    rmin2= 0
    rcylmax2= 10**22
    rmax2= 10**22
    rhozero=1.
    rhodens=2.
    shape=2
    _p_z_grid=None
    case=1
    boundary=0.25
    zmin = 0.0

    def __init__(self):
        pass
    def create(self,nx,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000',rhodiff=2.,case=1,nlowdens=2.0,smooth=1,beta=0.0,bdir=1):
        self.rhodens=self.rhozero*rhodiff
        self.dICdensouter=self.rhodens
        self.smooth=smooth
        self.case=case
        if smooth==0:
                self.dICdensprofile=2
        #przero = 0.3*self.rhodens/self.gamma
        gam1 = self.gamma-1
        # Gravity magnitude the hydrostatic pressure profile is integrated with.
        # MUST equal the run's body force: IC_createparamfile writes
        # bBodyForce=1 with dBodyForceConst=1.0, and pkdBodyForce applies
        # a_z = -dBodyForceConst*sign(z) (constant gravity toward the midplane;
        # mti integrates its profile with the same g=1.0). This was 0.5 --
        # only HALF the run's gravity was supported, so the slabs free-fell
        # and over-compressed (bug found + fixed 2026-06-09).
        grav=1.0
        C=4.*np.pi
        dvz=0.025
        print('Setup for Rayleigh-Taylor...')
        self.dxbound=0.5
        deltax = self.dxbound/nx
        deltadens = deltax*(self.rhozero/self.rhodens)**(1./3.)
        dx=self.dxbound/2.
        if case==1:
            dz=dx*2+self.dICdensR*(nlowdens-1.0)
        else:
            dz=dx*2+self.dICdensR*0.5*(rhodiff-1.0)+self.dICdensR*(nlowdens-1.0)
            print("MYDZ",dz)
        dy=2*np.sqrt(6)*deltax
        distribute.setbound(self,-dx,dx,-dy,dy,-dz,dz)
        self.boundary=dz-self.dICdensR
        forceboundary=self.boundary-self.dICdensR*nlowdens;
        print("FORCEBOUDNARY",forceboundary,"BOUNDARY",self.boundary)
        if distri==0:
            smooth=0
            self.smooth=0
            if case == 1:
                distribute.setcloseddist(self,-dx,dx,-dy,dy,self.boundary,dz,deltadens)
                distribute.setcloseddist(self,-dx,dx,-dy,dy,-self.boundary,self.boundary,deltax)
                distribute.setcloseddist(self,-dx,dx,-dy,dy,-dz,-self.boundary,deltadens)
            else:
                distribute.setcloseddist(self,-dx,dx,-dy,dy,self.boundary,dz,deltadens)
                distribute.setcloseddist(self,-dx,dx,-dy,dy,-dz,self.boundary,deltax)
        elif distri==1:
            if case == 1:
                distribute.setrandomdist(self,-dx,dx,-dy,dy,self.boundary,dz,deltadens)
                distribute.setrandomdist(self,-dx,dx,-dy,dy,-self.boundary,self.boundary,deltax)
                distribute.setrandomdist(self,-dx,dx,-dy,dy,-dz,-self.boundary,deltadens)
            else:
                distribute.setrandomdist(self,-dx,dx,-dy,dy,self.boundary,dz,deltadens)
                distribute.setrandomdist(self,-dx,dx,-dy,dy,-dz,self.boundary,deltax)
        else:
            tgdata,tddata,tsdata,data_header,time=tip.readtipsy(entry);
            self.x=tgdata[:,1]
            self.y=tgdata[:,2]
            self.z=tgdata[:,3]
            self.rho=tgdata[:,7]
   #     if(self.smooth==1):
   #         distribute.setcloseddist(self,-dx,dx,-dy,dy,-dz,dz,deltamean)
   #         fixdens.set_density_profile(self,-dy,dy,1,icoord=1)
        self.npart=len(self.x)
        self.ngas=self.npart
        if self.case ==1:
            totmass = self.dxbound*self.dybound*(2*self.dICdensR*self.rhodens+self.rhozero*(self.dzbound-2*self.dICdensR))
        else:
            totmass = self.dxbound*self.dybound*(self.dICdensR*self.rhodens+self.rhozero*(self.dzbound-self.dICdensR))

        self.mass = [totmass/self.npart]*self.npart

        self.z = [z_val - forceboundary for z_val in self.z]
        self.zmin = -dz - forceboundary
        self.dzbound = self.dzbound+forceboundary*2.0
        self.boundary = self.dICdensR*nlowdens
        if vm==1:
            print("Varying masses -> ONLY FOR RAND + RELAXING TO GLASS");
            for i in range(0,self.npart,2):
                dmrat=0.25*self.mass[i]*np.random.rand()
                self.mass[i]=self.mass[i]+dmrat
                self.mass[i+1]=self.mass[i+1]-dmrat
        self.h=smth.getsmooth(self,200)
        print('npart = ',self.npart,' particle mass = ',self.mass[0])
        zmax = self.boundary + self.dICdensR
        if self.case==1:
            zmin=0.0
        else:
            zmin=self.zmin
        self.compute_pressure_profile(grav,zmin,zmax)
        # Magnetic field from the plasma beta AT THE INTERFACE (beta=0 -> hydro):
        # B0 = sqrt(2 p(z_iface) / beta), uniform; bdir picks the geometry --
        # 1 = x (horizontal, the default: perpendicular to gravity, the classic
        # magnetic-RT configuration), 2 = y (horizontal, along the thin slab),
        # 3 = z (vertical, along gravity).
        if beta > 0.0:
            B0 = np.sqrt(2.0*self.p_at_z(self.boundary)/beta)
        else:
            B0 = 0.0
        bvec = {1: (B0, 0.0, 0.0), 2: (0.0, B0, 0.0),
                3: (0.0, 0.0, B0)}[int(bdir)]
        print('beta = ', beta, ' B0 = ', B0, ' direction = ', int(bdir))
        for i in range(self.npart):
            self.vx.append(0.0)
            self.vy.append(0.0)
            if(self.case==1):
                zterm=abs(self.z[i])
            else:
                zterm=self.z[i]
            if zterm < self.boundary+0.2 and zterm > self.boundary-0.2 :
                # dvz, C are from your code (amplitude / wavenumber)
                vz_pert = dvz * (1.0 + np.cos(2*C*(self.x[i]))) * (1.0 + np.cos(5*np.pi*(zterm-self.boundary)))
                self.vz.append(vz_pert)
            else:
                self.vz.append(0.0)
           
            self.Bx.append(bvec[0])
            self.By.append(bvec[1])
            self.Bz.append(bvec[2])
            p_i = self.p_at_z(self.z[i])  # uses hydrostatic integral8
            u_i = p_i / ( gam1 * self.rho[i] )
            #u_i = (przero - grav * self.getrhoi(i) * (abs(self.z[i])-self.dICdensR)) / (gam1 * self.getrhoi(i))
            self.u.append(u_i)
            #if distri==2:
            #    self.u.append((przero-grav*self.getrhoi(i)*np.abs(self.z[i]))/(gam1*self.rho[i]))
            #else:
            #    self.u.append((przero-grav*self.getrhoi(i)*np.abs(self.z[i]))/(gam1*self.getrhoi(i)))
        if distri==0 or distri==2:
            self.dzbound = self.dzbound*2
    def rampfunc(self,z):
        rampf=0.5*(np.tanh((z+self.dICdensR)/self.dICdensRsmooth)-np.tanh((z-self.dICdensR)/self.dICdensRsmooth))
        return rampf

    def getrho(self, z):
        """Density as function of z (position). Uses ramp to blend."""
        if(self.case==1):
                zterm=abs(z)
        else:
                zterm=z
        if self.smooth == 1:
            r = self.rhodens+ self.rampfunc(z)*(self.rhozero - self.rhodens)
            return r
        else:
            return self.rhozero if zterm <= self.boundary else self.rhodens

    def getrhoi(self, i):
        """Wrapper for index-based access."""
        return self.getrho(self.z[i])

    def compute_pressure_profile(self, grav, zmin, zmax, nz=4001):
        print("ZMIN",zmin,"ZMAX",zmax,"DZBOUND",self.dzbound)

        zgrid = np.linspace(zmin, zmax, nz)
        a_grid = -grav * np.sign(zgrid)
        rho_grid = np.array([self.getrho(z) for z in zgrid])
    
        # trapezoidal cumulative integral: p(z) = p0 + integral_0^z rho(z') * a(z') dz'
        integrand = rho_grid * a_grid
        dp = np.cumsum(0.5*(integrand[1:] + integrand[:-1]) * np.diff(zgrid))
        
        p0x = 1E-6
        p0 = -np.min(dp)+p0x
        dp = np.insert(dp, 0, 0.0)
        p_grid = p0 + dp
        
        # store interpolator (simple numpy interp)
        self._p_z_grid = (zgrid, p_grid)
    def p_at_z(self, z):
        if self.case==1:
            zabs = abs(z)
        else:
            zabs = z
        zgrid, p_grid = self._p_z_grid
        return np.interp(zabs, zgrid, p_grid)
