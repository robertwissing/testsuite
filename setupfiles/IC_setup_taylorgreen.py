#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 15 23:18:49 2018

@author: robertwi
"""
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
from numba import njit


class setup_taylorgreen(object):
    dICdensRsmooth = 0.1
    dICdensprofile = 1
    dICdensdir = 1
    dICdensR = 1.0 
    dICdensinner = 1.0 # calculated depending on mass
    dICdensouter = 1.0
    rhoit=0
    npart=0
    ngas=0
    ndark=0
    nstar=0
    ndim=3
    time=0
    gamma=5./3.
    deltastep=0.01
    periodic=1
    nsteps=200
    dmsolunit=1.0
    dkpcunit=1.0
    adi=0
    grav=0
    cosmo=0
    molweight=1.0
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
    Bpsi=[]
    soft=[]
    tform=[]
    metals=[]
    dxbound=[]
    dybound=[]
    dzbound=[]
    totvol=[]
    rcylmin2= 0
    rmin2= 0
    rcylmax2= 10**22
    rmax2= 10**22
    rhozero=1.
    shape=2

    
    def __init__(self):
        pass
    def create(self,nx,v0=0.1,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000'):
        przero = 1
        gam1 = self.gamma-1
        print('Setup for Taylor–Green vortex.. Isothermal')
        self.rhozero==1.0
        self.dxbound=1.
        deltax = self.dxbound/nx
        dx=self.dxbound/2.
        dy=dx
        dz=2*np.sqrt(6)*deltax
        #dz=dx/4
        distribute.setbound(self,-dx,dx,-dy,dy,-dz,dz)
        #distribute.setrandomdist(self,-dx,dx,-dy,dy,-dz,dz,deltax)
        distribute.setcloseddist(self,-dx,dx,-dy,dy,-dz,dz,deltax)
        #distribute.setcubicdist(self,-dx,dx,-dy,dy,-dz,dz,deltax)
        self.npart=len(self.x)
        self.ngas=self.npart
        totmass = self.dxbound*self.dybound*self.dzbound*self.rhozero
        self.mass = [totmass/self.npart]*self.npart
        self.h=smth.getsmooth(self,200)
        print('npart = ',self.npart,' particle mass = ',self.mass[0],' total mass = ',totmass)
  
        for i in range(self.npart):
            cx=2*np.pi*self.x[i]
            cy=2*np.pi*self.y[i]
            self.vx.append(v0*np.sin(cx)*np.cos(cy))
            self.vy.append(-v0*np.sin(cy)*np.cos(cx))
            self.vz.append(0.)
            self.Bx.append(0.)
            self.By.append(0.)
            self.Bz.append(0.)
            self.u.append(1.0)
            self.rho[i]=self.getrhoi(i)

    def getrho(self,x):
        return self.rhozero
    
    def getrhoi(self,i):
        return self.rhozero
    
