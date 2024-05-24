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


class setup_mhdisowave(object):
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
    grav=0
    cosmo=0
    molweight=1.0
    gamma=5./3.
    periodic=1
    deltastep=0.002
    nsteps=200
    dmsolunit=1.0
    dkpcunit=1.0
    adi=1
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
    def create(self,nx,beta=0.1,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000'):
        przero = 1.0
        gam1 = self.gamma-1
        B0=np.sqrt(2.0*przero/beta)
        C=4.*np.pi
        print('Setup forisolated mhd wave..')
        self.rhozero==1.0
        self.dxbound=4.
        deltax = self.dxbound/nx
        dx=self.dxbound/2.
        dy=dx*0.5
        dz=2*np.sqrt(6)*deltax
        distribute.setbound(self,-dx,dx,-dy,dy,-dz,dz)
        distribute.setcloseddist(self,-dx,dx,-dy,dy,-dz,dz,deltax)
        self.npart=len(self.x)
        self.ngas=self.npart
        totmass = self.dxbound*self.dybound*self.dzbound*self.rhozero
        self.mass = [totmass/self.npart]*self.npart
        self.h=smth.getsmooth(self,200)
        print('npart = ',self.npart,' particle mass = ',self.mass[0],' total mass = ',totmass)
  
        for i in range(self.npart):
            vxi = 0.01*np.exp(-(self.x[i]/(3.0*self.h[i]))**2)
            self.vx.append(vxi)
            self.vy.append(0.0)
            self.vz.append(0.0)
            self.Bx.append(B0)
            self.By.append(0.)
            self.Bz.append(0.)
            self.u.append(przero/(gam1*self.getrhoi(i)))
            self.rho[i]=self.getrhoi(i)

    def getrho(self,x):
        return self.rhozero
    
    def getrhoi(self,i):
        return self.rhozero