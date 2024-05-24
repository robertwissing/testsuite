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


class setup_rt(object):
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
    
    def __init__(self):
        pass
    def create(self,nx,rhodiff=2.):
        self.rhodens=self.rhozero*rhodiff
        przero = self.rhodens/self.gamma
        gam1 = self.gamma-1
        grav=0.5
        C=4.*np.pi
        B0=1/np.sqrt(C)
        dvy=0.025
        print('Setup for Rayleigh-Taylor...')
        
        self.dxbound=1.
        deltax = self.dxbound/nx
        dx=self.dxbound/4.
        dy=dx*2.
        dz=2*np.sqrt(6)*deltax
        distribute.setbound(self,-dx,dx,-dy,dy,-dz,dz)
        distribute.setcloseddist(self,-dx,dx,-dy,dy,-dz,dz,deltax)
        fixdens.set_density_profile(self,-dy,dy,1,icoord=1)
        
        self.npart=len(self.x)
        self.ngas=self.npart
        totmass = self.dxbound*self.dybound*self.dzbound*(self.rhozero+self.rhodens)/2.
        self.mass = [totmass/self.npart]*self.npart
        self.h=smth.getsmooth(self,200)
        print('npart = ',self.npart,' particle mass = ',self.mass[0])
  
        for i in range(self.npart):
            self.vx.append(0.)
            if(-0.2<self.y[i]<0.2):    
                self.vy.append(dvy*(1+np.cos(2*C*(self.x[i]+0.75)))*(1+np.cos(5*np.pi*self.y[i])))
            else:
                self.vy.append(0.)
            self.vz.append(0.)
            self.Bx.append(B0)
            self.By.append(0.)
            self.Bz.append(0.)
            self.u.append((przero-grav*self.getrhoi(i)*self.y[i])/(gam1*self.getrhoi(i)))
            
    def getrho(self,x):
        delta=0.025
        rampf = 1./(1. + np.exp(-2.*(x)/delta))
        
        return self.rhozero + rampf*(self.rhodens - self.rhozero)
    
    def getrhoi(self,i):
        return self.rhozero + self.rampfunc(i)*(self.rhodens - self.rhozero)  
    
    def rampfunc(self,i):
        delta=0.025
        
        rampf = 1./(1. + np.exp(-2.*(self.y[i])/delta))
        return rampf

@njit
def rhofunc(x,c):
        delta=0.05
        fac1 = (1. - 1./(1. + np.exp(2.*(x+0.25)/delta)))
        
        rampf = fac1*fac2
        return c[0] + rampf*(c[1] - c[0])