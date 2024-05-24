#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 01:13:46 2019

@author: robertwi
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 15 23:18:49 2018

@author: robertwi
"""
## Setup for ORZAG-TANG
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
from numba import njit
import random as random


class setup_khslab(object):
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
    periodic=1
    deltastep=0.02
    nsteps=200
    dmsolunit=1.0
    dkpcunit=1.0
    adi=1
    choice=1
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
    rhodens=2.
    shape=2
    smooth=0
    
    def __init__(self):
        pass
    def create(self,nx,choice,rhodiff=2.,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000'):
        self.rhodens=self.rhozero*rhodiff
        przero = 2.5
        gam1 = self.gamma-1
        print('Setup for Kevin-Helmholtz..')
        
        self.dxbound=1.0
        deltax = self.dxbound/nx
        deltadens = deltax*(self.rhozero/self.rhodens)**(1./3.)
        dx=self.dxbound/2.
        dz=2*np.sqrt(6)*deltax
        if(choice==1):    
            dy=dx*0.25
        else:
            dy=dx*0.5
            deltax=deltadens
        distribute.setbound(self,-dx,dx,-dy,dy,-dz,dz)
        distribute.setcloseddist(self,-dx,dx,-dy,dy,-dz,dz,deltax)
        self.npart=len(self.x)
        self.ngas=self.npart
        density=self.getrhoi(0)
        totmass = self.dxbound*self.dybound*self.dzbound*density
        self.mass = [totmass/self.npart]*self.npart
        self.h=smth.getsmooth(self,200)
        print('npart = ',self.npart,' particle mass = ',self.mass[0])
  
        for i in range(self.npart): 
            rand = random.uniform(-1.0*np.sqrt(przero/density),1.0*np.sqrt(przero/density))
            self.vx.append(rand)
            self.vy.append(rand)
            self.vz.append(rand)
            self.Bx.append(0.)
            self.By.append(0.)
            self.Bz.append(0.)
            self.u.append(przero/(gam1*self.getrhoi(i)))
            
    def getrho(self,x):
            if(self.choice==2):
                return self.rhodens
            else:
                return self.rhozero
    
    def getrhoi(self,i):
            if(self.choice==2):
                return self.rhodens
            else:
                return self.rhozero