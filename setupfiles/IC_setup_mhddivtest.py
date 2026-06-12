#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 17:02:08 2018

@author: robertwi
"""

import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
import random
from numba import njit
import readtipsy as tip

class setup_mhddivtest(object):
    dICdensRsmooth = 0.01
    dICdensprofile = 3
    dICdensdir = 1
    dICdensR = 1.0 
    dICdensinner = 1.0 # calculated depending on mass
    dICdensouter = 1.0
    inflow=0
    rhoit=0
    npart=0
    rhoit=0
    ns=64
    ngas=0
    ndark=0
    nstar=0
    ndim=3
    time=0
    grav=0
    H0=2.894405
    cosmo=0
    molweight=1.0
    gamma=5./3.
    periodic=1
    deltastep=0.1
    freqout=10
    grav=0
    nsteps=100
    dmsolunit=1.0
    dkpcunit=1.0
    adi=1
    shape=1
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
    rhodens=1.
    smooth=0.
    def __init__(self):
        pass
    
    def create(self,nx,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000',rhodiff=1.0,smooth=0.0,square=1,cosmo=0.0):
         #--setup parameters
         przero = 6.0
         self.rhozero = 1.0
         self.rhodens=self.rhozero*rhodiff
         self.dICdensinner=self.rhozero
         self.dICdensouter=self.rhodens
         self.smooth=smooth
         gam1=self.gamma-1
         uzero=przero/(gam1*self.rhozero)
         if int(smooth)==0:
            self.dICdensprofile=2
         r0=1./(np.sqrt(8))
         r02=r0*r0
         #Bz0=0.01*1./np.sqrt(np.pi*4)
         Bz0=1./np.sqrt(np.pi*4)
        #boundaries    
         dz=1.0
         dy=1.0
         dx=1.0
         self.dxbound = (dx+dx)
         deltax = self.dxbound/nx
         deltadens = deltax*(self.rhozero/self.rhodens)**(1./3.)
         if square==1:
             dz=2*np.sqrt(6)*deltax
         self.dICdensR = dx*0.5 #offset density contrast
         print(self.dICdensR,dx,-dx)
         distribute.setbound(self,-dx,dx,-dy,dy,-dz,dz)
         if distri==0:
             distribute.setcubicdist(self,-self.dICdensR,self.dICdensR,-dy,dy,-dz,dz,deltax)
             distribute.setcubicdist(self,self.dICdensR,dx,-dy,dy,-dz,dz,deltadens)
             distribute.setcubicdist(self,-dx,-self.dICdensR,-dy,dy,-dz,dz,deltadens)
         elif distri==1:
             distribute.setrandomdist(self,-self.dICdensR,self.dICdensR,-dy,dy,-dz,dz,deltax)
             distribute.setrandomdist(self,self.dICdensR,dx,-dy,dy,-dz,dz,deltadens)
             distribute.setrandomdist(self,-dx,-self.dICdensR,-dy,dy,-dz,dz,deltadens)
         else:
             tgdata,tddata,tsdata,data_header,time=tip.readtipsy(entry);
             self.x=tgdata[:,1]
             self.y=tgdata[:,2]
             self.z=tgdata[:,3]
             self.rho=tgdata[:,7]    
         
         self.npart=len(self.x)
         self.ngas=self.npart    
         totmass = 0.5*(self.rhozero+self.rhodens)*self.dxbound*self.dybound*self.dzbound
         self.mass = [totmass/self.npart]*self.npart
         if vm==1:
             print("Varying masses -> ONLY FOR RAND + RELAXING TO GLASS");              
             for i in range(0,self.npart,2):
                 dmrat=0.25*self.mass[i]*np.random.rand()
                 self.mass[i]=self.mass[i]+dmrat
                 self.mass[i+1]=self.mass[i+1]-dmrat
         print('npart = ',self.npart)
         self.h=smth.getsmooth(self,200)
         for i in range(self.npart):
             if square==1:
                r2 = self.x[i] * self.x[i] + self.y[i] * self.y[i]
             else:
                r2 = self.x[i] * self.x[i] + self.y[i] * self.y[i] + self.z[i] * self.z[i]
             Bx0=0.0;
             rat2=r2/r02
             if(rat2<1.0):
                Bx0=Bz0*(rat2*rat2*rat2*rat2-2*rat2*rat2+1)
             self.Bx.append(Bx0)
             self.By.append(0.)
             self.Bz.append(Bz0)
             self.vx.append(0.0)
             self.vy.append(0.0)
             self.vz.append(0.0)
             if distri==2:
                self.u.append(przero/(gam1*self.rho[i]))
             else:
                self.u.append(przero/(gam1*self.getrhoi(i)))

    def rampfunc(self,x):
        rampf=0.5*(np.tanh((x+self.dICdensR)/self.dICdensRsmooth)-np.tanh((x-self.dICdensR)/self.dICdensRsmooth))
        return rampf

    def getrho(self, x):
        xterm=abs(x)
        """Density as function of z (position). Uses ramp to blend."""
        if self.smooth == 1:
            r = self.rhodens + self.rampfunc(x)*(self.rhozero - self.rhodens)
            return r
        else:
            return self.rhozero if xterm <= self.boundary else self.rhodens

    def getrhoi(self, i):
        """Wrapper for index-based access."""
        return self.getrho(self.x[i])
