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
    dICdensRsmooth = 0.1
    dICdensprofile = 1
    dICdensdir = 1
    dICdensR = 1.0 
    dICdensinner = 1.0 # calculated depending on mass
    dICdensouter = 1.0
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
    cosmo=1
    molweight=1.0
    gamma=5./3.
    periodic=1
    deltastep=0.0005
    freqout=10
    grav=0
    nsteps=100
    dmsolunit=1.0
    dkpcunit=1000.0
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
    wk=0.
    x0=0.
    xdot0=0.
    rhozero=1.
    drho=0.
    Bzero=1.
    dB=0.
    a=0.01
    kappa=1.
      
    def __init__(self):
        pass
    
    def create(self,nx,vzero=0.0,Lboxinkpc=1.0,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000'):
        
        
         #--setup parameters
         self.dkpcunit=Lboxinkpc
         przero = 6.0
         self.rhozero = 1.0
         gam1=self.gamma-1
         uzero=przero/(gam1*self.rhozero)
         r0=1./(4*np.sqrt(8))
         r02=r0*r0
         #r02=0.01*0.01
         Bz0=1./np.sqrt(np.pi*4)
         #boundaries    
         dz=0.5
         dy=0.5
         dx=0.5
         cs=np.sqrt(self.gamma*przero)
         vx0=0.0
         
         distribute.setbound(self,-dx,dx,-dy,dy,-dz,dz)
         deltax = self.dxbound/nx
         if distri==0:
             distribute.setcloseddist(self,-dx,dx,-dy,dy,-dz,dz,deltax)
         elif distri==1:
             distribute.setrandomdist(self,-dx,dx,-dy,dy,-dz,dz,deltax)
         else:
             tgdata,tddata,tsdata,data_header,time=tip.readtipsy(entry);
             self.x=tgdata[:,1]
             self.y=tgdata[:,2]
             self.z=tgdata[:,3]
             self.rho=tgdata[:,7]    
         
         self.npart=len(self.x)
         self.ngas=self.npart    
         totmass = self.rhozero*self.dxbound*self.dybound*self.dzbound
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
             #if(displace==1):   
                 #self.displacepart(i,deltaz*1.9)
             r2 = self.x[i] * self.x[i] + self.y[i] * self.y[i] + self.z[i] * self.z[i]
             Bx0=0.0;
             rat2=r2/r02
             if(rat2<1.0):
                Bx0=Bz0*(rat2*rat2*rat2*rat2-2*rat2*rat2+1)
             self.Bx.append(Bx0)
             self.By.append(0.)
             self.Bz.append(Bz0)
             self.vx.append(vx0)
             self.vy.append(vx0)
             self.vz.append(vx0)
             self.u.append(uzero)

             
    def getrhoi(self,i):
            return self.rhozero
        
    def displacepart(self,i,rs):
        r1=random.random()
        r2=random.random()
        r3=random.random()
        rabs=np.sqrt(r1*r1+r2*r2+r3*r3)
        self.x[i] = self.x[i] + rs*r1/rabs
        self.y[i] = self.y[i] + rs*r2/rabs
        self.z[i] = self.z[i] + rs*r3/rabs
        if (self.x[i]>self.dxbound/2):
            self.x[i]=self.x[i]-self.dxbound;
        if (self.x[i]<-self.dxbound/2):
            self.x[i]=self.x[i]+self.dxbound;
        if (self.y[i]>self.dybound/2):
            self.y[i]=self.y[i]-self.dybound;
        if (self.y[i]<-self.dybound/2):
            self.y[i]=self.y[i]+self.dybound;
        if (self.z[i]>self.dzbound/2):
            self.z[i]=self.z[i]-self.dzbound;
        if (self.z[i]<-self.dzbound/2):
            self.z[i]=self.z[i]+self.dzbound;