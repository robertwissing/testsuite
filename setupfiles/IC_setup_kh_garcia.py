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


class setup_kh_garcia(object):
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
    cosmo=0
    time=0
    gamma=5./3.
    periodic=1
    deltastep=0.02
    nsteps=300
    dmsolunit=1.0
    dkpcunit=1.0
    adi=1
    grav=0
    molweight=1
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
    rhodisk=2.
    rhomean=1.5
    shape=2
    smooth=0
    v1=-0.5
    v2=0.5
    rloop=0.5
    
    def __init__(self):
        pass
    def create(self,nx,rhodiff=2.,smooth=0,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000'):
        self.rhodisk=self.rhozero*rhodiff
        self.smooth=smooth
        v1 = -0.5
        v2 = 0.5
        przero = 2.5
        gam1 = self.gamma-1
        B0=0.1
        C=4.*np.pi
        print('Setup for Kevin-Helmholtz..')
        self.rhomean=(self.rhozero+self.rhodisk)/2.
        self.dxbound=2.
        deltax=self.dxbound/nx
        deltadisk=deltax*(self.rhozero/self.rhodisk)**(1./3.)
        #deltamean = self.dxbound/nx
        #deltax = deltamean*(self.rhomean/self.rhozero)**(1./3.)
        #print(deltax,deltamean,deltadens)
        dx=self.dxbound/2.
        dy=dx
        dz=dx
        distribute.setbound(self,-dx,dx,-dy,dy,-dz,dz)
        costheta = self.dxbound/np.sqrt(self.dxbound**2 + self.dybound**2)
        sintheta = self.dybound/np.sqrt(self.dxbound**2 + self.dybound**2)
        deltadisk = deltax*(self.rhozero/self.rhodisk)**(1./3.)
        if(self.smooth==0):
            self.rcylmin2 = self.rloop**2
            distribute.setcloseddist(self,-dx,dx,-dy,dy,-dz,dz,deltax)
            self.rcylmin2 = 0.0
            self.rcylmax2 = self.rloop**2
            distribute.setcloseddist(self,-dx,dx,-dy,dy,-dz,dz,deltadisk)
        else:
             s#elf.rcylmax2 = 1.0
             distribute.setcloseddist(self,-dx,dx,-dy,dy,-dz,dz,deltax)
             fixdens.set_density_profile(self,0.0,rbig,rhofunc=True,icoord=0,igeom=2)
             totmass = self.rhomean*self.dxbound*self.dybound*self.dzbound
        #distribute2.setcloseddist(self,-dx,dx,-dy,dy,-dz,dz,deltadisk,np.sqrt(rhodiff))
        totmass = (self.rhozero*self.dxbound*self.dybound + (self.rhodisk-self.rhozero)*np.pi*self.rloop**2)*self.dzbound
     
        
        self.npart=len(self.x)
        self.ngas=self.npart
        #totmass = self.dxbound*self.dybound*self.dzbound*self.rhomean
        self.mass = [totmass/self.npart]*self.npart
        self.h=smth.getsmooth(self,200)
        print('npart = ',self.npart,' particle mass = ',self.mass[0],' total mass = ',totmass)
  
        for i in range(self.npart):
            r=np.sqrt(self.x[i]**2+self.y[i]**2)
            #vr=0.01*np.sin(C*z);
            vr=0.05*np.exp(-np.abs(r-0.5)/0.1)*np.sin(C*self.z[i]);
            #vr=1.0
            phi=np.arctan2(self.y[i],self.x[i])
            self.vx.append(vr*np.cos(phi))
            self.vy.append(vr*np.sin(phi))
            self.vz.append(self.getveli(i))
            self.Bx.append(0.)
            self.By.append(0.)
            self.Bz.append(B0)
            self.u.append(przero/(gam1*self.getrhoi(i)))
            self.rho[i]=self.getrhoi(i)
            
    def getrho(self,x,y):
        r=np.sqrt(x**2+y**2)
        if(self.smooth==1):
            delta=0.0025
            fac1 = (1. - 1./(1. + np.exp(2.*(r+self.rloop)/delta)))
            fac2 = (1. - 1./(1. + np.exp(2.*(self.rloop-r)/delta)))
            rampf = fac1*fac2
            return self.rhozero + rampf*(self.rhodisk - self.rhozero)
        else:
            if(r<=self.rloop):
                return self.rhodisk
            else:
                return self.rhozero

    def getveli(self,i):
        r=np.sqrt(self.x[i]**2+self.y[i]**2)
        if(self.smooth==1): 
            return self.v1 + self.rampfunc(i)*(self.v2 - self.v1)  
        else:
            if(r<=self.rloop):
                return self.v2
            else:
                return self.v1
    
    
    def getrhoi(self,i):
        r=np.sqrt(self.x[i]**2+self.y[i]**2)
        if(self.smooth==1): 
            return self.rhozero + self.rampfunc(i)*(self.rhodisk - self.rhozero)  
        else:
            if(r<=self.rloop):
                return self.rhodisk
            else:
                return self.rhozero
            
    
    def rampfunc(self,i):
        r=np.sqrt(self.x[i]**2+self.y[i]**2)
        delta=0.0025
        fac1 = (1. - 1./(1. + np.exp(2.*(r+self.rloop)/delta)))
        fac2 = (1. - 1./(1. + np.exp(2.*(self.rloop-r)/delta)))
        rampf = fac1*fac2
        return rampf