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
import readtipsy as tip

class setup_kh(object):
    dICdensRsmooth = 0.025
    dICdensprofile = 3
    dICdensdir = 2
    dICdensR = 0.25
    dICdensinner = 2.0 # calculated depending on mass
    dICdensouter = 1.0
    rhoit=0
    ns=64
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
    freqout=20
    dmsolunit=1.0
    dkpcunit=1.0
    adi=1
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
    rhodens=2.
    rhomean=1.5
    shape=2
    smooth=0
    v1=-0.5
    v2=0.5
    
    def __init__(self):
        pass
    def create(self,nx,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000',rhodiff=2.,mach=np.sqrt(6)/5,smooth=1,B0=0.0,Bdir=1):
        self.rhodens=self.rhozero*rhodiff
        self.dICdensinner=self.rhodens
        self.smooth=smooth
        if smooth==0:
            self.dICdensprofile=2
        v1 = -0.5
        v2 = -v1
        gam1 = self.gamma-1
        self.przero = 1/(self.gamma*mach**2)
        lambdawave=0.5
        C=2.*np.pi/lambdawave
        deltavy=0.01
        print('Setup for Kevin-Helmholtz..')
        self.rhomean=(self.rhozero+self.rhodens)/2.
        self.dxbound=1.
        deltax = self.dxbound/nx
        deltadens = deltax*(self.rhozero/self.rhodens)**(1./3.)
        dx=self.dxbound/2.
        dy=dx
        dz=2*np.sqrt(6)*deltax
        distribute.setbound(self,-dx,dx,-dy,dy,-dz,dz)
        if distri==0:
            distribute.setcloseddist(self,-dx,dx,-dy,-dy*0.5,-dz,dz,deltax)
            distribute.setcloseddist(self,-dx,dx,dy*0.5,dy,-dz,dz,deltax)
            distribute.setcloseddist(self,-dx,dx,-dy*0.5,dy*0.5,-dz,dz,deltadens)
        elif distri==1:
            distribute.setrandomdist(self,-dx,dx,-dy,-dy*0.5,-dz,dz,deltax)
            distribute.setrandomdist(self,-dx,dx,dy*0.5,dy,-dz,dz,deltax)
            distribute.setrandomdist(self,-dx,dx,-dy*0.5,dy*0.5,-dz,dz,deltadens)
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
        totmass = self.dxbound*self.dybound*self.dzbound*self.rhomean
        self.mass = [totmass/self.npart]*self.npart
        if vm==1:
            print("Varying masses -> ONLY FOR RAND + RELAXING TO GLASS");
            for i in range(0,self.npart,2):
                dmrat=0.25*self.mass[i]*np.random.rand()
                self.mass[i]=self.mass[i]+dmrat
                self.mass[i+1]=self.mass[i+1]-dmrat
        self.h=smth.getsmooth(self,200) #Redundant
        print('npart = ',self.npart,' particle mass = ',self.mass[0],' total mass = ',totmass)

        Bx0=By0=Bz0=0.0;
        if Bdir == 1:
            Bx0=B0;
        if Bdir == 2:
            By0=B0;
        if Bdir == 3:
            Bz0=B0;
        else:
            print('choose valid direction')
            exit
        
        for i in range(self.npart):                
            self.vx.append(self.getveli(i))
            self.vy.append(deltavy*np.sin(C*self.x[i]))
            self.vz.append(0.)
            self.Bx.append(Bx0)
            self.By.append(By0)
            self.Bz.append(Bz0)
            if distri==2:
                self.u.append(self.przero/(gam1*self.rho[i]))
            else:
                self.u.append(self.przero/(gam1*self.getrhoi(i)))

    def getrho(self,x):
        if(self.smooth==1):
            
            delta=self.dICdensRsmooth
            fac1 = (1. - 1./(1. + np.exp(2.*(x+self.dICdensR)/delta)))
            fac2 = (1. - 1./(1. + np.exp(2.*(self.dICdensR-x)/delta)))
        
            rampf = fac1*fac2
            return self.rhozero + rampf*(self.rhodens - self.rhozero)
        else:
            if(abs(x)<=self.dICdensR):
                return self.rhodens
            else:
                return self.rhozero

    def getveli(self,i):
         if(self.smooth==1): 
            return self.v1 + self.rampfunc(i)*(self.v2 - self.v1)  
         else:
            if(abs(self.y[i])<=self.dICdensR):
                return self.v2
            else:
                return self.v1
    
    
    def getrhoi(self,i):
        if(self.smooth==1): 
            return self.rhozero + self.rampfunc(i)*(self.rhodens - self.rhozero)  
        else:
            if(abs(self.y[i])<=self.dICdensR):
                return self.rhodens
            else:
                return self.rhozero
    
    def rampfunc(self,i):
        rampf=0.5*(np.tanh((self.y[i]+self.dICdensR)/self.dICdensRsmooth)-np.tanh((self.y[i]-self.dICdensR)/self.dICdensRsmooth))
        return rampf
