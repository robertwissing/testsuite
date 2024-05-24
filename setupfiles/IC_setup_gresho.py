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
import readtipsy as tip

class setup_gresho(object):
    dICdensRsmooth = 0.1
    dICdensprofile = 1
    dICdensdir = 1
    dICdensR = 1.0 
    dICdensinner = 1.0 # calculated depending on mass
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
    deltastep=0.01
    freqout=100
    nsteps=100
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
    shape=2
    przero=5.0
    
    def __init__(self):
        pass
    def create(self,nx,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000',mach=np.sqrt(3/25)):
        self.przero = 1/(self.gamma*mach**2)
        self.deltastep=self.deltastep*(mach/np.sqrt(3/25))
        gam1 = self.gamma-1
        C=4.*np.pi
        print('Setup for Gresho vortex..')
        self.rhozero==1.0
        self.dxbound=1.
        deltax = self.dxbound/nx
        dx=self.dxbound/2.
        dy=dx
        dz=2*np.sqrt(6)*deltax
        distribute.setbound(self,-dx,dx,-dy,dy,-dz,dz)
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
        totmass = self.dxbound*self.dybound*self.dzbound*self.rhozero
        self.mass = [totmass/self.npart]*self.npart
        if vm==1:
            print("Varying masses -> ONLY FOR RAND + RELAXING TO GLASS");              
            for i in range(0,self.npart,2):
                dmrat=0.25*self.mass[i]*np.random.rand()
                self.mass[i]=self.mass[i]+dmrat
                self.mass[i+1]=self.mass[i+1]-dmrat
        self.h=smth.getsmooth(self,200)
        print('npart = ',self.npart,' particle mass = ',self.mass[0],' total mass = ',totmass)
  
        for i in range(self.npart):
            vxi,vyi,vzi = self.getveli(i)
            self.vx.append(vxi)
            self.vy.append(vyi)
            self.vz.append(vzi)
            self.Bx.append(0.)
            self.By.append(0.)
            self.Bz.append(0.)
            self.u.append(self.getpressi(i)/(gam1*self.getrhoi(i)))
            self.rho[i]=self.getrhoi(i)

    def getrho(self,x):
        return self.rhozero

    def getveli(self,i):
        rcyl=np.sqrt(self.x[i]*self.x[i]+self.y[i]*self.y[i])
        if(rcyl<=0.2):
            vazi=5*rcyl
        elif(rcyl<=0.4):
            vazi=2.0-5.0*rcyl
        else:
            vazi=0.0
        vx,vy,vz = self.transform_coordinate(0.0,vazi,0.0,self.x[i],self.y[i],self.z[i],"cyl","cart")
        return vx,vy,vz
    
    def getpressi(self,i):
        rcyl=np.sqrt(self.x[i]*self.x[i]+self.y[i]*self.y[i])
        if(rcyl<=0.2):
            P=self.przero + 12.5*rcyl**2
        elif(rcyl<=0.4):
            P=self.przero+4.0+12.5*rcyl**2 - 20.0*rcyl + 4*np.log(5*rcyl)
        else:
            P=self.przero-2.0+4*np.log(2)
        return P
    
    def getrhoi(self,i):
        return self.rhozero
    
    def transform_coordinate(self,v1,v2,v3,x,y,z,fromcord,tocord):
        
        if(fromcord=="cyl"):
            r=np.sqrt(x**2+y**2)
            if(tocord=="cart"):
                ux=v1*x/r-v2*y/r
                uy=v1*y/r+v2*x/r 
                uz=v3
                
        return ux,uy,uz
