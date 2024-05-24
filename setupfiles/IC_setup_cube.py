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
import readtipsy as tip
class setup_cube(object):
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
    deltastep=0.03
    nsteps=400
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
    dxbound=[]
    dybound=[]
    dzbound=[]
    totvol=[]
    rcylmin2= 0
    rmin2= 0
    rcylmax2= 10**22
    rmax2= 10**22
    rhozero=0.
    shape=1
      
    def __init__(self):
        pass
    
    def create(self,nx,displace=0,vzero=0.0,rotated=0,zvel=0.0,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000'):        
         #--setup parameters
         przero = 1.0
         self.rhozero = 1.0
         gam1=self.gamma-1
         uzero=przero/(gam1*self.rhozero)
    
    
         #boundaries    
         zmin=-0.5
         ymin=-0.5
         xmin=-0.5
         dz=-zmin
         dy=-ymin
         dx=-xmin
         distribute.setbound(self,xmin,-xmin,ymin,-ymin,zmin,-zmin)
         deltax = self.dxbound/nx
         if rotated==1:
             costheta = self.dxbound/np.sqrt(self.dxbound**2 + self.dybound**2)
             sintheta = self.dybound/np.sqrt(self.dxbound**2 + self.dybound**2)
         else:
             costheta = 1
             sintheta = 0
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
         deltaz=deltax*np.sqrt(6.)/3.
         print(self.h/deltax,deltaz)
         for i in range(self.npart):
             if(displace==1):   
                 self.displacepart(i,deltaz*1.9)
             self.Bx.append(0.)
             self.By.append(0.)
             self.Bz.append(0.)
             self.vx.append(vzero*costheta)
             self.vy.append(vzero*sintheta)
             self.vz.append(zvel*0.1*vzero)
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