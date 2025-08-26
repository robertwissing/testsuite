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
    kmag = 1.0
    def __init__(self):
        pass
    def create(self,nx,distri=0,vm=0.0,entry='datafiles/alfvenwave128_preglass.00000',mach=np.sqrt(3/25),beta=10e6):
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
            self.mass=tgdata[:,0]
        self.npart=len(self.x)
        self.ngas=self.npart
        totmass = self.dxbound*self.dybound*self.dzbound*self.rhozero
        if distri != 2:
            self.mass = [totmass/self.npart]*self.npart
            
            if vm > 1.0 or vm < -1.0:
                print("Varying masses -> ONLY FOR RAND + RELAXING TO GLASS");              
                for i in range(0,self.npart,2):
                    #max mass - min mass * random
                    #dmrat=(vm*self.mass[i]-(1.0/vm)*self.mass[i])*np.random.rand()
                    total_mass = self.mass[i] + self.mass[i+1]
                    # Calculate allowable mass range for each particle to respect max_ratio
                    min_mass = total_mass / (1 + np.abs(vm))  # Minimum allowed mass for either particle
                    max_mass = total_mass - min_mass         # Maximum allowed mass 
                    if (vm > 1.0):
                        # Randomly redistribute mass within [min_mass, max_mass]
                        m1 = np.random.uniform(min_mass, max_mass)
                    else:
                        m1 = max_mass
                    m2 = total_mass - m1

                    self.mass[i]=m1
                    self.mass[i+1]=m2
        self.h=smth.getsmooth(self,200)
        print('npart = ',self.npart,' particle mass = ',self.mass[0],' total mass = ',totmass)
        self.kmag = np.sqrt(self.getpress(0.2)/( (1+beta)*0.5*self.getBform(0.2)**2 ))
        
        for i in range(self.npart):
            vxi,vyi,vzi = self.getveli(i)
            self.vx.append(vxi)
            self.vy.append(vyi)
            self.vz.append(vzi)
            Bxi,Byi,Bzi = self.getBi(i)
            Ptot = self.getpressi(i)
            Pmag = (Bxi**2+Byi**2+Bzi**2)*0.5
            Pth = Ptot-Pmag
            if Pth < 0.0:
                print(Pth,rcyl)
            self.Bx.append(Bxi)
            self.By.append(Byi)
            self.Bz.append(Bzi)
            self.u.append(Pth/(gam1*self.getrhoi(i)))
            self.rho[i]=self.getrhoi(i)

    def getrho(self,x):
        return self.rhozero

    def getBform(self,rcyl):
        select=3
        if(select==1):
            Bform =np.exp(-5*(rcyl-0.2)**2)*np.sin(2.5*np.pi*rcyl)
            if (rcyl>=0.4):
                Bform = 0.0
        if(select==2):
            Bform = rcyl**2*(rcyl-0.4)**2
            if (rcyl>=0.4):
                Bform = 0.0
        if(select==3):
            if(rcyl<=0.2):
                Bform=5*rcyl
            elif(rcyl<=0.4):
                Bform=2.0-5.0*rcyl
            else:
                Bform=0.0
        return Bform; 

    def getBi(self,i):
        rcyl=np.sqrt(self.x[i]*self.x[i]+self.y[i]*self.y[i])
        Bazi = self.kmag*self.getBform(rcyl)
        Bx,By,Bz = self.transform_coordinate(0.0,Bazi,0.0,self.x[i],self.y[i],self.z[i],"cyl","cart")

        return Bx,By,Bz

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

    def getpress(self,rcyl):
        if(rcyl<=0.2):
            P=self.przero + 12.5*rcyl**2
        elif(rcyl<=0.4):
            P=self.przero+4.0+12.5*rcyl**2 - 20.0*rcyl + 4*np.log(5*rcyl)
        else:
            P=self.przero-2.0+4*np.log(2)
        return P

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
