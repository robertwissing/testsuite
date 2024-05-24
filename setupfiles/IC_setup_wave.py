#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 21:30:46 2018

@author: robertwi
"""

## Setup for ORZAG-TANG
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
class setup_wave(object):
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
    deltastep=0.1
    grav=0
    nsteps=100
    dmsolunit=1.0
    dkpcunit=1.0
    adi=1
    shape=2
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
    
    def __init__(self):
        pass
    def create(self,nx,rotate):
        if(rotate==1):           
            sina = 2./3.
            sinb = 2./np.sqrt(5.)
            cosa = np.sqrt(5.)/3.
            cosb = 1./np.sqrt(5.)
        else:
            sina = 0.
            sinb = 0.
            cosa = 1.
            cosb = 1.
        
        runit = (cosa*cosb,cosa*sinb,sina)
        self.drho = 10**(-3)
        #self.rhozero  = 1.-self.drho-0.000005
        self.rhozero = 1.
        cs = 1.
        dx = 0.5
        wavelength=2*dx
        deltax = wavelength/nx
        # try to give y boundary that is a multiple of 6 particle spacings in the low density part
        fac = 12.*(int((1.-np.finfo(float).eps)*2./6.) + 1)
        dy = fac*deltax*np.sqrt(0.75)
        dz = fac*deltax*np.sqrt(6.)/3.
        distribute.setbound(self,-dx,dx,-dy,dy,-dz,dz)
        self.wk = 2.*np.pi/wavelength
        self.x0 = (-dx,-dy,-dz)
        self.xdot0 = np.dot(self.x0,runit)
        distribute.setcloseddist(self,-dx,dx,-dy,dy,-dz,dz,deltax)
        if self.drho != 0:
            fixdens.set_density_profile(self,-dx,dx,1)
        
        #uuzero = 3./2.*cs**2
        #przero = cs**2*self.rhozero
        
        gam1 = self.gamma - 1.
        uuzero = cs**2/(self.gamma*gam1)
        przero = gam1*self.rhozero*uuzero
        
        self.npart=len(self.x)
        self.ngas=self.npart
        totmass = self.rhozero*self.dxbound*self.dybound*self.dzbound
        self.mass = [totmass/self.npart]*self.npart
        self.h=smth.getsmooth(self,200)
        
        print('npart = ',self.npart,' particle mass = ',self.mass[0], 'wk',self.wk,'x0',self.xdot0)
  
        for i in range(self.npart):
            xyz=(self.x[i],self.y[i],self.z[i])
            x1    = np.dot(xyz,runit)
            sinx1 = np.sin(self.wk*(x1-self.xdot0))
            vvec = self.drho*np.array([sinx1,0.,0.])
            vxyz=self.transform_vec(vvec,sina,sinb,cosa,cosb)
            self.vx.append(vxyz[0])
            self.vy.append(vxyz[1])
            self.vz.append(vxyz[2])
            self.Bx.append(0.)
            self.By.append(0.)
            self.Bz.append(0.)
            uui=uuzero + przero/self.rhozero*self.drho*sinx1
            self.u.append(uui)
             
    def getrho(self,x):
        return (self.rhozero + self.drho*np.sin(self.wk*(x - self.xdot0)))
    def getrhoi(self,i):
        return (self.rhozero + self.drho*np.sin(self.wk*(self.x[i] - self.xdot0))) 
    def getamplitudes():
        return
    
    #------------------------------------------------------
    #+
    #  transform vectors from rotated to unrotated coords
    #+
    #------------------------------------------------------
    def transform_vec(self,xvec,sina,sinb,cosa,cosb):
         x = [0,0,0]
         x[0] = xvec[0]*cosa*cosb - xvec[1]*sinb - xvec[2]*sina*cosb
         x[1] = xvec[0]*cosa*sinb + xvec[1]*cosb - xvec[2]*sina*sinb
         x[2] = xvec[0]*sina + xvec[2]*cosa
         return x