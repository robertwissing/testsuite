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
class setup_cosmowave(object):
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
    H0=2.894405
    cosmo=1
    molweight=1.0
    gamma=5./3.
    periodic=1
    deltastep=0.0005
    grav=0
    nsteps=100
    dmsolunit=1.0
    dkpcunit=1000.0
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
    Bzero=1.
    dB=0.
    a=0.
    kappa=1.
    def __init__(self):
        pass
    def create(self,nx,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000',trav=1,rotate=0,Lboxinkpc=1,zi=1.0):
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
        self.drho = 0
        #self.rhozero  = 1.-self.drho-0.000005
        self.dkpcunit=Lboxinkpc
        ascale=1./(1.+zi)
        """
        #dimensionless amplitude Au=10**-6
        #wavenumber kL=2pi k=2pi
        #z=127 ai=1/128
        #du(ai)=dA*Va cos(kx) for alfven waves
        #du(ai)=dA*Vs cos(kx) for compressible
        #The travelling waves are initialized using the analytic solutions at
        #a = ai . The amplitudes shown in the figures have been divided by
        #the value of Au such that the velocity amplitude appears as 1 at a =
        #ai
        We use this solution to initialize a simulation with A = π and
        show the amplitude evolution of the wave in Fig. 1.
        #omega_A=pi
        #omega_A=2pi
        show δBc /Bc multiplied with a1/4 and δu/VA
        multiplied with a3/4 as a function of ln a/ai 
        
        omega_S/V_S=k/H_0
        omega_A/V_A=k/H_0
        
        V_A=B_c/sqrt(rho_c)
        V_S=sqrt(gamma*rho_c0/rho_c)
        V_G=sqrt(4piGrho_c)/k
        omega_x=kV_x/H_0
        K=sqrt(omega_A**2-1/16)
        phi=Kln(a/a_i) phi(a_i)=0   
        
        du(x,a_i)=V_A*A_u*cos(phi)*cos(kx)
        du(x,a_i)=H_0/k*omega_A*A_u*cos(kx)
        B_c(a_i)=(H_0/k)*omega_A*sqrt(rho_c)
        The critical density of the Universe in simulation units is 1.
        
        dTheta2          = 0.725
        dHubble0         = 2.894405
        dOmega0          = 1.0
        dLambda          = 0.0
        dRedTo   = 0.0
        bComove          = 1
        """
        dx = 0.5
        wavelength=2*dx
        k=2.*np.pi/wavelength
        deltax = wavelength/nx
        A_u = 10**(-1)

        omega_A = np.pi
        omega_S = np.pi
        omega_G = 0.0
        #omega_G = np.pi/2
        
        self.gamma=4./3.

        V_A=(self.H0/k)*omega_A
        V_S=(self.H0/k)*omega_S
        V_G=(self.H0/k)*omega_G
        print("V_a ",V_A,"V_S ",V_S,"V_G ",V_G)
        
        sigma=np.sqrt(omega_A**2+omega_S**2-omega_G**2)
        
        kappa=np.sqrt(sigma**2-1./16.)
        self.kappa=kappa
        if trav==0:
            """ Standing wave """
            a=1/128 
        else:
            """ Moving wave """
            n=3
            a=np.exp(-2*n*np.pi/kappa)
        self.a=a
        phi=kappa*np.log(1) #log(a/a_i)
        if omega_G > 0:
            self.rhozero=((k*V_G)**2)/(4*np.pi)
        print("rhozero", self.rhozero)
        
        valf=V_A*a**(-0.5)
        Bzero=V_A*np.sqrt(self.rhozero)*a**(-0.5)
        gam1 = self.gamma - 1.
        
        print("Bzero",Bzero)
        
        """ Incompressible """
        #uzero=V_A*A_u*a**(-1)
        #cs=uzero*10
        
        """ Compressible """
        uzero=V_S*A_u*a**(-1)
        cs=V_S*a**((-3*gam1)/2)
        print('cs',cs,'cs**2',cs**2,'Va',valf)
        print('uzero',uzero)
        
        uuzero = cs**2/(self.gamma*gam1)
        pzero=uuzero*self.rhozero*gam1
        print("uuzeroth",uuzero,"pzero",pzero)
        #uuzero=vA**2  B**2/(rho*self.gamma*gam1)
        dy=dx
        dz=dx
        distribute.setbound(self,-dx,dx,-dy,dy,-dz,dz)
        self.wk = k
        self.x0 = (-dx,-dy,-dz)
        self.xdot0 = np.dot(self.x0,runit)
        #distribute.setcloseddist(self,-dx,dx,-dy,dy,-dz,dz,deltax)
        distribute.setcubicdist(self,-dx,dx,-dy,dy,-dz,dz,deltax)
        #dbrho=A_u*omega_S*np.sqrt(a)/(4*sigma**2)*4*kappa
        self.drho= A_u*omega_S*np.sqrt(a)/(4*sigma**2)
        print("DRHO",self.drho,"wk",self.wk,"xdot0",self.xdot0,"kappa",kappa)
        if self.drho != 0:
            print("FIXING THE DENSITY HERE")
            fixdens.set_density_profile(self,-dx,dx,1)
        

        #Bzero=V_A*np.sqrt(self.rhozero)*a**(-0.5)
        #c_s=V_S*a**(-0.5)
        #P_th=rho*uzero*gam1
        #P_th=rho*cs**2/(gamma)
        #Pmag=B**2/2=rho*vA**2/2          P_th/P_mag = (cs**2/vA**2)*(2/gamma)
        

        self.npart=len(self.x)
        self.ngas=self.npart
        totmass = self.rhozero*self.dxbound*self.dybound*self.dzbound
        self.mass = [totmass/self.npart]*self.npart
        self.h=smth.getsmooth(self,200)
        
        print('npart = ',self.npart,' particle mass = ',self.mass[0], 'wk',self.wk,'x0',self.xdot0)
  
        for i in range(self.npart):
            xyz=(self.x[i],self.y[i],self.z[i])
            x1    = np.dot(xyz,runit)
            cosx1 = np.cos(self.wk*(x1-self.xdot0))
            #vvec = uzero*np.array([0.,cosx1,0.])
            vvec = uzero*np.array([cosx1,0.,0.]) #comp
            vxyz=self.transform_vec(vvec,sina,sinb,cosa,cosb)
            self.vx.append(vxyz[0])
            self.vy.append(vxyz[1])
            self.vz.append(vxyz[2])
            if trav==1:
                #self.Bz.append(Bzero+Bzero*self.drho*(4*self.kappa*np.cos(self.wk*(self.x[i] - self.xdot0))-np.sin(self.wk*(self.x[i] - self.xdot0))))
                #self.Bz.append(Bzero+Bzero*self.drho*(4*self.kappa*np.cos(self.wk*(self.x[i] - self.xdot0))-np.sin(self.wk*(self.x[i] - self.xdot0))))
                self.Bz.append(Bzero+Bzero*(self.rho[i]-self.rhozero)/self.rhozero)
                #self.Bz.append(Bzero)
                #self.Bz.append(Bzero-Bzero*A_u*np.sqrt(a)/(4*np.pi*omega_A)*(4*self.kappa*np.cos(self.wk*(self.x[i] - self.xdot0))-np.sin(self.wk*(self.x[i] - self.xdot0)))) 
            else:
                self.Bz.append(Bzero)
            self.By.append(0.)
            self.Bx.append(0.)
            #print(pzero, cs**2,(self.rho[i]-self.rhozero),gam1,"rho",self.rho[i])
            uui=(pzero+cs**2*(self.rho[i]-self.rhozero))/(gam1*self.rho[i])
            self.u.append(uui)

             
    def getrho(self,x):
        return (self.rhozero + self.drho*(4*self.kappa*np.cos(self.wk*(x - self.xdot0))-np.sin(self.wk*(x - self.xdot0))))

    def getrhoi(self,i):
        return (self.rhozero + self.drho*(4*self.kappa*np.cos(self.wk*(self.x[i] - self.xdot0))-np.sin(self.wk*(self.x[i] - self.xdot0))))

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
