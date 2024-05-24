## Setup for ORZAG-TANG
import numpy as np
import random 
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
class setup_mri(object):
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
    gamma=1.00001
    shape=2
    periodic=1
    grav=1
    deltastep=1.0
    nsteps=200
    dmsolunit=1.0
    dkpcunit=1.0
    adi=0
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
    pot=[]
    dxbound=[]
    dybound=[]
    dzbound=[]
    totvol=[]
    rcylmin2= 0
    rmin2= 0
    rcylmax2= 10**22
    rmax2= 10**22
    rhozero=1.0
    stratified = 1
    adiabatic = 0
    Binvert = 0
    netflux = 0
    H = 1.0
    glass = 1
    def __init__(self):
        pass
    def create(self,nx,case,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000'):
         #--setup parameters
         q=3/2
         orbvel=1.0
         m=4
         #B0=np.sqrt(15)/(8*np.pi*m)
         self.rhozero = 1.0
         #boundaries    
         #dz = 6*deltax
 #        deltax = 1.0/nx
 #        deltadisk = deltax*(rhozero/rhodisk)**(1./3.)
 #        nxdisk=(2*rdisk)/deltadisk+(1.0-2*rdisk)
         self.H=1.0
         xmin=-self.H/2
         if case==1: # thin unstratified nonetflux
             dz=2*np.sqrt(6)/(nx)
             ymin=xmin
             self.stratified = 0
             self.Binazi = 0
             self.netflux = 0
             beta=400
         elif case==2: # regular unstratified nonetflux
             dz=-xmin
             ymin=-np.pi*dz
             #ymin=-dz
             self.stratified = 0
             self.Binazi = 0
             self.netflux = 0
             beta=400
         elif case==3: # tall unstratified nonetflux
             dz=-4.0*xmin
             ymin=-dz
             self.stratified = 0
             self.Binazi = 0
             self.netflux = 0
             beta=400
         elif case==4: #stratified netflux
             dz=(self.H/2)*4
             xmin=-np.sqrt(2)*(self.H/2)
             ymin=4.0*xmin
             self.stratified = 1
             self.Binazi = 1
             self.netflux = 1
             beta=400
         elif case==5: #unstratified netflux
             dz=-xmin
             ymin=6.28*xmin
             self.stratified = 0
             self.Binazi = 0
             self.netflux = 1
             beta=400
         elif case==6: # regular unstratified nonetflux
             dz=-xmin
             ymin=-np.pi*dz
             self.stratified = 0
             self.Binazi = 1
             self.netflux = 0
             beta=400
         elif case==7:
             dz=-xmin
             ymin=-dz
             self.stratified = 0
             self.Binazi = 0
             self.netflux = 0
             beta=400
         elif case==8: #unstratified netflux tall
             dz=-4.0*xmin
             ymin=-dz
             self.stratified = 0
             self.Binazi = 0
             self.netflux = 1
             beta=400
             
         if self.adiabatic==1:
                 cs=np.sqrt(self.gamma)
                 gam1=self.gamma-1
                 przero=self.rhozero
                 uzero=przero/(gam1*self.rhozero)
                 B0=np.sqrt(2*przero/beta)
         else:
             cs=1.0
             gam1=self.gamma-1
             przero=1.0
             uzero=przero/(gam1*self.rhozero)
             print("UZERO",uzero)
             B0=np.sqrt(2*przero/beta)
             
         distribute.setbound(self,xmin,-xmin,ymin,-ymin,-dz,dz)
         deltax = self.dxbound/nx
         distribute.setcloseddist(self,xmin,-xmin,ymin,-ymin,-dz,dz,deltax)
         self.npart=len(self.x)
         totmass = self.rhozero*self.dxbound*self.dybound*self.dzbound
         print('npart = ',self.npart)
         print("totmass ",totmass)
         self.mass = [totmass/self.npart]*self.npart
         print('npart = ',self.npart, 'mass = ', self.mass[1])
         if self.stratified==1:
             fixdens.set_density_profile(self,-dz,dz,1,icoord=2)
         
         self.npart=len(self.x)
         self.ngas=self.npart
         totmass = self.rhozero*self.dxbound*self.dybound*self.dzbound
         self.mass = [totmass/self.npart]*self.npart
         print('npart = ',self.npart, 'mass = ', self.mass[1])
         self.h=smth.getsmooth(self,200)
         print(self.h/deltax)

         for i in range(self.npart):
             if self.glass == 1:
                 self.Bx.append(0.)
                 self.By.append(0.)
                 self.Bz.append(0.)
                 rand= random.uniform(-0.2*cs,0.2*cs)
                 rand2= random.uniform(-0.2*cs,0.2*cs)
                 rand3= random.uniform(-0.2*cs,0.2*cs)
                 self.vx.append(rand)
                 self.vy.append(rand2)
                 #self.vy.append(-q*orbvel*self.x[i])
                 self.vz.append(rand3)
                 self.u.append(uzero)
             else:
                 if self.netflux == 1:
                     Bterm=B0
                 else:                    
                     Bterm=B0*np.sin(2*np.pi*self.x[i]/self.H)
                 self.Bx.append(0.)
                 if self.Binazi == 1:
                     self.By.append(Bterm)
                     self.Bz.append(0.) 
                 else:                                      
                     self.By.append(0.)
                     self.Bz.append(Bterm)
                 if self.stratified==1:
                     self.u.append(cs/(gam1))
                 else:
                     self.u.append(uzero)
                 rand= random.uniform(-5.3,5.3)
                 rand2= random.uniform(-5.3,5.3)
                 self.vx.append(rand)
                 self.vy.append(-q*orbvel*self.x[i])
                 self.vz.append(rand2)

             
    def getrhoi(self,i):
        if(self.stratified==1):
            return self.rhozero*np.exp(-self.z[i]**2/(2*self.H**2))
        else:
            return self.rhozero
        
    def getrho(self,x):
        if(self.stratified==1):
            return self.rhozero*np.exp(-x**2/(2*self.H**2))
        else:
            return self.rhozero