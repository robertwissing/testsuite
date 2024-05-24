## Setup for ORZAG-TANG
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
class setup_mhdbalsaravortex(object):
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
    deltastep=0.1
    nsteps=200
    dmsolunit=1.0
    dkpcunit=1.0
    grav=0
    cosmo=0
    molweight=1
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
    rhodisk=2.0
    rhomean=1.5
    rloop=0.3
    smooth=0
    
    def __init__(self):
        pass
    def create(self,nx,u=1.0,k=1.0,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000'):
        
         #--setup parameters
         vzero = np.sqrt(2.)
         przero = 1.0
         self.rhozero = 1.0
         C=1.0/(2*np.pi)
         gam1=self.gamma-1
         uzero=przero/(gam1*self.rhozero)
         print('Setup for MHD balsara vortex problem')
         self.dxbound=20.0
         deltax = self.dxbound/nx
         dx=self.dxbound/2.
         dy=dx
         dz=2*np.sqrt(6)*deltax
         distribute.setbound(self,-dx,dx,-dy,dy,-dz,dz)
         costheta = self.dxbound/np.sqrt(self.dxbound**2 + self.dybound**2)
         sintheta = self.dybound/np.sqrt(self.dxbound**2 + self.dybound**2)
         deltadisk = deltax*(self.rhozero/self.rhodisk)**(1./3.)
         distribute.setcloseddist(self,-dx,dx,-dy,dy,-dz,dz,deltax)
         totmass = self.rhozero*self.dxbound*self.dybound*self.dzbound
         self.npart=len(self.x)
         self.ngas=self.npart 
         self.mass = [totmass/self.npart]*self.npart
         print('npart = ',self.npart)
         self.h=smth.getsmooth(self,200)
    
         for i in range(self.npart):            
             radius2 = self.x[i]**2 + self.y[i]**2
             rsq1 = 1.0 - radius2
             self.Bx.append(-self.y[i]*C*u*np.exp(0.5*rsq1))
             self.By.append(self.x[i]*C*u*np.exp(0.5*rsq1))
             self.Bz.append(0.)
             self.vx.append(1-self.y[i]*C*k*np.exp(0.5*rsq1))
             self.vy.append(1+self.x[i]*C*k*np.exp(0.5*rsq1))
             self.vz.append(1.0)
             self.u.append(1.5*(1.0+np.exp(rsq1)/(32*np.pi**3)*(u*u*rsq1-0.5*k*k)))

             
    def getrhoi(self,i):
        radius2 = self.x[i]**2 + self.y[i]**2
        if (radius2 < self.rloop**2):
            return self.rhodisk
        else:
            return self.rhozero
        
    def getrho(self,r):
        radius2 = r**2
        if (radius2 < self.rloop**2):
            return self.rhodisk
        else:
            return self.rhozero