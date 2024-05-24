## Setup for ORZAG-TANG
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
class setup_mhdloopshear(object):
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
    rloop=0.3
    y0=0.0
    
    def __init__(self):
        pass
    def create(self,nx,rhodiff,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000'):
        
        
         #--setup parameters
         cs = 0.6
         q=3/2
         orbvel=1.0
         przero = 1.0
         self.rhozero = 1.0
         self.rhodisk = self.rhozero*rhodiff
         Azero = 1.e-3
         self.rloop = 0.3
         gam1=self.gamma-1
         uzero=przero/(gam1*self.rhozero)
         print('Setup for MHD current loop advection problem')
         print(' Azero        = ',Azero,', rotation speed = ',cs,
                    ' radius of loop   = ',self.rloop,', density = ',self.rhozero,', pressure = ',przero)
    
    
    
         #boundaries    
         #dz = 6*deltax
 #        deltax = 1.0/nx
 #        deltadisk = deltax*(rhozero/rhodisk)**(1./3.)
 #        nxdisk=(2*rdisk)/deltadisk+(1.0-2*rdisk)
         dz=4*np.sqrt(6)/(nx)
         ymin=-2
         xmin=-0.5
         self.y0=1
         distribute.setbound(self,xmin,-xmin,ymin,-ymin,-dz,dz)
         deltax = self.dxbound/nx
         costheta = self.dxbound/np.sqrt(self.dxbound**2 + self.dybound**2)
         sintheta = self.dybound/np.sqrt(self.dxbound**2 + self.dybound**2)
         deltadisk = deltax*(self.rhozero/self.rhodisk)**(1./3.)
         self.rcylmin2 = self.rloop**2
         distribute.setcloseddist(self,xmin,-xmin,ymin,-ymin,-dz,dz,deltax)
         self.rcylmin2 = 0.0
         self.rcylmax2 = self.rloop**2
         distribute.setcloseddist(self,xmin,-xmin,ymin,-ymin,-dz,dz,deltadisk)
         self.npart=len(self.x)
         self.ngas=self.npart    
         totmass = (self.rhozero*self.dxbound*self.dybound + self.rhodisk*np.pi*self.rloop**2)*self.dzbound
         self.mass = [totmass/self.npart]*self.npart
         print('npart = ',self.npart)
         self.h=smth.getsmooth(self,200)
         print(self.h/deltax)
         
         for i in range(self.npart):            
             radius2 = self.x[i]**2 + (self.y[i]-self.y0)**2
             if (radius2 <= self.rloop**2):
                 self.Bx.append(Azero*(-(self.y[i]-self.y0)/np.sqrt(radius2)))
                 self.By.append(Azero*(self.x[i]/np.sqrt(radius2)))
                 self.u.append(przero/(gam1*self.rhozero))
             else:
                 self.Bx.append(0.)
                 self.By.append(0.)
                 self.u.append(przero/(gam1*self.rhodisk))
             self.Bz.append(0.)
             self.vx.append(cs)
             self.vy.append(-q*orbvel*self.x[i])
             self.vz.append(0.)

             
    def getrhoi(self,i):
        radius2 = self.x[i]**2 + (self.y[i]-self.y0)**2
        if (radius2 <= self.rloop**2):
            return self.rhodisk
        else:
            return self.rhozero