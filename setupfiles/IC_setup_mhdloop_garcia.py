IC_setup_mhdloop_garcia.py## Setup for ORZAG-TANG
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
class setup_mhdloop_garcia(object):
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
    cosmo=0
    grav=0
    molweight=1
    gamma=5./3.
    periodic=1
    deltastep=0.1
    nsteps=400
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
    rhomean=1.5
    rloop=0.3
    smooth=0
    
    def __init__(self):
        pass
    def create(self,nx,rhodiff=1.0,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000'):
        
        
         #--setup parameters
         vzero = 2.0

         przero = 1.0
         self.rhozero = 1.0
         self.rhodisk = self.rhozero*rhodiff
         Azero = 1.e-3
         self.rloop = 0.3
         gam1=self.gamma-1
         uzero=przero/(gam1*self.rhozero)
         print('Setup for MHD current loop advection problem')
         print(' Azero        = ',Azero,', rotation speed = ',vzero,
                    ' radius of loop   = ',self.rloop,', density = ',self.rhozero,', pressure = ',przero)
         self.dxbound=2.0
         self.rhomean=(self.rhozero+self.rhodisk)/2.
         deltax = self.dxbound/nx
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
             #distribute2.setcloseddist(self,-dx,dx,-dy,dy,-dz,dz,deltadisk,np.sqrt(rhodiff))
             totmass = (self.rhozero*self.dxbound*self.dybound + (self.rhodisk-self.rhozero)*np.pi*self.rloop**2)*self.dzbound
         else:
             #self.rcylmax2 = 1.0
             distribute.setcloseddist(self,-dx,dx,-dy,dy,-dz,dz,deltax)
             fixdens.set_density_profile(self,0.0,rbig,rhofunc=True,icoord=0,igeom=2)
             totmass = self.rhomean*self.dxbound*self.dybound*self.dzbound
         self.npart=len(self.x)
         self.ngas=self.npart 
         self.mass = [totmass/self.npart]*self.npart
         print('npart = ',self.npart)
         self.h=smth.getsmooth(self,200)
    
         for i in range(self.npart):            
             radius2 = self.x[i]**2 + self.y[i]**2 + self.z[i]**2
             rcyl2 = self.x[i]**2 + self.y[i]**2
             if (radius2 < self.rloop**2):
                 self.Bx.append(Azero*(-self.y[i]/np.sqrt(rcyl2)))
                 self.By.append(Azero*(self.x[i]/np.sqrt(rcyl2)))
                 self.u.append(przero/(gam1*self.rhodisk))
                 self.Bz.append(Azero*0.0)
             else:
                 self.Bx.append(0.)
                 self.By.append(0.)
                 self.u.append(przero/(gam1*self.rhozero))
                 self.Bz.append(0.)
             self.vx.append(vzero*costheta*0.0)
             self.vy.append(vzero*sintheta*0.0)
             self.vz.append(1.0*vzero)

             
    def getrhoi(self,i):
        radius2 = self.x[i]**2 + self.y[i]**2 + self.z[i]**2
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