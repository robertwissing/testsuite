## Setup for ORZAG-TANG
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
class setup_mhdblast_garcia(object):
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
    gamma=1.666666667
    periodic=1
    deltastep=0.0002
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
    dxbound=[]
    dybound=[]
    dzbound=[]
    totvol=[]
    rcylmin2= 0
    rmin2= 0
    rcylmax2= 10**22
    rmax2= 10**22
    rhozero= 1.0
    
    def __init__(self):
        pass
    def create(self,nx,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000'):
        
         gam1=self.gamma-1
    
         #--setup parameters
         self.rhozero=1.0
         energytot=5.0
         width=0.1
         Bzero=10
         ener0=energytot/np.pi**(3.0/2.0)/1.0/width**3
         #boundaries
         dz=1.0
         hfact=5
         ymin=-1.0
         xmin=-1.0
         distribute.setbound(self,xmin,-xmin,ymin,-ymin,-dz,dz)
         deltax = self.dxbound/nx
         distribute.setcloseddist(self,xmin,-xmin,ymin,-ymin,-dz,dz,deltax)
         self.npart=len(self.x)
         self.ngas=self.npart    
         totmass = self.rhozero*self.dxbound*self.dybound*self.dzbound
         self.mass = [totmass/self.npart]*self.npart
         print('npart = ',self.npart,' particle mass = ',self.mass[1])
         toten=0.
         self.h=smth.getsmooth(self,200)
         print(self.h/deltax)
         for i in range(self.npart):
             rrr = np.sqrt(self.x[i]**2 + self.y[i]**2+self.z[i]**2)
             self.vx.append(0.)
             self.vy.append(0.)
             self.vz.append(0.)
             self.u.append(ener0*np.exp(-(rrr**2/width**2))+1.0)
             self.Bx.append(0.)
             self.By.append(Bzero)
             self.Bz.append(0.)
    def getrhoi(self,i):
        return self.rhozero

 