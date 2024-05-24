## Setup for ORZAG-TANG
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
from numba import njit
import readtipsy as tip

class setup_mhdblast(object):
    dICdensRsmooth = 0.1
    dICdensprofile = 1
    dICdensdir = 1
    dICdensR = 1.0 
    dICdensinner = 1.0 # calculated depending on mass
    dICdensouter = 1.0
    npart=0
    rhoit=0
    ns=64
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
    freqout=20
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
         Pout=1
         betaout=0.02
         Pin=100
         betain=2
         Rin=0.125
         Bzero=np.sqrt(2*Pin/betain)
         ublast=1
         #boundaries    
         dz=0.5
         hfact=5
         dy=0.5
         dx=0.5
         distribute.setbound(self,-dx,dx,-dy,dy,-dz,dz)
         deltax = self.dxbound/nx
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
         print('npart = ',self.npart,' particle mass = ',self.mass[1])
         self.h=smth.getsmooth(self,200)
         print(self.h/deltax)
         for i in range(self.npart):
             radius2 = self.x[i]**2 + self.y[i]**2+self.z[i]**2
             self.vx.append(0.)
             self.vy.append(0.)
             self.vz.append(0.)
             if (radius2 < Rin**2):
                self.u.append(Pin/(gam1*self.rhozero))
             else:
                self.u.append(Pout/(gam1*self.rhozero))
             self.Bx.append(Bzero*sqrt(0.5))
             self.By.append(0.)
             self.Bz.append(Bzero*sqrt(0.5))
             #toten=toten + self.mass[i]*self.u[i]
         #self.u=[x * (ublast/toten) for x in self.u]
         
    def getrhoi(self,i):
        return self.rhozero

 
