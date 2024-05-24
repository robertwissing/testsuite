## Setup for ORZAG-TANG
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
class setup_mhdpinch_garcia(object):
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
    gamma=1.66666666667
    periodic=1
    deltastep=0.002
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
         Bzero=3.0
         #boundaries
         dz=0.32
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
            self.vx.append(0.)
            self.vy.append(0.)
            self.vz.append(-1.0)
            xcen=self.x[i]
            ycen=self.y[i]
            xcena=abs(self.x[i])
            ycena=abs(self.y[i])
            rrr=np.sqrt(xcen**2+ycen**2)
            if(xcen > 0.0 and ycen > 0.0):
                phiaz=np.arccos(xcena/rrr)
            if(xcen < 0.0 and ycen > 0.0):
                phiaz=np.pi-np.arccos(xcena/rrr)
            if(xcen < 0.0 and ycen < 0.0):
                phiaz=np.pi+np.arccos(xcena/rrr)
            if(xcen > 0.0 and ycen < 0.0):
                phiaz=-np.arccos(xcena/rrr)
            rrr=np.exp(-(rrr-0.5)**2/0.01)
            self.Bx.append(-Bzero*np.sin(phiaz)*rrr)
            self.By.append(Bzero*np.cos(phiaz)*rrr)
            self.Bz.append(0.)
            self.u.append(1.0/gam1)
            if(rrr==0.0):
                self.Bx.append(0.)
                self.By.append(0.)
                self.Bz.append(0.)
             
             
             
         
        
    def getrhoi(self,i):
        return self.rhozero

 