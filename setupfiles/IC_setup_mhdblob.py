## Setup for ORZAG-TANG
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
class setup_mhdblob(object):
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
    shape = 0
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
    rhoin = 10.0
    rhoout = 1.0
    Rin = 0.1
    
    def __init__(self):
        pass
    def create(self,nx,Bzero,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000'):
         gam1=self.gamma-1
         przero = 1.0
         denszero  = 1.0
         denscloud = 10.0
         velmachno = 2.7
         spsound = np.sqrt(self.gamma*przero/self.rhoout)
         vzero = velmachno*spsound
         tcrush = 2.*self.Rin*np.sqrt(self.rhoin/self.rhoout)/vzero
         taukh = 1.6*tcrush
         print("tcrush: ", tcrush)
         print("taukh: ", taukh)

         #boundaries
         zmin=-5*self.Rin
         ymin=-5*self.Rin
         xmin=-5*self.Rin
         xmax=35*self.Rin
         distribute.setbound(self,xmin,xmax,ymin,-ymin,zmin,-zmin)
         deltax = self.dxbound/nx
         deltasphere = deltax*(self.rhoout/self.rhoin)**(1./3.)
         self.rmin2=self.Rin**2
         distribute.setcloseddist(self,xmin,xmax,ymin,-ymin,zmin,-zmin,deltax)
         self.rmin2=0
         self.rmax2=self.Rin**2
         distribute.setcloseddist(self,xmin,xmax,ymin,-ymin,zmin,-zmin,deltasphere)
         self.npart=len(self.x)
         self.ngas=self.npart
         totmass = (self.rhoout*self.dxbound*self.dybound*self.dzbound + (self.rhoin-self.rhoout)*4/3*np.pi*self.Rin**3)
         self.mass = [totmass/self.npart]*self.npart
         print('npart = ',self.npart,' particle mass = ',self.mass[1], ' totmass = ', totmass, ' sphere mass in solar mass = ', self.rhoin*4/3*np.pi*self.Rin**3)
         
         for i in range(self.npart):
             radius2 = self.x[i]**2 + self.y[i]**2+self.z[i]**2
             #velvec in lengthunit/s
             if (radius2 < self.Rin**2):
                #self.u.append(self.getP(self.rhoin)/(gam1*self.rhoin))
                self.vx.append(0.0)
                self.u.append(przero/(gam1*self.rhoin))
             else:
                #self.u.append(self.getP(self.rhoout)/(gam1*self.rhoout))
                self.vx.append(vzero)
                self.u.append(przero/(gam1*self.rhoout))
             self.vy.append(0.0)
             self.vz.append(0.0)
             self.Bx.append(0.)
             self.By.append(Bzero)
             self.Bz.append(0.)
             self.x[i]=self.x[i]-(xmax+xmin)/2
             
         self.h=smth.getsmooth(self,200)
         
    def getrhoi(self,i):
        radius2=self.x[i]**2 + self.y[i]**2+self.z[i]**2
        if (radius2 < self.Rin**2):
            return self.rhoin
        else:
            return self.rhoout