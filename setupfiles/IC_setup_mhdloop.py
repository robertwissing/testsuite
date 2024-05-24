## Setup for ORZAG-TANG
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
import readtipsy as tip
class setup_mhdloop(object):
    dICdensRsmooth = 0.003
    dICdensprofile = 2
    dICdensdir = 5
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
    freqout=10
    ns=64
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
    rloop=0.3
    smooth=0
    
    def __init__(self):
        pass
    def create(self,nx,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000',rhodiff=1.0):
                
         #--setup parameters
         vzero = np.sqrt(5.)
         przero = 1.0
         self.rhozero = 1.0
         self.rhodisk = self.rhozero*rhodiff
         self.dICdensinner = self.rhodisk
         self.dICdensouter = self.rhozero
         Azero = 1.e-3
         self.rloop = 0.3
         self.dICdensR = 0.3
         gam1=self.gamma-1
         uzero=przero/(gam1*self.rhozero)
         print('Setup for MHD current loop advection problem')
         print(' Azero        = ',Azero,', rotation speed = ',vzero,
                    ' radius of loop   = ',self.rloop,', density = ',self.rhozero,', pressure = ',przero)
         self.dxbound=2.0
         deltax = self.dxbound/nx
         dx=self.dxbound/2.
         dy=dx/2.0
         dz=2*np.sqrt(6)*deltax
         distribute.setbound(self,-dx,dx,-dy,dy,-dz,dz)
         costheta = self.dxbound/np.sqrt(self.dxbound**2 + self.dybound**2)
         sintheta = self.dybound/np.sqrt(self.dxbound**2 + self.dybound**2)
         deltadisk = deltax*(self.rhozero/self.rhodisk)**(1./3.)

         if distri==0:
             self.rcylmin2 = self.rloop**2
             distribute.setcloseddist(self,-dx,dx,-dy,dy,-dz,dz,deltax)
             self.rcylmin2 = 0.0
             self.rcylmax2 = self.rloop**2
             distribute.setcloseddist(self,-dx,dx,-dy,dy,-dz,dz,deltadisk)
         elif distri==1:
             self.rcylmin2 = self.rloop**2
             distribute.setrandomdist(self,-dx,dx,-dy,dy,-dz,dz,deltax)
             self.rcylmin2 = 0.0
             self.rcylmax2 = self.rloop**2
             distribute.setrandomdist(self,-dx,dx,-dy,dy,-dz,dz,deltadisk)
         else:
            tgdata,tddata,tsdata,data_header,time=tip.readtipsy(entry);
            self.x=tgdata[:,1]
            self.y=tgdata[:,2]
            self.z=tgdata[:,3]
            self.rho=tgdata[:,7]


         
         totmass = (self.rhozero*self.dxbound*self.dybound + (self.rhodisk-self.rhozero)*np.pi*self.rloop**2)*self.dzbound
         self.npart=len(self.x)
         self.ngas=self.npart 
         self.mass = [totmass/self.npart]*self.npart

         if vm==1:
            print("Varying masses -> ONLY FOR RAND + RELAXING TO GLASS");
            for i in range(0,self.npart,2):
                dmrat=0.25*self.mass[i]*np.random.rand()
                self.mass[i]=self.mass[i]+dmrat
                self.mass[i+1]=self.mass[i+1]-dmrat
         
         print('npart = ',self.npart)
         self.h=smth.getsmooth(self,200)
    
         for i in range(self.npart):            
             radius2 = self.x[i]**2 + self.y[i]**2
             if (radius2 < self.rloop**2):
                 self.Bx.append(Azero*(-self.y[i]/np.sqrt(radius2)))
                 self.By.append(Azero*(self.x[i]/np.sqrt(radius2)))
                 self.u.append(przero/(gam1*self.rhodisk))
             else:
                 self.Bx.append(0.)
                 self.By.append(0.)
                 self.u.append(przero/(gam1*self.rhozero))
             self.Bz.append(0.)
             self.vx.append(vzero*costheta)
             self.vy.append(vzero*sintheta)
             self.vz.append(0.1*vzero)

             
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
