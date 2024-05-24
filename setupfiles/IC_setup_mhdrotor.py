## Setup for ORZAG-TANG
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
import readtipsy as tip
class setup_mhdrotor(object):
    dICdensRsmooth = 0.001
    dICdensprofile = 3
    dICdensdir = 5
    dICdensR = 0.1
    dICdensinner = 1.0 # calculated depending on mass
    dICdensouter = 1.0
    rhoit=0
    npart=0
    ngas=0
    ndark=0
    nstar=0
    ndim=3
    time=0
    freqout=10
    ns=64
    gamma=1.4
    periodic=1
    deltastep=0.005
    nsteps=100
    dmsolunit=1.0
    dkpcunit=1.0
    molweight=1.0
    cosmo=0
    grav=0
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
    rhozero=0.
    rhodisk=0.
    rdisk=0.
    
    def __init__(self):
        pass
    def create(self,nx,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000',rhodisk=10.0):
        
        
         #--setup parameters
         const = np.sqrt(4.*np.pi)
         self.rdisk = 0.1           # radius of the initial disk
         vzero = 2.0             # rotation speed of initial disk
         Bzero = 5.0/const      # uniform field in bx direction
         przero = 1.0             # initial pressure
         self.rhozero = 1.0            # ambient density
         self.rhodisk = rhodisk          # density of rotating disk
         self.dICdensinner = self.rhodisk
         self.dICdensouter = self.rhozero
         self.dICdensRsmooth= self.rdisk*0.01
     
         print('Setup for MHD rotor problem')
         print(' radius of disk        = ',self.rdisk,', rotation speed = ',vzero,
                    ' initial B   = ',Bzero,', density = ',self.rhozero,', pressure = ',przero)
    
    
    
         dz=2*np.sqrt(6)/(nx)
         dy=0.5
         dx=0.5
         distribute.setbound(self,-dx,dx,-dy,dy,-dz,dz)
         deltax = self.dxbound/nx
         deltadisk = deltax*(self.rhozero/self.rhodisk)**(1./3.)

         if distri==0:
             self.rcylmin2=self.rdisk**2
             distribute.setcloseddist(self,-dx,dx,-dy,dy,-dz,dz,deltax)
             self.rcylmin2=0
             self.rcylmax2=self.rdisk**2
             distribute.setcloseddist(self,-dx,dx,-dy,dy,-dz,dz,deltadisk)
         elif distri==1:
            self.rcylmin2=self.rdisk**2
            distribute.setrandomdist(self,-dx,dx,-dy,dy,-dz,dz,deltax)
            self.rcylmin2=0
            self.rcylmax2=self.rdisk**2
            distribute.setcloseddist(self,-dx,dx,-dy,dy,-dz,dz,deltadisk)
         else:
            tgdata,tddata,tsdata,data_header,time=tip.readtipsy(entry);
            self.x=tgdata[:,1]
            self.y=tgdata[:,2]
            self.z=tgdata[:,3]
            self.rho=tgdata[:,7]
            
         self.npart=len(self.x)
         self.ngas=self.npart    
         totmass = (self.rhozero*self.dxbound*self.dybound + self.rhodisk*np.pi*self.rdisk**2)*self.dzbound
         self.mass = [totmass/self.npart]*self.npart
         if vm==1:
            print("Varying masses -> ONLY FOR RAND + RELAXING TO GLASS");
            for i in range(0,self.npart,2):
                dmrat=0.25*self.mass[i]*np.random.rand()
                self.mass[i]=self.mass[i]+dmrat
                self.mass[i+1]=self.mass[i+1]-dmrat
         print('npart = ',self.npart)
         self.h=smth.getsmooth(self,200)
         print(self.h/deltax)
    
         for i in range(self.npart):            
             radius2 = self.x[i]**2 + self.y[i]**2
             self.vz.append(0.)
             if (radius2 <= self.rdisk**2):
                self.vx.append(-vzero*self.y[i]/self.rdisk)
                self.vy.append(vzero*self.x[i]/self.rdisk)
             else:
                self.vx.append(0.)
                self.vy.append(0.)
             if distri==2:
                 self.u.append(przero/((self.gamma - 1.)*self.rho[i]))
             else:
                 self.u.append(przero/((self.gamma - 1.)*self.getrhoi(i)))
             self.Bx.append(Bzero)
             self.By.append(0.)
             self.Bz.append(0.)
             
    def getrhoi(self,i):
        radius2 = self.x[i]**2 + self.y[i]**2
        if (radius2 <= self.rdisk**2):
            return self.rhodisk
        else:
            return self.rhozero
        
