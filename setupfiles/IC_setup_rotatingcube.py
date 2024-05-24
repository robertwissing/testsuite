## Setup for ORZAG-TANG
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
class setup_rotatingcube(object):
    npart=0
    ngas=0
    ndark=0
    nstar=0
    ndim=3
    time=0
    gamma=5./3.
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
    rhocube=0.
    rcube=0.
    
    def __init__(self):
        pass
    def create(self,nx,vin):
        
        
         #--setup parameters
         const = np.sqrt(4.*np.pi)
         self.rcube = 0.25           # length of innercube
         vzero = vin          # rotation speed of initial disk
         Bzero = 5.0/const      # uniform field in bx direction
         przero = 3.75             # initial pressure
         self.rhozero = 7./4.            # ambient density
         self.rhocube = 7.0           # density of rotating disk
     
         print('Setup for MHD rotor problem')
         print(' radius of disk        = ',self.rdisk,', rotation speed = ',vzero,
                    ' initial B   = ',Bzero,', density = ',self.rhozero,', pressure = ',przero)
    
    
    
         #boundaries    
         #dz = 6*deltax
 #        deltax = 1.0/nx
 #        deltadisk = deltax*(rhozero/rhodisk)**(1./3.)
 #        nxdisk=(2*rdisk)/deltadisk+(1.0-2*rdisk)
         dz=0.5
         ymin=-0.5
         xmin=-0.5
         distribute.setbound(self,xmin,-xmin,ymin,-ymin,-dz,dz)
         deltax = self.dxbound/nx
         deltacube = deltax*(self.rhozero/self.rhocube)**(1./3.)
         self.rmin2=self.rcube
         distribute.setcloseddist(self,xmin,0.5,ymin,0.5,-dz,dz,deltax)
         self.rmin2=0
         self.rmax2=self.rcube
         distribute.setcloseddist(self,xmin,0.5,ymin,0.5,-dz,dz,deltacube)
         self.npart=len(self.x)
         self.ngas=self.npart
         totmass = (self.rhozero*self.dxbound*self.dybound*self.dzbound+(self.rhocube-self.rhozero)*self.rcube**3)
         self.mass = [totmass/self.npart]*self.npart
         print('npart = ',self.npart)
         self.h=smth.getsmooth(self,200)
         print(self.h/deltax)
    
         for i in range(self.npart):
             self.vz.append(0.)
             self.vx.append(-vzero*self.y[i]/self.dxbound)
             self.vy.append(vzero*self.x[i]/self.dybound)
             self.u.append(przero/((self.gamma - 1.)*self.getrhoi(i)))
             self.Bx.append(0.)
             self.By.append(0.)
             self.Bz.append(0.)
             
    def getrhoi(self,i):
        if (D.x[i] <= D.rcube and D.y[i] <= D.rcube and D.z[i] <= D.rcube):
            return self.rhocube
        else:
            return self.rhozero