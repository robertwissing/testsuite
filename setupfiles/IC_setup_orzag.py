## Setup for ORZAG-TANG
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
import readtipsy as tip
class setup_orzag(object):
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
    shape=0
    ndim=3
    time=0
    cosmo=0
    freqout=25
    ns=64
    gamma=5./3.
    periodic=1
    grav=0
    deltastep=0.01
    nsteps=250
    dmsolunit=1.0
    dkpcunit=1.0
    molweight=1.0
    adi=1
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
    rhozero= 1.0
    
    def __init__(self):
        pass
    def create(self,nx,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000'):
    #--setup parameters
     const = 4.*np.pi
     betazero = 10./3.
     machzero = 1.0
     vzero = 1.0
     Bzero = 1.0/np.sqrt(const)
     przero = 0.5*Bzero**2*betazero
     rhozero = self.gamma*przero*machzero
     self.rhozero = rhozero
     self.dICdensinner=self.rhozero
     self.dICdensouter=self.rhozero
     gam1 = self.gamma - 1.
     uuzero = przero/(gam1*rhozero)

     print('Setup for 3D Orszag-Tang vortex problem...')
     print(' beta        = ',betazero,', mach number = ',machzero,
                ' initial B   = ',Bzero,', density = ',rhozero,', pressure = ',przero)

     self.dxbound=1.

     deltax = self.dxbound/nx
     dz=2*np.sqrt(6)*deltax
     hfact=5
     dy=0.5
     dx=0.5
     distribute.setbound(self,-dx,dx,-dy,dy,-dz,dz)

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
     totmass = self.dxbound*self.dybound*self.dzbound*self.rhozero
     self.mass = [totmass/self.npart]*self.npart
     if vm==1:
        print("Varying masses -> ONLY FOR RAND + RELAXING TO GLASS");
        for i in range(0,self.npart,2):
            dmrat=0.25*self.mass[i]*np.random.rand()
            self.mass[i]=self.mass[i]+dmrat
            self.mass[i+1]=self.mass[i+1]-dmrat
            
     print('npart = ',self.npart,' particle mass = ',self.mass[0])
     self.h=smth.getsmooth(self,200)
     for i in range(self.npart):
         self.vx.append(-vzero*np.sin(2.*np.pi*(self.y[i]+dy)))
         self.vy.append(vzero*np.sin(2.*np.pi*(self.x[i]+dx)))
         self.vz.append(0.01*vzero)
         if(self.x[i]<-0.498 and i==0):
             print("i",i,"x:",self.x[i],"vy:",self.vy[i])
         self.u.append(uuzero)
         self.Bz.append(0.)
         self.Bx.append(-Bzero*np.sin(2.*np.pi*(self.y[i]+dy))) 
         self.By.append(Bzero*np.sin(4.*np.pi*(self.x[i]+dx)))

             
    def getrhoi(self,i):
            return self.rhozero
