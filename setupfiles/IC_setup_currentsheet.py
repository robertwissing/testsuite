## Setup for the MHD CURRENT SHEET test
## Gardiner & Stone (2005) JCP 205, 509; Hopkins & Raives (2016) MNRAS 455, 51.
## Periodic 2D box with two antiparallel current sheets. A small velocity
## perturbation drives reconnection at the field reversals -> a robustness /
## divergence-handling test (complements orzag / mhdloop / mhddivtest).
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
import readtipsy as tip
class setup_currentsheet(object):
    dICdensRsmooth = 0.1
    dICdensprofile = 1
    dICdensdir = 1
    dICdensR = 1.0
    dICdensinner = 1.0 # calculated depending on mass
    dICdensouter = 1.0
    inflow=0
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
    nsteps=500   # end time = nsteps*deltastep = 5.0 (Hopkins & Raives 2016)
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
    def create(self,nx,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000',beta=0.1,v0=0.1):
    #--setup parameters
    # beta : plasma beta -> sets the (uniform) gas pressure przero = beta*B0^2/2
    # v0   : amplitude of the sinusoidal vx perturbation that drives reconnection
     Bzero = 1.0
     przero = beta*0.5*Bzero**2
     rhozero = 1.0
     self.rhozero = rhozero
     self.dICdensinner=self.rhozero
     self.dICdensouter=self.rhozero
     gam1 = self.gamma - 1.
     uuzero = przero/(gam1*rhozero)

     print('Setup for 2D MHD current sheet problem...')
     print(' beta        = ',beta,', v0 = ',v0,
                ' initial B   = ',Bzero,', density = ',rhozero,', pressure = ',przero)

     self.dxbound=1.

     deltax = self.dxbound/nx
     dz=2*np.sqrt(6)*deltax
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
         # Sinusoidal shear perturbation (one period across the box height)
         self.vx.append(v0*np.sin(2.*np.pi*(self.y[i]+dy)))
         self.vy.append(0.)
         self.vz.append(0.)
         self.u.append(uuzero)
         self.Bx.append(0.)
         self.Bz.append(0.)
         # Two antiparallel current sheets at x = -0.25 and x = +0.25:
         # By = -B0 in the inner half (|x|<0.25), +B0 in the outer regions.
         if abs(self.x[i]) < 0.25:
             self.By.append(-Bzero)
         else:
             self.By.append(Bzero)


    def getrhoi(self,i):
            return self.rhozero
