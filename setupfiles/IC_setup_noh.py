## Setup for the NOH problem
## Noh (1987) JCP 72, 78.
## Cold (P->0) uniform gas with a uniform radial inflow velocity |v|=v0 directed
## at the origin. An infinitely strong, self-similar shock forms at the centre and
## propagates outward at speed D = 0.5*(gamma-1)*v0 (= 1/3 for gamma=5/3, v0=1).
## Exact post-shock (stationary) density:
##   1D planar      rho2 = rho0 * (gamma+1)/(gamma-1)        = 4
##   2D cylindrical rho2 = rho0 * ((gamma+1)/(gamma-1))**2   = 16
##   3D spherical   rho2 = rho0 * ((gamma+1)/(gamma-1))**3   = 64
## A stringent test of shock capturing, wall heating and symmetry preservation.
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
import readtipsy as tip
class setup_noh(object):
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
    freqout=50
    ns=64
    gamma=5./3.
    periodic=1
    grav=0
    deltastep=0.001
    nsteps=600   # end time = nsteps*deltastep = 1.0 (shock radius D*t = 1/3)
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
    def create(self,nx,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000',is3d=1,pini=1e-6,v0=1.0):
    #--setup parameters
    # is3d : 0 = 2D cylindrical (thin slab, rho2=16); 1 = 3D spherical (rho2=64)
    # pini : (tiny) initial pressure -> cold gas
    # v0   : magnitude of the uniform radial inflow velocity (directed at origin)
     self.rhozero = 1.0
     self.dICdensinner=self.rhozero
     self.dICdensouter=self.rhozero
     gam1 = self.gamma - 1.
     uuzero = pini/(gam1*self.rhozero)

     print('Setup for Noh problem...')
     print(' is3d = ',is3d,', v0 = ',v0,', density = ',self.rhozero,', pressure = ',pini)

     self.dxbound=1.

     deltax = self.dxbound/nx
     dx=1.0
     dy=1.0
     if is3d==1:
         dz=1.0
     else:
         dz=2*np.sqrt(6)*deltax
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
     self.h=smth.getsmooth(self,self.ns)
     for i in range(self.npart):
         if is3d==1:
             radius = np.sqrt(self.x[i]**2 + self.y[i]**2 + self.z[i]**2)
         else:
             radius = np.sqrt(self.x[i]**2 + self.y[i]**2)
         if radius > 0.0:
             self.vx.append(-v0*self.x[i]/radius)
             self.vy.append(-v0*self.y[i]/radius)
             if is3d==1:
                 self.vz.append(-v0*self.z[i]/radius)
             else:
                 self.vz.append(0.)
         else:
             self.vx.append(0.)
             self.vy.append(0.)
             self.vz.append(0.)
         self.u.append(uuzero)
         self.Bx.append(0.)
         self.By.append(0.)
         self.Bz.append(0.)


    def getrhoi(self,i):
            return self.rhozero
