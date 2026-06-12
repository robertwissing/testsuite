## Setup for the ISENTROPIC (YEE) VORTEX
## Yee, Vinokur & Djomehri (2000) JCP 162, 33; widely used (e.g. Springel 2010,
## Hopkins 2015) as a SMOOTH, steady, exact solution of the Euler equations to
## measure numerical dissipation and the formal order of convergence.
##
## A vortex perturbation is superimposed on a uniform background (rho=P=T=1):
##   du = -(beta_v/2pi) * y * exp((1-r^2)/2)
##   dv =  (beta_v/2pi) * x * exp((1-r^2)/2)
##   dT = -(gamma-1) beta_v^2/(8 gamma pi^2) * exp(1-r^2)
##   T  = 1 + dT ,   rho = T^(1/(gamma-1)) ,   P = rho^gamma  (entropy = const)
## The centrifugal force balances the pressure gradient -> steady state; with a
## uniform mean flow (vadv,vadv) the vortex simply advects unchanged. After one
## box-crossing it should return to the initial state (Galilean-invariance /
## low-dissipation test).
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
import readtipsy as tip
class setup_isentropicvortex(object):
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
    gamma=1.4
    periodic=1
    grav=0
    deltastep=0.01
    nsteps=1000   # end time = nsteps*deltastep = 10.0 (one box crossing for vadv=1)
    dmsolunit=1.0
    dkpcunit=1.0
    molweight=1.0
    adi=1
    beta_v=5.0
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
    def create(self,nx,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000',beta_v=5.0,vadv=1.0):
    #--setup parameters
    # beta_v : vortex strength (Yee standard = 5)
    # vadv   : uniform mean-flow velocity applied to both x and y (advection; 0 = stationary)
     self.beta_v = beta_v
     self.rhozero = 1.0
     self.dICdensinner=self.rhozero
     self.dICdensouter=self.rhozero
     gam1 = self.gamma - 1.

     print('Setup for isentropic (Yee) vortex...')
     print(' beta_v = ',beta_v,', vadv = ',vadv,', gamma = ',self.gamma)

     # Domain [-5,5]^2 so the vortex perturbation has decayed to ~0 at the
     # periodic boundary (exp((1-25)/2) ~ 6e-6).
     self.dxbound=10.

     dx=5.0
     dy=5.0
     deltax = self.dxbound/nx
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
     # Equal-mass particles sample the (non-uniform) density via the glass relaxation;
     # mean density is rhozero so total mass = rhozero * box volume.
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
     dTfac = gam1*beta_v**2/(8.*self.gamma*np.pi**2)
     velfac = beta_v/(2.*np.pi)
     for i in range(self.npart):
         r2 = self.x[i]**2 + self.y[i]**2
         expf = np.exp(0.5*(1.-r2))
         T = 1. - dTfac*np.exp(1.-r2)
         # u = P/((gamma-1) rho) = T/(gamma-1) for this isentropic background
         self.u.append(T/gam1)
         self.vx.append(vadv - velfac*self.y[i]*expf)
         self.vy.append(vadv + velfac*self.x[i]*expf)
         self.vz.append(0.)
         self.Bx.append(0.)
         self.By.append(0.)
         self.Bz.append(0.)


    def getrhoi(self,i):
        # rho(r) = T(r)^(1/(gamma-1)) with T the isentropic-vortex temperature profile
        gam1 = self.gamma - 1.
        r2 = self.x[i]**2 + self.y[i]**2
        T = 1. - gam1*self.beta_v**2/(8.*self.gamma*np.pi**2)*np.exp(1.-r2)
        return T**(1./gam1)
