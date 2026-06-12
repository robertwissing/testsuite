## Setup for HYDROSTATIC SQUARE test
## Saitoh & Makino 2013; Hopkins 2013, 2015; Hu et al. 2014.
## High-density square in pressure equilibrium with a low-density medium.
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
import readtipsy as tip

class setup_square(object):
    dICdensRsmooth = 0.01
    dICdensprofile = 2
    dICdensdir = 0
    dICdensR = 0.25
    dICdensinner = 4.0
    dICdensouter = 1.0
    rhoit=0
    inflow=0
    ns=64
    npart=0
    ngas=0
    ndark=0
    nstar=0
    ndim=3
    time=0
    gamma=5./3.
    periodic=1
    deltastep=0.02
    nsteps=200
    freqout=20
    dmsolunit=1.0
    dkpcunit=1.0
    adi=1
    grav=0
    cosmo=0
    molweight=1.0
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
    dxbound=[]
    dybound=[]
    dzbound=[]
    totvol=[]
    rcylmin2= 0
    rmin2= 0
    rcylmax2= 10**22
    rmax2= 10**22
    rhozero=1.
    rhodens=4.
    przero=2.5
    shape=0

    def __init__(self):
        pass

    def create(self,nx,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000',rhodiff=4.0,przero=2.5,gamma=5./3.,is3d=0):
        is3d = int(is3d)
        self.gamma   = gamma
        self.przero  = przero
        self.rhodens = self.rhozero*rhodiff
        self.dICdensinner = self.rhodens
        self.dICdensouter = self.rhozero
        gam1 = self.gamma - 1.

        mode = '3D cube (Sandnes 2025)' if is3d else '2D square thin-slab (Saitoh-Makino 2013)'
        print('Setup for hydrostatic-square test [',mode,']...')
        print(' rho_in = ',self.rhodens,', rho_out = ',self.rhozero,
              ', chi = ',rhodiff,', P = ',self.przero,', gamma = ',self.gamma)

        self.dxbound = 1.
        deltax = self.dxbound/nx
        deltadens = deltax*(self.rhozero/self.rhodens)**(1./3.)
        dx = self.dxbound/2.
        dy = dx
        if is3d:
            dz = dx
        else:
            dz = 2*np.sqrt(6)*deltax
        R = self.dICdensR
        distribute.setbound(self,-dx,dx,-dy,dy,-dz,dz)
        # shape=1 makes checkboundary cube the inner/outer region via rmin2/rmax2:
        # (xi^2<rmax2 AND yi^2<rmax2 AND zi^2<rmax2) AND (rmin2 <= xi^2 OR yi^2 OR zi^2).
        # In 3D mode this is a literal inner-cube carve. In 2D-slab mode the z-test
        # collapses to "always pass" only while dz <= R (else inner particles get clipped in z).
        self.shape = 1
        if (not is3d) and dz > R:
            raise ValueError(
                "IC_setup_square: slab dz={:.4f} exceeds inner half-side R={:.4f}; "
                "the cubic-shell carve will clip inner particles in z. "
                "Use nx >= {:d} or pass is3d=1.".format(dz, R, int(np.ceil(2*np.sqrt(6)/R))))

        if distri==0 or distri==1:
            setdist = distribute.setcloseddist if distri==0 else distribute.setrandomdist
            # Outer (low-density) lattice over the full box, excluding the inner cube
            self.rmin2 = R*R
            self.rmax2 = (dx*dx)*4   # generous: never clips the [-dx,dx] box
            setdist(self,-dx,dx,-dy,dy,-dz,dz,deltax)
            # Inner (high-density) lattice confined to the inner cube/square
            self.rmin2 = 0.0
            self.rmax2 = R*R
            setdist(self,-R,R,-R,R,-dz,dz,deltadens)
            # Restore non-carving defaults for downstream code
            self.rmin2 = 0.0
            self.rmax2 = 10**22
        else:
            tgdata,tddata,tsdata,data_header,time=tip.readtipsy(entry);
            self.x=tgdata[:,1]
            self.y=tgdata[:,2]
            self.z=tgdata[:,3]
            self.rho=tgdata[:,7]

        self.npart = len(self.x)
        self.ngas  = self.npart
        if is3d:
            V_inner = (2*R)**3
            V_outer = self.dxbound*self.dybound*self.dzbound - V_inner
            totmass = V_inner*self.rhodens + V_outer*self.rhozero
        else:
            A_inner = (2*R)**2
            A_outer = self.dxbound*self.dybound - A_inner
            totmass = (A_inner*self.rhodens + A_outer*self.rhozero)*self.dzbound
        self.mass = [totmass/self.npart]*self.npart
        if vm==1:
            print("Varying masses -> ONLY FOR RAND + RELAXING TO GLASS");
            for i in range(0,self.npart,2):
                dmrat=0.25*self.mass[i]*np.random.rand()
                self.mass[i]=self.mass[i]+dmrat
                self.mass[i+1]=self.mass[i+1]-dmrat
        self.h = smth.getsmooth(self,200)
        print('npart = ',self.npart,' particle mass = ',self.mass[0],' total mass = ',totmass)

        for i in range(self.npart):
            self.vx.append(0.)
            self.vy.append(0.)
            self.vz.append(0.)
            self.Bx.append(0.)
            self.By.append(0.)
            self.Bz.append(0.)
            if distri==2:
                self.u.append(self.przero/(gam1*self.rho[i]))
            else:
                self.u.append(self.przero/(gam1*self.getrhoi(i)))

    def getrhoi(self,i):
        # Inner region is a cube in all 3 axes; in 2D-slab mode |z|<dz<=R holds
        # for every slab particle, so the z-test is a no-op there.
        if (abs(self.x[i])<=self.dICdensR
            and abs(self.y[i])<=self.dICdensR
            and abs(self.z[i])<=self.dICdensR):
            return self.rhodens
        else:
            return self.rhozero
