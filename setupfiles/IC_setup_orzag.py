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
    def create(self,nx,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000',dim3=0):
     # dim3 != 0 selects the true 3D Orszag-Tang vortex (Helzel et al. 2011);
     # the default (dim3 == 0) keeps the 2D vortex on a thin z-slab.
     if int(dim3) != 0:
         self.create3d(nx,distri=distri,vm=vm,entry=entry)
         return
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

    def create3d(self,nx,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000'):
     #--True 3D Orszag-Tang vortex.
     #  References:
     #    C. Helzel, J. A. Rossmanith & B. Taetz, J. Comput. Phys. 230 (2011) 3803
     #      (arXiv:1007.2606) -- original 3D extension of the OT vortex.
     #    X. Tu et al., arXiv:2202.03761, "Meshless methods for MHD with vector
     #      potential" (the SPH-relevant reference).
     #    M. Rossazza et al. (2025).
     #  Normalisation: we follow the Tu/Rossazza UNIT-BOX [0,1]^3 form, which is
     #  the same B0=1/sqrt(4pi), beta=10/3 scaling as this code's 2D OT vortex
     #  (k=2*pi). Helzel's equivalent [0,2*pi]^3, rho=gamma^2 scaling is the same
     #  problem stretched ~2*pi in length, so its eddy-turnover time is ~2*pi
     #  longer -- using the unit box keeps the evolution time consistent with the
     #  2D test and with the t=1 figures in Tu (2022)/Rossazza (2025).
     #  Initial conditions on the periodic cube [0,1]^3 (eps = 0.2):
     #    rho = 25/(36*pi),  p = 5/(12*pi),  gamma = 5/3
     #    v = (-(1+eps sin 2pi z) sin 2pi y, (1+eps sin 2pi z) sin 2pi x,
     #          eps sin 2pi z)
     #    B = (-B0 sin 2pi y, B0 sin 4pi x, 0),  B0 = 1/sqrt(4*pi)
     const = 4.*np.pi
     betazero = 10./3.
     machzero = 1.0
     vzero = 1.0
     eps = 0.2
     Bzero = 1.0/np.sqrt(const)
     przero = 0.5*Bzero**2*betazero
     rhozero = self.gamma*przero*machzero
     self.rhozero = rhozero
     self.dICdensinner = rhozero
     self.dICdensouter = rhozero
     gam1 = self.gamma - 1.
     uuzero = przero/(gam1*rhozero)

     print('Setup for the true 3D Orszag-Tang vortex problem (Tu 2022 / Helzel et al. 2011)...')
     print(' beta        = ',betazero,', mach number = ',machzero,
                ', initial B   = ',Bzero,', density = ',rhozero,
                ', pressure = ',przero,', eps = ',eps)

     half = 0.5
     deltax = (2.*half)/nx
     distribute.setbound(self,-half,half,-half,half,-half,half)

     if distri==0:
         distribute.setcloseddist(self,-half,half,-half,half,-half,half,deltax)
     elif distri==1:
         distribute.setrandomdist(self,-half,half,-half,half,-half,half,deltax)
     else:
         tgdata,tddata,tsdata,data_header,time=tip.readtipsy(entry);
         self.x=tgdata[:,1]
         self.y=tgdata[:,2]
         self.z=tgdata[:,3]
         self.rho=tgdata[:,7]
     self.npart=len(self.x)
     self.ngas=self.npart
     totmass = self.dxbound*self.dybound*self.dzbound*rhozero
     self.mass = [totmass/self.npart]*self.npart
     if vm==1:
        print("Varying masses -> ONLY FOR RAND + RELAXING TO GLASS");
        for i in range(0,self.npart,2):
            dmrat=0.25*self.mass[i]*np.random.rand()
            self.mass[i]=self.mass[i]+dmrat
            self.mass[i+1]=self.mass[i+1]-dmrat

     print('npart = ',self.npart,' particle mass = ',self.mass[0])
     self.h=smth.getsmooth(self,200)
     twopi = 2.*np.pi
     for i in range(self.npart):
         # shift box-centred coordinates into [0,1] for the field pattern
         xx = self.x[i] + half
         yy = self.y[i] + half
         zz = self.z[i] + half
         envz = 1. + eps*np.sin(twopi*zz)
         self.vx.append(-vzero*envz*np.sin(twopi*yy))
         self.vy.append( vzero*envz*np.sin(twopi*xx))
         self.vz.append( vzero*eps*np.sin(twopi*zz))
         self.u.append(uuzero)
         self.Bx.append(-Bzero*np.sin(twopi*yy))
         self.By.append( Bzero*np.sin(2.*twopi*xx))
         self.Bz.append(0.)

    def getrhoi(self,i):
            return self.rhozero
