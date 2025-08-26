
## Setup for ORZAG-TANG
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
from numba import njit
import readtipsy as tip
class setup_areablob(object):
    #dICdensRsmooth = 0.001
    #dICdensprofile = 3
    #dICdensdir = 4
    #dICdensR = 0.015
    #dICdensinner = 360.0
    #dICdensouter = 1.0
    dICdensinner = 1000.0
    dICdensouter = 1.0
    dICdensprofile = 3
    dICdensdir = 4
    dICdensRsmooth = 0.0075
    dICdensR = 0.05
    npart=0
    rhoit=0
    ns=64
    ngas=0
    ndark=0
    nstar=0
    ndim=3
    time=0
    gamma=5./3.
    periodic=1
    deltastep=0.0001
    freqout=10
    nsteps=1100
    dmsolunit=1.0
    dkpcunit=1.0
    adi=0
    grav=1
    shape = 0
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
    pot=[]
    dxbound=[]
    dybound=[]
    dzbound=[]
    totvol=[]
    rcylmin2= 0
    rmin2= 0
    rcylmax2= 10**22
    rmax2= 10**22
    rhoin = 1.0
    rhoout = 1.0
    Rin = 0.015

    def __init__(self):
        pass
    def create(self,nx,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000',rhodiff=360.0,multcloud=1):
         """Create particle distribution for MHDCollapse.

         Parameters:
         - nx: Number of particles aalong x dimension.
         - mu: .
         - distri: Distribution type (0: closed, 1: random, 2: read from file).
         - vm: Varying masses (1 for varying masses).
         - entry: File entry for reading particle data (used when distri=2).
         """ 
         gam1=self.gamma-1
         #--setup parameters
         self.Rin=self.dICdensR
         self.rhoin= rhodiff #msol/pc3
         self.rhoout=self.rhoin/rhodiff

         self.dICdensinner=self.rhoin
         self.dICdensouter=self.rhoout
         self.dICdensRsmooth=0.0
         # boundaries
         fac=2.0

         dz = 0.5*fac
         dy = (0.5/multcloud)
         dx = (0.5/multcloud)
         print(multcloud)
         distribute.setbound(self,-dx,dx,-dy,dy,-dz,dz)
         deltax = (self.dxbound)/nx
         deltasphere = deltax*(self.rhoout/self.rhoin)**(1./3.)
         self.rmin2=self.Rin**2
         
         if distri==0:
             distribute.setcloseddist(self,-dx,dx,-dy,dy,-dz,dz,deltax)
         elif distri==1:
             distribute.setrandomdist(self,-dx,dx,-dy,dy,-dz,dz,deltax)
             
         npartout=len(self.x)
         self.rmin2=0
         self.rmax2=self.Rin**2
         
         if distri==0:
             distribute.setcloseddist(self,-dx,dx,-dy,dy,-dz,dz,deltasphere)
         elif distri==1:
             distribute.setrandomdist(self,-dx,dx,-dy,dy,-dz,dz,deltasphere)
         else:
             tgdata,tddata,tsdata,data_header,time=tip.readtipsy(entry);
             self.x=tgdata[:,1]
             self.y=tgdata[:,2]
             self.z=tgdata[:,3]
             self.rho=tgdata[:,7]
        
         if multcloud > 1 and not distri==1:
            term = (np.arange(multcloud) - (multcloud-1)/2)
            x_offsets = term * self.dxbound
            y_offsets = term * self.dybound
            xx, yy = np.meshgrid(x_offsets, y_offsets)
            offsets = np.column_stack((xx.ravel(), yy.ravel()))
            # Initialize lists to hold all replicated points
            all_x = []
            all_y = []
            all_z = []
            all_rho = []
            # Create the replicas
            for dx, dy in offsets:
                all_x.append(self.x + dx)
                all_y.append(self.y + dy)
                all_z.append(self.z)  # z stays the same
                all_rho.append(self.rho)  # z stays the same
            # Combine all replicas into single arrays
            self.x = np.concatenate(all_x)
            self.y = np.concatenate(all_y)
            self.z = np.concatenate(all_z)
            self.rho = np.concatenate(all_rho)
            self.dxbound *= multcloud;
            self.dybound *= multcloud;
         self.npart=len(self.x)
         self.ngas=self.npart
         totmass = (self.rhoout*self.dxbound*self.dybound*self.dzbound + (multcloud**2)*(self.rhoin-self.rhoout)*4/3*np.pi*self.Rin**3)
         self.mass = [totmass/self.npart]*self.npart
         print(self.dxbound,self.dybound,self.dzbound,self.npart)
         if vm==1:
             print("Varying masses -> ONLY FOR RAND + RELAXING TO GLASS");              
             for i in range(0,self.npart,2):
                 dmrat=0.25*self.mass[i]*np.random.rand()
                 self.mass[i]=self.mass[i]+dmrat
                 self.mass[i+1]=self.mass[i+1]-dmrat

         
         for i in range(self.npart):
             self.vx.append(0.)
             self.vy.append(0.)
             self.vz.append(0.)
             if (np.abs(self.z[i]) > 0.75):
                 self.u.append(1E-10)
             else:
                 self.u.append(1.0/self.rho[i])
             self.Bx.append(0.)
             self.By.append(0.)
             self.Bz.append(0.)
             self.h.append((self.mass[i]/self.rho[i])**0.3333333)
        
    def getrhoi(self,i):
        radius2=self.x[i]**2 + self.y[i]**2+self.z[i]**2
        rpartp = np.sqrt(radius2)
        dens =  self.rhoout + (self.rhoin-self.rhoout)*1.0/(1.0+np.exp(-(self.dICdensR-rpartp)/self.dICdensRsmooth));
        return dens
    
    def getP(self,rho):
        gtosol = 5.0279933*10**(-34)*self.dmsolunit
        kmtokpc= 3.24077929*10**(-17)*self.dkpcunit
        cmtokpc= 3.24077929*10**(-22)*self.dkpcunit
        rhocrit=10**(-14)*gtosol/(cmtokpc**3)
        print(rhocrit,rho)
        return (0.2*kmtokpc)**2*rho*np.sqrt(1+(rho/rhocrit)**(4/3))

 
