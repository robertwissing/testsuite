
## Setup for ORZAG-TANG
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
from numba import njit
import readtipsy as tip
class setup_mhdcollapse(object):
    dICdensRsmooth = 0.001
    dICdensprofile = 3
    dICdensdir = 4
    dICdensR = 0.015
    dICdensinner = 360.0
    dICdensouter = 1.0
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
    dmsolunit=10**3
    dkpcunit=10**(-3)
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
    kpctounit=10**3
    msoltounit=10**(-3)
    
    def __init__(self):
        pass
    def create(self,nx,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000',mu=10,rhodiff=360.0,EkoverEp=0.045):
         """Create particle distribution for MHDCollapse.

         Parameters:
         - nx: Number of particles aalong x dimension.
         - mu: .
         - distri: Distribution type (0: closed, 1: random, 2: read from file).
         - vm: Varying masses (1 for varying masses).
         - entry: File entry for reading particle data (used when distri=2).
         """ 
         gam1=self.gamma-1
         kpctounit=self.kpctounit
         msoltounit=self.msoltounit
         #--setup parameters
         self.Rin=0.015*10**(-3)*kpctounit # pc
         Mtot=1*msoltounit # solar mass
         self.rhoin= Mtot/(4/3*np.pi*self.Rin**3) #msol/pc3
         self.rhoout=self.rhoin/rhodiff

         self.dICdensinner=self.rhoin
         self.dICdensouter=self.rhoout

         P = 4.7*10**5 #period yrs
         yrstosec=31556926
         w=(2*np.pi/(P*yrstosec))*(EkoverEp/0.045)**2 # angular momentum in rad/sec
         Bzero=610/mu*10**(-6) #gauss
         
         # boundaries
         dz=0.075*10**(-3)*kpctounit #pc
         dy=0.075*10**(-3)*kpctounit
         dx=0.075*10**(-3)*kpctounit
         distribute.setbound(self,-dx,dx,-dy,dy,-dz,dz)
         deltax = self.dxbound/nx
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
             
         self.npart=len(self.x)
         self.ngas=self.npart
         totmass = (self.rhoout*self.dxbound*self.dybound*self.dzbound + (self.rhoin-self.rhoout)*4/3*np.pi*self.Rin**3)
         self.mass = [totmass/self.npart]*self.npart
         
         if vm==1:
             print("Varying masses -> ONLY FOR RAND + RELAXING TO GLASS");              
             for i in range(0,self.npart,2):
                 dmrat=0.25*self.mass[i]*np.random.rand()
                 self.mass[i]=self.mass[i]+dmrat
                 self.mass[i+1]=self.mass[i+1]-dmrat
                 
         sphmass = (self.npart-npartout)*self.mass[1]/msoltounit
         print('npart = ',self.npart,'npartcloud = ',self.npart-npartout,' n ' , (self.npart-npartout)**(1./3.),' particle mass = ',self.mass[1], ' totmass = ', totmass, ' sphere mass in solar mass = ',sphmass )
         kpctoau=206264806/kpctounit
         rhosm=(1.0/(5.0/kpctoau)**3)*(self.npart**2/436396**2)
         hsm2=1.0/rhosm**(1./3.)
         rhoc=10**-12
         msoltog=1.989e+33
         hsm=(((np.mean(self.mass)/msoltounit)*msoltog/rhoc)**(1./3.))*3.24078e-22*kpctounit
         print(hsm*kpctoau/kpctounit,hsm2*kpctoau/kpctounit)
         print(hsm,hsm2)
         
         #hsm=(5.0/kpctoau)
         for i in range(self.npart):
             radius2 = self.x[i]**2 + self.y[i]**2+self.z[i]**2
             radvec=(self.x[i],self.y[i],self.z[i])
             wvec=(0,0,w)
             velvec=np.cross(wvec,radvec)
             #velvec in lengthunit/s
             if (radius2 < self.Rin**2):
                #self.u.append(self.getP(self.rhoin)/(gam1*self.rhoin))
                self.vx.append(velvec[0])
                self.vy.append(velvec[1])
                self.vz.append(velvec[2])
             else:
                #self.u.append(self.getP(self.rhoout)/(gam1*self.rhoout))
                self.vx.append(0.)
                self.vy.append(0.)
                self.vz.append(0.)
             self.u.append(1.0)
             self.Bx.append(0.)
             self.By.append(0.)
             self.Bz.append(Bzero)
             self.h.append(hsm)

             
         gtosol = 5.0279933*10**(-34)*self.msoltounit
         kmtokpc= 3.24077929*10**(-17)*self.kpctounit
         cmtokpc=3.24077929*10**(-22)*self.kpctounit
         print(10**(-14)*gtosol/cmtokpc**3)
         print((self.rhoin)/(10**(-14)*gtosol/cmtokpc**3))
         print((self.rhoout)/(10**(-14)*gtosol/cmtokpc**3))
         print(self.getP(self.rhoout)/(gam1*self.rhoout))
        
         
         print(velvec[0])
         print(0.2*kmtokpc*14877722399232.906)
         print("CLOUD DENSITY g/cm**3 ", self.rhoin*cmtokpc**3/gtosol)
         
    def getrhoi(self,i):
        radius2=self.x[i]**2 + self.y[i]**2+self.z[i]**2
        if (radius2 < self.Rin**2):
            return self.rhoin
        else:
            return self.rhoout
    
    def getP(self,rho):
        gtosol = 5.0279933*10**(-34)*self.msoltounit
        kmtokpc= 3.24077929*10**(-17)*self.kpctounit
        cmtokpc= 3.24077929*10**(-22)*self.kpctounit
        rhocrit=10**(-14)*gtosol/(cmtokpc**3)
        print(rhocrit,rho)
        return (0.2*kmtokpc)**2*rho*np.sqrt(1+(rho/rhocrit)**(4/3))

 
