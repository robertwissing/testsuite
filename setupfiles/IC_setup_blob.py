## Setup for ORZAG-TANG
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
import readtipsy as tip
class setup_blob(object):
    dICdensRsmooth = 0.005
    dICdensprofile = 3
    dICdensdir = 4
    dICdensR = 0.1 
    dICdensinner = 10.0 # calculated depending on mass
    dICdensouter = 1.0
    molweight=0.59259
    cosmo = 0
    periodic = 1
    ns=64
    deltastep=0.004
    adi = 1
    inflow = 1
    cooling = 0
    dxInflow=100.0
    dxOutflow=-100.0
    dyOutflowRight =100.0
    dyOutflowLeft = -100.0
    dMInflow=1.0
    dTInflow=1.0
    dVelInflow=1.0
    dAccInflow=0.0
    dCloudDensity=100.0
    grav = 0
    nsteps=800
    freqout=5
    rhoit=0
    dmsolunit=1.0
    dkpcunit=1.0
    npart=0
    ngas=0
    ndark=0
    nstar=0
    ndim=3
    time=0
    gamma=5./3.
    shape = 0
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
    metals=[]
    dxbound=[]
    dybound=[]
    dzbound=[]
    totvol=[]
    rcylmin2= 0
    rmin2= 0
    rcylmax2= 10**22
    rmax2= 10**22
    rhoin = 10.0
    rhoout = 1.0
    Rin = 0.1
    vzero=1.0
    windref=0
    przero=1
    rhophyswind = 1E-26
    Tcloud = 4000


    def __init__(self):
        pass
    def create(self,nx,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000',rhodiff=10.,beta=-1, mach=2.7, inflow=0, cooling=0, windref=0,rhophyswind=26,Tcloud=4000):
         if inflow == 0:
             self.inflow = 0
         self.windref = windref
         self.Tcloud = Tcloud
         self.rhophyswind = 10**(-rhophyswind)
         self.dCloudDensity = rhodiff
         gam1=self.gamma-1
         self.rhoout  = 1.0
         self.rhoin = rhodiff
         self.dICdensinner = rhodiff
         self.setphysunits()
         przero=self.przero
         Bzero = 0.0
         machno = mach
         csound = np.sqrt(self.gamma*przero/self.rhoout)
         vzero = machno*csound
         self.vzero=machno*csound
         tcrush = self.Rin*np.sqrt(self.rhoin/self.rhoout)/vzero
         taukh = 1.6*tcrush
         self.deltastep=(20*tcrush/self.nsteps)
         #boundaries
         zmin=-5*self.Rin
         if cooling==1:
            zmax=55*self.Rin
         else:
            zmax=15*self.Rin
         if inflow==0:
            zmax=35*self.Rin
         ymin=-5*self.Rin
         ymax=-ymin
         xmin=-5*self.Rin
         xmax=-xmin

         distribute.setbound(self,xmin,xmax,ymin,ymax,zmin,zmax)
         deltax = self.dxbound/nx
         deltasphere = deltax*(self.rhoout/self.rhoin)**(1./3.)
         self.rmin2=self.Rin**2
         if distri==0:
             distribute.setcubicdist(self,xmin,xmax,ymin,ymax,zmin,zmax,deltax)
         elif distri==1:
             print("Doing random dist")
             distribute.setrandomdist(self,xmin,xmax,ymin,ymax,zmin,zmax,deltax)
             print("Finished random dist")
         self.rmin2=0
         self.rmax2=self.Rin**2
         if distri==0:
            distribute.setcubicdist(self,xmin,xmax,ymin,ymax,zmin,zmax,deltasphere)
         elif distri==1:
            distribute.setrandomdist(self,xmin,xmax,ymin,ymax,zmin,zmax,deltasphere)
         else:
            tgdata,tddata,tsdata,data_header,time=tip.readtipsy(entry);
            self.x=tgdata[:,1]
            self.y=tgdata[:,2]
            self.z=tgdata[:,3]
            self.rho=tgdata[:,7]
            m_tot = np.sum(tgdata[:,0])
            x_com = np.sum(tgdata[:,0] * self.x) / m_tot
            y_com = np.sum(tgdata[:,0] * self.y) / m_tot
            z_com = np.sum(tgdata[:,0] * self.z) / m_tot
            self.x -= x_com
            self.y -= y_com
            self.z -= z_com
            self.wrapbound()

         self.npart=len(self.x)
         self.ngas=self.npart
         masscloud = self.rhoin*4/3*np.pi*self.Rin**3
         totmass = (self.rhoout*self.dxbound*self.dybound*self.dzbound + masscloud - self.rhoout*4/3*np.pi*self.Rin**3)
         self.mass = [totmass/self.npart]*self.npart
         hfact3 = 3.0*64.0/(np.pi*2.0*2.0*2.0*4.0);
         hestimate = (hfact3*self.mass[0])**(1./3.)
         if(self.inflow == 1):
             self.dxInflow = -self.dzbound/2 + hestimate*2.0
             self.dxOutflow = self.dzbound/2 - hestimate*2.0
             if (hestimate*2.0 > 0.25):
                 print("H is very large compared to cloud position")
                 exit
             self.dyOutflowLeft = -self.dybound/2
             self.dyOutflowRight = self.dybound/2
             self.dMInflow = np.mean(self.mass)
             self.dTInflow = self.Tcloud*self.dCloudDensity
             self.dVelInflow = vzero
             self.dAccInflow = 0.0
             self.dCloudPosition = -self.dzbound/2 + 0.5 
         tcross = self.dzbound/vzero
         print("csound: ", csound)
         print("crossing time: ",tcross)
         print("crushing time: ", tcrush)
         print("tcross/tcrush: ", tcross/tcrush)
         print("masscloud: ", masscloud)

         print("taukh: ", taukh)

         print('npart = ',self.npart,' particle mass = ',self.mass[1], ' totmass = ', totmass, ' sphere mass in solar mass = ', self.rhoin*4/3*np.pi*self.Rin**3)
         
         for i in range(self.npart):
             radius2 = self.x[i]**2 + self.y[i]**2+self.z[i]**2
             #velvec in lengthunit/s
             if (self.rho[i] > self.rhoout * 1.5):
                self.vz.append(self.getveli(i))
                self.u.append(przero/(gam1*self.rhoin))
             else:
                self.vz.append(self.getveli(i))
                self.u.append(przero/(gam1*self.rhoout))
             if distri == 2:
                 self.u[i] = (przero/(gam1*self.rho[i]))
             self.vx.append(0.0)
             self.vy.append(0.0)
             self.Bx.append(0.)
             self.By.append(Bzero)
             self.Bz.append(0.)
             self.metals.append(1E-6)
             if distri != 1:
                self.z[i]=self.z[i]-(zmax+zmin)/2
         print("wrapping")
         self.wrapbound()
         self.h=[0.1]*self.npart
         
    def getrhoi(self,i):
        radius2=self.x[i]**2 + self.y[i]**2+self.z[i]**2
        if (radius2 < self.Rin**2):
            return self.rhoin
        else:
            return self.rhoout

    def getveli(self,i):
        radius2=self.x[i]**2 + self.y[i]**2+self.z[i]**2
        r=np.sqrt(radius2)
        vout=self.vzero*self.rhoout
        vinner=0.0
        v = vout+(vinner-vout)*1.0/(1.0+np.exp(-(self.dICdensR-r)/self.dICdensRsmooth));
        v = v/self.rho[i];
        if self.windref == 1:
            v -= self.vzero
            v = -v
        if self.windref == 2:
            masscloud = self.rhoin*4/3*np.pi*self.Rin**3
            totmass = (self.rhoout*self.dxbound*self.dybound*self.dzbound + masscloud - self.rhoout*4/3*np.pi*self.Rin**3)
            masswind = totmass-masscloud
            v -= self.vzero*masswind/totmass
            v = -v
        if self.windref == 3:
            v -= self.vzero*0.5
            v = -v
        return v;

    def wrapbound(self):
        self.x=np.array(self.x)
        self.y=np.array(self.y)
        self.z=np.array(self.z)
        self.x=np.mod(self.x + self.dxbound / 2.0, self.dxbound) - self.dxbound / 2.0
        self.y=np.mod(self.y + self.dybound / 2.0, self.dybound) - self.dybound / 2.0
        self.z=np.mod(self.z + self.dzbound / 2.0, self.dzbound) - self.dzbound / 2.0

    def setphysunits(self):
            MSOLG=1.99e33
            gtoMsol=1/MSOLG
            KBOLTZ=1.38e-16
            MHYDR=1.67e-24
            KPCCM=3.085678e21
            cmtoKpc=1/KPCCM
            GCGS=6.67e-8
            self.dmsolunit=self.rhophyswind*pow(self.dkpcunit*KPCCM,3.0)/MSOLG 
            self.przero = self.dCloudDensity*self.Tcloud*self.dkpcunit*KPCCM*KBOLTZ/(self.molweight*self.dmsolunit*MSOLG*MHYDR*GCGS)
