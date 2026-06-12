## Setup for ORZAG-TANG
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
import readtipsy as tip
import IC_denstable as denstable
class setup_evrard(object):
    dICdensRsmooth = 0.1
    dICdensprofile = 8
    dICdensdir = 4
    dICdensR = 1.0 
    dICdensinner = 1.0 # calculated depending on mass
    dICdensouter = 1.0
    inflow=0
    rhoit=0
    npart=0
    ngas=0
    ndark=0
    nstar=0
    ndim=3
    ns=64
    time=0
    grav=1
    cosmo=0
    molweight=1.0
    gamma=1.6666666666667
    periodic=0
    deltastep=0.01
    freqout=1
    nsteps=300
    dmsolunit=1.0
    dkpcunit=1.0
    adi=1
    shape=2
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
    rhozero= 1.0
    M=1.0
    R=1.0
    eps = 0.02
    Rtab = 1.0
    def __init__(self):
        pass
    def create(self,nx,distri=0,vm=0,entry='glassreadyfile',beta=1):
        
         gam1=self.gamma-1
         self.rhozero=1.0
         term1 = (self.R**2)/2 - self.eps*self.R + self.eps**2 * np.log((self.R + self.eps)/self.eps)
         Mtrue = (2*self.M / self.R**2) * term1
         G = 1.0
         uzero = 0.05 * G * self.M / self.R
         Bzero = 0.0;
         #boundaries    
         dz=self.R+self.R*0.25
         dy=self.R+self.R*0.25
         dx=self.R+self.R*0.25
         distribute.setbound(self,-dx,dx,-dy,dy,-dz,dz)
         deltax = self.dxbound/nx
         print('DISTRI: ',distri)
         if distri==0:
             distribute.setcloseddist(self,-dx,dx,-dy,dy,-dz,dz,deltax)
         elif distri==1:
            distribute.setrandomdist(self,-dx,dx,-dy,dy,-dz,dz,deltax)            
         else:
            tgdatain,tddata,tsdata,data_header,time=tip.readtipsy(entry);
            Rpart = np.sqrt(tgdatain[:,1]*tgdatain[:,1]+tgdatain[:,2]*tgdatain[:,2]+tgdatain[:,3]*tgdatain[:,3])
            mapp = (Rpart < self.R)
            #    (tgdatain[:,1]**2 + tgdatain[:,2]**2 < rdisk**2)
            tgdata=tgdatain[mapp,:]
            self.x=tgdata[:,1]
            self.y=tgdata[:,2]
            self.z=tgdata[:,3]
            self.rho=tgdata[:,7]
            self.mass=tgdata[:,0]

         
            
         self.npart=len(self.x)
         self.ngas=self.npart
         if distri==1:
             totmass = denstable.compute_total_mass(self,512)
             print("Numerically integrated mass(512**3) " , totmass)
             # Generate 1D radial density table
             ntab=10000
             dimtab=1
             self.Rtab = np.sqrt(self.dxbound*self.dxbound+self.dybound*self.dybound+self.dzbound*self.dzbound)
             gridtab, rhotab = denstable.generate_density_table(self,n=ntab,dim=dimtab)
             denstable.write_density_table_xdr(self,"densitytable_xdr",gridtab,rhotab,n=ntab,dim=dimtab)
         else:
             totmass = self.M
             
         self.mass = [totmass/self.npart]*self.npart
         print('npart = ',self.npart,' particle mass = ',self.mass[1])
         #self.h=smth.getsmooth(self,64)
         self.h = [self.eps]*self.npart
         
         for i in range(self.npart):
             self.vx.append(0.)
             self.vy.append(0.)
             self.vz.append(0.)
             self.u.append(uzero) #isothermal initially
             self.Bx.append(0.)
             self.By.append(0.)
             self.Bz.append(Bzero)

    def getrhoi(self,i):
        r=np.sqrt(self.x[i]**2+self.y[i]**2+self.z[i]**2)
        return (self.M/(2*np.pi*self.R**2*(r+self.eps)))

    def getrho(self,r):
        return (self.M/(2*np.pi*self.R**2*(r+self.eps)))
    
    def getrho_vec(self,r_vals):
        return np.where(
            r_vals <= self.R,
            self.M / (2 * np.pi * self.R**2 * (r_vals + self.eps)),
            self.M / (2 * np.pi * self.R**2 * (r_vals + self.eps) * 1.0)
        )
    
