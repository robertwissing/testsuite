## Setup for ORZAG-TANG
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
import readtipsy as tip
import IC_denstable as denstable
from scipy.integrate import simps, odeint  # Add simps to imports
from scipy.interpolate import interp1d
#from createplanet_2025_b import get_inputs, run_model
from createplanet_2025_b import run_or_load

class setup_planet(object):
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
    deltastep=0.1
    freqout=10
    nsteps=3000
    dmsolunit=1.0
    dkpcunit=1.0
    adi=2
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
    Mat1=[]
    Mat2=[]
    Mat3=[]
    Mat4=[]
    Mat5=[]
    dxbound=[]
    dybound=[]
    dzbound=[]
    totvol=[]
    rcylmin2= 0
    rmin2= 0
    rcylmax2= 10**22
    rmax2= 10**22
    M=1.0
    R=1.0
    Rtab = 1.0
    eps = 0.01
    materialfile = ""
   
    def __init__(self):
        pass
    def create(self,nx,distri=0,vm=0,entry='glassreadyfile',planet='Earth_PES095_2layer.dat'):
         v0 = 0.0
         self.materialfile = planet+".dat_material";
         #inputs = get_inputs( planet )
         #r_values,rho_values,T_values,M,R = run_model(inputs)
         r_values,rho_values,T_values,Mat_values,M,R = run_or_load(planet, force=False)
         print(r_values,rho_values,T_values,M,R)
         ## We want stuff in units of Earth mass and Earth Radius
         M_earth = 5.972*10**24 #kg
         R_earth = 6378*1000 # m
         M_sol = 1.989*10**30 #kg
         R_kpc = 3.086*10**19 #m
         #we set M_earth and R_earth as code units.
         #1 Code unit =  1 Msol unit * M_earth/M_sol
         self.dmsolunit = M_earth/M_sol
         self.dkpcunit = R_earth/R_kpc
         Rho_earth = M_earth/R_earth**3 
         M = M / M_earth
         R = R / R_earth
         r_values = r_values / R_earth
         rho_values = rho_values / Rho_earth
         ## T_values is in kelvin
         unique, counts = np.unique(r_values, return_counts=True)
         duplicates = unique[counts > 1]
         for dup in duplicates:
            idx = np.where(r_values == dup)[0]
            idx = idx[1]
            r_values = np.delete(r_values, idx)
            rho_values = np.delete(rho_values, idx)
            T_values = np.delete(T_values, idx)
            Mat_values = np.delete(Mat_values, idx)

         print(r_values[-1],r_values[-2])
         print(duplicates)
         self.R = R
         self.M = M
         #Get r_values and densities and temperature from planet file.
         # Create an interpolator for theta(r)
         mindens=rho_values[-1]*0.5
         print(r_values)
         print(rho_values)
         self.density_interp = interp = interp1d(
              r_values, 
              rho_values,
              bounds_error=False,
              fill_value=(mindens),
              kind='linear'
          )
         self.temperature_interp = interp = interp1d(
              r_values,
              T_values,
              bounds_error=False,
              fill_value=(T_values[-1]),
              kind='linear'
          )
         self.material_interp = interp = interp1d(
              r_values,
              Mat_values,
              bounds_error=False,
              fill_value=(Mat_values[-1]),
              kind='linear'
          )

         mass_integrand = 4 * np.pi * r_values**2 * rho_values
         M3 = simps(mass_integrand, r_values)  # Numerical integration
         print("my mass",self.M,M3,"rmax",r_values[-1])
         Bzero = 0.0;
         #boundaries    
         dz=self.R*1.25
         dy=self.R*1.25
         dx=self.R*1.25
         distribute.setbound(self,-dx,dx,-dy,dy,-dz,dz)
         deltax = 1.0/nx ##deltax based on earth radius (1.0)
         print('DISTRI: ',distri)
         if distri==0:
             print("DOES NOT WORK WITH LATTICE SETUP AT THE MOMENT (add stretch mapping)")
             exit
             distribute.setcloseddist(self,-dx,dx,-dy,dy,-dz,dz,deltax)
         elif distri==1:
            distribute.setrandomdist(self,-dx,dx,-dy,dy,-dz,dz,deltax)            
         else:
            tgdatain,tddata,tsdata,data_header,time=tip.readtipsy(entry);
            Rpart = np.sqrt(tgdatain[:,1]*tgdatain[:,1]+tgdatain[:,2]*tgdatain[:,2]+tgdatain[:,3]*tgdatain[:,3])
            mapp = (Rpart < self.R)
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
         
         self.Mat1 = np.zeros(self.npart)
         self.Mat2 = np.zeros(self.npart)
         self.Mat3 = np.zeros(self.npart)
         self.Mat4 = np.zeros(self.npart)
         self.Mat5 = np.zeros(self.npart)
         
         for i in range(self.npart):
             self.vx.append(v0*self.x[i])
             self.vy.append(v0*self.y[i])
             self.vz.append(v0*self.z[i])
             self.u.append( self.gettemperaturei(i)  )
             mymaterial = self.getmateriali(i)
             if mymaterial < 0.5:
                 self.Mat1[i] = 1.0
             elif mymaterial < 1.5:
                 self.Mat2[i] = 1.0
             elif mymaterial < 2.5:
                 self.Mat3[i] = 1.0
             elif mymaterial < 3.5:
                 self.Mat4[i] = 1.0
             elif mymaterial < 4.5:
                 self.Mat5[i] = 1.0
             self.Bx.append(0.)
             self.By.append(0.)
             self.Bz.append(Bzero)

    def gettemperaturei(self,i):
        r=np.sqrt(self.x[i]**2+self.y[i]**2+self.z[i]**2)
        return self.temperature_interp(r)

    def getmateriali(self,i):
        r=np.sqrt(self.x[i]**2+self.y[i]**2+self.z[i]**2)
        return self.material_interp(r)

    def getrhoi(self,i):
        r=np.sqrt(self.x[i]**2+self.y[i]**2+self.z[i]**2)
        return self.density_interp(r)

    def getrho(self,r):
        return self.density_interp(r)
    
    def getrho_vec(self,r_vals):
        return self.density_interp(r_vals)
    
