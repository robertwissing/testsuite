## Setup for ORZAG-TANG
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
import readtipsy as tip
import IC_denstable as denstable
from scipy.integrate import simps, odeint  # Add simps to imports
from scipy.interpolate import interp1d

class setup_polytrope(object):
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
    eps = 0.01
    Rtab = 1.0
    n = 1/(gamma-1.0)
    rho_c = 1.0
    K = 1.0
    alpha = np.sqrt((n + 1) * K * rho_c**(1.0/n - 1.0) / (4 * np.pi))
    theta_crit = 0.005
    
    def __init__(self):
        pass
    def create(self,nx,distri=0,vm=0,entry='glassreadyfile',gamma=5./3.):
         v0=0.0
         self.gamma = gamma
         self.n = n = 1.0 / (gamma - 1.0)
         R=self.R
         M=self.M
         # Step 1: Solve dimensionless Lane-Emden equation
         xi_values, theta_values, dtheta_values = self.solve_lane_emden()
         xi_1 = xi_values[-1]
         dtheta_xi1 = dtheta_values[-1]  # Surface derivative (negative)
        
         # Physical constants
         G = 1.0  # Gravitational constant
        
         # Step 2: Compute alpha
         self.alpha = R / xi_1
        
         # Step 3: Compute rho_c from mass condition
         self.rho_c = M / (-4 * np.pi * self.alpha**3 * xi_1**2 * dtheta_xi1)
        
         # Step 4: Compute K
         self.K = (self.alpha**2 * 4 * np.pi * G) / ((n + 1) * self.rho_c**((1 - n)/n))
         print("MY K and n     self.n*self.K*(self.rho[i]**(1/self.n)   ", self.K, self.n);      
         # Step 5: Verify and compute profile
         r_values = self.alpha * xi_values
         densities = self.rho_c * theta_values**n
        
         # Create an interpolator for theta(r)
         mindens=densities[-1]*0.5
         self.density_interp = interp = interp1d(
              r_values, 
              densities,
              bounds_error=False,
              fill_value=(mindens),
              kind='cubic'
          )
         
         mass_integrand = 4 * np.pi * r_values**2 * densities
         M3 = simps(mass_integrand, r_values)  # Numerical integration
         surface_xi = xi_values[-1]
         surface_dtheta = dtheta_values[-1]
         M2 = -4 * np.pi * self.alpha**3 * self.rho_c * surface_xi**2 * surface_dtheta
         print("thet",theta_values[-1]);
         print("my mass",self.M,M2,M3,"rmax",r_values[-1])
         
         Bzero = 0.0;
         #boundaries    
         dz=self.R+self.R*0.25
         dy=self.R+self.R*0.25
         dx=self.R+self.R*0.25
         distribute.setbound(self,-dx,dx,-dy,dy,-dz,dz)
         deltax = self.dxbound/nx
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
             self.vx.append(v0*self.x[i])
             self.vy.append(v0*self.y[i])
             self.vz.append(v0*self.z[i])
             self.u.append(self.n*self.K*(self.rho[i]**(1/self.n)))
             self.Bx.append(0.)
             self.By.append(0.)
             self.Bz.append(Bzero)

    def getrhoi(self,i):
        r=np.sqrt(self.x[i]**2+self.y[i]**2+self.z[i]**2)
        return self.density_interp(r)

    def getrho(self,r):
        return self.density_interp(r)
    
    def getrho_vec(self,r_vals):
        return self.density_interp(r_vals)
    
    def solve_lane_emden(self, xi_max=10, num_points=10000):
        """Solves Lane-Emden equation and returns xi, theta, and dtheta"""
        xi_start = 0.0
        xi = np.linspace(xi_start, xi_max, num_points)
        theta0 = 1.0 - (self.n * xi_start**2) / 6.0
        dtheta0 = - (self.n * xi_start) / 3.0
        y0 = [theta0, dtheta0]
        
        sol = odeint(self._lane_emden_ode, y0, xi, args=(self.n,))
        theta = sol[:, 0]
        dtheta = sol[:, 1]  # dtheta/dxi
       
        surface_mask = (theta**self.n <= self.theta_crit) | np.isnan(theta)
        if np.any(surface_mask):
            surface_idx = np.argmax(surface_mask)
            xi = xi[:surface_idx]
            theta = theta[:surface_idx]
            dtheta = dtheta[:surface_idx]
        return xi, theta, dtheta  # Return dtheta for analytical mass

    def _lane_emden_ode(self, y, xi, n):
        theta, dtheta = y
        dydxi = [dtheta, -2*dtheta/xi - theta**n] if xi > 1e-6 else [dtheta, -theta**n]
        return dydxi
