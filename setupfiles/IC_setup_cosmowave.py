#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 21:30:46 2018

@author: robertwi
"""

## Setup for ORZAG-TANG
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
class setup_cosmowave(object):
    dICdensRsmooth = 0.1
    dICdensprofile = 1
    dICdensdir = 1
    dICdensR = 1.0 
    dICdensinner = 1.0 # calculated depending on mass
    dICdensouter = 1.0
    rhoit=0
    npart=0
    ngas=0
    ndark=0
    nstar=0
    ndim=3
    time=0
    grav=0
    H0=2.894405
    cosmo=1
    molweight=1.0
    gamma=5./3.
    periodic=1
    deltastep=0.0005
    grav=0
    nsteps=100
    dmsolunit=1.0
    dkpcunit=1000.0
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
    _rho_z_grid=None
    wk=0.
    x0=0.
    xdot0=0.
    rhozero=1.
    drho=0.
    Bzero=1.
    dB=0.
    a=0.
    kappa=1.
    def __init__(self):
        pass
    def create(self,nx,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000',case=3,factor_A=1,factor_G=0):
        # case selects one of the unique test setups from Berlok 2022 (MNRAS 515 3492):
        #   0 = standing Alfven wave           (paper Sec. 5.1.1, eqs. 75, 76)
        #   1 = traveling Alfven wave          (paper Sec. 5.1.2, eqs. 79, 80; n=2)
        #   2 = standing compressible wave     (paper Sec. 5.2.1/2/3, eqs. 81, 82)
        #   3 = traveling compressible gamma=4/3 (paper Sec. 5.2.1, eqs. 83, 84; n=3)
        #   4 = divergence-cleaning test in EdS (this work):
        #       Linearized Tricco et al. 2016 eqs. 23-24 with v=0, transformed to
        #       comoving (B = B_c/a^2, gradient picks up 1/a). The induction eq.
        #       acquires a +2(adot/a) dB_c Hubble drag; substituting dB_c = a^2 g(a)
        #       gives a g'' + (1/2) g' + Omega_h^2 g = 0, which under u = 2*Omega_h*sqrt(a)
        #       reduces to g'' + g = 0. Hyperbolic-only (tau -> inf) closed form:
        #         dB_c (x,a) = (a/a_i)^2 * dB_c_i * cos[2*Omega_h*(sqrt(a)-sqrt(a_i))] * sin(k(x-x0))
        #         dB_phys(x,a) =          dB_phys_i * cos[2*Omega_h*(sqrt(a)-sqrt(a_i))] * sin(k(x-x0))
        #       i.e. physical divergence amplitude is bounded; comoving grows as (a/a_i)^2.
        #       Omega_h = k * c_h / H_0 using the code's run-time c_h.
        #       IC: B0 along z, dBx = A_u*B0*sin(k(x-x0)), psi(a_i)=0, v=0, drho=0.
        # factor_A, factor_G scale the paper-default frequencies:
        #   omega_A = factor_A * pi    (0 = MHD off; 1 = paper default)
        #   omega_G = factor_G * pi/2  (0 = self-gravity off; 1 = paper default)
        # omega_S stays at pi (paper baseline).
        self.drho = 0
        #self.rhozero  = 1.-self.drho-0.000005
        """
        #dimensionless amplitude Au=10**-6
        #wavenumber kL=2pi k=2pi
        #z=127 ai=1/128
        #du(ai)=dA*Va cos(kx) for alfven waves
        #du(ai)=dA*Vs cos(kx) for compressible
        #The travelling waves are initialized using the analytic solutions at
        #a = ai . The amplitudes shown in the figures have been divided by
        #the value of Au such that the velocity amplitude appears as 1 at a =
        #ai
        We use this solution to initialize a simulation with A = π and
        show the amplitude evolution of the wave in Fig. 1.
        #omega_A=pi
        #omega_A=2pi
        show δBc /Bc multiplied with a1/4 and δu/VA
        multiplied with a3/4 as a function of ln a/ai 
        
        omega_S/V_S=k/H_0
        omega_A/V_A=k/H_0
        
        V_A=B_c/sqrt(rho_c)
        V_S=sqrt(gamma*rho_c0/rho_c)
        V_G=sqrt(4piGrho_c)/k
        omega_x=kV_x/H_0
        K=sqrt(omega_A**2-1/16)
        phi=Kln(a/a_i) phi(a_i)=0   
        
        du(x,a_i)=V_A*A_u*cos(phi)*cos(kx)
        du(x,a_i)=H_0/k*omega_A*A_u*cos(kx)
        B_c(a_i)=(H_0/k)*omega_A*sqrt(rho_c)
        The critical density of the Universe in simulation units is 1.
        
        dTheta2          = 0.725
        dHubble0         = 2.894405
        dOmega0          = 1.0
        dLambda          = 0.0
        dRedTo   = 0.0
        bComove          = 1
        """
        dx = 0.5
        wavelength=2*dx
        k=2.*np.pi/wavelength
        deltax = wavelength/nx
        A_u = 10**(-1)

        omega_A = factor_A * np.pi
        omega_S = np.pi
        omega_G = factor_G * np.pi/2
        
        self.gamma=4./3.

        V_A=(self.H0/k)*omega_A
        V_S=(self.H0/k)*omega_S
        V_G=(self.H0/k)*omega_G
        print("V_a ",V_A,"V_S ",V_S,"V_G ",V_G)
        
        sigma=np.sqrt(omega_A**2+omega_S**2-omega_G**2)

        kappa=np.sqrt(sigma**2-1./16.)
        self.kappa=kappa
        kappa_A=np.sqrt(omega_A**2-1./16.)
        self.kappa_A=kappa_A
        alfven = case in (0,1)
        if case == 0:
            """ standing Alfven (paper 5.1.1) """
            a = 1./128.
        elif case == 1:
            """ traveling Alfven (paper 5.1.2), n=2 """
            n = 2
            a = np.exp(-2*n*np.pi/kappa_A)
        elif case == 2:
            """ standing compressible (paper 5.2.1/2/3) """
            a = 1./128.
        elif case == 4:
            """ divergence-cleaning test (this work) """
            a = 1./128.
        else:
            """ traveling compressible gamma=4/3 (paper 5.2.1), n=3 """
            n = 3
            a = np.exp(-2*n*np.pi/kappa)
        self.a=a
        phi=kappa*np.log(1) #log(a/a_i)
        if omega_G > 0:
            self.rhozero=((k*V_G)**2)/(4*np.pi)
        print("rhozero", self.rhozero)
        
        valf=V_A*a**(-0.5)
        Bzero=V_A*np.sqrt(self.rhozero)*a**(-0.5)
        gam1 = self.gamma - 1.
        
        print("Bzero",Bzero)
        
        """ Incompressible """
        #uzero=V_A*A_u*a**(-1)
        #cs=uzero*10
        
        """ Compressible / Alfven velocity amplitude """
        uzero=(V_A if alfven else V_S)*A_u*a**(-1)
        cs=V_S*a**((-3*gam1)/2)
        print('cs',cs,'cs**2',cs**2,'Va',valf)
        print('uzero',uzero)
        
        uuzero = cs**2/(self.gamma*gam1)
        pzero=uuzero*self.rhozero*gam1
        print("uuzeroth",uuzero,"pzero",pzero)
        #uuzero=vA**2  B**2/(rho*self.gamma*gam1)
        dy=dx
        dz=dx
        distribute.setbound(self,-dx,dx,-dy,dy,-dz,dz)
        self.wk = k
        self.x0 = (-dx,-dy,-dz)
        self.xdot0 = -dx
        #distribute.setcloseddist(self,-dx,dx,-dy,dy,-dz,dz,deltax)
        distribute.setcubicdist(self,-dx,dx,-dy,dy,-dz,dz,deltax)
        #dbrho=A_u*omega_S*np.sqrt(a)/(4*sigma**2)*4*kappa
        if case == 3:
            self.drho = A_u*omega_S*np.sqrt(a)/(4*sigma**2)
        else:
            # standing waves (0,2) and traveling Alfven (1) have no density perturbation at a_i
            self.drho = 0.
        print("DRHO",self.drho,"wk",self.wk,"xdot0",self.xdot0,"kappa",kappa)
        if self.drho != 0:
            print("FIXING THE DENSITY HERE")
            xgrid = np.linspace(-dx, dx, 20001)
            self._rho_z_grid = (xgrid, self.getrho(xgrid))
            fixdens.set_density_profile(self,-dx,dx,1)
        

        #Bzero=V_A*np.sqrt(self.rhozero)*a**(-0.5)
        #c_s=V_S*a**(-0.5)
        #P_th=rho*uzero*gam1
        #P_th=rho*cs**2/(gamma)
        #Pmag=B**2/2=rho*vA**2/2          P_th/P_mag = (cs**2/vA**2)*(2/gamma)
        

        self.npart=len(self.x)
        self.ngas=self.npart
        totmass = self.rhozero*self.dxbound*self.dybound*self.dzbound
        self.mass = [totmass/self.npart]*self.npart
        self.h=smth.getsmooth(self,200)
        
        print('npart = ',self.npart,' particle mass = ',self.mass[0], 'wk',self.wk,'x0',self.xdot0)
  
        for i in range(self.npart):
            cosx1 = np.cos(self.wk*(self.x[i]-self.xdot0))
            sinx1 = np.sin(self.wk*(self.x[i]-self.xdot0))
            if alfven:
                # Alfven: dv perpendicular to B0 (B0 along x-hat, dv along y-hat)
                self.vx.append(0.)
                self.vy.append(uzero*cosx1)
                self.vz.append(0.)
            elif case == 4:
                # divergence-cleaning test: no velocity perturbation
                self.vx.append(0.)
                self.vy.append(0.)
                self.vz.append(0.)
            else:
                # compressible: dv parallel to k (along x-hat)
                self.vx.append(uzero*cosx1)
                self.vy.append(0.)
                self.vz.append(0.)
            if case == 0:
                # standing Alfven: B0 along x-hat, dB(a_i)=0
                self.Bx.append(Bzero)
                self.By.append(0.)
                self.Bz.append(0.)
            elif case == 1:
                # traveling Alfven eigenmode (paper eq. 79):
                # dBc/Bc = A_u*sqrt(a_i)/(4*OmegaA) * [sin(kx-psi) - 4*kappa_A*cos(kx-psi)], psi=0 at a_i
                dBy = Bzero*A_u*np.sqrt(a)/(4.*omega_A)*(sinx1 - 4.*kappa_A*cosx1)
                self.Bx.append(Bzero)
                self.By.append(dBy)
                self.Bz.append(0.)
            elif case == 2:
                # standing compressible: B0 along z-hat, dB(a_i)=0
                self.Bx.append(0.)
                self.By.append(0.)
                self.Bz.append(Bzero)
            elif case == 4:
                # divergence-cleaning IC: B0 along z-hat, dBx = A_u*B0*sin(k(x-x0))
                # gives div(Bc) = A_u*B0*k*cos(k(x-x0)) at a=a_i, psi(a_i)=0
                self.Bx.append(Bzero*A_u*sinx1)
                self.By.append(0.)
                self.Bz.append(Bzero)
            else:
                # case 3: traveling compressible gamma=4/3 (existing behavior)
                # uses dBc/Bc = drho/rho simplifying assumption (paper eq. 37)
                self.Bx.append(0.)
                self.By.append(0.)
                self.Bz.append(Bzero+Bzero*(self.rho[i]-self.rhozero)/self.rhozero)
            uui=(pzero+cs**2*(self.rho[i]-self.rhozero))/(gam1*self.rho[i])
            self.u.append(uui)

             
    def getrho(self,x):
        return (self.rhozero + self.drho*(4*self.kappa*np.cos(self.wk*(x - self.xdot0))-np.sin(self.wk*(x - self.xdot0))))

    def getrhoi(self,i):
        return (self.rhozero + self.drho*(4*self.kappa*np.cos(self.wk*(self.x[i] - self.xdot0))-np.sin(self.wk*(self.x[i] - self.xdot0))))

    def getamplitudes():
        return
    
    #------------------------------------------------------
    #+
    #  transform vectors from rotated to unrotated coords
    #+
    #------------------------------------------------------
    def transform_vec(self,xvec,sina,sinb,cosa,cosb):
         x = [0,0,0]
         x[0] = xvec[0]*cosa*cosb - xvec[1]*sinb - xvec[2]*sina*cosb
         x[1] = xvec[0]*cosa*sinb + xvec[1]*cosb - xvec[2]*sina*sinb
         x[2] = xvec[0]*sina + xvec[2]*cosa
         return x
