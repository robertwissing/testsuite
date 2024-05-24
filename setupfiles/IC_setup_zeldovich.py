## Setup for ORZAG-TANG
import numpy as np
import random 
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
class setup_zeldovich(object):
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
    cosmo=1
    gamma=1.666666667
    #molweight=1.0
    molweight=1.0
    #dMeanMolWeight  = 0.59259
    shape=2
    periodic=1
    grav=1
    deltastep=0.002
    nsteps=400
    dmsolunit=1.0
    dkpcunit=1000.0
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
    A=0.0
    A2=0.0
    k=0.0
    rcylmax2= 10**22
    rmax2= 10**22
    rhozero=1.0
    stratified=0
    def __init__(self):
        pass
    def create(self,nx,Lboxinkpc,zi,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000'):
         self.dkpcunit=Lboxinkpc
         ascale=1./(1.+zi)
         #--setup parameters
         self.rhozero = 1.0

         przero=6.2*10**(-8)

         beta=10**-6
         L=1.0
         zc=2.0
         #0.09523809523 
         #0.09523809523 * 4.58257569496 
         #A2 is 0.43643578047
         self.A=(1.+zc)/(1.+zi)
         self.A2=self.A*np.sqrt((1.+zi))
         A3peak=0.65/(1.0+zi)
         gam1=self.gamma-1
         Ti=1000.0*self.A*self.A
         przero=gam1*Ti/394471218
         xmin=-L/2.0
         ymin=xmin
         zmin=xmin
         self.k=2*np.pi/L
         cs=np.sqrt(self.gamma)
         H0=2.894405
         B0=np.sqrt(2.0*przero/beta)
         uzero=przero/(gam1*self.rhozero)
         print("pressure ",przero, "beta ",beta ,'B0 ',B0,'uzero',uzero)
         distribute.setbound(self,xmin,-xmin,ymin,-ymin,zmin,-zmin)
         deltax = self.dxbound/nx
         distribute.setcubicdist(self,xmin,-xmin,ymin,-ymin,zmin,-zmin,deltax)
         
         self.npart=len(self.x)
         totmass = self.rhozero*self.dxbound*self.dybound*self.dzbound
         print('npart = ',self.npart)
         print("totmass ",totmass)
         self.mass = [totmass/self.npart]*self.npart
         print('npart = ',self.npart, 'mass = ', self.mass[1])
         self.npart=len(self.x)
         self.ngas=self.npart
         totmass = self.rhozero*self.dxbound*self.dybound*self.dzbound
         self.mass = [totmass/self.npart]*self.npart
         print('npart = ',self.npart, 'mass = ', self.mass[1])
         #self.h=smth.getsmooth(self,64)
         self.h=0.01
         print(self.h/deltax)
         for i in range(self.npart):
                     self.Bx.append(0.)
                     self.By.append(B0)
                     self.Bz.append(0.) 
                     #self.rho[i]=self.rhozero/(1.-self.A*np.cos(self.k*self.x[i]))
                     #self.u.append(Ti*(self.rho[i]/self.rhozero)**(2./3.))
                     self.u.append(uzero)
                     #self.vx.append(A3peak*np.sin(self.k*self.x[i]))
                     self.vx.append(-H0*self.A2*np.sin(self.k*self.x[i])/self.k) #used before
                     #self.vx.append(-H0*np.sin(self.k*self.x[i])/self.k)
                     #self.vx.append(-(1+zi)*H0*np.sin(self.k*self.x[i])/self.k)
                     #self.vx.append(0.0)
                     self.vy.append(0.0)
                     self.vz.append(0.0)
                     self.x[i]=self.x[i]-self.A*np.sin(self.k*self.x[i])/self.k #used before
             
        
    def getrhoi(self,i):
        if(self.stratified==1):
            return self.rhozero/(1-self.A*np.cos(self.k*self.x[i]))
        else:
            return self.rhozero
        
    def getrho(self,x):
        if(self.stratified==1):
            return self.rhozero/(1-self.A*np.cos(self.k*x))
        else:
            return self.rhozero