## Setup for ORZAG-TANG
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
import readtipsy as tip
class setup_alfven(object):
    dICdensRsmooth = 0.1
    dICdensprofile = 1
    dICdensdir = 1
    dICdensR = 1.0 
    dICdensinner = 1.0 # calculated depending on mass
    dICdensouter = 1.0
    rhoit=0
    ns=64
    npart=0
    ngas=0
    ndark=0
    nstar=0
    ndim=3
    time=0
    gamma=5./3.
    periodic=1
    deltastep=0.01
    nsteps=1000
    freqout=200
    dmsolunit=1.0
    dkpcunit=1.0
    molweight=1
    cosmo=0
    grav=0
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
    wk=0.
    x0=0.
    xdot0=0.
    rhozero=1.
    drho=0.
    shape=0.
    
    def __init__(self):
        pass
    def create(self,nx,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000',rotate=1):
        if(rotate==1):           
            sina = 2./3.
            sinb = 2./np.sqrt(5.)
            cosa = np.sqrt(5.)/3.
            cosb = 1./np.sqrt(5.)
            igeom = "igeom_rotated"
        else:
            sina = 0.
            sinb = 0.
            cosa = 1.
            cosb = 1.
            igeom = "igeom_cartesian"
        
        runit = (cosa*cosb,cosa*sinb,sina)
        self.rhozero  = 1.
        rvec     = 0.
        du = 0.
        self.drho = 0.
        dv = 0.
        dB = 0.
        gam1 = self.gamma - 1.
        ampl   = 0.1
        przero = 0.1
        Bzero  = 1.0 #in x
        vzero  = 0.
        uuzero = przero/(gam1*self.rhozero)
        wavelength = 1.
        vwave = 1.
    
        #prim_to_cons(rhozero,vzero,Bzero,uuzero,q0)
    
        print('Setup for 3D Alfenwave problem..')
        print('LATTICE SETUP REQUIRES rho=Kh relation otherwise will have small oscillations')
        print('sina= ',sina,', sinb= ',sinb, 'density= ', self.rhozero, 'B= ', Bzero, 'pressure = ',przero, 'period =', wavelength/vwave) 
    
        dx=1.5
        dy=0.75
        dz=0.75
        distribute.setbound(self,-dx,dx,-dy,dy,-dz,dz)
        self.wk = 2.*np.pi/wavelength
        self.x0 = (-dx,-dy,-dz)
        self.xdot0 = np.dot(self.x0,runit)
        deltax=self.dxbound/nx
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
            
        if self.drho != 0:   
            fixdens.set_density_profile(self,-dx,dx,1)
        self.npart=len(self.x)
        self.ngas=self.npart
        totmass = self.rhozero*self.dxbound*self.dybound*self.dzbound
        self.mass = [totmass/self.npart]*self.npart
        if vm==1:
            print("Varying masses -> ONLY FOR RAND + RELAXING TO GLASS");              
            for i in range(0,self.npart,2):
                dmrat=0.25*self.mass[i]*np.random.rand()
                self.mass[i]=self.mass[i]+dmrat
                self.mass[i+1]=self.mass[i+1]-dmrat

        
        self.h=smth.getsmooth(self,200)
        print('npart = ',self.npart,' particle mass = ',self.mass[0])
  
        for i in range(self.npart):
            xyz=(self.x[i],self.y[i],self.z[i])
            x1    = np.dot(xyz,runit)
            sinx1 = np.sin(self.wk*(x1-self.xdot0))
            cosx1 = np.cos(self.wk*(x1-self.xdot0))
            vvec = ampl*np.array([0.,sinx1,cosx1])
            Bvec = ampl*np.array([Bzero/ampl,sinx1,cosx1])
            #Avec =  np.array([0.0,ampl*sinx1/self.wk,-ampl*cosx1/self.wk])
            #Bvec = Avec
            uui  = uuzero
            vxyz=self.transform_vec(vvec,sina,sinb,cosa,cosb)
            Bxyz=self.transform_vec(Bvec,sina,sinb,cosa,cosb)
            self.vx.append(vxyz[0])
            self.vy.append(vxyz[1])
            self.vz.append(vxyz[2])
            self.Bx.append(Bxyz[0])
            self.By.append(Bxyz[1])
            self.Bz.append(Bxyz[2])
            self.u.append(uuzero)
             
    def getrho(self,x):
        return (self.rhozero + self.drho*np.cos(self.wk*(x - self.xdot0)))
    def getrhoi(self,i):
        return (self.rhozero + self.drho*np.cos(self.wk*(self.x[i] - self.xdot0)))   
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
