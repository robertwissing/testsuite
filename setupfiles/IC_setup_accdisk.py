## Setup for ORZAG-TANG
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
import readtipsy as tip

class setup_accdisk(object):
    dICdensRsmooth = 0.1
    dICdensprofile = 7
    dICdensdir = 5
    dICdensR = 10.0 
    dICdensinner = 360.0 # calculated depending on mass
    dICdensouter = 1E-12
    rhoit=0
    ns=58
    npart=0
    ngas=0
    ndark=0
    nstar=0
    ndim=3
    time=0
    gamma=5./3.
    periodic=1
    deltastep=0.62831853071  #1/10 of an orbit
    freqout=2
    nsteps=1000 #100 orbits
    dmsolunit=1.0
    dkpcunit=4.84814E-9
    adi=0
    grav=1
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
    dxbound=[]
    dybound=[]
    dzbound=[]
    totvol=[]
    rcylmin2= 0
    rmin2= 0
    rcylmax2= 10**22
    rmax2= 10**22
    shape=2
    rhozero = 1.0
    massdisk= 1E-6
    massrat = 1E-3
    rdisk = 10.0
    q = 0.5
    hr = 0.1
    rorb = 1.0
    onlyoneBH = 0
    
    def __init__(self):
        pass
    def create(self,nx,distri=0,vm=0,entry=' ',hr=0.1,massrat=1E-3,onlyoneBH=0,q=0.5,rorb=1.0,rdisk=10.0,rinner=0.5):
        
         gam1=self.gamma-1
         self.onlyoneBH=onlyoneBH
         self.dICdensRsmooth = hr
         self.dICdensR = rdisk;
         self.hr = hr;
         self.massrat = massrat;
         self.q = q;
         self.rorb = rorb;
         self.rdisk = rdisk;
         if(q < 2.0):
             self.dICdensinner = self.massdisk/2*pi*(rdisk**(2-q)-rinner**(2-q))/(2-q); #calculate normalization density from simple surface density R^(-qdisk)  
         elif(q==2.0):
             self.dICdensinner =self.massdisk/2*pi*(log(rdisk)-log(rinner));
         else:
             self.dICdensinner = self.massdisk/2*pi*(rinner**(2-q)-rdisk**(2-q))/(q-2);
         #--setup parameters
         #self.masszero=1.0
         #cs0 is same as scaleheight h/r and = 1/mach vorb=r**(-q)
         #P=cs0*density*r**(-q)
         #just cube, icgenerator will do rest
         dz=1.0
         dy=self.rdisk
         dx=self.rdisk
         distribute.setbound(self,-dx,dx,-dy,dy,-dz,dz)
         deltax = self.dxbound/nx
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
            self.periodic=0
         self.npart=len(self.x)
         self.ngas=self.npart
         totmass = self.massdisk
         print('npart = ',self.npart,' total disk mass = ',totmass)  
         self.mass = [totmass/self.npart]*self.npart
         print(' particle mass = ',self.mass[1])
         self.dxbound = self.rdisk*4
         self.dybound =	self.rdisk*4
         self.dzbound =	self.rdisk*2
         
         for i in range(self.npart):
             rad = np.sqrt(self.x[i]**2 + self.y[i]**2+self.z[i]**2)
             vxi,vyi,vzi = self.getveli(i)
             self.vx.append(vxi)
             self.vy.append(vyi)
             self.vz.append(vzi)
             self.u.append(self.hr*rad**(-self.q)/(gam1))
             self.Bx.append(0.)
             self.By.append(0.)
             self.Bz.append(0.)
             self.h.append(hr*0.6*1E-6) ## 1E-6 times smaller than secondary BH
             
    def getrhoi(self,i):
        return self.rhozero
    
    def getveli(self,i):
        rcyl=np.sqrt(self.x[i]*self.x[i]+self.y[i]*self.y[i])
        vazi=rcyl**(-self.q)
        vx,vy,vz = self.transform_coordinate(0.0,vazi,0.0,self.x[i],self.y[i],self.z[i],"cyl","cart")
        return vx,vy,vz
    
    def transform_coordinate(self,v1,v2,v3,x,y,z,fromcord,tocord):
        
        if(fromcord=="cyl"):
            r=np.sqrt(x**2+y**2)
            if(tocord=="cart"):
                ux=v1*x/r-v2*y/r
                uy=v1*y/r+v2*x/r 
                uz=v3
                
        return ux,uy,uz
