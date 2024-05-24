import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
import readtipsy as tip
class setup_sedov(object):
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
    rhoit=0
    ndim=3
    ns=64
    time=0
    grav=0
    cosmo=0
    molweight=1.0
    gamma=1.666666667
    periodic=1
    deltastep=0.0002
    freqout=20
    nsteps=200
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
    
    def __init__(self):
        pass


    #Select 0 we define and inner and outer pressure and an inner radius of blast. And corresponding inner plasmabeta for mhd.
    #Select 1 and 2 we choose an blast energy and give it to a few particles in centre.
    #In select 1 we have different magnetic strength in the blast and outside blast(need divergence cleaning) and select 2 we have a uniform field.
    def create(self,nx,distri=0,vm=0,entry='glassreadyfile',select=1,betain=0,ublast=10.0,Rin=0.125,Pin=100,Pout=1):
        
         gam1=self.gamma-1
         self.rhozero=1.0
         #--setup parameters
         if select==1 or select==2:
             print("You picked select 1 or 2")
             Bin=Bout=0.0
         else:
             print("You picked select 0")
             if betain==0:
                 Bin=Bout=0.0
             else:
                 Bin=np.sqrt(2*Pin/betain)
                 Bout=Bin
         
         #boundaries    
         dz=0.5
         dy=0.5
         dx=0.5
         distribute.setbound(self,-dx,dx,-dy,dy,-dz,dz)
         deltax = self.dxbound/nx
         print('DISTRI: ',distri)
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
            self.mass=tgdata[:,0]
         self.npart=len(self.x)
         self.ngas=self.npart
         if distri==0 or distri==1:
             totmass = self.rhozero*self.dxbound*self.dybound*self.dzbound
             self.mass = [totmass/self.npart]*self.npart


         print('npart = ',self.npart,' particle mass = ',self.mass[1])
         self.h=smth.getsmooth(self,64)

         if select == 1 or select == 2:
             ublast=ublast/self.mass[0]
             Rin=np.mean(self.h)
         print('Rin',Rin)
         print(self.h/deltax)
         n=0
         n2=0
         sum1=0
         sum2=0
         Uin=0.0
         Uout=0.0
         for i in range(self.npart):
             radius2 = self.x[i]**2 + self.y[i]**2+self.z[i]**2
             if (radius2 < Rin**2):
                 n=n+1
             else:
                 n2=n2+1
         if select == 1 or select == 2:
             Uin=ublast/n
             Uout=2*10**-5/n2
             Pin=Uin*(gam1*self.rhozero)
             Pout=Uout*(gam1*self.rhozero)
             if betain==0:
                 Bzero=0.0
             else:
                 if select == 1:
                     Emagin=Uin/betain
                     Emagout=Uout/betain
                     Bin=np.sqrt(2*Emagin*self.rhozero)
                     Bout=np.sqrt(2*Emagout*self.rhozero)
                 else:
                     Bin=np.sqrt(2*Pin/(betain))
                     Bout=np.sqrt(2*Pin/(betain))

         print("Bin",Bin,"Bout",Bout);
            
         for i in range(self.npart):
             radius2 = self.x[i]**2 + self.y[i]**2+self.z[i]**2
             self.vx.append(0.)
             self.vy.append(0.)
             self.vz.append(0.)
             if (radius2 < Rin**2):
                 sum1=sum1+Uin
                 self.u.append(Pin/(gam1*self.rhozero))
                 Bzero=Bin
             else:
                 sum2=sum2+Uout
                 self.u.append(Pout/(gam1*self.rhozero))
                 Bzero=Bout
             self.Bx.append(Bzero*np.sqrt(0.5))
             self.By.append(0.)
             self.Bz.append(Bzero*np.sqrt(0.5))
         print("ETH: ", np.sum(self.u),"inside ",sum1,"outside ",sum2," ",n,n2)

    def getrhoi(self,i):
        return self.rhozero
