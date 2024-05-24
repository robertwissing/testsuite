## Setup for ORZAG-TANG
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
from numba import jit
import readtipsy as tip
class setup_shocktube(object):
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
    periodic=1
    freqout=10
    ns=64
    deltastep=0.01
    nsteps=200
    dmsolunit=1.0
    dkpcunit=1.0
    adi=1
    grav=0
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
    dxbound=[]
    dybound=[]
    dzbound=[]
    totvol=[]
    rcylmin2= 0
    rmin2= 0
    rcylmax2= 10**22
    rmax2= 10**22
    rightstate=[]
    leftstate=[]
    polyk=0.1
    gamma=5./3.
    xshock=0.0
    shape=0
    tmax=0.2
    extend=2.0
    
    def __init__(self):
        pass
    def create(self,nx,distri=0,vm=0,entry='datafiles/alfvenwave128_preglass.00000',choice=1):
        #Get shock parameters
         [xleft,xright]=self.choose_shock(choice)
         dxleft = -xleft/nx
         self.xshock = 0.5*(xleft + xright)

         self.dICdensprofile = 2
         self.dICdensdir = 6
         self.dICdensR = self.xshock
         self.dICdensinner = self.leftstate[0]
         self.dICdensouter = self.rightstate[0]
    
         # adjust boundaries to allow space for boundary particles and inflow
         [xleft,xright,ymin,ymax,zmin,zmax,dxright]=self.adjust_shock_boundaries(xleft,xright,dxleft)
        
        
         #3 set domain boundaries - use a longer box in x direction since
         # it is not actually periodic in x
         distribute.setbound(self,self.extend*xleft,self.extend*xright,ymin,ymax,zmin,zmax)
         
         # set limits of the different domains
         xminleft  = np.array([self.extend*xleft,ymin,zmin])
         xmaxleft  = np.array([self.extend*xright,ymax,zmax])
         xminright = np.array([self.extend*xleft,ymin,zmin])
         xmaxright = np.array([self.extend*xright,ymax,zmax])  
             
         #then divide the x axis into two halves at xshock
         xmaxleft[0]  = self.xshock
         xminright[0] = self.xshock
         if(abs(xminleft[0])>abs(xmaxright[0])):
             xmaxright[0]=-xminleft[0]
         else:
             xminleft[0]=-xmaxright[0]
         print(dxleft*nx,xmaxleft[0]-xminleft[0],dxleft,nx)


         print("distri",distri);
         if distri==0:
             distrifunc=distribute.setcloseddist;
         elif distri==1:
             distrifunc=distribute.setrandomdist;

         if distri==0 or distri==1:
             distrifunc(self,xminleft[0],xmaxleft[0],xminleft[1],xmaxleft[1],xminleft[2],xmaxleft[2],dxleft)  # set left half
             # set particle mass
             volume           = np.prod(xmaxleft-xminleft)
             totmass          = volume*self.leftstate[0]
             npart=len(self.x)
             masspart = totmass/npart
             # (closepacked) now adjust spacing on right hand side to get correct density given the particle mass
             volume  = np.prod(xmaxright - xminright)
             totmass2 = volume*self.rightstate[0]
             [ny,nz]=distribute.get_ny_nz_closepacked(dxright,xminright[1],xmaxright[1],xminright[2],xmaxright[2])
             dxright = (xmaxright[0] - xminright[0])/((totmass2/masspart)/(ny*nz))
             distrifunc(self,xminright[0], xmaxright[0], xminright[1], xmaxright[1],xminright[2], xmaxright[2], dxright) # set right half
             self.npart=len(self.x)
             self.ngas=self.npart
             self.mass = [(totmass+totmass2)/self.npart]*self.npart
             if vm==1:
                 print("Varying masses -> ONLY FOR RAND + RELAXING TO GLASS");
                 for i in range(0,self.npart,2):
                     dmrat=0.25*self.mass[i]*np.random.rand()
                     self.mass[i]=self.mass[i]+dmrat
                     self.mass[i+1]=self.mass[i+1]-dmrat
            
         if distri==2:
             tgdata,tddata,tsdata,data_header,time=tip.readtipsy(entry);
             self.x=tgdata[:,1]
             self.y=tgdata[:,2]
             self.z=tgdata[:,3]
             self.rho=tgdata[:,7]
             self.mass=tgdata[:,0]
             self.npart=len(self.x)
             self.ngas=self.npart
             
         self.h=smth.getsmooth(self,200)
         rhozero = (self.leftstate[0]*np.prod(xmaxleft-xminleft) + self.rightstate[0]*np.prod(xmaxright - xminright)) + np.prod(xmaxright - xminleft)
         
         gam1 = self.gamma - 1.
         if (np.abs(gam1) > 1.e-3):
            uuleft  = self.leftstate[1]/(gam1*self.leftstate[0])
            uuright = self.rightstate[1]/(gam1*self.rightstate[0])
         else:
            uuleft  = 3.*self.leftstate[1]/(2.*self.leftstate[0])
            uuright = 3.*self.rightstate[1]/(2.*self.rightstate[0])
        
        
         for i in range(self.npart):
            delta = self.x[i] - self.xshock
            lim=self.dxbound*0.4
            if (delta > 0.):
                if(self.x[i]>lim):
                   self.vx.append(self.rightstate[2]*np.exp(-8*(self.x[i]-lim)))
                   self.vy.append(self.rightstate[3]*np.exp(-8*(self.x[i]-lim)))
                   self.vz.append(self.rightstate[4]*np.exp(-8*(self.x[i]-lim)))
                else:
                   self.vx.append(self.rightstate[2])
                   self.vy.append(self.rightstate[3])
                   self.vz.append(self.rightstate[4])
                self.u.append(uuright)
                self.Bx.append(self.rightstate[5])
                self.By.append(self.rightstate[6])
                self.Bz.append(self.rightstate[7])
            else:
                if(self.x[i]<-lim):
                   self.vx.append(self.leftstate[2]*np.exp(8*(lim+self.x[i])))
                   self.vy.append(self.leftstate[3]*np.exp(8*(lim+self.x[i])))
                   self.vz.append(self.leftstate[4]*np.exp(8*(lim+self.x[i])))
                else:
                   self.vx.append(self.leftstate[2])
                   self.vy.append(self.leftstate[3])
                   self.vz.append(self.leftstate[4])
                self.u.append(uuleft)
                self.Bx.append(self.leftstate[5])
                self.By.append(self.leftstate[6])
                self.Bz.append(self.leftstate[7])
        
        
    def adjust_shock_boundaries(self,xleft,xright,dxleft):
            vxleft=self.leftstate[2]
            vxright=self.rightstate[2]
            densleft=self.leftstate[0]
            densright=self.rightstate[0]
            tmax=self.tmax
            dxright = dxleft*(densleft/densright)**(1./3) # NB: dxright here is only approximate
             # try to give y boundary that is a multiple of 6 particle spacings in the low density part
            fac = -8*(int(1.99*2./6.) + 1)*max(dxleft,dxright)
            print("FACTOR",fac)
            ymin   = fac*np.sqrt(0.75)
            zmin   = fac*np.sqrt(6.)/3.
            ymax  = -ymin
            zmax  = -zmin
            return xleft,xright,ymin,ymax,zmin,zmax,dxright
        
        #-----------------------------------------------------------------------
        #+
        #  Choose which shock tube problem to set up
        #+
        #-----------------------------------------------------------------------
    def choose_shock (self,choice):
        
        #--set default file output parameters
            xleft  = -0.500
            yleft  = 0.0
            zleft  = 0.0
            xright = 0.0
            yright = 0.0
            zright = 0.0
            const  = np.sqrt(4.*np.pi)
            
            #--list of shocks
            
            # shocks(:) = 'none'
            # shocks(1) = 'Sod shock'
            # shocks(2) = 'Ryu 1a'
            # shocks(3) = 'Ryu 1b'
            # shocks(4) = 'Ryu 2a'
            # shocks(5) = 'Ryu 2b'
            # shocks(6) = 'Brio-Wu (Ryu 5a)'
            # shocks(7) = 'C-shock'
            # shocks(8) = 'Steady shock'
            #[dens,press,vx,vy,vz,bx,by,bz]
            if(choice==1):
                #--Sod shock
                shocktype = "Sod shock"
                self.extend = 3.0
                self.gamma      = 5./3.
                self.deltastep=0.002
                self.nsteps=100
                self.leftstate  = [1.000,1.0,0.,0.,0.,0.,0.,0.]
                self.rightstate = [0.125,0.1,0.,0.,0.,0.,0.,0.]
            elif(choice==2):
                #--Ryu et al. shock 1a
                shocktype = "Ryu et al. shock 1a"
                self.extend = 5.0
                self.deltastep=0.0008
                self.nsteps=100
                self.gamma      =   5./3.
                self.leftstate  = [1.,20.,10.,0.,0.,5./const,5./const,0.]
                self.rightstate = [1.,1.,-10.,0.,0.,5./const,5./const,0.]
            elif(choice==3):
                #--Ryu et al. shock 1b
                shocktype = "Ryu et al. shock 1b"
                self.extend = 2.0
                self.gamma      =  5./3.
                self.deltastep=0.0003
                self.nsteps=200
                self.leftstate  = [1.0,1. ,0.,0.,0.,5./const,5./const,0.]
                self.rightstate = [0.1,10.,0.,0.,0.,5./const,2./const,0.]
            elif(choice==4):
                #--Ryu et al. shock 2a
                shocktype  = "Ryu et al. shock 2a"
                self.extend = 2.0
                self.deltastep=0.002
                self.nsteps=100
                self.gamma      = 5./3.
                self.leftstate  = [1.08,0.95,1.2,0.01,0.5,2./const,3.6/const,2./const]
                self.rightstate = [1.  ,1.  ,0. ,0.  ,0. ,2./const,4.0/const,2./const]
            elif(choice==5):
                #--Ryu et al. shock 2b
                shocktype = "Ryu et al. shock 2b"
                self.extend = 2.0
                self.gamma      = 5./3.
                self.deltastep=0.00035
                self.nsteps=100
                self.leftstate  = [1.0,1. ,0.,0.,0.,3./const,6./const,0.]
                self.rightstate = [0.1,10.,0.,2.,1.,3./const,1./const,0.]
            elif(choice==6):
                #--Brio-Wu shock
                shocktype = "Brio/Wu (Ryu/Jones shock 5a)"
                self.extend = 2.0
                self.gamma      = 2.0
                self.deltastep=0.001
                self.nsteps=100
                self.leftstate  = [1.000,1.0,0.,0.,0.,0.75, 1.,0.]
                self.rightstate = [0.125,0.1,0.,0.,0.,0.75,-1.,0.]
            elif(choice==7):
                #--C-shock
                shocktype = "C-shock"
                self.extend = 2.0
                self.deltastep=0.0008
                self.nsteps=1000
                self.gamma      =  1.0
                self.polyk      =  0.01
                self.leftstate  = [1.,0.006, 4.45,0.,0.,1./sqrt(2.),1./sqrt(2.),0.]
                self.rightstate = [1.,0.006,-4.45,0.,0.,1./sqrt(2.),1./sqrt(2.),0.]
                xleft      = -4.00
            elif(choice==8):
                #--Steady shock (Falle 2003)
                #nx         = 512
                self.extend = 2.0
                self.deltastep=0.0008
                self.nsteps=1000
                self.polyk      = 0.01
                self.gamma      =  1.0
                self.leftstate  = [1.7942,0.017942,-0.9759,-0.6561,0.,1.,1.74885,0.]
                self.rightstate = [1.    ,0.01    ,-1.7510, 0.    ,0.,1.,0.6    ,0.]
                xleft      = -2.0
                #Toth 2000 shocktube
            elif(choice==9):
                self.extend = 2.0
                self.gamma      = 5./3.
                self.tmax = 0.08
                self.deltastep=0.0008
                self.nsteps=100
                self.leftstate  = [1.,20.,10.,0.,0.,5./const,5/const, 0.]
                self.rightstate = [1.,1.,-10.,0.,0.,5./const,5./const, 0.]
                #MHD discontinuity shocktube
            elif(choice==10):
                self.extend = 2.0
                self.gamma      = 5./3.
                self.leftstate  = [1.,1.,0.,0.,0.,20./const,20./const, 0.]
                self.rightstate = [1.,1.,0.,0.,0.,5./const,5./const, 0.]
            elif(choice==11):
                self.extend = 2.0
                self.gamma = 1.4
                self.leftstate  = [1.,1000.,0.,0.,0.,0.,0., 0.]
                self.rightstate = [1.,0.1,0.,0.,0.,0.,0., 0.]
                
            if (abs(xright)  < np.finfo(float).eps):
                  xright  = -xleft
            return xleft,xright
         
    def getrhoi(self,i):
        
        if (self.x[i]<self.xshock):
            return self.leftstate[0]
        else:
            return self.rightstate[0]
                           
