## Setup for ORZAG-TANG
import numpy as np
import IC_distribute as distribute
import IC_fixdens as fixdens
import IC_smoothlength as smth
from numba import jit
class setup_mhdslabcond(object):
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
    periodic=1
    deltastep=0.01
    nsteps=200
    dmsolunit=1.0
    dkpcunit=1.0
    adi=1
    xshock=0.0
    shape=0
    tmax=0.2
    
    def __init__(self):
        pass
    def create(self,nx,choice):
        #--general parameters
         ti=2
        #Get shock parameters
         [xleft,xright]=self.choose_shock(choice)
         dxleft = -xleft/nx
         self.xshock = 0.5*(xleft + xright)
    
         #
         # adjust boundaries to allow space for boundary particles and inflow
         #
         [xleft,xright,ymin,ymax,zmin,zmax,dxright]=self.adjust_shock_boundaries(xleft,xright,dxleft)
        
        
         #3 set domain boundaries - use a longer box in x direction since
         # it is not actually periodic in x
         distribute.setbound(self,ti*xleft,ti*xright,ymin,ymax,zmin,zmax)
         
         # set limits of the different domains
         xminleft  = np.array([ti*xleft,ymin,zmin])
         xmaxleft  = np.array([ti*xright,ymax,zmax])
         xminright = np.array([ti*xleft,ymin,zmin])
         xmaxright = np.array([ti*xright,ymax,zmax])
         print(xminright)
         
         # setup the particles
         
#         if (np.abs(self.leftstate[0]-self.rightstate[0]) > np.finfo(float).eps):
             
             

         #then divide the x axis into two halves at xshock
         if(1):
            xmaxleft[0]  = self.xshock
            xminright[0] = self.xshock
            if(abs(xminleft[0])>abs(xmaxright[0])):
                xmaxright[0]=-xminleft[0]
            else:
                xminleft[0]=-xmaxright[0]
            #write(iprint,'(1x,3(a,es16.8))') 'Setup_shock: left half  ',xminleft(1), ' to ',xmaxleft(1), ' with dx_left  = ',dxleft
            print(dxleft*nx,xmaxleft[0]-xminleft[0],dxleft,nx)
            # set up a uniform lattice
            setcloseddist=distribute.setcloseddist
            #setcloseddist=jit(distribute.setcloseddist)
            setcloseddist(self,xminleft[0],xmaxleft[0],xminleft[1],xmaxleft[1],xminleft[2],xmaxleft[2],dxleft)  # set left half
        
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
            # now set up box for right half
            #write(iprint,'(1x,3(a,es16.8))') 'Setup_shock: right half ',xminright(1),' to ',xmaxright(1),' with dx_right = ',dxright
            setcloseddist(self,xminright[0], xmaxright[0], xminright[1], xmaxright[1],xminright[2], xmaxright[2], dxright) # set right half
            self.npart=len(self.x)
            self.ngas=self.npart
            self.mass = [(totmass+totmass2)/self.npart]*self.npart
            self.h=smth.getsmooth(self,200)
            # define rhozero as average density; required for certain simulations (e.g. non-ideal MHD with constant resistivity)
            rhozero = (self.leftstate[0]*np.prod(xmaxleft-xminleft) + self.rightstate[0]*np.prod(xmaxright - xminright)) + np.prod(xmaxright - xminleft)
 #        else:
            
            # set all of volume if densities are equal
            #write(iprint,'(3(a,es16.8))') 'Setup_shock: one density  ',xminleft(1), ' to ',xmaxright(1), ' with dx  = ',dxleft
 #           distribute.setcloseddist(self,xminleft[0],xmaxleft[0],xminleft[1], xmaxleft[1],xminleft[2],xmaxleft[2],dxleft)
  #          self.npart=len(self.x)
   #         volume           = np.inner(xmaxleft-xminleft,xmaxright-xminright)
    #        rhozero          = self.leftstate[0]
     #       dxright          = dxleft
      #      totmass = rhozero*volume
       #     self.npart=len(self.x)
        #    self.ngas=self.npart
         #   self.mass = [totmass/self.npart]*self.npart
          #  self.h=smth.getsmooth(self,200)
         
         gam1 = self.gamma - 1.
         print(gam1)      
        
         for i in range(self.npart):
            delta = self.x[i] - self.xshock
            lim=1.0
            if (delta > 0.):
                if(self.x[i]>lim):
                   self.vx.append(self.rightstate[2]*np.exp(-8*(self.x[i]-lim)))
                   self.vy.append(self.rightstate[3]*np.exp(-8*(self.x[i]-lim)))
                   self.vz.append(self.rightstate[4]*np.exp(-8*(self.x[i]-lim)))
                   #print(self.vx[i])
                else:
                   self.vx.append(self.rightstate[2])
                   self.vy.append(self.rightstate[3])
                   self.vz.append(self.rightstate[4])
                self.u.append(self.rightstate[1])
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
                self.u.append(self.leftstate[1])
                self.Bx.append(self.leftstate[5])
                self.By.append(self.leftstate[6])
                self.Bz.append(self.leftstate[7])
        
        
        #-----------------------------------------------------------------------
        #+
        #  Adjust the shock boundaries to allow for inflow/outflow
        #+
        #-----------------------------------------------------------------------
    def adjust_shock_boundaries(self,xleft,xright,dxleft):
            vxleft=self.leftstate[2]
            vxright=self.rightstate[2]
            densleft=self.leftstate[0]
            densright=self.rightstate[0]
            tmax=self.tmax
            #if (vxleft > 0):
            #    xleft = xleft - vxleft*tmax
            dxright = dxleft*(densleft/densright)**(1./3) # NB: dxright here is only approximate
            #if (vxright < 0):
            #    xright = xright - vxright*tmax
             # try to give y boundary that is a multiple of 6 particle spacings in the low density part
            fac = -6*(int(1.99*2./6.) + 1)*max(dxleft,dxright)
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
            #nx     = 256
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
            #[dens,u,vx,vy,vz,bx,by,bz]
            if(choice==1):
                #--Isotropic
                #shocktype = "Isotropic"
                self.gamma      = 5./3.
                self.leftstate  = [1.,1.0,0.,0.,0.,0.,0.,0.]
                self.rightstate = [1.,2.0,0.,0.,0.,0.,0.,0.]
            if(choice==2):
                #--Anisotropic parallel
                #shocktype = "Anisotropic parallel"
                self.gamma      = 5./3.
                self.leftstate  = [1.,1.0,0.,0.,0.,1.0,0.,0.]
                self.rightstate = [1.,2.0,0.,0.,0.,1.0,0.,0.]
            if(choice==3):
                #--Anisotropic 1/sqrt2 angle
                #shocktype = "Anisotropic"
                self.gamma      = 5./3.
                self.leftstate  = [1.,1.0,0.,0.,0.,1.0,1.0,0.]
                self.rightstate = [1.,2.0,0.,0.,0.,1.0,1.0,0.]
            if(choice==4):
                #--Anisotropic transverse
                #shocktype = "Anisotropic transverse"
                self.gamma      = 5./3.
                self.leftstate  = [1.,1.0,0.,0.,0.,0.0,1.0,0.]
                self.rightstate = [1.,2.0,0.,0.,0.,0.0,1.0,0.]
            if(choice==5):
                #--Anisotropic small angl
                #shocktype = "Anisotropic transverse"
                self.gamma      = 5./3.
                self.leftstate  = [1.,1.0,0.,0.,0.,0.1,0.9,0.]
                self.rightstate = [1.,2.0,0.,0.,0.,0.1,0.9,0.]
               
        
            if (abs(xright)  < np.finfo(float).eps):
                  xright  = -xleft
            return xleft,xright
         
    def getrhoi(self,i):
        
        if (self.x[i]<self.xshock):
            return self.leftstate[0]
        else:
            return self.rightstate[0]
             