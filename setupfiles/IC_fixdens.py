import numpy as np
from numba import njit

ngrid = 2048 #r number of points used when integrating rho to get mass
maxits = 500  # max number of iterations
maxits_nr =200  # max iterations with Newton-Raphson
tol = 1.e-12  # tolerance on iterations
igeom_cartesian = 1
igeom_cylindrical = 2
igeom_spherical = 3

#--------------------------------------------------------------------------
#+
#  Subroutine to implement the stretch mapping procedure
#
#  IN/OUT:
#    xyzh : particle coordinates and smoothing length
#
#  IN:
  #  np : number of particles
#    cmin, cmax         : range in the coordinate to apply transformation
#  start (optional)   : only consider particles between start and np
#    geom  (optional)   : geometry in which stretch mapping is to be performed
#                         1 - cartesian
#                          2 - cylindrical
       #                   3 - spherical
      #                    4 - toroidal
     #                    (if not specified, assumed to be cartesian)
    #coord (optional)   : coordinate direction in which stretch mapping is to be performed
         #                (if not specified, assumed to be the first coordinate)
    #rhofunc (optional) : function containing the desired density function rho(r) or rho(x)
    #rhotab  (optional) : array of tabulated density profile
   # ctab    (optional) : array of tabulated coordinate positions for density bins in table

  #The function rhofunc is assumed to be a real function with a single argument:

 # real function rhofunc(r)
   #real, intent(in) :: r

  # rhofunc = 1./r**2

 # end function rhofunc

 # If the ctab array is not present, the table rhotab(:) is assumed to
#  contain density values on equally spaced bins between cmin and cmax


"""
Function to perform strecth mapping and get density profile for IC.
Subroutine given above. Based on Daniel Price Fortran code.
The given density function needs to be set in:
get_mass(x,xmin) or
get_mass_rcyl(x,xmin)
get_mass_r
Was made more general before with D.getrho(ri), but for now it is done by writing functions
in this script, to allow for numba to be used and for faster calculations.
"""


def set_density_profile(D,xmin,xmax,rhofunc=None,rhotab=None,xtab=None,start=0,igeom=1,icoord=0):
        is_r    = is_rspherical(igeom,icoord)
        is_rcyl = is_rcylindrical(igeom,icoord)
        #get total mass integrated along coordinate direction to normalise the density profile
        #  (we assume total mass is 1 for both particles and the desired profile,
        #   i.e. that the same total mass is in both)
        if (rhotab != None):
           nt = len(rhotab)

           # if no coordinate table is passed, create a table
           # of equally spaced points between xmin and xmax

           if (xtab != None):
              if (len(xtab) < nt):
                  print('set_density_profile','coordinate table different size to density table')
                  exit
              xtable[1:nt] = xtab[1:nt]
           else:
              xtable[1:nt]=np.linspace(xmin,xmax,nt)
           # calculate mass in coord
           if (is_r):
              get_mass_tab_r(masstab,rhotab,xtable)
           elif (is_rcyl):
              get_mass_tab_rcyl(masstab,rhotab,xtable)
           else:
              get_mass_tab(masstab,rhotab,xtable)

           totmass    = yinterp(masstab,xtable[1:nt],xmax)
           rho_at_min = yinterp(rhotab,xtable[1:nt],xmin)

        else:

           if (is_r):
              totmass = get_mass_r(D,xmax,xmin)
           elif (is_rcyl):
              totmass = get_mass_rcyl(D,xmax,xmin)
           else:
              totmass = get_mass(xmax,xmin)

           rho_at_min = D.getrho(xmin)
           print("rhoatmin ",rho_at_min)

        if (is_r):
           rhozero = totmass/(4./3.*np.pi*(xmax**3 - xmin**3))
        elif (is_rcyl):
           rhozero = totmass/(np.pi*(xmax**2 - xmin**2))
        else:
            rhozero = totmass/(xmax-xmin)
        print("Total mass ",totmass," rhozero ",rhozero)

        x=np.array([0,0,0], dtype=np.float64)
        xt=np.array([0,0,0], dtype=np.float64)
        for i in range(start,len(D.x)):
           x[0] = D.x[i]
           x[1] = D.y[i]
           x[2] = D.z[i]
           hi = 1
           if (igeom > 1):
              xt[0],xt[1],xt[2] = coord_transform(x[0],x[1],x[2],igeom_cartesian,igeom)
              xold = xt[icoord]
           else:
              xt = x
              xold = x[icoord]
           if (is_r):
              fracmassold = 4./3.*np.pi*rhozero*(xold**3 - xmin**3)
           elif (is_rcyl):
              fracmassold = np.pi*rhozero*(xold**2 - xmin**2)
           else:
              fracmassold = rhozero*(xold - xmin)
           if (xold > xmin):  # if x=0 do nothing
              its   = 0
              xprev = 0.0
              xi    = xold   # starting guess
              # calc func to determine if tol is met
              # tablefunction still fortran
              if (rhotab != None):
                 func  = yinterp(masstab,xtable[1:nt],xi)
              else:
                 if (is_r):
                    func  = get_mass_r(D,xi,xmin)
                 elif (is_rcyl):
                    func  = get_mass_rcyl(D,xi,xmin)
                 else:
                    func  = get_mass(xi,xmin)

              func = func - fracmassold
              xminbisect = xmin
              xmaxbisect = xmax
              bisect = False  # use Newton-Raphson by default
              while (np.abs(func/totmass) > tol and its < maxits):
                 xprev = xi
                 its   = its + 1
                 # tablefunction still fortran
                 if (rhotab != None):
                    func  = yinterp(masstab,xtable[1:nt],xi) - fracmassold
                    if (is_r):
                       dfunc = 4.*np.pi*xi**2*yinterp(rhotab,xtable[1:nt],xi)
                    elif (is_rcyl):
                       dfunc = 2.*np.pi*xi*yinterp(rhotab,xtable[1:nt],xi)
                    else:
                       dfunc = yinterp(rhotab,xtable[1:nt],xi)
                 else:
                    if (is_r):
                       func  = get_mass_r(D,xi,xmin) - fracmassold
                       dfunc = 4.*np.pi*xi**2*D.getrho(xi)
                    elif (is_rcyl):
                       func  = get_mass_rcyl(D,xi,xmin) - fracmassold
                       dfunc = 2.*np.pi*xi*D.getrho(xi)
                    else:
                       func  = get_mass(xi,xmin) - fracmassold
                       dfunc = D.getrho(xi)
                 if (bisect):
                    if (func > 0.):
                       xmaxbisect = xi
                    else:
                       xminbisect = xi
                    xi = 0.5*(xminbisect + xmaxbisect)
                 else:
                    # Newton-Raphson
                    xi = xprev - func/dfunc
                    # do not allow N-R to take big jumps
                    if(xprev != 0.0):
                        if (abs(xi) < 0.9*abs(xprev)):
                           xi = 0.9*xprev
                        elif (abs(xi) > 1.1*abs(xprev)):
                           xi = 1.1*xprev
                    # Use bisection if Newton-Raphson is failing
                    if (xi > xmax or xi < xmin or its > maxits_nr):
                       bisect = True
                       xi = 0.5*(xminbisect + xmaxbisect)


# tablefunction still fortran
              if (rhotab!=None):
                 rhoi = yinterp(rhotab,xtable[1:nt],xi)
              else:
                 rhoi = D.getrho(xi)

              xt[icoord] = xi
              x[0],x[1],x[2] = coord_transform(xt[0],xt[1],xt[2],igeom,igeom_cartesian)
              D.x[i]= x[0]
              D.y[i] = x[1]
              D.z[i] = x[2]
              D.rho[i] = rhoi
              if (its >= maxits):
                  print('set_density_profile','Stretch mapping not converged')
                  print(i,its,maxits)



        print(' >>>>>> done')


#--------------------------------------------------------------
#+
#  query function for whether we have spherical r direction
#+
#--------------------------------------------------------------
def is_rspherical(igeom,icoord):
     if (igeom==igeom_spherical and icoord==0):
         return True
     else:
         return False


#--------------------------------------------------------------
#+
#  query function for whether we have cylindrical r direction
#+
#--------------------------------------------------------------
def is_rcylindrical(igeom,icoord):
     if (igeom==igeom_cylindrical and icoord==0):
         return True
     else:
         return False


#--------------------------------------------------------------
#+
#  Integrate to get total mass along the coordinate direction
#+
#--------------------------------------------------------------
#
# mass integrated along spherical radius

def get_mass_r(D,r,rmin):

     dr = (r - rmin)/ngrid
     dmprev     = 0.
     get_mass_r = 0.
     for i in range(1,ngrid):
        ri         = rmin + i*dr
        dmi        = ri*ri*D.getrho(ri)*dr
        get_mass_r = get_mass_r + 0.5*(dmi + dmprev) # trapezoidal rule
        dmprev     = dmi

     get_mass_r = 4.*np.pi*get_mass_r
     return get_mass_r
#
# mass integrated along cylindrical radius
#
def get_mass_rcyl(D,rcyl,rmin):
     dr = (rcyl - rmin)/ngrid
     dmprev   = 0.
     get_mass_rcyl = 0.
     for i in range(1,ngrid):
        ri            = rmin + i*dr
        dmi           = ri*D.getrho(ri)*dr
        get_mass_rcyl = get_mass_rcyl + 0.5*(dmi + dmprev) # trapezoidal rule
        dmprev        = dmi

     get_mass_rcyl = 2.*np.pi*get_mass_rcyl
     return get_mass_rcyl
#
# mass integrated along cartesian direction
#
# For faster iterations function rhoxi defined here by hand(maybe should create table of function instead or find better way).
@njit
def get_mass(x,xmin):
     dx = (x - xmin)/ngrid
     dmprev   = 0.
     get_mass = 0.
     delta=0.0025
     for i in range(1,ngrid):
         xi       = xmin + i*dx
         #rampf = 1./(1. + np.exp(-2.*(xi)/delta))

         #fac1 = (1. - 1./(1. + np.exp(2.*(xi+0.25)/delta)))
         #fac2 = (1. - 1./(1. + np.exp(2.*(0.25-xi)/delta)))
         #rampf = fac1*fac2
         #rhoxi= 1. + rampf*(2. - 1.)

         #rhoxi = 1.0*np.exp(-xi**2/(2.0))
         #rhoxi= 1.0
         #rhoxi=D.getrho(xi)
         rhoxi=1.0 + 0.00047535759060996836*(4*4.435843640411451*np.cos(6.283185307179586*(xi+0.5))-np.sin(6.283185307179586*(xi+0.5)))
         #rhoxi = 1.0 + 10**(-3)*np.sin(6.283185307179586*(xi + 0.5))
         #rhoxi = np.exp(-xi**2/(2.0))
         #dmi      = D.getrho(xi)*dx
         dmi = rhoxi*dx
         get_mass = get_mass + 0.5*(dmi + dmprev) # trapezoidal rule
         dmprev   = dmi
     return get_mass


#------------------------------------
#+
#  Same as above, but fills a table
#+
#------------------------------------
#
# version that integrates along spherical radius

def get_mass_tab_r(masstab,rhotab,rtab):


     masstab[0] = 0.
     dmprev     = 0.
     rprev      = rtab[0]
     for i in range(1,len(rhotab)):
        ri     = rtab(i)
        dr     = ri - rprev
        dmi    = ri*ri*rhotab(i)*dr
        masstab[i] = masstab[i-1] + 0.5*(dmi + dmprev) # trapezoidal rule
        dmprev = dmi
        rprev  = ri

     masstab = 4.*np.pi*masstab


# version that integrates along cylindrical radius

def get_mass_tab_rcyl(masstab,rhotab,rtab):


     masstab[0] = 0.
     dmprev     = 0.
     rprev      = rtab[0]
     for i in range(1,len(rhotab)):
        ri     = rtab[i]
        dr     = ri - rprev
        dmi    = ri*rhotab(i)*dr
        masstab[i] = masstab[i-1] + 0.5*(dmi + dmprev) # trapezoidal rule
        dmprev = dmi
        rprev  = ri

     masstab = 2.*np.pi*masstab


# version that integrates along a cartesian direction

def get_mass_tab(masstab,rhotab,xtab):

     masstab[0] = 0.
     dmprev     = 0.
     xprev      = xtab[0]
     for i in range(1,len(rhotab)):
        xi     = xtab[i]
        dx     = xi - xprev
        dmi    = rhotab[i]*dx
        masstab[i] = masstab[i-1] + 0.5*(dmi + dmprev) # trapezoidal rule
        dmprev = dmi
        xprev  = xi

def coord_transform(xin,yin,zin,itypein,itypeout):
    #input output same
    if(itypein==itypeout):
        xout = xin
        yout = yin
        zout = zin
    #input cyl output cart r,phi,z -> x,y,z
    if(itypein==2 and itypeout==1):
        xout = xin*np.cos(yin)
        yout = xin*np.sin(yin)
        zout = zin
    #input cart output cyl
    if (itypein==1 and itypeout==2):
        xout = np.sqrt(xin*xin+yin*yin)
        yout = np.arctan2(yin,xin) # phi
        zout = zin # z
    return xout,yout,zout
