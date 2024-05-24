#!/usr/bin/env python3
# convert makedisk files to gasoline binary
 ##with magneticfield
from setup import setup
import numpy as np
import smoothlength as smth
from writeascii import writeascii, convertfromcodeunits
import readtipsy as tip


D = setup()

file_gas = 'gas.dat'
file_halo = 'halo.dat'
file_disk = 'disk.dat'
file_bulge = 'bulge.dat'
    
fG = open(file_gas, 'r');
Dgas=np.loadtxt(fG)
#Dgas=Dgas[1:5,:]
D.ngas=len(Dgas)
fG.close()

fH = open(file_halo, 'r');
Ddark=np.loadtxt(fH)
#Ddark=Ddark[1:5,:]
D.ndark=len(Ddark)
fH.close()

fB = open(file_bulge, 'r');
fD = open(file_disk,'r');
Dbulge = np.loadtxt(fB)
Ddisk = np.loadtxt(fD)
#Ddisk=Ddisk[1:5,:]
#Dbulge=Dbulge[1:5,:]
Dstar = np.vstack((Dbulge,Ddisk))
D.nstar = len(Dstar)
fB.close()
fD.close()
D.npart=D.ngas+D.ndark+D.nstar
D.time = 0.0
temperature = 10**4
D.ndim = 3
x=np.hstack((Dgas[:,0],Ddark[:,0],Dstar[:,0]))
y=np.hstack((Dgas[:,1],Ddark[:,1],Dstar[:,1]))
z=np.hstack((Dgas[:,2],Ddark[:,2],Dstar[:,2]))
vx=np.hstack((Dgas[:,3],Ddark[:,3],Dstar[:,3]))
vy=np.hstack((Dgas[:,4],Ddark[:,4],Dstar[:,4]))
vz=np.hstack((Dgas[:,5],Ddark[:,5],Dstar[:,5]))
mass=np.hstack((Dgas[:,6],Ddark[:,6],Dstar[:,6]))
u=np.hstack([temperature]*D.ngas)
rho=np.hstack([1.0]*D.ngas)
h=smth.getsmooth2(Dgas[:,6],[1.0]*D.ngas,200)

mgas=Dgas[1,6]
mgastot=sum(Dgas[:,6])
#cmtopc=3.24077929*10**-19
#mh=8.4152384*10**-58
#rhohthres=(10**-3/cmtopc**3)*mh
#jeans=2*(64*mgas/((4*np.pi/3)*rhohthres))**(1./3.)
#softparam=jeans/4
softparam=80*10**-3 #kpc
softdark=np.hstack([softparam]*D.ndark)
softstar=np.hstack([softparam]*D.nstar)
metalsgas=np.hstack((([0.02041]*D.ngas)))
metalsstar=np.hstack((([0.02041]*D.nstar)))
tform=[10**10]*D.nstar
potgas=[0.0]*D.ngas
potdark=[0.0]*D.ndark
potstar=[0.0]*D.nstar
kmtokpc=3.24077929 * 10**-17
vx=np.multiply(vx,kmtokpc)
vy=np.multiply(vy,kmtokpc)
vz=np.multiply(vz,kmtokpc)
convertfromcodeunits(D,1.0,10**9, 0,1,1)


tipsygasdata = np.transpose(np.vstack((Dgas[:,6],Dgas[:,0],Dgas[:,1],Dgas[:,2],Dgas[:,3],Dgas[:,4],Dgas[:,5],rho,u,h,metalsgas,potgas)))
tipsydarkdata = np.transpose(np.vstack((Ddark[:,6],Ddark[:,0],Ddark[:,1],Ddark[:,2],Ddark[:,3],Ddark[:,4],Ddark[:,5],softdark,potdark)))
tipsystardata = np.transpose(np.vstack((Dstar[:,6],Dstar[:,0],Dstar[:,1],Dstar[:,2],Dstar[:,3],Dstar[:,4],Dstar[:,5],metalsstar,tform,softstar,potstar)))
data_header = np.array([D.npart, D.ndim, D.ngas, D.ndark, D.nstar, 0])
file_string = "agoralowres.00000"
tip.writetipsy(tipsygasdata,tipsydarkdata,tipsystardata,file_string,data_header)

Bx=np.hstack(([0.0]*D.ngas))
By=np.hstack(([0.0]*D.ngas))
Bz=np.hstack(([10**-12]*D.ngas))
B=np.transpose(np.array([Bx,By,Bz]))
tip.writealltipsyaux(B,file_string)