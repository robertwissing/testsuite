#!/usr/bin/env python3.10
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 16:35:25 2023

@author: robertwi
"""

import h5py as h5
import numpy as np
import readtipsy as tip
import matplotlib.pylab as plt
f = h5.File("ic_imridisc.hdf5", "r")
#f = h5.File("disc_single.hdf5", "r")

#header = f['Header']
parttype0 = f['PartType0']
parttype5 = f['PartType5']


coord=parttype0['Coordinates']
vel=parttype0['Velocities']

#GAS PARTICLES
x=coord[:,0]
y=coord[:,1]
z=coord[:,2]
vx=vel[:,0]
vy=vel[:,1]
vz=vel[:,2]
m=np.array(parttype0['Masses'])
rho=np.array(parttype0['Density'])
u=np.array(parttype0['InternalEnergy'])
h=np.array(parttype0['SmoothingLength'])
pot=np.array(parttype0['Potential'])
cs=np.array(parttype0['SoundSpeed'])
P=np.array(parttype0['Pressure'])

r=np.sqrt(x*x+y*y+z*z)

q=0.5
cs0=cs*r**q
#plt.plot(r,cs*r,'.')
#plt.plot(r,cs**2*rho,r,P,'.')
plt.plot(r,cs0,'.')
adasD=asd2
massgas=np.sum(m)
#m=m*(1e-8/massgas)
#massgas2=np.sum(m)


#SINK PARTICLES
coordS=parttype5['Coordinates']
velS=parttype5['Velocities']
xS=coordS[:,0]
yS=coordS[:,1]
zS=coordS[:,2]
vxS=velS[:,0]
vyS=velS[:,1]
vzS=velS[:,2]
mS=np.array(parttype5['Masses'])
Srad=np.array(parttype5['SinkRadius'])
BHAL=np.array(parttype5['BH_AccretionLength'])
potS=np.array(parttype5['Potential'])

npartgas=len(x)
npartsink=len(xS)
npart=npartgas+npartsink
tgdata=np.zeros((npartgas,12))
tsdata=np.zeros((npartsink,11))
xh=h+np.min(h)-h
xu=u-u+max(u)
tgdata[:,0]=m
tgdata[:,1]=x
tgdata[:,2]=y
tgdata[:,3]=z
tgdata[:,4]=vx
tgdata[:,5]=vy
tgdata[:,6]=vz
tgdata[:,7]=rho
tgdata[:,8]=xu
#tgdata[:,9]=xh
tgdata[:,9]=0.001
tgdata[:,11]=pot

tsdata[:,0]=mS
tsdata[:,1]=xS
tsdata[:,2]=yS
tsdata[:,3]=zS
tsdata[:,4]=vxS
tsdata[:,5]=vyS
tsdata[:,6]=vzS
tsdata[:,9]=0.001
tsdata[:,10]=potS


massbh=np.sum(mS)



data_header=np.array([npart, 3, npartgas, 0, npartsink,0])

file_string='accretiondiskstd'

tip.writetipsy(tgdata,[],tsdata,file_string,data_header,time=np.array([0.0]))


print("GAS PARTICLES: ")
datasetNames = [n for n in parttype0.keys()]
for n in datasetNames:
    print(n)
    
    
print("\n SINK PARTICLES: ")
datasetNames = [n for n in parttype5.keys()]
for n in datasetNames:
    print(n)

print("HEADER: ")
# Get and print list of datasets within the H5 file
datasetNames = [n for n in f.keys()]
for n in datasetNames:
    print(n)
    
    
    
