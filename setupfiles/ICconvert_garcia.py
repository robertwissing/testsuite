#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 17:32:16 2020

@author: robertwi
"""

import sys
sys.path.insert(0, '$MY_ANALYSIS_PATH')
import numpy as np
from matplotlib import pyplot as plt
import random as random
import readtipsy as tip
import os
import fnmatch
from shutil import copyfile
from collections import OrderedDict

"""
Script to convert columntxt files from Garcia Senz to SPLASH for the Axisymmetric MHD paper.
"""

file="./KH_AxisSPHYNX_t2.8.dat"
file_string="KHaxist28"
file_stringcyl="cloudaximu10"
#Cloudmu_inftyt1_1e12
#Cloudmu10t1_1e12
#Cloudmu20t1_1e12
test=np.loadtxt(file)
nres=12
z=np.linspace(-1,1,nres)
timy=0.0

#columns are:   r,z,massparticle,h,rho,pres,vr, vz,Br,Bz
npart=len(test[:,0])
X=np.zeros((npart*nres,12))
#for kh x,y,mass,h,density
for i in range(nres):
    zarr=np.array(npart*[z[i]])
    densarr=np.array(npart*[1.0])
    X[npart*i:npart*(i+1),0]=test[:,2]
    X[npart*i:npart*(i+1),1]=test[:,0]
    X[npart*i:npart*(i+1),2]=test[:,1]
    X[npart*i:npart*(i+1),3]=zarr
    X[npart*i:npart*(i+1),7]=test[:,4]
    X[npart*i:npart*(i+1),8]=test[:,4]
    X[npart*i:npart*(i+1),9]=test[:,3]
    X[npart*i:npart*(i+1),4]=test[:,4]
    X[npart*i:npart*(i+1),5]=test[:,4]
    X[npart*i:npart*(i+1),6]=test[:,4]
    X[npart*i:npart*(i+1),10]=test[:,4]
    X[npart*i:npart*(i+1),11]=test[:,4]

#X[:,8]=X[:,7]

data_header=np.array([npart*nres, 3, npart*nres, 0, 0, 0])
ascale=timy
time=np.array([ascale])
tip.writetipsy(X,[],[],file_string,data_header,time)

#
#r1=test[:,0]
#phi1=np.linspace(0,2*np.pi,nres)
#z1=test[:,1]
#
#X=np.zeros((npart*(nres-1),12))
#for i in range(nres-1):
#    phiarr=np.array(npart*[phi1[i]])
#    X[npart*i:npart*(i+1),0]=test[:,2]
#    X[npart*i:npart*(i+1),1]=test[:,0]*np.cos(phi1[i])
#    X[npart*i:npart*(i+1),2]=test[:,0]*np.sin(phi1[i])
#    X[npart*i:npart*(i+1),3]=test[:,1]
#    X[npart*i:npart*(i+1),7]=test[:,4]
#    X[npart*i:npart*(i+1),9]=test[:,3]
#
#X[:,4]=X[:,0]
#X[:,5]=X[:,0]
#X[:,6]=X[:,0]
#X[:,8]=X[:,0]
#X[:,10]=X[:,0]
#X[:,11]=X[:,0]
#
#data_header=np.array([npart*(nres-1), 3, npart*(nres-1), 0, 0, 0])
#ascale=timy
#time=np.array([ascale])
#tip.writetipsy(X,[],[],file_stringcyl,data_header,time)
