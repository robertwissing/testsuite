#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 22:27:53 2019

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
tgdata,tddata,tsdata,data_header,time,label=tip.readalltipsy('mhdcollapse200mu10gd.00760');
N=data_header[0]
ngas=data_header[2]
ndark=data_header[3]
nstar=data_header[4]

rad=tgdata[:,1]*tgdata[:,1]+tgdata[:,2]*tgdata[:,2]+tgdata[:,3]*tgdata[:,3]
dens=tgdata[:,7]
mapp=((rad<0.000064) & (dens<50000.0))

tgdatanew=tgdata[mapp,:]
#tgdatanew=tgdata[:,:]
#fieldarr=np.transpose(np.array([ tgdatanew[:,7],tgdatanew[:,7],tgdatanew[:,4],tgdatanew[:,5],tgdatanew[:,6],tgdatanew[:,12],tgdatanew[:,13],tgdatanew[:,14] ]))
#fieldarr=np.transpose(np.array([ tgdatanew[:,7],tgdatanew[:,12],tgdatanew[:,13],tgdatanew[:,14] ]))
fieldarr=np.transpose(np.array([ tgdatanew[:,7],tgdatanew[:,7]]  ))
N=len(tgdatanew)
data_header=np.array([data_header[1],N])
tip.writeprepASCII(tgdatanew[:,1:4],'qhullParticles.txt',data_header)
tip.writeprepASCII(fieldarr,'qhullFields.txt',data_header)
