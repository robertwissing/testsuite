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


file="/Users/robertwi/OneDrive - Lund University/PHD/MHD/initial/LOW_512XLARGEMHDVERTDENSWend64FBSB64AB05.00850_splash.grid"
griddata,[nx,ny,nz],time=tip.readgrid(file)
asd=griddata[0,:,:,:]

for i in range(len(griddata[:,0,0,0])):
#for i in range(1):
    test=str(i)
    newfile='LOW_512XLARGEMHDVERTDENSHIGH' + str(i) + '.splash_col' + '.raw'
    #if(i==0):
    #    griddatanew = np.log10(griddata[i,:,:,:])
    #else:
    griddatanew = griddata[i,:,:,:]
    griddatanew2 = np.transpose(griddatanew)
    tip.writevtkgrid(griddatanew2,newfile)

