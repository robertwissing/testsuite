#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 15:24:38 2018

@author: robertwi
"""
import numpy as np

def getsmooth(D,N):
    m=np.asarray(D.mass)
    rho=np.asarray(D.rho)
    hfact3=N*3/(2**3*4*np.pi)
    D.h=(hfact3*m/rho)**(1./3.)
    return D.h

def getsmooth2(Dmass,Drho,N):
    m=np.asarray(Dmass)
    rho=np.asarray(Drho)
    hfact3=N*3/(2**3*4*np.pi)
    h=(hfact3*m/rho)**(1./3.)
    return h

#def getdensity(x,y,z,m):
 #   npart=len(x)
  #  r=np.sqrt(x**2+y**2+z**2)
   # for i in range (1,npart):
    #    for j in range (i,npart):
     #       rij=abs(r(i)-r(j))
      #      if 2*h(i) < rij:
       #         pair_i(niac+1)=i
                #pair_j(niac)=
        #        neigh(i)=neigh(i)+1
                #neigh(j)=neight(j)+1
         #   if 2*h(j) < rij:
          #      pair_j(niac+1)=
           #     neigh(j)=neight(j)+1
            #niac=niac+1