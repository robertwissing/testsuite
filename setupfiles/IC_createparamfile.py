#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 11:37:47 2020

@author: robertwi
"""

"""
Script that generates a parameter file for given IC.
"""


def createparamfile(filename,dxP,dyP,dzP,periodic,delta,nsteps,dmsol,dkpc,out,adi,molw,gamma,grav,cosmo,rhoit,ns,glasssetup,dICdensRsmooth,dICdensprofile,dICdensdir,dICdensR,dICdensinner,dICdensouter):
    f = open('datafiles/' + filename+'.param', 'w')
    f.write('achInFile \t = '+ filename + '.00000' +' \n')
    f.write('achOutName \t = '+ filename + ' \n')
    if(periodic==1):
        f.write('bPeriodic \t = 1 \n')
        f.write('dxPeriod \t = '+ str(dxP) + ' \n')
        f.write('dyPeriod \t = '+ str(dyP) + ' \n')
        f.write('dzPeriod \t = '+ str(dzP) + ' \n')
    else:
        f.write('bPeriodic \t = 0 \n')
    f.write('dDelta \t = '+ str(delta) + ' \n')
    f.write('nSteps \t = '+ str(nsteps) + ' \n')
    f.write('dMsolUnit \t = '+ str(dmsol) + ' \n')
    f.write('dKpcUnit \t = '+ str(dkpc) + ' \n')
    f.write('iCheckInterval \t = '+ str(round(out*10)) + ' \n')
    f.write('iLogInterval \t = 1 \n')
    f.write('iOutInterval \t = '+ str(out) + ' \n')
    f.write('dMeanMolWeight \t = '+ str(molw) + ' \n')
    if(adi==1):
        f.write('bGasAdiabatic \t = 1 \n')
        f.write('dConstGamma \t = '+ str(gamma) + ' \n')
    else:
        f.write('bGasIsothermal \t = 1 \n')
        f.write('dConstGamma \t = 1.0000000001  \n')
    if(grav==1):
        f.write('bDoGravity \t = 1 \n')
        f.write('dhMinOverSoft \t = 0.1 \n')
    else:
        f.write('bDoGravity \t = 0 \n')

    f.write('#bLogTiming \t = 1 \n')
    f.write('bStandard \t = 1 \n')
    f.write('bVDetails \t = 2 \n')
    if(rhoit==1):
        f.write('nSmooth \t = ' + str(ns+20) + ' \n')
        f.write('nSmoothMean \t = ' + str(ns) + ' \n')
    else:
        f.write('nSmooth \t = ' + str(ns) + ' \n')
    f.write('#dConstAlpha \t = 2 \n')
    f.write('#dConstBeta \t = 4 \n')
    f.write('dExtraStore \t = 20.0 \n')
    f.write('#dThermalDiffusionCoeff \t = 0.03 \n')
    f.write('dEta \t = 0.175 \n')
    f.write('dEtaCourant \t = 0.4 \n')
    f.write('dEtaDiffusion \t = 0.25 \n')
    f.write('iMaxRung \t = 29 \n')
    f.write('iBinaryOutput \t = 1 \n')
    f.write('dTheta \t = 0.6 \n')
    f.write('bViscosityLimiter \t = 0 \n')
    f.write('bViscosityLimitdt \t = 0 \n')
    if(cosmo==1):
        f.write('dTheta2 \t = 0.8 \n')
        f.write('dHubble0 \t = 2.894405 \n')
        f.write('dOmega0 \t = 1.0 \n')
        f.write('dLambda \t = 0.0 \n')
        f.write('dRedTo \t = 0.0 \n')
        f.write('bComove \t = 1 \n')
        f.write('dFracNoDomainDecomp \t = 0.1 \n')
        f.write('nTruncateRung \t = 128 \n')
    else:
        f.write('bEwald \t = 0 \n')
    if(glasssetup==1):
        f.write('dICdensRsmooth \t = ' + str(dICdensRsmooth) + ' \n')
        f.write('dICdensprofile \t = ' + str(dICdensprofile) + ' \n')
        f.write('dICdensdir \t = ' + str(dICdensdir) + ' \n')
        f.write('dICdensR \t = ' + str(dICdensR) + ' \n')
        f.write('dICdensinner \t = ' + str(dICdensinner) + ' \n')
        f.write('dICdensouter \t = ' + str(dICdensouter) + ' \n')
    f.close()
