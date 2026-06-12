#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 07:28:57 2022

@author: robertwi
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 22:06:33 2019

@author: robertwi
"""

import pynbody
import matplotlib.pylab as plt
import matplotlib
import pylab
import pynbody; from pynbody.analysis import profile;
import pynbody.plot as pp
import pynbody.plot.sph as sph
from pynbody.derived import array
import numpy as np
import matplotlib.gridspec as gridspec
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams.update({'font.size': 18})

@pynbody.derived_array
def BField(sim) :
    BFieldi=np.empty_like(sim['vel'])
    CodetoGauss=0.00019134614504210336*10**6
    BFieldi[:,0]=sim['BFieldx']*CodetoGauss
    BFieldi[:,1]=sim['BFieldy']*CodetoGauss
    BFieldi[:,2]=sim['BFieldz']*CodetoGauss
    BFieldi.sim = sim
    BFieldi.units = 'Mpc'
    return BFieldi

@pynbody.derived_array
def BFieldAbs(sim) :
    CodetoGauss=0.00019134614504210336*10**6
    BFieldAbsi=np.sqrt(sim['BFieldx']*sim['BFieldx']+sim['BFieldy']*sim['BFieldy']+sim['BFieldz']*sim['BFieldz'])*CodetoGauss
    BFieldAbsi.sim = sim
    BFieldAbsi.units = 'Mpc^2'
    return BFieldAbsi

#@pynbody.derived_array
#def betatherm(sim) :
#    BFieldAbsi=np.sqrt(sim['BFieldx']*sim['BFieldx']+sim['BFieldy']*sim['BFieldy']+sim['BFieldz']*sim['BFieldz'])
#    BFieldAbsi.sim = sim
#    return BFieldAbsi
#

def titles(title,fig,axs,T1,T2,T3,T4,T5):
    axs[0].set_title(T1)
    axs[1].set_title(T2)
    axs[2].set_title(T3)
    axs[3].set_title(T4)
    axs[4].set_title(T5)
    
def titles2(title,fig,axs,T1,T2,T3,T4):
    axs[0].set_title(T1)
    axs[1].set_title(T2)
    axs[2].set_title(T3)
    axs[3].set_title(T4)
    
def titles3(title,fig,axs,T1,T2,T3):
    axs[0].set_title(T1)
    axs[1].set_title(T2)
    axs[2].set_title(T3)
    
def titles4(title,fig,axs,T1,T2):
    axs[0].set_title(T1)
    axs[1].set_title(T2)
    
def fieldsetup(s):
    pg = pynbody.analysis.profile.Profile(s.g,rmax=20,nbins=1000,type='lin')
    
    rxy=np.array([0.0])
    rxy=np.append(rxy,pg['rxy'])
    rxy=np.array(rxy)
    ttt3=np.copy(s.g['vphi'])
    npx=len(rxy)
    for j in range(npx-1):
        mapp = (s.g['rxy']>rxy[j]) & (s.g['rxy']<rxy[j+1]) 
        ttt3[mapp]=ttt3[mapp]-pg['vphi'][j]
    s.g['vphi4']=ttt3
    s.g['vturb']=np.sqrt(s.g['vrxy']*s.g['vrxy']+s.g['vphi4']*s.g['vphi4']+s.g['vz']*s.g['vz'])
    s.g['uDotBcleandiss']=np.sqrt(s.g['DivB']*s.g['DivB'])
    s.g['Pvel'] = 0.5*s.g['rho']*(s.g['vturb']*s.g['vturb'])

    CodetoGauss=0.00019134614504210336*10**6
    s.g['vxtu'] = s.g['vx']-s.g['v_mean'][:,0]
    s.g['vytu'] = s.g['vy']-s.g['v_mean'][:,1]
    s.g['vztu'] = s.g['vz']-s.g['v_mean'][:,2]
    s.g['vtu'] = np.sqrt(s.g['vxtu']**2+s.g['vytu']**2+s.g['vztu']**2)
    s.g['Babs'] = np.sqrt(s.g['BFieldx']**2+s.g['BFieldy']**2+s.g['BFieldz']**2)
    s.g['Curlabs'] = np.sqrt(s.g['CurlBx']**2+s.g['CurlBy']**2+s.g['CurlBz']**2)
    s.g['Curlnorm'] = s.g['Curlabs']/s.g['Babs']
    s.g['Pmag'] = 0.5*s.g['Babs']*s.g['Babs']
    s.g['mfrac'] = s.g['MassHot']/s.g['mass']
    s.g['Pth'] = 0.6666666667*s.g['rho']*(s.g['uHot']*s.g['mfrac']+s.g['u']*(1.0-s.g['mfrac']))
    s.g['betath'] = s.g['Pth']/s.g['Pmag']
    s.g['betavel'] = s.g['Pvel']/s.g['Pmag']
    s.g['cs'] = np.sqrt(1.6666666667*s.g['Pth']/(s.g['rho']))
    s.g['lambdaJ'] = (s.g['cs']*np.sqrt(np.pi/s.g['rho']))/s.g['smooth']
    s.g['mach'] = np.sqrt(s.g['vturb']*s.g['vturb'])/s.g['cs']
    s.g['mach2'] = np.sqrt((s.g['vx']*s.g['vx']+s.g['vy']*s.g['vy']+s.g['vz']*s.g['vz']))/s.g['cs']
    s.g['divBerr'] = s.g['smoothlength']*np.sqrt(s.g['DivB']*s.g['DivB'])/s.g['Babs']
    s.g['Br']    =(s.g['pos'] * (s.g['BField']/CodetoGauss)).sum(axis=1) / s.g['r']
    s.g['B2']    =((s.g['BField']/CodetoGauss) ** 2).sum(axis=1)
    s.g['Bt']    =np.sqrt(s.g['B2'] - s.g['Br'] ** 2) 
    s.g['Bx']    = s.g['BFieldx']*CodetoGauss
    s.g['By']    = s.g['BFieldy']*CodetoGauss
    s.g['Bz']    = s.g['BFieldz']*CodetoGauss
    s#.g['Brxy']    =(s.g['pos'][:, 0:2] * s.g['BField'][:, 0:2]/CodetoGauss).sum(axis=1) / s.g['rxy']
    s.g['Brxy'] = ((s.g['x'] * s.g['BFieldx'] + s.g['y'] * s.g['BFieldy']) / s.g['rxy'])*CodetoGauss
    s.g['Bphi']    = (s.g['x'] * s.g['BFieldy'] - s.g['y'] * s.g['BFieldx']) / s.g['rxy']*CodetoGauss
    s.g['Bphi'][np.where(s.g['Bphi'] != s.g['Bphi'])] = 0
    s.g['MwRphi'] = s.g['Brxy']*s.g['Bphi']
    s.g['MwRz'] = s.g['Brxy']*s.g['Bz']
    s.g['MwZphi'] = s.g['Bz']*s.g['Bphi']  
    return s

def fieldsetup2(s):
    pg = pynbody.analysis.profile.Profile(s.g,rmax=20,nbins=1000,type='lin')
    
    rxy=np.array([0.0])
    rxy=np.append(rxy,pg['rxy'])
    rxy=np.array(rxy)
    ttt3=np.copy(s.g['vphi'])
    npx=len(rxy)
    for j in range(npx-1):
        mapp = (s.g['rxy']>rxy[j]) & (s.g['rxy']<rxy[j+1]) 
        ttt3[mapp]=ttt3[mapp]-pg['vphi'][j]
    s.g['vphi4']=ttt3
    s.g['vturb']=np.sqrt(s.g['vrxy']*s.g['vrxy']+s.g['vphi4']*s.g['vphi4']+s.g['vz']*s.g['vz'])
    s.g['uDotBcleandiss']=np.sqrt(s.g['DivB']*s.g['DivB'])
    s.g['Pvel'] = 0.5*s.g['rho']*(s.g['vturb']*s.g['vturb'])
 
    CodetoGauss=0.00019134614504210336*10**6
    s.g['vxtu'] = s.g['vx']-s.g['v_mean'][:,0]
    s.g['vytu'] = s.g['vy']-s.g['v_mean'][:,1]
    s.g['vztu'] = s.g['vz']-s.g['v_mean'][:,2]
    s.g['vtu'] = np.sqrt(s.g['vxtu']**2+s.g['vytu']**2+s.g['vztu']**2)
    s.g['Babs'] = np.sqrt(s.g['BFieldx']**2+s.g['BFieldy']**2+s.g['BFieldz']**2)
    s.g['Curlabs'] = np.sqrt(s.g['CurlBx']**2+s.g['CurlBy']**2+s.g['CurlBz']**2)
    s.g['Curlnorm'] = s.g['Curlabs']/s.g['Babs']
    s.g['Pmag'] = 0.5*s.g['Babs']*s.g['Babs']
    s.g['Pth'] = 0.6666666667*s.g['rho']*(s.g['temp']*4.86263e-06)
    s.g['betath'] = s.g['Pth']/s.g['Pmag']
    s.g['betavel'] = s.g['Pvel']/s.g['Pmag']
    s.g['cs'] = np.sqrt(1.6666666667*s.g['Pth']/(s.g['rho']))
    s.g['lambdaJ'] = (s.g['cs']*np.sqrt(np.pi/s.g['rho']))/s.g['smooth']
    s.g['mach'] = np.sqrt(s.g['vturb']*s.g['vturb'])/s.g['cs']
    s.g['mach2'] = np.sqrt((s.g['vx']*s.g['vx']+s.g['vy']*s.g['vy']+s.g['vz']*s.g['vz']))/s.g['cs']
    s.g['divBerr'] = s.g['smoothlength']*np.sqrt(s.g['DivB']*s.g['DivB'])/s.g['Babs']
    s.g['Br']    =(s.g['pos'] * (s.g['BField']/CodetoGauss)).sum(axis=1) / s.g['r']
    s.g['B2']    =((s.g['BField']/CodetoGauss) ** 2).sum(axis=1)
    s.g['Bt']    =np.sqrt(s.g['B2'] - s.g['Br'] ** 2)
    s.g['Bx']    = s.g['BFieldx']*CodetoGauss
    s.g['By']    = s.g['BFieldy']*CodetoGauss
    s.g['Bz']    = s.g['BFieldz']*CodetoGauss
    #s.g['Brxy']    =((s.g['pos'][:, 0:2] * s.g['BField'][:, 0:2]/CodetoGauss).sum(axis=1) / s.g['rxy'])*CodetoGauss
    s.g['Brxy'] = ((s.g['x'] * s.g['BFieldx'] + s.g['y'] * s.g['BFieldy']) / s.g['rxy'])*CodetoGauss
    s.g['Bphi']    = ((s.g['x'] * s.g['BFieldy'] - s.g['y'] * s.g['BFieldx']) / s.g['rxy'])*CodetoGauss
    s.g['Bphi'][np.where(s.g['Bphi'] != s.g['Bphi'])] = 0
    s.g['MwRphi'] = s.g['Brxy']*s.g['Bphi']
    s.g['MwRz'] = s.g['Brxy']*s.g['Bz']
    s.g['MwZphi'] = s.g['Bz']*s.g['Bphi']  
    
    return s





#s1 = pynbody.load('./LOW_8XLARGEMHDVERTDENSWend64FBBWweakNEWTIMESTEP.02000')
#s1 = pynbody.load('./LOW_8XLARGEMHDVERTDENSWend64FBSB64weak.01980')
#s1 = pynbody.load('./LOW_8XLARGEMHDVERTDENSWend64FBBWweak.01000')
#s1 = pynbody.load('./LOW_8XLARGEMHDVERTDENSWend64FBSB64weak.01000')
#s2 = pynbody.load('./LOW_8XLARGEMHDVERTDENSWend64FBSB64J14weak.01000')
#s3 = pynbody.load('./LOW_8XLARGEMHDVERTDENSWend64FBSB64J56weak.01000')
#s4 = pynbody.load('./LOW_8XLARGEMHDVERTDENSWend64FBSB64J174weak.01000')

s1 = pynbody.load('./LOW_64XLARGEMHDVERTDENSWend64FBSB64BW.01820')
#s1 = pynbody.load('./FEEDBACKSTRENGTH/LOW_64XLARGEMHDVERTDENSWend64FBSB64BW.01820')
#s1 = pynbody.load('./FEEDBACKSTRENGTH/LOW_64XLARGEMHDVERTDENSWend64FBSB64BW.00250')
s1b = pynbody.load('./LOW_64XLARGEMHDVERTDENSWend64FBSB1.01600')
s2 = pynbody.load('./LOW_64XLARGEMHDVERTDENSWend64FB05SB64.01820')
s3 = pynbody.load('./LOW_64XLARGEMHDVERTDENSWend64FBSB64.01750')
s4 = pynbody.load('./LOW_64XLARGEMHDVERTDENSWend64FB2SB64.01820')
#
#pp.sfh(s1b)
#pp.sbprofile(s1b)
pp.sfh(s1,trange=[0,2],nbins=500,clear=True)
pp.sfh(s1b,trange=[0,2],nbins=500)
pp.sfh(s2,trange=[0,2],nbins=500)
pp.sfh(s3,trange=[0,2],nbins=500)
pp.sfh(s4,trange=[0,2],nbins=500)
plt.legend(['BW','FB1NB1','FB05NB64','FB1NB64','FB2NB64'])

#pp.sfh(s2,linestyle='dashed',color='k')
asd=Asd2


#s1 = pynbody.load('./JEANSFILESFB/LOW_8XLARGEMHDVERTDENSWend64FBSB64weak.03000')
#s2 = pynbody.load('./FEEDBACKSTRENGTH/LOW_8XLARGEMHDVERTDENSWend64FBBWweakNEWTIMESTEP.03000')
#s3 = pynbody.load('./JEANSFILESFB/LOW_8XLARGEMHDVERTDENSWend64FBSB64J56weak.02000')
#s4 = pynbody.load('./JEANSFILESFB/LOW_8XLARGEMHDVERTDENSWend64FBSB64J174weak.02000')

#s1 = pynbody.load('./FEEDBACKSTRENGTH/LOW_8XLARGEMHDVERTDENSWend64FBSB64weak.03000')
#s2 = pynbody.load('./FEEDBACKSTRENGTH/LOW_8XLARGEMHDVERTDENSWend64FBBWFB2weak.02000')
#s3 = pynbody.load('./FEEDBACKSTRENGTH/LOW_8XLARGEMHDVERTDENSWend64FBSB64FB05weak.03000')
#s4 = pynbody.load('./FEEDBACKSTRENGTH/LOW_8XLARGEMHDVERTDENSWend64FBSB64FB2weak.03000')

#s1 = pynbody.load('./FEEDBACKSTRENGTH/LOW_8XLARGEMHDVERTDENSWend64FBSB64weak.03000')
#s2 = pynbody.load('./FEEDBACKSTRENGTH/LOW_8XLARGEMHDVERTDENSWend64FBBWFB2weak.02000')
#s3 = pynbody.load('./FEEDBACKSTRENGTH/LOW_8XLARGEMHDVERTDENSWend64FBSB64FB05weak.03000')
#s4 = pynbody.load('./FEEDBACKSTRENGTH/LOW_8XLARGEMHDVERTDENSWend64FBSB64FB2weak.03000')

#s1 = pynbody.load('./FIELDSTRENGTH/LOW_8XLARGEMHDVERTDENSWend64FBSB64beta01.00005')
#s2 = pynbody.load('./FIELDSTRENGTH/LOW_8XLARGEMHDVERTDENSWend64FBSB64beta01.00040')
#s3 = pynbody.load('./FIELDSTRENGTH/LOW_8XLARGEMHDVERTDENSWend64FBSB64beta01.00400')
#s1 = pynbody.load('./DIFFSTRENGTH/LOW_8XLARGEMHDVERTDENSWend64FBSB64weakAR025.02000')
#s2 = pynbody.load('./DIFFSTRENGTH/LOW_64XLARGEMHDVERTDENSWend64FBSB64AB025.02000')
#s4 = pynbody.load('./DIFFSTRENGTH/LOW_8XLARGEMHDVERTDENSWend64FBSB64weakAR025.02000')
#s4 = pynbody.load('./DIFFSTRENGTH/LOW_64XLARGEMHDVERTDENSWend64FBSB64AB025.02000')

#s2 = pynbody.load('./FEEDBACKSTRENGTH/LOW_8XLARGEMHDVERTDENSWend64FBSB1weak.03000')
#s1 = pynbody.load('./FEEDBACKSTRENGTH/LOW_8XLARGEMHDVERTDENSWend64FBBWweakNEWTIMESTEP.03000')
#s3 = pynbody.load('./FEEDBACKSTRENGTH/LOW_8XLARGEMHDVERTDENSWend64FBSB200weak.03000')
#s4 = pynbody.load('./FIELDSTRENGTH/LOW_8XLARGEMHDVERTDENSWend64FBSB64weakAR025.02000')
    


#s5 = pynbody.load('./LOW_8XLARGEMHDVERTDENSWend64FBSB64FB05weak.01000') #1
#s6 = pynbody.load('./LOW_8XLARGEMHDVERTDENSWend64FBSB64FB2weak.01000') #200
#s5 = pynbody.load('./LOW_8XLARGEMHDVERTDENSWend64FBSB1weak.01000')
#s6 = pynbody.load('./LOW_8XLARGEMHDVERTDENSWend64FBSB200weak.01000')

#s1 = pynbody.load('./LOW_8XLARGEMHDVERTDENSWend64FBBWweak.01000')
#s2 = pynbody.load('./LOW_8XLARGEMHDVERTDENSWend64FBBW64J14weak.01000')
#s3 = pynbody.load('./LOW_8XLARGEMHDVERTDENSWend64FBBW64J56weak.01000')
#s4 = pynbody.load('./LOW_8XLARGEMHDVERTDENSWend64FBBW64J174weak.01000')
#s1 = pynbody.load('./LOW_8XLARGEMHDVERTDENSWend64FBSB64weak.01000')



#s5 = pynbody.load('./LOW_64XFBBWLARGEMHDVERTWend64SBFB200.01400')
#s5 = pynbody.load('./LOW_8XLARGEMHDVERTDENSWend64FBBW64J174weak.01000')
#s5 = pynbody.load('./LOW_64XLARGEMHDVERTDENSWend64FBSB64J56.00990')
#s5 = pynbody.load('JEANSSB/LOW_8XLARGEMHDVERTDENSWend64FBSB64J56weak.01000')
#s5 = pynbody.load('FBruns/LOW_8XLARGEMHDVERTWend64FBSB200.02000')                


if(s1['mass'].units==s2['mass'].units==s3['mass'].units):
    print('OK GOOD TO GO')
    print(s1['mass'].units)
else:
    print('NOPP')
    print(s1['mass'].units)
    print(s2['mass'].units)
    print(s3['mass'].units)
    print(s4['mass'].units)
    print(s1b['mass'].units)
    asd=asd22222
    
T1='density'
T2='$|B|$'
T3='$B_{r}$'
T5='$B_z$'
T4='$B_{\phi}$'
time=''

pynbody.analysis.angmom.faceon(s1)
#pynbody.analysis.angmom.faceon(s1b)
pynbody.analysis.angmom.faceon(s2)
pynbody.analysis.angmom.faceon(s3)
#pynbody.analysis.angmom.faceon(s4)

s1=fieldsetup2(s1)
#s1b=fieldsetup(s1b)
s2=fieldsetup(s2)
s3=fieldsetup(s3)
#s4=fieldsetup(s4)






cmapx="RdBu"
cmapx2="inferno"
cmapx3="jet"
cmap3x="jet"

s1.g["BField"]



lintest=[-1e-3,1e-3]
linme=1e-5

rabBvminx=1e-12*10**6
rabBvmaxx=1e-4*10**6

densrabBvminx=5e-7
densrabBvmaxx=2e-1

phirabBvminx=-1e-4*10**6
phirabBvmaxx=1e-4*10**6

zabBvminx=1e-12
zabBvmaxx=1e-5

mxz=35


f2, axs2 = plt.subplots(1,5,figsize=(20,5),sharex=True,sharey=True,squeeze=True)
#daamg2=sph.velocity_image(s1.g, width=mxz, mode='stream',vector_color='white',av_z="rho", qty="rho", cmap = cmapx2,denoise=False,subplot=axs2[0],show_cbar=True,vmin=densrabBvminx,vmax=densrabBvmaxx, clear=False,ret_im=True)
daamg=sph.image(s2.g, width=mxz, qty="BFieldAbs",av_z="rho", cmap = cmapx,denoise=False,subplot=axs2[1],show_cbar=True,vmin=rabBvminx,vmax=rabBvmaxx, clear=False,ret_im=True)
daamg3=sph.image(s2.g, width=mxz, qty="Brxy",av_z="rho", cmap = cmapx3,denoise=False,subplot=axs2[2],show_cbar=True,vmin=phirabBvminx,vmax=phirabBvmaxx, clear=False,ret_im=True,linthresh=linme)
daamg4=sph.image(s2.g, width=mxz, qty="Bphi",av_z="rho", cmap = cmapx3,denoise=False,subplot=axs2[3],show_cbar=True,vmin=phirabBvminx,vmax=phirabBvmaxx, clear=False,ret_im=True,linthresh=linme)
daamg5=sph.image(s2.g, width=mxz, qty="Bz",av_z="rho", cmap = cmapx3,denoise=False,subplot=axs2[4],show_cbar=True,vmin=phirabBvminx,vmax=phirabBvmaxx, clear=False,ret_im=True,linthresh=linme)
daamg2=sph.image(s2.g, width=mxz, qty="rho",units="g cm^-2", cmap = cmapx2,denoise=False,subplot=axs2[0],show_cbar=True,vmin=densrabBvminx,vmax=densrabBvmaxx, clear=False,ret_im=True)
titles('BField'+time,f2,axs2,T1,T2,T3,T4,T5)
cbar_ax = f2.add_axes([0.30,0.18,0.1,0.02])
cbar = f2.colorbar(daamg, cax=cbar_ax,orientation='horizontal',ticks=[1e-6,1e-3, 1e-0])
cbar_ax = f2.add_axes([0.14,0.18,0.1,0.02])
cbar = f2.colorbar(daamg2, cax=cbar_ax,orientation='horizontal',ticks=[1e-6, 1e-4,1e-2])
cbar_ax = f2.add_axes([0.46,0.18,0.1,0.02])
cbar = f2.colorbar(daamg3, cax=cbar_ax,orientation='horizontal',ticks=[-1e-0, 0,1e-0])
cbar_ax = f2.add_axes([0.62,0.18,0.1,0.02])
cbar = f2.colorbar(daamg4, cax=cbar_ax,orientation='horizontal',ticks=[-1e-0,0,1e-0])
cbar_ax = f2.add_axes([0.78,0.18,0.1,0.02])
cbar = f2.colorbar(daamg5, cax=cbar_ax,orientation='horizontal',ticks=[-1e-0,0,1e-0])

f2.savefig('bfieldfront.png')

pynbody.analysis.angmom.sideon(s1)

f2, axs2 = plt.subplots(1,5,figsize=(20,5),sharex=True,sharey=True,squeeze=True)
#daamg2=sph.velocity_image(s1.g, width=mxz, mode='stream',vector_color='white',av_z="rho", qty="rho", cmap = cmapx2,denoise=False,subplot=axs2[0],show_cbar=True,vmin=densrabBvminx,vmax=densrabBvmaxx, clear=False,ret_im=True)
daamg=sph.velocity_image(s1.g, width=mxz,key_x=0.15, key_y=0.45, mode='stream',vector_resolution=120,density=1.0,vector_qty='BField',vector_color='black', qty="BFieldAbs",av_z="rho", cmap = cmapx,denoise=False,subplot=axs2[1],show_cbar=True,vmin=rabBvminx,vmax=rabBvmaxx, clear=False,ret_im=True)

#daamg=sph.image(s1.g, width=mxz, qty="BFieldAbs",av_z="rho", cmap = cmapx,denoise=False,subplot=axs2[1],show_cbar=True,vmin=rabBvminx,vmax=rabBvmaxx, clear=False,ret_im=True)
daamg3=sph.image(s1.g, width=mxz, qty="Bx",av_z="rho", cmap = cmapx3,denoise=False,subplot=axs2[2],show_cbar=True,vmin=phirabBvminx,vmax=phirabBvmaxx, clear=False,ret_im=True,linthresh=linme)
daamg4=sph.image(s1.g, width=mxz, qty="By",av_z="rho", cmap = cmapx3,denoise=False,subplot=axs2[3],show_cbar=True,vmin=phirabBvminx,vmax=phirabBvmaxx, clear=False,ret_im=True,linthresh=linme)
daamg5=sph.image(s1.g, width=mxz, qty="Bz",av_z="rho", cmap = cmapx3,denoise=False,subplot=axs2[4],show_cbar=True,vmin=phirabBvminx,vmax=phirabBvmaxx, clear=False,ret_im=True,linthresh=linme)
daamg2=sph.image(s1.g, width=mxz, qty="rho",units="g cm^-2", cmap = cmapx2,denoise=False,subplot=axs2[0],show_cbar=True,vmin=densrabBvminx,vmax=densrabBvmaxx, clear=False,ret_im=True)
titles('BField'+time,f2,axs2,T1,T2,T3,T4,T5)
cbar_ax = f2.add_axes([0.30,0.18,0.1,0.02])
cbar = f2.colorbar(daamg, cax=cbar_ax,orientation='horizontal',ticks=[1e-6,1e-3, 1e-0])
cbar_ax = f2.add_axes([0.14,0.18,0.1,0.02])
cbar = f2.colorbar(daamg2, cax=cbar_ax,orientation='horizontal',ticks=[1e-6, 1e-4,1e-2])
cbar_ax = f2.add_axes([0.46,0.18,0.1,0.02])
cbar = f2.colorbar(daamg3, cax=cbar_ax,orientation='horizontal',ticks=[-1e-0, 0,1e-0])
cbar_ax = f2.add_axes([0.62,0.18,0.1,0.02])
cbar = f2.colorbar(daamg4, cax=cbar_ax,orientation='horizontal',ticks=[-1e-0,0,1e-0])
cbar_ax = f2.add_axes([0.78,0.18,0.1,0.02])
cbar = f2.colorbar(daamg5, cax=cbar_ax,orientation='horizontal',ticks=[-1e-0,0,1e-0])


#T1='BW $0.5\epsilon_{SN}$'
#T2='BW $2\epsilon_{SN}$'
#T3='SB $0.5\epsilon_{SN}$'
#T4='SB $2\epsilon_{SN}$'
#T1=r'$t=0$'
#T2=r'$t=40Myr$'
#T3=r'$t=400Myr$'
#T4=r'$\alpha_B=0.25$ $t=2Gyr$'
#T1=r'$N_{FB}=1$'
#T2=r'$N_{FB}=64$'
#T3=r'$N_{FB}=200$'
#T4=r'$\alpha_B=0.25$ $t=2Gyr$'

#T1=r'Low $\alpha_B=0.25$'
#T2=r'Medium $\alpha_B=0.25$'
#T2=r'$N_{FB}=64$'

T1='BW $1\epsilon_{SN}$'
T2='SB $1\epsilon_{SN} \\ N_{FB}=1$'
T3='SB $1\epsilon_{SN} \\ N_{FB}=200$'
T4='SB $1\epsilon_{SN} \\ N_{FB}=64$'
T5='SB $2\epsilon_{SN} \\ N_{FB}=64$'

pynbody.analysis.angmom.sideon(s1)
pynbody.analysis.angmom.sideon(s2)
pynbody.analysis.angmom.sideon(s3)
#pynbody.analysis.angmom.sideon(s4)
f2, axs2 = plt.subplots(2,3,figsize=(15,10),sharex=True,sharey=True,squeeze=True)
#daamg2=sph.velocity_image(s1.g, width=mxz, mode='stream',vector_color='white',av_z="rho", qty="rho", cmap = cmapx2,denoise=False,subplot=axs2[0],show_cbar=True,vmin=densrabBvminx,vmax=densrabBvmaxx, clear=False,ret_im=True)
sph.image(s1.g, width=mxz, qty="BFieldAbs",av_z="rho", cmap = cmapx,denoise=False,subplot=axs2[1,0],show_cbar=True,vmin=rabBvminx,vmax=rabBvmaxx, clear=False,ret_im=True)
sph.image(s2.g, width=mxz, qty="BFieldAbs",av_z="rho", cmap = cmapx,denoise=False,subplot=axs2[1,1],show_cbar=True,vmin=rabBvminx,vmax=rabBvmaxx, clear=False,ret_im=True)
daamg2=sph.image(s3.g, width=mxz, qty="BFieldAbs",av_z="rho", cmap = cmapx,denoise=False,subplot=axs2[1,2],show_cbar=True,vmin=rabBvminx,vmax=rabBvmaxx, clear=False,ret_im=True)
#daamg2=sph.image(s4.g, width=mxz, qty="BFieldAbs",av_z="rho", cmap = cmapx,denoise=False,subplot=axs2[4],show_cbar=True,vmin=rabBvminx,vmax=rabBvmaxx, clear=False,ret_im=True)
pynbody.analysis.angmom.faceon(s1)
pynbody.analysis.angmom.faceon(s2)
pynbody.analysis.angmom.faceon(s3)
sph.image(s1.g, width=mxz, qty="BFieldAbs",av_z="rho", cmap = cmapx,denoise=False,subplot=axs2[0,0],show_cbar=True,vmin=rabBvminx,vmax=rabBvmaxx, clear=False,ret_im=True)
sph.image(s2.g, width=mxz, qty="BFieldAbs",av_z="rho", cmap = cmapx,denoise=False,subplot=axs2[0,1],show_cbar=True,vmin=rabBvminx,vmax=rabBvmaxx, clear=False,ret_im=True)
daamg2=sph.image(s3.g, width=mxz, qty="BFieldAbs",av_z="rho", cmap = cmapx,denoise=False,subplot=axs2[0,2],show_cbar=True,vmin=rabBvminx,vmax=rabBvmaxx, clear=False,ret_im=True)

titles3('BField'+time,f2,axs2,T1,T2,T3)
cbar_ax = f2.add_axes([0.91,0.25,0.02,0.50])
cbar = f2.colorbar(daamg2, cax=cbar_ax)


f2, axs2 = plt.subplots(1,3,figsize=(15,5),sharex=True,sharey=True,squeeze=True)
#daamg2=sph.velocity_image(s1.g, width=mxz, mode='stream',vector_color='white',av_z="rho", qty="rho", cmap = cmapx2,denoise=False,subplot=axs2[0],show_cbar=True,vmin=densrabBvminx,vmax=densrabBvmaxx, clear=False,ret_im=True)
sph.image(s1.g, width=mxz, qty="rho",units="g cm^-2", cmap = cmapx2,denoise=False,subplot=axs2[0],show_cbar=True,vmin=densrabBvminx,vmax=densrabBvmaxx, clear=False,ret_im=True)
#sph.image(s1b.g, width=mxz, qty="rho",units="g cm^-2", cmap = cmapx2,denoise=False,subplot=axs2[1],show_cbar=True,vmin=densrabBvminx,vmax=densrabBvmaxx, clear=False,ret_im=True)
sph.image(s2.g, width=mxz, qty="rho",units="g cm^-2", cmap = cmapx2,denoise=False,subplot=axs2[1],show_cbar=True,vmin=densrabBvminx,vmax=densrabBvmaxx, clear=False,ret_im=True)
daamg3=sph.image(s3.g, width=mxz, qty="rho",units="g cm^-2", cmap = cmapx2,denoise=False,subplot=axs2[2],show_cbar=True,vmin=densrabBvminx,vmax=densrabBvmaxx, clear=False,ret_im=True)
#sph.image(s4.g, width=mxz, qty="rho",units="g cm^-2", cmap = cmapx2,denoise=False,subplot=axs2[3],show_cbar=True,vmin=densrabBvminx,vmax=densrabBvmaxx, clear=False,ret_im=True)

titles3('BField'+time,f2,axs2,T1,T2,T3)
cbar_ax = f2.add_axes([0.91,0.25,0.02,0.50])
cbar = f2.colorbar(daamg3, cax=cbar_ax)


f2, axs2 = plt.subplots(1,4,figsize=(20,5),sharex=True,sharey=True,squeeze=True)
#daamg2=sph.velocity_image(s1.g, width=mxz, mode='stream',vector_color='white',av_z="rho", qty="rho", cmap = cmapx2,denoise=False,subplot=axs2[0],show_cbar=True,vmin=densrabBvminx,vmax=densrabBvmaxx, clear=False,ret_im=True)
sph.image(s1.g, width=mxz, qty="Bz",av_z="rho", cmap = cmapx3,denoise=False,subplot=axs2[0],show_cbar=True,vmin=phirabBvminx,vmax=phirabBvmaxx, clear=False,ret_im=True)
#sph.image(s1b.g, width=mxz, qty="Bz",av_z="rho", cmap = cmapx3,denoise=False,subplot=axs2[1],show_cbar=True,vmin=phirabBvminx,vmax=phirabBvmaxx, clear=False,ret_im=True)
sph.image(s2.g, width=mxz, qty="Bz",av_z="rho", cmap = cmapx3,denoise=False,subplot=axs2[2],show_cbar=True,vmin=phirabBvminx,vmax=phirabBvmaxx, clear=False,ret_im=True)
sph.image(s3.g, width=mxz, qty="Bz",av_z="rho", cmap = cmapx3,denoise=False,subplot=axs2[3],show_cbar=True,vmin=phirabBvminx,vmax=phirabBvmaxx, clear=False,ret_im=True)
daamg4=sph.image(s4.g, width=mxz, qty="Bz",av_z="rho", cmap = cmapx3,denoise=False,subplot=axs2[4],show_cbar=True,vmin=phirabBvminx,vmax=phirabBvmaxx, clear=False,ret_im=True)

titles('BField'+time,f2,axs2,T1,T2,T3,T4)
#titles2('BField'+time,f2,axs2,T1,T2,T3,T4)
cbar_ax = f2.add_axes([0.91,0.25,0.02,0.50])
cbar = f2.colorbar(daamg4, cax=cbar_ax,ticks=[-1e-0,0,1e-0])

