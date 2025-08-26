#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 22:27:53 2019

@author: robertwi
"""
import numpy as np
from matplotlib import pyplot as plt
import random as random
import readtipsy as tip
import os
import fnmatch
import sys
from shutil import copyfile
import pynbody
from SPH_general import *
from SPH_render import *

def match_columns_with_zeros(X_base_cols, X_other):
    """
    Ensures X_other has the same number of columns as X_base_cols.
    If X_other has fewer columns, append zero-filled columns.
    Returns the updated X_other.
    """
    target_cols = X_base_cols
    other_cols = X_other.shape[1]
    if other_cols < target_cols:
        # Need to add (target_cols - other_cols) columns of zeros
        diff = target_cols - other_cols
        zero_cols = np.zeros((X_other.shape[0], diff))
        X_other = np.hstack((X_other, zero_cols))
    # If it already matches or has more columns (unexpected?), do whatever fits your logic.
    return X_other

def setQuantity(entry,output,choice=0):
    tgdata,tddata,tsdata,data_header,time=tip.readtipsy(entry);
    X, labels = tip.readalltipsyaux(entry)
    N=data_header[0]
    ndim=data_header[1]
    ngas=data_header[2]
    ndark=data_header[3]
    nstar=data_header[4]
    print(N,ndim,ngas,ndark,nstar,time)
    if(X.size!=0):
        current_rows, ncols = X.shape
        if current_rows < N:
                    # We need more rows to match ndark
                    missing_rows = N - current_rows
                    zero_rows = np.zeros((missing_rows, ncols))
                    # Append these rows to X
                    X = np.concatenate((X, zero_rows), axis=0)
        elif current_rows > N:
                    # If X has too many rows, decide whether to keep the first ndark or handle it differently
                    X = X[:N, :]

    print(X.shape)
    Xgas=X[:ngas]
    Xdark=X[ngas:ngas+ndark]
    Xstar=X[ngas+ndark:N]
    



    print("What quantity do you want to modify? \n")
    print("1: Velocity \n")
    print("2: Softening \n")
    print("3: Tform \n")  
    print("4: Temperature \n")
    print("5: Mass \n")
    print("6: Metals \n")
    print("7: Center particles based on potential \n")
    print("8: Remove particles \n")
    print("9: Change direction \n")
    print("10: Time \n")
    print("11: Boverrho \n")
    print("12: Magnetic field removal \n")
    print("13: Impact (clone) \n")
    print("14: Spin \n")
    print("15: Add particle \n")
    print("16: Output BField vector file \n")
    print("17: Smooth BField \n")
    print("18: Merge particles \n")
    print("19: Split particles \n")
    choice = input("Option: ")
    choice = int(choice)
    print(choice)
    if (choice==1):
        coorsh=input("1: constant vel 2: shear vel 3: kh vel")
        coorsh=int(coorsh)
        if(coorsh==1):
            direc=int(input("1:x 2:y 3:z 4: r"))
            velval=float(input("Velocity value (applied to all particles): "))
            if(direc < 4):
                tgdata[:,3+direc]=velval
                tddata[:,3+direc]=velval
                tsdata[:,3+direc]=velval
            else:
                 mappr1= tgdata[:,1] < 0.0 
                 mappr2= tgdata[:,1] >= 0.0
                 mappr3= tgdata[:,2] < 0.0
                 mappr4= tgdata[:,2] >= 0.0
                 velval2=velval*(tgdata[:,1]*tgdata[:,1]+tgdata[:,2]*tgdata[:,2])
                 tgdata[mappr1,4]=-velval2[mappr1]
                 tgdata[mappr2,4]=velval2[mappr2]
                 tgdata[mappr3,5]=-velval2[mappr3]
                 tgdata[mappr4,5]=velval2[mappr4]
        elif(coorsh==2):
            velval=float(input("Orbital velocity, (q=3/2): "))
            tgdata[:,5]=-1.5*velval*tgdata[:,1]
            tddata[:,5]=-1.5*velval*tddata[:,1]
            tsdata[:,5]=-1.5*velval*tsdata[:,1]
            setzerovel=int(input("1: Set zero vel in other dir?: "))
            if setzerovel==1:
                tgdata[:,4]=0.0*tgdata[:,1]
                tddata[:,4]=0.0*tddata[:,1]
                tsdata[:,4]=0.0*tsdata[:,1]
                tgdata[:,6]=0.0*tgdata[:,1]
                tddata[:,6]=0.0*tddata[:,1]
                tsdata[:,6]=0.0*tsdata[:,1]
        elif(coorsh==3):
            delta=0.05
            fac1 = (1. - 1./(1. + np.exp(2.*(tgdata[:,1]+0.25)/delta)))
            fac2 = (1. - 1./(1. + np.exp(2.*(0.25-tgdata[:,1])/delta)))
            rampf = fac1*fac2
            velval=float(input("Velocity diff: "))
            print(len(fac1),len(rampf))
            tgdata[:,5]=1.0+rampf*velval
            #tddata[:,5]=1.0+rampf*velval
            #tsdata[:,5]=1.0+rampf*velval
            setzerovel=int(input("1: Set zero vel in other dir?: "))
            if setzerovel==1:
                tgdata[:,4]=0.0*tgdata[:,1]
                tddata[:,4]=0.0*tddata[:,1]
                tsdata[:,4]=0.0*tsdata[:,1]
                tgdata[:,6]=0.0*tgdata[:,1]
                tddata[:,6]=0.0*tddata[:,1]
                tsdata[:,6]=0.0*tsdata[:,1]                   
        elif(coorsh==4):
            tgdata[:,4]=-tgdata[:,2]
            tgdata[:,5]=tgdata[:,1]
    if (choice==2):
        softval=input("Softening value: ")
        setzerovel=int(input("1: Set stars+gas 1 or darkmatter 2 or gas 3 or specific star 4 "))
        if setzerovel==1:
            tsdata[:,9]=softval
            tgdata[:,9]=softval
        if setzerovel==2:
            tddata[:,7]=softval
        if setzerovel==3:
            tgdata[:,9]=softval
        if setzerovel==4:
            order=int(input("order: "))
            tsdata[order,9]=softval
        
    if (choice==3):
        tformval=input("Tform value: ")
        mapp = tsdata[:,8] == 0.5
        tsdata[mapp,8]=tformval
        
    if (choice==4):
        tempval=input("Temperature value: ")
        tgdata[:,8]=tempval

    if (choice==5):
        massval=input("Mass value: ")
        tsdata[:,0]=massval
        tgdata[:,0]=massval
        tddata[:,0]=massval
        
    if (choice==6):
        metalval=input("Metal value: ")
        tsdata[:,7]=metalval
        tgdata[:,10]=metalval
    
    if (choice==7):
        minpotg=min(tgdata[:,11])
        minpotidx=np.argmin(tgdata[:,11])
        rmin=tgdata[minpotidx,1:4]
        print(rmin)
        tddata[:,1:4]=tddata[:,1:4]-rmin
        tsdata[:,1:4]=tsdata[:,1:4]-rmin
        tgdata[:,1:4]=tgdata[:,1:4]-rmin
    if (choice==8):
         print("Which particle type should be removed? \n")
         remtype=input("1: gas 2: DM 3: star  \n")
         remtype=int(remtype)
         print("Which quantity should removal be based on? \n")
         qtype=input("1: r 2: y 3: z 4: rcyl 5: nr of particle  \n")
         qtype=int(qtype)
         if(qtype == 5):
            val=input("value: ")
            val=int(val)
         else:
            dirval=input("1: below value or 2: above value: \n")
            dirval=int(dirval)
            val=input("value: ")
            val=float(val)
         
        # Load auxdata
         if (X.size != 0):
             Xg=X[0:ngas,:]
             Xd=X[ngas:ngas+ndark, :]
             Xs=X[ngas+ndark:ngas+ndark+nstar, :]
         if (remtype==1):
             tdata=np.copy(tgdata)
             mapp = np.ones(ngas,dtype=int)
         elif(remtype==2):
             tdata=np.copy(tddata)
             mapp = np.ones(ndark,dtype=int)
         else:
             tdata=np.copy(tsdata)
             mapp = np.ones(nstar,dtype=int)
         if(qtype == 1):
             if(dirval == 2):
                 mapp=np.sqrt(tdata[:,1]*tdata[:,1]+tdata[:,2]*tdata[:,2]+tdata[:,3]*tdata[:,3]) < val
             else:
                 mapp=np.sqrt(tdata[:,1]*tdata[:,1]+tdata[:,2]*tdata[:,2]+tdata[:,3]*tdata[:,3]) > val
         if(qtype == 2):
             if(dirval == 2):
                 mapp=np.sqrt(tdata[:,1]*tdata[:,1]) < val
             else:
                  mapp=np.sqrt(tdata[:,1]*tdata[:,1]) > val
         if(qtype == 3):
             if(dirval == 2):
                 mapp=np.sqrt(tdata[:,3]*tdata[:,3]) < val
             else:
                  mapp=np.sqrt(tdata[:,3]*tdata[:,3]) > val
         if(qtype == 4):
             if(dirval == 2):
                 mapp=np.sqrt(tdata[:,1]*tdata[:,1]+tdata[:,2]*tdata[:,2]) < val
             else:
                 mapp=np.sqrt(tdata[:,1]*tdata[:,1]+tdata[:,2]*tdata[:,2]) > val
         nnew=np.sum(mapp)
         if (remtype==1):
             if(qtype == 5):
                 nnew=nnew-1
                 tgdata=np.delete(tdata, val, axis=0)
             else:
                 tgdata=tdata[mapp,:]
                 Xg=Xg[mapp,:]
             ndiff=ngas-nnew
             data_header[2]=nnew
         elif(remtype==2):
             if(qtype == 5):
                 nnew=nnew-1
                 tddata=np.delete(tdata, val, axis=0)
             else:
                tddata=tdata[mapp,:]
                Xd=Xd[mapp,:]
             ndiff=ndark-nnew
             data_header[3]=nnew
         else:
             if(qtype == 5):
                 nnew=nnew-1
                 tsdata=np.delete(tdata, val, axis=0)
             else:
                tsdata=tdata[mapp,:]
                Xs=Xs[mapp,:]
             ndiff=nstar-nnew
             data_header[4]=nnew
         if (X.size != 0):
            X = np.concatenate((Xg, Xd, Xs))
         data_header[0]=(ngas+ndark+nstar)-ndiff
         print("Value:",data_header[0],data_header[2],data_header[3],data_header[4],nnew,ndiff)
    if (choice==9):
        coorsh=input("1: Replace this direction")
        coorsh=int(coorsh)
        coorsh2=input("2: With this direction")
        coorsh2=int(coorsh2)
        tgdatatemp=np.copy(tgdata)
        tddatatemp=np.copy(tddata)
        tsdatatemp=np.copy(tsdata)
        tgdata[:,coorsh]=tgdatatemp[:,coorsh2]
        tddata[:,coorsh]=tddatatemp[:,coorsh2]
        tsdata[:,coorsh]=tsdatatemp[:,coorsh2]
        tgdata[:,coorsh2]=tgdatatemp[:,coorsh]
        tddata[:,coorsh2]=tddatatemp[:,coorsh]
        tsdata[:,coorsh2]=tsdatatemp[:,coorsh]

    if (choice==10):
        newtime = input("Set time: ")
        time=np.array([float(newtime)])

    if (choice)==11:
        Bx,lab=tip.readtipsyaux(entry,'BFieldx')
        By,lab2=tip.readtipsyaux(entry,'BFieldy')
        Bz,lab3=tip.readtipsyaux(entry,'BFieldz')
        X2 = np.hstack([tgdata[:,7],([1.0]*ndark),([1.0]*nstar)])
        print(np.shape(X2))
        Bx2=Bx/X2
        By2=By/X2
        Bz2=Bz/X2
        print(np.shape(Bx2))
        X3=np.transpose(np.array([Bx2,By2,Bz2]))
        print(np.shape(X3))
    if (choice)==12:
        minpotg=min(tgdata[:,11])
        minpotidx=np.argmin(tgdata[:,11])
        rmin=tgdata[minpotidx,1:4]
        tddata[:,1:4]=tddata[:,1:4]-rmin
        tsdata[:,1:4]=tsdata[:,1:4]-rmin
        tgdata[:,1:4]=tgdata[:,1:4]-rmin
        
        tdata=tgdata
        dirval=input("1: below value or 2: above value: \n")
        dirval=int(dirval)
        val=input("value: ")
        val=float(val)
        X2 = np.hstack([tgdata[:,3],tddata[:,3],tsdata[:,3]])
        if(dirval == 2):
                       mapp=np.sqrt(X2*X2) < val
        else:
                       mapp=np.sqrt(X2*X2) > val
        Bx,lab=tip.readtipsyaux(entry,'BFieldx')
        By,lab2=tip.readtipsyaux(entry,'BFieldy')
        Bz,lab3=tip.readtipsyaux(entry,'BFieldz')
        Bx[mapp]=0.0
        By[mapp]=0.0
        Bz[mapp]=0.0
        X3=np.transpose(np.array([Bx,By,Bz]))
    if (choice==13):
        # Load data
        f = pynbody.load(entry)
        f.properties['z'] = 1
        
        # Align the object face-on and calculate virial radius and spin parameter
        pynbody.analysis.angmom.faceon(f)
        Vr = pynbody.analysis.halo.virial_radius(f, overden=178)
        print('Vr: ', Vr, 'at z: ', f.properties['z'])
        spin = pynbody.analysis.angmom.spin_parameter(f)
        print('Spin: ', spin)
        
        # Filter based on the virial radius
        coorsh = float(input("1: Object size or impact distance/2 (Rvir=50): "))
        R=Vr
        mapp = tddata[:, 1]**2 + tddata[:, 2]**2 + tddata[:, 3]**2 < R**2
        mapp2 = tsdata[:, 1]**2 + tsdata[:, 2]**2 + tsdata[:, 3]**2 < R**2
        mapp3 = tgdata[:, 1]**2 + tgdata[:, 2]**2 + tgdata[:, 3]**2 < R**2
        M = (np.sum(tddata[mapp, 0]) + np.sum(tsdata[mapp2, 0]) + np.sum(tgdata[mapp3, 0]))
        dens = M / (4 / 3 * np.pi * R**3)
        
        # Critical density calculations
        densc = pynbody.analysis.cosmology.rho_crit(f, unit="g cm**-3")
        densc2 = pynbody.analysis.cosmology.rho_crit(f, z=3, unit="g cm**-3")
        print("mvir", M, dens * 1.57446e-26, 'rhoc', densc, 'rhoc2', densc2, dens * 6.77331e-23 / 178, R, M)
        
        #coorsh4=M

        # Impact setup - gather user inputs
        #coorsh = float(input("1: Object size or impact distance/2 (Rvir=50): "))
        coorsh2 = int(input("2: Impact geometry (along plane or above plane (y or z)): "))
        coorsh3 = float(input("3: Impact Parameter (sin angle): "))
        coorsh4 = float(input("4: Object mass, Virial mass: "))
        coorsh5 = float(input("5: Tilt Impactor around x-axis (angle in degrees): "))
        coorsh6 = float(input("6: Tilt Impactor around y or z-axis (angle in degrees): "))

        
        # --- Check for Magnetic Fields ---
        # Find the indices for BFieldx, BFieldy, BFieldz if they exist in labels
        BField_indices = {}
        for i, label in enumerate(labels):
            if 'BFieldx' in label:
                BField_indices['x'] = i
            elif 'BFieldy' in label:
                BField_indices['y'] = i
            elif 'BFieldz' in label:
                BField_indices['z'] = i
        # Copy data for impactor and target
        Ximp = np.copy(X)
        tgdataimp = np.copy(tgdata)
        tddataimp = np.copy(tddata)
        tsdataimp = np.copy(tsdata)
        Xtar = np.copy(X)
        tgdatatar = np.copy(tgdata)
        tddatatar = np.copy(tddata)
        tsdatatar = np.copy(tsdata)

        # --- Apply rotations before adjusting the impactor's position ---
        # 1. Tilt around the x-axis using coorsh5
        angle_x_rad = np.radians(coorsh5)
        cos_x = np.cos(angle_x_rad)
        sin_x = np.sin(angle_x_rad)
        
        tgdataimp[:, [2, 3]] = tgdataimp[:, [2, 3]] @ np.array([[cos_x, -sin_x], [sin_x, cos_x]])
        tddataimp[:, [2, 3]] = tddataimp[:, [2, 3]] @ np.array([[cos_x, -sin_x], [sin_x, cos_x]])
        tsdataimp[:, [2, 3]] = tsdataimp[:, [2, 3]] @ np.array([[cos_x, -sin_x], [sin_x, cos_x]])

        tgdataimp[:, [5, 6]] = tgdataimp[:, [5, 6]] @ np.array([[cos_x, -sin_x], [sin_x, cos_x]])
        tddataimp[:, [5, 6]] = tddataimp[:, [5, 6]] @ np.array([[cos_x, -sin_x], [sin_x, cos_x]])
        tsdataimp[:, [5, 6]] = tsdataimp[:, [5, 6]] @ np.array([[cos_x, -sin_x], [sin_x, cos_x]])
        
        # Rotate magnetic field components if they exist
        if BField_indices:
            Ximp[:, [BField_indices['y'], BField_indices['z']]] = Ximp[:, [BField_indices['y'], BField_indices['z']]] @ np.array([[cos_x, -sin_x], [sin_x, cos_x]])
        
        # 2. Tilt around the selected axis (y or z) using coorsh6
        angle_other_rad = np.radians(coorsh6)
        cos_other = np.cos(angle_other_rad)
        sin_other = np.sin(angle_other_rad)
        
        if coorsh2 == 2:  # Tilt around the y-axis
            tgdataimp[:, [1, 3]] = tgdataimp[:, [1, 3]] @ np.array([[cos_other, -sin_other], [sin_other, cos_other]])
            tddataimp[:, [1, 3]] = tddataimp[:, [1, 3]] @ np.array([[cos_other, -sin_other], [sin_other, cos_other]])
            tsdataimp[:, [1, 3]] = tsdataimp[:, [1, 3]] @ np.array([[cos_other, -sin_other], [sin_other, cos_other]])

            tgdataimp[:, [4, 6]] = tgdataimp[:, [4, 6]] @ np.array([[cos_other, -sin_other], [sin_other, cos_other]])
            tddataimp[:, [4, 6]] = tddataimp[:, [4, 6]] @ np.array([[cos_other, -sin_other], [sin_other, cos_other]])
            tsdataimp[:, [4, 6]] = tsdataimp[:, [4, 6]] @ np.array([[cos_other, -sin_other], [sin_other, cos_other]])

            # Rotate magnetic field components around y-axis
            if BField_indices:
                Ximp[:, [BField_indices['x'], BField_indices['z']]] = Ximp[:, [BField_indices['x'], BField_indices['z']]] @ np.array([[cos_other, -sin_other], [sin_other, cos_other]])


        elif coorsh2 == 3:  # Tilt around the z-axis
            tgdataimp[:, [1, 2]] = tgdataimp[:, [1, 2]] @ np.array([[cos_other, -sin_other], [sin_other, cos_other]])
            tddataimp[:, [1, 2]] = tddataimp[:, [1, 2]] @ np.array([[cos_other, -sin_other], [sin_other, cos_other]])
            tsdataimp[:, [1, 2]] = tsdataimp[:, [1, 2]] @ np.array([[cos_other, -sin_other], [sin_other, cos_other]])
            
            tgdataimp[:, [4, 5]] = tgdataimp[:, [4, 5]] @ np.array([[cos_other, -sin_other], [sin_other, cos_other]])
            tddataimp[:, [4, 5]] = tddataimp[:, [4, 5]] @ np.array([[cos_other, -sin_other], [sin_other, cos_other]])
            tsdataimp[:, [4, 5]] = tsdataimp[:, [4, 5]] @ np.array([[cos_other, -sin_other], [sin_other, cos_other]])
            # Rotate magnetic field components around z-axis
            if BField_indices:
                Ximp[:, [BField_indices['x'], BField_indices['y']]] = Ximp[:, [BField_indices['x'], BField_indices['y']]] @ np.array([[cos_other, -sin_other], [sin_other, cos_other]])

           
        # Adjust position for impactor and target based on inputs
        Rvir = coorsh
        tgdataimp[:, 1] += 2 * Rvir
        tddataimp[:, 1] += 2 * Rvir
        tsdataimp[:, 1] += 2 * Rvir
        tgdataimp[:, coorsh2] += (2 * Rvir) * coorsh3
        tddataimp[:, coorsh2] += (2 * Rvir) * coorsh3
        tsdataimp[:, coorsh2] += (2 * Rvir) * coorsh3
        
        # Virial mass and velocity
        Mvir = coorsh4
        G = 1.0  # Gravitational constant
        vesc = np.sqrt((2 * G * (Mvir + Mvir)) / (2 * Rvir))
        vinf = 0.0
        vimp = np.sqrt(vesc**2 + vinf**2)
        vx1 = (Mvir / (Mvir + Mvir)) * vimp
        print("VESC IS ",vesc," velocity is ", vx1)
        # Set velocities for impactor and target
        tgdataimp[:, 4] -= vx1
        tddataimp[:, 4] -= vx1
        tsdataimp[:, 4] -= vx1
        tgdatatar[:, 4] += vx1
        tddatatar[:, 4] += vx1
        tsdatatar[:, 4] += vx1
        
        # Concatenate impactor and target data
        Xg = np.concatenate((Xtar[0:ngas, :], Ximp[0:ngas, :]))
        Xd = np.concatenate((Xtar[ngas:ngas+ndark, :], Ximp[ngas:ngas+ndark, :]))
        Xs = np.concatenate((Xtar[ngas+ndark:ngas+ndark+nstar, :], Ximp[ngas+ndark:ngas+ndark+nstar, :]))
        X = np.concatenate((Xg, Xd, Xs))
        print(X.shape)
        tgdata = np.concatenate((tgdatatar, tgdataimp))
        tsdata = np.concatenate((tsdatatar, tsdataimp))
        tddata = np.concatenate((tddatatar, tddataimp))
        
        # Update the header information
        ngas = len(tgdata[:, 1])
        ndark = len(tddata[:, 1])
        nstar = len(tsdata[:, 1])
        data_header[0] = ngas + ndark + nstar
        data_header[2] = ngas
        data_header[3] = ndark
        data_header[4] = nstar
        

        
    if (choice==14):
        w=[[0.0]*ngas,[0.0]*ngas,[0.0]*ngas]
        r2=tgdata[:,1]**2+tgdata[:,2]**2+tgdata[:,3]**2
        w0 = float ( input("spin value: ") )
        direction = int ( input("What direction: 1=x 2=y 3=z ") ) - 1
        posneg = float ( input("positive or negative?: -1=neg 1=pos") )
        wval=w0*np.exp(-r2)
        w[direction] = posneg*wval
        wx=np.hstack((w[0],([0.0]*ndark),([0.0]*nstar)))
        wy=np.hstack((w[1],([0.0]*ndark),([0.0]*nstar)))
        wz=np.hstack((w[2],([0.0]*ndark),([0.0]*nstar)))
        wfinal=np.transpose(np.array([wx,wy,wz]))
        tip.writealltipsyaux(wfinal,output,case="spin")
        
    if (choice==15):
        acc = int ( input("LISA primordial + secondary? = 1 : ") )
        if (acc==1):
            hr = float ( input("h/r ?  : ") )
            mass = 1.0
            soft = 0.5
            sink = -1.0
            print(tsdata)
            tsdata=np.array([mass,0.0,0.0,0.0,0.0,0.0,0.0,0.0,sink,soft,0.0])
            soft = hr*0.6
            mass = mass*0.0001
            tsdata=np.array([tsdata,np.array([mass,-1.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,soft,0.0])])
            data_header[4]=data_header[4]+2;
            data_header[0]=data_header[0]+2;
        else:
            typ = int ( input("What kind of particle do you want to add? 1: gas 2: darkmatter 3: star : ") )
            pos = np.fromstring(  input("Position enter like this: x,y,z default is 0,0,0"), dtype=float, sep=',' )
            vel = np.fromstring(  input("Velocity enter like this: vx,vy,vz or default is 0,0,0"), dtype=float, sep=',' )
            mass = float ( input("Mass of particle") )
            soft = float ( input("Softening of particle") )
            data_header[typ+1]=data_header[typ+1]+1;
            data_header[0]=data_header[0]+1;
            ## mass 1.0 soft 0.5
            ## 0.6*0.1 = 0.06 in the alignment run and 0.6*0.03 = 0.018
            tsdata=np.append(tsdata,np.array([mass,pos[0],pos[1],pos[2],vel[0],vel[1],vel[2],0.0,0.0,soft,0.0]), axis=1)

    if (choice==16):
        Bx,lab=tip.readtipsyaux(entry,'BFieldx')
        By,lab2=tip.readtipsyaux(entry,'BFieldy')
        Bz,lab3=tip.readtipsyaux(entry,'BFieldz')
        X2=np.transpose(np.array([Bx,By,Bz]))
        tip.writealltipsyaux(X2,output,case="mhd")
    if (choice==18):
        coorsh = int(input("1: What type of particle do you want to split/merge? : 1: gas, 2: dark, 3:star \n "))
        # Generate initial array
        # Set values based on input particle arrays:
        # gas: m,x,y,z,vx,vy,vz,rho,soft,T,metals,potential
        # dark: m,x,y,z,vx,vy,vz,soft,potential
        # star: m,x,y,z,vx,vy,vz,metals,tform,soft,potential
        # particles array structure for merging and splitting:
        # m,x,y,z,vx,vy,vz,rho,soft,T,metals,potential,Bx,By,Bz,wx,wy,wz,momenti,tform
        
        if(coorsh==1):
            tdata=np.copy(tgdata)
            particles = np.zeros((ngas, 21))
            particles[:,0:12]=tdata[:,0:12]
            print("test",len(tdata[:,0]),tdata.size/20)
            Xres=Xgas
        elif(coorsh==2):
            tdata=np.copy(tddata)
            particles = np.zeros((ndark, 21))
            particles[:,0:7]=tdata[:,0:7]
            particles[:,9]=tdata[:,7]
            particles[:,11]=tdata[:,8]
            Xres=Xdark
        elif(coorsh==3):
            tdata=np.copy(tsdata)   
            particles = np.zeros((nstar, 21))
            particles[:,0:7]=tdata[:,0:7]
            particles[:,10]=tdata[:,7]
            particles[:,19]=tdata[:,8]
            particles[:,9]=tdata[:,9]
            particles[:,11]=tdata[:,10]
            Xres=Xstar
        else:
            print("INVALID CHOICE")
            exit
       

        find_labels=["BFieldx", "BFieldy", "BFieldz", "wx", "wy", "wz"]
        auxarray = np.zeros((particles.shape[0], len(find_labels)))
        # Map labels to indices
        label_to_index = {label: i for i, label in enumerate(labels)}

        # Iterate over find_labels to update BField and wfield
        for i, label in enumerate(find_labels):
            if label in label_to_index:
                if i < 3:  # First three are BFieldx, BFieldy, BFieldz
                    auxarray[:, i] = Xres[:, label_to_index[label]]
        
        particles[:,12:18] = auxarray


        coorsh2 = int(input("Split/merge on condition: 1: none 2: mass 3: h \n "))
        
        coorshms = int(input("Do you want to split, merge or both? both means that you split and then merge. 1:split 2:merge 3:merge+split  \n "))

        scond=None
        mcond=None
        if(coorsh2==1):
            print("No condition")
        else:
            if(coorsh2==2):
                print("Mass condition, if only splitting then enter target max mass, if merging then enter target min mass. If both then you will set max and min.")
                if(coorshms==1 or coorshms==3):
                    maxmass=float(input("Target maximum mass \n "))
                    scond = lambda p: p[:, 0] > maxmass
                if(coorshms==2 or coorshms==3):
                    minmass=float(input("Target minimum mass \n "))
                    mcond = lambda p: p[:, 0] < minmass
            else:
                print("smoothing length condition")
                if(coorshms==1 or coorshms==3):
                    maxh=float(input("Target maximum h \n "))
                    scond = lambda p: p[:, 8] > maxh
                if(coorshms==2 or coorshms==3):
                    minh=float(input("Target minimum h \n "))
                    mcond = lambda p: p[:, 8] < minh
         
        #scond = lambda p: (p[:, 1] > 0.0) & (p[:, 2] > 0.0)
        #mcond = lambda p: (p[:, 1] > 0.0) & (p[:, 2] > 0.0)
        
        coorsh4 = float(input("Periodic or not? 1: True else: False "))
        if(coorsh4==1):
            periodic=True
            Lx, Ly, Lz = map(float, input("Enter box sizes (x, y, z) separated by spaces: ").split())
            box_size=np.array([Lx,Ly,Lz])
        else:
            periodic=False
            box_size=np.array([1.0,1.0,1.0])
        
        nt = int(input("How many times do you wanna split ? \n "))
        
        if(coorsh2==1):
            print("No condition, limit nsmooth to 3+2^nt ")
            #nsmooth=int(3+2**nt)
            nsmooth=64
        else:
            nsmooth=int(input(" nsmooth (typical 64) "))
        

        print(particles.size)
        tree,neighbors_list,distances,vec,particles =initial_setup(particles,box_size,nsmooth,periodic_box=periodic,h_init=None)
        print(particles.size)
        if(coorshms==1):
            result_particles = split_particles_full(particles, tree, neighbors_list, distances, vec ,periodic_box=periodic,box_size=box_size,condition=scond,split_rounds=nt)
        if(coorshms==2):
            result_particles, unmerged_particles = merge_particles_full(particles, tree, neighbors_list, distances, vec ,periodic_box=periodic,box_size=box_size,condition=mcond,merge_rounds=nt)
        if(coorshms==3):
            for _ in range(nt):
              particles_merged, unmerged_particles = merge_particles_full(particles, tree, neighbors_list, distances, vec ,periodic_box=periodic,box_size=box_size,condition=mcond,merge_rounds=1)
              tree_merged,neighbors_list_merged,distances_merged,vec_merged,particles_merged =initial_setup(particles_merged,box_size=box_size,nsmooth=nsmooth,periodic_box=periodic,h_init=particles_merged[:,20])
              particles = split_particles_full(particles_merged, tree_merged, neighbors_list_merged, distances_merged, vec_merged ,periodic_box=periodic,box_size=box_size,condition=scond,split_rounds=1)
              tree,neighbors_list,distances,vec,particles =initial_setup(particles,box_size,nsmooth=nsmooth,periodic_box=periodic,h_init=particles[:, 20])
            result_particles=particles
        tree,neighbors_list,distances,vec,result_particles =initial_setup(result_particles,box_size,nsmooth=64,periodic_box=periodic,h_init=result_particles[:, 20])
        
        if(coorsh==1):
            tgdata = np.zeros((result_particles.shape[0], tgdata.shape[1]))
            tgdata[:,0:12]=result_particles[:,0:12]
            ngas=len(tgdata[:,0])
            

            ## all aux values of gas set to zero first
            Xgas = np.zeros((result_particles.shape[0], len(labels)))
            BField=result_particles[:,12:15]
            wField=result_particles[:,15:18]
            
            find_labels=["BFieldx", "BFieldy", "BFieldz", "wx", "wy", "wz"]            
            # Map labels to indices
            label_to_index = {label: i for i, label in enumerate(labels)}

            # Iterate over find_labels to update BField and wfield
            for i, label in enumerate(find_labels):
                if label in label_to_index:
                     if i < 3:  # First three are BFieldx, BFieldy, BFieldz
                         Xgas[:, label_to_index[label]] = BField[:, i]
                     else:  # Next three are wx, wy, wz
                         Xgas[:, label_to_index[label]] = wField[:, i - 3]            
                else:
                    # Label not found, so we add it to labels and Xgas
                    labels.append(label)
                    new_index = len(labels) - 1
                    label_to_index[label] = new_index
                    #zeroarrstar = np.zeros((nstar, 1))
                    #zeroarrdark = np.zeros((ndark, 1))
                    if i < 3:
                        new_data = BField[:, i]   # shape (nParticles,)
                    else:
                        new_data = wField[:, i - 3]  # shape (nParticles,)
                     # Expand Xgas by adding a new column
                     # new_data must be shape (nParticles,)
                    print(Xgas.shape)
                    Xgas = np.column_stack((Xgas, new_data))
                    print(Xgas.shape)
                    #Xstar = np.column_stack((Xstar, zeroarrstar))
                    #Xdark = np.column_stack((Xdark, zeroarrdark))
        elif coorsh == 2:
            tddata = np.zeros((result_particles.shape[0], tddata.shape[1]))
            tddata[:,0:7]=result_particles[:,0:7]
            tddata[:,7]=result_particles[:,9]
            tddata[:,8]=result_particles[:,10]
            ndark=len(tddata[:,0])
            ## all aux values of dark set to zero first
            Xdark = np.zeros((result_particles.shape[0], len(labels)))

        elif(coorsh==3):
            tsdata = np.zeros((result_particles.shape[0], tsdata.shape[1]))
            tsdata[:,0:7]=result_particles[:,0:7]
            tsdata[:,7]=result_particles[:,10]
            tsdata[:,8]=result_particles[:,19]
            tsdata[:,9]=result_particles[:,9]
            tsdata[:,10]=result_particles[:,11]
            nstar=len(tsdata[:,0])
            ## all aux values of star set to zero first
            Xstar = np.zeros((result_particles.shape[0], len(labels) ))

        else:
            print("INVALID CHOICE")
            exit  

        data_header[0] = ngas + ndark + nstar
        data_header[2] = ngas
        data_header[3] = ndark
        data_header[4] = nstar
        ncols_final = Xgas.shape[1]
        Xstar = match_columns_with_zeros(ncols_final, Xstar)
        Xdark = match_columns_with_zeros(ncols_final, Xdark)
        print(Xgas.shape,Xdark.shape,Xstar.shape, X.shape)
        X = np.concatenate((Xgas, Xdark, Xstar))

        #use labels to find corresponding indice for BField,Spin,BClean. If it does not exist for Spin(wx,wy,wz) generate arrays
        # Find label that equal BFieldx,BFieldy,BFieldz,BField,BClean,wx,wy,wz
    if (choice==19):
        tdata=np.copy(tgdata)
        result_particles = np.zeros((ngas, 21))
        result_particles[:,0:12]=tdata[:,0:12]
        print(ngas)
        periodic=False
        box_size=np.array([1.0,1.0,1.0])
        tree,neighbors_list,distances,vec,result_particles =initial_setup(result_particles,box_size,nsmooth=64,periodic_box=periodic,h_init=None)
        tgdata = np.zeros((result_particles.shape[0], tgdata.shape[1]))
        tgdata[:,0:12]=result_particles[:,0:12]
        ngas=len(tgdata[:,0])
        data_header[0] = ngas + ndark + nstar
        data_header[2] = ngas
        data_header[3] = ndark
        data_header[4] = nstar

        print(ngas)

    tip.writetipsy(tgdata,tddata,tsdata,output,data_header,time)
    tip.writealltipsyauxB(X, output, labels)
    return choice

entry = input("File name: ")
output = input("Output name: ")
choice=setQuantity(entry,output)
