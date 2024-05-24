#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 13:31:40 2018

@author: robertwi
"""

import numpy as np
import random as random
import glob, os
import fnmatch
import sys

# Define data types
di = np.dtype('>i4')
df = np.dtype('>f')
dd = np.dtype('>d')


intauxfiles = ["iord","igasorder"]
banauxfiles = ["AHF","grid"]
vecauxfiles = ["BField","CurlB"]


def readtipsy(file_string,endian=1):
    fB=open(file_string,'rb')
    time=np.fromfile(fB,dtype=dd,count=1)
    data_header=np.fromfile(fB,dtype=di,count=6)
    tipsygdata=np.fromfile(fB,dtype=df,count=data_header[2]*12)
    tipsygdata=np.reshape(tipsygdata,(-1,12))
    tipsyddata=np.fromfile(fB,dtype=df,count=data_header[3]*9)
    tipsyddata=np.reshape(tipsyddata,(-1,9))
    tipsysdata=np.fromfile(fB,dtype=df,count=data_header[4]*11)
    tipsysdata=np.reshape(tipsysdata,(-1,11))
    fB.close()
    return tipsygdata,tipsyddata,tipsysdata,data_header,time
    
def readtipsyaux(file_stringaux,label,binary=1,endian=1):
    file_stringaux=file_stringaux+'.'+label
    fB=open(file_stringaux,'rb')
    npart=np.fromfile(fB,dtype=di,count=1)
    print(file_stringaux,npart)
    if(label in intauxfiles):
        X=np.fromfile(fB,dtype=di)
    else:    
        X=np.fromfile(fB,dtype=df)
    fB.close()
    X=fixendian(X,endian)
    return X,label

def readtipsyauxB(file_stringaux,label,ngas,ndark,nstar,binary=1,endian=1):
    file_stringaux=file_stringaux+'.'+label
    fB=open(file_stringaux,'rb')
    npart=np.fromfile(fB,dtype=di,count=1)
    print(file_stringaux,npart)
    if(label in intauxfiles):
        Xg=np.fromfile(fB,dtype=di,count=ngas)
        Xd=np.fromfile(fB,dtype=di,count=ndark)
        Xs=np.fromfile(fB,dtype=di,count=nstar)
    else:
        Xg=np.fromfile(fB,dtype=df,count=ngas)
        Xd=np.fromfile(fB,dtype=df,count=ndark)
        Xs=np.fromfile(fB,dtype=df,count=nstar)
    fB.close()
    Xg=fixendian(Xg,endian)
    Xd=fixendian(Xd,endian)
    Xs=fixendian(Xs,endian)
    return Xg,Xd,Xs,label

def readtipsyVecaux(file_stringaux,label,binary=1,endian=1):
    file_stringaux=file_stringaux+'.'+label
    fB=open(file_stringaux,'rb')
    npart=np.fromfile(fB,dtype=di,count=1)
    print(file_stringaux,npart)
    if(label in intauxfiles):
        X=np.fromfile(fB,dtype=di)
    else:    
        X=np.fromfile(fB,dtype=df)
    fB.close()
    X=fixendian(X,endian)
    return X,label
    
def readalltipsyaux(file_string,binary=1,endian=1):
    file_stringaux=glob.glob(file_string + ".*")
    file_stringaux=sorted(file_stringaux, key=str.lower)
    j=0
    vec=0
    while j <len(file_stringaux):
        print(file_string)
        print(file_stringaux[j].replace(file_string + '.',''))
        if any(substring == file_stringaux[j].replace(file_string + '.','') for substring in vecauxfiles):
            print("Vector field: " + file_stringaux[j])
            vec = vec + 2
        elif all(substring not in file_stringaux[j].replace(file_string + '.','') for substring in banauxfiles):
            pass
        else:
            print("Did not include: " + file_stringaux[j])
            del file_stringaux[j]
            j=j-1
        j=j+1
    #label=file_stringaux.copy()
    #label=file_stringaux[:]
    label=[]
    naux=len(file_stringaux)
    if naux == 0:
        X=[]
        label=[]
    else: 
        k=0
        for i in range (naux):
                print(file_stringaux[i])
                vecread=0
                if any(substring == file_stringaux[i].replace(file_string + '.','') for substring in vecauxfiles):
                    print("Vector field: " + file_stringaux[i])
                    label.append(file_stringaux[i].replace(file_string + '.','')+ "x")
                    label.append(file_stringaux[i].replace(file_string + '.','')+ "y")
                    label.append(file_stringaux[i].replace(file_string + '.','') + "z")
                    vecread=1
                else:
                    label.append(file_stringaux[i].replace(file_string + '.',''))
                    
                fB=open(file_stringaux[i],'rb')
                npart=np.fromfile(fB,dtype=di,count=1)
                print(npart,naux)
                if(i==0):
                    X=np.zeros((npart[0],naux+vec))
                if(label[i] in intauxfiles):
                    X[:,k]=np.fromfile(fB,dtype=di)
                else:
                    if vecread == 1:
                        X[:,k:k+3]=(np.fromfile(fB,dtype=df)).reshape((-1, 3))
                        k=k+2
                    else:
                        X[:,k]=np.fromfile(fB,dtype=df)
                k=k+1
                fB.close()
        X=fixendian(X,endian)
    return X,label

def readalltipsy(file_string,binary=1,endian=1):
    tipsygdata,tipsyddata,tipsysdata,data_header,time=readtipsy(file_string)
    tipsydataaux,label=readalltipsyaux(file_string)
    ng=len(tipsygdata[:,0])
    nd=len(tipsyddata[:,0])
    ns=len(tipsysdata[:,0])
    tipsygdata = np.column_stack((tipsygdata,tipsydataaux[0:ng,:]))
    tipsyddata = np.column_stack((tipsyddata,tipsydataaux[ng:(ng+nd),:]))
    tipsysdata = np.column_stack((tipsysdata,tipsydataaux[(ng+nd):(ng+nd+ns),:]))
    tipsygdata=fixendian(tipsygdata,endian)
    tipsyddata=fixendian(tipsyddata,endian)
    tipsysdata=fixendian(tipsysdata,endian)
    return tipsygdata,tipsyddata,tipsysdata,data_header,time,label
    
def writetipsyaux(X,label,file_string,binary=1,endian=1):  
    file_stringB=file_string + "." + label
    fB=open(file_stringB,'w');
    fB=open(file_stringB,'a');
    B=df.type(X)
    npart=len(X)
    bnpart=di.type(npart)
    if(binary==0):
        fB.write(str(npart))
        fB.write("\n")
        np.savetxt(fB,B);
        fB.close()
    else:
        B=fixendian(B)
        bnpart=bnpart.byteswap()
        bnpart.tofile(fB)
        B.tofile(fB)
        fB.close()
        fB=open(file_stringB,'rb')     
        fB.close()
        
def writetipsy(tipsygdata,tipsyddata,tipsysdata,file_string,data_header,time=np.array([0.0]),endian=1):
    fbin=open(file_string,'wb');
    fbin=open(file_string,'ab');
    tipsygdata=df.type(tipsygdata)
    tipsyddata=df.type(tipsyddata)
    tipsysdata=df.type(tipsysdata)
    data_header=di.type(data_header)
    time=fixendian(time)
    data_header=fixendian(data_header,endian)
    tipsygdata=fixendian(tipsygdata,endian)
    tipsyddata=fixendian(tipsyddata,endian)
    tipsysdata=fixendian(tipsysdata,endian)
    time.tofile(fbin)
    data_header.tofile(fbin)
    tipsygdata.tofile(fbin)
    tipsyddata.tofile(fbin)
    tipsysdata.tofile(fbin)
    fbin.close()

def writealltipsyaux(X,file_string,case="mhd",binary=1,endian=1):
    if case=="mhd":
        writetipsyaux(X[:,0],"BFieldx",file_string)
        writetipsyaux(X[:,1],"BFieldy",file_string)
        writetipsyaux(X[:,2],"BFieldz",file_string)
        writetipsyVecaux(X,"BField",file_string)
    if case=="spin":
        writetipsyaux(X[:,0],"wx",file_string)
        writetipsyaux(X[:,1],"wy",file_string)
        writetipsyaux(X[:,2],"wz",file_string)

def writealltipsyauxB(X,file_string,label,binary=1,endian=1):
    i=0
    for lab in label:
        writetipsyaux(X[:,i],lab,file_string,binary=1,endian=1)
        i=i+1

def writeprepASCII(tipsydata,file_string,data_header,time=np.array([0.0]),endian=1):
    fbin=open(file_string,'wb');
    fbin=open(file_string,'ab');
    tipsydata=df.type(tipsydata)
    data_header=di.type(data_header)
    time=fixendian(time)
    data_header=fixendian(data_header,endian)
    tipsydata=fixendian(tipsydata,endian)
    np.savetxt(fbin,data_header, fmt='%i',delimiter=' ');
    np.savetxt(fbin,tipsydata, fmt='%f',delimiter=' ');
    fbin.close()
    
    
def writetipsyVecaux(X,label,file_string,binary=1,endian=1):  
        file_stringB=file_string + "." + label
        fB=open(file_stringB,'w');
        fB=open(file_stringB,'a');
        B=df.type(np.dstack((X[:,0],X[:,1],X[:,2])))
        npart=len(X[:,0])
        bnpart=di.type(npart)
        if(binary==0):
            B=df.type(np.hstack((X[:,0],X[:,1],X[:,2])))
            fB.write(str(npart))
            fB.write("\n")
            np.savetxt(fB,B);
            fB.close()
        else:
            B=B.byteswap().newbyteorder()
            bnpart=bnpart.byteswap()
            bnpart.tofile(fB)
            B.tofile(fB)
            fB.close()
            fB=open(file_stringB,'rb')     
            fB.close()
            
def createtipsysplash(file_string,binary=1,endian=1,centering=0):
    tipsygdata,tipsyddata,tipsysdata,data_header,time,label=readalltipsy(file_string)
    file_stringsplash=file_string + "_splash"
    ##Center to gas potential
    if (centering==1):
        print("centering");
        minpotg=min(tipsygdata[:,11])
        minpotidx=np.argmin(tipsygdata[:,11])
        rmin=tipsygdata[minpotidx,1:4]
        vmin=tipsygdata[minpotidx,4:7]
        npart = len(tipsygdata[:,0])+len(tipsyddata[:,0])+len(tipsysdata[:,0])
        M = (np.sum(tipsygdata[:,0])+np.sum(tipsyddata[:,0])+np.sum(tipsysdata[:,0]))
        print(npart)
        cmx = (np.sum(tipsygdata[:,1]*tipsygdata[:,0])+np.sum(tipsyddata[:,1]*tipsyddata[:,0])+np.sum(tipsysdata[:,1]*tipsysdata[:,0]))/M
        cmy = (np.sum(tipsygdata[:,2]*tipsygdata[:,0])+np.sum(tipsyddata[:,2]*tipsyddata[:,0])+np.sum(tipsysdata[:,2]*tipsysdata[:,0]))/M
        cmz = (np.sum(tipsygdata[:,3]*tipsygdata[:,0])+np.sum(tipsyddata[:,3]*tipsyddata[:,0])+np.sum(tipsysdata[:,3]*tipsysdata[:,0]))/M
        vcmx = (np.sum(tipsygdata[:,4]*tipsygdata[:,0])+np.sum(tipsyddata[:,4]*tipsyddata[:,0])+np.sum(tipsysdata[:,4]*tipsysdata[:,0]))/M
        vcmy = (np.sum(tipsygdata[:,5]*tipsygdata[:,0])+np.sum(tipsyddata[:,5]*tipsyddata[:,0])+np.sum(tipsysdata[:,5]*tipsysdata[:,0]))/M
        vcmz = (np.sum(tipsygdata[:,6]*tipsygdata[:,0])+np.sum(tipsyddata[:,6]*tipsyddata[:,0])+np.sum(tipsysdata[:,6]*tipsysdata[:,0]))/M
        print("npart ",npart,"cmx ",cmx,"cmy ",cmy,"cmz ",cmz,"vcmx ",vcmx,"vcmy ",vcmy,"vcmz ",vcmz )
        print("mass",tipsygdata[:,0])
        #cmx=0.0;cmy=0.0;cmz=0.0;vcmx=0.0;vcmy=0.0;vcmz=0.0;
        vcmx=0.0;vcmy=0.0;vcmz=0.0;
        tipsyddata[:,1]=tipsyddata[:,1]-cmx
        tipsyddata[:,2]=tipsyddata[:,2]-cmy
        tipsyddata[:,3]=tipsyddata[:,3]-cmz
        tipsyddata[:,4]=tipsyddata[:,4]-vcmx
        tipsyddata[:,5]=tipsyddata[:,5]-vcmy
        tipsyddata[:,6]=tipsyddata[:,6]-vcmz
        tipsygdata[:,1]=tipsygdata[:,1]-cmx
        tipsygdata[:,2]=tipsygdata[:,2]-cmy
        tipsygdata[:,3]=tipsygdata[:,3]-cmz
        tipsygdata[:,4]=tipsygdata[:,4]-vcmx
        tipsygdata[:,5]=tipsygdata[:,5]-vcmy
        tipsygdata[:,6]=tipsygdata[:,6]-vcmz
        tipsysdata[:,1]=tipsysdata[:,1]-cmx
        tipsysdata[:,2]=tipsysdata[:,2]-cmy
        tipsysdata[:,3]=tipsysdata[:,3]-cmz
        tipsysdata[:,4]=tipsysdata[:,4]-vcmx
        tipsysdata[:,5]=tipsysdata[:,5]-vcmy
        tipsysdata[:,6]=tipsysdata[:,6]-vcmz
    if (centering==2):
        print("centering");
        minpotg=min(tipsygdata[:,11])
        minpotidx=np.argmin(tipsygdata[:,11])
        rmin=tipsygdata[minpotidx,1:4]
        vmin=tipsygdata[minpotidx,4:7]
        tipsyddata[:,1:4]=tipsyddata[:,1:4]-rmin
        tipsysdata[:,1:4]=tipsysdata[:,1:4]-rmin
        tipsygdata[:,1:4]=tipsygdata[:,1:4]-rmin
        tipsyddata[:,4:7]=tipsyddata[:,4:7]
        tipsysdata[:,4:7]=tipsysdata[:,4:7]
        tipsygdata[:,4:7]=tipsygdata[:,4:7]

    #tispyddata=np.array([])
    #tispysdata=np.array([])
    #data_header[0]=data_header[0]-data_header[3]-data_header[4]
    #data_header[3]=0
    #data_header[4]=0
    writetipsy(tipsygdata,tipsyddata,tipsysdata,file_stringsplash,data_header,time)
    file_stringlabel=file_stringsplash + ".label"
    fB=open(file_stringlabel,'w');
    fB=open(file_stringlabel,'a');
    fB.write(str(len(label)))
    fB.write("\n")
    fB.write("\n".join(label))
    fB.close()
    return tipsygdata,tipsyddata,tipsysdata
    
def readgrid(file_string):
    fB=open(file_string,'rb')
    headid=np.fromfile(fB,dtype=di,count=1)
    cell=np.fromfile(fB,dtype=di,count=3)
    ntot=cell[0]*cell[1]*cell[2]
    ncol=np.fromfile(fB,dtype=di,count=1)
    ncol=ncol[0]
    time=np.fromfile(fB,dtype=df,count=1)
    headid=np.fromfile(fB,dtype=di,count=1)
    tipsydata=np.zeros((ntot,ncol),dtype=df)
    
    for i in range(0,ncol):
        bodid=np.fromfile(fB,dtype=di,count=1)
        tipsydata[:,i]=np.fromfile(fB,dtype=df,count=ntot)
        bodid=np.fromfile(fB,dtype=di,count=1)
    nx, ny, nz = (cell[0], cell[1], cell[2])
    griddata=np.zeros((len(tipsydata[0,:]),nx,ny,nz),dtype=df)
    for j in range(0,len(tipsydata[0,:])):
        datai=np.reshape(tipsydata[:,j],(nz,ny,nx))
        datai=np.swapaxes(datai,0,2)
        griddata[j,:,:,:]=datai

    return griddata,cell,time

def getfiles(endtag="[0-9]"):
    filenames=[]
    listOfFiles = os.listdir('.')
    for entryinit in listOfFiles:
        if fnmatch.fnmatch(entryinit, "*.*"+endtag):
            filenames.append(entryinit)
    if(len(filenames)==0):
        print("No file ending with " + endtag + " in directory")
        sys.exit()
    return filenames

def getgridall(listOfFiles,replace=0):
        tipsydata,cell,time =readgrid(listOfFiles[0])
        nx, ny, nz = (cell[0], cell[1], cell[2])
        if(os.path.isfile("GridDataAllDumps.npy")):      
            data=np.load("GridDataAllDumps.npy")
            print(ny,len(data[0,0,0,:,0]))
        if(os.path.isfile("GridDataAllDumps.npy" or "TimeData.npy") and len(data[0,0,:,0,0])==nx and len(data[0,0,0,:,0]) == ny and len(data[0,0,0,0,:]) == nz and replace==0 and len(data[0,:,0,0,0]) == len(listOfFiles)):
            data=np.load("GridDataAllDumps.npy")
            timedata=np.load("TimeData.npy")
        else:
            nx, ny, nz = (cell[0], cell[1], cell[2])
            nfiles=len(listOfFiles)
            data=np.zeros((len(tipsydata[:,0,0,0]),nfiles,nx,ny,nz),dtype=df)
            timedata=np.zeros(nfiles,dtype=df)
            i=0
            for entry in listOfFiles:
                tipsydata,cell,time=readgrid(entry)
                nx, ny, nz = (cell[0], cell[1], cell[2])
                print(nx,ny,nz,entry)
                data[:,i,:,:,:]=tipsydata
                timedata[i]=time
                i=i+1
            np.save("GridDataAllDumps.npy",data)
            np.save("TimeData.npy",timedata)
                
        return data,timedata
            
def fixendian(arr,endian=1):
    if(arr.dtype.byteorder == "=" or arr.dtype.byteorder == "<" and endian==1):
        arr=arr.byteswap().newbyteorder()
    elif(arr.dtype.byteorder == ">" and endian==0):
        arr=arr.byteswap().newbyteorder()
    return arr

def convertfromcodeunits(D2,dKpcUnit,dMsolUnit,Tin,vin,Bin):
    MSOLG=1.99e33
    gtoMsol=1/MSOLG
    KBOLTZ=1.38e-16
    MHYDR=1.67e-24
    KPCCM=3.085678e21
    cmtoKpc=1/KPCCM
    GCGS=6.67e-8
    
    #permiability
    #mu=kg*m/(s^2*A^2) A=sqrt((kg*m)/(mu))*s 10^-4*sqrt(kg*mu)/(sqrt(m)*s^3)
    #G=10^-4kg/(A*s^2)
    #INPUTS ARE CODE UNITS IN GASOLINE - 
    #CODE UNIT SET BY dKPC dMsol
    # if dKPC 1 dMsol 1 convert all inputs kpc and msol
    
    dConstGamma = D2.gamma
    dMeanMolWeight = 1.
    dGasConst = dKpcUnit*KPCCM*KBOLTZ/(MHYDR*GCGS*dMsolUnit*MSOLG);
    dErgPerGmUnit = GCGS*dMsolUnit*MSOLG/(dKpcUnit*KPCCM);
    dGmPerCcUnit = (dMsolUnit*MSOLG)/(dKpcUnit*KPCCM)**3
    dSecUnit= np.sqrt(1/(dGmPerCcUnit*GCGS))
    ##
    if(Tin==0):
        duTFac = 1
    else:
        duTFac = (dConstGamma - 1)*dMeanMolWeight/dGasConst 

    if(vin==0):
        vFac = 1.0
    else:        
        vFac = dSecUnit

    if(Bin==0):    
        Bfac=1.0
        CodetoGauss=1.0
    else:
        MU0CGS=4*np.pi
        GausstoCode=np.sqrt((gtoMsol/dMsolUnit)/((cmtoKpc/dKpcUnit)*MU0CGS))*dSecUnit
        Bfac=GausstoCode
        CodetoGauss=1/GausstoCode
    print("BFAC: ", Bfac, "vFac: ", vFac, "duTFac", duTFac, "dSecunit ", dSecUnit, "dBunit ", CodetoGauss)
   
    D2.vx=np.multiply(D2.vx,vFac)
    D2.vy=np.multiply(D2.vy,vFac)  
    D2.vz=np.multiply(D2.vz,vFac)
    print(D2.vz)
    D2.u=np.multiply(D2.u,duTFac)
    D2.Bx=np.multiply(D2.Bx,Bfac)
    D2.By=np.multiply(D2.By,Bfac)
    D2.Bz=np.multiply(D2.Bz,Bfac)
    return D2
