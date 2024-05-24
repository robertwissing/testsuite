#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 15:27:28 2017

@author: robertwi
"""
import subprocess
import numpy as np

import struct
import pickle
import math

def writeascii(file_string,D,mhd):
    asciihere=0
    ##with magneticfield
    file_stringB=file_string+"_B"
    di =np.dtype('>i4')
    df=np.dtype('>f')
    dd=np.dtype('>d')
    bpad=di.type(0)
    bndark=di.type(int(D.npart-D.ngas-D.nstar))
    btime=dd.type(D.time)
    bnpart=di.type(D.npart)

    data_1=np.column_stack((int(D.npart),int(D.ngas),int(D.nstar)));
    data_out = np.stack((int(D.ndim),float(D.time)));
    data_out2 = np.hstack((D.mass,D.x,D.y,D.z,D.vx,D.vy,D.vz,D.soft,D.rho,D.u,D.h,D.metals,D.tform,D.pot));
    data_outB = df.type(np.hstack((D.mass,D.x,D.y,D.z,D.vx,D.vy,D.vz,D.soft,D.rho,D.u,D.h,D.metals,D.tform,D.pot,D.Bx,D.By,D.Bz,D.Bpsi)));

    data_header = di.type(np.hstack((D.npart,D.ndim,D.ngas,bndark,D.nstar,bpad)))
    data_outorg = df.type(np.dstack((D.mass,D.x,D.y,D.z,D.vx,D.vy,D.vz,D.rho,D.u,D.h,D.metals,D.pot)));
    data_outB2 = df.type(np.dstack((D.mass,D.x,D.y,D.z,D.vx,D.vy,D.vz,D.rho,D.u,D.h,D.metals,D.pot,D.Bx,D.By,D.Bz,D.Bpsi)))
    data_headerlittle=data_header

    btime=btime.byteswap()
    data_header=data_header.byteswap().newbyteorder()
    data_outorg=data_outorg.byteswap().newbyteorder()
    asda=Asddas
    if(mhd==2):
        file_stringB=file_string+'bin'
        fB=open(file_stringB,'w');
        fB=open(file_stringB,'ab');
        btime.tofile(fB)
        data_header.tofile(fB)
        data_outB2.tofile(fB)
        fB.close()
        fB=open(file_stringB,'rb')
        #test1=np.fromfile(fB,dd.type,1)
        #test2=np.fromfile(fB,di.type,6)
        #test3=np.fromfile(fB,df.type,3)
        fB.close()
        print('hej',test1,test2)
        print('float',test3,D.mass[0])
        exefile = './bin2ascii '
        extra='<'+ file_string + 'bin> ' + file_string + 'ascii'
        args=exefile+extra
        print(args)
        cmd=subprocess.Popen(args,shell=True,stdout=subprocess.PIPE ,stderr=subprocess.PIPE)
        cmd_out, cmd_err = cmd.communicate()
        print(cmd_out)
        print(cmd_err)
        exefile = './totipstd '
        extra='<'+ file_string + 'bin> ' + file_string + 'std'
        args=exefile+extra
        cmd=subprocess.Popen(args,shell=True,stdout=subprocess.PIPE ,stderr=subprocess.PIPE)
        cmd_out, cmd_err = cmd.communicate()
        print(cmd_out)
        print(cmd_err)
    elif(mhd==1):
        fB=open(file_stringB,'w');
        fB=open(file_stringB,'ab');
        np.savetxt(fB,data_1, fmt='%i',delimiter=', ');
        np.savetxt(fB, data_out, fmt='%i');
        np.savetxt(fB,data_outB);
        fB.close()


        fileextra=file_string+'_ex'
        fex=open(fileextra,'w');
        fex=open(fileextra,'ab');
        data_outex = np.hstack((D.dxbound, D.dybound, D.dzbound));
        np.savetxt(fex,data_outex, fmt='%1.8e',delimiter=', ');
        fex.close()


        exefile = './ascii2bin '
        extra='<'+ file_stringB + '> ' + file_string + 'bin'
        args=exefile+extra
        cmd=subprocess.Popen(args,shell=True,stdout=subprocess.PIPE ,stderr=subprocess.PIPE)
        cmd_out, cmd_err = cmd.communicate()
        print(cmd_out)
        print(cmd_err)

        fB=open(file_string + 'bin','rb')

        test1=np.fromfile(fB,np.float64,1)
        test2=np.fromfile(fB,np.intc,6)
        test3=np.fromfile(fB,np.float32,20)
        fB.close()
        print('hej',test1,test2)
        print('float',test3,D.mass[0])



        exefile = './bin2ascii '
        extra='<'+ file_string + 'bin> ' + file_string + 'ascii'
        args=exefile+extra
        cmd=subprocess.Popen(args,shell=True,stdout=subprocess.PIPE ,stderr=subprocess.PIPE)
        cmd_out, cmd_err = cmd.communicate()
        print(cmd_out)
        print(cmd_err)
        exefile = './totipstd '
        extra='<'+ file_string + 'bin> ' + file_string + 'std'
        args=exefile+extra
        cmd=subprocess.Popen(args,shell=True,stdout=subprocess.PIPE ,stderr=subprocess.PIPE)
        cmd_out, cmd_err = cmd.communicate()
        print(cmd_out)
        print(cmd_err)
    else:
        file_stringbin=file_string+'std'
        fbin=open(file_stringbin,'wb');
        fbin=open(file_stringbin,'ab');
        print(data_header.dtype.byteorder)
        print(data_header.tobytes()==data_headerlittle.tobytes())
        btime.tofile(fbin)
        data_header.tofile(fbin)
        data_outorg.tofile(fbin)
        fbin.close()

        fbin=open(file_stringbin,'rb')
        test1=np.fromfile(fbin,dtype=dd,count=1)
        test2=np.fromfile(fbin,dtype=di,count=6)
        test3=np.fromfile(fbin,dtype=df,count=4)
        fbin.close()
        print('hej',test1,test2)
        print('float',test3,D.mass[0])

        fileextra=file_string+'_ex'
        fex=open(fileextra,'w');
        fex=open(fileextra,'ab');
        data_outex = np.hstack((D.dxbound, D.dybound, D.dzbound));
        np.savetxt(fex,data_outex, fmt='%1.8e',delimiter=', ');
        fex.close()


        file_stringB=file_string + "std.BField"
        fB=open(file_stringB,'w');
        fB=open(file_stringB,'a');
        B=df.type(np.dstack((D.Bx,D.By,D.Bz)))
        print(D.Bx[1],D.By[1],D.Bz[1])
        if(asciihere==1):
            B=df.type(np.hstack((D.Bx,D.By,D.Bz)))
            fB.write(str(D.npart))
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
            test1=np.fromfile(fB,dtype=di,count=1)
            test2=np.fromfile(fB,dtype=df,count=4)
            fB.close()

        file_stringB=file_string + "std.BFieldx"
        fB=open(file_stringB,'w');
        fB=open(file_stringB,'a');
        B=df.type(D.Bx)
        if(asciihere==1):
            fB.write(str(D.npart))
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

        file_stringB=file_string + "std.BFieldy"
        fB=open(file_stringB,'w');
        fB=open(file_stringB,'a');
        B=df.type(D.By)
        if(asciihere==1):
            fB.write(str(D.npart))
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

        file_stringB=file_string + "std.BFieldz"
        fB=open(file_stringB,'w');
        fB=open(file_stringB,'a');
        B=df.type(D.Bz)
        if(asciihere==1):
            fB.write(str(D.npart))
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

        #exefile = './totipstd '
        #extra='<'+ file_string + 'orgbin> ' + file_string + 'std'
        #args=exefile+extra
        #cmd=subprocess.Popen(args,shell=True,stdout=subprocess.PIPE ,stderr=subprocess.PIPE)
        #cmd_out, cmd_err = cmd.communicate()
        #print(cmd_out)
        #print(cmd_err)

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
    dMeanMolWeight = D2.molweight
    dGasConst = dKpcUnit*KPCCM*KBOLTZ/(MHYDR*GCGS*dMsolUnit*MSOLG);
    print("gas and mol",dGasConst," ",dMeanMolWeight,"gamma",dConstGamma)
    dErgPerGmUnit = GCGS*dMsolUnit*MSOLG/(dKpcUnit*KPCCM);
    dGmPerCcUnit = (dMsolUnit*MSOLG)/(dKpcUnit*KPCCM)**3
    dSecUnit= np.sqrt(1/(dGmPerCcUnit*GCGS))
    ##
    if(Tin==0):
        duTFac = 1.0
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


def calculatecosmoMsol(h0,dKpcUnit):
    MSOLG=1.9891e33
    gtoMsol=1/MSOLG
    KBOLTZ=1.38e-16
    MHYDR=1.67e-24
    kpctocm=3.085678e21
    cmtoKpc=1/kpctocm
    GCGS=6.67259e-8

    H0code = np.sqrt(8.*np.pi/3.);
    H0phys = h0*100. * 1.0e5/(1.0e3 * kpctocm);
    dcode_to_phys = dKpcUnit * kpctocm;
    vcode_to_phys = H0phys/H0code * dcode_to_phys;
    mcode_to_phys = (1.0/GCGS)*dcode_to_phys*(vcode_to_phys**2.0);
    rhocode_to_phys = (1.0/GCGS)/dcode_to_phys**2.0 * vcode_to_phys**2.0;
    Msolunit=mcode_to_phys/MSOLG
    print( "dcode to phys", dcode_to_phys)
    print("vcode to phys in km/s",vcode_to_phys/1.0e5)
    print("Msolunit", Msolunit)
    print( "dsecunit", dcode_to_phys/vcode_to_phys)
    print( "density unit in g/cc (comoving)", rhocode_to_phys)
    return Msolunit
