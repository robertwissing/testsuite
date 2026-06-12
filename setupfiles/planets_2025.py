# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 04:01:33 2016

@author: Khyzath
"""
import sys
import scipy.integrate as integrate
import scipy.constants as constants
import numpy as np
import scipy.optimize as optimize
from varpoly import varpoly

# Materials (edited)
# (rho0,B0,dBdP0,AA,ZZ,tmin,delta,a_exp,c_exp)
iron = [7744,166e9,5.1,55,26,1.8,2,4.92,5./3.,1.3] # a_exp and c_exp uncertain
outer_core = [6920,115e9,5.4,55,26,1.8,2,4.445,5./3.,0.8]
forsterite = [3380,130e9,4.2,36,18,2.2,13,5.615,5./3.,0.64] #Fitted to vapor curve for forsterite

material_list=[]

def planet(name):
    if name=="Earth_PES095_2layer":
        file_string = name + '.dat'
        p_00 = 3.54e11 # Central pressure
        T_c = 5400 #Central core temperature
        r_end = [0.33800E+07, 9.9e9] #end of each layer, last is dummy variable
        T_cond = [0,2000,0] # If there is a temperature jump between two layers, the starting temperature will be this for the given layer(if > 0).
        EOSobj_1 = varpoly();
        EOSobj_1.creatematerial(*iron)
        EOSobj_2 = varpoly();
        EOSobj_2.creatematerial(*forsterite)
        material_list.append(EOSobj_1)
        material_list.append(EOSobj_2)
        return [file_string, p_00, T_c, r_end, T_cond, material_list]
    if name=="Earth_PES095_2layer_warm":
        file_string = name + '.dat'
        p_00 = 3.54e11 # Central pressure
        T_c = 20000 #Central core temperature
        r_end = [0.33800E+07, 9.9e9] #end of each layer, last is dummy variable
        T_cond = [0,0,0] # If there is a temperature jump between two layers, the starting temperature will be this for the given layer(if > 0).
        EOSobj_1 = varpoly();
        EOSobj_1.creatematerial(*iron)
        EOSobj_2 = varpoly();
        EOSobj_2.creatematerial(*forsterite)
        material_list.append(EOSobj_1)
        material_list.append(EOSobj_2)
        return [file_string, p_00, T_c, r_end, T_cond, material_list]
    elif name=="testplanet_1layer":
        file_string = 'testplanet_1layer.dat'
        p_00 = 3.54e11 # Central pressure
        T_c = 5400 #Central core temperature
        r_end = [9.9e9] #end of each layer, last is dummy variable
        T_cond = [0,0] # If there is a temperature jump between two layers, the starting temperature will be this for the given layer(if > 0).
        EOSobj_1 = varpoly();
        EOSobj_1.creatematerial(*forsterite)
        material_list.append(EOSobj_1)
        return [file_string, p_00, T_c, r_end, T_cond, material_list]
    elif name=="testplanet_cold_2layer":
        file_string = 'testplanet_cold_2layer.dat'
        p_00 = 3.54e11 # Central pressure
        T_c = 100 #Central core temperature
        r_end = [0.33800E+07, 9.9e9] #end of each layer, last is dummy variable
        T_cond = [0,2000,0] # If there is a temperature jump between two layers, the starting temperature will be this for the given layer(if > 0).
        EOSobj_1 = varpoly();
        EOSobj_1.creatematerial(*iron)
        EOSobj_2 = varpoly();
        EOSobj_2.creatematerial(*forsterite)
        material_list.append(EOSobj_1)
        material_list.append(EOSobj_2)
        return [file_string, p_00, T_c, r_end, T_cond, material_list]
    elif name=="testplanet_cold_1layer":
        file_string = 'testplanet_cold_1layer.dat'
        p_00 = 3.54e11 # Central pressure
        T_c = 100 #Central core temperature
        r_end = [9.9e9] #end of each layer, last is dummy variable
        T_cond = [0,0] # If there is a temperature jump between two layers, the starting temperature will be this for the given layer(if > 0).
        EOSobj_1 = varpoly();
        EOSobj_1.creatematerial(*forsterite)
        material_list.append(EOSobj_1)
        return [file_string, p_00, T_c, r_end, T_cond, material_list]
    elif name=="theia13ms":
        file_string = name + '.dat'
        p_00 = 0.719e11
        T_c = 1800
        r_end = [0.19E+07,0.400E+07, 9.9e9]
        T_cond = [0,0,0,0]
        EOSobj_1 = varpoly();
        EOSobj_1.creatematerial(*iron)
        EOSobj_2 = varpoly();
        EOSobj_2.creatematerial(*forsterite)
        material_list.append(EOSobj_1)
        material_list.append(EOSobj_2)
        return [file_string, p_00, T_c, r_end, T_cond, material_list]
    else:
        print("No planet")
    
