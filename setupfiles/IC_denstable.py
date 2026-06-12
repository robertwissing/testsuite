#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 24 11:30:33 2025

@author: robertwi
"""

import xdrlib
import numpy as np


def generate_density_table(D, dim=1, n=1000):
    if dim == 1:
        xtab = np.linspace(0, D.Rtab, n)
        densities = np.array([D.getrho(x) for x in xtab])
        print(xtab,densities)
        return (xtab,), densities
    
    elif dim == 2:
        xtab = np.linspace(-D.Rtab, D.Rtab, n)
        ytab = np.linspace(-D.Rtab, D.Rtab, n)
        xx, yy = np.meshgrid(xtab, ytab, indexing='ij')
        densities = np.zeros((n, n))
        
        for i in range(n):
            for j in range(n):
                r = np.sqrt(xx[i,j]**2 + yy[i,j]**2)
                densities[i,j] = D.getrho(r)
                
        return (xtab, ytab), densities
    
    elif dim == 3:
        xtab = np.linspace(-D.Rtab, D.Rtab, n)
        ytab = np.linspace(-D.Rtab, D.Rtab, n)
        ztab = np.linspace(-D.Rtab, D.Rtab, n)
        densities = np.zeros((n, n, n))
        
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    r = np.sqrt(xtab[i]**2 + ytab[j]**2 + ztab[k]**2)
                    densities[i,j,k] = D.getrho(r)
        return (xtab, ytab, ztab), densities



def write_density_table_xdr(D, filename, grid, densities , dim=1, n=1000):

    with open(filename, 'wb') as f:
        # Write header: dimension, points per dim, domain size
        header = f"{dim}\n{n}\n{D.R:.15g}\n"
        f.write(header.encode('ascii'))
        
        # Create XDR packer
        packer = xdrlib.Packer()
        
        # Pack grid axes
        for axis in grid:
            packer.pack_array(axis.astype(np.float32), packer.pack_float)
        
        # Pack densities (flatten in C-order)
        flat_densities = densities.ravel().astype(np.float32)
        packer.pack_array(flat_densities, packer.pack_float)
        
        f.write(packer.get_buffer())
        
        
        
def read_density_table_xdr(filename):
    """
    Read density table from XDR file
    
    Parameters:
    filename : str
        Input file name
        
    Returns:
    tuple: (dim, n, R, grid, densities)
        dim : int
            Dimensionality (1, 2, or 3)
        n : int
            Points per dimension
        R : float
            Domain radius
        grid : tuple of arrays
            Coordinate axes (1-3 arrays)
        densities : ndarray
            Density values with shape (n,) for 1D, (n,n) for 2D, or (n,n,n) for 3D
    """
    with open(filename, 'rb') as f:
        # Read header
        header_lines = []
        for _ in range(3):
            line = b''
            while True:
                char = f.read(1)
                if char == b'\n' or char == b'':
                    break
                line += char
            header_lines.append(line.decode('ascii').strip())
        
        dim = int(header_lines[0])
        n = int(header_lines[1])
        R = float(header_lines[2])
        
        # Read XDR data
        data = f.read()
        unpacker = xdrlib.Unpacker(data)
        
        # Read grid axes
        grid = []
        for _ in range(dim):
            arr = unpacker.unpack_farray(n, unpacker.unpack_float)
            grid.append(np.array(arr))
        
        # Read densities (flattened array)
        total_points = n ** dim
        flat_densities = unpacker.unpack_farray(total_points, unpacker.unpack_float)
        
        # Reshape to multidimensional array
        if dim == 1:
            densities = np.array(flat_densities)
        elif dim == 2:
            densities = np.array(flat_densities).reshape((n, n))
        elif dim == 3:
            densities = np.array(flat_densities).reshape((n, n, n))
        
        return dim, n, R, tuple(grid), densities

def read_density_table_xdr(filename):
    with open(filename, 'rb') as f:
        # Read header
        dim = int(f.readline().decode().strip())
        n = int(f.readline().decode().strip())
        R = float(f.readline().decode().strip())
        
        # Read remaining data
        data = f.read()
        unpacker = xdrlib.Unpacker(data)
        
        # Read grid axes
        grid = []
        for _ in range(dim):
            # First read array length (XDR includes this)
            length = unpacker.unpack_uint()
            if length != n:
                raise ValueError(f"Expected axis length {n}, got {length}")
                
            # Read axis values
            axis = np.zeros(n, dtype=np.float32)
            for i in range(n):
                axis[i] = unpacker.unpack_float()
            grid.append(axis)
        
        # Read densities
        total_points = n ** dim
        length = unpacker.unpack_uint()
        if length != total_points:
            raise ValueError(f"Expected {total_points} density values, got {length}")
        
        flat_densities = np.zeros(total_points, dtype=np.float32)
        for i in range(total_points):
            flat_densities[i] = unpacker.unpack_float()
        
        # Reshape to multidimensional array
        if dim == 1:
            densities = flat_densities
        elif dim == 2:
            densities = flat_densities.reshape((n, n))
        elif dim == 3:
            densities = flat_densities.reshape((n, n, n))
        
        return dim, n, R, tuple(grid), densities


def compute_total_mass(D, N):
        V=D.dxbound*D.dybound*D.dzbound;
        dV = V*(1/N)**3  # Volume per particle
        dx=D.dxbound/N
        dy=D.dybound/N
        dz=D.dzbound/N
        
        # Coordinates for particle centers
        x = np.linspace(-D.dxbound/2 + dx/2, D.dxbound/2 - dx/2, N)
        y = np.linspace(-D.dybound/2 + dy/2, D.dybound/2 - dy/2, N)
        z = np.linspace(-D.dzbound/2 + dz/2, D.dzbound/2 - dz/2, N)
        
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
        R_vals = np.sqrt(X**2 + Y**2 + Z**2)
        
        # Vectorized density calculation
        rho_vals = D.getrho_vec(R_vals)
        
        total_mass = np.sum(rho_vals) * dV
        return total_mass

#dim, n, R, grid, densities = read_density_table_xdr("densitytable_xdr")
#=grid[0]
