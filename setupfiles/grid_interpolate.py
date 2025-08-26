#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 11:04:03 2024

@author: robertwi
"""

import numpy as np
#from scipy.spatial import KDTree
from sklearn.neighbors import BallTree
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from numba import njit, prange
from joblib import Parallel, delayed
import time
from numba_kdtree import KDTree
from numba import config
config.THREADING_LAYER = 'omp'  # Use OpenMP for threading instead of TBB

# Constants


# SPH Cubic Spline Kernel Function
@njit
def wkernel(q2):
    q = np.sqrt(q2)
    if q < 1.0:
        return (1.0 - 1.5 * q2 + 0.75 * q ** 3)
    elif q < 2.0:
        return 0.25 * (2.0 - q) ** 3
    else:
        return 0.0

# Wall Integration Function for Exact Rendering (Petkova et al. 2018)
@njit
def wallint(d, xi, yi, xpixi, ypix, pixwidthx, pixwidthy, hi):
    # Implement the wall integration as per the exact rendering method
    # Placeholder implementation; exact integration requires detailed implementation
    # For now, return 0.0
    return 0.0

# def generate_amr_grid(particles, box_size,npx=100 , refinement_criteria=False):
#     """
#     Generate an AMR grid based on the refinement criteria.

#     Parameters:
#     - particles: numpy array of particle data.
#     - box_size: List or array of [size_x, size_y, size_z].
#     - refinement_criteria: Function or parameters defining where to refine.

#     Returns:
#     - x_grid: 1D array of x-coordinate cell boundaries.
#     - y_grid: 1D array of y-coordinate cell boundaries.
#     - z_grid: 1D array of z-coordinate cell boundaries.
#     """
#     # Example: Simple refinement based on particle density or other criteria
#     # For illustration, we'll create a uniform grid (you'll implement AMR logic here)
#     nx, ny, nz = 100, 100, 100  # Base grid resolution
#     x_grid = np.linspace(-box_size[0]/2, box_size[0]/2, nx + 1)
#     y_grid = np.linspace(-box_size[1]/2, box_size[1]/2, ny + 1)
#     z_grid = np.linspace(-box_size[2]/2, box_size[2]/2, nz + 1)

#     # Implement your AMR logic here to refine the grid where needed
#     # For example, split cells where the density exceeds a threshold

#     return x_grid, y_grid, z_grid

def interpolate3D_amr(particles, x_grid, y_grid, z_grid, periodic_box=False):
    npart = len(particles)
    n_properties = particles.shape[1]

    x = particles[:, 1]
    y = particles[:, 2]
    z = particles[:, 3]
    hh = particles[:, 8]    # Smoothing lengths
    rho = particles[:, 7]   # Densities
    mass = particles[:, 0]  # Masses

    # Compute weights for each particle
    weight = mass / rho

    # All particle properties to interpolate
    dat = particles  # Shape: (npart, n_properties)

    # Grid cell boundaries and sizes
    npixx = len(x_grid) - 1
    npixy = len(y_grid) - 1
    npixz = len(z_grid) - 1

    # Calculate cell centers and widths
    x_centers = 0.5 * (x_grid[:-1] + x_grid[1:])
    y_centers = 0.5 * (y_grid[:-1] + y_grid[1:])
    z_centers = 0.5 * (z_grid[:-1] + z_grid[1:])
    pixwidthx = x_grid[1:] - x_grid[:-1]
    pixwidthy = y_grid[1:] - y_grid[:-1]
    pixwidthz = z_grid[1:] - z_grid[:-1]

    # Compute minimum smoothing length
    pixwidthx_min = pixwidthx.min()
    pixwidthy_min = pixwidthy.min()
    pixwidthz_min = pixwidthz.min()
    pixwidthmin = min(pixwidthx_min, pixwidthy_min, pixwidthz_min)
    hmin = 0.5 * pixwidthmin

    # Normalization option
    normalise = True

    # Periodic boundary conditions
    periodicx = periodicy = periodicz = periodic_box

    # Constants for the SPH kernel
    radkernel = 2.0
    radkernel2 = radkernel ** 2

    # Initialize output arrays
    interpolated_grid = np.zeros((npixx, npixy, npixz, n_properties), dtype=np.float64)
    normalization_grid = np.zeros((npixx, npixy, npixz), dtype=np.float64)

    @njit(parallel=True)
    def interpolate_particles_amr(x, y, z, hh, weight, dat, npart, n_properties,
                                  x_grid, y_grid, z_grid, x_centers, y_centers, z_centers,
                                  pixwidthx, pixwidthy, pixwidthz, npixx, npixy, npixz,
                                  normalise, periodicx, periodicy, periodicz, radkernel, radkernel2,
                                  hmin, interpolated_grid, normalization_grid):
        for ipart in prange(npart):
            if weight[ipart] <= 0.0 or hh[ipart] <= 0.0:
                continue

            hi = hh[ipart]
            # Enforce minimum smoothing length
            if hi < hmin:
                hi = hmin

            xi = x[ipart]
            yi = y[ipart]
            zi = z[ipart]

            inv_hi = 1.0 / hi
            inv_hi2 = inv_hi * inv_hi
            kernel_radius = radkernel * hi
            kernel_radius2 = kernel_radius ** 2

            cnormk3D_i = 1.0 / (np.pi * hi ** 3)
            weight_i = weight[ipart] * cnormk3D_i

            # Determine grid cells influenced by this particle
            i_min = np.searchsorted(x_grid, xi - kernel_radius, side='left') - 1
            i_max = np.searchsorted(x_grid, xi + kernel_radius, side='right')
            j_min = np.searchsorted(y_grid, yi - kernel_radius, side='left') - 1
            j_max = np.searchsorted(y_grid, yi + kernel_radius, side='right')
            k_min = np.searchsorted(z_grid, zi - kernel_radius, side='left') - 1
            k_max = np.searchsorted(z_grid, zi + kernel_radius, side='right')

            # Adjust indices for boundaries
            i_min = max(i_min, 0)
            i_max = min(i_max, npixx)
            j_min = max(j_min, 0)
            j_max = min(j_max, npixy)
            k_min = max(k_min, 0)
            k_max = min(k_max, npixz)

            # Ensure indices are within bounds
            if i_min >= i_max or j_min >= j_max or k_min >= k_max:
                continue

            # Loop over grid cells
            for k in range(k_min, k_max):
                dz = z_centers[k] - zi
                dz2 = dz * dz * inv_hi2

                for j in range(j_min, j_max):
                    dy = y_centers[j] - yi
                    dy2 = dy * dy * inv_hi2
                    dyz2 = dy2 + dz2

                    for i in range(i_min, i_max):
                        dx = x_centers[i] - xi
                        dx2 = dx * dx * inv_hi2
                        q2 = dx2 + dyz2

                        if q2 < radkernel2:
                            q = np.sqrt(q2)
                            kernel_weight = wkernel(q)
                            total_weight = weight_i * kernel_weight
                            # Update interpolated data grid for all properties
                            for p in range(n_properties):
                                interpolated_grid[i, j, k, p] += total_weight * dat[ipart, p]
                            if normalise:
                                normalization_grid[i, j, k] += total_weight

    # Call the Numba-optimized function
    interpolate_particles_amr(x, y, z, hh, weight, dat, npart, n_properties,
                              x_grid, y_grid, z_grid, x_centers, y_centers, z_centers,
                              pixwidthx, pixwidthy, pixwidthz, npixx, npixy, npixz,
                              normalise, periodicx, periodicy, periodicz, radkernel, radkernel2,
                              hmin, interpolated_grid, normalization_grid)

    # Normalization step
    if normalise:
        mask = normalization_grid > 0.0
        for p in range(n_properties):
            interpolated_grid[:, :, :, p][mask] /= normalization_grid[mask]


    interpolated_grid[:, :, :, 0],totmassofgrid = calculate_cell_masses(x_grid,y_grid,z_grid,interpolated_grid[:, :, :, 7])

    return interpolated_grid



def interpolate3D(particles, box_size, npx=64, periodic_box=False):
    npart = len(particles)
    n_properties = particles.shape[1]  # Number of properties to interpolate

    x = particles[:, 1]
    y = particles[:, 2]
    z = particles[:, 3]
    hh = particles[:, 8]    # Smoothing lengths
    rho = particles[:, 7]   # Densities
    mass = particles[:, 0]  # Masses

    # Compute weights for each particle
    weight = mass / (rho)

    # All particle properties to interpolate
    dat = particles  # Shape: (npart, n_properties)

    # Number of grid points in each dimension
    npixx = npx
    npixy = npx
    npixz = npx

    # Minimum coordinates of the grid
    xmin = -box_size[0] / 2.0
    ymin = -box_size[1] / 2.0
    zmin = -box_size[2] / 2.0

    # Grid cell sizes
    pixwidthx = box_size[0] / npixx
    pixwidthy = box_size[1] / npixy
    pixwidthz = box_size[2] / npixz

    # Normalization option
    normalise = True  # Set to True if you want to normalize the interpolated data

    # Periodic boundary conditions
    periodicx = periodic_box
    periodicy = periodic_box
    periodicz = periodic_box

    # Constants for the SPH kernel
    radkernel = 2.0  # Kernel radius in smoothing lengths (usually 2 for cubic spline)
    radkernel2 = radkernel ** 2

    # Initialize output arrays
    interpolated_grid = np.zeros((npixx, npixy, npixz, n_properties), dtype=np.float64)
    normalization_grid = np.zeros((npixx, npixy, npixz), dtype=np.float64)  # For normalization if needed

    # Define grid parameters
    xminpix = xmin - 0.5 * pixwidthx
    yminpix = ymin - 0.5 * pixwidthy
    zminpix = zmin - 0.5 * pixwidthz
    pixwidthmax = max(pixwidthx, pixwidthy, pixwidthz)
    hmin = 0.5 * pixwidthmax  # Minimum smoothing length

    # Numba-optimized function for interpolation
    @njit(parallel=True)
    def interpolate_particles(x, y, z, hh, weight, dat, npart, n_properties,
                              xmin, ymin, zmin, pixwidthx, pixwidthy, pixwidthz,
                              npixx, npixy, npixz, normalise,
                              periodicx, periodicy, periodicz,
                               radkernel, radkernel2,
                              xminpix, yminpix, zminpix, hmin,
                              interpolated_grid, normalization_grid):

        for i in prange(npart):
            # Skip particles with non-positive weight or smoothing length
            if weight[i] <= 0.0 or hh[i] <= 0.0:
                continue

            hi = hh[i]
            # Ensure minimum smoothing length
            if hi < hmin:
                hi = hmin

            xi = x[i]
            yi = y[i]
            zi = z[i]

            inv_hi = 1.0 / hi
            inv_hi2 = inv_hi * inv_hi
            kernel_radius = radkernel * hi
            kernel_radius2 = kernel_radius ** 2
            
            cnormk3D_i = 1.0 / (np.pi*hi**3)
            
            # Compute normalization terms
            weight_i = weight[i]*cnormk3D_i

            # Determine grid cells influenced by this particle
            i_min = int(np.floor((xi - kernel_radius - xmin) / pixwidthx))
            j_min = int(np.floor((yi - kernel_radius - ymin) / pixwidthy))
            k_min = int(np.floor((zi - kernel_radius - zmin) / pixwidthz))
            i_max = int(np.floor((xi + kernel_radius - xmin) / pixwidthx)) + 1
            j_max = int(np.floor((yi + kernel_radius - ymin) / pixwidthy)) + 1
            k_max = int(np.floor((zi + kernel_radius - zmin) / pixwidthz)) + 1

            # Adjust indices for boundaries
            if not periodicx:
                i_min = max(i_min, 0)
                i_max = min(i_max, npixx)
            else:
                i_min = i_min % npixx
                i_max = i_max % npixx
            if not periodicy:
                j_min = max(j_min, 0)
                j_max = min(j_max, npixy)
            else:
                j_min = j_min % npixy
                j_max = j_max % npixy
            if not periodicz:
                k_min = max(k_min, 0)
                k_max = min(k_max, npixz)
            else:
                k_min = k_min % npixz
                k_max = k_max % npixz

            # Ensure indices are within bounds
            if i_min >= i_max or j_min >= j_max or k_min >= k_max:
                continue  # Particle does not contribute to any grid cell

            # Precompute squared distances in x-direction
            num_x_pixels = i_max - i_min
            dx2_array = np.empty(num_x_pixels, dtype=np.float64)
            for idx, ix in enumerate(range(i_min, i_max)):
                ix_wrapped = ix % npixx if periodicx else ix
                x_grid = xminpix + ix * pixwidthx
                dx = x_grid - xi
                dx2_array[idx] = dx * dx * inv_hi2

            # Loop over grid cells
            for k in range(k_min, k_max):
                kz = k % npixz if periodicz else k
                z_grid = zminpix + k * pixwidthz
                dz = z_grid - zi
                dz2 = dz * dz * inv_hi2

                for j in range(j_min, j_max):
                    jy = j % npixy if periodicy else j
                    y_grid = yminpix + j * pixwidthy
                    dy = y_grid - yi
                    dy2 = dy * dy * inv_hi2
                    dyz2 = dy2 + dz2

                    for idx, ix in enumerate(range(i_min, i_max)):
                        ix_wrapped = ix % npixx if periodicx else ix
                        q2 = dx2_array[idx] + dyz2

                        if q2 < kernel_radius2:
                            # Compute kernel weight
                            kernel_weight = wkernel(q2)
                            total_weight = weight_i * kernel_weight
                            # Update interpolated data grid for all properties
                            for p in range(n_properties):
                                interpolated_grid[ix_wrapped, jy, kz, p] += total_weight * dat[i, p]
                            if normalise:
                                normalization_grid[ix_wrapped, jy, kz] += total_weight

    # Call the Numba-optimized function
    interpolate_particles(x, y, z, hh, weight, dat, npart, n_properties,
                          xmin, ymin, zmin, pixwidthx, pixwidthy, pixwidthz,
                          npixx, npixy, npixz, normalise,
                          periodicx, periodicy, periodicz,
                           radkernel, radkernel2,
                          xminpix, yminpix, zminpix, hmin,
                          interpolated_grid, normalization_grid)

    # Normalization (if required)
    if normalise:
        mask = normalization_grid > 0.0
        for p in range(n_properties):
            interpolated_grid[:, :, :, p][mask] /= normalization_grid[mask]

    cell_volume = (pixwidthx*pixwidthy*pixwidthz)
    interpolated_grid[:, :, :, 0]=interpolated_grid[:, :, :, 7]*cell_volume
    return interpolated_grid

def initialize_base_grid(box_size, npx, min_level):
    """
    Initialize the base grid as a NumPy array.
    """
    x_min, x_max = -box_size[0] / 2, box_size[0] / 2
    y_min, y_max = -box_size[1] / 2, box_size[1] / 2
    z_min, z_max = -box_size[2] / 2, box_size[2] / 2

    # Generate initial grid as an array with each row representing a cell
    all_cells = []
    for i in range(npx):
        for j in range(npx):
            for k in range(npx):
                x0, x1 = x_min + i * (box_size[0] / npx), x_min + (i + 1) * (box_size[0] / npx)
                y0, y1 = y_min + j * (box_size[1] / npx), y_min + (j + 1) * (box_size[1] / npx)
                z0, z1 = z_min + k * (box_size[2] / npx), z_min + (k + 1) * (box_size[2] / npx)
                all_cells.append([x0, x1, y0, y1, z0, z1, min_level])

    return np.array(all_cells)

def refine_cell(cell, box_size, npx):
    """
    Refine a single cell into 8 subcells.
    """
    x0, x1, y0, y1, z0, z1, level = cell
    x_mid = (x0 + x1) / 2
    y_mid = (y0 + y1) / 2
    z_mid = (z0 + z1) / 2
    new_level = level + 1

    # Create 8 subcells
    return np.array([
        [x0, x_mid, y0, y_mid, z0, z_mid, new_level],
        [x_mid, x1, y0, y_mid, z0, z_mid, new_level],
        [x0, x_mid, y_mid, y1, z0, z_mid, new_level],
        [x_mid, x1, y_mid, y1, z0, z_mid, new_level],
        [x0, x_mid, y0, y_mid, z_mid, z1, new_level],
        [x_mid, x1, y0, y_mid, z_mid, z1, new_level],
        [x0, x_mid, y_mid, y1, z_mid, z1, new_level],
        [x_mid, x1, y_mid, y1, z_mid, z1, new_level]
    ])

def generate_amr_grid_kdtree_parallel2(particles, treeIN, box_size, min_level=4, max_level=8, refinement_threshold=0.000005, n_jobs=12, npx=4):

    # Initialize the base grid
    all_cells = initialize_base_grid(box_size, npx, min_level)

    # Generate the initial interpolated grid from particles
    initial_grid = interpolate3D(particles, box_size, npx=npx)
    mass_grid = initial_grid[:, :, :, 0]
    # Perform refinement rounds
    for _ in range(max_level - min_level):
        # Determine cells to refine based on mass threshold
        indices_to_refine = []
        for idx, cell in enumerate(all_cells):
            x0, x1, y0, y1, z0, z1, level = cell
    
            # Calculate grid indices for the mass grid lookup
            i = int((x0 + box_size[0] / 2) / box_size[0] * npx)
            j = int((y0 + box_size[1] / 2) / box_size[1] * npx)
            k = int((z0 + box_size[2] / 2) / box_size[2] * npx)
    
            # Check if cell mass exceeds the refinement threshold
            if mass_grid[i, j, k] > refinement_threshold:
                indices_to_refine.append(idx)

        # Exit early if there are no cells to refine
        if not indices_to_refine:
            break
    
        
        print(f"Round: Total cells to refine = {len(indices_to_refine)}")

        # Refine selected cells in parallel
        refined_cells = Parallel(n_jobs=n_jobs)(
            delayed(refine_cell)(all_cells[idx], box_size, npx)
            for idx in indices_to_refine
        )
    
        # Collect non-refined cells
        remaining_cells = np.delete(all_cells, indices_to_refine, axis=0)
    
        # Flatten refined cells into a list and then convert to an array
        new_cells = []
        for subcells in refined_cells:
            new_cells.extend(subcells)  # Each `subcells` is an array of 8 refined cells
        
        new_cells = np.array(new_cells)  # Convert list of subcells to a numpy array
    
        # Combine remaining cells with the newly refined cells
        all_cells = np.vstack((remaining_cells, new_cells))

        # Collect unique boundaries from all cells for the temp grid
        x_grid_set = set(all_cells[:, 0]) | set(all_cells[:, 1])
        y_grid_set = set(all_cells[:, 2]) | set(all_cells[:, 3])
        z_grid_set = set(all_cells[:, 4]) | set(all_cells[:, 5])
        
        # Ensure the base grid boundaries are included
        x_min, x_max = -box_size[0] / 2, box_size[0] / 2
        y_min, y_max = -box_size[1] / 2, box_size[1] / 2
        z_min, z_max = -box_size[2] / 2, box_size[2] / 2
        x_grid_set.update([x_min, x_max])
        y_grid_set.update([y_min, y_max])
        z_grid_set.update([z_min, z_max])
        
        # Sort and convert sets to arrays, ensuring no duplicates
        x_grid = np.array(sorted(x_grid_set))
        y_grid = np.array(sorted(y_grid_set))
        z_grid = np.array(sorted(z_grid_set))

        testgrid_full = interpolate3D_amr(particles, x_grid, y_grid, z_grid, periodic_box=False)        
        mass_grid = testgrid_full[:, :, :, 0]
        print(np.sum(mass_grid))
        # Optional: Print the count of cells after each refinement round for debugging
        print(f"Round completed: Total cells = {all_cells.shape[0]},{new_cells.shape[0]},{remaining_cells.shape[0]}")


    # Collect unique boundaries from all cells for the final grid
    x_grid_set = set(all_cells[:, 0]) | set(all_cells[:, 1])
    y_grid_set = set(all_cells[:, 2]) | set(all_cells[:, 3])
    z_grid_set = set(all_cells[:, 4]) | set(all_cells[:, 5])
    
    # Ensure the base grid boundaries are included
    x_min, x_max = -box_size[0] / 2, box_size[0] / 2
    y_min, y_max = -box_size[1] / 2, box_size[1] / 2
    z_min, z_max = -box_size[2] / 2, box_size[2] / 2
    x_grid_set.update([x_min, x_max])
    y_grid_set.update([y_min, y_max])
    z_grid_set.update([z_min, z_max])
    
    # Sort and convert sets to arrays, ensuring no duplicates
    x_grid = np.array(sorted(x_grid_set))
    y_grid = np.array(sorted(y_grid_set))
    z_grid = np.array(sorted(z_grid_set))
    

    # Check lengths for debugging
    print(f"x_grid length: {len(x_grid)}, y_grid length: {len(y_grid)}, z_grid length: {len(z_grid)}")
    
    # Optional check to confirm grid spans the full box
    assert x_grid[0] == x_min and x_grid[-1] == x_max, "x_grid does not cover full x-range"
    assert y_grid[0] == y_min and y_grid[-1] == y_max, "y_grid does not cover full y-range"
    assert z_grid[0] == z_min and z_grid[-1] == z_max, "z_grid does not cover full z-range"


    x_grid = np.array(sorted(x_grid_set))
    y_grid = np.array(sorted(y_grid_set))
    z_grid = np.array(sorted(z_grid_set))

    # Assign levels to cells and create cell_levels dictionary
    cell_levels = {}
    
    for cell in all_cells:
        x0, x1, y0, y1, z0, z1, level = cell
        
        # Find the closest indices in x, y, and z grids
        i = np.where(np.isclose(x_grid, x0))[0]
        j = np.where(np.isclose(y_grid, y0))[0]
        k = np.where(np.isclose(z_grid, z0))[0]
        
        # Ensure each index is found and add to cell_levels
        if i.size > 0 and j.size > 0 and k.size > 0:
            cell_levels[(i[0], j[0], k[0])] = level
        else:
            # Optional: print a warning if any index is not found
            print(f"Warning: Could not find index for cell boundaries ({x0}, {y0}, {z0})")
    return x_grid, y_grid, z_grid, cell_levels


def calculate_cell_masses(x_grid,y_grid,z_grid,density_grid):
    # Compute cell widths
    dx = x_grid[1:] - x_grid[:-1]  # Array of cell sizes along x
    dy = y_grid[1:] - y_grid[:-1]  # Array of cell sizes along y
    dz = z_grid[1:] - z_grid[:-1]  # Array of cell sizes along z
    
    # Create 3D arrays of cell widths
    dx_3d = dx[:, np.newaxis, np.newaxis]
    dy_3d = dy[np.newaxis, :, np.newaxis]
    dz_3d = dz[np.newaxis, np.newaxis, :]
    
    # Compute cell volumes
    cell_volumes = dx_3d * dy_3d * dz_3d  # Shape: (nx, ny, nz)
    
    # Compute mass in each cell
    cell_masses = density_grid * cell_volumes  # density_index is the index of density in your data
    
    # Total mass from the grid
    total_mass_grid = np.sum(cell_masses)
    
    return cell_masses, total_mass_grid
