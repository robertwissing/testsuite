#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 11:08:32 2024

@author: robertwi
"""

import numpy as np
import time
from SPH_methods import *
from grid_interpolate import *
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import Axes3D

import readtipsy as tip

from SPH_constants import *


#njobs=os.cpu_count()
#print(f"Number of available CPUs (logical cores): {njobs}")


def generate_particles_lattice(N, box_size, ):

    # Random positions in the box
    # particles[:, 1] = 1.0*np.random.rand(N) * box_size[0] - box_size[0] / 2.0  # Column 1 is 'x'
    # particles[:, 2] = 1.0*np.random.rand(N) * box_size[1] - box_size[1] / 2.0  # Column 2 is 'y'
    # particles[:, 3] = 1.0*np.random.rand(N) * box_size[2] - box_size[2] / 2.0  # Column 3 is 'z'
    
    ##FOR LATTICE
    # N = int(np.floor(N / 4) * 4)
    # n_cells = N // 4  # number of FCC unit cells
    # n_side = int(round(n_cells ** (1 / 3)))  # approximate side length in cells
    # correct_n_cells = n_side**3
    # if correct_n_cells != n_cells:
    #     # Update n_cells and N
    #     n_cells = correct_n_cells
    #     N = 4 * n_cells




    # Create a 2D array to store all particle properties
    # Columns: mass (0), x (1), y (2), z (3), vx (4), vy (5), vz (6), Bx (7), By (8), Bz (9)
    particles = np.zeros((N, PARTICLE_ARRAY_SIZE),dtype=PRECISION)  # 10 columns for the different particle properties
    
    # Set masses to be equal and sum to Mtot
    particles[:, 0] = Mtot / N  # Column 0 is 'mass'
    



    #FOR LATTICE
    # a = box_size[0] / n_side

    # positions = []
    # for i in range(n_side):
    #     for j in range(n_side):
    #         for k in range(n_side):
    #             # The corner of the current cell
    #             corner_x = i * a - box_size[0] / 2.0
    #             corner_y = j * a - box_size[1] / 2.0
    #             corner_z = k * a - box_size[2] / 2.0

    #             # In an FCC cell, there are 4 atoms:
    #             # corner
    #             positions.append([corner_x, corner_y, corner_z])
    #             # face-center (x,y)
    #             positions.append([corner_x + a/2, corner_y + a/2, corner_z])
    #             # face-center (y,z)
    #             positions.append([corner_x + a/2, corner_y, corner_z + a/2])
    #             # face-center (x,z)
    #             positions.append([corner_x, corner_y + a/2, corner_z + a/2])

    # # Convert to NumPy array (now exactly N rows, because 4 * n_side^3 = N)
    # positions = np.array(positions, dtype=np.float32)
    # # No slicing needed, because we matched the logic exactly.
    # # If there's a mismatch, check your n_side**3 vs. n_cells logic.

    # # Assign positions
    # particles[:, 1] = positions[:, 0]  # x
    # particles[:, 2] = positions[:, 1]  # y
    # particles[:, 3] = positions[:, 2]  # z

    
    
    
    
    
    
    
    mean = 0.0
    std_dev = box_size[0] / 6.0
    particles[:, 1] = np.random.normal(mean, std_dev, N)
    particles[:, 2] = np.random.normal(mean, std_dev, N)
    particles[:, 3] = np.random.normal(mean, std_dev, N)
    particles[:, 1] = np.clip(particles[:, 1], -box_size[0]/2.0, box_size[0]/2.0)
    particles[:, 2] = np.clip(particles[:, 2], -box_size[1]/2.0, box_size[1]/2.0)
    particles[:, 3] = np.clip(particles[:, 3], -box_size[2]/2.0, box_size[2]/2.0)

    
    # Generate random velocities
    particles[:, 4] = np.random.randn(N)  # Column 4 is 'vx'
    particles[:, 5] = np.random.randn(N)  # Column 5 is 'vy'
    particles[:, 6] = np.random.randn(N)  # Column 6 is 'vz'
    
    # Normalize velocity to achieve desired total kinetic energy
    initial_kinetic_energy = 0.5 * np.sum(particles[:, 0] * (particles[:, 4]**2 + particles[:, 5]**2 + particles[:, 6]**2))
    velocity_scaling = np.sqrt(E_kin / initial_kinetic_energy)
    particles[:, 4] *= velocity_scaling  # Scale 'vx'
    particles[:, 5] *= velocity_scaling  # Scale 'vy'
    particles[:, 6] *= velocity_scaling  # Scale 'vz'
    
    # Volume of the entire box and per particle (using the product of box_size dimensions)
    box_volume = np.prod(box_size)
    particle_volume = box_volume / N
    
    # Density is 7 column, we can put up an estimate
    particles[:, 7] = particles[:, 0]/particle_volume
    # thermal  is 8 column, leave at zero
    # soft is 9 column, leave at zero
    #particles[:, 9] = 0.0+100*np.abs(1.5+particles[:,1])
    particles[:, 8] = 100.0

    # metals is 10 column, leave at zero
    # potential is 11 column, leave at zero

        
    
    # Generate random magnetic fields
    particles[:, 12] = np.random.randn(N)  # Column 9 is 'Bx'
    particles[:, 13] = np.random.randn(N)  # Column 10 is 'By'
    particles[:, 14] = np.random.randn(N)  # Column 11 is 'Bz'
    

    
    # Normalize magnetic field to achieve desired total magnetic energy (considering volume)
    initial_magnetic_energy = 0.5 * np.sum(particle_volume * (particles[:, 12]**2 + particles[:, 13]**2 + particles[:, 14]**2))
    magnetic_scaling = np.sqrt(E_mag / initial_magnetic_energy)
    particles[:, 12] *= magnetic_scaling  # Scale 'Bx'
    particles[:, 13] *= magnetic_scaling  # Scale 'By'
    particles[:, 14] *= magnetic_scaling  # Scale 'Bz'
    
    # rest is spin and momenti  15 16 17 18 leave at zero  
    
    return particles


def generate_particles(N, box_size, Mtot=1.0, E_kin=1.0, E_mag=1.0):

    # Random positions in the box
    # particles[:, 1] = 1.0*np.random.rand(N) * box_size[0] - box_size[0] / 2.0  # Column 1 is 'x'
    # particles[:, 2] = 1.0*np.random.rand(N) * box_size[1] - box_size[1] / 2.0  # Column 2 is 'y'
    # particles[:, 3] = 1.0*np.random.rand(N) * box_size[2] - box_size[2] / 2.0  # Column 3 is 'z'
    
    ##FOR LATTICE
    # N = int(np.floor(N / 4) * 4)
    # n_cells = N // 4  # number of FCC unit cells
    # n_side = int(round(n_cells ** (1 / 3)))  # approximate side length in cells
    # correct_n_cells = n_side**3
    # if correct_n_cells != n_cells:
    #     # Update n_cells and N
    #     n_cells = correct_n_cells
    #     N = 4 * n_cells




    # Create a 2D array to store all particle properties
    # Columns: mass (0), x (1), y (2), z (3), vx (4), vy (5), vz (6), Bx (7), By (8), Bz (9)
    particles = np.zeros((N, PARTICLE_ARRAY_SIZE),dtype=PRECISION)  # 10 columns for the different particle properties
    
    # Set masses to be equal and sum to Mtot
    particles[:, 0] = Mtot / N  # Column 0 is 'mass'
    



    #FOR LATTICE
    # a = box_size[0] / n_side

    # positions = []
    # for i in range(n_side):
    #     for j in range(n_side):
    #         for k in range(n_side):
    #             # The corner of the current cell
    #             corner_x = i * a - box_size[0] / 2.0
    #             corner_y = j * a - box_size[1] / 2.0
    #             corner_z = k * a - box_size[2] / 2.0

    #             # In an FCC cell, there are 4 atoms:
    #             # corner
    #             positions.append([corner_x, corner_y, corner_z])
    #             # face-center (x,y)
    #             positions.append([corner_x + a/2, corner_y + a/2, corner_z])
    #             # face-center (y,z)
    #             positions.append([corner_x + a/2, corner_y, corner_z + a/2])
    #             # face-center (x,z)
    #             positions.append([corner_x, corner_y + a/2, corner_z + a/2])

    # # Convert to NumPy array (now exactly N rows, because 4 * n_side^3 = N)
    # positions = np.array(positions, dtype=np.float32)
    # # No slicing needed, because we matched the logic exactly.
    # # If there's a mismatch, check your n_side**3 vs. n_cells logic.

    # # Assign positions
    # particles[:, 1] = positions[:, 0]  # x
    # particles[:, 2] = positions[:, 1]  # y
    # particles[:, 3] = positions[:, 2]  # z

    
    
    
    
    
    
    
    mean = 0.0
    std_dev = box_size[0] / 6.0
    particles[:, 1] = np.random.normal(mean, std_dev, N)
    particles[:, 2] = np.random.normal(mean, std_dev, N)
    particles[:, 3] = np.random.normal(mean, std_dev, N)
    particles[:, 1] = np.clip(particles[:, 1], -box_size[0]/2.0, box_size[0]/2.0)
    particles[:, 2] = np.clip(particles[:, 2], -box_size[1]/2.0, box_size[1]/2.0)
    particles[:, 3] = np.clip(particles[:, 3], -box_size[2]/2.0, box_size[2]/2.0)

    
    # Generate random velocities
    particles[:, 4] = np.random.randn(N)  # Column 4 is 'vx'
    particles[:, 5] = np.random.randn(N)  # Column 5 is 'vy'
    particles[:, 6] = np.random.randn(N)  # Column 6 is 'vz'
    
    # Normalize velocity to achieve desired total kinetic energy
    initial_kinetic_energy = 0.5 * np.sum(particles[:, 0] * (particles[:, 4]**2 + particles[:, 5]**2 + particles[:, 6]**2))
    velocity_scaling = np.sqrt(E_kin / initial_kinetic_energy)
    particles[:, 4] *= velocity_scaling  # Scale 'vx'
    particles[:, 5] *= velocity_scaling  # Scale 'vy'
    particles[:, 6] *= velocity_scaling  # Scale 'vz'
    
    # Volume of the entire box and per particle (using the product of box_size dimensions)
    box_volume = np.prod(box_size)
    particle_volume = box_volume / N
    
    # Density is 7 column, we can put up an estimate
    particles[:, 7] = particles[:, 0]/particle_volume
    # softening length is 8 column, leave at zero
    # thermal energy is 9 column, leave at zero
    #particles[:, 9] = 0.0+100*np.abs(1.5+particles[:,1])
    particles[:, 8] = 100.0

    # metals is 10 column, leave at zero
    # potential is 11 column, leave at zero

        
    
    # Generate random magnetic fields
    particles[:, 12] = np.random.randn(N)  # Column 9 is 'Bx'
    particles[:, 13] = np.random.randn(N)  # Column 10 is 'By'
    particles[:, 14] = np.random.randn(N)  # Column 11 is 'Bz'
    

    
    # Normalize magnetic field to achieve desired total magnetic energy (considering volume)
    initial_magnetic_energy = 0.5 * np.sum(particle_volume * (particles[:, 12]**2 + particles[:, 13]**2 + particles[:, 14]**2))
    magnetic_scaling = np.sqrt(E_mag / initial_magnetic_energy)
    particles[:, 12] *= magnetic_scaling  # Scale 'Bx'
    particles[:, 13] *= magnetic_scaling  # Scale 'By'
    particles[:, 14] *= magnetic_scaling  # Scale 'Bz'
    
    # rest is spin and momenti  15 16 17 18 leave at zero  
    
    return particles

def calculate_angular_momentum_COM(particles):
    # Extract masses, positions, and velocities from the particle array
    masses = particles[:, 0]
    positions = particles[:, 1:4]
    velocities = particles[:, 4:7]

    # Calculate the total mass of the system
    total_mass = np.sum(masses)
    
    # Calculate the center of mass (COM)
    com = np.sum(masses[:, np.newaxis] * positions, axis=0) / total_mass
    
    # Calculate positions relative to the center of mass
    relative_positions = positions - com
    
    # Calculate angular momentum with respect to the COM
    angular_momentum_trans = np.sum(np.cross(relative_positions, masses[:, np.newaxis] * velocities), axis=0)
    #angular_momentum_trans = np.sum(np.cross(relative_positions, velocities), axis=0)
    
    
    angular_momentum_trans = np.nan_to_num(angular_momentum_trans) 
    return angular_momentum_trans, com

def calculate_conservative_properties(particles):
    N = len(particles)
    
    # Extract properties from particle array
    masses = particles[:, 0]
    positions = particles[:, 1:4]
    velocities = particles[:, 4:7]
    spin = particles[:, 15:18] * particles[:, 0][:, np.newaxis]
    #spin = particles[:, 15:18]
    momenti = particles[:, 18]
    magnetic_fields = particles[:, 12:15]
    thermal_energy = particles[:, 8]
    volume = particles[:, 0]/particles[:, 7]
    
    # Calculate Kinetic Energy
    translational_energy = 0.5 * np.sum(masses * np.sum(velocities**2, axis=1))
    
    
    # Calculate Rotational Energy, skipping cases where momenti == 0
    valid_momenti = momenti > 0  # Boolean mask for non-zero momenti
    rotational_energy = 0.5 * np.sum(np.sum(spin[valid_momenti]**2, axis=1) / momenti[valid_momenti])
    
    # Calculate Magnetic Energy (assuming uniform particle volume)
    magnetic_energy = 0.5 * np.sum(np.sum(volume[:, np.newaxis]*magnetic_fields**2, axis=1))  # ME = 0.5 * B^2 * V
    
    # Calculate Total Thermal Energy
    total_thermal_energy = np.sum(masses*thermal_energy)  # Sum of thermal energy in all particles
    
    # Calculate Linear Momentum (P = m * v)
    linear_momentum = np.sum(masses[:, np.newaxis] * velocities, axis=0)  # Sum of (m * v) for all particles
    
    # Calculate Angular Momentum (L = r x (m * v))
    angular_momentum_trans, com = calculate_angular_momentum_COM(particles)
    
    # Calculate Angular Momentum spin
    angular_momentum_spin = np.sum(spin, axis=0)
    
    angular_momentum = angular_momentum_trans+angular_momentum_spin
    
    # Total volume
    tot_volume = np.sum(volume)
    
    kinetic_energy=translational_energy + rotational_energy

    
    # Return all calculated quantities
    return {
        'kinetic_energy': kinetic_energy,
        'translational_energy': translational_energy,
        'rotational_energy': rotational_energy,
        'magnetic_energy': magnetic_energy,
        'thermal_energy': total_thermal_energy,
        'linear_momentum': linear_momentum,
        'angular_momentum_trans_rat': np.nan_to_num(np.abs(angular_momentum_trans)/np.abs(angular_momentum)),
        'angular_momentum_spin_rat': np.nan_to_num(np.abs(angular_momentum_spin)/np.abs(angular_momentum)),
        'trans_angular_momentum': angular_momentum_trans,
        'spin_angular_momentum': angular_momentum_spin,
        'total_angular_momentum': angular_momentum,
        'total_angular_momentum/N': angular_momentum/N,
        'com': com,
        'volume': tot_volume
    }


def initial_setup(particles,box_size,nsmooth,periodic_box,h_init=None,gradient=False):
       particles_initial = particles.copy()  # or np.copy(particles)
       particles = remove_exact_duplicates(particles_initial)
       print("Removed exact position duplicates(why are they here?):", len(particles_initial)-len(particles))

       N=len(particles[:,0])
       if(nsmooth > N):
           nsmooth=N
       #Step 0: if we do not have smoothing length, do a very rough estimate of it
       if h_init is None:
           h_init=estimate_smoothing_length(particles[:,1:4], box_size=box_size, nSmooth=nsmooth);
       else:
           h_init=np.max(h_init)
       # Step 1: Time building the tree
       start_time = time.time()
       tree = build_tree(particles[:,1:4],box_size=box_size,h=2*h_init,periodic_box=periodic_box)
       tree_build_time = time.time() - start_time
       print(f"Time to build kd tree: {tree_build_time:.4f} seconds")

       # Step 2: Time calculating the sif unmerged_particles.size > 0 and unmerged_particles.ndim > 1:printmoothing length and neighbor list (Staged Querying)
       start_time = time.time()
       h_smooth, neighbors_list, distances, vec = calculate_smoothing_length_staged(
            tree, particles[:,1:4],box_size=box_size, periodic_box=periodic_box, final_k=nsmooth)
       smoothing_length_time_staged = time.time() - start_time
       print(f"Time to calculate smoothing lengths and neighbors (Staged) with kd tree: {smoothing_length_time_staged:.4f} seconds")


       start_time = time.time()
       density = calculate_density(particles[:,1:4], particles[:,0], h_smooth, neighbors_list, distances, box_size=box_size, periodic_box=periodic_box)
       density_calculation_time = time.time() - start_time
       print(f"Time to calculate density with kd tree: {density_calculation_time:.4f} seconds")
       particles[:,7]=density
       particles[:,20]=h_smooth
       
       if gradient==True:
           
           particles = EOS(particles)
           
           start_time = time.time()
           Q0,Q1x,Q1y,Q1z,E0x,E0y,E0z,dvxdt,dvydt,dvzdt,dudt = calculate_gradients(particles, neighbors_list, distances,vec, box_size=box_size, periodic_box=periodic_box)
           gradient_calculation_time = time.time() - start_time
           print(f"Time to calculate gradients with kd tree: {gradient_calculation_time:.4f} seconds")

           # Calculate mean values
           mean_Q0 = np.mean(Q0)
           mean_Q1x = np.mean(Q1x)
           mean_Q1y = np.mean(Q1y)
           mean_Q1z = np.mean(Q1z)
           mean_E0x = np.mean(E0x)
           mean_E0y = np.mean(E0y)
           mean_E0z = np.mean(E0z)
           mean_dvxdt = np.mean(dvxdt)
           mean_dvydt = np.mean(dvydt)
           mean_dvzdt = np.mean(dvzdt)
           
           print(dvxdt)
           
           # Print only the means
           print(f"Mean of Q0: {mean_Q0:.4f}")
           print(f"Mean of Q1x: {mean_Q1x:.4f}")
           print(f"Mean of Q1y: {mean_Q1y:.4f}")
           print(f"Mean of Q1z: {mean_Q1z:.4f}")
           print(f"Mean of E0x: {mean_E0x:.4f}")
           print(f"Mean of E0y: {mean_E0y:.4f}")
           print(f"Mean of E0z: {mean_E0z:.4f}")
           print(f"Mean of dvxdt: {mean_dvxdt:.4f}")
           print(f"Mean of dvydt: {mean_dvydt:.4f}")
           print(f"Mean of dvzdt: {mean_dvzdt:.4f}")
           #print("distances!!",distances,"nsmooth len distances",nsmooth,len(distances))
          
       properties = calculate_conservative_properties(particles)
       print("Translational Energy:", properties['translational_energy'])
       print("Rotational Energy:", properties['rotational_energy'])
       print("Kinetic Energy:", properties['kinetic_energy'])
       print("Magnetic Energy:", properties['magnetic_energy'])
       print("Thermal Energy:", properties['thermal_energy'])
       print("Linear Momentum:", properties['linear_momentum'])
       print("Angular Momentum:", properties['total_angular_momentum'])
       print("Trans Angular Momentum:", properties['trans_angular_momentum'])
       print("Spin Angular Momentum:", properties['spin_angular_momentum'])
       print("Angular Momentum/N:", properties['total_angular_momentum/N'])
       print("Angular Momentum Translational/ratio:", properties['angular_momentum_trans_rat'])
       print("Angular Momentum Error/Spin/ratio:", properties['angular_momentum_trans_rat'])
       print("Centre of mass", properties['com'])
       print("Total Volume:", properties['volume'])
       
       return tree,neighbors_list,distances,vec,particles


def merge_particles_full(particles, tree, neighbors_list, distances, vec, periodic_box, box_size, condition=None, merge_rounds=1):
    
    h=particles[:,20]
    nsmooth=len(neighbors_list[0,:])
    merged_particles=particles
    unmerged_particles=None
    # Loop over the number of merge rounds
    
    for i in range(merge_rounds):
        print(f"Starting merge round {i + 1}")
        
        merge_condition=map_all_from_condition(particles,condition)
        
        
        # Step 4: Time merging particles
        start_time = time.time()
        merged_particles, unmerged_particles = merge_particles(
            particles, tree, neighbors_list, distances, vec, box_size=box_size, periodic_box=periodic_box,condition=condition)
        
        if merged_particles.size == 0:  # Checks if merged_particles is empty
            merged_particles = unmerged_particles
            unmerged_particles = []
            break
        if len(unmerged_particles) > 0:
            print("Final number of unmerged particles ",len(unmerged_particles))
            merged_particles = addunmerged(merged_particles, unmerged_particles)
        
        
        
        merge_time = time.time() - start_time
        print(f"Time to merge particles with kd tree in round {i + 1}: {merge_time:.4f} seconds")
        
        # Update tree and other parameters for the next merge round
        if i < merge_rounds - 1:
            tree_merged, neighbors_list_merged, distances_merged, vec_merged, merged_particles = initial_setup(
                merged_particles, box_size=box_size, nsmooth=nsmooth, h_init=h, periodic_box=periodic_box)
            # Use the merged data as the input for the next round
            particles, tree, neighbors_list, distances, vec = merged_particles, tree_merged, neighbors_list_merged, distances_merged, vec_merged

    # Return final merged and unmerged particles after all rounds
    return merged_particles, unmerged_particles

def check_for_nan_values(array, label):
    nan_mask = np.isnan(array)
    if nan_mask.any():
        print(f"[DEBUG] NaN values found in {label}.")
        nan_indices = np.where(nan_mask.any(axis=1))[0]  # Find indices of rows with NaNs
        print(f"NaN found at rows (showing up to 10): {nan_indices[:10]}")
        print(f"Sample rows with NaN values in {label}:\n", array[nan_indices[:10]])
    else:
        print(f"[DEBUG] No NaN values in {label}.")


def split_particles_full(particles, tree, neighbors_list, distances, vec, periodic_box, box_size, condition=None, split_rounds=1):
    # Loop over the number of merge rounds
    split_particles=particles
    nsmooth=len(neighbors_list[0,:])
    kmove=0
    ktar=1
    for i in range(split_rounds):
        print(f"Starting split round {i + 1}")
        h=split_particles[:,20]
        
        split_target=map_all_from_condition(particles,condition)
        
        if particles[split_target].size == 0:
            print("Condition hit for all particles! !! !")
            break
        
        #Step check for nan values in particles
        check_for_nan_values(particles, "particles")  
        # Step 4: Time merging particles
        print("length of particles to be split",split_target.size)
        print("number of particles:",len(particles[:,0]))
        start_time = time.time()
        split_particles = split_particless(
            particles, neighbors_list,distances,vec,box_size=box_size, periodic_box=periodic_box,method="nn_ort",kin=ktar,split_target=split_target)
        split_time = time.time() - start_time
        print(f"Time to split particles with kd tree in round {i + 1}: {split_time:.4f} seconds")
        if ktar < nsmooth-kmove:
            ktar += kmove
            #kmove = ktar

        check_for_nan_values(split_particles, "split_particles")        
        print("split based on neighbour",ktar)
        # Update tree and other parameters for the next merge round
        if i < split_rounds - 1:
            tree_split, neighbors_list_split, distances_split, vec_split,  split_particles = initial_setup(
                split_particles, box_size=box_size, nsmooth=nsmooth, h_init=h, periodic_box=periodic_box
            )
            # Use the split data as the input for the next round
            particles, tree, neighbors_list, distances, vec = split_particles, tree_split, neighbors_list_split, distances_split, vec_split
    # Return final split particles after all rounds
    return split_particles

def mergeandsplit(particles, tree, neighbors_list, distances, vec, periodic_box, box_size, nsmooth,scondition=None ,mcondition=None, mergesplit_rounds=1):
    for _ in range(mergesplit_rounds):
        # Perform merging
        particles, unmerged_particles = merge_particles_full(
            particles, tree, neighbors_list, distances, vec, periodic_box=periodic_box, box_size=box_size, condition=mcondition,merge_rounds=1
        )
        
        # Rebuild tree and neighbor lists after merging
        tree, neighbors_list, distances, vec, particles = initial_setup(
            particles, box_size=box_size, nsmooth=nsmooth, periodic_box=periodic_box, h_init=particles[:, 20]
        )
        
        # Perform splitting
        particles = split_particles_full(
            particles, tree, neighbors_list, distances, vec, periodic_box=periodic_box, box_size=box_size, condition=scondition, split_rounds=1
        )
        
        # Rebuild tree and neighbor lists after splitting
        tree, neighbors_list, distances, vec, particles = initial_setup(
            particles, box_size=box_size, nsmooth=nsmooth, periodic_box=periodic_box, h_init=particles[:, 20]
        )

    return particles, tree, neighbors_list, distances, vec
