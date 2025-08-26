#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 16:10:30 2024

@author: robertwi
"""

import numpy as np
import time
from SPH_methods import *
from grid_render import *
from grid_interpolate import *
from SPH_general import *
from SPH_render import *
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import Axes3D

from SPH_constants import *

import readtipsy as tip

def compute_timestep_criterion(particles, factor_acc=0.3, factor_cfl=0.4, factor_u=0.25):
    """
    Compute an SPH-like time step criterion based on:
      1) CFL condition: dt_CFL ~ (h / (c + alpha * c_sig)) 
      2) Force condition: dt_force ~ sqrt(h / |a|)
    Then take the minimum over all particles and multiply by a safety factor.

    Parameters
    ----------
    particles   : ndarray
        The array of particles (positions, velocities, densities, sound speeds, etc.)
    cfl_factor  : float
        Safety factor (like 0.3 or 0.4).
    alpha       : float
        Artificial viscosity factor (typical range 0.5–1.0+).

    Returns
    -------
    dt_crit : float
        The recommended maximum timestep to maintain stability.
    """
    # Example placeholders:
    h = particles[:, 20]          # smoothing length
    ax = particles[:, 12]          # acceleration x
    ay = particles[:, 13]          # acceleration y
    az = particles[:, 14]          # acceleration z
    u = particles[:, 8]          # u
    dudt = particles[:, 15]          # dudt
    a_magnitude = np.sqrt(ax**2 + ay**2 + az**2) + 1e-20  # avoid division by zero

    Poverrho=particles[:, 11]/particles[:, 7] 
    # If you store sound speed in an array or compute from EOS:
    c_sound = np.sqrt((5./3.)*Poverrho);

    # Some approximation of a signal speed for artificial viscosity 
    # (often c_sig ~ c_sound + k * |v|, or some function of local velocities).
    # For simplicity, let's just set c_sig = c_sound here.
    c_sig = c_sound

    # 1) dt from CFL condition
    dt_cfl_array = h*0.5 / (c_sound + ALPHA * c_sig + 1e-20)
    dt_cfl = factor_cfl*dt_cfl_array.min()

    # 2) dt from force/acceleration condition
    dt_force_array = np.sqrt(h*0.5 / (a_magnitude + 1e-20))
    dt_force = factor_acc*dt_force_array.min()
    
    # 3) dt from expansion
    dt_u=1e20
#    if dudt < 0.0:
    dt_u_array = u/(np.abs(dudt)+1e-20)
    dt_u = factor_u*dt_u_array.min()

    # Combine by taking the minimum
    dt_crit = min(dt_cfl, dt_force, dt_u)


    return dt_crit


def leapfrog_step(particles,box_size, nsmooth, periodic_box,dt):
    # Step 1: Half kick (update velocities by half timestep)
    particles[:, 4:7] += 0.5 * dt * particles[:, 12:15]
    # Half-step internal energy update
    particles[:, 8] += 0.5 * dt * particles[:, 15]
    
    
    
    # Step 2: Drift (update positions by full timestep)
    particles[:, 1:4] += dt * particles[:, 4:7]
    
    # Predict value of v and u at t+1
    v_old = particles[:, 4:7].copy()
    u_old = particles[:, 8].copy()
    particles[:, 4:7] += 0.5 * dt * particles[:, 12:15]
    
    particles[:, 8] += 0.5 * dt * particles[:, 15]
    
    

    # Apply periodic boundary conditions if necessary
    if periodic_box:
        particles[:, 1:4],particles[:, 4:7] = apply_boundary_condition_all(particles[:, 1:4], box_size,velocities=particles[:, 4:7])

    
    particles=compute_forces(particles,box_size, nsmooth, periodic_box)
    
    
    # Step 4: Second half kick (final velocity update)
    # use the old values and evolve with new velocities
    particles[:, 4:7] = v_old + 0.5 * dt * particles[:, 12:15]
    # Internal energy
    
    particles[:, 8] = u_old + 0.5 * dt * particles[:, 15]
    

    
    return particles
    
    
def compute_forces(particles,box_size, nsmooth, periodic_box):
    particles=compute_sph_forces(particles,box_size, nsmooth, periodic_box)
    return particles
    
def compute_sph_forces(particles,box_size, nsmooth, periodic_box):
    # Step 3: Recalculate accelerations after drift
    start_time = time.time()
    tree = build_tree(particles[:, 1:4], box_size=box_size, h=particles[:, 20], periodic_box=periodic_box)
    tree_build_time = time.time() - start_time
    print(f"Time to build kd tree: {tree_build_time:.4f} seconds")
    
    
    # ## DEBUG PARTICLES REP
    # particles_rep=np.copy(tree[2])

    # particles_rep[:,2]=particles_rep[:,1]
    # particles_rep[:,1]=particles_rep[:,0]


    # colors = ['blue', 'red', 'red','orange','purple']
    # plot_density_vs_smoothing_length(particles_rep,particles, colors=colors,density_index=1,smoothing_length_index=2,log_x=False,log_y=False)


    start_time = time.time()
    h_smooth, neighbors_list, distances, vec = calculate_smoothing_length_staged(
        tree, particles[:, 1:4], box_size=box_size, periodic_box=periodic_box, final_k=nsmooth
    )
    smoothing_length_time = time.time() - start_time
    print(f"Time to calculate smoothing lengths and neighbors (Staged): {smoothing_length_time:.4f} seconds")

    start_time = time.time()
    density = calculate_density(
        particles[:, 1:4], particles[:, 0], h_smooth, neighbors_list, distances, 
        box_size=box_size, periodic_box=periodic_box
    )
    density_calculation_time = time.time() - start_time
    print(f"Time to calculate density: {density_calculation_time:.4f} seconds")
    particles[:, 7] = density            
    particles[:, 20] = h_smooth
    
    start_time = time.time()
    particles = EOS(particles)
    eos_time = time.time() - start_time
    print(f"Time to apply EOS: {eos_time:.4f} seconds")

    start_time = time.time()
    Q0, Q1x, Q1y, Q1z, E0x, E0y, E0z, dvxdt, dvydt, dvzdt,dudt = calculate_gradients(
        particles, neighbors_list, distances, vec, box_size=box_size, periodic_box=periodic_box
    )
    gradient_calculation_time = time.time() - start_time
    print(f"Time to calculate gradients: {gradient_calculation_time:.4f} seconds")
    
    #dvydt -= 9.8
    #print(np.sum(dvxdt),np.sum(dvydt),np.sum(dvxdt))
    
    
    particles[:, 12] = dvxdt
    particles[:, 13] = dvydt
    particles[:, 14] = dvzdt
    particles[:, 15] = dudt
    return particles
   

def implicit_midpoint_step(particles, box_size, nsmooth, periodic_box, dt, tol=1e-3, max_iter=4):
    """
    Perform a single implicit midpoint integration step for SPH particles.
    This method is symplectic and time-reversible.

    Parameters:
    - particles: ndarray, particle array with positions, velocities, etc.
    - box_size: ndarray, simulation box dimensions
    - nsmooth: int, number of smoothing neighbors
    - periodic_box: bool, whether to apply periodic boundaries
    - dt: float, timestep
    - tol: float, tolerance for convergence
    - max_iter: int, maximum number of iterations for convergence

    Returns:
    - particles: ndarray, updated particle properties
    """
    # Store initial values for iterative process
    x_old = particles[:, 1:4].copy()
    v_old = particles[:, 4:7].copy()
    u_old = particles[:, 8].copy()
    
    # Predictor step (initial guess for midpoint)
    x_mid = x_old + 0.5 * dt * v_old
    v_mid = v_old + 0.5 * dt * particles[:, 12:15]  # Use initial acceleration
    u_mid = u_old + 0.5 * dt * particles[:, 15]
    
    converged = False
    it=0
    for _ in range(max_iter):
        # Apply boundary conditions to midpoint positions
        if periodic_box:
            x_mid, v_mid = apply_boundary_condition(x_mid, box_size, velocities=v_mid)

        # Recalculate forces at the midpoint
        particles[:, 1:4] = x_mid
        particles[:, 4:7] = v_mid
        particles[:, 8] = u_mid
        particles = compute_forces(particles, box_size, nsmooth, periodic_box)
        
        # Update midpoint velocities based on midpoint accelerations
        v_mid_new = v_old + 0.5 * dt * particles[:, 12:15]
        u_mid_new = u_old + 0.5 * dt * particles[:, 15]
        x_mid_new = x_old + 0.5 * dt * v_mid_new
        
        it += 1
        # Check convergence
        if np.allclose(v_mid, v_mid_new, atol=tol) and np.allclose(x_mid, x_mid_new, atol=tol):
            converged = True
            print("iterations",it)
            break
        
        # Update guesses for next iteration
        x_mid = x_mid_new
        v_mid = v_mid_new
        u_mid = u_mid_new
    
    if not converged:
        print("WARNING: Implicit Midpoint method did not converge within max_iter.")    
    
    # Final update to positions and velocities
    particles[:, 1:4] = x_old + dt * v_mid
    particles[:, 4:7] = v_old + dt * particles[:, 12:15]
    particles[:, 8] = u_old + dt * particles[:, 15]
                                               
    return particles

def PEFRL_sub_step(particles,box_size, nsmooth, periodic_box, dt,coef):
    particles[:, 4:7] += coef * dt * particles[:, 12:15]
    particles[:, 8] += coef * dt * particles[:, 15]
    #drift
    particles[:, 1:4] += coef * dt * particles[:, 4:7]
    # Apply periodic boundary conditions if necessary
    if periodic_box:
        particles[:, 1:4],particles[:, 4:7] = apply_boundary_condition_all(particles[:, 1:4], box_size,velocities=particles[:, 4:7])

    particles = compute_forces(particles, box_size, nsmooth, periodic_box)
    return particles

def PEFRL_step(particles,box_size, nsmooth, periodic_box, dt):
    # Coefficients
    c1 =  0.1786178958448091
    c2 = -0.2123418310626054
    c3 =  0.5353258025961350
    c4 =  c2
    c5 =  c1
    particles=PEFRL_sub_step(particles,box_size, nsmooth, periodic_box, dt,c1)
    particles=PEFRL_sub_step(particles,box_size, nsmooth, periodic_box, dt,c2)
    particles=PEFRL_sub_step(particles,box_size, nsmooth, periodic_box, dt,c3)
    particles=PEFRL_sub_step(particles,box_size, nsmooth, periodic_box, dt,c4)
    particles=PEFRL_sub_step(particles,box_size, nsmooth, periodic_box, dt,c5)

    return particles
    

def time_step(particles, box_size, nsmooth, periodic_box, dt_old=None,time_in=0.0):
    """
    Perform a single time step for SPH evolution using Kick-Drift-Kick (KDK) integration.

    Parameters:
    - particles: ndarray, shape (N, 7), containing mass, positions, and velocities
    - box_size: ndarray, shape (3,), size of the box (for periodic boundary)
    - nsmooth: int, number of smoothing neighbors
    - periodic_box: bool, whether to apply periodic boundary conditions
    - dt: float, timestep for integration

    Returns:
    - particles: ndarray, updated particles with new positions and velocities
    """
    start_time = time.time()
    if time==0.0:
        particles = EOS(particles)
        particles = compute_sph_forces(particles,box_size, nsmooth, periodic_box)

    dt_crit = compute_timestep_criterion(particles)
    if dt_crit < dt_old:
        print(" WARNING dt",dt_old,"dt_crit",dt_crit)
    # Override or reduce the input dt
    dt = min(dt_old, dt_crit)
    
    remain_t=dt_old-(time_in % dt_old)
    if dt > remain_t and remain_t > 1E-8:
        dt=remain_t
        
    eos_time = time.time() - start_time
    print(f"Time to apply timestepcritertion: {eos_time:.4f} seconds")
    
    particles=leapfrog_step(particles,box_size, nsmooth, periodic_box,dt)
    #particles=PEFRL_step(particles,box_size, nsmooth, periodic_box,dt)
    #particles=implicit_midpoint_step(particles,box_size, nsmooth, periodic_box,dt)
    
    return particles,dt


def evolve_sph(particles, box_size, dt, t_final, nsmooth, periodic_box):
    """
    Evolve the SPH system over multiple timesteps using Kick-Drift-Kick (KDK) integration.

    Parameters:
    - particles: ndarray, initial particles with mass, positions, and velocities
    - box_size: ndarray, size of the box
    - dt: float, timestep
    - nsteps: int, number of timesteps to evolve
    - nsmooth: int, number of smoothing neighbors
    - periodic_box: bool, whether to apply periodic boundary

    Returns:
    - particles: ndarray, evolved particles with new positions and velocities
    """
    labels = ['Original Particles', 'split  Particles', 'merge Particles','orange','purp']
    colors = ['blue', 'green', 'red','orange','purple']
    t=0.0
    step=0
    while t < t_final:
        dt_old=np.copy(dt)
        if step==0:            
            dt=dt*0.000001
            
            
        start_time0 = time.time()
        particles,dt_new = time_step(
            particles, box_size=box_size, nsmooth=nsmooth, periodic_box=periodic_box, dt_old=dt,time_in=t
        )
        step_time = time.time() - start_time0
        print(f"Step time: {step_time:.4f} seconds")
        
        
        start_time0 = time.time()
        properties = calculate_conservative_properties(particles)
        #print("Linear Momentum:", properties['linear_momentum'])
        #print("Centre of mass", properties['com'])
        #print("Kinetic Energy:", properties['kinetic_energy'])
        #print("Thermal Energy:", properties['thermal_energy'])
        #print("Total Energy:", properties['thermal_energy']+properties['kinetic_energy'])
        print("Time:", t)
        print("Step:", step)
        t += dt_new
        dt=dt_old
        step += 1
        # npxi=100
        if (t % dt) < 1E-8:
            print("out")
            #plot_density_vs_smoothing_length(particles, labels=labels, colors=colors,density_index=1,smoothing_length_index=7,log_x=False,log_y=False)
            plot_azimuthal_velocity(particles,t,properties['kinetic_energy'],properties['thermal_energy'])
        
        plot_time = time.time() - start_time0
        print(f"Plot time: {plot_time:.4f} seconds")
        # testgrid_full = interpolate3D(particles, box_size,npx=npxi,periodic_box=periodic_box)
        # plot_projection(testgrid_full[:, :, :, 7], box_size, axis='z', method='mean')
        
    
    return particles

def plot_azimuthal_velocity(particles, time,Ekin,Eth):
    Etot=Ekin+Eth
    # Extract positions and velocities
    x, y, z = particles[:, 1], particles[:, 2], particles[:, 3]
    vx, vy, vz = particles[:, 4], particles[:, 5], particles[:, 6]

    # Calculate cylindrical radius (R)
    R = np.sqrt(x**2 + y**2)

    # Calculate azimuthal velocity (v_phi)
    v_phi = (-y * vx + x * vy) / R  # Cross product in cylindrical coordinates

    # Plotting azimuthal velocity vs cylindrical radius
    plt.figure(figsize=(8, 6))
    plt.scatter(R, v_phi, s=2, alpha=0.7)
    plt.xlabel('Cylindrical Radius (R)')
    plt.ylabel('Azimuthal Velocity (v_phi)')
    plt.title(f'Azimuthal Velocity vs Cylindrical Radius\nTime = {time:.2f} {Etot:.6f} {Ekin:.6f} {Eth:.6f}')
    plt.grid(True)
    plt.show()
