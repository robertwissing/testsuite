import numpy as np
#from scipy.spatial import KDTree
from sklearn.neighbors import BallTree
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from numba import njit, prange,get_thread_id,get_num_threads
from joblib import Parallel, delayed
import time
import sys
import os
from itertools import product

# Add the folder containing numba_kdtree to the Python module search path
##odule_path = "/mn/stornext/u3/robertwi/analysis/numba-kdtree"  # Replace with the actual path
##sys.path.append(os.path.abspath(module_path))

from numba_kdtree import KDTree
from numba import config
from collections import defaultdict

from SPH_constants import *



@njit
def apply_1d_periodic(coord, L):
    """
    Map 'coord' into [-L/2, L/2] using modulo arithmetic.
    """
    return np.mod(coord + L/2, L) - L/2

@njit
def apply_1d_reflecting(coord, L):
    """
    Reflect 'coord' if it crosses +/- L/2.
    """
    half_L = L / 2.0
    if coord < -half_L:
        coord = -half_L + (-half_L - coord)  # Reflect position
    elif coord > half_L:
        coord = half_L - (coord - half_L)    # Reflect position
    return coord

@njit
def apply_1d_reflecting_vel(coord,vel, L):
    """
    Reflect 'coord' if it crosses +/- L/2.
    """
    half_L = L / 2.0
    if coord < -half_L:
        vel = -vel
    elif coord > half_L:
        vel = -vel
    return vel


@njit
def apply_boundary_condition(positions, box_size, boundary_types=np.array([PERIODIC, PERIODIC, PERIODIC]),velocities=None):
    for i in range(3):
        btype = boundary_types[i]
        if btype == PERIODIC:
            positions[i] = apply_1d_periodic(positions[i], box_size[i])
        elif btype == REFLECTING:
            if velocities is not None and velocities.size > 0:
                velocities[i] = apply_1d_reflecting_vel(positions[i],velocities[i], box_size[i])  
            positions[i] = apply_1d_reflecting(positions[i], box_size[i])   
        else:
            raise ValueError("Unknown boundary type code.")
    return positions,velocities

@njit(parallel=True)
def apply_boundary_condition_all(positions, box_size, boundary_types=np.array([PERIODIC, PERIODIC, PERIODIC]),velocities=None):
    """
    Apply boundary conditions to all particles.
    
    Parameters:
    - positions: 2D array of shape (N, 3) representing N particles' x, y, z positions.
    - box_size: 1D array of length 3 representing the size of the simulation box in each dimension.
    - boundary_types: 1D array of length 3 representing boundary types for x, y, z respectively.
    """
    N = len(positions)
    for i in prange(N):
        apply_boundary_condition(positions[i,:], box_size, boundary_types,velocities[i,:])
    
    return positions,velocities


@njit
def create_periodic_replicas_1d_numba(particles, box_size, h, dim):
    """
    Creates periodic replicas of particles near the boundaries in one dimension using Numba for acceleration.

    Parameters:
    - particles (np.ndarray): Array of particle positions with shape (N, 3).
    - box_size (np.ndarray): Size of the simulation box in each dimension (Lx, Ly, Lz).
    - h (float): Margin distance to determine proximity to boundaries.
    - dim (int): Dimension (0, 1, or 2) along which to generate periodic replicas.

    Returns:
    - Tuple[np.ndarray, np.ndarray]:
        - Array of replica particles with shape (M, 3), where M ≤ 2N.
        - Array of indices mapping replicas to their original particles with shape (M,).
    """
    N = particles.shape[0]
    replicas_min = 0
    replicas_max = 0

    # First pass: Count the number of replicas needed near min and max boundaries
    for i in range(N):
        if particles[i, dim] < (-0.5 * box_size[dim] + h):
            replicas_min += 1
        if particles[i, dim] > (0.5 * box_size[dim] - h):
            replicas_max += 1

    total_replicas = replicas_min + replicas_max

    if total_replicas == 0:
        # Return empty arrays if no replicas are needed
        return np.empty((0, 3), dtype=particles.dtype), np.empty(0, dtype=np.int64)

    # Initialize arrays for replicas and their corresponding indices
    replicas = np.zeros((total_replicas, 3), dtype=particles.dtype)
    replica_indices = np.zeros(total_replicas, dtype=np.int64)

    count = 0  # Counter for the number of replicas added

    # Second pass: Generate replicas and record their original indices
    for i in range(N):
        if particles[i, dim] < (-0.5 * box_size[dim] + h):
            # Create a replica shifted positively by box_size in the specified dimension
            replicas[count, :] = particles[i, :]
            replicas[count, dim] += box_size[dim]
            replica_indices[count] = i
            count += 1

        if particles[i, dim] > (0.5 * box_size[dim] - h):
            # Create a replica shifted negatively by box_size in the specified dimension
            replicas[count, :] = particles[i, :]
            replicas[count, dim] -= box_size[dim]
            replica_indices[count] = i
            count += 1

    return replicas, replica_indices


@njit
def create_reflected_replicas_1d_numba(particles, box_size, h, dim):
    """
    Creates reflected replicas of particles near the boundaries in one dimension using Numba for acceleration.

    Parameters:
    - particles (np.ndarray): Array of particle positions with shape (N, 3).
    - box_size (np.ndarray): Size of the simulation box in each dimension (Lx, Ly, Lz).
    - h (float): Margin distance to determine proximity to boundaries.
    - dim (int): Dimension (0, 1, or 2) along which to generate reflected replicas.

    Returns:
    - Tuple[np.ndarray, np.ndarray]:
        - Array of replica particles with shape (M, 3), where M ≤ 2N.
        - Array of indices mapping replicas to their original particles with shape (M,).
    """
    N = particles.shape[0]
    replicas_min = 0
    replicas_max = 0

    # Define boundaries
    lower_bound = -0.5 * box_size[dim]
    upper_bound = 0.5 * box_size[dim]

    # First pass: Count the number of replicas needed near min and max boundaries
    for i in range(N):
        if particles[i, dim] < (lower_bound + h):
            replicas_min += 1
        if particles[i, dim] > (upper_bound - h):
            replicas_max += 1

    total_replicas = replicas_min + replicas_max

    if total_replicas == 0:
        # Return empty arrays if no replicas are needed
        return np.empty((0, 3), dtype=particles.dtype), np.empty(0, dtype=np.int64)

    # Initialize arrays for replicas and their corresponding indices
    replicas = np.zeros((total_replicas, 3), dtype=particles.dtype)
    replica_indices = np.zeros(total_replicas, dtype=np.int64)

    count = 0  # Counter for the number of replicas added

    # Second pass: Generate reflected replicas and record their original indices
    for i in range(N):
        # Check lower boundary
        if particles[i, dim] < (lower_bound + h):
            replicas[count, :] = particles[i, :]
            # Reflect across the lower boundary
            replicas[count, dim] = 2 * lower_bound - particles[i, dim]
            replica_indices[count] = i
            count += 1

        # Check upper boundary
        if particles[i, dim] > (upper_bound - h):
            replicas[count, :] = particles[i, :]
            # Reflect across the upper boundary
            replicas[count, dim] = 2 * upper_bound - particles[i, dim]
            replica_indices[count] = i
            count += 1

    return replicas, replica_indices


def create_boundary_replicas_diff(
    particles,
    box_size,
    h,
    boundary_types=np.array([PERIODIC, PERIODIC, PERIODIC])):
    """
    Creates boundary replicas for particles based on specified boundary conditions per dimension.

    Parameters:
    - particles (np.ndarray): Array of particle positions with shape (N, 3).
    - box_size (Tuple[float, float, float]): Size of the simulation box in each dimension (Lx, Ly, Lz).
    - h (float): Margin distance to determine proximity to boundaries.
    - boundary_types (np.ndarray): Array indicating boundary types for each dimension (0: PERIODIC, 1: REFLECTING).

    Returns:
    - Tuple[np.ndarray, np.ndarray]:
        - Combined array of original and replica particles.
        - Array of indices mapping replicas to their original particles.
    """
    particles = np.asarray(particles)
    box_size = np.asarray(box_size)
    boundary_types = np.asarray(boundary_types)
    h=np.max(h)

    # Lists to store per-dimension replicas
    per_dim_replicas = []
    per_dim_indices = []

    replica_count=0

    start_time = time.time()
    # Iterate through each dimension and generate replicas based on boundary type
    for dim in range(3):
        if boundary_types[dim] == PERIODIC:
            replicas, indices = create_periodic_replicas_1d_numba(particles, box_size, h, dim)
            if replicas.size > 0:
                per_dim_replicas.append(replicas)
                per_dim_indices.append(indices)
        elif boundary_types[dim] == REFLECTING:
            replicas, indices = create_reflected_replicas_1d_numba(particles, box_size, h, dim)
            if replicas.size > 0:
                per_dim_replicas.append(replicas)
                per_dim_indices.append(indices)
        else:
            raise ValueError(f"Unknown boundary type code: {boundary_types[dim]} for dimension {dim}.")
            
    print("time1",time.time()-start_time)
    start_time = time.time() 
    if not per_dim_replicas:
        # No replicas needed
        return particles.copy(), np.arange(len(particles))

    # To handle multi-dimensional replicas, generate all possible combinations of per-dimension replicas
    # First, group replicas by their originating particles
    from collections import defaultdict

    replica_dict = defaultdict(dict)  # {particle_index: {dim: [replicas]}}

    for dim, (replicas, indices) in enumerate(zip(per_dim_replicas, per_dim_indices)):
        for replica, idx in zip(replicas, indices):
            replica_dict[idx][dim] = replica

    # Now, for each particle that has replicas in multiple dimensions, create combined replicas
    combined_replicas = []
    combined_indices = []

    print("time2",time.time()-start_time)
    start_time = time.time() 
    
    for idx, dim_replicas in replica_dict.items():
        # List to store operations: original or per-dimension replica
        # Start with the original position
        original = particles[idx]

        # Collect shift/reflection options per dimension
        options = []
        for dim in range(3):
            if dim in dim_replicas:
                # Include the shifted/reflected replica and the original position
                # This allows combinations where some dimensions are shifted/reflected and others are not
                options.append([original.copy(), dim_replicas[dim]])
            else:
                # Only the original position is available
                options.append([original.copy()])
        # Generate all combinations excluding the all-original combination
        all_combinations = list(product(*options))
        # Exclude the combination where all dimensions are original

        all_combinations = [
            comb for comb in all_combinations
            if not all(np.array_equal(comb[dim], original) for dim in range(3))
        ]

        # Add each unique combination as a replica
        for comb in all_combinations:
            # Create a new replica by selecting the appropriate coordinate in each dimension
            replica = original.copy()
            for dim in range(3):
                if not np.array_equal(comb[dim], original[dim]):
                    replica[dim] = comb[dim][dim]
            combined_replicas.append(replica)
            combined_indices.append(idx)
            replica_count += 1
    print(f"Total replicas created: {replica_count}")

    print("time3",time.time()-start_time)
    if combined_replicas:
        combined_replicas = np.array(combined_replicas)
        combined_indices = np.array(combined_indices)
        print(f"Total replicas created: {len(combined_replicas)}")
        print(f"Total particles: {len(particles)}")
        # Combine original particles with replicas
        combined_particles = np.vstack([particles, combined_replicas])
        # Map replicas to original particles
        combined_indices = np.concatenate([np.arange(len(particles)), combined_indices])
        return combined_particles, combined_indices
    else:
        # No replicas generated beyond per-dimension
        return particles.copy(), np.arange(len(particles))


def create_boundary_replicas(
    particles,
    box_size,
    h,
    boundary_types=np.array([PERIODIC, PERIODIC, PERIODIC])
):
    """
    Create boundary replicas for particles based on specified boundary conditions.

    Parameters:
    - particles: Array or list of particle positions.
    - box_size: Size of the simulation box.
    - h: Smoothing length
    - boundary_types: NumPy array specifying boundary conditions for each dimension (0 for PERIODIC, 1 for REFLECTING).

    Returns:
    - Replicated particles based on boundary conditions.
    """
    # Define the default periodic and reflecting arrays
    periodic = np.array([PERIODIC, PERIODIC, PERIODIC])
    reflecting = np.array([REFLECTING, REFLECTING, REFLECTING])
    lattice = np.array([LATTICE, LATTICE, LATTICE])
    if np.array_equal(boundary_types, periodic):
        return create_periodic_replicas_old(particles, box_size, h)
        #return create_boundary_replicas_diff(particles, box_size, h, boundary_types=boundary_types)
    elif np.array_equal(boundary_types, reflecting):
        return create_reflecting_replicas_old(particles, box_size, h)
        #return create_boundary_replicas_diff(particles, box_size, h, boundary_types=boundary_types)
    elif np.array_equal(boundary_types, lattice):
        return create_lattice_replicas_old(particles, box_size, h)
    else:
        return create_boundary_replicas_diff(particles, box_size, h, boundary_types=boundary_types)


def create_periodic_replicas(particles, box_size, h):
    """
    Creates periodic replicas of particles near the boundaries of a simulation box.

    Parameters:
    - particles (np.ndarray): Array of particle positions with shape (N, 3).
    - box_size (Tuple[float, float, float]): Size of the simulation box in each dimension (Lx, Ly, Lz).
    - h (float or np.ndarray): Margin distance to determine proximity to boundaries.
                                  Can be a scalar or array-like with three elements for each dimension.

    Returns:
    - Tuple[np.ndarray, np.ndarray]:
        - Combined array of original and replica particles.
        - Array of indices mapping replicas to their original particles.
    """
    particles = np.asarray(particles)
    box_size = np.asarray(box_size)
    
    h=np.max(h)*1.0
    # Identify particles near each boundary
    near_min = particles < (-0.5 * box_size + h)
    near_max = particles > (0.5 * box_size - h)
    
    # Combined condition: particles near any boundary
    boundary_condition = near_min.any(axis=1) | near_max.any(axis=1)
    boundary_particles = particles[boundary_condition]
    boundary_indices = np.nonzero(boundary_condition)[0]
    
    replicas = []
    replica_indices = []
    
    # Iterate over each boundary particle to determine necessary shifts
    for i, p in enumerate(boundary_particles):
        shift_options = []
        for dim in range(3):
            if p[dim] < (-0.5 * box_size[dim] + h):
                # Particle is near the negative boundary in this dimension; shift positively
                shift_options.append([1, 0])  # Shift by +1 or no shift
            elif p[dim] > (0.5 * box_size[dim] - h):
                # Particle is near the positive boundary in this dimension; shift negatively
                shift_options.append([-1, 0])  # Shift by -1 or no shift
            else:
                # Particle is not near any boundary in this dimension; no shift
                shift_options.append([0])
        
        # Generate all possible non-zero shift combinations
        all_shift_combinations = list(product(*shift_options))
        # Exclude the (0, 0, 0) shift to avoid duplicating original particles
        all_shift_combinations = [shift for shift in all_shift_combinations if any(s != 0 for s in shift)]
        
        # Create replicas based on shift combinations
        for shift in all_shift_combinations:
            shift_vector = np.array(shift) * box_size
            replica = p + shift_vector
            replicas.append(replica)
            replica_indices.append(boundary_indices[i])
    
    if replicas:
        replicas = np.array(replicas)
        replica_indices = np.array(replica_indices)

        # Combine original particles with replicas
        combined_particles = np.vstack([particles, replicas])
        combined_indices = np.concatenate([np.arange(len(particles)), replica_indices])
        return combined_particles, combined_indices
    else:
        # No replicas needed
        return particles, np.arange(len(particles))

def create_reflecting_replicas(particles, box_size, h):
    # Margin to determine if particles are near the boundaries
    margin = h
    boundary_particles = []

    # Identify particles near each boundary (adjusted for [-0.5, 0.5] range)
    near_x_min = particles[:,0] < (-box_size[0] / 2 + margin)
    near_x_max = particles[:,0] > (box_size[0] / 2 - margin)
    near_y_min = particles[:,1] < (-box_size[1] / 2 + margin)
    near_y_max = particles[:,1] > (box_size[1] / 2 - margin)
    near_z_min = particles[:,2] < (-box_size[2] / 2 + margin)
    near_z_max = particles[:,2] > (box_size[2] / 2 - margin)

    # Combine all conditions to select boundary particles
    boundary_condition = near_x_min | near_x_max | near_y_min | near_y_max | near_z_min | near_z_max
    boundary_particles = particles[boundary_condition]
    boundary_indices = np.arange(len(particles))[boundary_condition]  # Save original indices

    # Reflect particles that are near the boundaries
    reflected_particles = []
    reflected_indices = []
    
    for i, particle in enumerate(boundary_particles):
        reflected_particle = apply_reflecting_boundary(particle.copy(), box_size)
        reflected_particles.append(reflected_particle)
        reflected_indices.append(boundary_indices[i])

    # Combine original and reflected particles
    if len(reflected_particles) > 0:
        reflected_particles = np.array(reflected_particles)
        reflected_indices = np.array(reflected_indices)

        # Return the original particles along with their reflected replicas
        return np.concatenate([particles, reflected_particles]), np.concatenate([np.arange(len(particles)), reflected_indices])
    else:
        return particles, np.arange(len(particles))  # No reflected particles, return original


def create_reflecting_replicas_old(particles, box_size, h):
    # Margin to determine if particles are near the boundaries
    margin = h
    boundary_particles = []
    
    half_box=box_size / 2

    # Identify particles near each boundary (adjusted for [-0.5, 0.5] range)
    near_x_min = particles[:,0] < (-box_size[0] / 2 + margin)
    near_x_max = particles[:,0] > (box_size[0] / 2 - margin)
    near_y_min = particles[:,1] < (-box_size[1] / 2 + margin)
    near_y_max = particles[:,1] > (box_size[1] / 2 - margin)
    near_z_min = particles[:,2] < (-box_size[2] / 2 + margin)
    near_z_max = particles[:,2] > (box_size[2] / 2 - margin)

    # Combine all conditions to select boundary particles
    boundary_condition = near_x_min | near_x_max | near_y_min | near_y_max | near_z_min | near_z_max
    boundary_particles = particles[boundary_condition]
    boundary_indices = np.arange(len(particles))[boundary_condition]  # Save original indices

    shifts = [-1, 0, 1]
    replicas = []
    replica_indices = []  # List to store the indices of the replicas

    # Create replicas by shifting only boundary particles
    for shift_x in shifts:
        for shift_y in shifts:
            for shift_z in shifts:
                if shift_x == 0 and shift_y == 0 and shift_z == 0:
                    continue  # Skip the original particles
                    
                mirrored = boundary_particles.copy()
                # Mirror across the original box boundaries
                if shift_x != 0:  # Reflect in x if box was shifted in x
                    mirrored[:, 0] = 2 * (shift_x * half_box[0]) - boundary_particles[:, 0]
                
                if shift_y != 0:  # Reflect in y if box was shifted in y
                    mirrored[:, 1] = 2 * (shift_y * half_box[1]) - boundary_particles[:, 1]
                
                if shift_z != 0:  # Reflect in z if box was shifted in z
                    mirrored[:, 2] = 2 * (shift_z * half_box[2]) - boundary_particles[:, 2]                
                

                replicas.append(mirrored)
                replica_indices.append(boundary_indices)  # Keep track of the original indices


    near_x_min = replicas[:,0] < (-box_size[0] / 2 - margin)
    near_x_max = replicas[:,0] > (box_size[0] / 2 + margin)
    near_y_min = replicas[:,1] < (-box_size[1] / 2 - margin)
    near_y_max = replicas[:,1] > (box_size[1] / 2 + margin)
    near_z_min = replicas[:,2] < (-box_size[2] / 2 - margin)
    near_z_max = replicas[:,2] > (box_size[2] / 2 + margin)
    boundary_condition = near_x_min | near_x_max | near_y_min | near_y_max | near_z_min | near_z_max
    boundary_particles = replicas[~boundary_condition]
    replica_indices = replica_indices[~boundary_condition]

    if len(replicas) > 0:
        all_replicas = np.concatenate(replicas)
        all_indices = np.concatenate(replica_indices)
        # Return the original particles, replicas, and corresponding original indices
        return np.concatenate([particles, all_replicas]), np.concatenate([np.arange(len(particles)), all_indices])
    else:
        return particles, np.arange(len(particles))  # No replicas, just return original particles and their indices


def create_periodic_replicas_old(particles, box_size, h):
    # Margin to determine if particles are near the boundaries
    margin = h
    particles = particles
    indices = np.arange(len(particles))  # Save original indices


    shifts = [-1, 0, 1]
    replicas = []
    replica_indices = []  # List to store the indices of the replicas

    # Create replicas by shifting only boundary particles
    for shift_x in shifts:
        for shift_y in shifts:
            for shift_z in shifts:
                if shift_x == 0 and shift_y == 0 and shift_z == 0:
                    continue  # Skip the original particles
                shift_vector = np.array([shift_x, shift_y, shift_z]) * box_size
                replica = particles + shift_vector
                replicas.append(replica)
                replica_indices.append(indices)  # Keep track of the original indices




    if len(replicas) > 0:
        all_replicas = np.concatenate(replicas)
        all_indices = np.concatenate(replica_indices)
        near_x_min = all_replicas[:,0] < (-box_size[0] / 2 - margin)
        near_x_max = all_replicas[:,0] > (box_size[0] / 2 + margin)
        near_y_min = all_replicas[:,1] < (-box_size[1] / 2 - margin)
        near_y_max = all_replicas[:,1] > (box_size[1] / 2 + margin)
        near_z_min = all_replicas[:,2] < (-box_size[2] / 2 - margin)
        near_z_max = all_replicas[:,2] > (box_size[2] / 2 + margin)
        boundary_condition = near_x_min | near_x_max | near_y_min | near_y_max | near_z_min | near_z_max
        all_replicas = all_replicas[~boundary_condition]
        all_indices = all_indices[~boundary_condition]
        # Return the original particles, replicas, and corresponding original indices
        return np.concatenate([particles, all_replicas]), np.concatenate([np.arange(len(particles)), all_indices])
    else:
        return particles, np.arange(len(particles))  # No replicas, just return original particles and their indices


# Function to build either a KDTree or a BallTree from the given particles

def build_tree(particles, box_size,h, periodic_box):
    tree_type='kd'
    particles_rep=np.copy(particles)
    indices = np.arange(len(particles))
    leaf_size=24
    start_time = time.time()
    if periodic_box and box_size is not None:
        particles_rep,indices= create_boundary_replicas(particles, box_size,np.max(h))
        
    # Choose the type of tree to build
    if tree_type == 'kd':
        # Using scipy's KDTree
        return KDTree(particles_rep,leafsize=leaf_size), indices, particles_rep,tree_type
    elif tree_type == 'ball':
        # Using scikit-learn's BallTree
        return BallTree(particles_rep,leaf_size=leaf_size, metric='euclidean'), indices, particles_rep,tree_type
    


@njit(parallel=True)
def mytree_query2(particles, k, tree):
    N = len(particles)
    h = np.zeros(N, dtype=PRECISION)
    distances = np.zeros((N, k), dtype=PRECISION)
    neighbors = np.zeros((N, k), dtype=np.int64)
    
    for i in prange(N):  # Parallel loop over all particles
        dist, neigh, nsmooth = tree.query(particles[i], k=k)
        distances[i, :] = dist
        neighbors[i, :] = neigh
        h[i]=distances[i,-1]
        
    return distances, neighbors, h



def mytree_query_ball(particle,rball,tree,tree_type):
    if tree_type == 'kd':
        neighbors = tree.query_radius(particle, r=rball)
    if tree_type == 'ball':
        neighbors = tree.query_radius(particle, r=rball)

    return neighbors


@njit(parallel=True)
def find_zero_distance_mask(distances):
    """
    Returns a boolean mask indicating which particles
    have >=2 neighbors at zero distance.
    """
    N = distances.shape[0]
    mask = np.zeros(N, dtype=np.bool_)

    for i in prange(N):
        zero_distance_count = np.sum(distances[i, :] == 0)
        if zero_distance_count >= 2:
            mask[i] = True

    return mask

def remove_exact_duplicates(particles):
    """
    Removes exact-position duplicates from `particles`.

    Assumes that:
    - The first 3 columns in `particles` are x, y, z.
    - If two particles have exactly the same (x,y,z), keep the first one, remove the others.
    """
    # unique_positions: Each unique (x,y,z) combination
    # idx_unique: The index of the *first* occurrence of each unique position
    rounded_positions = np.round(particles[:, 1:4], decimals=8)
    unique_positions, idx_unique = np.unique(rounded_positions, axis=0, return_index=True)
    
    # Sort idx_unique to keep the order in which we found them
    # (np.unique may not return idx_unique in ascending order)
    idx_unique_sorted = np.sort(idx_unique)
    
    # Build a new array with only the unique particles
    particles_unique = particles[idx_unique_sorted]
    
    return particles_unique


def estimate_smoothing_length(particles, box_size, nSmooth):
    # Number of particles
    N_particles = len(particles)
    
    # Volume of the box (assuming a 3D box)
    V_box = np.prod(box_size)
    
    # Estimate number density
    n_density = N_particles / V_box
    
    # Rough estimate of the smoothing length based on nSmooth neighbors
    h_estimate = (3 * nSmooth / (4 * np.pi * n_density)) ** (1/3)
    
    return h_estimate



    
#wendland
@njit
def kernel(r,H):
    h=0.5*H
    sigma=(21/16.)/(np.pi*h**3)
    au= r/H
    ak=0.0
    if au <= 0.0:
        ak=0.95693359
    else:
        ak = 1-au;                                                  
        ak = ak*ak;                                                 
        ak = ak*ak;                                                 
        ak = ak*(1+4*au);
    return sigma*ak

@njit
def Dkernel(r,H):
    h=0.5*H
    sigma=(21/16)/(np.pi*h**5)
    q=r/H
    if (q <= 0.0):     
        adk=0.0
    else:       
        adk = 1-q;                                     
        adk = (-20./4.)*adk*adk*adk;     
    return sigma*adk


@njit
def gaussian_kernel(r, h):
    """A Gaussian SPH kernel function."""
    q = r / h
    return np.exp(-q**2) / (h**3 * np.pi**1.5)


def read_box_size_from_param(filename_param):
    box_size = np.array([1.0, 1.0, 1.0])  # Default box size if not set in file
    
    # Open the parameter file and extract the necessary values
    with open(filename_param, 'r') as file:
        for line in file:
            if 'dxPeriod' in line:
                dxPeriod = float(line.split('=')[-1].strip())
            elif 'dyPeriod' in line:
                dyPeriod = float(line.split('=')[-1].strip())
            elif 'dzPeriod' in line:
                dzPeriod = float(line.split('=')[-1].strip())
    
    # Set the box size based on the read values
    box_size = np.array([dxPeriod, dyPeriod, dzPeriod])
    
    return box_size

@njit
def rotation_matrix(axis, theta):
    """
    Compute the rotation matrix for a counterclockwise rotation around the given axis by theta radians.
    """
    axis = axis / np.linalg.norm(axis)
    a = np.cos(theta / 2)
    b, c, d = -axis * np.sin(theta / 2)
    return np.array([
        [a*a + b*b - c*c - d*d, 2*(b*c + a*d),       2*(b*d - a*c)],
        [2*(b*c - a*d),       a*a + c*c - b*b - d*d, 2*(c*d + a*b)],
        [2*(b*d + a*c),       2*(c*d - a*b),       a*a + d*d - b*b - c*c]
    ])

# Function to split a particle into two daughter particles with periodic boundary conditions
@njit
def split_particle_periodic(particle,nearest_distance,vectors_to_neighbor,ort_neigh_vector, smoothing_length, box_size,periodic_box, method="nn_ort",dist=0.4):

    if method == "random":
        # Random displacement within smoothing length
        displacement =  np.random.randn(3)  # Random 3D direction
        displacement /= np.linalg.norm(displacement)  # Normalize direction
        displacement *= np.random.uniform(0, smoothing_length) # Scale by random distance within smoothing length
    
    elif method == "nn_rand":
        # Random displacement within nearest neighbour k
        displacement =  np.random.randn(3)  # Random 3D direction
        displacement /= np.linalg.norm(displacement)  # Normalize direction
        displacement *= np.random.uniform(0, nearest_distance) # Scale by random distance within smoothing length
    elif method == "h_shell":
        # Random displacement within nearest neighbour k
        displacement =  np.random.randn(3)  # Random 3D direction
        displacement /= np.linalg.norm(displacement)  # Normalize direction
        displacement *= 3*smoothing_length
        
    elif method == "nn_ort":
        # Use nearest neighbor information to split along a specific direction
        dx, dy, dz = vectors_to_neighbor
        dr = nearest_distance
        drdh=dr/smoothing_length
        
        unit_vector_to_neighbor=vectors_to_neighbor/dr
        
        orthogonal_vector = np.cross(unit_vector_to_neighbor, ort_neigh_vector)
        norm_orthogonal = np.linalg.norm(orthogonal_vector)
        epsilon = 1e-8
        if norm_orthogonal < epsilon:
            # Vectors are parallel, define an alternative orthogonal vector
            # Choose an arbitrary vector orthogonal to unit_vector_to_neighbor
            abs_unit = np.abs(unit_vector_to_neighbor)
            if abs_unit[0] < abs_unit[1] and abs_unit[0] < abs_unit[2]:
                temp_vec = np.array([1.0, 0.0, 0.0], dtype=PRECISION)
            elif abs_unit[1] < abs_unit[2]:
                temp_vec = np.array([0.0, 1.0, 0.0], dtype=PRECISION)
            else:
                temp_vec = np.array([0.0, 0.0, 1.0], dtype=PRECISION)
            orthogonal_vector = np.cross(unit_vector_to_neighbor, temp_vec)

        orthogonal_vector /= np.linalg.norm(orthogonal_vector)  # Normalize to get a unit vector
        

        #print("myort vec in split particle",orthogonal_vector,ort_neigh_vector)

        #to keep balance between locality and not increasing density too much when splitting several times.
        #dist = (1.0 / (2.25*drdh))
        #dist =  1/np.sqrt(2) # for equal lengt
        #dist = 1/np.sqrt(2)
        #dist = 0.4
            
        displacement=orthogonal_vector*dr*dist

    else:
        raise ValueError("Invalid method or missing nearest neighbor/nearest neighbor distance for orthogonal method.")
   


    particle_1 = particle + displacement
    particle_2 = particle - displacement

    has_nan_1 = np.any(np.isnan(particle_1))
    has_nan_2 = np.any(np.isnan(particle_2))
    if has_nan_1 or has_nan_2:
        print("NaN detected in particle_1 or particle_2.")
        print(dr,norm_orthogonal,particle,displacement,particle_1,particle_2,orthogonal_vector,unit_vector_to_neighbor,temp_vec)

    if periodic_box:
        particle_1, _ = apply_boundary_condition(particle_1, box_size)
        particle_2, _ = apply_boundary_condition(particle_2, box_size)
    
    return particle_1, particle_2, orthogonal_vector

# Function to split all particles in a given array, considering periodic boundary conditions
@njit
def split_particless(particles, neighbors_list,nearest_distances,vectors_to_neighbors, box_size, periodic_box,method="nn_ort",kin=1,split_target=None):
    N = len(particles)
    
    
    # Count particles to be split and allocate space for new particles
    num_splits = np.sum(split_target) if split_target is not None else N
    new_particles = np.zeros((num_splits * 2, PARTICLE_ARRAY_SIZE), dtype=PRECISION)
    
    # Initialize new particle array to store all split particles
    #new_particles = np.zeros((N * 2, 19), dtype=np.float32)
    
      
    parent_masses=particles[:,0]
    parent_particles_pos=particles[:,1:4]
    parent_particles_vel=particles[:,4:7]
    parent_density=particles[:,7]
    parent_h=particles[:,20]
    parent_u=particles[:,8]
    parent_soft=particles[:,9]
    parent_metals=particles[:,10]
    parent_B=particles[:,12:15]
    parent_spin=particles[:,15:18]
    parent_momenti=particles[:,18]
    parent_tform=particles[:,19]
    parent_volume=parent_masses/parent_density
    
    my_vector = np.zeros((N, 3), dtype=PRECISION) + 100.0  # Initialize with a default value
    input_vector = np.zeros((1, 3), dtype=PRECISION)
    split_idx = 0  # Tracks position in new_particles array
    
    count=0
    counth=0
    
    for i in prange(N):
        if split_target is not None and not split_target[i]:
            continue  # Skip particles not marked for splitting
        dist=0.4
        k=kin
        nsmooth=len(neighbors_list[i])
        nearest_distanceoverhmin=1/(0.5*nsmooth)**(1/3)
        #print(nsmooth,nearest_distanceoverhmin)
        havg = parent_h[i]
        while k < nsmooth - 1:
             nearest_distanceoverh = nearest_distances[i,k]/havg
             #print(nearest_distanceoverh)
             #print(nearest_distanceoverh)
             #if(parent_particles_pos[i,0] < -2.584 and parent_particles_pos[i,0]>-2.588 and parent_particles_pos[i,1] < -1.444 and parent_particles_pos[i,1]>-1.447):
             #    print("k, h and distanceh ",k,havg,nearest_distanceoverh,nearest_distances[i,:],nearest_distances[i,:]/havg)
             #1/(N/2.0)**(1/3) Is about minimum distance for totally relaxed glass
             #Ensure that splitting length is atleast around this length
             
             if(nearest_distanceoverh<nearest_distanceoverhmin):
                  dist = 0.4
                  k += 1
                  #print("smaller than limit h (1/(N)**(1/3)) ",nearest_distanceoverh,k,nearest_distances[i,k],i,neighbors_list[i,k],vectors_to_neighbors[i,k])  
             else:   
                 break

        particle_to_split = parent_particles_pos[i]
        nearest_neighbor = neighbors_list[i,k]
        nearest_distance = nearest_distances[i,k]
        vectors_to_neighbor=vectors_to_neighbors[i,k]
        havg = parent_h[i]

        # Split the particle into two, considering periodic boundary conditions
        if my_vector[nearest_neighbor,0]==100:
            #generate two unit vectors to nearest neighbour and to second nearest neighbour
            #input_vector=vectors_to_neighbors[i,k+1]
            input_vector=vectors_to_neighbors[i, k+1].astype(PRECISION)
            #input_angle = np.arctan2(dy, dx)+np.pi  # Angle between 0 and 2π
            
        else:
            input_vector=my_vector[nearest_neighbor].astype(PRECISION)
            if neighbors_list[nearest_neighbor,k] != i:
                input_vector=vectors_to_neighbors[i,k+1].astype(PRECISION)




        #print(f"vectors_to_neighbors dtype: {vectors_to_neighbors.dtype}")
        if(nearest_distance==0):
            p1_pos, p2_pos, vector = split_particle_periodic(particle_to_split,nearest_distance,vectors_to_neighbor,input_vector, havg, box_size=box_size,periodic_box=periodic_box,method="h_shell",dist=dist)
        else:
            p1_pos, p2_pos, vector = split_particle_periodic(particle_to_split,nearest_distance,vectors_to_neighbor,input_vector, havg, box_size=box_size,periodic_box=periodic_box,method=method,dist=dist)
            my_vector[i]=vector
                
            # Rotate orthogonal vector by half a turn in respect to angle of neighbour
            #input_angle=0.5 * np.pi - my_angle[nearest_neighbor]
         
        if k>kin:
            count+=1
        # Store the new positions and masses
        new_particles[split_idx, 1:4] = p1_pos
        new_particles[split_idx + 1, 1:4] = p2_pos
        new_particles[split_idx, 0] = parent_masses[i] / 2  # Halve the mass for each daughter particle
        new_particles[split_idx + 1, 0] = parent_masses[i] / 2
        
       # To retain angular momentum and kinetic energy, we need to set two of v_omega,L_daugther or I_daugther
       # So we have opportunity to transfer some spin to translational through v_omega.
       # If we just half the angular momentum spin, we need to adjust I such that Ekin is conserved.
       # So two methods 
       #Here we split inertia and spin in two, conserving both. 
       
        
        # Split velocity, conserve momentum
        new_particles[split_idx, 4:7] = parent_particles_vel[i]
        new_particles[split_idx + 1, 4:7] = parent_particles_vel[i]
        
        # Split internal energy (u), distribute equally or adjust as needed
        
        new_particles[split_idx, 7] = parent_density[i]
        new_particles[split_idx + 1, 7] = parent_density[i]
        
        new_particles[split_idx, 20] = parent_h[i]
        new_particles[split_idx + 1, 20] = parent_h[i]
        
        new_particles[split_idx, 8] = parent_u[i] / 2
        new_particles[split_idx + 1, 8] = parent_u[i] / 2
        
        new_particles[split_idx, 9] = parent_soft[i]
        new_particles[split_idx + 1, 9] = parent_soft[i]
        
        new_particles[split_idx, 10] = parent_metals[i]
        new_particles[split_idx + 1, 10] = parent_metals[i]
        
        # Split magnetic field, conserve flux
        new_particles[split_idx, 12:15] = parent_B[i] / 2
        new_particles[split_idx + 1, 12:15] = parent_B[i] / 2
        
        # Split spin (angular momentum), adjust based on the splitting method
        new_particles[split_idx, 15:18] = parent_spin[i]
        new_particles[split_idx + 1, 15:18] = parent_spin[i]
        
        # Split moment of inertia based on new sizes
        new_particles[split_idx, 18] = parent_momenti[i] / 2
        new_particles[split_idx + 1, 18] = parent_momenti[i] / 2
        
        # Same tform as before
        new_particles[split_idx, 19] = parent_tform[i]
        new_particles[split_idx + 1, 19] = parent_tform[i]
        # Increment split index for the next pair of split particles
        split_idx += 2
        
        counth += k
    
    # Filter out the particles that were not split
    unsplit_particles = particles[~split_target]
    
    # Combine unsplit particles with the newly generated split particles
    combined_particles = np.vstack((unsplit_particles, new_particles))
    print("THIS IS HOW MANY WAS BELOW nearest_distanceoverhmin and mean k ", count,counth/N)
    return combined_particles

@njit(parallel=True)
def sort_neighbours(N,particles,prt_rep,prtidx,neighbors_initial,distances_initial,periodic_box,final_k,h,neighbor_list,distances_list,vectors_list):
     
    # Populate arrays based on neighbor queries
    for j in prange(N):
        neighbors = neighbors_initial[j, :]

        vectors_to_neighbors = prt_rep[neighbors] - particles[j]

        if periodic_box:
            neighbors = prtidx[neighbors]

        neighbor_list[j, :len(neighbors)] = neighbors
        vectors_list[j, :len(neighbors), :] = vectors_to_neighbors
    return h,neighbor_list,distances_list,vectors_list



def calculate_smoothing_length_staged(treeIN, particles, box_size, periodic_box, final_k=2):
    N = len(particles)
    # Adjust initial_k and final_k
    final_k = max(final_k - 1, 1)
    # Unpack tree information
    tree, prtidx, prt_rep, tree_type = treeIN
    
    # Initialize output arrays
    h = np.zeros(N)
    neighbor_list = -np.ones((N, final_k + 1), dtype=np.int64)
    distances_list = np.zeros((N, final_k + 1), dtype=PRECISION)
    vectors_list = np.zeros((N, final_k + 1, 3),dtype=PRECISION)
    
    # Perform tree query for initial distances and neighbors
    distances_list, neighbors_list,h = mytree_query2(particles, final_k + 1, tree)
    
    distances_list = np.array(distances_list,dtype=PRECISION)
    neighbors_list = np.array(neighbors_list,dtype=np.int64)
  
    h,neighbors_list,distances_list,vectors_list = sort_neighbours(N,particles,prt_rep,prtidx,neighbors_list,distances_list,periodic_box,final_k,h,neighbor_list,distances_list,vectors_list)
    
    # Return computed smoothing lengths and neighbor information
    return h, neighbors_list, distances_list, vectors_list



@njit(parallel=True)
def calculate_density(particles, masses, h, neighbors_list, distances_list, box_size, periodic_box):
    N = len(particles)
    density = np.zeros(N)
    
    for i in prange(N):  # Use parallel range for multi-threading
        neighbors = neighbors_list[i]
        distances = distances_list[i]  # Get the precomputed distances for this particle
        for ni, j in enumerate(neighbors):
            if j == -1:  # Skip invalid padded entries
                continue
            
            # Use the precomputed distance from the smoothing function
            r = distances[ni]  # Distance already includes periodic adjustments (if periodic_box was enabled)
            
            density[i] += masses[j] * kernel(r,h[i])  # Include mass in density calculation
    
    return density

@njit(parallel=True)
def EOS(particles, gamma=5./3.):
    N = len(particles)
    for i in prange(N):
        rho = particles[i, 7]  # density
        u = particles[i, 8]  # either T or u, depending on your setup
        # Adiabatic EoS: P = rho * x_val * (gamma - 1)
        P = rho * u * (gamma - 1.0)
        # Store in column 10
        particles[i, 11] = P
    
    return particles


@njit(parallel=True)
def calculate_gradients(particles, neighbors_list, distances_list, vectors_list, box_size, periodic_box):
    N = len(particles)
    m = particles[:, 0]
    h = particles[:, 20]
    density = particles[:, 7]
    P = particles[:, 11]
    vx = particles[:,4]
    vy = particles[:,5]
    vz = particles[:,6]
    Poverrho=particles[:, 11]/particles[:, 7] 
    c_sound = np.sqrt((5./3.)*Poverrho);

    Q0 = np.zeros(N)
    Q1x = np.zeros(N)
    Q1y = np.zeros(N)
    Q1z = np.zeros(N)
    
    E0x = np.zeros(N)
    E0y = np.zeros(N)
    E0z = np.zeros(N)
    
    dvxdt = np.zeros(N)
    dvydt = np.zeros(N)
    dvzdt = np.zeros(N)
    dudt = np.zeros(N)

    num_threads = get_num_threads()  # Dynamically get available threads

    # Buffer sized by number of threads to prevent race conditions
    buffer_dvxdt = np.zeros((N, num_threads))
    buffer_dvydt = np.zeros((N, num_threads))
    buffer_dvzdt = np.zeros((N, num_threads))
    buffer_dudt = np.zeros((N, num_threads))

    for i in prange(N):
        neighbors = neighbors_list[i]
        distances = distances_list[i]
        vectors = vectors_list[i]

        for ni, j in enumerate(neighbors):
            if j == -1:
                continue
            
            r = distances[ni]
            dx, dy, dz = vectors[ni]
            
            

            vx_temp=vx[j]
            vy_temp=vy[j]
            vz_temp=vz[j]
            

            test_fac=1.0
            
            # for reflecting
            # vx_temp=apply_1d_reflecting_vel(x[i]+dx,vx_temp, box_size[0])
            # vy_temp=apply_1d_reflecting_vel(y[i]+dy,vy_temp, box_size[1])
            # vz_temp=apply_1d_reflecting_vel(z[i]+dz,vz_temp, box_size[2])
            
            # 
            # if vx_temp != vx[j] or vy_temp != vy[j] or vz_temp != vz[j]:
            #     test_fac=-1.0


            
            dvx = vx_temp-vx[i]
            dvy = vy_temp-vy[i]
            dvz = vz_temp-vz[i]
            

            
            # Compute the dot product of velocity difference and position vector
            dvdotdr = dvx * dx + dvy * dy + dvz * dz

            # Artificial viscosity only acts when particles are approaching
            havg=0.5*0.5*(h[i]+h[j])
            if dvdotdr < 0:
                # Compute the viscosity coefficient
                absmu = -(havg * dvdotdr) / (r**2 + ETA * h[i]**2)
                # Monaghan's artificial viscosity
                visc = ALPHA*(c_sound[i]+c_sound[j]) + BETA * 2 * absmu
                visc = visc*absmu/(density[i]+density[j])
            else:
                visc = 0.0
                
            
            W = kernel(r, h[i])
            dW = Dkernel(r, h[i])
            
            Va = m[i] / density[i]
            Vb = m[j] / density[j]
            Pi = P[i]
            Pj = P[j]

            Q0[i] += Vb * W
            Q1x[i] += Vb * dx * W
            Q1y[i] += Vb * dy * W
            Q1z[i] += Vb * dz * W
            E0x[i]     += 2.0 * Vb * dx * dW
            E0y[i]     += 2.0 * Vb * dy * dW
            E0z[i]     += 2.0 * Vb * dz * dW

            dvxdt[i] += 0.5 * Vb * (Pi + Pj) * dW * (dx / density[i])
            dvydt[i] += 0.5 * Vb * (Pi + Pj) * dW * (dy / density[i])
            dvzdt[i] += 0.5 * Vb * (Pi + Pj) * dW * (dz / density[i])
            dudt[i] += 0.5 * Vb * (Pi)/density[i] * dvdotdr * dW
            
            dvxdt[i] += 0.5 * m[j] * (visc+visc) * dW * dx
            dvydt[i] += 0.5 * m[j] * (visc+visc) * dW * dy
            dvzdt[i] += 0.5 * m[j] * (visc+visc) * dW * dz
            dudt[i] += 0.5 * m[j] * visc * dvdotdr * dW
            

            tid = get_thread_id()  # Thread ID
            buffer_dvxdt[j, tid] -= test_fac*0.5 * Va * (Pi + Pj) * dW * (dx / density[j])
            buffer_dvydt[j, tid] -= test_fac*0.5 * Va * (Pi + Pj) * dW * (dy / density[j])
            buffer_dvzdt[j, tid] -= test_fac*0.5 * Va * (Pi + Pj) * dW * (dz / density[j])
            buffer_dudt[j, tid] += 0.5 * Va * (Pj)/density[j] * dvdotdr * dW
            
            buffer_dvxdt[j, tid] -= test_fac*0.5 * m[i] * (visc+visc) * dW * dx
            buffer_dvydt[j, tid] -= test_fac*0.5 * m[i] * (visc+visc) * dW * dy
            buffer_dvzdt[j, tid] -= test_fac*0.5 * m[i] * (visc+visc) * dW * dz
            buffer_dudt[j, tid] += 0.5 * m[i] * visc * dvdotdr * dW
    
    # Apply buffer to avoid race conditions
    for j in prange(N):
        dvxdt[j] += np.sum(buffer_dvxdt[j])
        dvydt[j] += np.sum(buffer_dvydt[j])
        dvzdt[j] += np.sum(buffer_dvzdt[j])
        dudt[j] += np.sum(buffer_dudt[j])

    return Q0,Q1x,Q1y,Q1z,E0x,E0y,E0z,dvxdt,dvydt,dvzdt,dudt

@njit
def calculate_angular_momentum_cm(pos1, pos2, vel1, vel2, mass1, mass2):
    # Calculate the center of mass
    center_of_mass = (mass1 * pos1 + mass2 * pos2) / (mass1 + mass2)

    # Position relative to the center of mass
    rel_pos1 = pos1 - center_of_mass
    rel_pos2 = pos2 - center_of_mass

    # Angular momentum
    L1 = np.cross(rel_pos1, mass1 * vel1)
    L2 = np.cross(rel_pos2, mass2 * vel2)

    # Total angular momentum
    total_L = L1 + L2

    return total_L

@njit
def merge_my_neighbor(unmerged_particles,merged_particles,nearest_neighbors,nearest_distance,vectors_to_neighbors,h_list,box_size,periodic_box,merge_condition=None):
    N = len(unmerged_particles)
    merge_target = -np.ones(N, dtype=np.int64)  # Initialize with -1 (no target)
    
    unmerged_masses=unmerged_particles[:,0]
    unmerged_particles_pos=unmerged_particles[:,1:4]
    unmerged_particles_vel=unmerged_particles[:,4:7]
    unmerged_density=unmerged_particles[:,7]
    unmerged_h=unmerged_particles[:,20]
    unmerged_u=unmerged_particles[:,8]
    unmerged_soft=unmerged_particles[:,9]
    unmerged_metals=unmerged_particles[:,10]
    unmerged_B=unmerged_particles[:,12:15]
    unmerged_spin=unmerged_particles[:,15:18]
    unmerged_momenti=unmerged_particles[:,18]
    unmerged_tform=unmerged_particles[:,19]
    unmerged_volume=unmerged_masses/unmerged_density
    
    merged = np.full(N, False)
    closest_distance = np.full(N, np.inf)
    closest_vectors = np.zeros((N, 3))  # Store displacement vectors

    nomergethisround = True
    alone = np.zeros(N, dtype=np.int64)
    alone_count = 0
    
    still_unmerged = np.zeros((N, unmerged_particles.shape[1]))
    still_unmerged_h = np.zeros(N)
    still_unmerged_count = 0
    
    merged_particles_array = np.zeros((N // 2, PARTICLE_ARRAY_SIZE))
    merged_count = 0
    
    
    # Process the nearest neighbor pairs
    for i in range(N):
        if merge_condition is not None and not merge_condition[i]:
            continue  # Skip particles not marked for merging
        k = 1  # Initialize the index for neighbors
        while k < len(nearest_neighbors[i, :]):
            nnd = nearest_distance[i, k]
            nnv = vectors_to_neighbors[i, k]
            nn = nearest_neighbors[i, k]
            ## if nearest is  further away then smoothing length do not merge
            if nnd > ((h_list[i]+h_list[nn])):
                alone[alone_count] = i
                alone_count += 1
                break
            if nnd <= closest_distance[i]:
                # Check if closer particle is found first
                if nnd <= closest_distance[nn]:                    
                    if merge_target[i] > -1:
                         merged[merge_target[i]] = False
                         merge_target[merge_target[i]] = -1
                         closest_distance[merge_target[i]] = np.inf
                    if merge_target[nn] > -1:
                         merged[merge_target[nn]] = False
                         merge_target[merge_target[nn]] = -1
                         closest_distance[merge_target[nn]] = np.inf
                         
                    # These are the closest neighbors at the moment
                    closest_distance[i] = nnd
                    closest_distance[nn] = nnd
                    closest_vectors[nn] = -nnv
                    closest_vectors[i] = nnv
                    merge_target[nn] = i
                    merge_target[i] = nn
                    merged[i] = True
                    merged[nn] = True
                    nomergethisround=False
                    break
            k += 1  # Move to the next neighbor
            

    # Merging particles based on targets
    for i in range(N):
        if merged[i]:
            nn = merge_target[i]
            if nn == -1:
                continue  # No valid merge target
            
            # Sum masses after merging
            new_mass = unmerged_masses[i] + unmerged_masses[nn]
            # Use the displacement vector to place the new particle at the midpoint
            new_particle_pos = (unmerged_masses[i] * unmerged_particles_pos[i] + 
                    unmerged_masses[nn] * (unmerged_particles_pos[i] + closest_vectors[i])) / new_mass

            # Wrap the new particle within the periodic box if necessary
            if periodic_box:
                new_particle_pos,_ = apply_boundary_condition(new_particle_pos, box_size)

            
            # Merge other quantities
            # Start with velocity, set to conserve momentum
            new_particle_vel = (unmerged_masses[i]*unmerged_particles_vel[i] + unmerged_masses[nn]*unmerged_particles_vel[nn])/new_mass
            # Calculate the difference in kinetic energy
            Ekin_diff = 0.5 * (unmerged_masses[i] * np.sum(unmerged_particles_vel[i]**2) +
                   unmerged_masses[nn] * np.sum(unmerged_particles_vel[nn]**2)) - 0.5 * new_mass * np.sum(new_particle_vel**2)
            # Save angular momentum in spin
            # Calculate the angular momentum relative to the center of mass
            #new_spin = new angular momentum
            LM=calculate_angular_momentum_cm(unmerged_particles_pos[i], unmerged_particles_pos[nn],
                                            unmerged_particles_vel[i], unmerged_particles_vel[nn],
                                            unmerged_masses[i], unmerged_masses[nn])
            
            new_spin = (LM + unmerged_spin[i]*unmerged_masses[i] + unmerged_spin[nn]*unmerged_masses[nn])/new_mass
            
            if unmerged_momenti[i] > 0.0:
                Ekin_diff += 0.5*np.sum(unmerged_spin[i]**2)/unmerged_momenti[i]
            if unmerged_momenti[nn] > 0.0:
                Ekin_diff += 0.5*np.sum(unmerged_spin[nn]**2)/unmerged_momenti[nn]
                
            if Ekin_diff > 0:
                new_momenti = np.linalg.norm(new_spin)**2 / (2 * Ekin_diff)
            else:
                new_momenti = 0  # Or set to an alternative value if zero or negative is invalid
            
            # Calculate BField of merged particle (Here we choose to be flux conservative, 
            #another way is to be energy conservative, though we can add the missing mag energy to thermal)
            #Here we add the total particle BFlux to the new particle, which can then be normalized after the merger with new particle volume.
            #As we do not know the particle volume yet, we have to also distribute any leftover magnetic energy later.
            #Could this perhaps be put in the cleaning field btw ? 
            #new_particle_BV =unmerged_B[i]*unmerged_volume[i] + unmerged_B[nn]*unmerged_volume[nn]
            new_particle_BV = 0.5 * (unmerged_B[i] + unmerged_B[nn])
            
            # Distribute leftover energy to thermal and average the thermal energy
            u_ex=0.0
            new_particle_u=unmerged_u[i]+unmerged_u[nn]+u_ex
            new_particle_soft=0.5*(unmerged_soft[i]+unmerged_soft[nn])
            new_particle_metals=0.5*(unmerged_metals[i]+unmerged_metals[nn])
            
            # Create new particle array with all merged quantities
            merged_particles_array[merged_count, 0] = new_mass
            merged_particles_array[merged_count, 1:4] = new_particle_pos
            merged_particles_array[merged_count, 4:7] = new_particle_vel
            merged_particles_array[merged_count, 7] = 0.0
            merged_particles_array[merged_count, 8] = new_particle_u
            merged_particles_array[merged_count, 9] = new_particle_soft
            merged_particles_array[merged_count, 10] = new_particle_metals
            
            merged_particles_array[merged_count, 12:15] = new_particle_BV
            merged_particles_array[merged_count, 15:18] = new_spin
            merged_particles_array[merged_count, 18] = new_momenti
            merged_particles_array[merged_count, 19] = (unmerged_tform[i]+unmerged_tform[nn])*0.5
            merged_count += 1
            
            # Remove merge_target from neighbour to avoid double merge
            merge_target[nn] = -1
            
        else:
            # Add unmerged particles to the list for the next round
            still_unmerged[still_unmerged_count] = unmerged_particles[i]
            still_unmerged_h[still_unmerged_count] = h_list[i]
            still_unmerged_count += 1


    # Convert lists to arrays for the next round
    unmerged_particles = still_unmerged[:still_unmerged_count]
    h_list = still_unmerged_h[:still_unmerged_count]
    if nomergethisround==False:
        alone=alone[:alone_count]

    return unmerged_particles,h_list,merged_particles_array[:merged_count],alone
    


def map_all_from_condition(particles,condition):
    """Applies a condition to filter multiple arrays based on the condition applied to h."""
    if condition==None:
        mask = particles[:, 0] > -1.0;
    else:
        mask = condition(particles)
    return mask

# Function to merge particles using the strategy described
def merge_particles(particles, treeIN, neighbors_list, distances_list,dvec_list, box_size, periodic_box=False,condition=None):
    h_list=particles[:,20]
    unmerged_particles = particles
    #unmerged_masses = particles[0]
    iteration = 0  # To keep track of iterations
    merged_particles = np.empty((0, PARTICLE_ARRAY_SIZE))
    alone = set()
    tree=treeIN[0]
    prtidx=treeIN[1]
    prt_rep=treeIN[2]
    tree_type=treeIN[3]
    
    

    while len(unmerged_particles)-len(alone) > 1:  # Continue merging until 1 or no particle is left unmerged
        N = len(unmerged_particles)
        print(N,len(alone))
        if(iteration==0):
            distances=distances_list
            indices=neighbors_list
            dvec = dvec_list
        else:
            tree_temp = build_tree(unmerged_particles[:,1:4], box_size, 2*h_list , periodic_box)
            tree=tree_temp[0]
            prtidx=tree_temp[1]
            h_temp, indices, distances, dvec = calculate_smoothing_length_staged(
                tree_temp, unmerged_particles[:,1:4],box_size=box_size, periodic_box=periodic_box, final_k=2)
        
        nearest_neighbors = prtidx[indices[:, :]]
        nearest_distance = distances[:, :]
        vectors_to_neighbors = dvec[:, :]
        
        
        #only include neighbours of merge_condition particles in unmerged_particles.
        merge_condition=map_all_from_condition(unmerged_particles,condition)
        
        
        unmerged_particles,h_list,merged_particles_temp,alone=merge_my_neighbor(unmerged_particles,merged_particles,nearest_neighbors,nearest_distance,vectors_to_neighbors,h_list,box_size,periodic_box,merge_condition=merge_condition)
        # Add merged_particles_temp to merged_particles
        merged_particles = np.vstack((merged_particles, merged_particles_temp))
        
        # Increment the iteration counter
        iteration += 1
    
    return np.array(merged_particles), np.array(unmerged_particles)

@njit(parallel=True)
def addunmerged(merged_particles, unmerged_particles, method="add"):
    """
    This function adds unmerged particles to merged particles using one of two methods:
    1. 'normalize': Transform unmerged particles into merged particles and normalize so that the mass of all particles is the same.
    2. 'add': Simply add the unmerged particles to the merged array.
    
    """
    # Extract the mass column (column 0) from both merged and unmerged particles
    merged_masses = merged_particles[:, 0]
    unmerged_masses = unmerged_particles[:, 0]

    if method == "normalize":
        # Method 1: Normalize the mass of all particles
        total_mass = np.sum(merged_masses) + np.sum(unmerged_masses)
        total_particles = len(merged_particles) + len(unmerged_particles)
        normalized_mass = total_mass / total_particles
        
        # Add unmerged particles to merged_particles
        merged_particles = np.concatenate((merged_particles, unmerged_particles))
        
        # Normalize the mass of all particles to the new normalized mass
        merged_particles[:, 0] = normalized_mass  # Set all masses to the normalized mass

    elif method == "add":
        # Method 2: Simply add unmerged particles to merged particles
        merged_particles = np.concatenate((merged_particles, unmerged_particles))
        
    else:
        raise ValueError("Invalid method. Choose 'normalize' or 'add'.")
    
    return merged_particles
