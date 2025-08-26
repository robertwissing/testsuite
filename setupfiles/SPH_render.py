#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 11:12:42 2024

@author: robertwi
"""


import matplotlib.pyplot as plt
import numpy as np


def plot_particle_scatter_and_projections(*particle_datasets, labels=None, colors=None, marker_sizes=None, sfac=7.0, unmerged_particles=None):
    """
    Plots 3D scatter plots and 2D projections (XY, XZ, YZ planes) for multiple particle datasets.

    Parameters:
    - particle_datasets: variable number of ndarrays, each containing particle data.
    - labels: list of strings, labels for each dataset.
    - colors: list of colors for each dataset.
    - marker_sizes: list of floats, scaling factors for marker sizes of each dataset.
    - sfac: float, base scatter size factor (default is 7.0).
    - unmerged_particles: ndarray or None, data for unmerged particles (optional).
    """
    num_datasets = len(particle_datasets)
    N = max(len(particles) for particles in particle_datasets)

    if N < 1000:
        fig = plt.figure(figsize=(16, 12))
        ax_3d = fig.add_subplot(221, projection='3d')

        # Default labels, colors, and marker sizes if none provided
        if labels is None:
            labels = [f'Dataset {i+1}' for i in range(num_datasets)]
        if colors is None:
            colors = ['green', 'red', 'blue', 'orange', 'purple', 'cyan', 'magenta']
        if marker_sizes is None:
            marker_sizes = [sfac * (i+1) for i in range(num_datasets)]

        # Plot each dataset in 3D
        for idx, particles in enumerate(particle_datasets):
            label = labels[idx]
            color = colors[idx % len(colors)]
            size = marker_sizes[idx % len(marker_sizes)]
            ax_3d.scatter(particles[:, 1], particles[:, 2], particles[:, 3],
                          s=size, c=color, marker='o', label=label)

        # Plot unmerged particles if provided
        if unmerged_particles is not None and unmerged_particles.size > 0 and unmerged_particles.ndim > 1:
            ax_3d.scatter(unmerged_particles[:, 1], unmerged_particles[:, 2], unmerged_particles[:, 3],
                          s=sfac * 5, c='orange', marker='o', label="Un-Merged Particles")

        ax_3d.set_title("3D Scatter Plot of Particle Positions")
        ax_3d.set_xlabel('X')
        ax_3d.set_ylabel('Y')
        ax_3d.set_zlabel('Z')
        ax_3d.legend()

        # 2D Projections
        # XY Plane
        ax_xy = fig.add_subplot(222)
        for idx, particles in enumerate(particle_datasets):
            label = labels[idx]
            color = colors[idx % len(colors)]
            size = marker_sizes[idx % len(marker_sizes)]
            ax_xy.scatter(particles[:, 2], particles[:, 1], s=size, c=color, marker='o', label=label)
        if unmerged_particles is not None and unmerged_particles.size > 0 and unmerged_particles.ndim > 1:
            ax_xy.scatter(unmerged_particles[:, 2], unmerged_particles[:, 1],
                          s=sfac * 5, c='orange', marker='o', label="Un-Merged Particles")
        ax_xy.set_title("2D Projection: XY Plane")
        ax_xy.set_xlabel('Y')
        ax_xy.set_ylabel('X')
        ax_xy.legend()

        # XZ Plane
        ax_xz = fig.add_subplot(223)
        for idx, particles in enumerate(particle_datasets):
            label = labels[idx]
            color = colors[idx % len(colors)]
            size = marker_sizes[idx % len(marker_sizes)]
            ax_xz.scatter(particles[:, 1], particles[:, 3], s=size, c=color, marker='o', label=label)
        if unmerged_particles is not None and unmerged_particles.size > 0 and unmerged_particles.ndim > 1:
            ax_xz.scatter(unmerged_particles[:, 1], unmerged_particles[:, 3],
                          s=sfac * 5, c='orange', marker='o', label="Un-Merged Particles")
        ax_xz.set_title("2D Projection: XZ Plane")
        ax_xz.set_xlabel('X')
        ax_xz.set_ylabel('Z')
        ax_xz.legend()

        # YZ Plane
        ax_yz = fig.add_subplot(224)
        for idx, particles in enumerate(particle_datasets):
            label = labels[idx]
            color = colors[idx % len(colors)]
            size = marker_sizes[idx % len(marker_sizes)]
            ax_yz.scatter(particles[:, 2], particles[:, 3], s=size, c=color, marker='o', label=label)
        if unmerged_particles is not None and unmerged_particles.size > 0 and unmerged_particles.ndim > 1:
            ax_yz.scatter(unmerged_particles[:, 2], unmerged_particles[:, 3],
                          s=sfac * 5, c='orange', marker='o', label="Un-Merged Particles")
        ax_yz.set_title("2D Projection: YZ Plane")
        ax_yz.set_xlabel('Y')
        ax_yz.set_ylabel('Z')
        ax_yz.legend()

        plt.tight_layout()
        plt.show()


def plot_density_scatter_projections(*particle_datasets, labels=None, colors=None, sfac=0.5,x_min=None, x_max=None, y_min=None, y_max=None, z_min=None, z_max=None):
    """
    Plots 2D scatter projections (XY and XZ planes) of particle datasets,
    colored by log density and sized by smoothing length.

    Parameters:
    - particle_datasets: variable number of ndarrays, each containing particle data.
    - labels: list of strings, labels for each dataset.
    - colors: list of colors for each dataset's markers.
    - sfac: float, scaling factor for marker sizes.
    """
    num_datasets = len(particle_datasets)
    fig, ax = plt.subplots(2, num_datasets, figsize=(6 * num_datasets, 12))

    # Default labels and colors if none provided
    if labels is None:
        labels = [f'Dataset {i+1}' for i in range(num_datasets)]
    if colors is None:
        colors = [plt.cm.tab10(i) for i in range(num_datasets)]

    # Compute global min and max densities for consistent color scaling
    all_log_densities = []
    all_positions = []
    for particles in particle_datasets:
        all_log_densities.append(np.log10(particles[:, 9]))
        all_positions.append(particles[:, 1:4])  # X, Y, Z positions

    min_density = np.min(np.hstack(all_log_densities))
    max_density = np.max(np.hstack(all_log_densities))

    # Compute global bounds for the axes
    all_positions = np.vstack(all_positions)
    x_min = x_min if x_min is not None else np.min(all_positions[:, 0])
    x_max = x_max if x_max is not None else np.max(all_positions[:, 0])
    y_min = y_min if y_min is not None else np.min(all_positions[:, 1])
    y_max = y_max if y_max is not None else np.max(all_positions[:, 1])
    z_min = z_min if z_min is not None else np.min(all_positions[:, 2])
    z_max = z_max if z_max is not None else np.max(all_positions[:, 2])

    for idx, particles in enumerate(particle_datasets):
        label = labels[idx]
        color = colors[idx % len(colors)]
        log_density = np.log10(particles[:, 9])
        marker_sizes = particles[:, 8] * 50 * sfac  # Adjust marker size by smoothing length

        # XY Plane
        scatter_xy = ax[0, idx].scatter(
            particles[:, 1], particles[:, 2],
            c=log_density, s=marker_sizes,
            cmap='plasma', vmin=min_density, vmax=max_density
        )
        ax[0, idx].set_title(f'{label} (XY plane)')
        ax[0, idx].set_xlabel('X')
        ax[0, idx].set_ylabel('Y')
        ax[0, idx].set_xlim(x_min, x_max)
        ax[0, idx].set_ylim(y_min, y_max)
        fig.colorbar(scatter_xy, ax=ax[0, idx], label='Log Density')

        # XZ Plane
        scatter_xz = ax[1, idx].scatter(
            particles[:, 1], particles[:, 3],
            c=log_density, s=marker_sizes,
            cmap='plasma', vmin=min_density, vmax=max_density
        )
        ax[1, idx].set_title(f'{label} (XZ plane)')
        ax[1, idx].set_xlabel('X')
        ax[1, idx].set_ylabel('Z')
        ax[1, idx].set_xlim(x_min, x_max)
        ax[1, idx].set_ylim(z_min, z_max)
        fig.colorbar(scatter_xz, ax=ax[1, idx], label='Log Density')

    plt.tight_layout()
    plt.show()

def plot_density_vs_smoothing_length(*particle_datasets, 
                                     labels=None, 
                                     colors=None, 
                                     density_index=7, 
                                     smoothing_length_index=8, 
                                     log_x=True, 
                                     log_y=True):
    """
    Plots Smoothing Length vs Density for multiple particle datasets.

    Parameters:
    - particle_datasets: variable number of ndarrays, each containing particle data.
    - labels: list of strings, labels for each dataset.
    - colors: list of colors for each dataset.
    - density_index: int, index of the density field in the dataset.
    - smoothing_length_index: int, index of the smoothing length field in the dataset.
    - log_x: bool, whether to apply logarithm to the x-axis (density).
    - log_y: bool, whether to apply logarithm to the y-axis (smoothing length).
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    # Default colors if none provided
    if colors is None:
        colors = plt.cm.get_cmap('tab10').colors

    for idx, particles in enumerate(particle_datasets):
        x_data = np.log10(particles[:, density_index]) if log_x else particles[:, density_index]
        y_data = np.log10(particles[:, smoothing_length_index]) if log_y else particles[:, smoothing_length_index]
        label = labels[idx] if labels and idx < len(labels) else f"Dataset {idx+1}"
        color = colors[idx % len(colors)]
        ax.scatter(x_data, y_data, s=1, c=[color], alpha=0.7, label=label)

    ax.set_title('Smoothing Length vs Density')
    x_label = 'Log Density' if log_x else 'Density'
    y_label = 'Log Smoothing Length' if log_y else 'Smoothing Length'
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.legend()
    plt.show()


def plot_radial_distribution(*particle_datasets, labels=None, colors=None):
    """
    Plots Radial Distance vs Smoothing Length for multiple particle datasets.

    Parameters:
    - particle_datasets: variable number of ndarrays, each containing particle data.
    - labels: list of strings, labels for each dataset.
    - colors: list of colors for each dataset.
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    # Default colors if none provided
    if colors is None:
        colors = plt.cm.get_cmap('tab10').colors

    for idx, particles in enumerate(particle_datasets):
        radial_distance = np.sqrt(np.sum(particles[:, 1:4] ** 2, axis=1))
        log_smoothing_length = np.log10(particles[:, 8])
        label = labels[idx] if labels and idx < len(labels) else f"Dataset {idx+1}"
        color = colors[idx % len(colors)]
        ax.scatter(radial_distance, log_smoothing_length, s=50, c=[color], alpha=0.7, label=label)

    ax.set_title('Radial Distance vs Log Smoothing Length')
    ax.set_xlabel('Radial Distance')
    ax.set_ylabel('Log Smoothing Length')
    ax.legend()
    plt.show()