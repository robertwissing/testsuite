import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.io as pio
import numpy as np
import matplotlib.pyplot as plt

def plot_slice(testgrid, x_grid, y_grid, z_grid, cell_levels, axis='z', index=None, cmap='rainbow'):
    """
    Plots a 2D slice of testgrid along the specified axis at the given index, or at z=0 by default, with AMR grid overlay.
    
    Parameters:
    - testgrid: 3D numpy array of interpolated data.
    - x_grid, y_grid, z_grid: 1D arrays specifying the grid coordinates in each dimension.
    - cell_levels: Dictionary with (i, j, k) tuples as keys and refinement levels as values.
    - axis: Axis along which to take the slice ('x', 'y', or 'z').
    - index: Index of the slice along the specified axis. If None, the closest slice to zero is used for 'z'.
    - cmap: Colormap for the plot.
    """
    # Set index to the closest slice at z=0 if axis is 'z'
    if axis == 'z':
        if index is None:
            index = np.where(np.isclose(z_grid, 0))[0][0]  # Find index where z=0
        data_slice = testgrid[:, :, index]
        X, Y = np.meshgrid(x_grid, y_grid, indexing='ij')
        slice_coord = z_grid[index]
        xlabel, ylabel = 'x', 'y'
        plot_cells = [(i, j, level) for (i, j, k), level in cell_levels.items() if k == index]
    elif axis == 'y':
        if index is None:
            index = np.where(np.isclose(y_grid, 0))[0][0]  # Find index where y=0
        data_slice = testgrid[:, index, :]
        X, Y = np.meshgrid(x_grid, z_grid, indexing='ij')
        slice_coord = y_grid[index]
        xlabel, ylabel = 'x', 'z'
        plot_cells = [(i, k, level) for (i, j, k), level in cell_levels.items() if j == index]
    elif axis == 'x':
        if index is None:
            index = np.where(np.isclose(x_grid, 0))[0][0]  # Find index where x=0
        data_slice = testgrid[index, :, :]
        X, Y = np.meshgrid(y_grid, z_grid, indexing='ij')
        slice_coord = x_grid[index]
        xlabel, ylabel = 'y', 'z'
        plot_cells = [(j, k, level) for (i, j, k), level in cell_levels.items() if i == index]
    else:
        raise ValueError("Axis must be 'x', 'y', or 'z'")

    # Plot the data slice
    plt.figure(figsize=(8, 6))
    plt.pcolormesh(X, Y, data_slice.T, shading='auto', cmap=cmap)
    plt.colorbar(label='Interpolated Data')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(f'Slice at {axis} = {slice_coord:.2f}')

    # Overlay AMR grid boundaries with color representing refinement level
    for (i, j, level) in plot_cells:
        # Get boundaries for each cell in the AMR grid
        if axis == 'z':
            x0, x1 = x_grid[i], x_grid[i + 1]
            y0, y1 = y_grid[j], y_grid[j + 1]
        elif axis == 'y':
            x0, x1 = x_grid[i], x_grid[i + 1]
            y0, y1 = z_grid[j], z_grid[j + 1]
        elif axis == 'x':
            x0, x1 = y_grid[i], y_grid[i + 1]
            y0, y1 = z_grid[j], z_grid[j + 1]

        # Adjust color or line thickness based on the refinement level
        plt.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], 
                 color='black', linewidth=0.5 + 0.2 * level)  # Adjust line width by level

    plt.show()

def plot_multiple_slices(testgrid, box_size, axis='z', indices=None, cmap='viridis'):
    """
    Plots multiple 2D slices of testgrid along the specified axis.

    Parameters:
    - testgrid: 3D numpy array of interpolated data.
    - box_size: List or array of [size_x, size_y, size_z].
    - axis: Axis along which to take the slices ('x', 'y', or 'z').
    - indices: List of indices of the slices. If None, three equally spaced slices are used.
    - cmap: Colormap for the plots.
    """
    import matplotlib.pyplot as plt
    npixx, npixy, npixz = testgrid.shape
    xmin, ymin, zmin = -box_size[0]/2, -box_size[1]/2, -box_size[2]/2
    xmax, ymax, zmax = box_size[0]/2, box_size[1]/2, box_size[2]/2

    # Create coordinate grids
    x_grid = np.linspace(xmin, xmax, npixx)
    y_grid = np.linspace(ymin, ymax, npixy)
    z_grid = np.linspace(zmin, zmax, npixz)

    if indices is None:
        if axis == 'x':
            indices = [npixx // 4, npixx // 2, 3 * npixx // 4]
        elif axis == 'y':
            indices = [npixy // 4, npixy // 2, 3 * npixy // 4]
        elif axis == 'z':
            indices = [npixz // 4, npixz // 2, 3 * npixz // 4]
        else:
            raise ValueError("Axis must be 'x', 'y', or 'z'")

    n_slices = len(indices)
    fig, axes = plt.subplots(1, n_slices, figsize=(6 * n_slices, 6))

    # Ensure axes is iterable even if there's only one subplot
    if n_slices == 1:
        axes = [axes]

    for ax, index in zip(axes, indices):
        if axis == 'z':
            data_slice = testgrid[:, :, index]
            X, Y = np.meshgrid(x_grid, y_grid, indexing='ij')
            slice_coord = z_grid[index]
            xlabel, ylabel = 'x', 'y'
            data_slice_to_plot = data_slice.T  # Transpose needed
        elif axis == 'y':
            data_slice = testgrid[:, index, :]
            X, Y = np.meshgrid(x_grid, z_grid, indexing='ij')
            slice_coord = y_grid[index]
            xlabel, ylabel = 'x', 'z'
            data_slice_to_plot = data_slice  # No transpose
        elif axis == 'x':
            data_slice = testgrid[index, :, :]
            X, Y = np.meshgrid(y_grid, z_grid, indexing='ij')
            slice_coord = x_grid[index]
            xlabel, ylabel = 'y', 'z'
            data_slice_to_plot = data_slice.T  # Transpose needed
        else:
            raise ValueError("Axis must be 'x', 'y', or 'z'")

        # Check if shapes match
        if data_slice_to_plot.shape != X.shape:
            print(f"Mismatch in shapes: data_slice {data_slice_to_plot.shape}, X {X.shape}")
        
        im = ax.pcolormesh(X, Y, data_slice_to_plot, shading='auto', cmap=cmap)
        ax.set_title(f'{axis} = {slice_coord:.2f}')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        fig.colorbar(im, ax=ax, label='Interpolated Data')

    plt.tight_layout()
    plt.show()


def plot_projection(testgrid, box_size, axis='z', method='sum', cmap='viridis'):
    """
    Creates and plots a 2D projection of testgrid by summing or averaging along the specified axis.
    
    Parameters:
    - testgrid: 3D numpy array of interpolated data.
    - box_size: List or array of [size_x, size_y, size_z].
    - axis: Axis along which to project ('x', 'y', or 'z').
    - method: 'sum' or 'mean' to define the projection method.
    - cmap: Colormap for the plot.
    """
    npixx, npixy, npixz = testgrid.shape
    xmin, ymin, zmin = -box_size[0]/2, -box_size[1]/2, -box_size[2]/2
    xmax, ymax, zmax = box_size[0]/2, box_size[1]/2, box_size[2]/2
    
    # Create coordinate grids
    x_grid = np.linspace(xmin, xmax, npixx)
    y_grid = np.linspace(ymin, ymax, npixy)
    z_grid = np.linspace(zmin, zmax, npixz)
    
    if axis == 'z':
        if method == 'sum':
            projection = np.sum(testgrid, axis=2)
        elif method == 'mean':
            projection = np.mean(testgrid, axis=2)
        X, Y = np.meshgrid(x_grid, y_grid, indexing='ij')
        xlabel, ylabel = 'x', 'y'
    elif axis == 'y':
        if method == 'sum':
            projection = np.sum(testgrid, axis=1)
        elif method == 'mean':
            projection = np.mean(testgrid, axis=1)
        X, Y = np.meshgrid(x_grid, z_grid, indexing='ij')
        xlabel, ylabel = 'x', 'z'
    elif axis == 'x':
        if method == 'sum':
            projection = np.sum(testgrid, axis=0)
        elif method == 'mean':
            projection = np.mean(testgrid, axis=0)
        X, Y = np.meshgrid(y_grid, z_grid, indexing='ij')
        xlabel, ylabel = 'y', 'z'
    else:
        raise ValueError("Axis must be 'x', 'y', or 'z'")
    
    plt.figure(figsize=(8, 6))
    plt.pcolormesh(X, Y, projection.T, shading='auto', cmap=cmap)
    plt.colorbar(label='Projected Data')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(f'Projection along {axis}-axis ({method})')
    plt.show()

import plotly.graph_objects as go
import plotly.io as pio

def plot_3d_volume(testgrid, box_size):
    """
    Visualizes testgrid using 3D volume rendering with Plotly.
    
    Parameters:
    - testgrid: 3D numpy array of interpolated data.
    - box_size: List or array of [size_x, size_y, size_z].
    """
    npixx, npixy, npixz = testgrid.shape
    xmin, ymin, zmin = -box_size[0]/2, -box_size[1]/2, -box_size[2]/2
    xmax, ymax, zmax = box_size[0]/2, box_size[1]/2, box_size[2]/2

    # Check for NaNs or Infs
    if np.isnan(testgrid).any() or np.isinf(testgrid).any():
        print("Warning: testgrid contains NaNs or Infs. Replacing them with zeros.")
        testgrid = np.nan_to_num(testgrid, nan=0.0, posinf=0.0, neginf=0.0)

    # Normalize the data
    data_min = testgrid.min()
    data_max = testgrid.max()
    if data_max - data_min == 0:
        print("Warning: testgrid has zero variance.")
        data_normalized = np.zeros_like(testgrid)
    else:
        data_normalized = (testgrid - data_min) / (data_max - data_min)

    # Create coordinate grids
    x_grid = np.linspace(xmin, xmax, npixx)
    y_grid = np.linspace(ymin, ymax, npixy)
    z_grid = np.linspace(zmin, zmax, npixz)
    X, Y, Z = np.meshgrid(x_grid, y_grid, z_grid, indexing='ij')

    # Flatten the arrays
    X_flat = X.flatten()
    Y_flat = Y.flatten()
    Z_flat = Z.flatten()
    data_flat = data_normalized.flatten()

    # Configure the renderer
    pio.renderers.default = 'browser'  # or 'notebook', 'notebook_connected', etc.

    # Create the figure
    fig = go.Figure(data=go.Volume(
        x=X_flat,
        y=Y_flat,
        z=Z_flat,
        value=data_flat,
        isomin=0.1,
        isomax=1.0,
        opacity=0.1,
        surface_count=20,
        colorscale='Viridis'
    ))

    fig.update_layout(scene=dict(
        xaxis_title='X',
        yaxis_title='Y',
        zaxis_title='Z'
    ))

    # Try to display the figure
    try:
        fig.show()
    except Exception as e:
        print(f"An error occurred: {e}")
        # Save to HTML as a fallback
        fig.write_html('volume_plot.html')
        print("Figure saved as 'volume_plot.html'. Open this file in a web browser to view the plot.")


# def visualize_amr_grid_2d(x_grid, y_grid, particles, z_slice=None, delta_z=None):
#     """
#     Visualize the AMR grid in the XY plane at a given Z slice.

#     Parameters:
#     - x_grid, y_grid: Arrays of cell boundaries along X and Y axes.
#     - particles: Numpy array of particle data.
#     - z_slice: The Z value at which to take the slice. If None, uses the middle of the Z range.
#     - delta_z: Thickness of the slice. If None, defaults to 1% of the Z range.
#     """
#     if z_slice is None:
#         z_slice = 0.0  # Middle of the Z range

#     if delta_z is None:
#         delta_z = (particles[:, 3].max() - particles[:, 3].min()) * 0.01  # 1% of the Z range

#     # Extract particles near the z_slice
#     mask = np.abs(particles[:, 3] - z_slice) < delta_z
#     particles_slice = particles[mask]

#     fig, ax = plt.subplots(figsize=(10, 10))

#     # Plot vertical grid lines (X)
#     for x in x_grid:
#         ax.axvline(x, color='black', linewidth=0.5)

#     # Plot horizontal grid lines (Y)
#     for y in y_grid:
#         ax.axhline(y, color='black', linewidth=0.5)

#     # Plot particles in the slice
#     ax.scatter(particles_slice[:, 1], particles_slice[:, 2], s=1, color='red', alpha=0.5)

#     ax.set_xlim(x_grid[0], x_grid[-1])
#     ax.set_ylim(y_grid[0], y_grid[-1])
#     ax.set_xlabel('X')
#     ax.set_ylabel('Y')
#     ax.set_title(f'AMR Grid Visualization at Z = {z_slice}')
#     plt.show()


# def visualize_amr_grid_3d(x_grid, y_grid, z_grid, max_cells=1000):
#     """
#     Visualize the AMR grid in 3D using Plotly.

#     Parameters:
#     - x_grid, y_grid, z_grid: Arrays of cell boundaries along X, Y, Z axes.
#     - max_cells: Maximum number of cells to plot to avoid performance issues.
#     """
#     # Create lists to store the cubes
#     cubes = []

#     # Initialize counters
#     cell_count = 0

#     # Iterate over cells
#     for i in range(len(x_grid) - 1):
#         for j in range(len(y_grid) - 1):
#             for k in range(len(z_grid) - 1):
#                 # Check if we've reached the maximum number of cells to plot
#                 if cell_count >= max_cells:
#                     break

#                 # Cell corners
#                 x0, x1 = x_grid[i], x_grid[i + 1]
#                 y0, y1 = y_grid[j], y_grid[j + 1]
#                 z0, z1 = z_grid[k], z_grid[k + 1]
                
#                 # Configure the renderer
#                 pio.renderers.default = 'browser'  # or 'notebook', 'notebook_connected', etc.

#                 # Create the cube
#                 cube = go.Mesh3d(
#                     x=[x0, x1, x1, x0, x0, x1, x1, x0],
#                     y=[y0, y0, y1, y1, y0, y0, y1, y1],
#                     z=[z0, z0, z0, z0, z1, z1, z1, z1],
#                     i=[0, 0, 0, 4, 4, 5, 7, 6, 6, 5, 2, 1],
#                     j=[1, 2, 3, 5, 7, 6, 2, 5, 4, 6, 3, 0],
#                     k=[2, 3, 1, 6, 5, 7, 5, 4, 7, 7, 0, 4],
#                     opacity=0.1,
#                     color='blue',
#                     flatshading=True,
#                     showscale=False
#                 )
#                 cubes.append(cube)
#                 cell_count += 1

#     # Create the figure
#     fig = go.Figure(data=cubes)
#     fig.update_layout(scene=dict(
#         xaxis=dict(title='X'),
#         yaxis=dict(title='Y'),
#         zaxis=dict(title='Z'),
#         aspectmode='data'
#     ))
#     fig.show()

def visualize_amr_grid_2d(x_grid, y_grid, z_grid, particles, cell_levels, z_slice=0.0, delta_z=None):
    """
    Visualize the AMR grid in the XY plane at a given Z slice, showing refinement levels.
    
    Parameters:
    - x_grid, y_grid, z_grid: 1D arrays defining the grid boundaries in each direction.
    - particles: Array of particle data with coordinates in columns.
    - cell_levels: Dictionary with (i, j, k) as keys and refinement levels as values.
    - z_slice: The z-coordinate of the slice to visualize.
    - delta_z: Thickness range around z_slice for including particles.
    """
    # Set delta_z if not provided
    if delta_z is None:
        delta_z = (particles[:, 3].max() - particles[:, 3].min()) * 0.01  # 1% of the Z range

    # Extract particles near the z_slice
    mask = np.abs(particles[:, 3] - z_slice) < delta_z
    particles_slice = particles[mask]

    # Determine the closest k index to z_slice
    k_slice = np.abs(z_grid - z_slice).argmin()
    if z_grid[k_slice] != z_slice:
        print(f"Closest z-grid slice to z={z_slice} is z={z_grid[k_slice]} at index {k_slice}")

    # Define the number of cells in x and y dimensions
    nx = len(x_grid) - 1
    ny = len(y_grid) - 1
    refinement_array = np.zeros((ny, nx))  # Note: Rows = y, Columns = x

    # Populate the refinement_array based on cell_levels
    for (i, j, k), level in cell_levels.items():
        if k == k_slice and 0 <= i < nx and 0 <= j < ny:
            refinement_array[j, i] = level  # Note: refinement_array[j, i] aligns with imshow

    # Plot the AMR grid
    plt.figure(figsize=(10, 8))
    plt.imshow(refinement_array, origin='lower', extent=[x_grid[0], x_grid[-1], y_grid[0], y_grid[-1]],
               cmap='viridis', interpolation='none')
    plt.colorbar(label='Refinement Level')
    plt.scatter(particles_slice[:, 1], particles_slice[:, 2], s=10, color='red', alpha=1.0)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(f'AMR Grid Refinement Levels at Z = {z_slice}')
    plt.show()
    

def compare_mass_and_volume(testgrid_full, particles, box_size):
    """
    Compares the total mass and volume from particles and grid data.

    Parameters:
    - testgrid_full: 4D ndarray, grid data where the first dimension represents mass,
      and the eighth dimension represents density.
    - box_size: list or array-like, the size of the simulation box in each dimension [x, y, z].
    - particles: 2D ndarray, particle data where the first column is mass and the eighth column is density.

    The function calculates:
    - Total mass from the grid and particles.
    - Total volume from the grid (both total box volume and calculated from mass/density).
    - Total volume from the particles calculated from mass/density.
    - Differences and ratios between the mass and volume calculations.
    """
    # Total mass from the grid
    totmassofgrid = np.sum(testgrid_full[:, :, :, 0])

    # Create a mask for valid density values in the grid
    valid_mask = (testgrid_full[:, :, :, 7] > 0) & np.isfinite(testgrid_full[:, :, :, 7])

    # Calculate the total volume from the grid using valid entries
    totvolofgrid2 = np.sum(testgrid_full[:, :, :, 0][valid_mask] / testgrid_full[:, :, :, 7][valid_mask])

    # Total volume of the grid (box volume)
    totvolofgrid = box_size[0] * box_size[1] * box_size[2]

    # Total mass from particles
    total_particle_mass = np.sum(particles[:, 0])

    # Total volume from particles calculated from mass/density
    total_particle_vol = np.sum(particles[:, 0] / particles[:, 7])

    # Printing the results
    print(f"Total mass from particles: {total_particle_mass}")
    print(f"Total mass from grid: {totmassofgrid}")
    print(f"Mass Difference (particles - grid): {total_particle_mass - totmassofgrid}")
    print(f"Mass Ratio (particles/grid): {total_particle_mass / totmassofgrid}\n")

    print(f"Total volume from particles: {total_particle_vol}")
    print(f"Total volume from grid (box volume): {totvolofgrid}")
    print(f"Total volume from grid (mass/density): {totvolofgrid2}")
    print(f"Volume Difference (particles - grid box volume): {total_particle_vol - totvolofgrid}")
    print(f"Volume Ratio (particles/grid box volume): {total_particle_vol / totvolofgrid}\n")


