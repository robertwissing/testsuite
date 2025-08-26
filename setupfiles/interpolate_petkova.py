## Script to conserve mass when going from particles to grid, translated to python from Petkova 2018 C script (vertex_integral) function

import math
import random
import pyvoro
from numba import njit
from itertools import product
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.colors import LogNorm  # For log-scale plotting

@njit
def compute_logs(u, tol=1e-12):
    # Handle edge cases with a tolerance
    if u >= 1.0 - tol:
        # As u approaches 1 from below: log(1+u) - log(1-u) → ∞
        return 1000.0
    elif u <= -1.0 + tol:
        # As u approaches -1 from above: log(1+u) - log(1-u) → -∞
        return -1000.0
    else:
        # Use atanh for numerical stability: 2*atanh(u) = log((1+u)/(1-u))
        return 2.0 * math.atanh(u)

@njit
def vertex_integral(phi, r0, R_0, h):
    if r0 == 0.0 or R_0 == 0.0 or phi == 0.0:
        return 0.0

    r03 = r0 ** 3
    r0h = r0 / h
    r0h2 = r0h ** 2
    r0h3 = r0h ** 3
    r0h_2 = (h / r0) ** 2
    r0h_3 = (h / r0) ** 3

    B1, B2, B3 = 0.0, 0.0, 0.0

    if r0 >= 2.0 * h:
        B3 = (h ** 3) / 4.0
    elif r0 > h:
        B3 = (r03 / 4.0) * (-4.0/3.0 + r0h - 0.3*r0h2 + (1.0/30.0)*r0h3 - (1.0/15.0)*r0h_3 + (8.0/5.0)*r0h_2)
        B2 = (r03 / 4.0) * (-4.0/3.0 + r0h - 0.3*r0h2 + (1.0/30.0)*r0h3 - (1.0/15.0)*r0h_3)
    else:
        B3 = (r03 / 4.0) * (-2.0/3.0 + 0.3*r0h2 - 0.1*r0h3 + (7.0/5.0)*r0h_2)
        B2 = (r03 / 4.0) * (-2.0/3.0 + 0.3*r0h2 - 0.1*r0h3 - (1.0/5.0)*r0h_2)
        B1 = (r03 / 4.0) * (-2.0/3.0 + 0.3*r0h2 - 0.1*r0h3)

    a = R_0 / r0
    a2 = a ** 2
    linedist = math.sqrt(r0**2 + R_0**2)
    R = R_0 / math.cos(phi) #if math.cos(phi) != 0 else float('inf')
    r = math.sqrt(r0**2 + R**2) #if R != float('inf') else float('inf')

    full_int = 0.0
    D2 = 0.0
    D3 = 0.0

    if linedist < 1.0 * h:
        try:
            sqrt_val = math.sqrt(h**2 - r0**2)
            phi1 = math.acos(R_0 / sqrt_val) #if sqrt_val > 0 else 0.0
        except:
            phi1 = 0.0

        cosp = math.cos(phi1)
        cosp2 = cosp ** 2
        mu = cosp / a / math.sqrt(1 + cosp2/a2) #if a != 0 else 0.0
        tanp = math.tan(phi1) #if math.cos(phi1) != 0 else float('inf')

        I0 = phi1
        I_2 = phi1 + a2 * tanp
        I_4 = phi1 + 2 * a2 * tanp + (a2**2 / 3.0) * tanp * (2 + 1.0 / cosp2)

        u = math.sin(phi1) * math.sqrt(1 - mu**2) if (1 - mu**2) >= 0 else 0.0
        #logs = math.log(1 + u) - math.log(1 - u) if u not in (1.0, -1.0) else 0.0
        logs=compute_logs(u)
        I1 = math.atan(u / a) #if a != 0 else 0.0

        I_1 = (a / 2.0) * logs + I1
        I_3 = I_1 + a * (1 + a2) / 4.0 * (2 * u / (1 - u**2) + logs) if (1 - u**2) != 0 else I_1 + a * (1 + a2) / 4.0 * logs
        I_5 = I_3 + a * (1 + a2)**2 / 16.0 * ((10*u - 6*u**3)/(1 - u**2)**2 + 3 * logs) if (1 - u**2) != 0 else I_3+a * (1 + a2)**2 / 16.0 * 3*logs

        D2 = (-1.0/6 * I_2 + 0.25 * r0h * I_3 - 0.15 * r0h2 * I_4 + 
              (1.0/30) * r0h3 * I_5 - (1.0/60) * r0h_3 * I1 + 
              (B1 - B2) / r03 * I0)

        try:
            sqrt_val = math.sqrt(4*h**2 - r0**2)
            phi2 = math.acos(R_0 / sqrt_val)  #if sqrt_val > 0 else 0.0
        except:
            phi2 = 0.0

        cosp = math.cos(phi2)
        cosp2 = cosp ** 2
        mu = cosp / a / math.sqrt(1 + cosp2/a2) #if a != 0 else 0.0
        tanp = math.tan(phi2) #if math.cos(phi2) != 0 else float('inf')

        I0 = phi2
        I_2 = phi2 + a2 * tanp
        I_4 = phi2 + 2 * a2 * tanp + (a2**2 / 3.0) * tanp * (2 + 1.0 / cosp2)

        u = math.sin(phi2) * math.sqrt(1 - mu**2) if (1 - mu**2) >= 0 else 0.0
        #logs = math.log(1 + u) - math.log(1 - u) #if u not in (1.0, -1.0) else 0.0
        logs=compute_logs(u)
        I1 = math.atan(u / a) #if a != 0 else 0.0

        I_1 = (a / 2.0) * logs + I1
        I_3 = I_1 + a * (1 + a2) / 4.0 * (2 * u / (1 - u**2) + logs) if (1 - u**2) != 0 else I_1 + a * (1 + a2) / 4.0 * logs
        I_5 = I_3 + a * (1 + a2)**2 / 16.0 * ((10*u - 6*u**3)/(1 - u**2)**2 + 3 * logs) if (1 - u**2) != 0 else I_3+a * (1 + a2)**2 / 16.0 * 3*logs

        D3 = (1.0/3 * I_2 - 0.25 * r0h * I_3 + 0.075 * r0h2 * I_4 - 
              (1.0/120)*r0h3 * I_5 + (4.0/15)*r0h_3 * I1 + 
              (B2 - B3)/r03 * I0 + D2)

    elif linedist < 2.0 * h:
        try:
            sqrt_val = math.sqrt(4*h**2 - r0**2)
            phi2 = math.acos(R_0 / sqrt_val) #if sqrt_val > 0 else 0.0
        except:
            phi2 = 0.0

        cosp = math.cos(phi2)
        cosp2 = cosp ** 2
        mu = cosp / a / math.sqrt(1 + cosp2/a2) #if a != 0 else 0.0
        tanp = math.tan(phi2) #if math.cos(phi2) != 0 else float('inf')

        I0 = phi2
        I_2 = phi2 + a2 * tanp
        I_4 = phi2 + 2 * a2 * tanp + (a2**2 / 3.0) * tanp * (2 + 1.0 / cosp2)

        u = math.sin(phi2) * math.sqrt(1 - mu**2) #if (1 - mu**2) >= 0 else 0.0
        logs=compute_logs(u)
        I1 = math.atan(u / a) #if a != 0 else 0.0

        I_1 = (a / 2.0) * logs + I1
        I_3 = I_1 + a * (1 + a2) / 4.0 * (2 * u / (1 - u**2) + logs) if (1 - u**2) != 0 else I_1 + a * (1 + a2) / 4.0 * logs
        I_5 = I_3 + a * (1 + a2)**2 / 16.0 * ((10*u - 6*u**3)/(1 - u**2)**2 + 3 * logs) if (1 - u**2) != 0 else I_3+a * (1 + a2)**2 / 16.0 * 3*logs

        D3 = (1.0/3 * I_2 - 0.25 * r0h * I_3 + 0.075 * r0h2 * I_4 - 
              (1.0/120)*r0h3 * I_5 + (4.0/15)*r0h_3 * I1 + 
              (B2 - B3)/r03 * I0 + D2)

    cosp = math.cos(phi)
    cosp2 = cosp ** 2
    mu = cosp / a / math.sqrt(1 + cosp2/a2) #if a != 0 else 0.0
    tanp = math.tan(phi) #if math.cos(phi) != 0 else float('inf')

    I0 = phi
    I_2 = phi + a2 * tanp
    I_4 = phi + 2 * a2 * tanp + (a2**2 / 3.0) * tanp * (2 + 1.0 / cosp2)

    u = math.sin(phi) * math.sqrt(1 - mu**2) if (1 - mu**2) >= 0 else 0.0
    logs=compute_logs(u)
    I1 = math.atan(u / a) #if a != 0 else 0.0

    I_1 = (a / 2.0) * logs + I1
    I_3 = I_1 + a * (1 + a2) / 4.0 * (2 * u / (1 - u**2) + logs) if (1 - u**2) != 0 else I_1 + a * (1 + a2) / 4.0 * logs
    I_5 = I_3 + a * (1 + a2)**2 / 16.0 * ((10*u - 6*u**3)/(1 - u**2)**2 + 3 * logs) if (1 - u**2) != 0 else I_3+a * (1 + a2)**2 / 16.0 * 3*logs

  # Calculating the integral expression.

    if r < 1.0 * h:
        term1 = (1.0/6 * I_2 - 0.075 * r0h2 * I_4 + (1.0/40) * r0h3 * I_5 + (B1 / r03) * I0)
        full_int = (r0h3 / math.pi) * term1
    elif r < 2.0 * h:
        term1 = 0.25 * (4.0/3 * I_2 - r0h * I_3 + 0.3 * r0h2 * I_4 - 
                       (1.0/30)*r0h3 * I_5 + (1.0/15)*r0h_3 * I1) + (B2 / r03) * I0 + D2
        full_int = (r0h3 / math.pi) * term1
    else:
        term1 = -0.25 * r0h_3 * I1 + (B3 / r03) * I0 + D3
        full_int = (r0h3 / math.pi) * term1

    return full_int


def generate_grid(method,cell_num,x_min,x_max,y_min,y_max,z_min,z_max):
    
    if method=='grid':
        cell_num=int(cell_num**(1./3.))+1
        # Generate grid coordinates
        spacing=(x_max-x_min)/(cell_num-1)
        x_temp = np.arange(x_min + spacing/2, x_max-spacing/2+spacing*0.0001, spacing)  # Centers of cells
        y_temp = np.arange(y_min + spacing/2, y_max-spacing/2+spacing*0.0001, spacing)
        z_temp = np.arange(z_min + spacing/2, z_max-spacing/2+spacing*0.0001, spacing)
        print(x_temp,y_temp,z_temp)
        
        # Create 3D grid points
        xx, yy, zz = np.meshgrid(x_temp, y_temp, z_temp, indexing='ij')
        points = np.vstack([xx.ravel(), yy.ravel(), zz.ravel()]).T
        
        #points = points+np.random*spacing*0.0001 
        print(len(points))
    else:
        # Generate random points
        random.seed(42)  # For reproducibility
        spacing = 1.0;
        points = [
            [
               random.uniform(x_min, x_max),
                random.uniform(y_min, y_max),
                random.uniform(z_min, z_max)
            ] for _ in range(cell_num)
        ]
    
    return points,spacing

def generate_voronoi(points, bounds):
    return pyvoro.compute_voronoi(
        points,
        bounds,
        dispersion=1000.0,
        periodic=[False]*3,
        radii=[]
    )


def refine_point(original_point, original_spacing):
    """Generate 8 points at cube octant centers"""
    offset = original_spacing / 4
    return np.array([
        original_point + np.array([dx, dy, dz])
        for dx in [-offset, offset]
        for dy in [-offset, offset] 
        for dz in [-offset, offset]
    ])

def refine_points(original_points, selected_mask, original_spacing):
    """Refine points with proper parent tracking"""
    unselected_points = original_points[~selected_mask]
    selected_points = original_points[selected_mask]

    delta = original_spacing / 4
    offsets = np.array(list(product([-delta, delta], repeat=3)))
    
    if len(selected_points) > 0:
        # Refined points keep original 3D format
        refined = selected_points[:, np.newaxis, :] + offsets[np.newaxis, :, :]
        refined_points = refined.reshape(-1, 3)
        
        # Track parents separately
        parent_indices = np.repeat(np.where(selected_mask)[0], 8)
    else:
        refined_points = np.empty((0, 3))
        parent_indices = np.array([], dtype=int)
    
    # Return points (N,3) and parents (N,) separately
    return np.vstack([unselected_points, refined_points]), parent_indices








# Constants matching C++ code
x_min, x_max = -1.0, 1.0
y_min, y_max = -1.0, 1.0
z_min, z_max = -1.0, 1.0
cell_num = 2048
h = 0.5
x0, y0, z0 = 0.0, 0.0, 0.0
#radius = 0.0

points,spacing = generate_grid('grid2',cell_num,x_min,x_max,y_min,y_max,z_min,z_max)

selected_mask = np.zeros(len(points), dtype=bool)  # Initialize all False
#selected_mask = np.all(points == [0,0,0], axis=1)
#selected_mask = points[:, 0] > 0.0
selected_mask[5] = True
#selected_mask = np.all(points == [0,0,0], axis=1)
#refined_points = refine_points(original_points, selected_mask, spacing)
#print("Original points count:", len(original_points))
#print("Selected points count:", sum(selected_mask))
#print("Refined points count:", len(refined_points))

#points=refined_points
bounds = np.array([[x_min,x_max],[y_min,y_max],[z_min,z_max]])
voronoi = generate_voronoi(points,bounds)
#voronoi = hierarchical_voronoi_with_flat_surfaces(points, spacing, selected_mask)








def plot_voronoi_grid(voronoi_cells):
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Assign unique colors to cells using a colormap
    colors = plt.cm.tab20(np.linspace(0, 1, len(voronoi_cells)))
    
    for idx, cell in enumerate(voronoi_cells):
        vertices = np.array(cell['vertices'])
        faces = cell['faces']
        
        # Plot each face with cell-specific color
        for face in faces:
            vert_indices = face['vertices']
            if len(vert_indices) < 3:
                continue  # Skip degenerate faces
            
            face_verts = [vertices[i] for i in vert_indices]
            ax.add_collection3d(Poly3DCollection(
                [face_verts],
                alpha=0.25,  # Reduced transparency for clarity
                linewidths=0.8,
                edgecolor='k',
                facecolor=colors[idx]
            ))
    
    # Plot generators (original points)
    generators = np.array([cell['original'] for cell in voronoi_cells])
    ax.scatter(generators[:,0], generators[:,1], generators[:,2], 
               c='k', s=40, marker='o', depthshade=False, label='Generators')
    
    # Configure axes and view
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_zlim(z_min, z_max)
    ax.set_xlabel('X', fontsize=12)
    ax.set_ylabel('Y', fontsize=12)
    ax.set_zlabel('Z', fontsize=12)
    ax.view_init(elev=25, azim=45)  # Optimal viewing angle
    ax.set_title(f'Voronoi Grid ({len(voronoi_cells)} Cells)', fontsize=14)
    ax.legend()
    
    plt.tight_layout()
    plt.show()

# Generate the plot
plot_voronoi_grid(voronoi)






# Initialize lists to store data
#Mcell_values = []
#densities = []


def exact_interpolate(voronoi):
    Mtot = 0.0
    epsilon = 1e-9 # Use epsilon for floating-point comparisons
    Mface_values = []
    Mcell_values = []
    densities = []
    
    for k, cell in enumerate(voronoi):
        vertices = cell['vertices']
        faces = cell['faces']
        vol = cell['volume']
        x, y, z = cell['original']
        Msumface = 0.0
        # rp_val = x*x+y*y+z*z
        # if rp_val > 2*4*h*h:
        #     densities.append(0.0)
        #     Mcell_values.append(0.0)
        #     continue
    
        # Face/vertex processing
        for face in faces:
            Msum = 0.0
            fv_indices = face['vertices']
            if len(fv_indices) < 3:
                continue  # Skip faces with <3 vertices
    
            # Correct indentation: process valid faces
            fv = [
                [vertices[i][0] , vertices[i][1] , vertices[i][2] ]
                for i in fv_indices
            ]
    
            # Plane equation using first three vertices
            v1, v2, v3 = fv[0], fv[1], fv[2]
            A = (v2[1]-v1[1])*(v3[2]-v1[2]) - (v3[1]-v1[1])*(v2[2]-v1[2])
            B = (v2[2]-v1[2])*(v3[0]-v1[0]) - (v3[2]-v1[2])*(v2[0]-v1[0])
            C = (v2[0]-v1[0])*(v3[1]-v1[1]) - (v3[0]-v1[0])*(v2[1]-v1[1])
            D = - (A * v1[0] + B * v1[1] + C * v1[2])
    
            norm = math.sqrt(A**2 + B**2 + C**2)
            if norm < epsilon:
                    continue
    
            r0_val = (A * x0 + B * y0 + C * z0 + D) / norm
            r0_val = r0_val+epsilon;
            ar0 = abs(r0_val)
    
    
            xp = x0 - r0_val * A / norm
            yp = y0 - r0_val * B / norm
            zp = z0 - r0_val * C / norm
            
            # rp_val = xp*xp+yp*yp+zp*zp;
            # if rp_val > 4*h*h+epsilon:
            #     continue
    
            s2 = v1[0]*(v2[1]*v3[2] - v2[2]*v3[1]) + \
                 v1[1]*(v2[2]*v3[0] - v2[0]*v3[2]) + \
                 v1[2]*(v2[0]*v3[1] - v2[1]*v3[0])
            
    
            num_vertices = len(fv_indices)
            for i in range(num_vertices):
                current_idx = fv_indices[i]
                next_idx = fv_indices[(i + 1) % num_vertices]
                #print(fv_indices)
                #print("idx",current_idx,next_idx)
    
                v2_coords = vertices[current_idx]
                v3_coords = vertices[next_idx]
                x2, y2, z2 = v2_coords[0] , v2_coords[1] , v2_coords[2] 
                x3, y3, z3 = v3_coords[0] , v3_coords[1] , v3_coords[2] 
                line_vec = [x3 - x2, y3 - y2, z3 - z2]
                point_vec = [xp - x2, yp - y2, zp - z2]
    
                cross = [
                    point_vec[1] * line_vec[2] - point_vec[2] * line_vec[1],
                    point_vec[2] * line_vec[0] - point_vec[0] * line_vec[2],
                    point_vec[0] * line_vec[1] - point_vec[1] * line_vec[0]
                ]
                cross_mag = math.sqrt(sum(c**2 for c in cross))
                line_len = math.sqrt(sum(lv**2 for lv in line_vec))
                R_0_val = cross_mag / line_len #if line_len != 0 else 0.0
                
                #print(R_0_val)
    
                r12 = math.sqrt((xp - x2)**2 + (yp - y2)**2 + (zp - z2)**2)
                r13 = math.sqrt((xp - x3)**2 + (yp - y3)**2 + (zp - z3)**2)
                r23 = math.sqrt(line_len**2)
    
                if r12 == 0 or r23 == 0:
                    cosa = 0.0
                else:
                    cosa = ((x3 - x2)*(xp - x2) + (y3 - y2)*(yp - y2) + (z3 - z2)*(zp - z2))
                    cosa /= (r12 * r23)
                    cosa = max(min(cosa, 1.0), -1.0)
                    
                R_0_val = R_0_val + epsilon
                """
                phi1 = math.acos(R_0_val / r12) #if R_0_val < r12 and r12 != 0 else 0.0
                phi2 = math.acos(R_0_val / r13) #if R_0_val < r13 and r13 != 0 else 0.0
    
                s1 = xp*(y2*z3 - z2*y3) + yp*(z2*x3 - x2*z3) + zp*(x2*y3 - y2*x3)
                sign_M = -1 if (s1 * s2 * r0_val) <= 0 else 1
                epsilon=1E-10
                if (r12 * math.sin(phi1) >= r23) or (r13 * math.sin(phi2) >= r23):
                    if phi1 >= phi2:
                        integral = vertex_integral(phi1, ar0, R_0_val, h) - vertex_integral(phi2, ar0, R_0_val, h)
                    else:
                        integral = vertex_integral(phi2, ar0, R_0_val, h) - vertex_integral(phi1, ar0, R_0_val, h)
                else:
                    integral = vertex_integral(phi1, ar0, R_0_val, h) + vertex_integral(phi2, ar0, R_0_val, h)
                """
                # Safely compute phi1 and phi2 to avoid division by zero and invalid acos arguments
                phi1 = math.acos(R_0_val / r12) if (r12 != 0 and R_0_val < r12) else 0.0
                phi2 = math.acos(R_0_val / r13) if (r13 != 0 and R_0_val < r13) else 0.0
                
                s1 = xp * (y2 * z3 - z2 * y3) + yp * (z2 * x3 - x2 * z3) + zp * (x2 * y3 - y2 * x3)
                sign_M = -1 if (s1 * s2 * r0_val) <= 0 else 1
                
    
                # Check conditions with epsilon to handle numerical inaccuracies
                condition1 = (r12 * math.sin(phi1)) + epsilon >= r23
                condition2 = (r13 * math.sin(phi2)) + epsilon >= r23
                
                if condition1 or condition2:
                    # Determine the order of phi1 and phi2 for subtraction
                    if phi1 >= phi2:
                        integral = vertex_integral(phi1, ar0, R_0_val, h) - vertex_integral(phi2, ar0, R_0_val, h)
                    else:
                        integral = vertex_integral(phi2, ar0, R_0_val, h) - vertex_integral(phi1, ar0, R_0_val, h)
                else:
                    # Sum the integrals if neither condition is met
                    integral = vertex_integral(phi1, ar0, R_0_val, h) + vertex_integral(phi2, ar0, R_0_val, h)
                            
                M_contribution = sign_M * integral
                Msum += M_contribution  # Now accumulates mass
                
            #print(phi1,phi2,r12 * math.sin(phi1),r13 * math.sin(phi2),r23,round(ar0),R_0_val,h,M_contribution)
            # Store results
            Mface_values.append(Msum)
            Msumface += Msum;
            #print(i,Msum)
    
            #print("Mface",Msum)
        if(Msumface<0.0):
            Msumface = 0.0;
        if(abs(Msumface) < 1E-9):
            Msumface = 1E-10;
        Mface_values = []
        Mcell_values.append(Msumface)
        Mtot += Msumface
        density = abs(Msumface) / vol / 5.0
        densities.append(density)
            
        #print(f"Current cell mass: {Msumface} \n \n")
    print(f"Total mass: {Mtot} (should be close to 1 for h <= 0.5)")
    return Mcell_values,densities


Mcell_values,densities=exact_interpolate(voronoi)




def plot_voronoi_density(voronoi_cells, densities):
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Normalize density values
    norm = Normalize(vmin=min(densities), vmax=max(densities))
    cmap = plt.cm.viridis
    mappable = ScalarMappable(norm=norm, cmap=cmap)
    
    # Define opacity range (adjust alpha_min/max as needed)
    alpha_min = 0.0  # Minimum opacity (most transparent)
    alpha_max = 0.5  # Maximum opacity (least transparent)
    
    for idx, cell in enumerate(voronoi_cells):
        vertices = np.array(cell['vertices'])
        faces = cell['faces']
        density = densities[idx]
        
        # Get color and compute alpha
        color = cmap(norm(density))  # RGBA color from colormap
        alpha = alpha_min + (alpha_max - alpha_min) * norm(density)  # Scale alpha
        
        # Plot each face
        for face in faces:
            vert_indices = face['vertices']
            if len(vert_indices) < 3:
                continue
            
            face_verts = [vertices[i] for i in vert_indices]
            ax.add_collection3d(Poly3DCollection(
                [face_verts],
                alpha=alpha,  # Use density-dependent alpha
                linewidths=0.5,
                edgecolor='k',
                facecolor=color  # RGBA includes alpha, but we override it
            ))
    
    # Add colorbar and labels
    cbar = plt.colorbar(mappable, ax=ax, shrink=0.6)
    cbar.set_label('Density', fontsize=12)
    
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_zlim(z_min, z_max)
    ax.set_xlabel('X', fontsize=12)
    ax.set_ylabel('Y', fontsize=12)
    ax.set_zlabel('Z', fontsize=12)
    ax.view_init(elev=25, azim=45)
    ax.set_title('Voronoi Cell Density Distribution', fontsize=14)
    
    plt.tight_layout()
    plt.show()
    
# Generate the plot
plot_voronoi_density(voronoi, densities)




def plot_voronoi_xy_projection(voronoi_cells, densities=None):
    fig, ax = plt.subplots(figsize=(12, 8))
    
    if densities is not None:
        # Normalize densities for coloring
        norm = Normalize(vmin=min(densities), vmax=max(densities))
        cmap = plt.cm.viridis
        mappable = ScalarMappable(norm=norm, cmap=cmap)
    
    for idx, cell in enumerate(voronoi_cells):
        vertices = np.array(cell['vertices'])
        faces = cell['faces']
        
        # Project vertices to XY plane (ignore Z)
        xy_vertices = vertices[:, :2]
        
        # Assign color (if densities are provided)
        color = cmap(norm(densities[idx])) if densities is not None else 'skyblue'
        
        # Plot each face
        for face in faces:
            vert_indices = face['vertices']
            if len(vert_indices) < 3:
                continue  # Skip degenerate faces
            
            face_verts = [xy_vertices[i] for i in vert_indices]
            poly = plt.Polygon(face_verts, alpha=0.2, edgecolor='k', facecolor=color)
            ax.add_patch(poly)
    
    # Plot generators (original points)
    generators = np.array([cell['original'] for cell in voronoi_cells])
    ax.scatter(generators[:,0], generators[:,1], c='red', s=30, label='Generators')
    
    # Add colorbar if densities are provided
    if densities is not None:
        cbar = plt.colorbar(mappable, ax=ax, shrink=0.6)
        cbar.set_label('Density', fontsize=12)
    
    # Configure axes
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_xlabel('X', fontsize=12)
    ax.set_ylabel('Y', fontsize=12)
    ax.set_title('Voronoi Grid (XY Projection)', fontsize=14)
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.1)
    
    plt.tight_layout()
    plt.show()

# Example usage:
#plot_voronoi_xy_projection(voronoi, densities)  # With densities




def plot_column_density(voronoi_cells, densities=None, resolution=100, axis='z'):
    """
    Plot column density by projecting Voronoi cells onto a 2D plane.
    
    Parameters:
        voronoi_cells (list): List of Voronoi cell dictionaries.
        densities (list): Optional density values for each cell.
        resolution (int): Resolution of the output grid.
        axis (str): Projection axis ('x', 'y', or 'z').
    """
    # Extract all vertices to determine grid bounds
    all_vertices = np.concatenate([np.array(cell['vertices']) for cell in voronoi_cells])
    x_min, x_max = all_vertices[:, 0].min(), all_vertices[:, 0].max()
    y_min, y_max = all_vertices[:, 1].min(), all_vertices[:, 1].max()
    z_min, z_max = all_vertices[:, 2].min(), all_vertices[:, 2].max()

    # Define projection axes
    if axis == 'x':
        proj_axes = (1, 2)  # Project onto YZ plane
        extent = [y_min, y_max, z_min, z_max]
        xlabel, ylabel = 'Y', 'Z'
    elif axis == 'y':
        proj_axes = (0, 2)  # Project onto XZ plane
        extent = [x_min, x_max, z_min, z_max]
        xlabel, ylabel = 'X', 'Z'
    else:  # Default: project onto XY plane
        proj_axes = (0, 1)
        extent = [x_min, x_max, y_min, y_max]
        xlabel, ylabel = 'X', 'Y'

    # Initialize a 2D grid for column density
    grid = np.zeros((resolution, resolution))

    # Compute column density
    for idx, cell in enumerate(voronoi_cells):
        vertices = np.array(cell['vertices'])
        
        # Get projected vertices (drop the axis we're integrating over)
        proj_vertices = vertices[:, proj_axes]
        
        # Approximate cell area in projection
        from scipy.spatial import ConvexHull
        try:
            hull = ConvexHull(proj_vertices)
            area = hull.volume  # For 2D, volume = area
        except:
            area = 0.0  # Skip degenerate cells
        
        # Assign density (default: 1.0 if no densities provided)
        density = densities[idx] if densities is not None else 1.0
        
        # Rasterize the cell's contribution to the grid
        from matplotlib.path import Path
        cell_path = Path(proj_vertices[hull.vertices])
        x_grid = np.linspace(extent[0], extent[1], resolution)
        y_grid = np.linspace(extent[2], extent[3], resolution)
        xx, yy = np.meshgrid(x_grid, y_grid)
        points = np.vstack([xx.ravel(), yy.ravel()]).T
        mask = cell_path.contains_points(points).reshape(resolution, resolution)
        
        # Add density to grid
        grid[mask] += density * area

    # Plot the column density
    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(
        grid.T,  # Transpose to match axis order
        extent=extent,
        origin='lower',
        cmap='viridis',
        norm=LogNorm()  # Use log scale for better contrast
    )
    
    # Add colorbar and labels
    cbar = plt.colorbar(im, ax=ax, label='Column Density')
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(f'Column Density Projection (Along {axis.upper()}-Axis)', fontsize=14)
    
    plt.tight_layout()
    plt.show()

# Example usage:
plot_column_density(voronoi, densities, axis='z')  # XY projection
plot_column_density(voronoi, densities, axis='x')  # YZ projection
plot_column_density(voronoi, densities, axis='y')  # XZ projection


