from pprint import pprint
import numpy as np
import random
import matplotlib.pyplot as plt

COMPARTMENTS = "Compartments"

def get_regular_edges(voxels, periodic=False, spacing=1.0):
    """
    Generates list of edge dictionaries for neighbor relationships between regular cubic voxels
    """
    position_to_key = {tuple(v['position']): k for k, v in voxels.items()}
    positions = np.array([v['position'] for v in voxels.values()])

    # Determine bounds in each dimension
    x_vals, y_vals, z_vals = positions[:, 0], positions[:, 1], positions[:, 2]
    x_min, x_max = x_vals.min(), x_vals.max()
    y_min, y_max = y_vals.min(), y_vals.max()
    z_min, z_max = z_vals.min(), z_vals.max()

    # 6-connected neighbor offsets in 3D
    neighbor_offsets = [
        (spacing, 0, 0), (-spacing, 0, 0),
        (0, spacing, 0), (0, -spacing, 0),
        (0, 0, spacing), (0, 0, -spacing)
    ]

    edges = {}
    edge_id = 1
    seen_pairs = set()

    for key, voxel in voxels.items():
        x, y, z = voxel['position']
        for dx, dy, dz in neighbor_offsets:
            nx, ny, nz = x + dx, y + dy, z + dz
            wrapped = False
            if periodic:
                # Apply periodic wrapping
                if nx > x_max:
                    nx = x_min
                    wrapped = True
                if nx < x_min:
                    nx = x_max
                    wrapped = True
                if ny > y_max:
                    ny = y_min
                    wrapped = True
                if ny < y_min:
                    ny = y_max
                    wrapped = True
                if nz > z_max:
                    nz = z_min
                    wrapped = True
                if nz < z_min:
                    nz = z_max
                    wrapped = True

            neighbor_pos = (nx, ny, nz)
            neighbor_key = position_to_key.get(neighbor_pos)
            if neighbor_key and neighbor_key != key:
                edge_key = tuple(sorted([key, neighbor_key]))
                if edge_key not in seen_pairs:
                    edge_label = f"{edge_id}"
                    edges[edge_label] = {}
                    edges[edge_label]['neighbors'] = [f"{comp}" for comp in edge_key]
                    edges[edge_label]['surface_area'] = spacing ** 2
                    edges[edge_label]['periodic'] = wrapped
                    seen_pairs.add(edge_key)
                    edge_id += 1
    return edges

def generate_voxels(dims, spacing):
    """Creates a spec for shared environments in Euclidean Space

    Parameters:
        dims: list of int, number of voxels in each spatial dimension [x, y, z]
        spacing: float, spacing between neighboring voxels
    """

    voxels = {}
    voxel = 0
    for i in np.arange(spacing/2, spacing*dims[0], spacing):
        for j in np.arange(spacing/2, spacing*dims[1], spacing):
            if dims[2] != 0:
                for k in np.arange(spacing/2, spacing*dims[2], spacing):
                    voxels[f"{voxel}"] = {}
                    voxels[f"{voxel}"]["position"] = [float(i), float(j), float(k)]
                    voxel += 1
            else:
                voxels[f"{voxel}"] = {}
                k=0
                voxels[f"{voxel}"]["position"] = [float(i), float(j), float(k)]
                voxel += 1

    return voxels

def generate_shared_environments(voxels, spacing, substrates):
    """Generate random substrate concentrations"""
    volume = spacing ** 3
    for key in voxels.keys():
        voxels[key]['Shared Environment'] = {}
        voxels[key]['Shared Environment']['volume'] = volume
        voxels[key]['Shared Environment']['counts'] = {}
        voxels[key]['Shared Environment']['concentrations'] = {}
        for substrate in substrates:
            count = random.uniform(0, 10)
            voxels[key]['Shared Environment']['counts'][substrate] = count
            voxels[key]['Shared Environment']['concentrations'][substrate] = count/volume
    compartments = voxels
    return compartments

def detect_boundary_positions(compartments, num_dims = 3, spacing=1.0):
    """
    Determines which compartments lie on the boundaries of the 3D domain
    and which specific boundaries (x_min, x_max, y_min, etc.) they touch.

    compartments: dict of form {key: {'position': (x, y, z)}}
    num_dims: int, number of dimensions
    spacing: grid spacing (used for floating point tolerance)

    Returns:
        boundary_info: dict {key: [list of boundary labels]}
                       boundary labels are from {'x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max'}
    """
    import numpy as np

    # Get all positions
    positions = np.array([v['position'] for v in compartments.values()])
    x_vals, y_vals, z_vals = positions[:, 0], positions[:, 1], positions[:, 2]

    x_min, x_max = x_vals.min(), x_vals.max()
    y_min, y_max = y_vals.min(), y_vals.max()
    z_min, z_max = z_vals.min(), z_vals.max()

    tolerance = spacing / 10  # To handle floating-point rounding

    boundary_info = {}

    for key, comp in compartments.items():
        x, y, z = comp['position']
        boundaries = []
        if np.isclose(x, x_min, atol=tolerance): boundaries.append('x_min')
        if np.isclose(x, x_max, atol=tolerance): boundaries.append('x_max')
        if np.isclose(y, y_min, atol=tolerance): boundaries.append('y_min')
        if np.isclose(y, y_max, atol=tolerance): boundaries.append('y_max')
        if num_dims == 3:
            if np.isclose(z, z_min, atol=tolerance): boundaries.append('z_min')
            if np.isclose(z, z_max, atol=tolerance): boundaries.append('z_max')
        compartments[key]["boundaries"] = boundaries

    return compartments

def plot_concentrations_2d(compartments, molecule='glucose', timepoint=None, **kwargs):
    """
    Plots a heatmap of the specified molecule's concentration for each compartment.

    Parameters:
        compartments : dict, A dictionary of compartment data.
        molecule : str, The molecule whose concentration to plot (default: 'glucose').
        timepoints : float, timepoint for the plot
        **kwargs: Additional keyword arguments passed to plt.imshow().

    Returns:
        fig, ax: Matplotlib figure and axis objects.
    """
    # Extract positions and concentrations
    positions = []
    concentrations = []

    for comp in compartments.values():
        loc = comp['position']
        conc = comp['Shared Environment']['concentrations'].get(molecule, np.nan)
        positions.append((loc[0], loc[1]))
        concentrations.append(conc)

    # Build position maps
    positions = np.array(positions)
    concentrations = np.array(concentrations)
    x_coords = sorted(set(loc[0] for loc in positions))
    y_coords = sorted(set(loc[1] for loc in positions))
    x_map = {x: i for i, x in enumerate(x_coords)}
    y_map = {y: i for i, y in enumerate(y_coords)}

    # Create the grid
    grid = np.full((len(y_coords), len(x_coords)), np.nan)
    for (x, y), value in zip(positions, concentrations):
        xi = x_map[x]
        yi = y_map[y]
        grid[yi, xi] = value

    # Plot
    fig, ax = plt.subplots(figsize=(6, 5))
    extent = [min(x_coords) - 0.5, max(x_coords) + 0.5, min(y_coords) - 0.5, max(y_coords) + 0.5]
    im = ax.imshow(grid, origin='lower', extent=extent, **kwargs)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(f'{molecule.capitalize()} Concentration')
    if timepoint is not None:
        ax.set_title(f"t = {timepoint:.2f}")
    else:
        ax.set_title(f'{molecule.capitalize()} Concentration Heatmap')
    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    ax.set_xticks(x_coords)
    ax.set_yticks(y_coords)
    # ax.grid(True, linestyle='--', alpha=0.4)
    fig.tight_layout()
    return fig, ax

if __name__ == "__main__":
    compartments = generate_voxels(dims=[2, 2, 2], spacing=1)
    pprint(compartments)
    compartments2 = generate_voxels(dims=[3, 3, 1], spacing=1)
    pprint(compartments2)
    edges = get_regular_edges(compartments2, spacing=1)
    pprint(edges)
    edges = get_regular_edges(compartments2, periodic=True, spacing=1)
    print("Periodic Boundary Edges")
    pprint(edges)
    substrates = ["glucose", "acetate", "biomass"]
    compartments3 = generate_shared_environments(compartments2, spacing=1, substrates=substrates)
    pprint({"Compartments": compartments3})
    kwargs = []
    fig, ax = plot_concentrations_2d(compartments3, molecule='glucose', cmap='plasma', vmin=0, vmax=10)
    plt.show()

    pprint(detect_boundary_positions(compartments2, num_dims=2, spacing=1))