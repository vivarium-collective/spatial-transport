from pprint import pprint
import numpy as np
import pandas as pd
import seaborn as sns
import math
import random
import matplotlib.pyplot as plt
from cdFBA.utils import get_substrates, make_cdfba_composite

COMPARTMENTS = "Compartments"

def get_regular_edges(voxels, periodic=False, spacing=[1,1,1]):
    """
    Generates list of edge dictionaries for neighbor relationships between regular cubic voxels

    Parameters:
        voxels: dict, output from generate_voxels function
        periodic: bool, True if periodic boundary conditions are desired
        spacing: list, spatial spacing in all three spatial directions [x,y,z],
            must be same as used to generate voxels
    """
    position_to_key = {tuple(v['position']): k for k, v in voxels.items()}
    positions = np.array([v['position'] for v in voxels.values()])

    # Determine bounds in each dimension
    x_vals, y_vals, z_vals = positions[:, 0], positions[:, 1], positions[:, 2]
    x_min, x_max = x_vals.min(), x_vals.max()
    y_min, y_max = y_vals.min(), y_vals.max()
    z_min, z_max = z_vals.min(), z_vals.max()

    dx, dy, dz = spacing

    # 6-connected neighbor offsets in 3D (accounting for anisotropic spacing)
    neighbor_offsets = [
        (dx, 0, 0), (-dx, 0, 0),
        (0, dy, 0), (0, -dy, 0),
        (0, 0, dz), (0, 0, -dz)
    ]

    edges = {}
    edge_id = 1
    seen_pairs = set()

    for key, voxel in voxels.items():
        x, y, z = voxel['position']
        for offset in neighbor_offsets:
            dx, dy, dz = offset
            nx, ny, nz = x + dx, y + dy, z + dz
            wrapped = False

            if periodic:
                # Apply periodic wrapping for each dimension independently
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
                    # Surface area now reflects the appropriate orthogonal face area
                    if offset[0] != 0:
                        surface_area = spacing[1] * spacing[2]
                    elif offset[1] != 0:
                        surface_area = spacing[0] * spacing[2]
                    else:
                        surface_area = spacing[0] * spacing[1]
                    edges[edge_label]['surface_area'] = surface_area
                    edges[edge_label]['position'] = [(a + b) / 2 for a, b in zip([x, y, z], [nx, ny, nz])]
                    edges[edge_label]['periodic'] = wrapped
                    seen_pairs.add(edge_key)
                    edge_id += 1

    return edges

def generate_voxels(dims, spacing=[1,1,1]):
    """Creates a spec for shared environments in Euclidean Space

    Parameters:
        dims: list of int, number of voxels in each spatial dimension [x, y, z]
        spacing: list, spacing between neighboring voxels in each spatial dimension [x, y, z]
    """

    voxels = {}
    voxel = 0
    for i in np.arange(spacing[0]/2, spacing[0]*dims[0], spacing[0]):
        for j in np.arange(spacing[1]/2, spacing[1]*dims[1], spacing[1]):
            if dims[2] != 0:
                if spacing[2] != 0:
                    for k in np.arange(spacing[2]/2, spacing[2]*dims[2], spacing[2]):
                        voxels[f"{voxel}"] = {}
                        voxels[f"{voxel}"]["position"] = [round(float(i), 2), round(float(j), 2), round(float(k), 2)]
                        voxel += 1
                else:
                    raise ValueError("z-spacing cannot be 0")
            else:
                voxels[f"{voxel}"] = {}
                k=spacing[2]/2
                voxels[f"{voxel}"]["position"] = [round(float(i), 2), round(float(j), 2), round(float(k), 2)]
                voxel += 1
    return voxels

def generate_shared_environment(volume, substrates, species, sub_range=(0, 10), bio_range=(0, 0.1)):
    shared_environment = {'volume': volume, 'counts': {}, 'concentrations': {}}
    for substrate in substrates:
        count = random.uniform(sub_range[0], sub_range[1])
        shared_environment['counts'][substrate] = count
        shared_environment['concentrations'][substrate] = count/volume
    for species in species:
        biomass = random.uniform(bio_range[0], bio_range[1])
        shared_environment['counts'][species] = biomass
        shared_environment['concentrations'][species] = biomass/volume
    return shared_environment

def add_shared_environments(voxels, substrates, spacing = [1,1,1]):
    """Generate random substrate concentrations
    Parameters:
        voxels: dict, spatial voxels output from the generate_voxels function
        spacing: list, spacing between neighboring voxels in each spatial dimension [x, y, z]
        substrates: list, substrate names output from the generate_shared_environment function
    """
    volume = spacing[0] * spacing[1] * spacing[2]
    for key in voxels.keys():
        voxels[key]['Shared Environment'] = {}
        voxels[key]['Shared Environment']['volume'] = volume
        voxels[key]['Shared Environment']['counts'] = {}
        voxels[key]['Shared Environment']['concentrations'] = {}
        for substrate in substrates:
            count = random.uniform(0, 1)
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

#cdFBA Utility Functions
def generate_simple_cdfba_composite(voxels, model_dict, exchanges, volume, sub_range=(0, 10), bio_range=(0, 0.1)):
    substrates = []
    species_list = [species for species in model_dict.keys()]
    for species, model in model_dict.items():
        substrates += get_substrates(model_file=model, exchanges=exchanges)
    base_spec = make_cdfba_composite(model_dict=model_dict, exchanges=exchanges, volume=volume, interval=0.1)
    for id in voxels:
        spec = base_spec.copy()
        shared_environment = generate_shared_environment(volume=1, substrates=substrates, species=species_list, sub_range=sub_range, bio_range=bio_range)
        spec["Shared Environment"] = shared_environment
        voxels[id].update(spec)
    return voxels

def plot_concentrations_2d(compartments, molecule='glucose', timepoint=None, **kwargs):
    """
    Plots a heatmap of the specified molecule's concentration for each compartment using Seaborn.

    Parameters:
        compartments : dict
            Dictionary of compartment data.
        molecule : str
            Molecule to plot (default: 'glucose').
        timepoint : float
            Timepoint for the title.
        **kwargs : additional keyword arguments passed to sns.heatmap().

    Returns:
        fig, ax : Matplotlib figure and axis objects.
    """
    # Extract positions and concentrations
    positions = []
    concentrations = []

    for comp in compartments.values():
        loc = comp['position']
        conc = comp['Shared Environment']['concentrations'].get(molecule, np.nan)
        positions.append((loc[0], loc[1]))
        concentrations.append(conc)

    positions = np.array(positions)
    concentrations = np.array(concentrations)

    # Create a DataFrame with X, Y, and concentration
    df = pd.DataFrame({
        'x': positions[:, 0],
        'y': positions[:, 1],
        molecule: concentrations
    })

    # Pivot to make a 2D grid
    grid = df.pivot(index='y', columns='x', values=molecule)

    # Sort axes
    grid = grid.sort_index(ascending=True)  # y-axis
    grid = grid[sorted(grid.columns)]       # x-axis

    # Plot with seaborn
    fig, ax = plt.subplots(dpi=300)
    sns.heatmap(grid, ax=ax, cbar_kws={'label': f'{molecule.capitalize()} Concentration', 'orientation':'horizontal'}, **kwargs)

    if timepoint is not None:
        ax.set_title(f"t = {timepoint:.2f}")
    else:
        ax.set_title(f'{molecule.capitalize()} Concentration Heatmap')

    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')

    fig.tight_layout()
    return fig, ax


def plot_linear_kymograph(results, molecule='glucose', timepoints=None, axis='x', aspect=None, half_ticks=False, **kwargs):
    """plots linear kymograph of compartments along a linear axis for a specified molecule"""
    concentration_dict = {}
    index = ["x", "y", "z"].index(axis)
    compartments = results[0]["compartments"]
    positions = [comp["position"][index] for comp in compartments.values()]

    for timepoint in timepoints:
        timepoint = float("{:.2f}".format(timepoint))
        result = [r for r in results if math.isclose(r["global_time"], timepoint, rel_tol=1e-9)][0]
        comps = result["compartments"]
        concentrations = [
            comp["Shared Environment"]["concentrations"].get(molecule, np.nan)
            for comp in comps.values()
        ]
        concentration_dict.update({timepoint: concentrations})

    concentration_df = pd.DataFrame(concentration_dict)
    concentration_df.index = positions

    fig, ax = plt.subplots(dpi=300)
    sns.heatmap(
        concentration_df,
        ax=ax,
        cbar_kws={"label": f"{molecule.capitalize()} Concentration (mM)"},
        **kwargs
    )

    if aspect is not None:
        ax.set_aspect(aspect)
    if half_ticks:
        ax.set_xticks(ax.get_xticks()[::2])
    ax.set_ylabel(f"{axis}-position")
    ax.set_xlabel("Timepoint")
    ax.set_title("Concentration Kymograph")
    return fig, ax

if __name__ == "__main__":
    spacing = [1, 0.1, 1]
    dims = [2, 5, 1]
    # voxels = generate_voxels(dims=dims, spacing=spacing)
    # pprint(voxels)
    voxels2 = generate_voxels(dims=dims, spacing=spacing)
    pprint(voxels2)
    edges = get_regular_edges(voxels2, spacing=spacing)
    pprint(edges)
#     edges = get_regular_edges(voxels2, periodic=True, spacing=1)
#     print("Periodic Boundary Edges")
#     pprint(edges)
    substrates = ["glucose", "acetate", "biomass"]
    compartments = add_shared_environments(voxels2, substrates=substrates, spacing=spacing)
    pprint({"Compartments": compartments})
    fig, ax = plot_concentrations_2d(compartments, molecule='glucose', cmap='plasma', vmin=0, vmax=10)
    plt.show()
#     pprint(detect_boundary_positions(voxels2, num_dims=2, spacing=1))
#     pprint(generate_shared_environment(volume=1, substrates=["glucose", "acetate"], species=["boimass"]))