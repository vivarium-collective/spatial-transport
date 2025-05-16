from pprint import pprint
import numpy as np
import random
import matplotlib.pyplot as plt

COMPARTMENTS = "Compartments"

def get_regular_edges(compartments, spacing=1.0):
    """
    Generates list of edge dictionaries for neighbor relationships between regular cubic compartments
    """
    position_to_key = {tuple(v['position']): k for k, v in compartments.items()}

    # 6-connected neighbor offsets in 3D
    neighbor_offsets = [
        (spacing, 0, 0), (-spacing, 0, 0),
        (0, spacing, 0), (0, -spacing, 0),
        (0, 0, spacing), (0, 0, -spacing)
    ]

    edges = {}
    edge_id = 1
    seen_pairs = set()

    for key, voxel in compartments.items():
        x, y, z = voxel['position']
        for dx, dy, dz in neighbor_offsets:
            neighbor_pos = (x + dx, y + dy, z + dz)
            neighbor_key = position_to_key.get(neighbor_pos)
            if neighbor_key:
                edge_key = tuple(sorted([key, neighbor_key]))
                if edge_key not in seen_pairs:
                    edge_label = f"{edge_id}"
                    edges[edge_label] = {}
                    edges[edge_label]['neighbors'] = [f"{comp}" for comp in edge_key]
                    edges[edge_label]['surface_area'] = spacing ** 2
                    seen_pairs.add(edge_key)
                    edge_id += 1
    return edges

def generate_compartments(dims, spacing):
    """Creates a spec for shared environments in Euclidean Space

    Parameters:
        dims: list of int, number of compartments in each spatial dimension [x, y, z]
        spacing: float, spacing between neighboring compartments
    """

    compartments = {}
    compartment = 0
    for i in np.arange(spacing/2, spacing*dims[0], spacing):
        for j in np.arange(spacing/2, spacing*dims[1], spacing):
            if dims[2] != 0:
                for k in np.arange(spacing/2, spacing*dims[2], spacing):
                    compartments[f"{compartment}"] = {}
                    compartments[f"{compartment}"]["position"] = [float(i), float(j), float(k)]
                    compartment += 1
            else:
                compartments[f"{compartment}"] = {}
                k=0
                compartments[f"{compartment}"]["position"] = [float(i), float(j), float(k)]
                compartment += 1

    return compartments

def generate_shared_environments(compartments, spacing, substrates):
    """Generate random substrate concentrations"""
    volume = spacing ** 3
    for key in compartments.keys():
        compartments[key]['Shared Environment'] = {}
        compartments[key]['Shared Environment']['volume'] = volume
        compartments[key]['Shared Environment']['counts'] = {}
        compartments[key]['Shared Environment']['concentrations'] = {}
        for substrate in substrates:
            count = random.uniform(0, 10)
            compartments[key]['Shared Environment']['counts'][substrate] = count
            compartments[key]['Shared Environment']['concentrations'][substrate] = count/volume
    return compartments

def plot_concentrations_2d(compartments, molecule='glucose', **kwargs):
    """
    Plots a heatmap of the specified molecule's concentration for each compartment.

    Parameters:
        compartments (dict): A dictionary of compartment data.
        molecule (str): The molecule whose concentration to plot (default: 'glucose').
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

    ax.set_title(f'{molecule.capitalize()} Concentration Heatmap')
    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    ax.set_xticks(x_coords)
    ax.set_yticks(y_coords)
    # ax.grid(True, linestyle='--', alpha=0.4)
    fig.tight_layout()

    return fig, ax



if __name__ == "__main__":
    compartments = generate_compartments(dims=[2, 2, 2], spacing=1)
    pprint(compartments)
    compartments2 = generate_compartments(dims=[3, 3, 0], spacing=1)
    pprint(compartments2)
    edges = get_regular_edges(compartments, spacing=1)
    pprint(edges)
    substrates = ["glucose", "acetate", "biomass"]
    compartments3 = generate_shared_environments(compartments2, spacing=1, substrates=substrates)
    pprint({"Compartments": compartments3})
    kwargs = []
    fig, ax = plot_concentrations_2d(compartments3, molecule='glucose', cmap='plasma', vmin=0, vmax=10)
    plt.show()