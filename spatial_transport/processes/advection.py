from pprint import pprint
import io
import imageio.v2 as imageio
import matplotlib.pyplot as plt
import numpy as np
import math

from process_bigraph import Process, Composite, Step, ProcessTypes
from process_bigraph.emitter import emitter_from_wires, gather_emitter_results
from spatial_transport.utils import get_regular_edges, generate_voxels, add_shared_environments, plot_concentrations_2d, detect_boundary_positions

ADVECTION = "Advection"

class SimpleAdvection(Process):
    """
    This process performs a simple advection on a network of volume elements, along a fixed vector direction.
    """

    config_schema = {
        "spacing": "list[float]",
        "substrates": "list[string]",
        "advection": "list[float]", #advection velocity vector
        "boundary": "string", # default or periodic
    }

    def __init__(self, config, core):
        super().__init__(config, core)

        self.substrates = config['substrates']
        self.spacing = config['spacing']
        self.advection = np.array(config['advection'])
        self.boundary = config['boundary']

    def inputs(self):
        return {
            "compartments": "map[compartment]",
            "edges": "map[edge_type]",
        }

    def outputs(self):
        return {
            "compartments": "map[compartment]"
        }

    def update(self, inputs, interval):
        edges = inputs['edges']
        compartments = inputs['compartments']

        # position_to_key = {tuple(v['position']): k for k, v in compartments.items()}
        positions = np.array([v['position'] for v in compartments.values()])

        # Determine bounds in each dimension
        x_vals, y_vals, z_vals = positions[:, 0], positions[:, 1], positions[:, 2]
        x_min, x_max = x_vals.min(), x_vals.max()
        y_min, y_max = y_vals.min(), y_vals.max()
        z_min, z_max = z_vals.min(), z_vals.max()
        max_values = {'x': x_max, 'y': y_max, 'z': z_max}

        update = {
            compartment_id: {
                "Shared Environment": {
                    'counts': {
                        substrate: 0 for substrate in self.substrates
                    },
                }
            }
            for compartment_id in compartments.keys()}

        for edge_id, edge in edges.items():
            compartment1 = edge["neighbors"][0]
            compartment2 = edge["neighbors"][1]
            area = edge["surface_area"]
            conc1 = compartments[compartment1]['Shared Environment']['concentrations']
            conc2 = compartments[compartment2]['Shared Environment']['concentrations']
            pos1 = np.array(compartments[compartment1]["position"])
            pos2 = np.array(compartments[compartment2]["position"])
            if self.boundary == "default":
                normal1 = (pos2 - pos1) / np.linalg.norm(pos2 - pos1)
            if self.boundary == "periodic":
                if not edge["periodic"]:
                    normal1 = (pos2 - pos1)/np.linalg.norm(pos2 - pos1)
                if edge["periodic"]:
                    boundaries1 = compartments[compartment1]['boundaries']
                    boundaries2 = compartments[compartment2]['boundaries']
                    for i, dim in enumerate(['x', 'y', 'z']):
                        dim_max = max_values[dim]
                        if f"{dim}_max" in boundaries2 and f"{dim}_min" in boundaries1:
                            pos1[i] += (dim_max + self.spacing/2)
                        if f"{dim}_max" in boundaries1 and f"{dim}_min" in boundaries2:
                            pos2[i] += (dim_max + self.spacing/2)
                    normal1 = (pos2 - pos1) / np.linalg.norm(pos2 - pos1)
            advect = self.advection
            vn = np.dot(normal1, advect)
            for substrate in self.substrates:
                concentration1 = conc1[substrate]
                concentration2 = conc2[substrate]
                if vn > 0:
                    delta1 = -vn * concentration1 * area * interval
                else:
                    delta1 = -vn * concentration2 * area * interval
                update[edge["neighbors"][0]]["Shared Environment"]["counts"][substrate] += delta1
                update[edge["neighbors"][1]]["Shared Environment"]["counts"][substrate] += -delta1

        return {"compartments": update}

class DynamicAdvection(Process):
    """
    Performs dynamic advection on a network of volume elements. Advection direction is based on
    a vector space that assigns an advection direction at each volume edge. Reads a dictionary mapping
    edge_ids to the corresponding vector as input.
    """

    config_schema = {
        "spacing": "list[float]",
        "substrates": "list[string]",
        "boundary": "string",  # default or periodic
    }

    def __init__(self, config, core):
        super().__init__(config, core)

        self.substrates = config['substrates']
        self.spacing = config['spacing']
        self.boundary = config['boundary']

    def inputs(self):
        return {
            "compartments": "map[compartment]",
            "edges": "map[edge_type]",
            "advection": "map[set_list]",
        }

    def outputs(self):
        return {
            "compartments": "map[compartment]"
        }

    def update(self, inputs, interval):
        edges = inputs['edges']
        compartments = inputs['compartments']

        # position_to_key = {tuple(v['position']): k for k, v in compartments.items()}
        positions = np.array([v['position'] for v in compartments.values()])

        # Determine bounds in each dimension
        x_vals, y_vals, z_vals = positions[:, 0], positions[:, 1], positions[:, 2]
        x_min, x_max = x_vals.min(), x_vals.max()
        y_min, y_max = y_vals.min(), y_vals.max()
        z_min, z_max = z_vals.min(), z_vals.max()
        max_values = {'x': x_max, 'y': y_max, 'z': z_max}

        # Initialize update dictionary
        update = {
            compartment_id: {
                "Shared Environment": {
                    'counts': {
                        substrate: 0 for substrate in self.substrates
                    },
                }
            }
            for compartment_id in compartments.keys()}

        # Iterate through all edges to calculate mass transfer across them
        for edge_id, edge in edges.items():
            area = edge["surface_area"]
            compartment1 = edge["neighbors"][0]
            compartment2 = edge["neighbors"][1]
            conc1 = compartments[compartment1]['Shared Environment']['concentrations']
            conc2 = compartments[compartment2]['Shared Environment']['concentrations']
            pos1 = np.array(compartments[compartment1]["position"])
            pos2 = np.array(compartments[compartment2]["position"])
            # Check for periodic boundaries
            if self.boundary == "default":
                normal1 = (pos2 - pos1) / np.linalg.norm(pos2 - pos1)
            if self.boundary == "periodic":
                if not edge["periodic"]:
                    normal1 = (pos2 - pos1) / np.linalg.norm(pos2 - pos1)
                if edge["periodic"]:
                    boundaries1 = compartments[compartment1]['boundaries']
                    boundaries2 = compartments[compartment2]['boundaries']
                    # Shift boundary compartment position for normal calculation
                    for i, dim in enumerate(['x', 'y', 'z']):
                        dim_max = max_values[dim]
                        if f"{dim}_max" in boundaries2 and f"{dim}_min" in boundaries1:
                            pos1[i] += (dim_max + self.spacing / 2)
                        if f"{dim}_max" in boundaries1 and f"{dim}_min" in boundaries2:
                            pos2[i] += (dim_max + self.spacing / 2)
                    # Calculate edge normal
                    normal1 = (pos2 - pos1) / np.linalg.norm(pos2 - pos1)
            # Retrieve advection vector for edge
            advect = np.array(inputs['advection'][edge_id])
            # Calculate dot product of normal and advection
            vn = np.dot(normal1, advect)
            # Calculate direction and amount of mass transfer
            for substrate in self.substrates:
                concentration1 = conc1[substrate]
                concentration2 = conc2[substrate]
                if vn > 0:
                    delta1 = -vn * concentration1 * area * interval
                else:
                    delta1 = -vn * concentration2 * area * interval
                update[edge["neighbors"][0]]["Shared Environment"]["counts"][substrate] += delta1
                update[edge["neighbors"][1]]["Shared Environment"]["counts"][substrate] += -delta1

        return {"compartments": update}

class Peristalsis(Process):
    """This step calculates a vector field by assigning vectors to each edge in the simulation
    to simulate peristaltic fluid flow along the x-axis. A gaussian distribution is used to generate the
    wave pulses.
    """
    config_schema = {
        "amplitude": "float",
        "velocity": "float",
        "wavelength": "float",
        "period": "float",
        "direction": "list[float]",
    }

    def __init__(self, config, core):
        super().__init__(config, core)

        self.amplitude = config["amplitude"]
        self.velocity = config["velocity"]
        self.wavelength = config["wavelength"]
        self.period = config["period"]
        self.direction = config["direction"]

    def inputs(self):
        return {
            "edges": "map[edge_type]",
            "global_time": "float",
        }

    def outputs(self):
        return {
            "advection": "map[set_list]",
        }

    def update(self, inputs, interval):
        advection = {}
        t = inputs["global_time"]

        # Normalize direction vector once
        dx, dy, dz = self.direction
        norm = math.sqrt(dx * dx + dy * dy + dz * dz)
        dx, dy, dz = dx / norm, dy / norm, dz / norm
        direction = np.array([dx, dy, dz])

        # Precompute Gaussian denominator
        sigma = self.wavelength / 2.355
        denom = 2 * sigma * sigma

        # Precompute wave centers (1D array)
        n_waves = int(t // self.period) + 1
        base_shift = self.velocity * t
        interval_shift = self.velocity * self.period
        wave_centers = base_shift - np.arange(n_waves) * interval_shift

        # Collect edge positions into array (N, 3)
        edge_ids = list(inputs["edges"].keys())
        positions = np.array([inputs["edges"][eid]["position"] for eid in edge_ids])

        # Project positions along direction (N,)
        projections = positions @ direction

        # Compute all distances (N, n_waves) by broadcasting
        distances = projections[:, None] - wave_centers[None, :]

        # Gaussian envelopes (N, n_waves)
        envelopes = np.exp(-(distances ** 2) / denom)

        # Sum contributions over waves (N,)
        contribs = envelopes.sum(axis=1) * self.amplitude

        # Multiply by direction vector → (N, 3)
        velocities = contribs[:, None] * direction

        # Map back to dictionary
        for eid, vel in zip(edge_ids, velocities):
            advection[eid] = tuple(vel)

        return {"advection": advection}

def get_simple_advection_spec(spacing, substrates, advection, boundary, interval):
    return {
        "_type": "process",
        "address": "local:SimpleAdvection",
        "config": {
            "spacing": spacing,
            "substrates": substrates,
            "advection": advection,
            "boundary": boundary,
        },
        "inputs": {
            "compartments": ["Compartments"],
            "edges": ["Edges"]
        },
        "outputs": {
            "compartments": ["Compartments"],
        },
        "interval": interval
    }

def get_dynamic_advection_spec(spacing, substrates, boundary, interval):
    return {
        "_type": "process",
        "address": "local:DynamicAdvection",
        "config": {
            "spacing": spacing,
            "substrates": substrates,
            "boundary": boundary,
            "interval": interval,
        },
        "inputs": {
            "compartments": ["Compartments"],
            "edges": ["Edges"],
            "advection": [ADVECTION]
        },
        "outputs": {
            "compartments": ["Compartments"],
        },
        "interval": interval
    }

def get_peristalsis_spec(amplitude, velocity, wavelength, period, direction, interval):
    return {
        "_type": "process",
        "address": "local:Peristalsis",
        "config": {
            "amplitude": amplitude,
            "velocity": velocity,
            "wavelength": wavelength,
            "period": period,
            "direction": direction,
            "interval": interval,
        },
        "inputs": {
            "edges": ["Edges"],
            "global_time": ["global_time"],
        },
        "outputs": {
            "advection": [ADVECTION],
        }
    }

def run_simple_advection(core):
    spec = {}
    substrates = {
        "glucose": 0.06,
        "acetate": 0.12,
    }
    substrate_list = list(substrates.keys())
    advection = [0.5,0.5,0]
    spec["Simple Advection"] = get_simple_advection_spec(spacing=1, substrates=substrate_list, advection=advection, boundary="default", interval=0.1)
    comps = generate_voxels(dims=[10, 10, 0], spacing=1)
    comps = add_shared_environments(comps, spacing=1, substrates=substrates)
    comps = detect_boundary_positions(comps, num_dims=2, spacing=1)
    spec["Compartments"] = comps
    edges = get_regular_edges(comps, periodic=False, spacing=1)
    spec["Edges"] = edges
    # set emitter specs
    spec["emitter"] = emitter_from_wires({
        "global_time": ["global_time"],
        'compartments': ['Compartments'],
    })
    print("Show Specs")
    pprint(spec)
    sim = Composite(
        {
            "state": spec,
        },
        core=core
    )
    sim.run(20)
    results = gather_emitter_results(sim)[("emitter",)]
    frames = []
    for result in results:
        fig, ax = plot_concentrations_2d(result["compartments"], molecule='glucose', timepoint=result["global_time"], cmap='plasma', vmin=0, vmax=10)

        # Save fig to buffer
        buf = io.BytesIO()
        fig.savefig(buf, format='png')
        buf.seek(0)
        frames.append(imageio.imread(buf))
        plt.close(fig)
    imageio.mimsave('advection_plot.gif', frames, duration=1/60)

def run_dynamic_advection(core):
    spec = {}
    substrates = {
        "glucose": 0.06,
        "acetate": 0.12,
    }
    substrate_list = list(substrates.keys())
    advection = [0.5,0.5,0]
    spec["Dynamic Advection"] = get_dynamic_advection_spec(spacing=1, substrates=substrate_list, boundary="default", interval=0.1)
    comps = generate_voxels(dims=[24, 1, 0], spacing=1)
    comps = add_shared_environments(comps, spacing=1, substrates=substrates)
    comps = detect_boundary_positions(comps, num_dims=2, spacing=1)
    spec["Compartments"] = comps
    edges = get_regular_edges(comps, periodic=False, spacing=1)
    spec["Edges"] = edges
    spec["Peristalsis"] = get_peristalsis_spec(amplitude=6, velocity=3, wavelength=2, period=3, direction=[1, 0, 0], interval=0.1)
    spec[ADVECTION] = {edge_id:[0,0,0] for edge_id in edges}
    # set emitter specs
    spec["emitter"] = emitter_from_wires({
        "global_time": ["global_time"],
        'compartments': ['Compartments'],
    })
    print("Show Specs")
    pprint(spec)
    sim = Composite(
        {
            "state": spec,
        },
        core=core
    )
    sim.run(20)
    results = gather_emitter_results(sim)[("emitter",)]
    frames = []
    for result in results:
        fig, ax = plot_concentrations_2d(result["compartments"], molecule='glucose', timepoint=result["global_time"], cmap='plasma', vmin=0, vmax=10)

        # Save fig to buffer
        buf = io.BytesIO()
        fig.savefig(buf, format='png')
        buf.seek(0)
        frames.append(imageio.imread(buf))
        plt.close(fig)
    imageio.mimsave('advection_plot.gif', frames, duration=1/60)

if __name__ == "__main__":
    from spatial_transport import register_types
    # create the core object
    core = ProcessTypes()
    # register data types
    core = register_types(core)
    core.register_process("SimpleAdvection", SimpleAdvection)

    # run_simple_advection(core)
    run_dynamic_advection(core)