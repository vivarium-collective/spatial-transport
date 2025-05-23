from pprint import pprint
import io
import imageio.v2 as imageio
import matplotlib.pyplot as plt
import numpy as np

from process_bigraph import Process, Composite, ProcessTypes
from process_bigraph.emitter import emitter_from_wires, gather_emitter_results
from spatial_transport.utils import get_regular_edges, generate_compartments, generate_shared_environments, plot_concentrations_2d

class SimpleAdvection(Process):

    config_schema = {
        "spacing": "float",
        "substrates": "list[string]",
        "advection": "list[float]", #advection velocity vector
    }

    def __init__(self, config, core):

        self.substrates = config['substrates']
        self.spacing = config['spacing']
        self.advection = np.array(config['advection'])

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
            conc1 = compartments[edge["neighbors"][0]]['Shared Environment']['concentrations']
            conc2 = compartments[edge["neighbors"][1]]['Shared Environment']['concentrations']
            pos1 = np.array(compartments[compartment1]["position"])
            pos2 = np.array(compartments[compartment2]["position"])
            normal1 = (pos2 - pos1)/np.linalg.norm(pos2 - pos1)
            # normal2 = np.linalg.norm(pos1 - pos2)
            advect = self.advection
            vn = np.dot(normal1, advect)
            for substrate in self.substrates:
                concentration1 = conc1[substrate]
                concentration2 = conc2[substrate]
                if vn > 0:
                    delta1 = -vn * concentration1 * self.spacing ** 2 * interval
                else:
                    delta1 = -vn * concentration2 * self.spacing ** 2 * interval
                update[edge["neighbors"][0]]["Shared Environment"]["counts"][substrate] += delta1
                update[edge["neighbors"][1]]["Shared Environment"]["counts"][substrate] += -delta1

        return {"compartments": update}

def get_simple_advection_spec(spacing, substrates, advection, interval):
    return {
        "_type": "process",
        "address": "local:SimpleAdvection",
        "config": {
            "spacing": spacing,
            "substrates": substrates,
            "advection": advection,
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

def run_simple_advection(core):
    spec = {}
    substrates = {
        "glucose": 0.06,
        "acetate": 0.12,
    }
    substrate_list = list(substrates.keys())
    advection = [0.1,0.1,0]
    spec["Simple Advection"] = get_simple_advection_spec(spacing=1, substrates=substrate_list, advection=advection, interval = 0.1)
    comps = generate_compartments(dims=[5, 5, 0], spacing=1)
    comps = generate_shared_environments(comps, spacing=1, substrates=substrates)
    spec["Compartments"] = comps
    edges = get_regular_edges(comps, spacing=1)
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
        fig, ax = plot_concentrations_2d(result["compartments"], molecule='glucose', cmap='plasma', vmin=0, vmax=10)

        # Save fig to buffer
        buf = io.BytesIO()
        fig.savefig(buf, format='png')
        buf.seek(0)
        frames.append(imageio.imread(buf))
        plt.close(fig)
    imageio.mimsave('animated_plot.gif', frames, duration=1/60)

if __name__ == "__main__":
    from spatial_transport import register_types
    # create the core object
    core = ProcessTypes()
    # register data types
    core = register_types(core)
    core.register_process("SimpleAdvection", SimpleAdvection)

    run_simple_advection(core)