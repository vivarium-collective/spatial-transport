from pprint import pprint

from process_bigraph import Process, Composite, ProcessTypes
from process_bigraph.emitter import emitter_from_wires, gather_emitter_results
from spatial_transport.utils import get_regular_edges, generate_voxels, generate_shared_environments, plot_concentrations_2d
import io
import imageio.v2 as imageio
import matplotlib.pyplot as plt

#Diffusion Processes

class SimpleDiffusion(Process):
    """Simple diffusion between compartments"""
    config_schema = {
        "spacing": "float",
        "substrates": "map[float]",
        "dimension": "float", # "2D or 3D"
    }

    def __init__(self, config, core):
        super().__init__(config, core)

        self.substrates = config['substrates']
        self.spacing = config['spacing']

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
                        substrate: 0 for substrate in self.substrates.keys()
                    },
                }
            }
            for compartment_id in compartments.keys()}

        for edge_id, edge in edges.items():
            compartment1 = compartments[edge["neighbors"][0]]
            compartment2 = compartments[edge["neighbors"][1]]
            conc1 = compartments[edge["neighbors"][0]]['Shared Environment']['concentrations']
            conc2 = compartments[edge["neighbors"][1]]['Shared Environment']['concentrations']
            for substrate in self.substrates.keys():
                concentration1 = conc1[substrate]
                concentration2 = conc2[substrate]
                diffusivity = self.substrates[substrate]
                d_conc = -diffusivity * (concentration2 - concentration1) * self.spacing ** 2 * interval
                update[edge["neighbors"][0]]["Shared Environment"]["counts"][substrate] += -d_conc
                update[edge["neighbors"][1]]["Shared Environment"]["counts"][substrate] += d_conc
        return {"compartments": update}

def get_simple_diffusion_spec(spacing, substrates, interval):
    return {
        "_type": "process",
        "address": "local:SimpleDiffusion",
        "config": {
            "spacing": spacing,
            "substrates": substrates,
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

def run_simple_diffusion(core):
    spec = {}
    substrates = {
        "glucose": 0.06,
        "acetate": 0.12,
    }
    spec["Simple Diffusion"] = get_simple_diffusion_spec(spacing=1, substrates=substrates, interval = 0.1)
    comps = generate_voxels(dims=[5, 10, 0], spacing=1)
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
    imageio.mimsave('diffusion_plot.gif', frames, duration=1/60)
    counts = []
    for result in results:
        glucose = 0
        for id, comp in result['compartments'].items():
            glucose += comp["Shared Environment"]["counts"]["glucose"]
        counts.append(glucose)
    print(counts)

if __name__ == "__main__":
    from spatial_transport import register_types
    # create the core object
    core = ProcessTypes()
    # register data types
    core = register_types(core)
    core.register_process("SimpleDiffusion", SimpleDiffusion)
    run_simple_diffusion(core)