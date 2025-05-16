from pprint import pprint

from process_bigraph import Process, Step, Composite, ProcessTypes
from process_bigraph.emitter import gather_emitter_results
from spatial_transport.utils import get_regular_edges, generate_compartments, generate_shared_environments

class SimpleDiffusion(Process):
    """Simple diffusion between compartments"""
    config_schema = {
        "spacing": "float",
        "substrates": "map[float]",
    }

    def __init__(self, config, core):
        super().__init__(config, core)

        self.substrates = config['substrates']
        self.spacing = config['spacing']

    def inputs(self):
        return {
            "compartments": "compartments",
            "edges": "map[edge_type]",
        }

    def outputs(self):
        return {
            "compartments": "any"
        }

    def update(self, inputs, interval):
        edges = inputs['edges']
        compartments = inputs['compartments']

        update = {
            compartment_id: {
                'counts': {
                    substrate: 0 for substrate in self.substrates.keys()
                },
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
                d_conc = -diffusivity * ((concentration2 - concentration1)/self.spacing) * self.spacing **2 * interval
                update[edge["neighbors"][0]]["counts"][substrate] += -d_conc
                update[edge["neighbors"][1]]["counts"][substrate] += d_conc
        return update

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
    comps = generate_compartments(dims=[3, 3, 0], spacing=1)
    comps = generate_shared_environments(comps, spacing=1, substrates=substrates)
    spec["Compartments"] = comps
    edges = get_regular_edges(comps, spacing=1)
    spec["Edges"] = edges
    # set emitter specs
    spec["emitter"] = {
        "_type": "step",
        "address": "local:ram-emitter",
        "config": {
            "emit": {
                "compartments": "any",
                "global_time": "any",
            }
        },
        "inputs": {
            "compartments": ["Compartments"],
            "global_time": ["global_time"]
        }
    }
    print("Show Specs")
    pprint(spec)
    sim = Composite(
        {
            "state": spec,
        },
        core=core
    )
    sim.run(10)
    results = gather_emitter_results(sim)[("emitter",)]
    print("Show Results")
    pprint(results)

if __name__ == "__main__":
    from spatial_transport import register_types
    # create the core object
    core = ProcessTypes()
    # register data types
    core = register_types(core)
    core.register_process("SimpleDiffusion", SimpleDiffusion)

    run_simple_diffusion(core)