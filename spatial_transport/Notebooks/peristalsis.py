from pprint import pprint
import io
import imageio.v2 as imageio
import matplotlib.pyplot as plt
import matplotlib as mpl
import time

from process_bigraph import Composite, ProcessTypes
from process_bigraph.emitter import emitter_from_wires, gather_emitter_results

from cdFBA import register_types
from cdFBA.utils import get_substrates, make_cdfba_composite, get_injector_spec


from spatial_transport.processes.diffusion import get_simple_diffusion_spec
from spatial_transport.processes.advection import get_simple_advection_spec
from spatial_transport import register_types as register_types2
from spatial_transport.utils import generate_voxels, generate_shared_environment, get_regular_edges, generate_simple_cdfba_composite, plot_concentrations_2d
from spatial_transport.plot import plot_linear_conc


if __name__ == "__main__":
    exchanges = ['EX_lcts_e', 'EX_gal_e']
    all_species = ['EX_lcts_e', 'EX_gal_e', 'deltaGal', 'deltaLac']
    substrates = get_substrates(model_file="iAF1260", exchanges=exchanges)

    model_dict = {
        'deltaGal':'E_no_galE.xml',
        'deltaLac':'E_no_LCTStex.xml'
    }

    voxels = generate_voxels([10,1,0], spacing=1)
    edges = get_regular_edges(voxels, periodic=False, spacing=1)

    compartments = generate_simple_cdfba_composite(voxels, model_dict, exchanges, 1, sub_range=(5,5), bio_range=(0.05,0.05))
    for id, compartment in compartments.items():
        compartment["Shared Environment"]["counts"]["D-Galactose"] = 0
        compartment["Shared Environment"]["concentrations"]["D-Galactose"] = 0

    spec = {"Compartments": compartments, "Edges": edges}

    spec["Simple Advection"] = get_simple_advection_spec(
        spacing=1,
        substrates=substrates,
        advection=[0.25, 0, 0],
        boundary="default",
        interval=0.1
    )

    for species in model_dict.keys():
        substrates.append(species)

    substrates_dict = {
        "Lactose C12H22O11": 0.03,
        "D-Galactose": 0.06,
        "deltaGal": 0.001,
        "deltaLac": 0.001
    }

    spec["Simple Diffusion"] = get_simple_diffusion_spec(
        substrates=substrates_dict,
        interval=0.1
    )

    injector_config = {
        "injection_params": {
            'Lactose C12H22O11': {
                "amount": 5,
                "interval": 5,
            }
        }
    }

    spec["Compartments"]["0"]["Injector"] = get_injector_spec(injector_config, interval=0.1)
    spec["Compartments"]["0"]["Injector"]["inputs"]["global_time"] = ["..", "..", "global_time"]

    spec["emitter"] = emitter_from_wires({
        "global_time": ["global_time"],
        "compartments": [
            "Compartments",
            "*",
            "Shared Environment",
        ],
        "position": [
            "Compartments",
            "*",
            "position",
        ],
    })

    pprint(spec)

    core = ProcessTypes()
    core = register_types(core)
    core = register_types2(core)

    sim = Composite(
        {
            "state": spec,
        },
        core=core
    )
    start_time = time.time()
    sim.run(20)
    print(f"{time.time() - start_time} seconds")
    results = gather_emitter_results(sim)[("emitter",)]
    # pprint(results[0])
    print(len(results))

    time_range = [0, 200]
    interval = 2
    fig, ax = plot_linear_conc(results, substrate="Lactose C12H22O11", time_range=time_range, interval=interval)
    plt.show()