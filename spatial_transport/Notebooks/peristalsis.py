from pprint import pprint
import io
import imageio.v2 as imageio
import matplotlib.pyplot as plt
import matplotlib as mpl
import time

from process_bigraph import Composite, allocate_core
from process_bigraph.emitter import emitter_from_wires, gather_emitter_results

from cdFBA.utils import get_substrates, make_cdfba_composite, get_injector_spec

from spatial_transport.data_types import register_types
from spatial_transport.processes.diffusion import get_simple_diffusion_spec
from spatial_transport.processes.advection import (
    get_simple_advection_spec, get_peristalsis_spec, get_dynamic_advection_spec,
    get_advection_store, ADVECTION,
)
from spatial_transport.utils import generate_voxels, generate_shared_environment, get_regular_edges, generate_simple_cdfba_composite, plot_concentrations_2d
from spatial_transport.plot import plot_linear_conc


if __name__ == "__main__":
    exchanges = ['EX_lcts_e', 'EX_gal_e']
    all_species = ['EX_lcts_e', 'EX_gal_e', 'deltaGal', 'deltaLac']
    substrates = get_substrates(model_file="iAF1260", exchanges=exchanges)

    model_dict = {
        'deltaGal': 'E_no_galE.xml',
        'deltaLac': 'E_no_LCTStex.xml'
    }

    spacing = [1, 1, 1]
    voxels = generate_voxels([24, 1, 0], spacing=spacing)
    edges = get_regular_edges(voxels, periodic=False, spacing=spacing)

    compartments = generate_simple_cdfba_composite(voxels, model_dict, exchanges, 1, sub_range=(5, 5),
                                                   bio_range=(0.001, 0.001))
    for id, compartment in compartments.items():
        compartment["Shared Environment"]["counts"]["D-Galactose"] = 0
        compartment["Shared Environment"]["concentrations"]["D-Galactose"] = 0
        if id != "0":
            compartment["Shared Environment"]["counts"]["Lactose C12H22O11"] = 0
            compartment["Shared Environment"]["concentrations"]["Lactose C12H22O11"] = 0

    spec = {"Compartments": compartments, "Edges": edges}

    spec["Dynamic Advection"] = get_dynamic_advection_spec(
        spacing=spacing,
        substrates=substrates,
        boundary="default",
        interval=0.1
    )

    spec["Peristalsis"] = get_peristalsis_spec(
        amplitude=3,
        velocity=2,
        wavelength=2,
        period=5,
        direction=[1, 0, 0],
        interval=0.1
    )

    spec[ADVECTION] = get_advection_store(edges)

    for species in model_dict.keys():
        substrates.append(species)

    substrates_dict = {
        "Lactose C12H22O11": 0.06,
        "D-Galactose": 0.12,
        "deltaGal": 0.01,
        "deltaLac": 0.01
    }

    # spec["Simple Diffusion"] = get_simple_diffusion_spec(
    #     substrates=substrates_dict,
    #     interval=0.1
    # )

    injector_config = {
        "injection_params": {
            'Lactose C12H22O11': {
                "amount": 10,
                "interval": 5,
            }
        }
    }

    spec["Compartments"]["0"]["Injector"] = get_injector_spec(injector_config, interval=0.1)
    spec["Compartments"]["0"]["Injector"]["inputs"]["global_time"] = ["..", "..", "global_time"]

    spec["emitter"] = emitter_from_wires({
        "global_time": ["global_time"],
        'compartments': ['Compartments'],
    })

    pprint(spec)

    # register_types pulls in the cdFBA types and processes as well
    core = register_types(allocate_core())

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

    for substrate in substrates:
        frames = []
        for result in results:
            fig, ax = plot_concentrations_2d(result["compartments"], molecule=substrate, timepoint=result["global_time"],
                                             cmap='plasma', vmin=0, vmax=10)

            # Save fig to buffer
            buf = io.BytesIO()
            fig.savefig(buf, format='png')
            buf.seek(0)
            frames.append(imageio.imread(buf))
            plt.close(fig)
        imageio.mimsave(f'peristalsis_plot_{substrate}.gif', frames, duration=1 / 60)