from pprint import pprint
import random
import numpy as np

from process_bigraph import Process, Composite, ProcessTypes
from process_bigraph.emitter import emitter_from_wires, gather_emitter_results

from spatial_transport.processes.diffusion import SimpleDiffusion
from spatial_transport import register_types
from spatial_transport.processes.diffusion import get_simple_diffusion_spec

from tyssue import Sheet, SheetGeometry
from tyssue.draw import sheet_view
from tyssue import config

import io
import imageio.v2 as imageio
import matplotlib.pyplot as plt
from matplotlib import colormaps as cm
from matplotlib.colors import BoundaryNorm

def get_tyssue_edges(sheet):
    sheet.get_extra_indices()
    height = sheet.vert_df.loc[0]["basal_shift"]
    lengths = list(sheet.edge_df.loc[sheet.east_edges]["length"])
    faces = list(sheet.edge_df.loc[sheet.east_edges]["face"])
    opposites = list(sheet.edge_df.loc[sheet.east_edges]["opposite"])
    opposite_faces = list(sheet.edge_df.loc[opposites]["face"])
    neighbors = [[f"{faces[i]}", f"{opposite_faces[i]}"] for i in range(len(faces))]
    edges = {
        f"{i}": {
            "neighbors": neighbors[i],
            "surface_area": float(lengths[i]*height)
        }
        for i in range(len(faces))
    }
    return edges

def generate_tyssue_environments(sheet, substrates):
    comps = {
        f"{i}" : {
            "position": [float(sheet.face_df.loc[i][dim]) for dim in sheet.coords]
        }
        for i in range(len(sheet.face_df))
    }
    height = sheet.vert_df.loc[0]["basal_shift"]
    for key in comps.keys():
        area = sheet.face_df.loc[int(key)]["area"]
        volume = float(area*height)
        comps[key]['Shared Environment'] = {}
        comps[key]['Shared Environment']['volume'] = volume
        comps[key]['Shared Environment']['counts'] = {}
        comps[key]['Shared Environment']['concentrations'] = {}
        for substrate in substrates:
            concentration = random.uniform(0, 10)
            comps[key]['Shared Environment']['counts'][substrate] = concentration*volume
            comps[key]['Shared Environment']['concentrations'][substrate] = concentration
    compartments = comps
    return compartments

def run_tyssue_diffusion(core, sheet, substrates):
    edges = get_tyssue_edges(sheet)
    compartments = generate_tyssue_environments(sheet, substrates)
    spec = {}
    spec["Simple Diffusion"] = get_simple_diffusion_spec(substrates, interval=0.1)
    spec["Compartments"] = compartments
    spec["Edges"] = edges
    spec["emitter"] = emitter_from_wires({
        "global_time": ["global_time"],
        'compartments': ['Compartments'],
    })
    pprint(spec)
    sim = Composite(
        {
            "state": spec,
        },
        core=core
    )
    sim.run(20)
    results = gather_emitter_results(sim)[("emitter",)]
    return results

def static_sheet_video_2d(results, sheet, substrate, vmax):
    """
    Parameters:
        results: dict, Vivarium results
        sheet: Sheet, tyssue Sheet object
        substrate: str, substrate name
        vmax: float, maximum substrate concentration
    """
    frames = []
    boundaries = [float(i) for i in np.arange(0, vmax+1, vmax/100)]
    norm = BoundaryNorm(boundaries=boundaries, ncolors=256)
    for result in results:
        cmap = cm.get_cmap('plasma')
        colors = [result["compartments"][i]["Shared Environment"]["concentrations"][substrate] for i in result["compartments"].keys()]
        colors = cmap(norm(colors))
        draw_specs = config.draw.sheet_spec()
        draw_specs["face"]["visible"] = True
        draw_specs["face"]["color"] = colors
        draw_specs["face"]["alpha"] = 1
        fig, ax = sheet_view(sheet, coords = ["x", "y"], **draw_specs)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title(f"t={result["global_time"]:.2f}")
        cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation='vertical')
        cbar.set_label("Glucose Concentration (mmol/L)")
        # Save fig to buffer
        buf = io.BytesIO()
        fig.savefig(buf, format='png')
        buf.seek(0)
        frames.append(imageio.imread(buf))
        plt.close(fig)
        imageio.mimsave('tyssue_diffusion.gif', frames, duration=1 / 30, loop=0)

if __name__ == "__main__":
    sheet = Sheet.planar_sheet_3d("sheet", nx=20, ny=20, distx=1, disty=1, noise=0.1)
    sheet.sanitize(trim_borders=True)
    sheet = sheet.extract_bounding_box(x_boundary=(0.2, 10.2), y_boundary=(0.2, 10.2), coords=['x', 'y', 'z'])
    sheet.sanitize(trim_borders=True)
    SheetGeometry.update_all(sheet)
    substrates = {
        "glucose": 0.06,
        "acetate": 0.12,
    }
    # create the core object
    core = ProcessTypes()
    # register data types
    core = register_types(core)
    core.register_process("SimpleDiffusion", SimpleDiffusion)
    results = run_tyssue_diffusion(core, sheet, substrates)
    pprint(len(results))
    static_sheet_video_2d(results, sheet, "glucose", vmax=10)