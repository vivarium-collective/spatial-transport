"""Whole-composite runs combining several transport processes, and the cdFBA
composite-inside-a-compartment arrangement the peristalsis demo is built on.

Tests marked `cobra` build dFBA processes from cobra's bundled `textbook`
model; deselect them with ``-m "not cobra"``. The `bigg` marked test pulls a
genome-scale model over the network.
"""
import pytest
from process_bigraph.emitter import emitter_from_wires

from spatial_transport.processes.diffusion import get_simple_diffusion_spec
from spatial_transport.processes.advection import (
    ADVECTION, get_dynamic_advection_spec, get_peristalsis_spec, get_advection_store)
from spatial_transport.utils import (
    generate_voxels, get_regular_edges, generate_simple_cdfba_composite)
from tests.conftest import (
    COMPARTMENTS, EDGES, SHARED_ENVIRONMENT, SPACING,
    uniform_grid, set_counts, spatial_emitter, run_spec, total_counts, profile)

EXCHANGES = ["EX_glc__D_e", "EX_ac_e"]


def transport_spec(compartments, edges, substrates, interval=0.1):
    """Diffusion and peristaltic advection acting on the same field."""
    return {
        COMPARTMENTS: compartments,
        EDGES: edges,
        ADVECTION: get_advection_store(edges),
        "Simple Diffusion": get_simple_diffusion_spec(
            substrates={substrate: 0.05 for substrate in substrates},
            interval=interval),
        "Dynamic Advection": get_dynamic_advection_spec(
            spacing=SPACING, substrates=substrates, boundary="default",
            interval=interval),
        "Peristalsis": get_peristalsis_spec(
            amplitude=3.0, velocity=2.0, wavelength=2.0, period=5.0,
            direction=[1.0, 0.0, 0.0], interval=interval),
        "emitter": spatial_emitter(),
    }


def test_diffusion_and_advection_together_conserve_mass(core):
    compartments, edges = uniform_grid(dims=[10, 1, 0], counts={"glucose": 0.0})
    set_counts(compartments, "2", {"glucose": 10.0})

    _, results = run_spec(core, transport_spec(compartments, edges, ["glucose"]), 3)

    totals = [total_counts(result, "glucose") for result in results]
    assert all(total == pytest.approx(totals[0]) for total in totals)


def test_diffusion_and_advection_spread_material_further_than_either_alone(core):
    compartments, edges = uniform_grid(dims=[10, 1, 0], counts={"glucose": 0.0})
    set_counts(compartments, "2", {"glucose": 10.0})

    _, results = run_spec(core, transport_spec(compartments, edges, ["glucose"]), 3)

    occupied = sum(1 for value in profile(results[-1], "glucose") if value > 1e-9)
    assert occupied > 1, "material has left the compartment it started in"


def test_emitter_records_one_result_per_step(core):
    compartments, edges = uniform_grid(dims=[4, 1, 0])
    spec = {
        COMPARTMENTS: compartments,
        EDGES: edges,
        "Simple Diffusion": get_simple_diffusion_spec(
            substrates={"glucose": 0.05}, interval=1.0),
        "emitter": spatial_emitter(),
    }

    _, results = run_spec(core, spec, 10)

    assert len(results) == 11
    assert [result["global_time"] for result in results] == [
        float(step) for step in range(11)]


#=========================================
# cdFBA composites inside compartments
#=========================================
def cdfba_field(model_dict, dims=(4, 1, 0), volume=1.0):
    """A line of voxels, each holding its own cdFBA composite."""
    voxels = generate_voxels(dims=list(dims), spacing=SPACING)
    edges = get_regular_edges(voxels, periodic=False, spacing=SPACING)
    compartments = generate_simple_cdfba_composite(
        voxels, model_dict, EXCHANGES, volume,
        sub_range=(5, 5), bio_range=(0.001, 0.001))
    return compartments, edges


@pytest.mark.cobra
def test_dfba_composites_grow_inside_advected_compartments(core):
    """The spatial types and the cdFBA types have to coexist in one core."""
    model_dict = {"ecoli": "textbook"}
    compartments, edges = cdfba_field(model_dict)

    spec = transport_spec(compartments, edges, ["D-Glucose", "Acetate"])
    _, results = run_spec(core, spec, 1.0)

    def biomass(result):
        return sum(
            c[SHARED_ENVIRONMENT]["counts"]["ecoli"]
            for c in result["compartments"].values())

    assert biomass(results[-1]) > biomass(results[0]), "the culture should grow"


@pytest.mark.cobra
def test_transport_redistributes_glucose_between_dfba_compartments(core):
    model_dict = {"ecoli": "textbook"}
    compartments, edges = cdfba_field(model_dict, dims=(6, 1, 0))
    for key in list(compartments)[1:]:
        environment = compartments[key][SHARED_ENVIRONMENT]
        environment["counts"]["D-Glucose"] = 0.0
        environment["concentrations"]["D-Glucose"] = 0.0

    spec = transport_spec(compartments, edges, ["D-Glucose", "Acetate"])
    _, results = run_spec(core, spec, 1.0)

    downstream = profile(results[-1], "D-Glucose")[1:]
    assert any(value > 0.0 for value in downstream), (
        "glucose seeded only upstream must reach the other compartments")


@pytest.mark.cobra
def test_each_compartment_gets_its_own_species_stores():
    """`generate_simple_cdfba_composite` must not alias stores between voxels.

    It used to shallow-copy the base cdFBA spec, so every voxel shared one
    `Species` dict and one `dFBA Results` dict and a per-voxel edit landed in
    all of them.
    """
    compartments, _ = cdfba_field({"ecoli": "textbook"}, dims=(3, 1, 0))

    for store in ("Species", "dFBA Results", "Shared Environment"):
        assert compartments["0"][store] is not compartments["1"][store]

    compartments["0"]["dFBA Results"]["ecoli"]["Acetate"] = 99.0
    assert compartments["1"]["dFBA Results"]["ecoli"]["Acetate"] != 99.0


@pytest.mark.cobra
def test_a_nested_injector_only_affects_its_own_compartment(core):
    """A cdFBA Injector placed in one voxel, as `Notebooks/peristalsis.py` does.

    This works on its own, but adding any process with a `map[compartment]`
    port replicates the Injector link into every compartment and then drops its
    updates entirely, so nothing is ever injected. The homogenizing happens
    when the Map port projection resolves one value schema across every key;
    it is a framework-level issue, not something this package configures.
    """
    from cdFBA.utils import get_injector_spec

    compartments, edges = cdfba_field({"ecoli": "textbook"}, dims=(3, 1, 0))

    injector = get_injector_spec(
        {"injection_params": {"Acetate": {"amount": 10.0, "interval": 0.5}}},
        interval=0.1)
    injector["inputs"]["global_time"] = ["..", "..", "global_time"]
    compartments["0"]["Injector"] = injector

    spec = transport_spec(compartments, edges, ["D-Glucose", "Acetate"])
    sim, results = run_spec(core, spec, 1.0)

    injected_into = [
        key for key, compartment in sim.state[COMPARTMENTS].items()
        if "Injector" in compartment]
    acetate = [total_counts(result, "Acetate") for result in results]

    assert injected_into == ["0"], "the Injector belongs to one compartment only"
    assert max(acetate) > acetate[0], "the injector adds acetate to the field"


# xfail rather than skip so the day the framework stops homogenizing Map value
# schemas, this reports as an unexpected pass and the demo can be re-enabled.
test_a_nested_injector_only_affects_its_own_compartment = pytest.mark.xfail(
    reason="a map[compartment] port replicates per-compartment links across "
           "every key and drops their updates",
    strict=True,
)(test_a_nested_injector_only_affects_its_own_compartment)


@pytest.mark.bigg
def test_bigg_model_runs_in_a_spatial_field(core):
    """End-to-end run on the genome-scale model used by the original demo."""
    model_dict = {"E.coli": "iAF1260"}
    compartments, edges = cdfba_field(model_dict, dims=(3, 1, 0))

    spec = transport_spec(compartments, edges, ["D-Glucose", "Acetate"])
    _, results = run_spec(core, spec, 1.0)

    assert len(results) == 11
    totals = [total_counts(result, "D-Glucose") for result in results]
    assert totals[-1] < totals[0], "the culture should draw down glucose"
