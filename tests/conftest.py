"""Shared fixtures and helpers for the spatial-transport test suite.

Grids are built with `uniform_grid`, which seeds every compartment with the same
counts, so any change a process makes to the field is entirely attributable to
that process rather than to the random initial conditions `add_shared_environments`
produces.
"""
import pytest
from process_bigraph import Composite, allocate_core
from process_bigraph.emitter import emitter_from_wires, gather_emitter_results

from spatial_transport.data_types import register_types
from spatial_transport.utils import (
    generate_voxels, get_regular_edges, detect_boundary_positions)

COMPARTMENTS = "Compartments"
EDGES = "Edges"
SHARED_ENVIRONMENT = "Shared Environment"

SPACING = [1.0, 1.0, 1.0]


@pytest.fixture(scope="session")
def core():
    """A core with the spatial-transport (and cdFBA) types and processes registered."""
    return register_types(allocate_core())


def uniform_grid(dims, spacing=SPACING, counts=None, num_dims=2, periodic=False):
    """Build a `(compartments, edges)` pair with identical counts everywhere.

    `counts` maps substrate name to the count placed in *every* compartment;
    tests then override individual compartments to set up a gradient.
    """
    counts = counts or {"glucose": 1.0}
    volume = spacing[0] * spacing[1] * spacing[2]

    compartments = generate_voxels(dims=dims, spacing=spacing)
    for compartment in compartments.values():
        compartment[SHARED_ENVIRONMENT] = {
            "volume": volume,
            "counts": dict(counts),
            "concentrations": {key: value / volume for key, value in counts.items()},
        }
    compartments = detect_boundary_positions(
        compartments, num_dims=num_dims, spacing=spacing)
    edges = get_regular_edges(compartments, periodic=periodic, spacing=spacing)

    return compartments, edges


def set_counts(compartments, compartment_id, counts):
    """Set counts on one compartment, keeping concentrations consistent."""
    environment = compartments[compartment_id][SHARED_ENVIRONMENT]
    volume = environment["volume"]
    for substrate, count in counts.items():
        environment["counts"][substrate] = count
        environment["concentrations"][substrate] = count / volume
    return compartments


def spatial_emitter(**extra_wires):
    """An emitter wired to global time and the compartment field."""
    wires = {
        "global_time": ["global_time"],
        "compartments": [COMPARTMENTS],
    }
    wires.update(extra_wires)
    return emitter_from_wires(wires)


def run_spec(core, spec, duration):
    """Run a spec as a Composite and return `(sim, emitted_results)`."""
    sim = Composite({"state": spec}, core=core)
    sim.run(duration)
    return sim, gather_emitter_results(sim)[("emitter",)]


def total_counts(result, substrate):
    """Sum one substrate over every compartment in an emitted result."""
    return sum(
        compartment[SHARED_ENVIRONMENT]["counts"][substrate]
        for compartment in result["compartments"].values())


def profile(result, substrate, field="concentrations"):
    """Per-compartment values, ordered by ascending x position."""
    compartments = sorted(
        result["compartments"].values(), key=lambda c: c["position"][0])
    return [c[SHARED_ENVIRONMENT][field][substrate] for c in compartments]
