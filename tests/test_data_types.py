"""Tests for the spatial-transport schema types.

`volumetric` is not defined here — it is imported from cdFBA so both packages
apply counts/concentration updates identically; the tests below check that the
shared type arrives intact and that the spatial types layered on top of it
carry the fields the processes read.
"""
import pytest

from cdFBA.data_types import Volumetric
from spatial_transport.data_types import Volumetric as SpatialVolumetric

COMPARTMENT = {
    "Shared Environment": {
        "volume": 2.0,
        "counts": {"glucose": 10.0},
        "concentrations": {"glucose": 5.0},
    },
    "position": [1.0, 0.5, 0.5],
    "boundaries": ["x_min"],
}

EDGE = {
    "neighbors": ["0", "1"],
    "surface_area": 1.0,
    "position": [1.5, 0.5, 0.5],
    "periodic": True,
}


def realize(core, type_name, state):
    _, realized, _ = core.realize(type_name, state)
    return realized


def apply_update(core, type_name, state, update):
    schema, realized, _ = core.realize(type_name, state)
    applied = core.apply(schema, realized, update)
    return applied[0] if isinstance(applied, tuple) else applied


def test_registered_types_are_available(core):
    for type_name in ("volumetric", "edge_type", "compartment", "advection_vector"):
        assert core.access(type_name) is not None


def test_volumetric_is_the_cdfba_type(core):
    """The type is imported rather than redefined, so it must be the same class."""
    assert SpatialVolumetric is Volumetric
    assert isinstance(core.access("volumetric"), Volumetric)


def test_cdfba_types_are_registered_alongside_the_spatial_ones(core):
    """register_types chains cdFBA's, so a composite needs only one call."""
    for type_name in ("bounds", "threshold", "dfba_changes"):
        assert core.access(type_name) is not None
    for link in ("dFBA", "UpdateEnvironment", "Injector"):
        assert core.link_registry.get(link) is not None


def test_spatial_processes_are_registered(core):
    for link in ("SimpleDiffusion", "SimpleAdvection", "DynamicAdvection", "Peristalsis"):
        assert core.link_registry.get(link) is not None


def test_realize_keeps_all_compartment_fields(core):
    state = realize(core, "compartment", COMPARTMENT)

    assert state["position"] == [1.0, 0.5, 0.5]
    assert state["boundaries"] == ["x_min"]
    assert state["Shared Environment"]["counts"] == {"glucose": 10.0}
    assert state["Shared Environment"]["volume"] == 2.0


def test_realize_keeps_the_edge_periodic_flag(core):
    """`periodic` and `position` are read by the advection processes.

    Both are now declared on the type rather than riding along as undeclared
    keys, so they are typed (`periodic` really is a boolean) and the port
    projection is explicit about what a process may read.
    """
    state = realize(core, "edge_type", EDGE)

    assert state["periodic"] is True
    assert state["position"] == [1.5, 0.5, 0.5]
    assert state["neighbors"] == ["0", "1"]
    assert state["surface_area"] == 1.0


def test_compartment_counts_update_recomputes_concentrations(core):
    result = apply_update(
        core, "compartment", COMPARTMENT,
        {"Shared Environment": {"counts": {"glucose": 2.0}}})

    # counts apply additively, concentrations are derived from counts / volume
    assert result["Shared Environment"]["counts"]["glucose"] == 12.0
    assert result["Shared Environment"]["concentrations"]["glucose"] == 6.0


def test_compartment_position_survives_a_counts_update(core):
    result = apply_update(
        core, "compartment", COMPARTMENT,
        {"Shared Environment": {"counts": {"glucose": 1.0}}})

    assert result["position"] == [1.0, 0.5, 0.5]
    assert result["boundaries"] == ["x_min"]


def test_advection_vector_overwrites_instead_of_accumulating(core):
    """Peristalsis recomputes the whole field each step, so updates replace.

    A plain `list[float]` would concatenate the update onto the stored vector,
    growing it without bound; `overwrite` is what makes the store hold a
    3-vector after every step.
    """
    result = apply_update(core, "advection_vector", [1.0, 0.0, 0.0], [3.0, 2.0, 1.0])

    assert result == [3.0, 2.0, 1.0]


def test_advection_field_overwrites_per_edge(core):
    result = apply_update(
        core, "map[advection_vector]",
        {"1": [1.0, 0.0, 0.0], "2": [1.0, 0.0, 0.0]},
        {"1": [5.0, 0.0, 0.0]})

    assert result["1"] == [5.0, 0.0, 0.0]
    assert result["2"] == [1.0, 0.0, 0.0], "untouched edges keep their vector"


def test_plain_list_would_accumulate(core):
    """Contrast case pinning down why `advection_vector` wraps in `overwrite`."""
    assert apply_update(core, "list[float]", [1.0, 0.0, 0.0], [3.0, 2.0, 1.0]) == [
        1.0, 0.0, 0.0, 3.0, 2.0, 1.0]
