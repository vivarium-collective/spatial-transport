"""Tests for the SimpleAdvection, DynamicAdvection and Peristalsis processes."""
import pytest

from spatial_transport.processes.advection import (
    ADVECTION, get_simple_advection_spec, get_dynamic_advection_spec,
    get_peristalsis_spec, get_advection_store)
from tests.conftest import (
    COMPARTMENTS, EDGES, uniform_grid, set_counts, spatial_emitter, run_spec,
    total_counts, profile)


def simple_spec(compartments, edges, substrates, advection,
                boundary="default", spacing=None, interval=0.1):
    return {
        COMPARTMENTS: compartments,
        EDGES: edges,
        "Simple Advection": get_simple_advection_spec(
            spacing=spacing or [1.0, 1.0, 1.0],
            substrates=substrates,
            advection=advection,
            boundary=boundary,
            interval=interval),
        "emitter": spatial_emitter(),
    }


def peristalsis_spec(compartments, edges, substrates, boundary="default",
                     spacing=None, interval=0.1, amplitude=3.0, velocity=2.0,
                     wavelength=2.0, period=5.0, direction=None):
    return {
        COMPARTMENTS: compartments,
        EDGES: edges,
        ADVECTION: get_advection_store(edges),
        "Dynamic Advection": get_dynamic_advection_spec(
            spacing=spacing or [1.0, 1.0, 1.0],
            substrates=substrates,
            boundary=boundary,
            interval=interval),
        "Peristalsis": get_peristalsis_spec(
            amplitude=amplitude,
            velocity=velocity,
            wavelength=wavelength,
            period=period,
            direction=direction or [1.0, 0.0, 0.0],
            interval=interval),
        "emitter": spatial_emitter(advection=[ADVECTION]),
    }


def centre_of_mass(result, substrate):
    """Mass-weighted mean x position of a substrate across the field."""
    compartments = result["compartments"].values()
    total = sum(c["Shared Environment"]["counts"][substrate] for c in compartments)
    return sum(
        c["position"][0] * c["Shared Environment"]["counts"][substrate]
        for c in compartments) / total


#==================
# SimpleAdvection
#==================
def test_simple_advection_conserves_mass(core):
    compartments, edges = uniform_grid(dims=[6, 1, 0], counts={"glucose": 2.0})
    set_counts(compartments, "1", {"glucose": 8.0})

    _, results = run_spec(
        core, simple_spec(compartments, edges, ["glucose"], [0.5, 0.0, 0.0]), 2)

    totals = [total_counts(result, "glucose") for result in results]
    assert all(total == pytest.approx(totals[0]) for total in totals)


def test_advection_carries_material_along_the_velocity_vector(core):
    compartments, edges = uniform_grid(dims=[6, 1, 0], counts={"glucose": 0.0})
    set_counts(compartments, "1", {"glucose": 8.0})

    _, results = run_spec(
        core, simple_spec(compartments, edges, ["glucose"], [0.5, 0.0, 0.0]), 2)

    assert centre_of_mass(results[-1], "glucose") > centre_of_mass(results[0], "glucose")


def test_a_negative_velocity_reverses_the_transport(core):
    compartments, edges = uniform_grid(dims=[6, 1, 0], counts={"glucose": 0.0})
    set_counts(compartments, "4", {"glucose": 8.0})

    _, results = run_spec(
        core, simple_spec(compartments, edges, ["glucose"], [-0.5, 0.0, 0.0]), 2)

    assert centre_of_mass(results[-1], "glucose") < centre_of_mass(results[0], "glucose")


def test_zero_velocity_leaves_the_field_alone(core):
    compartments, edges = uniform_grid(dims=[4, 1, 0], counts={"glucose": 1.0})
    set_counts(compartments, "0", {"glucose": 6.0})

    _, results = run_spec(
        core, simple_spec(compartments, edges, ["glucose"], [0.0, 0.0, 0.0]), 2)

    assert profile(results[-1], "glucose") == pytest.approx(
        profile(results[0], "glucose"))


def test_a_velocity_orthogonal_to_the_grid_axis_moves_nothing(core):
    """A y-velocity on a line of voxels laid out along x has no face to cross."""
    compartments, edges = uniform_grid(dims=[4, 1, 0], counts={"glucose": 1.0})
    set_counts(compartments, "0", {"glucose": 6.0})

    _, results = run_spec(
        core, simple_spec(compartments, edges, ["glucose"], [0.0, 0.5, 0.0]), 2)

    assert profile(results[-1], "glucose") == pytest.approx(
        profile(results[0], "glucose"))


def test_advection_is_upwind(core):
    """Material crosses a face from whichever side the flow comes from.

    With flow in +x and everything concentrated downstream, the upstream
    compartment must stay empty rather than being fed backwards.
    """
    compartments, edges = uniform_grid(dims=[3, 1, 0], counts={"glucose": 0.0})
    set_counts(compartments, "2", {"glucose": 6.0})

    _, results = run_spec(
        core, simple_spec(compartments, edges, ["glucose"], [0.5, 0.0, 0.0]), 1)

    assert profile(results[-1], "glucose")[:2] == pytest.approx([0.0, 0.0])


def test_periodic_advection_conserves_mass_and_wraps_around(core):
    compartments, edges = uniform_grid(
        dims=[4, 1, 0], counts={"glucose": 0.0}, periodic=True)
    set_counts(compartments, "3", {"glucose": 8.0})

    _, results = run_spec(
        core,
        simple_spec(compartments, edges, ["glucose"], [0.5, 0.0, 0.0],
                    boundary="periodic"),
        1)

    totals = [total_counts(result, "glucose") for result in results]
    assert all(total == pytest.approx(totals[0]) for total in totals)

    # the last compartment feeds the first one across the wrapping edge
    assert profile(results[-1], "glucose")[0] > 0.0


#===============
# Peristalsis
#===============
def test_peristalsis_writes_a_vector_for_every_edge(core):
    compartments, edges = uniform_grid(dims=[6, 1, 0])

    _, results = run_spec(
        core, peristalsis_spec(compartments, edges, ["glucose"]), 1)

    advection = results[-1]["advection"]
    assert set(advection) == set(edges)
    assert all(len(vector) == 3 for vector in advection.values())


def test_peristalsis_vectors_stay_three_dimensional_over_many_steps(core):
    """Regression test for the `overwrite` wrapper on the advection store.

    Without it the list updates concatenate, so each edge's vector grows by
    three entries per step and the dot product against the face normal fails.
    """
    compartments, edges = uniform_grid(dims=[4, 1, 0])

    _, results = run_spec(
        core, peristalsis_spec(compartments, edges, ["glucose"]), 3)

    assert all(
        len(vector) == 3
        for result in results
        for vector in result["advection"].values())


def test_peristalsis_points_along_the_configured_direction(core):
    compartments, edges = uniform_grid(dims=[6, 1, 0])

    _, results = run_spec(
        core,
        peristalsis_spec(compartments, edges, ["glucose"], direction=[1.0, 0.0, 0.0]),
        2)

    vectors = [v for result in results for v in result["advection"].values()]
    assert all(v[1] == 0.0 and v[2] == 0.0 for v in vectors)
    assert any(v[0] > 0.0 for v in vectors), "the wave reaches the field"


def test_the_peristaltic_wave_travels_downstream(core):
    """The Gaussian pulse moves along +x, so its peak edge index increases."""
    compartments, edges = uniform_grid(dims=[12, 1, 0])

    _, results = run_spec(
        core,
        peristalsis_spec(compartments, edges, ["glucose"], velocity=2.0, period=20.0),
        4)

    def peak_position(result):
        edge_id = max(result["advection"], key=lambda k: result["advection"][k][0])
        return edges[edge_id]["position"][0]

    assert peak_position(results[-1]) > peak_position(results[len(results) // 4])


def test_amplitude_scales_the_advection_vectors(core):
    compartments, edges = uniform_grid(dims=[8, 1, 0])

    def peak_speed(amplitude):
        _, results = run_spec(
            core,
            peristalsis_spec(compartments, edges, ["glucose"], amplitude=amplitude),
            2)
        return max(v[0] for r in results for v in r["advection"].values())

    assert peak_speed(6.0) == pytest.approx(2 * peak_speed(3.0))


#====================
# DynamicAdvection
#====================
def test_dynamic_advection_conserves_mass(core):
    compartments, edges = uniform_grid(dims=[8, 1, 0], counts={"glucose": 2.0})

    _, results = run_spec(
        core, peristalsis_spec(compartments, edges, ["glucose"]), 3)

    totals = [total_counts(result, "glucose") for result in results]
    assert all(total == pytest.approx(totals[0]) for total in totals)


def test_dynamic_advection_transports_with_the_wave(core):
    compartments, edges = uniform_grid(dims=[10, 1, 0], counts={"glucose": 0.0})
    set_counts(compartments, "1", {"glucose": 10.0})

    _, results = run_spec(
        core, peristalsis_spec(compartments, edges, ["glucose"]), 3)

    assert centre_of_mass(results[-1], "glucose") > centre_of_mass(results[0], "glucose")


def test_a_zero_advection_field_moves_nothing(core):
    """With amplitude 0 the Peristalsis step writes zero vectors everywhere."""
    compartments, edges = uniform_grid(dims=[6, 1, 0], counts={"glucose": 1.0})
    set_counts(compartments, "0", {"glucose": 9.0})

    _, results = run_spec(
        core, peristalsis_spec(compartments, edges, ["glucose"], amplitude=0.0), 2)

    assert profile(results[-1], "glucose") == pytest.approx(
        profile(results[0], "glucose"))
