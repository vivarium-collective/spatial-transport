"""Tests for the SimpleDiffusion process."""
import pytest

from spatial_transport.processes.diffusion import get_simple_diffusion_spec
from tests.conftest import (
    COMPARTMENTS, EDGES, SHARED_ENVIRONMENT,
    uniform_grid, set_counts, spatial_emitter, run_spec, total_counts, profile)


def diffusion_spec(compartments, edges, substrates, interval=0.1):
    return {
        COMPARTMENTS: compartments,
        EDGES: edges,
        "Simple Diffusion": get_simple_diffusion_spec(
            substrates=substrates, interval=interval),
        "emitter": spatial_emitter(),
    }


def test_diffusion_conserves_mass(core):
    compartments, edges = uniform_grid(dims=[4, 1, 0])
    set_counts(compartments, "0", {"glucose": 10.0})

    _, results = run_spec(
        core, diffusion_spec(compartments, edges, {"glucose": 0.1}), 2)

    totals = [total_counts(result, "glucose") for result in results]
    assert all(total == pytest.approx(totals[0]) for total in totals)


def test_diffusion_flattens_a_gradient(core):
    compartments, edges = uniform_grid(dims=[4, 1, 0], counts={"glucose": 0.0})
    set_counts(compartments, "0", {"glucose": 8.0})

    _, results = run_spec(
        core, diffusion_spec(compartments, edges, {"glucose": 0.1}), 20)

    start, end = profile(results[0], "glucose"), profile(results[-1], "glucose")
    assert start == [8.0, 0.0, 0.0, 0.0]
    assert max(end) - min(end) < max(start) - min(start)
    assert end[0] > end[1] > end[2] > end[3], "the profile stays monotonic"


def test_diffusion_moves_material_down_the_gradient(core):
    compartments, edges = uniform_grid(dims=[2, 1, 0], counts={"glucose": 0.0})
    set_counts(compartments, "0", {"glucose": 4.0})

    _, results = run_spec(
        core, diffusion_spec(compartments, edges, {"glucose": 0.1}), 1)

    first, second = profile(results[-1], "glucose")
    assert first < 4.0 and second > 0.0


def test_a_uniform_field_does_not_change(core):
    compartments, edges = uniform_grid(dims=[3, 3, 0], counts={"glucose": 2.0})

    _, results = run_spec(
        core, diffusion_spec(compartments, edges, {"glucose": 0.5}), 5)

    assert profile(results[-1], "glucose") == pytest.approx(
        profile(results[0], "glucose"))


def test_higher_diffusivity_equilibrates_faster(core):
    def spread_after_run(diffusivity):
        compartments, edges = uniform_grid(dims=[4, 1, 0], counts={"glucose": 0.0})
        set_counts(compartments, "0", {"glucose": 8.0})
        _, results = run_spec(
            core, diffusion_spec(compartments, edges, {"glucose": diffusivity}), 5)
        values = profile(results[-1], "glucose")
        return max(values) - min(values)

    assert spread_after_run(0.2) < spread_after_run(0.05)


def test_substrates_diffuse_independently(core):
    compartments, edges = uniform_grid(
        dims=[4, 1, 0], counts={"glucose": 0.0, "acetate": 0.0})
    set_counts(compartments, "0", {"glucose": 8.0, "acetate": 8.0})

    _, results = run_spec(
        core,
        diffusion_spec(compartments, edges, {"glucose": 0.2, "acetate": 0.02}),
        5)

    glucose, acetate = profile(results[-1], "glucose"), profile(results[-1], "acetate")
    assert max(glucose) - min(glucose) < max(acetate) - min(acetate)


def test_a_substrate_left_out_of_the_config_is_untouched(core):
    compartments, edges = uniform_grid(
        dims=[3, 1, 0], counts={"glucose": 0.0, "inert": 1.0})
    set_counts(compartments, "0", {"glucose": 6.0, "inert": 5.0})

    _, results = run_spec(
        core, diffusion_spec(compartments, edges, {"glucose": 0.2}), 5)

    assert profile(results[-1], "inert") == pytest.approx([5.0, 1.0, 1.0])


def test_periodic_boundaries_conserve_mass_and_wrap(core):
    compartments, edges = uniform_grid(
        dims=[4, 1, 0], counts={"glucose": 0.0}, periodic=True)
    set_counts(compartments, "0", {"glucose": 8.0})

    _, results = run_spec(
        core, diffusion_spec(compartments, edges, {"glucose": 0.1}), 2)

    totals = [total_counts(result, "glucose") for result in results]
    assert all(total == pytest.approx(totals[0]) for total in totals)

    # on a ring the far compartment is a neighbour, so it fills as fast as
    # its mirror image on the other side
    values = profile(results[-1], "glucose")
    assert values[1] == pytest.approx(values[3])


def test_surface_area_scales_the_flux(core):
    """A wider face between the same two compartments moves more material.

    Both runs start from the same concentration gradient (counts are scaled by
    the voxel volume), so the only difference is the area of the shared face.
    The flux is a count, not a concentration: widening the face by 4x moves 4x
    the material but into a voxel that is also 4x larger, which would leave the
    concentration unchanged.
    """
    def transferred(spacing):
        compartments, edges = uniform_grid(
            dims=[2, 1, 0], spacing=spacing, counts={"glucose": 0.0})
        volume = spacing[0] * spacing[1] * spacing[2]
        set_counts(compartments, "0", {"glucose": 4.0 * volume})
        _, results = run_spec(
            core, diffusion_spec(compartments, edges, {"glucose": 0.05}), 1)
        return profile(results[-1], "glucose", field="counts")[1]

    wide, narrow = transferred([1.0, 4.0, 1.0]), transferred([1.0, 1.0, 1.0])
    assert wide == pytest.approx(4 * narrow)
