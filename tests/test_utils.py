"""Tests for the grid construction helpers in `spatial_transport.utils`."""
import pytest

from spatial_transport.utils import (
    generate_voxels, get_regular_edges, detect_boundary_positions,
    add_shared_environments, generate_shared_environment)


def positions(voxels):
    return [tuple(v["position"]) for v in voxels.values()]


def neighbor_pairs(edges):
    return {tuple(edge["neighbors"]) for edge in edges.values()}


#=================
# generate_voxels
#=================
def test_two_dimensional_grid_is_flat():
    voxels = generate_voxels(dims=[2, 2, 0], spacing=[1, 1, 1])

    assert len(voxels) == 4
    assert positions(voxels) == [
        (0.5, 0.5, 0.5), (0.5, 1.5, 0.5), (1.5, 0.5, 0.5), (1.5, 1.5, 0.5)]


def test_three_dimensional_grid_fills_the_volume():
    voxels = generate_voxels(dims=[2, 2, 2], spacing=[1, 1, 1])

    assert len(voxels) == 8
    assert (0.5, 0.5, 0.5) in positions(voxels)
    assert (1.5, 1.5, 1.5) in positions(voxels)


@pytest.mark.parametrize("dims,expected", [
    ([1, 1, 0], 1),
    ([3, 2, 0], 6),
    ([2, 3, 4], 24),
])
def test_voxel_count_matches_dimensions(dims, expected):
    assert len(generate_voxels(dims=dims, spacing=[1, 1, 1])) == expected


def test_voxels_are_centred_and_spaced_by_spacing():
    voxels = generate_voxels(dims=[3, 1, 0], spacing=[4, 1, 1])

    xs = sorted(position[0] for position in positions(voxels))
    # voxel centres sit half a step in, one full step apart
    assert xs == [2.0, 6.0, 10.0]


def test_anisotropic_spacing_is_applied_per_axis():
    voxels = generate_voxels(dims=[2, 2, 0], spacing=[2, 4, 1])

    assert sorted({position[0] for position in positions(voxels)}) == [1.0, 3.0]
    assert sorted({position[1] for position in positions(voxels)}) == [2.0, 6.0]


def test_zero_z_spacing_in_three_dimensions_is_rejected():
    with pytest.raises(ValueError, match="z-spacing"):
        generate_voxels(dims=[2, 2, 2], spacing=[1, 1, 0])


#===================
# get_regular_edges
#===================
def test_a_line_of_voxels_has_one_fewer_edge_than_voxels():
    voxels = generate_voxels(dims=[4, 1, 0], spacing=[1, 1, 1])
    edges = get_regular_edges(voxels, periodic=False, spacing=[1, 1, 1])

    assert len(edges) == 3
    assert neighbor_pairs(edges) == {("0", "1"), ("1", "2"), ("2", "3")}


def test_grid_edges_are_six_connected_and_deduplicated():
    voxels = generate_voxels(dims=[3, 3, 0], spacing=[1, 1, 1])
    edges = get_regular_edges(voxels, periodic=False, spacing=[1, 1, 1])

    # 3x3 grid: 2 internal edges per row x 3 rows, in each of 2 directions
    assert len(edges) == 12
    assert len(neighbor_pairs(edges)) == len(edges), "each pair appears once"


def test_edges_sit_midway_between_their_neighbors():
    voxels = generate_voxels(dims=[2, 1, 0], spacing=[2, 1, 1])
    edges = get_regular_edges(voxels, periodic=False, spacing=[2, 1, 1])

    (edge,) = edges.values()
    assert edge["position"] == [2.0, 0.5, 0.5]


def test_surface_area_uses_the_orthogonal_face():
    """An x-normal edge is bounded by the y and z spacing, not the x spacing."""
    spacing = [2, 3, 5]
    voxels = generate_voxels(dims=[2, 2, 0], spacing=spacing)
    edges = get_regular_edges(voxels, periodic=False, spacing=spacing)

    areas = set()
    for edge in edges.values():
        first, second = (voxels[key]["position"] for key in edge["neighbors"])
        axis = next(i for i in range(3) if first[i] != second[i])
        areas.add((axis, edge["surface_area"]))

    assert areas == {(0, 3 * 5), (1, 2 * 5)}


def test_non_periodic_edges_are_never_flagged_periodic():
    voxels = generate_voxels(dims=[4, 1, 0], spacing=[1, 1, 1])
    edges = get_regular_edges(voxels, periodic=False, spacing=[1, 1, 1])

    assert not any(edge["periodic"] for edge in edges.values())


def test_periodic_boundaries_add_a_wrapping_edge_per_axis():
    voxels = generate_voxels(dims=[4, 1, 0], spacing=[1, 1, 1])
    edges = get_regular_edges(voxels, periodic=True, spacing=[1, 1, 1])

    wrapped = {tuple(e["neighbors"]) for e in edges.values() if e["periodic"]}
    assert wrapped == {("0", "3")}
    assert len(edges) == 4, "the line closes into a ring"


def test_periodic_grid_wraps_both_axes():
    voxels = generate_voxels(dims=[3, 3, 0], spacing=[1, 1, 1])
    edges = get_regular_edges(voxels, periodic=True, spacing=[1, 1, 1])

    # 12 interior edges plus one wrap per row and per column
    assert len(edges) == 18
    assert sum(edge["periodic"] for edge in edges.values()) == 6


#============================
# detect_boundary_positions
#============================
def test_corner_compartments_touch_two_boundaries():
    voxels = generate_voxels(dims=[3, 3, 0], spacing=[1, 1, 1])
    voxels = detect_boundary_positions(voxels, num_dims=2, spacing=[1, 1, 1])

    boundaries = {key: v["boundaries"] for key, v in voxels.items()}
    assert boundaries["0"] == ["x_min", "y_min"]
    assert boundaries["8"] == ["x_max", "y_max"]
    assert boundaries["4"] == [], "the centre voxel touches nothing"


def test_two_dimensional_detection_ignores_the_z_axis():
    voxels = generate_voxels(dims=[2, 2, 0], spacing=[1, 1, 1])
    voxels = detect_boundary_positions(voxels, num_dims=2, spacing=[1, 1, 1])

    assert not any(
        b.startswith("z") for v in voxels.values() for b in v["boundaries"])


def test_three_dimensional_detection_labels_the_z_faces():
    voxels = generate_voxels(dims=[2, 2, 2], spacing=[1, 1, 1])
    voxels = detect_boundary_positions(voxels, num_dims=3, spacing=[1, 1, 1])

    assert voxels["0"]["boundaries"] == ["x_min", "y_min", "z_min"]


#==========================
# shared environment setup
#==========================
def test_shared_environment_volume_is_the_voxel_volume():
    spacing = [2, 3, 4]
    voxels = generate_voxels(dims=[2, 1, 0], spacing=spacing)
    compartments = add_shared_environments(
        voxels, substrates=["glucose"], spacing=spacing)

    for compartment in compartments.values():
        assert compartment["Shared Environment"]["volume"] == 24


def test_shared_environment_concentrations_match_counts_over_volume():
    spacing = [2, 1, 1]
    voxels = generate_voxels(dims=[3, 1, 0], spacing=spacing)
    compartments = add_shared_environments(
        voxels, substrates=["glucose", "acetate"], spacing=spacing)

    for compartment in compartments.values():
        environment = compartment["Shared Environment"]
        for substrate, count in environment["counts"].items():
            assert environment["concentrations"][substrate] == pytest.approx(
                count / environment["volume"])


def test_shared_environment_counts_stay_within_the_requested_range():
    voxels = generate_voxels(dims=[4, 1, 0], spacing=[1, 1, 1])
    compartments = add_shared_environments(
        voxels, substrates=["glucose"], spacing=[1, 1, 1], min_max=[2, 3])

    counts = [c["Shared Environment"]["counts"]["glucose"] for c in compartments.values()]
    assert all(2 <= count <= 3 for count in counts)


def test_generate_shared_environment_covers_substrates_and_species():
    environment = generate_shared_environment(
        volume=2.0, substrates=["glucose"], species=["E.coli"],
        sub_range=(4, 4), bio_range=(1, 1))

    assert environment["counts"] == {"glucose": 4.0, "E.coli": 1.0}
    assert environment["concentrations"] == {"glucose": 2.0, "E.coli": 0.5}
