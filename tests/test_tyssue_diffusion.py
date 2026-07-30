"""Tests for the tyssue-backed diffusion setup.

`tyssue` is an optional dependency, so the whole module is skipped when it is
not installed; deselect it explicitly with ``-m "not tyssue"``.
"""
import pytest

tyssue = pytest.importorskip("tyssue", reason="needs the optional `tyssue` dependency")

from tyssue import Sheet, SheetGeometry  # noqa: E402

from spatial_transport.processes.tyssue_diffusion import (  # noqa: E402
    get_tyssue_edges, generate_tyssue_environments, run_tyssue_diffusion)

pytestmark = pytest.mark.tyssue

SUBSTRATES = {"glucose": 0.06, "acetate": 0.12}


@pytest.fixture(scope="module")
def sheet():
    sheet = Sheet.planar_sheet_3d("sheet", nx=6, ny=6, distx=1, disty=1, noise=0.0)
    sheet.sanitize(trim_borders=True)
    SheetGeometry.update_all(sheet)
    return sheet


def test_every_face_becomes_a_compartment(sheet):
    compartments = generate_tyssue_environments(sheet, SUBSTRATES)

    assert set(compartments) == {f"{face}" for face in sheet.face_df.index}


def test_compartment_positions_come_from_the_face_centres(sheet):
    """Regression test: this indexed `face_df.loc[face_df.loc == i]`, which
    raises `KeyError` on current pandas rather than returning a centre."""
    compartments = generate_tyssue_environments(sheet, SUBSTRATES)

    face = sheet.face_df.index[0]
    expected = [float(sheet.face_df.loc[face, dim]) for dim in sheet.coords]
    assert compartments[f"{face}"]["position"] == expected


def test_compartment_volume_is_the_prism_over_the_face(sheet):
    compartments = generate_tyssue_environments(sheet, SUBSTRATES)

    height = sheet.vert_df.loc[0]["basal_shift"]
    for key, compartment in compartments.items():
        area = sheet.face_df.loc[int(key), "area"]
        assert compartment["Shared Environment"]["volume"] == pytest.approx(area * height)


def test_counts_and_concentrations_are_consistent(sheet):
    compartments = generate_tyssue_environments(sheet, SUBSTRATES)

    for compartment in compartments.values():
        environment = compartment["Shared Environment"]
        for substrate, count in environment["counts"].items():
            assert environment["concentrations"][substrate] == pytest.approx(
                count / environment["volume"])


def test_edges_join_faces_that_share_a_boundary(sheet):
    edges = get_tyssue_edges(sheet)
    faces = {f"{face}" for face in sheet.face_df.index}

    assert edges, "an interior sheet has neighbouring faces"
    for edge in edges.values():
        first, second = edge["neighbors"]
        assert first in faces and second in faces
        assert first != second
        assert edge["surface_area"] > 0


def test_diffusion_on_a_sheet_conserves_mass(core, sheet):
    results = run_tyssue_diffusion(core, sheet, SUBSTRATES)

    def total(result):
        return sum(
            compartment["Shared Environment"]["counts"]["glucose"]
            for compartment in result["compartments"].values())

    totals = [total(result) for result in results]
    assert all(value == pytest.approx(totals[0]) for value in totals)


def test_diffusion_on_a_sheet_reduces_the_spread(core, sheet):
    results = run_tyssue_diffusion(core, sheet, SUBSTRATES)

    def spread(result):
        values = [
            compartment["Shared Environment"]["concentrations"]["glucose"]
            for compartment in result["compartments"].values()]
        return max(values) - min(values)

    assert spread(results[-1]) < spread(results[0])
