import math

import numpy as np
import pytest

from mn_molecular_design_studio.md_interaction_layout import _crossing_count, distribute_labels, reduce_crossings


ANCHORS = np.array([[.582471, -.56718425], [-1.175, .12171825], [1.175, .28241125],
    [.6334895, .02868175], [-.2681115, -.37878475], [-.653911, -.31223275],
    [-.577066, .13325675], [.7399915, .40541025], [.357059, .11339925],
    [.9020005, .56302475], [.582471, -.56718425], [-.653911, -.31223275]])


@pytest.mark.parametrize("count", [1, 2, 4, 12, 16])
def test_equal_arc_spacing_keeps_nodes_on_the_ellipse_without_crowding_the_ends(count):
    anchors = np.array([(math.cos(index), .5 * math.sin(index)) for index in range(count)])
    positions = distribute_labels(anchors, x_radius=2.18, y_radius=1.4)
    assert len(positions) == count
    assert len({tuple(position) for position in positions}) == count
    np.testing.assert_allclose((positions[:, 0] / 2.18) ** 2 + (positions[:, 1] / 1.4) ** 2, 1)
    sample_angles = np.linspace(0, 2 * math.pi, 20001)
    points = np.column_stack((2.18 * np.cos(sample_angles), 1.4 * np.sin(sample_angles)))
    arc = np.concatenate(([0], np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))))
    angles = np.mod(np.arctan2(positions[:, 1] / 1.4, positions[:, 0] / 2.18), 2 * math.pi)
    fractions = np.sort(np.interp(angles, sample_angles, arc) / arc[-1])
    gaps = np.diff(np.append(fractions, fractions[0] + 1))
    np.testing.assert_allclose(gaps, 1 / count, atol=1e-6)


def test_mn_ligand_layout_preserves_attachment_order_and_is_deterministic():
    source = ANCHORS.copy()
    positions = distribute_labels(source, x_radius=2.18, y_radius=1.4)
    anchor_scale = np.maximum(np.max(np.abs(source), axis=0), 1)
    angles = np.mod(np.arctan2(source[:, 1] / anchor_scale[1], source[:, 0] / anchor_scale[0]), 2 * math.pi)
    order = np.argsort(angles, kind="stable")
    target_angles = np.mod(np.arctan2(positions[:, 1] / 1.4, positions[:, 0] / 2.18), 2 * math.pi)
    assert np.all(np.diff(target_angles[order]) >= 0)
    np.testing.assert_array_equal(source, ANCHORS)
    np.testing.assert_array_equal(positions, distribute_labels(source, x_radius=2.18, y_radius=1.4))


def test_crossing_refinement_removes_the_real_data_rank_order_crossings_without_changing_slots():
    count = len(ANCHORS)
    old = np.array([(2.18 * math.cos(math.pi * .9 - 2 * math.pi * index / count),
        1.4 * math.sin(math.pi * .9 - 2 * math.pi * index / count)) for index in range(count)])
    initial = distribute_labels(ANCHORS, x_radius=2.18, y_radius=1.4)
    before = initial.copy()
    refined = reduce_crossings(ANCHORS, initial)
    assert _crossing_count(ANCHORS, old) == 25
    assert _crossing_count(ANCHORS, initial) == 3
    assert _crossing_count(ANCHORS, refined) == 0
    assert {tuple(position) for position in initial} == {tuple(position) for position in refined}
    np.testing.assert_array_equal(initial, before)
    np.testing.assert_array_equal(refined, reduce_crossings(ANCHORS, initial))


def test_shared_atom_anchors_get_distinct_residue_positions_not_merged_scientific_records():
    anchors = np.repeat([[.5, -.5]], 4, axis=0)
    positions = reduce_crossings(anchors, distribute_labels(anchors, x_radius=2.18, y_radius=1.4))
    assert len({tuple(position) for position in positions}) == 4
    assert _crossing_count(anchors, positions) == 0


def test_empty_graphs_do_not_require_a_fake_anchor():
    assert distribute_labels([], x_radius=2, y_radius=1).shape == (0, 2)
    assert reduce_crossings([], []).shape == (0, 2)


@pytest.mark.parametrize("anchors", [[[0, float("nan")]], [[0, float("inf")]], [[0, 1, 2]]])
def test_invalid_anchors_fail_explicitly(anchors):
    with pytest.raises(ValueError):
        distribute_labels(anchors, x_radius=2, y_radius=1)
    with pytest.raises(ValueError):
        reduce_crossings(anchors, [[0, 1]])


@pytest.mark.parametrize("radii", [(0, 1), (-1, 1), (2, float("inf"))])
def test_invalid_radii_fail_explicitly(radii):
    with pytest.raises(ValueError):
        distribute_labels([[0, 1]], x_radius=radii[0], y_radius=radii[1])
