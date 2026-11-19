"""Deterministic presentation geometry for ligand interaction diagrams."""
from __future__ import annotations

import numpy as np


def distribute_labels(anchors, *, x_radius, y_radius):
    """Preserve attachment order with MN-Ligand's equal-perimeter ellipse layout."""
    coordinates = np.asarray(anchors, dtype=float)
    if coordinates.size == 0:
        return coordinates.reshape((-1, 2))
    if coordinates.ndim != 2 or coordinates.shape[1] != 2 or not np.isfinite(coordinates).all():
        raise ValueError("Interaction label anchors must be finite two-dimensional coordinates.")
    radii = np.asarray([x_radius, y_radius], dtype=float)
    if not np.isfinite(radii).all() or (radii <= 0).any():
        raise ValueError("Interaction label radii must be positive and finite.")
    anchor_scale = np.maximum(np.max(np.abs(coordinates), axis=0), 1.0)
    anchor_angles = np.mod(np.arctan2(coordinates[:, 1] / anchor_scale[1], coordinates[:, 0] / anchor_scale[0]), 2 * np.pi)
    order = np.argsort(anchor_angles, kind="stable")
    count = len(coordinates)
    sample_angles = np.linspace(0, 2 * np.pi, 4097)
    perimeter_points = np.column_stack((x_radius * np.cos(sample_angles), y_radius * np.sin(sample_angles)))
    cumulative = np.concatenate(([0.0], np.cumsum(np.linalg.norm(np.diff(perimeter_points, axis=0), axis=1))))
    fractions = cumulative / cumulative[-1]
    best_angles = None
    best_cost = float("inf")
    ordered_angles = anchor_angles[order]
    for phase in np.linspace(0, 1 / count, 129, endpoint=False):
        slot_fractions = np.mod(phase + np.arange(count) / count, 1)
        slot_angles = np.interp(slot_fractions, fractions, sample_angles)
        differences = np.angle(np.exp(1j * (slot_angles - ordered_angles)))
        cost = float(np.sum(differences ** 2))
        if cost < best_cost:
            best_cost, best_angles = cost, slot_angles
    positions = np.zeros_like(coordinates)
    positions[order, 0] = x_radius * np.cos(best_angles)
    positions[order, 1] = y_radius * np.sin(best_angles)
    return positions


def _crossing_count(anchors, positions):
    first, second = np.triu_indices(len(anchors), 1)
    start_first, start_second = anchors[first], anchors[second]
    end_first, end_second = positions[first], positions[second]

    def orientation(start, end, points):
        delta, offset = end - start, points - start
        return delta[:, 0] * offset[:, 1] - delta[:, 1] * offset[:, 0]

    side_first = orientation(start_first, end_first, start_second) * orientation(start_first, end_first, end_second)
    side_second = orientation(start_second, end_second, start_first) * orientation(start_second, end_second, end_first)
    return int(np.count_nonzero((side_first < -1e-12) & (side_second < -1e-12)))


def reduce_crossings(anchors, positions):
    """Refine ellipse-slot assignments without moving atoms or changing spacing."""
    anchors = np.asarray(anchors, dtype=float)
    positions = np.asarray(positions, dtype=float).copy()
    if not anchors.size:
        return positions.reshape((-1, 2))
    if anchors.shape != positions.shape or anchors.ndim != 2 or anchors.shape[1] != 2 or not np.isfinite(anchors).all() or not np.isfinite(positions).all():
        raise ValueError("Interaction routing requires matching finite two-dimensional coordinates.")

    def cost(candidate):
        return _crossing_count(anchors, candidate), round(float(np.sum((candidate - anchors) ** 2)), 12)

    score = cost(positions)
    for iteration in range(32):
        best, best_score = positions, score
        for first in range(len(positions)):
            for second in range(first + 1, len(positions)):
                candidate = positions.copy()
                candidate[[first, second]] = candidate[[second, first]]
                candidate_score = cost(candidate)
                if candidate_score < best_score:
                    best, best_score = candidate, candidate_score
        if best_score >= score:
            break
        positions, score = best, best_score
    return positions


def place_edge_percentages(figure, axis, labels, *, residue_positions, molecule_bounds):
    """Place labels on the outer edge sections, avoiding rendered label boxes."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.transforms import Bbox

    if not labels:
        return
    canvas = FigureCanvasAgg(figure)
    canvas.draw()
    renderer = canvas.get_renderer()
    edge_labels = {text for _, text in labels}
    fixed = [text.get_window_extent(renderer).padded(2) for text in axis.texts if text not in edge_labels]
    node_radius = 11 * figure.dpi / 72
    for position in residue_positions:
        horizontal, vertical = axis.transData.transform(position)
        fixed.append(Bbox.from_extents(horizontal - node_radius, vertical - node_radius,
            horizontal + node_radius, vertical + node_radius))
    molecule = Bbox(axis.transData.transform(molecule_bounds)).padded(3)
    frame = axis.get_window_extent(renderer)
    occupied = []
    curves = []
    for patch, text in labels:
        path = patch.get_path().transformed(patch.get_transform())
        curve = next(segment for segment, _ in path.iter_bezier() if segment.degree > 0)
        curves.append(curve)
    edge_points = [curve(np.linspace(.15, .95, 17)) for curve in curves]
    for label_index, ((patch, text), curve) in enumerate(zip(labels, curves)):
        bounds = text.get_window_extent(renderer)
        width, height = bounds.width + 4, bounds.height + 3
        best = None
        for fraction in (.78, .69, .87, .60, .93, .52):
            point = curve(fraction)
            tangent = curve(min(1, fraction + .01)) - curve(max(0, fraction - .01))
            length = np.linalg.norm(tangent)
            normal = np.array([-tangent[1], tangent[0]]) / max(length, 1e-9)
            for offset in (0, 5, -5, 10, -10):
                candidate = point + offset * normal
                box = Bbox.from_bounds(candidate[0] - width / 2, candidate[1] - height / 2, width, height)
                cost = 100000 * sum(box.overlaps(other) for other in occupied)
                cost += 10000 * sum(box.overlaps(other) for other in fixed)
                cost += 1000 * box.overlaps(molecule)
                cost += 10000 * (not frame.contains(box.x0, box.y0) or not frame.contains(box.x1, box.y1))
                cost += 30 * sum(any(box.contains(*point) for point in points)
                    for index, points in enumerate(edge_points) if index != label_index)
                cost += 10 * abs(fraction - .78) + abs(offset) / 10
                if best is None or cost < best[0]:
                    best = cost, candidate, box
        text.set_position(axis.transData.inverted().transform(best[1]))
        occupied.append(best[2])
