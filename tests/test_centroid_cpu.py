from pathlib import Path

import numpy as np
import pytest
from rdkit import DataStructs

from mn_molecular_design_studio.centroid_cpu import cluster_packed


def reference(fps, order, threshold, weights):
    vectors = [DataStructs.CreateFromBinaryText(row.tobytes()) for row in fps]
    representatives, populations = [], []
    for source in order:
        values = DataStructs.BulkTanimotoSimilarity(vectors[source], [vectors[i] for i in representatives])
        match = next((i for i, value in enumerate(values) if value >= threshold), None)
        if match is None:
            representatives.append(int(source))
            populations.append(0)
            match = len(representatives) - 1
        populations[match] += int(weights[source])
    return representatives, populations


@pytest.mark.parametrize("threshold", [0.01, 0.3, 0.6, 1.0])
@pytest.mark.parametrize("block_size,cores", [(1, 1), (17, 4), (1024, 4)])
def test_parallel_prefix_index_matches_rdkit_first_match(tmp_path, threshold, block_size, cores):
    rng = np.random.default_rng(42)
    base = rng.random((20, 256)) < 0.2
    bits = base[rng.integers(0, len(base), size=230)].copy()
    bits ^= rng.random(bits.shape) < 0.025
    bits[:2] = False  # Exercise the installed RDKit version's empty-vector rule.
    bits[2:4] = True
    fps = np.ascontiguousarray(np.packbits(bits, axis=1).view(np.uint64))
    order = np.asarray(rng.permutation(len(fps)), dtype=np.uint64)
    weights = rng.integers(1, 10, size=len(fps), dtype=np.uint64)
    expected, counts = reference(fps, order, threshold, weights)
    actual, populations, metrics = cluster_packed(
        fps, order, tmp_path, similarity_threshold=threshold, cpu_cores=cores,
        block_size=block_size, weights=weights,
    )
    assert actual.tolist() == expected
    assert populations.tolist() == counts
    assert sum(populations) == sum(weights)
    assert metrics["approximate"] is False


def test_threshold_boundary_and_cancellation(tmp_path):
    fps = np.array([[0b111], [0b11111], [0b11100000]], dtype=np.uint64)
    order = np.arange(3, dtype=np.uint64)
    reps, _, _ = cluster_packed(fps, order, tmp_path, similarity_threshold=0.6, cpu_cores=2)
    assert reps.tolist() == [0, 2]
    reps, _, _ = cluster_packed(fps, order, tmp_path, similarity_threshold=np.nextafter(0.6, 1), cpu_cores=2)
    assert reps.tolist() == [0, 1, 2]

    def cancel(rows, representatives):
        raise RuntimeError("cancelled by user")

    with pytest.raises(RuntimeError, match="cancelled by user"):
        cluster_packed(fps, order, tmp_path, similarity_threshold=0.6, cpu_cores=2,
                       progress_callback=cancel)
