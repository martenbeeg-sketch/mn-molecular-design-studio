import json
from pathlib import Path
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from mn_molecular_design_studio import centroid_chunks as chunks
from mn_molecular_design_studio.centroid_cpu import cluster_chunk_file
from mn_molecular_design_studio.external_libraries import (
    FEATURE_SCHEMA, FingerprintSpec, _fingerprint_batch, atomic_json, digest_file,
    feature_arrow_schema,
)


@pytest.fixture
def library(tmp_path):
    spec = FingerprintSpec()
    smiles = ["CCO", "CCN", "CCCC", "c1ccccc1", "c1ccncc1", "CC(=O)O",
              "C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O", "CCOC", "CCS"] * 5
    rows = _fingerprint_batch([(i + 1, str(i), None, smi) for i, smi in enumerate(smiles)], spec)
    table = pa.Table.from_pylist([
        {"compound_key": row["row_number"], "compound_id": row["compound_id"],
         "compound_name": row["compound_name"], "source_smiles": row["source_smiles"],
         "fingerprint": row["fingerprint"]} for row in rows
    ], schema=feature_arrow_schema(spec))
    path = tmp_path / "features.parquet"
    pq.write_table(table, path)
    manifest = tmp_path / "manifest.json"
    atomic_json(manifest, {"schema_version": FEATURE_SCHEMA, "status": "completed",
                          "library_id": "test-library", "fingerprint_spec": spec.metadata(),
                          "row_count": len(rows), "shards": [{"path": path.name,
                          "rows": len(rows), "sha256": digest_file(path)}]})
    return manifest, path


def test_balanced_sizes_and_cluster_count_proportional_quotas():
    sizes = chunks.balanced_chunk_sizes(100_000_000)
    assert len(sizes) == 334
    assert sum(sizes) == 100_000_000
    assert max(sizes) - min(sizes) == 1
    assert chunks.balanced_chunk_sizes(10, chunk_count=3) == [4, 3, 3]
    quotas = chunks.representative_quotas([1, 100, 100], 12)
    assert quotas[0] == 1 and sum(quotas) == 12
    assert abs(quotas[1] - quotas[2]) <= 1
    assert chunks.representative_quotas([100, 900], 100) == [10, 90]
    assert chunks.representative_quotas([1, 1000], 2) == [1, 1]
    assert chunks.representative_quotas([1, 1000], 1) == [0, 1]
    assert chunks.representative_quotas([1, 2], 100) == [1, 2]
    assert sum(chunks.representative_quotas([1000] * 324, 300_000)) == 300_000
    with pytest.raises(ValueError):
        chunks.balanced_chunk_sizes(5, chunk_count=6)


def test_shuffle_is_balanced_complete_and_reproducible(library, tmp_path):
    _, source = library
    outputs = []
    for name, seed in [("first", 42), ("repeat", 42), ("different", 43)]:
        files = chunks._shuffle_chunks([source], [17, 17, 16], tmp_path / name,
                                       seed=seed, cpu_cores=2, fingerprint_bytes=128,
                                       progress=lambda *args: None)
        outputs.append([pq.read_table(path)["compound_key"].to_pylist() for path in files])
        assert [len(values) for values in outputs[-1]] == [17, 17, 16]
        assert sorted(sum(outputs[-1], [])) == list(range(1, 51))
        assert not (tmp_path / name / "sort-spill").exists()
    assert outputs[0] == outputs[1]
    assert outputs[0] != outputs[2]


def test_cpu_fingerprints_keep_openvs_features_and_chirality(library):
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem

    _, source = library
    table = pq.read_table(source)
    for smiles, actual in zip(table["source_smiles"].to_pylist(), table["fingerprint"].to_pylist()):
        expected = AllChem.GetMorganFingerprintAsBitVect(
            Chem.MolFromSmiles(smiles), radius=2, nBits=1024,
            useFeatures=True, useChirality=True,
        )
        assert actual == DataStructs.BitVectToBinaryText(expected)
    assert table["fingerprint"][6].as_py() != table["fingerprint"][7].as_py()


@pytest.mark.parametrize("method", ["sampled", "hierarchical"])
def test_chunk_workflow_retains_correct_counts_and_removes_temporaries(library, tmp_path, monkeypatch, method):
    from mn_molecular_design_studio.cluster import distributed_library

    manifest, _ = library
    specs = [{"node_id": "local", "name": "local", "local": True, "cpu_cores": 2}]
    monkeypatch.setattr(distributed_library, "_node_specs", lambda nodes: specs)

    def local_round(inputs, specs, directory, threshold, ssh_user, progress):
        directory.mkdir(parents=True)
        outputs, metrics = [], []
        for index, path in enumerate(inputs):
            output = directory / f"chunk-{index:08d}.parquet"
            metrics.append(cluster_chunk_file(path, output, threshold, 2))
            outputs.append(output)
        summaries = [{"node_id": "local", "name": "local", "cpu_worker_count": 2,
                      "input_rows": sum(item["input_rows"] for item in metrics),
                      "row_count": sum(item["row_count"] for item in metrics),
                      "chunk_count": len(inputs), "status": "completed"}]
        return outputs, summaries, metrics

    monkeypatch.setattr(chunks, "_run_round", local_round)
    output = tmp_path / (method + "-output")
    result_path = chunks.prepare_chunked_centroids(
        manifest, output, worker_nodes=[{}], chunk_size=17,
        final_candidate_count=8, method=method, max_merge_rounds=0,
    )
    result = json.loads(result_path.read_text())
    table = pq.read_table(output / "centroids.parquet")
    assert result["row_count"] == table.num_rows
    assert result["fingerprint_spec"]["use_chirality"] is True
    assert result["fingerprint_spec"]["use_features"] is True
    assert set(path.name for path in output.iterdir()) == {"manifest.json", "centroids.parquet"}
    assert result["input_rows"] == 50
    assert sum(table["cluster_size"].to_pylist()) == (8 if method == "sampled" else 50)
    assert result["coverage_guarantee"] is False
    if method == "sampled":
        assert result["candidate_count"] == 8
        assert sum(result["chunk_representative_quotas"]) == 8
        assert result["cluster_size_scope"] == "sampled_chunk_representatives"
    else:
        assert result["stop_reason"] in {"one_final_chunk", "no_reduction_in_round"}


def test_failure_never_publishes_completed_manifest(library, tmp_path, monkeypatch):
    from mn_molecular_design_studio.cluster import distributed_library

    manifest, _ = library
    monkeypatch.setattr(distributed_library, "_node_specs", lambda nodes: [
        {"node_id": "local", "name": "local", "local": True, "cpu_cores": 1}])

    def fail(*args, **kwargs):
        raise RuntimeError("worker cancelled")

    monkeypatch.setattr(chunks, "_run_round", fail)
    output = tmp_path / "failure"
    with pytest.raises(RuntimeError, match="worker cancelled"):
        chunks.prepare_chunked_centroids(manifest, output, worker_nodes=[{}], chunk_size=17)
    assert not (output / "manifest.json").exists()
    assert not (output / ".centroid-work").exists()
    assert (output / "failure.json").exists()


def test_internal_worker_processes_multiple_chunks(library, tmp_path):
    from mn_molecular_design_studio.cluster.distributed_library import _Worker
    import mn_molecular_design_studio

    root = Path(mn_molecular_design_studio.__file__).parent
    _, source = library
    inputs = chunks._shuffle_chunks([source], [17, 17, 16], tmp_path / "inputs",
                                    seed=42, cpu_cores=2, fingerprint_bytes=128,
                                    progress=lambda *args: None)
    worker = _Worker({"node_id": "unit-test", "name": "unit-test", "local": True, "cpu_cores": 2},
                     worker_script=root / "cluster/library_centroid_worker.py",
                     library_module=root / "external_libraries.py", temp_root=tmp_path / "workers",
                     ssh_user="compute", worker_args=["--tasks", "3", "--threshold", "0.6"])
    try:
        for index, path in enumerate(inputs):
            worker.send([{"chunk": index, "input": str(path)}])
        worker.close_input()
        worker.wait()
        manifest = json.loads((worker.fetch() / "worker-manifest.json").read_text())
        assert manifest["parallel_chunks"] == 2
        assert manifest["cores_per_chunk"] == 1
        assert sum(item["input_rows"] for item in manifest["chunks"]) == 50
        assert sum(item["population"] for item in manifest["chunks"]) == 50
        assert len({item["path"] for item in manifest["chunks"]}) == 3
    finally:
        worker.abort()
