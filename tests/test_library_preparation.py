from __future__ import annotations

import json
from pathlib import Path

import pyarrow.parquet as pq
import pytest
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

from mn_molecular_design_studio.active_learning_settings import ActiveLearningSettings
from mn_molecular_design_studio.collection_policy import (
    automatic_collection_allowed, collection_scope, policy_for_run, set_run_policy,
)
from mn_molecular_design_studio.external_libraries import (
    LibraryCatalog, iter_feature_batches, prepare_fingerprints,
)
from mn_molecular_design_studio.library_centroids import prepare_centroids


@pytest.fixture
def library(tmp_path, monkeypatch):
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path / "app"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(tmp_path / "app"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(tmp_path / "app/workdir/runs"))
    source = tmp_path / "input.csv"
    source.write_text("vendor,structure\n0001,CCO\n0002,OCC\n0003,c1ccccc1\n0004,CC(=O)O\n0005,not_a_smiles\n")
    catalog = LibraryCatalog()
    record = catalog.register([source], name="Example", id_column="vendor", smiles_column="structure")
    return catalog, record


def test_metadata_registration_never_copies_molecules(library):
    catalog, record = library
    assert catalog.get(record["id"]) == record
    assert "fingerprint" not in json.dumps(record)
    assert len(catalog.list()) == 1
    with catalog.database.connection() as connection:
        assert connection.execute("SELECT COUNT(*) FROM external_library_records").fetchone()[0] == 1
        assert connection.execute("SELECT COUNT(*) FROM molecules").fetchone()[0] == 0
    assert catalog.register([record["sources"][0]["path"]], name="Again", id_column="vendor", smiles_column="structure")["id"] == record["id"]


def test_packed_features_preserve_ids_and_record_rejections(library, tmp_path):
    catalog, record = library
    manifest_path = prepare_fingerprints(record, tmp_path / "features", batch_size=2, shard_rows=2)
    manifest = json.loads(manifest_path.read_text())
    assert (manifest["input_rows"], manifest["row_count"], manifest["rejected_rows"]) == (5, 4, 1)
    rows = [row for batch in iter_feature_batches(manifest_path) for row in batch.to_pylist()]
    assert rows[0]["compound_id"] == "0001"
    assert len(rows[0]["fingerprint"]) == 128
    assert rows[0]["source_smiles"] == "CCO"
    assert rows[0]["compound_key"] != rows[1]["compound_key"]
    expected = AllChem.GetMorganFingerprintAsBitVect(
        Chem.MolFromSmiles("CCO"), radius=2, nBits=1024, useFeatures=True, useChirality=True
    )
    assert DataStructs.CreateFromBinaryText(rows[0]["fingerprint"]).ToBitString() == expected.ToBitString()
    assert pq.read_table(manifest_path.parent / "rejected.parquet").to_pylist()[0]["compound_id"] == "0005"
    catalog.attach_preparation(record["id"], manifest_path)
    assert len(catalog.get(record["id"])["preparations"]) == 1
    with catalog.database.connection() as connection:
        assert connection.execute("SELECT COUNT(*) FROM molecules").fetchone()[0] == 0


def test_source_changes_and_duplicate_ids_cannot_publish_completed_manifest(library, tmp_path):
    _, record = library
    source = Path(record["sources"][0]["path"])
    source.write_text("vendor,structure\nA,CCO\nA,CCC\n")
    with pytest.raises(ValueError, match="changed"):
        prepare_fingerprints(record, tmp_path / "changed")
    fresh = LibraryCatalog().register([source], name="Duplicates", id_column="vendor", smiles_column="structure")
    with pytest.raises(ValueError, match="Duplicate"):
        prepare_fingerprints(fresh, tmp_path / "duplicates")
    assert not (tmp_path / "duplicates/manifest.json").exists()
    assert (tmp_path / "duplicates/failure.json").is_file()


def test_centroids_accumulate_all_fingerprint_shards_and_keep_cluster_sizes(library, tmp_path):
    _, record = library
    # Two Parquet shards expose any logic that accidentally counts only the
    # final shard when validating the prepared feature row total.
    features = prepare_fingerprints(record, tmp_path / "features", batch_size=2, shard_rows=2)
    progress = []
    manifest_path = prepare_centroids(
        features, tmp_path / "centroids", similarity_threshold=1.0,
        progress_callback=lambda rows, representatives: progress.append((rows, representatives)),
    )
    manifest = json.loads(manifest_path.read_text())
    reps = pq.read_table(manifest_path.parent / "centroids.parquet").to_pylist()
    assert progress[-1] == (4, 3)
    assert manifest["input_rows"] == 4
    assert manifest["row_count"] == len(reps) == 3
    assert {row["compound_id"] for row in reps} == {"0001", "0003", "0004"}
    assert sum(r["cluster_size"] for r in reps) == 4
    assert not (manifest_path.parent / "assignments.parquet").exists()


def test_exact_clustering_processes_every_input_without_compound_caps(library, tmp_path):
    _, record = library
    features = prepare_fingerprints(record, tmp_path / "features")
    manifest_path = prepare_centroids(features, tmp_path / "all-compounds")
    manifest = json.loads(manifest_path.read_text())
    centroids = pq.read_table(manifest_path.parent / "centroids.parquet").to_pylist()
    assert manifest["input_rows"] == 4
    assert sum(row["cluster_size"] for row in centroids) == 4
    assert manifest["limits"] is None


def test_manifest_content_validation_rejects_changed_shards(library, tmp_path):
    _, record = library
    path = prepare_fingerprints(record, tmp_path / "features")
    shard = path.parent / json.loads(path.read_text())["shards"][0]["path"]
    with shard.open("ab") as handle:
        handle.write(b"changed")
    with pytest.raises(ValueError, match="changed"):
        list(iter_feature_batches(path))


def test_percentile_schedule_is_independent_of_zero_iteration_limit():
    settings = ActiveLearningSettings(max_iterations=0)
    assert settings.percentile(1) == 10
    assert settings.percentile(2) == pytest.approx(4.6415888336)
    assert settings.percentile(10) == settings.percentile(1000) == 0.01
    assert not settings.iteration_limit_reached(1_000_000)
    assert settings.convergence_count(10_000_000) == 1000
    assert ActiveLearningSettings(max_iterations=3).iteration_limit_reached(3)
    with pytest.raises(ValueError):
        ActiveLearningSettings(max_iterations=-1)
    with pytest.raises(ValueError):
        ActiveLearningSettings(binder_percent_target=float("nan"))


def _run(path):
    path.mkdir()
    (path / "metadata.json").write_text("{}")
    return path


def test_collection_policy_survives_refresh_and_parent_restrictions(tmp_path, monkeypatch):
    parent = _run(tmp_path / "parent")
    child = _run(tmp_path / "child")
    assert automatic_collection_allowed(child)
    set_run_policy(parent, "manual")
    set_run_policy(child, "automatic", parent_run_dir=parent)
    assert policy_for_run(child) == "manual"
    set_run_policy(child, "automatic")
    assert not automatic_collection_allowed(child)
    from mn_molecular_design_studio import workspace, molecular_database_ingest as ingest
    monkeypatch.setattr(workspace, "is_studio_job", lambda _: True)
    monkeypatch.setattr(ingest, "MolecularDatabase", lambda: pytest.fail("A restricted run must not open the molecular Database"))
    assert ingest.ingest_run(child, force=True) == 0
    ingest.record_completed_run(child)
    with pytest.raises(ValueError, match="Cycle"):
        set_run_policy(parent, "manual", parent_run_dir=child)


def test_job_creation_persists_collection_scope(tmp_path, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(tmp_path / "runs"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path / "app"))
    with collection_scope("manual"), collection_scope("automatic"):
        job = jobs.create_job("library-preparation", "library_preparation", "library_preparation", {})
    assert policy_for_run(job.run_dir) == "manual"


def test_complete_workflow_keeps_bulk_products_external(library, tmp_path, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.workflows.library_preparation import run_library_preparation
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
    from mn_molecular_design_studio.cluster import distributed_library
    # Exercise the real fingerprint products without needing a live scheduler
    # or starting a second scientific worker process in this contract test.
    monkeypatch.setattr(distributed_library, "_node_specs", lambda workers, **kwargs: [
        {**worker, "name": "test-local", "local": True} for worker in workers])
    def prepare_locally(record, output, workers, **kwargs):
        # The distributed writer accepts the empty directory staged by the
        # workflow; the serial writer creates that directory itself.
        output.rmdir()
        return prepare_fingerprints(record, output), None, None

    monkeypatch.setattr(distributed_library, "prepare_distributed_fingerprints", prepare_locally)
    _, record = library
    path = run_library_preparation([Path(record["sources"][0]["path"])], "From tile", str(tmp_path / "prepared"), id_column="vendor", smiles_column="structure")
    assert read_json(path / "result.json")["success"]
    assert read_json(path / "result.json")["metrics"]["compound_count"] == 4
    result = read_json(path / "result.json")
    assert "centroid_count" not in result["metrics"]
    assert not automatic_collection_allowed(path)
    assert (tmp_path / "prepared/fingerprints/manifest.json").is_file()
    with LibraryCatalog().database.connection() as connection:
        assert connection.execute("SELECT COUNT(*) FROM molecules").fetchone()[0] == 0


def test_standalone_fingerprint_worker_uses_local_scheduler_without_enrollment(tmp_path, monkeypatch):
    from mn_molecular_design_studio.cluster import config, distributed_library
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    config.save({"mode": "disabled", "name": "standalone"})
    monkeypatch.setattr(distributed_library, "_live_local_resources", lambda: (
        True, True, {"cpu_free": 6, "cpu_total": 8}))
    workers = distributed_library._node_specs([{"node_id": config.identity(), "cpu_cores": 4}])
    assert len(workers) == 1 and workers[0]["local"] and workers[0]["cpu_cores"] == 4
    with pytest.raises(ValueError, match="registered"):
        distributed_library._node_specs([{"node_id": "unregistered-remote", "cpu_cores": 4}])
    monkeypatch.setattr(distributed_library, "_live_local_resources", lambda: (False, True, {}))
    with pytest.raises(ValueError, match="scheduler"):
        distributed_library._node_specs([{"node_id": config.identity(), "cpu_cores": 4}])


def test_tile_queues_studio_workflow(library, tmp_path):
    from streamlit.testing.v1 import AppTest
    source = library[1]["sources"][0]["path"]
    app = AppTest.from_string('''
import streamlit as st
from mn_molecular_design_studio.library_preparation_ui import render_setup
def submit(schema, workflow, parameters, project, settings, label):
    st.session_state['submitted'] = dict(workflow=workflow, parameters=parameters, project=project, settings=settings)
render_setup({'schema_version': 1}, user={'role':'admin'}, project_selector=lambda key: 'project-1', submit=submit)
''').run()
    assert not app.exception
    app.text_input(key="library-prep-new-name").set_value("UI library")
    app.text_area(key="library-prep-new-files").set_value(source)
    app.text_input(key="library-prep-new-id").set_value("vendor")
    app.text_input(key="library-prep-new-smiles").set_value("structure")
    app.text_input(key="library-prep-new-output").set_value(str(tmp_path / "ui-output"))
    app.run()
    app.button(key="library-prep-new-submit").click().run()
    assert not app.exception
    submitted = app.session_state["submitted"]
    assert submitted["workflow"] == "library.prepare"
    assert submitted["parameters"]["source_files"] == [source]
    assert submitted["settings"]["_execution"]["placement"] == "local"


def test_library_workflow_bypasses_bulk_peer_staging(monkeypatch, tmp_path):
    from mn_molecular_design_studio.cluster import service, config, bundles
    from mn_molecular_design_studio.cluster.agent import Agent
    monkeypatch.setattr(config, "load", lambda: {"mode": "peer", "default_placement": "auto"})
    monkeypatch.setattr(Agent, "storage_online", lambda _: True)
    monkeypatch.setattr(bundles, "pack", lambda *a, **kw: pytest.fail("Bulk library must not enter a peer input bundle"))
    payload = {"workflow": "library.prepare", "execution": {"placement": "local"}}
    assert service.submit_from_api(payload, {"role": "admin"}, request_dir=tmp_path, token="", authorize_path=lambda _: None) is None
    with pytest.raises(PermissionError):
        service.submit_from_api(payload, {"role": "user"}, request_dir=tmp_path, token="", authorize_path=lambda _: None)
