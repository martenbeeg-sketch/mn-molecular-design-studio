"""Bulk jobs never index the loop; only a bounded final shortlist is promoted."""
from pathlib import Path
import json

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from mn_molecular_design_studio import collection_policy as policy
from mn_molecular_design_studio.bulk_screening_results import publish_final_selection, import_final_results
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio.engine.workflow_engine.core import jobs, workflow_queue
from mn_molecular_design_studio.engine.workflow_engine.workflows import ligand_tools, active_learning_controller as controller

USER = {"id": "a" * 32, "role": "admin"}


@pytest.fixture
def bulk(tmp_path, monkeypatch):
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(tmp_path / "runs"))
    run = tmp_path / "runs/campaign/final"
    data = run / "artifacts/active_learning"
    data.mkdir(parents=True)
    (run / "molecular_design_studio.json").write_text(json.dumps({"studio_owner_id": USER["id"]}))
    (run / "metadata.json").write_text(json.dumps({"status": "completed"}))
    (run / "result.json").write_text(json.dumps({"success": True}))
    rows = [dict(library_id="lib", compound_key=i, compound_id=f"original-{i}", compound_name=f"Name {i}",
                 smiles=s, score=-10 + i, score_field="test_score", score_units="REU",
                 first_round=1, best_round=2) for i, s in enumerate(["CCO", "CCN", "CCC", "CCCC", "COC"])]
    (data / "retained_poses").mkdir()
    pose = data / "retained_poses/top.pdbqt"
    pose.write_text("ROOT\nPOSE\n")
    rows[0]["pose_path"] = str(pose)
    pq.write_table(pa.Table.from_pylist(rows), data / "docking_scores.parquet", row_group_size=2)
    (data / "round_history.json").write_text(json.dumps({"score_field": "test_score", "score_units": "REU"}))
    return run, MolecularDatabase(tmp_path / "database")


def test_only_final_selected_rows_import_and_external_table_is_unchanged(bulk):
    run, db = bulk
    path = run / "artifacts/active_learning/docking_scores.parquet"
    original = path.read_bytes()
    result = publish_final_selection(run, USER, 2, database=db)
    assert result["imported_count"] == 2
    records = db.list_molecules(USER, group_id=result["group_id"])["items"]
    assert {r["metadata"]["compound_id"] for r in records} == {"original-0", "original-1"}
    assert {r["name"] for r in records} == {"Name 0", "Name 1"}
    assert len(records[0]["observations"]) == 1
    with db.connection() as c:
        assert c.execute("SELECT count(*) FROM molecules").fetchone()[0] == 2
        obs = [json.loads(r[0]) for r in c.execute("SELECT context FROM molecule_observations")]
        assert all(o["units"] == "REU" and o["saved_selection"] for o in obs)
        structure = c.execute("SELECT path,context FROM molecule_structures").fetchone()
        assert Path(structure[0]).read_text() == "ROOT\nPOSE\n"
        assert json.loads(structure[1])["portable"]
    assert path.read_bytes() == original
    assert pq.ParquetFile(path).metadata.num_rows == 5
    assert publish_final_selection(run, USER, 2, database=db)["group_id"] == result["group_id"]
    with db.connection() as c:
        assert c.execute("SELECT count(*) FROM molecule_observations").fetchone()[0] == 2


def test_import_can_be_disabled_completely_and_limit_is_hard(bulk):
    run, db = bulk
    assert publish_final_selection(run, USER, 0, database=db)["imported_count"] == 0
    for bad in (1001, -1, True, 1.5):
        with pytest.raises(ValueError, match="1,000"):
            publish_final_selection(run, USER, bad, database=db)
    with db.connection() as c:
        assert c.execute("SELECT count(*) FROM molecules").fetchone()[0] == 0
    manifest = json.loads((run / "artifacts/active_learning/external_results.json").read_text())
    assert manifest["compound_count"] == 5 and manifest["database_selection"]["group_id"] is None


def test_thousand_row_cap_does_not_read_or_import_the_rest(bulk):
    run, db = bulk
    rows = [dict(library_id="lib", compound_key=i, compound_id=f"ID-{i}",
                 smiles="C" + "".join("N" if (i >> bit) & 1 else "C" for bit in range(10)),
                 score=float(i), pose_path="") for i in range(1002)]
    path = run / "artifacts/active_learning/docking_scores.parquet"
    pq.write_table(pa.Table.from_pylist(rows), path, row_group_size=127)
    result = publish_final_selection(run, USER, 1000, database=db)
    # Some structures are symmetric; membership count may deduplicate them.
    assert 0 < result["imported_count"] <= 1000
    with db.connection() as c:
        assert c.execute("SELECT count(*) FROM molecule_sources").fetchone()[0] == 1000
        assert c.execute("SELECT count(*) FROM molecule_observations").fetchone()[0] == 1000
        assert not c.execute("SELECT 1 FROM molecule_sources WHERE json_extract(provenance,'$.compound_id') IN ('ID-1000','ID-1001')").fetchone()
    assert pq.ParquetFile(path).metadata.num_rows == 1002


def test_manual_promotion_rejects_a_running_loop(bulk, monkeypatch):
    run, db = bulk
    from mn_molecular_design_studio import molecular_database_results
    monkeypatch.setattr(molecular_database_results, "_authorized_run", lambda *a: (run,))
    (run / "metadata.json").write_text(json.dumps({"status": "running"}))
    with pytest.raises(ValueError, match="completed"):
        import_final_results(run, USER, 2, database=db)


@pytest.mark.parametrize("engine", ["vina", "gnina", "unidock_pro", "autodock_gpu", "rosetta_ligand"])
@pytest.mark.parametrize("enabled", [False, True])
def test_docking_policy_persisted_before_worker_dispatch(tmp_path, monkeypatch, engine, enabled):
    monkeypatch.setattr(jobs, "runs_root", lambda: tmp_path / "runs")
    def dispatched(run):
        assert policy.policy_for_run(run) == ("automatic" if enabled else "manual")
    monkeypatch.setattr(workflow_queue, "spawn_worker_for_run", dispatched)
    kwargs = dict(target_pdb=tmp_path / "target.pdb", compound_set_file=tmp_path / "compounds.json",
                  center_x=0., center_y=0., center_z=0., import_to_database=enabled)
    run = workflow_queue.enqueue_workflow_call(getattr(ligand_tools, "run_ligand_" + engine),
        task_group="ligand-docking", tool=engine, job_type="ligand_docking", kwargs=kwargs,
        cpu_cores=1, gpu_device="none")
    assert policy.automatic_collection_allowed(run) is enabled


def test_direct_docking_respects_import_switch(tmp_path, monkeypatch):
    monkeypatch.setattr(jobs, "runs_root", lambda: tmp_path)
    def fake(**kwargs):
        return jobs.create_job("ligand-docking", "ligand_docking", "vina", {}, {}).run_dir
    monkeypatch.setattr(ligand_tools, "_run_docking", fake)
    run = ligand_tools.run_ligand_vina(tmp_path / "p.pdb", tmp_path / "c.json", 0., 0., 0., import_to_database=False)
    assert policy.policy_for_run(run) == "manual"
    from mn_molecular_design_studio.molecular_database_ingest import ingest_run
    assert ingest_run(run, force=True) == 0


def test_local_active_learning_batch_forces_external_only(tmp_path, monkeypatch):
    monkeypatch.setattr(jobs, "runs_root", lambda: tmp_path / "runs")
    monkeypatch.setattr(workflow_queue, "spawn_worker_for_run", lambda run: None)
    from mn_molecular_design_studio import workspace
    monkeypatch.setattr(workspace, "mark_studio_job", lambda *a, **kw: None)
    monkeypatch.setattr(controller, "_wait_local_job", lambda run, **kw: {"status": "completed"})
    result = controller._dispatch_docking_shard(campaign=tmp_path / "campaign",
        worker={"local": True, "node_id": "local", "node_name": "Organizer", "cpu_cores": 1, "gpu_device": "0"},
        request={"parameters": dict(target_pdb=tmp_path / "p.pdb", compound_set_file=tmp_path / "c.json",
                 center_x=0., center_y=0., center_z=0., cpu_cores=1, gpu_device="0", import_to_database=True),
                 "owner_id": USER["id"], "project_id": "", "round_number": 1}, ssh_user="compute", cancel_check=lambda: False)
    assert policy.policy_for_run(result["worker_run_dir"]) == "manual"


def test_remote_active_learning_launcher_forces_external_only_before_dispatch(tmp_path, monkeypatch):
    from mn_molecular_design_studio.cluster import active_learning_worker
    from mn_molecular_design_studio import workspace
    monkeypatch.setattr(workspace, "ensure_studio_workspace", lambda: None)
    monkeypatch.setattr(workspace, "mark_studio_job", lambda *a, **kw: None)
    monkeypatch.setattr(jobs, "runs_root", lambda: tmp_path / "runs")
    def dispatch(run):
        assert policy.policy_for_run(run) == "manual"
        metadata = jobs.read_json(run / "metadata.json")
        jobs.write_json(run / "metadata.json", {**metadata, "status": "completed"})
        jobs.write_json(run / "result.json", {"success": True})
    monkeypatch.setattr(workflow_queue, "spawn_worker_for_run", dispatch)
    request = tmp_path / "request.json"
    request.write_text(json.dumps({"schema_version": "studio.active-learning-docking-shard.v1",
        "parameters": {"target_pdb": str(tmp_path / "p.pdb"), "compound_set_file": str(tmp_path / "c.json"),
        "center_x": 0., "center_y": 0., "center_z": 0., "cpu_cores": 1, "gpu_device": "0", "import_to_database": True},
        "owner_id": USER["id"], "project_id": "", "parent_run_dir": "/organizer/only/campaign",
        "round_number": 1, "node_id": "remote", "node_name": "Remote"}))
    result = active_learning_worker.run(request)
    assert result["status"] == "completed"
    assert policy.policy_for_run(result["run_dir"]) == "manual"
