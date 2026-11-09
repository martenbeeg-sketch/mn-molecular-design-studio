from __future__ import annotations

import copy
import hashlib
import json
import time

import pytest

from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio import scientific_lineage as catalog
from mn_molecular_design_studio.scientific_lineage_ingest import (
    artifact, snapshot_from_candidates, link_retained, add_md_study, chaperone_records,
)

ADMIN = {"id": "owner", "role": "admin"}


@pytest.fixture
def db(tmp_path, monkeypatch):
    for key in ("MN_MOLECULAR_DESIGN_STUDIO_HOME", "MN_MOLECULAR_DESIGN_STUDIO_APP_HOME"):
        monkeypatch.setenv(key, str(tmp_path / "app"))
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(tmp_path / "runs"))
    return MolecularDatabase(tmp_path / "db")


def job(run, tool="rfd4_proteina", **extra):
    return dict(run_id=run, task_group="design", owner_id="owner", project_id="",
                tool=tool, created_at="2026-11-09", updated_at="2026-11-09", **extra)


def candidate(identifier="c1", model=None, passing=True):
    c = {"candidate_id": identifier, "stage": "complex_refolding", "source_tool": "boltz2",
         "binder_sequence": "A", "binder_chains": ["B"], "target_chains": ["A"],
         "metrics": {"native_pass_filters": passing}, "parents": [], "raw_metadata": {"native_candidate_id": identifier}}
    if model is not None:
        c["model_id"] = model
    return c


def save(db, snapshot):
    catalog.put_snapshot(db, snapshot, author=snapshot["job"]["node_id"])


def test_distributed_parent_controller_receipt_and_same_coordinates_do_not_double_count_models(db, tmp_path):
    imported = db.import_records(ADMIN, [{"name": "binder", "components": [{"type": "protein", "chain": "B", "sequence": "A"}]}])
    molecule = imported["molecule_ids"][0]
    path = tmp_path / "model.pdb"
    path.write_text("ATOM      1  CA  ALA B   1       1.000   2.000   3.000  1.00 50.00           C\nEND\n")
    root = snapshot_from_candidates(job("pipeline", "visual_pipeline"), [candidate("collected")], node_id="node-a")
    controller = snapshot_from_candidates(job("controller", "bindcraft2_campaign", campaign_id="pipeline", pipeline_parent_node_id="node-a", campaign_step_index="1"), [candidate()], node_id="node-a")
    worker = snapshot_from_candidates(job("worker", "bindcraft2", cluster_job_id="worker-cluster"), [candidate()], node_id="node-b")
    controller["dependencies"].append({"source_ref": "cluster-job:worker-cluster", "target_ref": controller["job"]["ref"], "role": "execution_child_of", "status": "recorded"})
    root_group = db.create_group(ADMIN, "pipeline results")
    db.add_to_group(ADMIN, root_group, [molecule])
    db.attach_structure(root_group, molecule, {"path": str(path)}, context={"candidate_id": "collected", "task_group": "design", "run_id": "pipeline"})
    link_retained(db, root, local_node="node-a")
    for snapshot, node in [(controller, "node-a"), (worker, "node-b")]:
        group_id = imported["group_id"] if node == "node-a" else db.create_group(ADMIN, "worker")
        if node != "node-a":
            db.add_to_group(ADMIN, group_id, [molecule])
        structure = db.attach_structure(group_id, molecule, {"path": str(path)}, context={"candidate_id": "c1", "task_group": "design", "run_id": snapshot["job"]["run_id"]})
        link_retained(db, snapshot, local_node=node)
        assert snapshot["outputs"][0]["retained"][0]["structure_id"] == structure
    # Worker can arrive first, before its parent or its receipt.
    for snapshot in [worker, root, controller]:
        save(db, snapshot)
    page = catalog.campaigns(db, ADMIN)
    assert page["total"] == 1
    summary = page["items"][0]["summary"]
    assert summary["jobs"] == 3
    assert summary["indexed_occurrences"] == 3
    assert summary["authoritative_output_occurrences"] == 1
    assert summary["acceptance_counts"]["native"]["passed"] == 1
    assert isinstance(summary["counts"], list)
    detail = catalog.campaign_detail(db, ADMIN, page["items"][0]["ref"])
    assert len(detail["outputs"]) == 3
    assert len({r["ref"] for r in detail["outputs"]}) == 3
    assert len(catalog.retained_models(db, ADMIN, page["items"][0]["ref"], stage="native")["items"]) == 1
    assert len(catalog.retained_models(db, ADMIN, page["items"][0]["ref"], stage="native", include_execution_refs=True)["items"]) == 3
    save(db, controller)
    with db.connection() as c:
        assert c.execute("SELECT count(*) FROM molecules").fetchone()[0] == 1
        assert c.execute("SELECT count(*) FROM scientific_outputs").fetchone()[0] == 3


def test_models_variants_and_engines_remain_separate_with_same_target(db):
    target = {"biological_target": "PD-L1", "sequences_by_chain": {"A": "ACDE"}, "boundaries": "1-4"}
    first = snapshot_from_candidates(job("rfd", target_construct=target), [candidate(model="0"), candidate(model="1", passing=False)], node_id="a")
    second = snapshot_from_candidates(job("bc", "bindcraft2", target_construct=target), [candidate()], node_id="b")
    for s in (first, second):
        save(db, s)
    assert catalog.campaigns(db, ADMIN)["total"] == 2
    assert len({o["ref"] for o in first["outputs"]}) == 2
    assert first["outputs"][0]["acceptance"] != first["outputs"][1]["acceptance"]
    assert catalog.campaigns(db, ADMIN, limit=1)["next_offset"] == 1
    assert len(catalog.campaigns(db, ADMIN, offset=1)["items"]) == 1


def test_unknown_model_is_not_assigned_to_known_models_and_zero_is_preserved(db, tmp_path):
    imported = db.import_records(ADMIN, [{"name": "binder", "components": [{"type": "protein", "chain": "B", "sequence": "A"}]}])
    path = tmp_path / "model.pdb"
    path.write_text("ATOM      1  CA  ALA B   1       1.000   2.000   3.000  1.00 50.00           C\nEND\n")
    for model in (0, 1):
        db.attach_structure(imported["group_id"], imported["molecule_ids"][0], {"path": str(path)},
            context={"candidate_id": "c1", "task_group": "design", "run_id": "models", "model_index": model})
    s = snapshot_from_candidates(job("models"), [candidate(), candidate(model=0), candidate(model=1)], node_id="a")
    s["outputs"][0]["source_metadata"]["extra_provenance"] = {"checkpoint": "recorded-v1"}
    link_retained(db, s, local_node="a")
    unknown, zero, one = s["outputs"][:3]
    assert unknown["retained"] == []
    assert zero["model_id"] == "0" and len(zero["retained"]) == 1
    assert one["model_id"] == "1" and len(one["retained"]) == 1
    assert zero["retained"][0]["sha256"] == one["retained"][0]["sha256"]
    assert zero["retained"][0]["structure_id"] != one["retained"][0]["structure_id"]
    save(db, s)
    detail = catalog.campaign_detail(db, ADMIN, catalog.campaigns(db, ADMIN)["items"][0]["ref"])
    assert any(o.get("source_metadata", {}).get("extra_provenance") == {"checkpoint": "recorded-v1"} for o in detail["outputs"])


def test_campaign_executor_does_not_replace_recorded_scientific_generator():
    from mn_molecular_design_studio.complex_analysis_collection_ui import _database_candidate_view
    structure = {"id": "structure1", "context": {"tool": "visual_pipeline", "task_group": "campaign", "candidate_id": "c1"}}
    record = {"id": "mol1", "sources": []}
    observation = {"context": {"candidate": {"candidate_id": "c1", "source_tool": "bindcraft2"}}}
    assert _database_candidate_view(record, structure, [], [], [(observation, {"representation_id": "structure1"})])["source_tool"] == "bindcraft2"
    assert _database_candidate_view(record, structure, [], [], [])["source_tool"] == ""


def test_mixed_construct_job_cannot_appear_as_one_exact_construct(db):
    rows = [candidate("first"), candidate("second")]
    rows[0]["raw_metadata"]["target_sequences_by_chain"] = {"A": "ACDE"}
    rows[1]["raw_metadata"]["target_sequences_by_chain"] = {"A": "ACDEF"}
    s = snapshot_from_candidates(job("mixed-target"), rows, node_id="a")
    save(db, s)
    campaign = catalog.campaigns(db, ADMIN)["items"][0]
    assert campaign["target_key"] == ""
    assert len(campaign["summary"]["targets"][0]["recorded_constructs"]) == 2
    assert s["outputs"][0]["target_construct"] != s["outputs"][1]["target_construct"]


def test_authorization_and_immutable_source_ownership(db, monkeypatch):
    user = {"id": "reader", "role": "user"}
    monkeypatch.setattr("mn_molecular_design_studio.auth.accessible_project_ids", lambda u: [])
    s = snapshot_from_candidates(job("private"), [], node_id="a")
    save(db, s)
    assert catalog.campaigns(db, user)["total"] == 0
    with pytest.raises(PermissionError):
        catalog.campaign_detail(db, user, catalog.campaigns(db, ADMIN)["items"][0]["ref"])
    with pytest.raises(PermissionError):
        catalog.put_snapshot(db, s, author="b")
    changed = copy.deepcopy(s)
    changed["job"]["owner_id"] = "reader"
    with pytest.raises(PermissionError):
        save(db, changed)
    with pytest.raises(PermissionError):
        catalog.validate_snapshot(s, "a", dict(s["job"], project_id="other"))


def test_container_escape_missing_remote_and_owner_local_execution(db, tmp_path):
    root = tmp_path / "job"
    root.mkdir()
    external = tmp_path / "private.pdb"
    external.write_text("private")
    (root / "container.pdb").symlink_to(external)
    j = snapshot_from_candidates(job("job"), [], node_id="a")["job"]
    assert artifact("container.pdb", j, run_dir=root)["reference"] == ""
    assert artifact("missing.cif", j, run_dir=root)["availability"] == "missing"
    assert artifact("/work/private.pdb", j, run_dir=root)["availability"] == "unresolved"
    remote = artifact("cluster://b/design/remote/artifacts/model.cif", j)
    assert remote["availability"] == "owner_local_unverified"
    assert catalog.analysis_placement([remote], local_node="a")["node_id"] == "b"
    assert catalog.analysis_placement([remote], local_node="a")["status"] == "owner_local"
    assert catalog.analysis_placement([remote, dict(remote, node_id="c")], local_node="a")["status"] == "requires_explicit_staging"
    with pytest.raises(ValueError):
        artifact("cluster://b/design/run/../secret", j)
    for invalid in ("cluster:///design/run/model.cif", "runs://unexpected/design/run/model.cif"):
        with pytest.raises(ValueError):
            artifact(invalid, j)
    assert not external.is_relative_to(root)


def test_sequence_only_chaperone_native_history_and_parentage_without_dependency_import(db):
    records = [{"design_id": "child", "parent_design_id": "backbone", "artifact_kind": "sequence", "sequences": '{"B":"ACD"}', "structures": [], "metrics": '{"mpnn":{"score":1.2}}'},
               {"design_id": "prediction", "parent_design_id": "child", "artifact_kind": "structure", "sequences": '{"B":"ACD"}', "structures": [{"name": "rfd", "path": "backbone.cif"}, {"name": "af3", "store_path": "models.catcif:prediction"}], "metrics": {}}]
    snapshot = chaperone_records(job("chaperone", "chaperone"), records, node_id="a")
    assert snapshot["outputs"][0]["artifacts"] == []
    assert snapshot["outputs"][0]["sequences"] == {"B": "ACD"}
    assert [a["step"] for a in snapshot["outputs"][1]["artifacts"]] == ["rfd", "af3"]
    assert snapshot["outputs"][1]["model_id"] == ""
    assert snapshot["dependencies"][0]["evidence"]["native_parent_id"] == "backbone"
    save(db, snapshot)
    with db.connection() as c:
        assert c.execute("SELECT count(*) FROM molecules").fetchone()[0] == 0


def test_md_replicas_continuations_and_exact_trajectory_inputs(db):
    s = snapshot_from_candidates(job("md", "gromacs"), [], node_id="owner-node")
    s["job"]["kind"] = "dynamics"
    study = {"study_id": "legacy-study", "production_runs": [
        {"run_id": "r1-initial", "repeat_index": 1},
        {"run_id": "r1-continued", "repeat_index": 1, "continuation_of_run_id": "r1-initial"},
        {"run_id": "r2", "repeat_index": 2}],
        "mmgbsa_runs": [{"run_id": "endpoint", "source_trajectory_ref": "cluster://owner-node/ligand-md/md/trajectory.dcd",
                          "recorded_trajectory_sha256": "a" * 64, "parameters": {"stride": 10}}],
        "analysis_runs": [{"run_id": "occupancy", "source_trajectory_ref": "cluster://owner-node/ligand-md/md/other.dcd"}]}
    add_md_study(s, study)
    assert len([o for o in s["outputs"] if o["kind"] == "md_replica"]) == 2
    assert len([o for o in s["outputs"] if o["kind"] == "md_segment"]) == 3
    assert len([e for e in s["dependencies"] if e["role"] == "continues"]) == 1
    assert len([e for e in s["dependencies"] if e["role"] == "trajectory_analysis_input"]) == 2
    save(db, s)
    assert catalog.campaigns(db, ADMIN, kind="dynamics")["total"] == 1


def test_docking_chemistry_and_model_provenance_are_distinct_from_compound_identity(db):
    c = candidate("compound1-pose3", model="pose3")
    c["raw_metadata"].update(pose_identity={"id": "pose3", "rank": 3}, ligand_smiles="C[C@H](O)N", atom_mapping={"1": "C1"}, preparation={"protonation": "recorded"}, receptor_preparation={"sha256": "r"})
    s = snapshot_from_candidates(job("dock", "vina", job_type="ligand_docking"), [c], node_id="a")
    assert s["outputs"][0]["kind"] == "docking_pose"
    assert s["outputs"][0]["molecular_context"]["ligand_smiles"] == "C[C@H](O)N"
    save(db, s)
    with db.connection() as conn:
        assert conn.execute("SELECT count(*) FROM molecules").fetchone()[0] == 0


def test_existing_evaluations_keep_id_method_parameters_and_multi_input_dependencies(db):
    imported = db.import_records(ADMIN, [{"name": "binder", "components": [{"type": "protein", "chain": "B", "sequence": "A"}]}])
    molecule = imported["molecule_ids"][0]
    for namespace in ("md.mmgbsa.gb", "md.mmgbsa.pb"):
        db.add_observation(imported["group_id"], molecule, namespace, "analysis", {"energy": -3}, context={
            "analysis_inputs": [{"reference": "cluster://a/ligand-md/md/trajectory.dcd", "role": "trajectory"},
                                {"reference": "cluster://a/ligand-md/md/topology.pdb", "role": "topology"}],
            "settings": {"stride": 10}, "candidate_id": "binder"}, source_key=namespace)
    with db.connection() as c:
        from mn_molecular_design_studio.evaluations import rebuild_missing
        rebuild_missing(c)
        before = [json.loads(r[0])["id"] for r in c.execute("SELECT descriptor FROM molecule_evaluations")]
    s = snapshot_from_candidates(job("analysis", "md_analysis"), [], node_id="a")
    link_retained(db, s, local_node="a")
    assert {o["evaluation_id"] for o in s["outputs"]} == set(before)
    assert {o["descriptor"]["evaluator"] for o in s["outputs"]} == {"md.mmgbsa.gb", "md.mmgbsa.pb"}
    assert len(s["dependencies"]) == 4
    save(db, s)
    with db.connection() as c:
        assert c.execute("SELECT count(*) FROM molecule_evaluations").fetchone()[0] == 2


def test_no_cross_project_parent_merge_and_no_filename_guess(db):
    parent = snapshot_from_candidates(job("parent"), [], node_id="a")
    child = snapshot_from_candidates(dict(job("child"), project_id="private"), [], node_id="a")
    child["dependencies"] = [{"source_ref": child["job"]["ref"], "target_ref": parent["job"]["ref"], "role": "execution_child_of", "status": "recorded"}]
    save(db, parent)
    save(db, child)
    assert catalog.campaigns(db, ADMIN)["total"] == 2


def test_api_pages_and_unauthorized_detail(db, monkeypatch):
    from fastapi import FastAPI
    from fastapi.testclient import TestClient
    from mn_molecular_design_studio import molecular_database_api
    save(db, snapshot_from_candidates(job("private"), [candidate()], node_id="a"))
    monkeypatch.setattr(molecular_database_api, "MolecularDatabase", lambda: db)
    user = dict(ADMIN)
    app = FastAPI()
    molecular_database_api.install(app, lambda: user)
    client = TestClient(app)
    result = client.get("/api/v1/molecules/scientific/campaigns?limit=1").json()
    ref = result["items"][0]["ref"]
    assert client.get(f"/api/v1/molecules/scientific/campaigns/{ref}").json()["output_total"] == 1
    assert client.get(f"/api/v1/molecules/scientific/campaigns/{ref}/models").json()["items"] == []
    user.update(id="other", role="user")
    monkeypatch.setattr("mn_molecular_design_studio.auth.accessible_project_ids", lambda u: [])
    assert client.get(f"/api/v1/molecules/scientific/campaigns/{ref}").status_code == 404


def test_campaign_page_is_index_only_and_requires_explicit_earlier_acceptance(db, monkeypatch):
    from streamlit.testing.v1 import AppTest
    save(db, snapshot_from_candidates(job("native-cohort"), [candidate()], node_id="a"))
    monkeypatch.setattr("mn_molecular_design_studio.logical_runs.catalog", lambda *a, **kw: pytest.fail("No live Results discovery"))
    monkeypatch.setattr("mn_molecular_design_studio.scientific_lineage_ingest.build_snapshot", lambda *a, **kw: pytest.fail("No page-load indexing"))
    app = AppTest.from_string(f'''
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio.scientific_lineage_ui import render_campaign_browser
render_campaign_browser(MolecularDatabase({str(db.root)!r}), {ADMIN!r})
''', default_timeout=10).run()
    assert not app.exception
    ref = catalog.campaigns(db, ADMIN)["items"][0]["ref"]
    start = time.monotonic()
    app.selectbox(key="scientific-campaign").select(ref).run()
    assert time.monotonic() - start < 5
    assert not app.exception
    assert app.selectbox(key="scientific-stage-" + ref).value == "final"
    assert any("No final decision" in info.value for info in app.info)
    app.selectbox(key="scientific-stage-" + ref).select("native").run()
    assert not app.exception


def test_existing_job_replication_carries_lineage_and_tombstones_without_files(db, tmp_path, monkeypatch):
    from mn_molecular_design_studio.cluster import config, replication
    from mn_molecular_design_studio.cluster.store import Store
    a, b = "11111111-1111-4111-8111-111111111111", "22222222-2222-4222-8222-222222222222"
    peer_db = MolecularDatabase(tmp_path / "peer-db")
    active = {"node": a, "database": db}
    monkeypatch.setattr(config, "identity", lambda: active["node"])
    monkeypatch.setattr("mn_molecular_design_studio.molecular_database.MolecularDatabase", lambda: active["database"])
    monkeypatch.setattr(replication, "capture_identity", lambda: None)
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "peer-a"))
    snapshot = snapshot_from_candidates(job("source"), [candidate()], node_id=a)
    save(db, snapshot)
    cluster_id = config.job_id(a, "design/source")
    envelope = dict(snapshot["job"], cluster_job_id=cluster_id, local_reference="design/source")
    with Store().connect() as c:
        c.execute("INSERT INTO jobs (id,node_id,local_ref,owner_id,project_id,payload,created_at) VALUES (?,?,?,?,?,?,?)",
                  (cluster_id, a, "design/source", "owner", "", json.dumps(envelope), "2026-11-09"))
    sender = replication.Replica()
    sender.capture_jobs()
    with sender.db() as c:
        row = c.execute("SELECT key,kind,author,counter,value FROM records WHERE kind='job'").fetchone()
        record = dict(row, value=json.loads(row["value"]))
    assert record["value"]["scientific_lineage"]["job"]["cluster_job_id"] == cluster_id
    assert "content" not in json.dumps(record["value"]["scientific_lineage"])
    active.update(node=b, database=peer_db)
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "peer-b"))
    with Store().connect() as c:
        c.execute("INSERT INTO nodes(id,name,token_hash) VALUES (?,?,?)", (a, "source", "fixture"))
    receiver = replication.Replica()
    receiver.receive([record])
    receiver.receive([record])
    assert catalog.campaigns(peer_db, ADMIN)["total"] == 1
    reopened = MolecularDatabase(peer_db.root)
    ref = catalog.campaigns(reopened, ADMIN)["items"][0]["ref"]
    assert catalog.campaign_detail(reopened, ADMIN, ref)["outputs"][0]["candidate_id"] == "c1"
    # Existing transport accepts old job payloads without a lineage extension.
    old = copy.deepcopy(record)
    old["value"].pop("scientific_lineage")
    old["counter"] += 1
    receiver.receive([old])
    assert catalog.campaigns(peer_db, ADMIN)["total"] == 1
    tombstone = dict(record, value=None, counter=record["counter"] + 2)
    # Existing receive() intentionally refuses demoting an accepted job to a
    # placeholder/tombstone. The explicit owner-authorized removal path must
    # remove its derived catalog too; do not change that transport policy.
    receiver._apply("job", record["key"], None, a)
    assert catalog.campaigns(peer_db, ADMIN)["total"] == 0
    assert catalog.campaigns(db, ADMIN)["total"] == 1


def test_constructs_with_same_target_name_are_not_pooled(db):
    for run, seq in [("one", "ACD"), ("two", "ACDE")]:
        save(db, snapshot_from_candidates(job(run, target_construct={"biological_target": "Tau", "sequence": seq}), [candidate()], node_id="a"))
    rows = catalog.campaigns(db, ADMIN)["items"]
    assert len({r["target_key"] for r in rows}) == 2
    assert catalog.campaigns(db, ADMIN, target=rows[0]["target_key"])["total"] == 1


def test_snapshot_bound_and_invalid_source_cannot_be_published(db):
    s = snapshot_from_candidates(job("oversized"), [candidate()], node_id="a")
    s["oversized"] = "x" * catalog.MAX_BYTES
    with pytest.raises(ValueError, match="oversized"):
        save(db, s)
    with pytest.raises(ValueError):
        catalog.job_reference("a?unsafe", "design", "run")
    assert catalog.campaigns(db, ADMIN)["total"] == 0


def test_completion_indexes_results_only_failed_candidates_without_promotion(db, tmp_path, monkeypatch):
    from mn_molecular_design_studio import molecular_database_ingest as ingest
    run = tmp_path / "runs/design/failed"
    (run / "artifacts/normalized_candidates").mkdir(parents=True)
    (run / "metadata.json").write_text(json.dumps(job("failed", status="failed")))
    (run / "artifacts/normalized_candidates/candidates.jsonl").write_text(json.dumps(candidate(passing=False)) + "\n")
    monkeypatch.setattr(ingest, "_record_completed_run", lambda r: None)
    monkeypatch.setattr("mn_molecular_design_studio.workspace.is_studio_job", lambda r: True)
    monkeypatch.setattr("mn_molecular_design_studio.molecular_database.MolecularDatabase", lambda: db)
    monkeypatch.setattr("mn_molecular_design_studio.cluster.config.identity", lambda: "a")
    ingest.record_completed_run(run)
    assert catalog.campaigns(db, ADMIN)["total"] == 1
    with db.connection() as c:
        assert c.execute("SELECT count(*) FROM molecules").fetchone()[0] == 0


def test_controlled_backfill_is_staging_only_and_idempotent(db, tmp_path, monkeypatch, capsys):
    from mn_molecular_design_studio.scientific_lineage_ingest import main
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps([{"node_id": "a", "metadata": job("explicit"), "candidates": [candidate()]}]))
    destination = tmp_path / "backfill"
    monkeypatch.setattr("sys.argv", ["backfill", "--staging-database-root", str(destination), "--manifest", str(manifest)])
    main()
    main()
    assert catalog.campaigns(MolecularDatabase(destination), ADMIN)["total"] == 1
    with MolecularDatabase(destination).connection() as c:
        assert c.execute("SELECT count(*) FROM molecules").fetchone()[0] == 0
    monkeypatch.setattr("sys.argv", ["backfill", "--staging-database-root", str(tmp_path / "app/database"), "--manifest", str(manifest)])
    with pytest.raises(SystemExit):
        main()
