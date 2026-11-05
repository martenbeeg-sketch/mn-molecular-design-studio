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


def _seed_historical_automatic_selection(run, database, limit):
    """Build the pre-P0 automatic source rows without invoking the new publisher."""
    import socket
    from mn_molecular_design_studio.molecular_database import _id, _json

    run = Path(run).resolve()
    rows = pq.read_table(run / "artifacts/active_learning/docking_scores.parquet").to_pylist()[:limit]
    origin = {"type": "job", "collection_kind": "curated", "saved_selection": True,
              "final_screening_selection": True, "task_group": "campaign", "run_id": run.name,
              "run_dir": str(run), "external_scores": str(run / "artifacts/active_learning/docking_scores.parquet"),
              "source_node": socket.gethostname(), "explicit_selection": False}
    group_id = _id("screening-final-", str(run))
    database.create_group(USER, "Historical automatic shortlist", import_method="tool_output",
                          origin=origin, group_id=group_id)
    with database.connection() as connection:
        for rank, row in enumerate(rows, 1):
            source_key = (f"{socket.gethostname()}/campaign/{run.name}:final:"
                          f"{row.get('library_id')}:{row.get('compound_key', rank)}")
            record = {"name": row["compound_name"],
                      "components": [{"type": "small_molecule", "chain": "L", "smiles": row["smiles"]}],
                      "metadata": {"compound_id": row["compound_id"]}}
            molecule_id = database._put(
                connection, group_id, record, source_key=source_key,
                origin={**origin, "compound_id": row["compound_id"], "final_rank": rank},
                preserve_existing=True,
            )
            observation_context = {**origin, "compound_id": row["compound_id"], "units": "REU"}
            connection.execute(
                "INSERT INTO molecule_observations VALUES (?,?,?,?,?,?,?,?)",
                (_id("score-", source_key), molecule_id, group_id, "docking", run.name,
                 _json({"test_score": row["score"], "final_rank": rank}),
                 _json(observation_context), "historical-fixture"),
            )
    return group_id


def _replicate_group_entry(source, target, group_id, molecule_id, author):
    from mn_molecular_design_studio.molecular_database_peer import apply_record

    with source.connection() as connection:
        group = dict(connection.execute("SELECT * FROM molecule_groups WHERE id=?", (group_id,)).fetchone())
        molecule = dict(connection.execute("SELECT * FROM molecules WHERE id=?", (molecule_id,)).fetchone())
        member = dict(connection.execute(
            "SELECT * FROM molecule_members WHERE group_id=? AND molecule_id=?", (group_id, molecule_id),
        ).fetchone())
        observations = [dict(row) for row in connection.execute(
            "SELECT * FROM molecule_observations WHERE group_id=?", (group_id,),
        )]
        sources = [dict(row) for row in connection.execute(
            "SELECT * FROM molecule_sources WHERE group_id=?", (group_id,),
        )]
    apply_record("molecule_group", f"molecule_group:{author}:{group_id}",
                 {"node_id": author, "group": group}, author, database=target)
    apply_record("molecule_entry", f"molecule_entry:{author}:{group_id}:{molecule_id}",
                 {"node_id": author, "group_id": group_id, "molecule": molecule,
                  "member": member, "structures": [], "observations": observations,
                  "sources": sources}, author, database=target)


def test_only_final_selected_rows_import_and_external_table_is_unchanged(bulk):
    run, db = bulk
    path = run / "artifacts/active_learning/docking_scores.parquet"
    original = path.read_bytes()
    result = publish_final_selection(run, USER, 2, database=db, explicit_selection=True)
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
    assert publish_final_selection(run, USER, 2, database=db, explicit_selection=True)["group_id"] == result["group_id"]
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


def test_explicit_compound_ids_import_only_the_selection_and_replay_idempotently(bulk, monkeypatch):
    run, db = bulk
    from mn_molecular_design_studio import molecular_database_results
    monkeypatch.setattr(molecular_database_results, "_authorized_run", lambda selected, _user: (selected,))
    path = run / "artifacts/active_learning/docking_scores.parquet"
    original_scores = path.read_bytes()
    imported = db.import_records(
        USER,
        [{"name": "HY-L126 existing identity", "components": [
            {"type": "small_molecule", "chain": "L", "smiles": "CCN"},
        ]}],
        name="Intentional reusable import", import_method="user_uploads",
    )
    original_id = imported["molecule_ids"][0]

    automatic_group_id = _seed_historical_automatic_selection(run, db, 2)
    with db.connection() as connection:
        historical_observation = connection.execute(
            "SELECT id,molecule_id,group_id,metrics,context FROM molecule_observations "
            "WHERE group_id=? AND json_extract(context,'$.compound_id')='original-1'",
            (automatic_group_id,),
        ).fetchone()
        historical_observation = tuple(historical_observation)

    first = import_final_results(
        run, USER, 0, database=db, selected_compound_ids=["original-1", "original-3"],
    )
    assert first["explicit_selection"] is True
    assert first["selected_compound_ids"] == ["original-1", "original-3"]
    assert first["group_id"].startswith("screening-selected-")
    selected = db.list_molecules(USER, group_id=first["group_id"])["items"]
    assert {row["metadata"]["compound_id"] for row in selected} == {"original-1", "original-3"}
    reused = next(row for row in selected if row["metadata"]["compound_id"] == "original-1")
    assert reused["id"] == original_id
    with db.connection() as connection:
        assert connection.execute("SELECT count(*) FROM molecules").fetchone()[0] == 3
        assert connection.execute("SELECT count(*) FROM molecule_members WHERE group_id=?",
                                  (first["group_id"],)).fetchone()[0] == 2
        assert connection.execute("SELECT count(*) FROM molecule_observations").fetchone()[0] == 3
        assert tuple(connection.execute(
            "SELECT id,molecule_id,group_id,metrics,context FROM molecule_observations WHERE id=?",
            (historical_observation[0],),
        ).fetchone()) == historical_observation

    second = import_final_results(
        run, USER, 0, database=db, selected_compound_ids=["original-1", "original-3"],
    )
    assert second["group_id"] == first["group_id"]
    with db.connection() as connection:
        assert connection.execute("SELECT count(*) FROM molecules").fetchone()[0] == 3
        assert connection.execute("SELECT count(*) FROM molecule_members WHERE group_id=?",
                                  (first["group_id"],)).fetchone()[0] == 2
        assert connection.execute("SELECT count(*) FROM molecule_observations").fetchone()[0] == 3
    assert path.read_bytes() == original_scores
    assert pq.ParquetFile(path).metadata.num_rows == 5


def test_exact_compound_import_rejects_missing_or_ambiguous_ids_without_writes(bulk, monkeypatch):
    run, db = bulk
    from mn_molecular_design_studio import molecular_database_results
    monkeypatch.setattr(molecular_database_results, "_authorized_run", lambda selected, _user: (selected,))
    with pytest.raises(ValueError, match="not found"):
        import_final_results(run, USER, 0, database=db, selected_compound_ids=["not-in-results"])
    with pytest.raises(ValueError, match="either an exact"):
        publish_final_selection(run, USER, 1, database=db, selected_compound_ids=["original-1"])
    with db.connection() as connection:
        assert connection.execute("SELECT count(*) FROM molecules").fetchone()[0] == 0


def test_foreign_private_evaluation_cannot_suppress_local_selected_evaluation(bulk, tmp_path, monkeypatch):
    run, database = bulk
    from mn_molecular_design_studio import molecular_database_results
    from mn_molecular_design_studio.molecular_database import _id
    monkeypatch.setattr(molecular_database_results, "_authorized_run", lambda selected, _user: (selected,))
    source = MolecularDatabase(tmp_path / "foreign-source")
    foreign_user = {"id": "another-owner", "role": "admin"}
    imported = source.import_records(foreign_user, [{
        "name": "Foreign experiment identity",
        "components": [{"type": "small_molecule", "chain": "L", "smiles": "CCN"}],
    }], name="Private foreign experiment", origin={"type": "job", "collection_kind": "curated"})
    foreign_molecule_id, foreign_group_id = imported["molecule_ids"][0], imported["group_id"]
    source.add_observation(
        foreign_group_id, foreign_molecule_id, "docking", run.name, {"test_score": -99},
        context={"run_id": run.name, "compound_id": "original-1", "source_node": "foreign-workstation"},
    )
    _replicate_group_entry(source, database, foreign_group_id, foreign_molecule_id, "foreign-workstation")
    foreign_replica_group_id = _id("peer-", "foreign-workstation:" + foreign_group_id)
    local_user = {"id": USER["id"], "role": "user", "_cluster_project_ids": []}

    selected = import_final_results(
        run, local_user, 0, database=database, selected_compound_ids=["original-1"],
    )
    with database.connection() as connection:
        local_rows = connection.execute(
            "SELECT count(*) FROM molecule_observations WHERE group_id=?", (selected["group_id"],),
        ).fetchone()[0]
        foreign_row = connection.execute(
            "SELECT metrics,context FROM molecule_observations WHERE group_id=?",
            (foreign_replica_group_id,),
        ).fetchone()
    selected_row = database.list_molecules(local_user, group_id=selected["group_id"])["items"][0]
    assert local_rows == 1
    assert selected_row["observations"]
    assert all(row["run_id"] == run.name for row in selected_row["observations"])
    assert selected["group_id"] in {group["id"] for group in selected_row["groups"]}
    assert all(row["metrics"].get("test_score") != -99 for row in selected_row["observations"])
    assert json.loads(foreign_row[0]) == {"test_score": -99}
    assert foreign_replica_group_id not in {group["id"] for group in selected_row["groups"]}
    with pytest.raises(PermissionError):
        database.group(foreign_replica_group_id, local_user)


def test_same_run_evaluation_from_another_workstation_is_a_distinct_source(bulk, tmp_path, monkeypatch):
    run, database = bulk
    from mn_molecular_design_studio import molecular_database_results
    from mn_molecular_design_studio.molecular_database import _json
    monkeypatch.setattr(molecular_database_results, "_authorized_run", lambda selected, _user: (selected,))
    source = MolecularDatabase(tmp_path / "other-workstation-source")
    imported = source.import_records(USER, [{
        "name": "Other workstation result",
        "components": [{"type": "small_molecule", "chain": "L", "smiles": "CCN"}],
    }], name="Distinct workstation source", origin={
        "type": "job", "collection_kind": "curated", "run_id": run.name,
        "source_node": "other-workstation",
    })
    molecule_id, group_id = imported["molecule_ids"][0], imported["group_id"]
    other_source_key = f"other-workstation/campaign/{run.name}:final:lib:1"
    with source.connection() as connection:
        connection.execute(
            "INSERT INTO molecule_sources VALUES (?,?,?,?)",
            (other_source_key, molecule_id, group_id, _json({
                "run_id": run.name, "compound_id": "original-1", "source_node": "other-workstation",
            })),
        )
    source.add_observation(group_id, molecule_id, "docking", run.name, {"test_score": -8.5},
                           context={"run_id": run.name, "compound_id": "original-1",
                                    "source_node": "other-workstation"})
    _replicate_group_entry(source, database, group_id, molecule_id, "other-workstation")

    selected = import_final_results(
        run, USER, 0, database=database, selected_compound_ids=["original-1"],
    )
    with database.connection() as connection:
        local_rows = connection.execute(
            "SELECT count(*) FROM molecule_observations WHERE group_id=?", (selected["group_id"],),
        ).fetchone()[0]
        other_group = connection.execute(
            "SELECT id FROM molecule_groups WHERE json_extract(origin,'$._replica_group')=?",
            (group_id,),
        ).fetchone()[0]
        other_evaluation = connection.execute(
            "SELECT metrics,context FROM molecule_observations WHERE group_id=?", (other_group,),
        ).fetchone()
    assert local_rows == 1
    assert json.loads(other_evaluation[0]) == {"test_score": -8.5}
    assert json.loads(other_evaluation[1])["source_node"] == "other-workstation"


def test_thousand_row_cap_does_not_read_or_import_the_rest(bulk):
    run, db = bulk
    rows = [dict(library_id="lib", compound_key=i, compound_id=f"ID-{i}",
                 smiles="C" + "".join("N" if (i >> bit) & 1 else "C" for bit in range(10)),
                 score=float(i), pose_path="") for i in range(1002)]
    path = run / "artifacts/active_learning/docking_scores.parquet"
    pq.write_table(pa.Table.from_pylist(rows), path, row_group_size=127)
    result = publish_final_selection(run, USER, 1000, database=db, explicit_selection=True)
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


@pytest.mark.parametrize("engine", ["vina", "gnina", "unidock_pro", "autodock_gpu", "rosetta_ligand"])
def test_omitted_docking_import_choice_defaults_to_results_only(tmp_path, monkeypatch, engine):
    monkeypatch.setattr(jobs, "runs_root", lambda: tmp_path / "runs")
    monkeypatch.setattr(workflow_queue, "spawn_worker_for_run", lambda _run: None)
    run = workflow_queue.enqueue_workflow_call(
        getattr(ligand_tools, "run_ligand_" + engine),
        task_group="ligand-docking", tool=engine, job_type="ligand_docking",
        kwargs={"target_pdb": tmp_path / "target.pdb", "compound_set_file": tmp_path / "compounds.json",
                "center_x": 0., "center_y": 0., "center_z": 0.},
        cpu_cores=1, gpu_device="none",
    )
    assert policy.policy_for_run(run) == "manual"
    assert not policy.automatic_collection_allowed(run)


def test_direct_docking_respects_import_switch(tmp_path, monkeypatch):
    monkeypatch.setattr(jobs, "runs_root", lambda: tmp_path)
    def fake(**kwargs):
        return jobs.create_job("ligand-docking", "ligand_docking", "vina", {}, {}).run_dir
    monkeypatch.setattr(ligand_tools, "_run_docking", fake)
    run = ligand_tools.run_ligand_vina(tmp_path / "p.pdb", tmp_path / "c.json", 0., 0., 0., import_to_database=False)
    assert policy.policy_for_run(run) == "manual"
    from mn_molecular_design_studio.molecular_database_ingest import ingest_run
    assert ingest_run(run, force=True) == 0


def test_direct_docking_without_import_choice_defaults_to_results_only(tmp_path, monkeypatch):
    monkeypatch.setattr(jobs, "runs_root", lambda: tmp_path)
    def fake(**kwargs):
        return jobs.create_job("ligand-docking", "ligand_docking", "vina", {}, {}).run_dir
    monkeypatch.setattr(ligand_tools, "_run_docking", fake)
    run = ligand_tools.run_ligand_vina(tmp_path / "p.pdb", tmp_path / "c.json", 0., 0., 0.)
    assert policy.policy_for_run(run) == "manual"


def test_batch_api_request_and_aggregate_completion_stay_results_only_by_default(tmp_path, monkeypatch):
    for key, value in {
        "MN_MOLECULAR_DESIGN_STUDIO_HOME": tmp_path / "studio",
        "MN_MOLECULAR_DESIGN_STUDIO_APP_HOME": tmp_path / "studio",
        "MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR": tmp_path / "runs",
        "MN_MOLECULAR_DESIGN_STUDIO_RESULTS_PARENT": tmp_path / "runs",
        "MN_MOLECULAR_DESIGN_STUDIO_REFERENCE_DIR": tmp_path / "references",
        "MDS_CLUSTER_HOME": tmp_path / "cluster",
    }.items():
        monkeypatch.setenv(key, str(value))

    from mn_molecular_design_studio.engine.workflow_engine.services import workflow_automation as automation
    from mn_molecular_design_studio.engine.workflow_engine.core import candidates
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs as jobs_module
    from mn_molecular_design_studio.engine.workflow_engine.core import workflow_queue as queue
    from mn_molecular_design_studio.engine.workflow_engine.workflows import input_batch
    from mn_molecular_design_studio.engine.workflow_engine.workflows import ligand_tools as tools
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    from mn_molecular_design_studio.workspace import mark_studio_job

    target, compounds = tmp_path / "target.pdb", tmp_path / "compounds.json"
    target.write_text("fixture")
    from mn_molecular_design_studio.ligand_compounds import write_compound_set
    write_compound_set(compounds, name="Fixture compounds", compounds=[
        {"compound_id": "one", "compound_name": "Ethanol", "smiles": "CCO"},
    ], owner_id=USER["id"], project_id="", created_at="fixture")
    monkeypatch.setattr(automation, "_validate_workflow_inputs", lambda *args, **kwargs: None)
    monkeypatch.setattr(queue, "spawn_worker_for_run", lambda *_args, **_kwargs: None)

    request = tmp_path / "request.json"
    items = [
        {"target_pdb": str(target), "compound_set_file": str(compounds), "center_x": 0.,
         "center_y": 0., "center_z": 0., "exhaustiveness": value}
        for value in (8, 16)
    ]
    request.write_text(json.dumps({"schema_version": 1, "workflow": "ligand.docking.vina",
                                   "parameters": items,
                                   "resources": {"cpu_cores": 1, "gpu_device": "none"}}))
    submitted = automation.submit_workflow_request(request)
    parent = Path(submitted["job"]["run_dir"])
    assert policy.policy_for_run(parent) == "manual"
    mark_studio_job(submitted, owner_id=USER["id"], project_id="", allow_private_source=True)

    child_policies = []
    def fake_engine(**_kwargs):
        child = jobs_module.create_job("ligand-docking", "ligand_docking", "vina", {}, {}).run_dir
        child_policies.append(policy.policy_for_run(child))
        candidates.write_candidates(child, "vina", [{
            "candidate_id": "one", "stage": "analysis", "source_tool": "vina",
            "metrics": {"score": -7.5},
            "raw_metadata": {"compound_id": "one", "components": [
                {"chain": "L", "type": "small_molecule", "smiles": "CCO"},
            ]},
        }])
        jobs_module.finish_job(child, True, {})
        return child

    monkeypatch.setattr(tools, "_run_docking", fake_engine)
    ingestion_attempts = []
    def observe_ingestion(run_dir, **kwargs):
        ingestion_attempts.append(Path(run_dir))
        return 0
    import mn_molecular_design_studio.molecular_database_ingest as ingest_module
    monkeypatch.setattr(ingest_module, "ingest_run", observe_ingestion)
    with jobs_module.adopt_next_job_creation(jobs_module.JobPaths("ligand-docking", parent.name, parent)):
        completed_parent = input_batch.run_input_batch(
            "ligand.docking.vina", items, cpu_cores=1, gpu_device="none",
        )

    assert completed_parent == parent
    assert child_policies == ["manual", "manual"]
    assert policy.policy_for_run(parent) == "manual"
    assert ingestion_attempts == []
    database = MolecularDatabase()
    with database.connection() as connection:
        for table in ("molecules", "molecule_members", "molecule_sources", "molecule_observations"):
            assert connection.execute(f"SELECT count(*) FROM {table}").fetchone()[0] == 0


def test_direct_publisher_requires_explicit_positive_top_n_action(bulk):
    run, database = bulk
    with pytest.raises(ValueError, match="explicit curation action"):
        publish_final_selection(run, USER, 1, database=database)
    with database.connection() as connection:
        assert connection.execute("SELECT count(*) FROM molecule_sources").fetchone()[0] == 0


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
