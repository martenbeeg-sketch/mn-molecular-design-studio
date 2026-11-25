"""API-only journeys use isolated accounts, files and tiny scientific fixtures."""
from __future__ import annotations

import asyncio
import base64
import copy
import hashlib
import json
import threading
from pathlib import Path

import httpx
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from fastapi.testclient import TestClient


@pytest.fixture
def queued_vina_batch(api_data, monkeypatch):
    """Use authenticated API submission and the real worker with tiny Docker outputs."""
    from mn_molecular_design_studio import ligand_compounds
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs, workflow_queue
    from mn_molecular_design_studio.engine.workflow_engine.workflows import ligand_tools
    from mn_molecular_design_studio.molecular_database import MolecularDatabase

    actor, project = api_data["admin"], api_data["project"]["id"]
    target, compounds = api_data["root"] / "target.pdb", api_data["root"] / "compounds.json"
    target.write_text("ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C\nEND\n")
    ligand_compounds.write_compound_set(
        compounds, name="Owned library", compounds=[{"compound_id": "ethanol", "smiles": "CCO"}],
        owner_id=actor["id"], project_id=project, created_at="fixture",
    )
    database = MolecularDatabase()
    intentional = database.import_records(actor, [{
        "name": "Intentional ethanol", "components": [{"type": "small_molecule", "chain": "L", "smiles": "CCO"}],
    }], name="Intentional import", project_id=project)
    launches, scores, pending = [], {}, []
    monkeypatch.setattr(ligand_tools, "load_manifest", lambda _name: {"image": "fixture:no-container"})
    monkeypatch.setattr(workflow_queue, "spawn_worker_for_run", lambda run: pending.append(Path(run)))

    def docker(run, _manifest, _command, *, gpu_device):
        run = Path(run)
        launches.append(run)
        params = jobs.read_json(run / "input.json")["params"]
        staged = ligand_compounds.load_compound_set(run / "artifacts/input/compound_set.json")["compounds"]
        pose_dir = run / "artifacts/poses"
        pose_dir.mkdir(exist_ok=True)
        rows = []
        for index, compound in enumerate(staged):
            compound_id = compound["compound_id"]
            score = -float(params["exhaustiveness"]) - (index + 1) / 4
            scores[(run.name, compound_id)] = score
            (pose_dir / f"{compound_id}.pdbqt").write_text(
                f"REMARK fixture pose {run.name} {compound_id} {score}\nENDMDL\n"
            )
            rows.append(f"{compound_id}\tcompleted\t{score}\tartifacts/poses/{compound_id}.pdbqt\t\n")
        (run / "artifacts/docking_scores.csv").write_text(
            "compound_id\tstatus\tbest_score_kcal_mol\tpose_path\terror\n" + "".join(rows)
        )
        return 0

    monkeypatch.setattr(ligand_tools, "run_docker_job", docker)

    def body(choices, compatible=False):
        items = []
        for index, choice in enumerate(choices):
            item = dict(target_pdb=str(target), compound_set_file=str(compounds),
                        center_x=0., center_y=0., center_z=0.,
                        exhaustiveness=8 if compatible else 8 + index)
            if choice is not None:
                item["import_to_database"] = choice
            items.append(item)
        return {"schema_version": 1, "workflow": "ligand.docking.vina", "project_id": project,
                "parameters": items, "resources": {"cpu_cores": 1, "gpu_device": "none"}}

    def observations():
        with database.connection() as connection:
            return [dict(row) for row in connection.execute(
                "SELECT o.*,g.owner_id,g.project_id FROM molecule_observations o "
                "JOIN molecule_groups g ON g.id=o.group_id WHERE namespace='vina' ORDER BY o.id"
            )]

    return dict(data=api_data, database=database, actor=actor, project=project, body=body,
                launches=launches, scores=scores, pending=pending, observations=observations,
                canonical_id=intentional["molecule_ids"][0])


@pytest.mark.parametrize("choices", [
    (None, None), (False, False), (True, True), (True, False), (False, True),
    (True, None), (None, True),
])
@pytest.mark.parametrize("compatible", [False, True])
def test_vina_batch_worker_before_api_marker_retains_exact_imports(queued_vina_batch, monkeypatch, choices, compatible):
    """Start the real worker before API marking, then finish and retry its queued call."""
    from mn_molecular_design_studio.collection_policy import policy_for_run
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs, local_worker, workflow_queue
    from mn_molecular_design_studio.engine.workflow_engine.workflows import ligand_tools
    from mn_molecular_design_studio.molecular_database_ingest import record_completed_run

    case = queued_vina_batch
    entered_docker, release_docker = threading.Event(), threading.Event()
    worker_thread = None
    parents, errors, results = [], [], []
    fixture_docker = ligand_tools.run_docker_job

    def paused_docker(run, manifest, command, *, gpu_device):
        if not entered_docker.is_set():
            entered_docker.set()
            assert release_docker.wait(10), "API did not release the isolated worker"
        return fixture_docker(run, manifest, command, gpu_device=gpu_device)

    def worker(run):
        try:
            results.append(local_worker.run_worker_job(run))
        except BaseException as exc:
            errors.append(exc)

    def dispatch(run):
        nonlocal worker_thread
        run = Path(run)
        parents.append(run)
        assert not (run / "molecular_design_studio.json").exists()
        metadata = jobs.read_json(run / "metadata.json")
        assert metadata["studio_owner_id"] == case["actor"]["id"]
        assert metadata["studio_project_id"] == case["project"]
        # A fresh thread gives the worker its own ContextVars, as a subprocess
        # would; submission's parent policy must not constrain true children.
        worker_thread = threading.Thread(target=worker, args=(run,), daemon=True)
        worker_thread.start()
        assert entered_docker.wait(10), "Worker did not enter the Docker substitute"

    monkeypatch.setattr(ligand_tools, "run_docker_job", paused_docker)
    monkeypatch.setattr(workflow_queue, "spawn_worker_for_run", dispatch)
    try:
        response = case["data"]["client"].post("/api/v1/workflows", json=case["body"](choices, compatible))
        assert response.status_code == 202, response.text
        parent = parents[0]
        assert (parent / "molecular_design_studio.json").exists()
        assert worker_thread.is_alive()  # Marking completes before either child finishes.
    finally:
        release_docker.set()
        if worker_thread is not None:
            worker_thread.join(10)
    assert not worker_thread.is_alive()
    assert not errors, errors
    assert results == [0]
    should_coalesce = compatible and choices[0] == choices[1]
    assert len(case["launches"]) == (1 if should_coalesce else 2)
    launch_choices = [choices[0]] if should_coalesce else choices
    assert [policy_for_run(run) for run in case["launches"]] == [
        "automatic" if choice is True else "manual" for choice in launch_choices
    ]
    assert policy_for_run(parent) == ("automatic" if should_coalesce and choices[0] is True else "manual")
    observations = case["observations"]()
    assert len(observations) == sum(choice is True for choice in choices)
    for run in case["launches"]:
        marker = jobs.read_json(run / "molecular_design_studio.json")
        assert marker["studio_owner_id"] == case["actor"]["id"]
        assert marker["project_id"] == case["project"]

    for observation in observations:
        context = json.loads(observation["context"])
        compound_id = context["row"]["compound_id"]
        assert observation["molecule_id"] == case["canonical_id"]
        assert observation["owner_id"] == case["actor"]["id"]
        assert observation["project_id"] == case["project"]
        assert context["run_id"] == observation["run_id"]
        score = case["scores"][(observation["run_id"], compound_id)]
        assert json.loads(observation["metrics"])["best_score_kcal_mol"] == score
        launch = next(run for run in case["launches"] if run.name == observation["run_id"])
        assert float(context["row"]["best_score_kcal_mol"]) == score
        assert context["settings"]["exhaustiveness"] == jobs.read_json(launch / "input.json")["params"]["exhaustiveness"]
        with case["database"].connection() as connection:
            structure = connection.execute("SELECT * FROM molecule_structures WHERE id=?", (context["representation_id"],)).fetchone()
            sources = connection.execute(
                "SELECT source_key FROM molecule_sources WHERE group_id=? AND molecule_id=? "
                "AND json_extract(provenance,'$.run_id')=? AND json_extract(provenance,'$.compound_id')=?",
                (observation["group_id"], case["canonical_id"], observation["run_id"], compound_id),
            ).fetchall()
        pose_bytes = (launch / context["row"]["pose_path"]).read_bytes()
        assert Path(structure["path"]).read_bytes() == pose_bytes
        assert context["pose_digest"] == hashlib.sha256(pose_bytes).hexdigest()
        assert sources and all(observation["run_id"] in row["source_key"] for row in sources)

    launches_before_retry = list(case["launches"])
    retained_files = {path: path.read_bytes() for run in case["launches"]
                      for path in [run / "artifacts/docking_scores.csv", *run.glob("artifacts/poses/*.pdbqt")]}
    before = copy.deepcopy(observations)
    record_completed_run(parent)
    assert case["observations"]() == before
    assert local_worker.run_worker_job(parent) == 0  # Restore worker_request.json and retry the entire batch.
    assert case["observations"]() == before
    # The scalar coalesced adapter reruns deterministically; completed separate
    # children are resumed without executing Docker again.
    assert case["launches"] == (launches_before_retry * 2 if should_coalesce else launches_before_retry)
    assert all(path.read_bytes() == content for path, content in retained_files.items())


@pytest.mark.parametrize("key,value", [
    ("owner_id", "foreign-owner"), ("studio_owner_id", "foreign-owner"),
    ("project_id", "foreign-project"), ("studio_project_id", "foreign-project"),
    ("molecular_design_studio", {"studio_owner_id": "foreign-owner"}),
    ("studio_lineage", {"studio_owner_id": "foreign-owner"}),
    ("database_collection_policy", "automatic"),
])
def test_vina_batch_parameters_cannot_forge_trusted_lineage(queued_vina_batch, key, value):
    case = queued_vina_batch
    body = case["body"]((True, False))
    body["parameters"][0][key] = value
    response = case["data"]["client"].post("/api/v1/workflows", json=body)
    assert response.status_code == 422, response.text
    assert not case["pending"] and not case["launches"]
    assert not case["observations"]()


@pytest.mark.parametrize("admit", [False, True])
def test_vina_opt_in_admits_new_scored_ligand_and_pose(queued_vina_batch, monkeypatch, admit):
    from mn_molecular_design_studio import ligand_compounds
    from mn_molecular_design_studio.engine.workflow_engine.core import local_worker
    from mn_molecular_design_studio.engine.workflow_engine.workflows import ligand_tools

    case = queued_vina_batch
    body = case["body"]((admit,))
    body["parameters"] = body["parameters"][0]
    ligand_compounds.write_compound_set(
        Path(body["parameters"]["compound_set_file"]), name="Unregistered ligands",
        compounds=[{"compound_id":"propanol","smiles":"CCCO"},
                   {"compound_id":"butanol","smiles":"CCCCO"},
                   {"compound_id":"pentanol","smiles":"CCCCCO"}],
        owner_id=case["actor"]["id"], project_id=case["project"], created_at="fixture",
    )
    fixture_docker = ligand_tools.run_docker_job

    def docker(run, manifest, command, *, gpu_device):
        result = fixture_docker(run, manifest, command, gpu_device=gpu_device)
        path = Path(run) / "artifacts/docking_scores.csv"
        path.write_text(path.read_text().replace("butanol\tcompleted", "butanol\tfailed"))
        (Path(run) / "artifacts/poses/pentanol.pdbqt").unlink()
        return result

    monkeypatch.setattr(ligand_tools, "run_docker_job", docker)
    response = case["data"]["client"].post("/api/v1/workflows", json=body)
    assert response.status_code == 202, response.text
    assert local_worker.run_worker_job(case["pending"][0]) == 0
    observations = case["observations"]()
    assert len(observations) == int(admit)
    if admit:
        context = json.loads(observations[0]["context"])
        assert context["row"]["compound_id"] == "propanol"
        assert context["calculation_recorded"] and context["representation_id"]
        with case["database"].connection() as connection:
            stored = connection.execute("SELECT path FROM molecule_structures WHERE id=?",
                                        (context["representation_id"],)).fetchone()
        assert Path(stored["path"]).read_bytes() == (
            case["launches"][0] / "artifacts/poses/propanol.pdbqt").read_bytes()


@pytest.fixture
def api_data(tmp_path, monkeypatch):
    for name, value in {
        "MN_MOLECULAR_DESIGN_STUDIO_HOME": tmp_path / "studio",
        "MN_MOLECULAR_DESIGN_STUDIO_APP_HOME": tmp_path / "studio",
        "MN_MOLECULAR_DESIGN_STUDIO_RESULTS_PARENT": tmp_path / "results",
        "MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR": tmp_path / "results/runs",
        "MN_MOLECULAR_DESIGN_STUDIO_REFERENCE_DIR": tmp_path / "references",
        "MDS_CLUSTER_HOME": tmp_path / "cluster",
        "MN_MOLECULAR_DESIGN_STUDIO_AUTH_MODE": "production",
    }.items():
        monkeypatch.setenv(name, str(value))
    from mn_molecular_design_studio import auth, api, workspace, molecular_database_api
    from mn_molecular_design_studio.cluster import config
    workspace.save_settings({"results_directory": str(tmp_path / "results/runs"),
                             "reference_directory": str(tmp_path / "references")})
    workspace.ensure_studio_workspace()
    config.save({"mode": "disabled"})
    monkeypatch.setattr(molecular_database_api, "sync_workspace_async", lambda *a, **k: [])
    admin = auth.create_initial_admin("research_instructor", "Instructor", "instructor-password-42")
    project = workspace.create_project("API course", owner_id=admin["id"])
    student = auth.create_user("student", "Student", "temporary-password-42", actor=admin)
    auth.change_password(student["id"], "student-password-42", forced=True)
    auth.grant_project_access(project["id"], student["id"], actor=admin)
    other = auth.create_user("otherstudent", "Other", "temporary-password-42", actor=admin)
    auth.change_password(other["id"], "otherstudent-password-42", forced=True)
    admin_token = auth.authenticate("research_instructor", "instructor-password-42")[1]
    student_token = auth.authenticate("student", "student-password-42")[1]
    other_token = auth.authenticate("otherstudent", "otherstudent-password-42")[1]
    with TestClient(api.app) as client:
        client.headers["Authorization"] = "Bearer " + admin_token
        yield {"client": client, "app": api.app, "admin": admin, "student": student,
               "other": other, "project": project, "root": tmp_path,
               "admin_token": admin_token, "student_token": student_token, "other_token": other_token}


def headers(data, role):
    return {"Authorization": "Bearer " + data[role + "_token"]}


def import_compounds(data):
    response = data["client"].post("/api/v1/molecules/imports", json={
        "method": "smiles", "text": "CCO ethanol\nCCN ethylamine", "project_id": data["project"]["id"]})
    assert response.status_code == 201, response.text
    return response.json()


def test_database_and_results_reads_do_not_trigger_ingestion(api_data, monkeypatch):
    from mn_molecular_design_studio import molecular_database_api

    monkeypatch.setattr(molecular_database_api, "sync_workspace_async",
                        lambda *args, **kwargs: pytest.fail("ordinary reads must not trigger ingestion"))
    data, client = api_data, api_data["client"]
    imported = import_compounds(data)
    group_id = imported["group_id"]
    molecule_id = imported["molecule_ids"][0]
    assert client.get("/api/v1/molecules/groups").status_code == 200
    assert client.get("/api/v1/molecules/metrics").status_code == 200
    assert client.get("/api/v1/molecules").status_code == 200
    assert client.get(f"/api/v1/molecules/groups/{group_id}").status_code == 200
    assert client.get(f"/api/v1/molecules/{molecule_id}").status_code == 200

    # The Results endpoint lists accessible run data without asking the
    # background workspace refresh to ingest it as a side effect.
    screening_run(data)
    result = client.get("/api/v1/molecules/jobs/campaign/api-screening")
    assert result.status_code == 200, result.text


def screening_run(data, count=11):
    from mn_molecular_design_studio.runtime import runs_root
    from mn_molecular_design_studio.workspace import STUDIO_SOURCE_ID
    run = Path(runs_root()) / "campaign" / "api-screening"
    folder = run / "artifacts/active_learning"
    folder.mkdir(parents=True)
    (run / "molecular_design_studio.json").write_text(json.dumps({"application": STUDIO_SOURCE_ID,
        "studio_owner_id": data["admin"]["id"], "project_id": data["project"]["id"]}))
    (run / "metadata.json").write_text(json.dumps({"status": "completed", "database_collection_policy": "manual"}))
    (run / "result.json").write_text(json.dumps({"success": True}))
    pq.write_table(pa.Table.from_pylist([{"compound_id": f"c-{i}", "compound_name": f"Compound {i}",
        "compound_key": i, "library_id": "lib", "smiles": "C" * (i % 8 + 1), "score": -20. + i / 100}
        for i in range(count)]), folder / "docking_scores.parquet", row_group_size=3)
    (folder / "round_history.json").write_text(json.dumps({"rounds": [{"round": 1}], "score_units": "REU"}))
    return run


def test_private_upload_retains_extension_and_never_imports_compounds(api_data):
    from mn_molecular_design_studio import api_files
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    data, client = api_data, api_data["client"]
    response = client.post("/api/v1/files?filename=library.csv", content=b"id,smiles\na,CCO\n")
    assert response.status_code == 201, response.text
    record = response.json()
    path = api_files.resolve(data["admin"], record["id"])
    assert path.suffix == ".csv" and path.read_bytes() == b"id,smiles\na,CCO\n"
    assert client.get(record["download_url"]).content == path.read_bytes()
    assert client.get(record["download_url"], headers=headers(data, "student")).status_code == 404
    assert client.get(record["download_url"] + "/metadata", headers=headers(data, "student")).status_code == 404
    assert MolecularDatabase().list_molecules(data["admin"])["total"] == 0
    forged = client.post("/api/v1/workflows/validate", headers=headers(data, "student"), json={
        "project_id": data["project"]["id"], "schema_version": 1, "workflow": "library.prepare",
        "parameters": {"source_files": [record["reference"]]}})
    assert forged.status_code == 404


def test_library_registration_sharing_and_lineage_are_metadata_only(api_data):
    from mn_molecular_design_studio.external_libraries import LibraryCatalog, FEATURE_SCHEMA
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    data, client = api_data, api_data["client"]
    source = data["root"] / "large-library.csv"
    source.write_text("ID,SMILES\na,CCO\nb,CCN\n")
    request = {"paths": [str(source)], "name": "External library", "id_column": "ID", "smiles_column": "SMILES"}
    assert client.post("/api/v1/libraries", json=request, headers=headers(data, "student")).status_code == 403
    response = client.post("/api/v1/libraries", json=request)
    assert response.status_code == 201, response.text
    library = response.json()
    assert client.get("/api/v1/libraries", headers=headers(data, "student")).json()["items"] == []
    assert client.get("/api/v1/libraries/" + library["id"], headers=headers(data, "student")).status_code == 404
    shared = client.patch("/api/v1/libraries/" + library["id"] + "/access", json={"project_ids": [data["project"]["id"]]})
    assert shared.status_code == 200, shared.text
    manifest = data["root"] / "fingerprints.json"
    manifest.write_text(json.dumps({"schema_version": FEATURE_SCHEMA, "library_id": library["id"],
        "status": "completed", "row_count": 2, "created_at": "2026-10-26T00:00:00Z"}))
    LibraryCatalog().attach_preparation(library["id"], manifest)
    path = "/api/v1/libraries/" + library["id"] + "/preparations?kind=fingerprints"
    assert client.get(path, headers=headers(data, "student")).json()["total"] == 1
    assert client.get(path, headers=headers(data, "other")).status_code == 404
    assert client.post("/api/v1/libraries", json=request).json()["project_ids"] == [data["project"]["id"]]
    assert MolecularDatabase().list_molecules(data["admin"])["total"] == 0


def test_screening_pagination_filters_original_rank_and_streamed_download(api_data, monkeypatch):
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    data, client = api_data, api_data["client"]
    screening_run(data, count=4101)
    monkeypatch.setattr(pq, "read_table", lambda *a, **k: pytest.fail("Must use bounded Parquet batches"))
    path = "/api/v1/screening/campaign/api-screening"
    response = client.post(path + "/search", json={"offset": 2047, "limit": 5})
    assert response.status_code == 200, response.text
    page = response.json()
    assert [r["rank"] for r in page["items"]] == list(range(2048, 2053))
    assert page["total"] == 4101 and page["has_more"] and page["next_offset"] == 2052
    filtered = client.post(path + "/search", json={"compound_ids": ["c-2", "c-2000", "c-4100"], "offset": 1, "limit": 2}).json()
    assert [r["rank"] for r in filtered["items"]] == [2001, 4101]
    assert filtered["total"] is None and filtered["has_more"] is False
    assert client.post(path + "/search", json={"score_min": 10, "score_max": -10}).status_code == 422
    assert client.post(path + "/search", json={"limit": 1001}).status_code == 422
    assert client.post(path + "/search", json={}, headers=headers(data, "other")).status_code == 404
    assert client.get(path + "/files/scores").content.startswith(b"PAR1")
    assert client.get(path + "/files/secrets").status_code == 404
    assert client.get(path).json()["score_units"] == "REU"
    assert MolecularDatabase().list_molecules(data["admin"])["total"] == 0


def test_final_screening_import_is_explicit_bounded_and_owner_only(api_data):
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    data, client = api_data, api_data["client"]
    run = screening_run(data)
    path = "/api/v1/screening/campaign/api-screening/final-selection"
    assert client.post(path, json={"limit": 1001}).status_code == 422
    assert client.post(path, json={"limit": 2}, headers=headers(data, "student")).status_code == 404
    (run / "metadata.json").write_text(json.dumps({"status": "running"}))
    assert client.post(path, json={"limit": 2}).status_code == 422
    (run / "metadata.json").write_text(json.dumps({"status": "completed"}))
    result = client.post(path, json={"limit": 2})
    assert result.status_code == 200, result.text
    assert result.json()["imported_count"] == 2
    assert client.post(path, json={"limit": 2}).json()["imported_count"] == 2
    assert MolecularDatabase().list_molecules(data["admin"])["total"] == 2
    assert pq.ParquetFile(run / "artifacts/active_learning/docking_scores.parquet").metadata.num_rows == 11


def test_final_selection_api_without_a_choice_is_results_only(api_data):
    from mn_molecular_design_studio.molecular_database import MolecularDatabase

    data = api_data
    run = screening_run(data)
    response = data["client"].post("/api/v1/screening/campaign/api-screening/final-selection", json={})
    assert response.status_code == 200, response.text
    assert response.json()["requested_limit"] == 0
    assert response.json()["imported_count"] == 0
    assert response.json()["group_id"] is None
    with MolecularDatabase().connection() as connection:
        assert connection.execute("SELECT count(*) FROM molecule_sources").fetchone()[0] == 0
        assert connection.execute("SELECT count(*) FROM molecule_observations").fetchone()[0] == 0
    assert pq.ParquetFile(run / "artifacts/active_learning/docking_scores.parquet").metadata.num_rows == 11


def test_final_selection_api_accepts_exact_ids_and_replays_idempotently(api_data):
    from mn_molecular_design_studio.molecular_database import MolecularDatabase

    data = api_data
    run = screening_run(data)
    scores_path = run / "artifacts/active_learning/docking_scores.parquet"
    original_scores = scores_path.read_bytes()
    endpoint = "/api/v1/screening/campaign/api-screening/final-selection"
    body = {"limit": 0, "selected_compound_ids": ["c-1", "c-3"]}
    first = data["client"].post(endpoint, json=body)
    assert first.status_code == 200, first.text
    assert first.json()["selected_compound_ids"] == ["c-1", "c-3"]
    assert first.json()["imported_count"] == 2
    database = MolecularDatabase()
    items = database.list_molecules(data["admin"], group_id=first.json()["group_id"])["items"]
    assert {item["metadata"]["compound_id"] for item in items} == {"c-1", "c-3"}

    replay = data["client"].post(endpoint, json=body)
    assert replay.status_code == 200, replay.text
    assert replay.json()["group_id"] == first.json()["group_id"]
    with database.connection() as connection:
        assert connection.execute("SELECT count(*) FROM molecule_members").fetchone()[0] == 2
        assert connection.execute("SELECT count(*) FROM molecule_sources").fetchone()[0] == 2
        assert connection.execute("SELECT count(*) FROM molecule_observations").fetchone()[0] == 2
    assert scores_path.read_bytes() == original_scores

    for invalid in (
        {"limit": 0, "selected_compound_ids": []},
        {"limit": 1, "selected_compound_ids": ["c-1"]},
        {"limit": 0, "selected_compound_ids": ["missing"]},
    ):
        assert data["client"].post(endpoint, json=invalid).status_code == 422
    with database.connection() as connection:
        assert connection.execute("SELECT count(*) FROM molecule_members").fetchone()[0] == 2


@pytest.mark.parametrize("format,extension", [("CSV", "csv"), ("SMILES", "smi"), ("SDF", "sdf")])
def test_exports_are_private_downloadable_files(api_data, format, extension):
    data, client = api_data, api_data["client"]
    imported = import_compounds(data)
    response = client.post("/api/v1/molecules/exports", json={"molecule_ids": imported["molecule_ids"], "format": format})
    assert response.status_code == 201, response.text
    exported = response.json()
    assert exported["name"].endswith("." + extension)
    assert client.get(exported["download_url"]).status_code == 200
    assert client.get(exported["download_url"], headers=headers(data, "other")).status_code == 404
    assert client.post("/api/v1/molecules/exports", json={"molecule_ids": imported["molecule_ids"], "format": format}, headers=headers(data, "other")).status_code == 404


def test_reviewed_table_import_persists_report_and_original_source(api_data):
    data, client = api_data, api_data["client"]
    content = b"ID,SMILES,Name\na,CCO,Ethanol\nb,OCC,Duplicate\nc,invalid,Bad\n"
    upload = client.post("/api/v1/files?filename=molecules.csv", content=content).json()
    analyzed = client.post("/api/v1/molecules/imports/analyze", json={"file_id": upload["id"], "smiles_column": "SMILES", "id_column": "ID"})
    assert analyzed.status_code == 200, analyzed.text
    assert analyzed.json()["summary"]["duplicate_structure_count"] == 1
    imported = client.post("/api/v1/molecules/imports", json={"method": "table", "file_id": upload["id"],
        "review_small_molecules": True, "chain_columns": {"SMILES": {"type": "small_molecule", "chain": "L"}},
        "id_column": "ID", "name_column": "Name"})
    assert imported.status_code == 201, imported.text
    assert len(imported.json()["molecule_ids"]) == 1
    batches = client.get("/api/v1/molecules/groups/" + imported.json()["group_id"] + "/imports").json()["items"]
    assert len(batches) == 1
    report = "/api/v1/molecules/imports/" + batches[0]["id"]
    assert client.get(report + "/rows?limit=2&offset=1").json()["total"] == 3
    assert client.get(report + "/source").content == content
    assert client.get(report + "/source", headers=headers(data, "other")).status_code == 404


def test_group_role_crop_and_membership_operations_use_existing_identities(api_data):
    data, client = api_data, api_data["client"]
    imported = client.post("/api/v1/molecules/imports", json={"method": "records", "records": [{
        "name": "Two chains", "components": [{"type": "protein", "chain": "A", "sequence": "ACDEFG"},
                                              {"type": "protein", "chain": "B", "sequence": "KLMNPQ"}]}]}).json()
    molecule, group = imported["molecule_ids"][0], imported["group_id"]
    assert client.patch(f"/api/v1/molecules/groups/{group}/molecules/{molecule}/chains/A/role", json={"role": "target"}).status_code == 200
    result = client.post(f"/api/v1/molecules/{molecule}/crop", json={"chains": ["A"]})
    assert result.status_code == 201, result.text
    cropped = client.get("/api/v1/molecules/" + result.json()["molecule_ids"][0]).json()
    assert [c["chain"] for c in cropped["components"]] == ["A"]
    assert client.patch(f"/api/v1/molecules/groups/{group}/project", json={"project_id": data["project"]["id"]}).status_code == 200
    assert client.post(f"/api/v1/molecules/groups/{group}/members/remove", json={"molecule_ids": [molecule]}, headers=headers(data, "student")).status_code == 404
    assert client.post(f"/api/v1/molecules/groups/{group}/members/remove", json={"molecule_ids": [molecule]}).status_code == 200
    # Removing one membership never removes another saved collection.
    assert client.get("/api/v1/molecules/" + cropped["id"]).status_code == 200


def test_account_administration_grants_disable_and_password_reset(api_data):
    data, client = api_data, api_data["client"]
    path = "/api/v1/admin/users"
    assert client.get(path, headers=headers(data, "student")).status_code == 403
    created = client.post(path, json={"username": "newstudent", "display_name": "New Student",
        "temporary_password": "temporary-password-42", "project_ids": [data["project"]["id"]]})
    assert created.status_code == 201, created.text
    account = created.json()
    assert "temporary-password-42" not in created.text and "password_hash" not in account
    assert account["project_ids"] == [data["project"]["id"]] and account["must_change_password"]
    member = path + "/" + account["id"]
    assert client.put(member + "/projects", json={"project_ids": ["not-a-project"]}).status_code == 422
    assert client.get(member).json()["project_ids"] == account["project_ids"]
    assert client.put(member + "/projects", json={"project_ids": []}).json()["project_ids"] == []
    assert client.post(member + "/password", json={"temporary_password": "changed-password-42"}).json()["sessions_revoked"]
    assert client.patch(member, json={"disabled": True}).json()["disabled"]
    assert client.patch(path + "/" + data["admin"]["id"], json={"disabled": True}).status_code == 422
    assert client.patch(member, json={"disabled": "yes"}).status_code == 422


def test_codex_pairing_writes_as_the_signed_in_account_and_revokes_separately(api_data):
    from concurrent.futures import ThreadPoolExecutor

    data, client = api_data, api_data["client"]
    pair = "/api/v1/auth/mcp/pair"
    redeem = "/api/v1/auth/mcp/redeem"
    disconnect = "/api/v1/auth/mcp/disconnect"
    assert client.post(pair, headers={"Authorization": ""}).status_code == 401
    response = client.post(pair)
    assert response.status_code == 200 and response.headers["Cache-Control"] == "no-store"
    code = response.json()["code"]
    with ThreadPoolExecutor(max_workers=2) as pool:
        attempts = list(pool.map(lambda _: client.post(redeem, json={"code": code}), range(2)))
    assert sorted(row.status_code for row in attempts) == [200, 401]
    token = next(row.json()["access_token"] for row in attempts if row.status_code == 200)
    paired = {"Authorization": "Bearer " + token}
    identity = client.get("/api/v1/auth/me", headers=paired).json()
    assert identity["id"] == data["admin"]["id"]
    assert identity["mcp_session"] is True and identity["session_expires_at"] is None
    project = client.post("/api/v1/projects", json={"name": "Created via Codex"}, headers=paired)
    assert project.status_code == 201, project.text
    assert project.json()["created_by"] == data["admin"]["id"]
    assert client.get("/api/v1/admin/users", headers=paired).status_code == 403
    unrestricted = client.post("/api/v1/workflows/validate", json={"project_id": data["project"]["id"]},
                               params={"base_dir": str(data["root"])}, headers=paired)
    assert unrestricted.status_code == 403
    assert client.post(pair, headers=paired).status_code == 403
    assert client.post(disconnect, headers=paired).status_code == 403
    pending = client.post(pair).json()["code"]
    assert client.post(disconnect, headers=headers(data, "student")).status_code == 200
    assert client.get("/api/v1/auth/me", headers=paired).status_code == 200
    assert client.post(disconnect).status_code == 200
    assert client.get("/api/v1/auth/me", headers=paired).status_code == 401
    assert client.post(redeem, json={"code": pending}).status_code == 401
    assert client.get("/api/v1/auth/me").status_code == 200


def test_codex_pairing_expires_but_redeemed_access_survives_browser_session(api_data, monkeypatch):
    from datetime import timedelta
    from mn_molecular_design_studio import auth

    data, client = api_data, api_data["client"]
    pair = "/api/v1/auth/mcp/pair"
    redeem = "/api/v1/auth/mcp/redeem"
    code = client.post(pair).json()["code"]
    original_now = auth._now()
    with monkeypatch.context() as clock:
        clock.setattr(auth, "_now", lambda: original_now + timedelta(minutes=6))
        assert client.post(redeem, json={"code": code}).status_code == 401
    code = client.post(pair).json()["code"]
    token = client.post(redeem, json={"code": code}).json()["access_token"]
    paired = {"Authorization": "Bearer " + token}
    unredeemed = client.post(pair).json()["code"]
    with monkeypatch.context() as clock:
        clock.setattr(auth, "_now", lambda: original_now + timedelta(days=30))
        assert client.get("/api/v1/auth/me").status_code == 401
        assert client.get("/api/v1/auth/me", headers=paired).status_code == 200
    auth.revoke_session(data["admin_token"])
    assert client.post(redeem, json={"code": unredeemed}).status_code == 401
    assert client.get("/api/v1/auth/me", headers=paired).status_code == 200
    auth.reset_admin_password("research_instructor", "new-password-42")
    assert client.get("/api/v1/auth/me", headers=paired).status_code == 401
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_AUTH_MODE", "development")
    assert client.post(pair).status_code == 403
    assert client.post(redeem, json={"code": code}).status_code == 403


def test_capabilities_exclude_node_routes_and_nonadmin_management(api_data):
    client = api_data["client"]
    response = client.get("/api/v1/capabilities", headers=headers(api_data, "student"))
    assert response.status_code == 200
    assert response.json()["final_database_import_limit"] == 0
    paths = [r["path"] for r in response.json()["operations"]]
    assert "/api/v1/libraries" in paths
    assert not any("/cluster/node/" in p or "/admin/" in p for p in paths)
    assert client.get("/api/v1/auth/me").json()["development_mode"] is False
    operation = client.get("/api/v1/operations/schema", params={"path": "/api/v1/molecules/exports", "method": "POST"}, headers=headers(api_data, "student"))
    assert operation.status_code == 200 and "ExportRequest" in operation.json()["components"]["schemas"]
    assert client.get("/api/v1/operations/schema", params={"path": "/api/v1/admin/users", "method": "POST"}, headers=headers(api_data, "student")).status_code == 404
    assert client.get("/api/v1/operations/schema", params={"path": "/api/v1/cluster/node/claim", "method": "POST"}).status_code == 404


def test_local_request_retries_never_submit_twice(api_data, monkeypatch):
    from mn_molecular_design_studio import api
    client = api_data["client"]
    calls = []
    monkeypatch.setattr(api, "_cluster_submission", lambda *a, **k: None)
    monkeypatch.setattr(api, "_guard_source_job", lambda *a, **k: None)
    monkeypatch.setattr(api, "submit_workflow_request", lambda path: calls.append(path) or {"job": {"run_id": "once", "status": "queued"}})
    payload = {"project_id": api_data["project"]["id"], "schema_version": 1, "workflow": "analysis.ipsae", "parameters": {}}
    key = {"Idempotency-Key": "logical-submission"}
    first = client.post("/api/v1/workflows", json=payload, headers=key)
    second = client.post("/api/v1/workflows", json=payload, headers=key)
    assert first.status_code == second.status_code == 202, first.text
    assert first.json() == second.json() and len(calls) == 1
    assert client.post("/api/v1/workflows", json={**payload, "parameters": {"changed": True}}, headers=key).status_code == 409
    status = client.get("/api/v1/submissions/status", params={"kind": "workflow", "key": "logical-submission"})
    assert status.json()["status"] == "completed" and status.json()["response"] == first.json()
    assert client.get("/api/v1/submissions/status", params={"kind": "workflow", "key": "logical-submission"}, headers=headers(api_data, "student")).status_code == 404


def test_artifact_streaming_stays_inside_the_authorized_run(api_data, monkeypatch):
    from mn_molecular_design_studio import api
    root = api_data["root"] / "artifact-job"
    root.mkdir()
    (root / "scores.txt").write_text("0123456789")
    outside = api_data["root"] / "private.txt"
    outside.write_text("private")
    (root / "escape").symlink_to(outside)
    monkeypatch.setattr(api, "_resolve_api_job", lambda reference, user: {"run_dir": str(root)})
    client = api_data["client"]
    path = "/api/v1/jobs/run/artifact"
    response = client.get(path, params={"path": "scores.txt"}, headers={"Range": "bytes=2-4"})
    assert response.status_code == 206 and response.content == b"234"
    for name in ("../private.txt", str(outside), "escape"):
        assert client.get(path, params={"path": name}).status_code == 404


def test_file_cleanup_cannot_remove_another_users_file(api_data):
    data, client = api_data, api_data["client"]
    record = client.post("/api/v1/files?filename=input.smi", content=b"CCO ethanol").json()
    assert client.delete(record["download_url"], headers=headers(data, "other")).status_code == 404
    assert client.delete(record["download_url"]).status_code == 200
    assert client.get(record["download_url"]).status_code == 404


def test_uncertain_local_submission_is_not_automatically_repeated(api_data, monkeypatch):
    from mn_molecular_design_studio import api
    client, calls = api_data["client"], []
    monkeypatch.setattr(api, "_cluster_submission", lambda *a, **k: None)
    monkeypatch.setattr(api, "_guard_source_job", lambda *a, **k: None)

    def failed(path):
        calls.append(path)
        raise OSError("Lost the scheduler reply")

    monkeypatch.setattr(api, "submit_workflow_request", failed)
    payload = {"project_id": api_data["project"]["id"], "schema_version": 1, "workflow": "analysis.ipsae", "parameters": {}}
    key = {"Idempotency-Key": "uncertain-request"}
    assert client.post("/api/v1/workflows", json=payload, headers=key).status_code == 422
    assert client.post("/api/v1/workflows", json=payload, headers=key).status_code == 409
    assert len(calls) == 1


@pytest.mark.parametrize("recorded", [True, False])
def test_comparison_versions_gate_preview_and_save(api_data, recorded):
    from mn_molecular_design_studio.molecular_database import MolecularDatabase

    client, user = api_data["client"], api_data["admin"]
    database = MolecularDatabase()
    imported = database.import_records(user, [
        {"name": "Fixture A", "components": [{"type": "small_molecule", "chain": "L", "smiles": "CCO"}]},
        {"name": "Fixture B", "components": [{"type": "small_molecule", "chain": "L", "smiles": "CCN"}]},
    ])
    group_id, molecule_ids = imported["group_id"], imported["molecule_ids"]
    digest = "sha256:" + "d" * 64
    for index, molecule_id in enumerate(molecule_ids):
        database.add_observation(group_id, molecule_id, "fixture", f"run-{index}", {"measurement": index + 1},
            context={"target_id": "fixture-target", "target_structure_digest": "a" * 64, "target_chains": ["A"],
                     "settings": {"fixture_setting": 3}, "image_digest": digest if recorded or index == 0 else ""})
    original = database.evaluation_rows(user, molecule_ids=molecule_ids)
    selection = {"molecule_ids": molecule_ids}
    settings = {"fields": ["fixture.measurement"], "directions": {"fixture.measurement": "higher"}}
    endpoint = "/api/v1/molecules/comparisons"
    preview = client.post(endpoint + "/preview", json={"selection": selection, "settings": settings})
    saved = client.post(endpoint, json={**selection, "settings": settings, "name": "Evaluation fixture"})
    if recorded:
        assert preview.status_code == 200, preview.text
        assert saved.status_code == 201, saved.text
        assert [row["Comparison rank"] for row in preview.json()["items"]] == [1, 2]
    else:
        assert preview.status_code == saved.status_code == 422
        assert preview.json()["detail"] == saved.json()["detail"]
        assert "unrecorded execution versions" in preview.json()["detail"]
    descriptive = client.post(endpoint + "/preview", json={"selection": selection})
    assert descriptive.status_code == 200, descriptive.text
    assert all(row["Comparison rank"] is None for row in descriptive.json()["items"])
    assert {row["fixture.measurement"] for row in descriptive.json()["items"]} == {1., 2.}
    assert database.evaluation_rows(user, molecule_ids=molecule_ids) == original


def test_saved_pipeline_revisions_are_accessible_only_in_the_original_scope(api_data):
    from mn_molecular_design_studio.pipeline_graph import new_graph
    data, client = api_data, api_data["client"]
    graph = new_graph()
    saved = client.post("/api/v1/pipelines/save", json={"project_id": data["project"]["id"], "graph": graph, "version": True})
    assert saved.status_code == 200, saved.text
    path = f"/api/v1/pipelines/{graph['id']}/versions"
    query = {"project_id": data["project"]["id"]}
    versions = client.get(path, params=query).json()["items"]
    assert len(versions) == 1
    version_path = path + "/" + versions[0]["revision"]
    assert client.get(version_path, params=query).json()["graph"]["id"] == graph["id"]
    assert client.get(version_path, params=query, headers=headers(data, "other")).status_code == 404
    assert client.get(version_path, params=query, headers=headers(data, "student")).status_code == 404


def test_pipeline_validates_a_database_selection_in_its_own_project(api_data):
    from mn_molecular_design_studio import molecule_imports, workspace
    from mn_molecular_design_studio.capability_pipeline import build_capability_graph
    from mn_molecular_design_studio.database_selections import freeze_protein_selection
    from mn_molecular_design_studio.molecular_database import MolecularDatabase

    data = api_data
    project = data["project"]["id"]
    student = data["student"]
    database = MolecularDatabase()
    imported = database.import_records(student, molecule_imports.parse_sequences("GAGA"), project_id=project)
    source = freeze_protein_selection(database, student, imported["molecule_ids"], project)
    graph = build_capability_graph("sequence_qc", source="molecule_source")
    node = next(node for node in graph["nodes"] if node["data"]["kind"] == "molecule_source")
    node["data"]["params"]["source_run_dir"] = str(source)
    headers = {"Authorization": "Bearer " + data["student_token"]}
    request = {"graph": graph, "project_id": project}

    response = data["client"].post("/api/v1/pipelines/validate", json=request, headers=headers)
    assert response.status_code == 200, response.text
    assert response.json()["ready_to_run"] is True

    other_project = workspace.create_project("Other API course", owner_id=data["admin"]["id"])["id"]
    response = data["client"].post("/api/v1/pipelines/validate", json={**request, "project_id": other_project})
    assert response.status_code == 422
    assert "must belong to the selected project" in response.json()["detail"]


def test_pipeline_save_requires_the_opened_revision(api_data):
    from mn_molecular_design_studio.pipeline_graph import new_graph

    client = api_data["client"]
    project_id = api_data["project"]["id"]
    path = "/api/v1/pipelines/save"
    payload = {"project_id": project_id, "graph": new_graph()}
    opened = client.post(path, json=payload)
    assert opened.status_code == 200, opened.text

    baseline = opened.json()["revision"]
    changed = dict(payload, graph=dict(payload["graph"], name="Newer edit"), expected_revision=baseline)
    result = client.post(path, json=changed)
    assert result.status_code == 200, result.text
    assert client.post(path, json=changed).json() == result.json()

    stale = dict(changed, graph=dict(payload["graph"], name="Stale edit"))
    assert client.post(path, json=stale).status_code == 409
    stale.pop("expected_revision")
    assert client.post(path, json=stale).status_code == 409
    assert client.post(path, json=dict(stale, expected_revision="invalid")).status_code == 422
    [draft] = client.get("/api/v1/pipelines", params={"project_id": project_id}).json()
    assert draft["revision"] == result.json()["revision"]
    assert draft["graph"]["name"] == "Newer edit"


def test_md_starter_and_remote_validation_report_deferred_workstation_checks(api_data, monkeypatch):
    from mn_molecular_design_studio.auth import user_input_root
    from mn_molecular_design_studio.cluster import config

    client = api_data["client"]
    project = api_data["project"]["id"]
    prefix = "/api/v1/pipelines"
    assert "md_analysis" in {item["id"] for item in client.get(prefix + "/templates").json()["items"]}
    draft = client.post(prefix + "/drafts", json={"template": "md_analysis", "project_id": project})
    assert draft.status_code == 200, draft.text
    graph = draft.json()["graph"]
    assert graph["template_id"] == "md_analysis"
    assert [node["data"]["kind"] for node in graph["nodes"]] == ["md_source", "md_analysis", "results"]
    payload = {"graph": graph, "project_id": project}
    local = client.post(prefix + "/validate", json=payload)
    assert local.status_code == 422
    config.save({"mode": "peer", "default_placement": "auto", "coordinator_url": "http://127.0.0.1:18766"})
    remote = client.post(prefix + "/validate", json=payload)
    assert remote.status_code == 422 and local.json() == remote.json()
    graph["nodes"][0]["data"]["params"].update(
        source_run_dir=config.artifact_uri(config.identity(), "md", "saved-fixture"),
        option_id="segment-1", topology_relative_path="artifacts/topology.pdb")
    ready = client.post(prefix + "/validate", json=payload)
    assert ready.status_code == 200, ready.text
    assert ready.json()["readiness"] == {"graph": "valid", "input_selectors": "valid", "source_files": "deferred", "workstation": "deferred"}
    assert not ready.json()["ready_to_run"]
    collision = user_input_root(api_data["admin"]["id"]).parent / "requests" / "artifacts" / "topology.pdb"
    collision.parent.mkdir(parents=True)
    collision.write_text("Not the topology of the saved run.")
    from mn_molecular_design_studio.cluster import transport
    def submit(path, envelope, **kwargs):
        assert path == "/api/v1/cluster/submissions"
        assert envelope["bundle_id"] is None
        selectors = envelope["payload"]["graph"]["nodes"][0]["data"]["params"]
        assert selectors["topology_relative_path"] == "artifacts/topology.pdb"
        assert selectors["source_run_dir"].startswith("cluster://")
        return {"job": {"cluster_job_id": "md-fixture"}}

    monkeypatch.setattr(transport, "request", submit)
    queued = client.post(prefix + "/run", json=payload)
    assert queued.status_code == 202, queued.text
    assert queued.json()["job"]["cluster_job_id"] == "md-fixture"


def test_rejected_pipeline_replays_original_error_but_lost_dispatch_stays_uncertain(api_data, monkeypatch):
    from mn_molecular_design_studio import api
    from mn_molecular_design_studio.cluster import config, service
    from mn_molecular_design_studio.pipeline_graph import new_graph

    client = api_data["client"]
    payload = {"project_id": api_data["project"]["id"], "graph": new_graph("md_analysis")}
    address = "/api/v1/pipelines/run"
    for mode in ("disabled", "peer"):
        config.save({"mode": mode, "default_placement": "local" if mode == "disabled" else "auto",
                     "coordinator_url": "http://127.0.0.1:18766"})
        key = {"Idempotency-Key": "missing-source-" + mode}
        first = client.post(address, json=payload, headers=key)
        replay = client.post(address, json=payload, headers=key)
        assert first.status_code == replay.status_code == 422
        assert first.headers["X-Studio-Submission-State"] == replay.headers["X-Studio-Submission-State"] == "rejected"
        assert first.json() == replay.json()
        status = client.get("/api/v1/submissions/status", params={"kind": "pipeline", "key": key["Idempotency-Key"]}).json()
        assert status["status"] == "rejected" and status["response"]["detail"] == first.json()["detail"]

    calls = []
    def lost_dispatch(*args, on_dispatch, **kwargs):
        on_dispatch()
        calls.append(True)
        raise OSError("Scheduler reply lost")

    monkeypatch.setattr(service, "submit_from_api", lost_dispatch)
    key = {"Idempotency-Key": "lost-reply"}
    lost = client.post(address, json=payload, headers=key)
    assert lost.status_code == 422 and "X-Studio-Submission-State" not in lost.headers
    assert client.post(address, json=payload, headers=key).status_code == 409
    assert calls == [True]
    status = client.get("/api/v1/submissions/status", params={"kind": "pipeline", "key": "lost-reply"}).json()
    assert status["status"] == "uncertain"


def test_archived_pipeline_history_keeps_its_owner_scope_and_rejects_writes(api_data):
    from mn_molecular_design_studio.pipeline_graph import new_graph
    from mn_molecular_design_studio.project_lifecycle import archive_project

    client = api_data["client"]
    project = api_data["project"]["id"]
    graph = new_graph("blank")
    payload = {"project_id": project, "graph": graph, "version": True}
    saved = client.post("/api/v1/pipelines/save", json=payload)
    assert saved.status_code == 200, saved.text
    archive_project(api_data["admin"], project, True)
    query = {"project_id": project}
    versions = f"/api/v1/pipelines/{graph['id']}/versions"
    assert client.get("/api/v1/pipelines", params=query).json()[0]["revision"] == saved.json()["revision"]
    assert len(client.get(versions, params=query).json()["items"]) == 1
    assert client.get(versions + "/" + saved.json()["revision"], params=query).json() == saved.json()
    assert client.get(versions, params=query, headers=headers(api_data, "other")).status_code == 404
    assert client.get("/api/v1/pipelines", params=query, headers=headers(api_data, "student")).json() == []
    assert client.post("/api/v1/pipelines/save", json=payload).status_code == 422
    assert client.post("/api/v1/pipelines/run", json={"project_id": project, "graph": graph}).status_code == 422


def test_run_results_requires_explicit_inspection_and_retries_failed_sources(api_data, monkeypatch):
    from types import SimpleNamespace
    from mn_molecular_design_studio import logical_runs, run_results
    from mn_molecular_design_studio.cluster import service, transport

    client = api_data["client"]
    query = {"project_id": api_data["project"]["id"]}
    source = {"status": "completed", "node_id": "fixture-node", "cluster_job_id": "fixture-job",
              "task_group": "analysis", "run_id": "fixture-run", "local_reference": "analysis/fixture-run",
              "updated_at": "fixture", "availability": "online"}
    run = {"id": "run:fixture", "outputs": [{"id": "results", "reference": "run:fixture", "sources": [source]}]}
    monkeypatch.setattr(logical_runs, "catalog", lambda *args, **kwargs: [run])
    monkeypatch.setattr(service, "lookup_job", lambda *args, **kwargs: source)
    clock = [100.0]
    requests = []
    monkeypatch.setattr(run_results, "time", SimpleNamespace(time=lambda: clock[0]))

    def request(path, payload=None, *, method="GET", token="", **kwargs):
        if method == "POST":
            requests.append(path)
            return {"command_id": f"inspection-{len(requests)}"}
        if path.endswith("inspection-1"):
            return {"status": "failed", "error": "Temporary fixture failure"}
        assert path.endswith("inspection-2")
        return {"status": "completed", "response": {"results": {"candidates": [{"candidate_id": "saved"}]}}}

    monkeypatch.setattr(transport, "request", request)
    address = "/api/v1/runs/run:fixture/results"
    initial = client.get(address, params=query).json()
    assert initial["inspection_required"]
    assert client.get(address, params={**query, "export": "csv"}).status_code == 409
    assert requests == []
    assert client.post(address, params=query, headers=headers(api_data, "other")).status_code == 403
    assert client.post(initial["inspection_url"]).json()["source_errors"][0]["retryable"]
    assert client.get(address, params=query).json()["source_errors"][0]["error"] == "Temporary fixture failure"
    assert client.post(address, params=query).json()["source_errors"][0]["retry_after_seconds"] == 2
    assert len(requests) == 1
    clock[0] += 3
    assert client.post(address, params=query).json()["complete"]
    assert client.get(address, params={**query, "export": "jsonl"}).json()["candidate_id"] == "saved"
    assert len(requests) == 2


def test_bound_pose_api_prepares_mapped_md_files_without_merging_molecules(api_data):
    from rdkit import Chem
    from rdkit.Geometry import Point3D
    from mn_molecular_design_studio import bound_docking
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    data, client = api_data, api_data["client"]
    receptor = data["root"] / "receptor.pdb"
    receptor.write_text("ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C\nEND\n")
    database = MolecularDatabase()
    target = bound_docking.register_receptor(database, data["admin"], receptor, project_id=data["project"]["id"])
    imported = database.import_records(data["admin"], [{"name": "Pose", "components": [{"type": "small_molecule", "chain": "L", "smiles": "CCO"}]}], project_id=data["project"]["id"])
    molecule = Chem.MolFromSmiles("CCO")
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    for i in range(molecule.GetNumAtoms()):
        conformer.SetAtomPosition(i, Point3D(5 + i, 2, 3))
    molecule.AddConformer(conformer)
    pose = data["root"] / "pose.sdf"
    with Chem.SDWriter(str(pose)) as writer:
        writer.write(molecule)
    identifier = imported["molecule_ids"][0]
    structure = database.attach_structure(imported["group_id"], identifier, {"path": str(pose), "label": "Pose"},
        context={**target, "database_structure_policy": "best_docking_pose"})
    before = database.list_molecules(data["admin"])["total"]
    prepared = client.post(f"/api/v1/molecules/{identifier}/bound-inputs", json={"structure_id": structure})
    assert prepared.status_code == 201, prepared.text
    result = prepared.json()
    assert result["workflow_inputs"]["complex_pdb"].startswith("studio-file:///")
    assert client.get(result["files"]["complex_pdb"]["download_url"]).content.count(b"HETATM") == 3
    assert client.get(result["files"]["ligand_sdf"]["download_url"]).status_code == 200
    assert database.list_molecules(data["admin"])["total"] == before == 2
    shared = client.post(f"/api/v1/molecules/{identifier}/bound-inputs", json={}, headers=headers(data, "student"))
    assert shared.status_code == 201, shared.text
    assert client.get(shared.json()["files"]["ligand_sdf"]["download_url"], headers=headers(data, "student")).status_code == 200
    assert client.post(f"/api/v1/molecules/{identifier}/bound-inputs", json={}, headers=headers(data, "other")).status_code == 404


def test_mcp_uses_rest_and_preserves_user_permissions(api_data):
    pytest.importorskip("mcp")
    from mcp import Client
    from mn_molecular_design_studio.api_client import StudioClient
    from mn_molecular_design_studio.mcp_server import create_server
    studio = StudioClient("http://studio.test", token=api_data["student_token"], transport=httpx.ASGITransport(app=api_data["app"]))

    async def check():
        async with Client(create_server(studio)) as client:
            tools = (await client.list_tools()).tools
            assert {"workflow_schema", "screening_results", "prepare_bound_pose", "run_pipeline"}.issubset({t.name for t in tools})
            schema = await client.call_tool("workflow_schema", {"workflow_id": "ligand.docking.unidock_pro"})
            assert not schema.is_error, schema.content
            assert "parameters_json_schema" in schema.structured_content
            identity = await client.call_tool("studio_identity", {})
            assert identity.structured_content["role"] == "user"
            rejected = await client.call_tool("studio_api_request", {"method": "GET", "path": "/api/v1/admin/users"})
            assert rejected.is_error
            blocked = await client.call_tool("studio_api_request", {"method": "POST", "path": "/api/v1/cluster/node/claim"})
            assert blocked.is_error
            staged = await client.call_tool("stage_text_input", {"filename": "ligands.smi", "text": "CCO ethanol\n"})
            assert staged.structured_content["reference"].startswith("studio-file:///")
    asyncio.run(check())


def test_mcp_http_never_falls_back_to_admin_or_development_identity(api_data, monkeypatch):
    pytest.importorskip("mcp")
    from mn_molecular_design_studio.api_client import StudioClient
    from mn_molecular_design_studio.mcp_server import http_app
    studio = StudioClient("http://studio.test", token=api_data["admin_token"], transport=httpx.ASGITransport(app=api_data["app"]))
    with TestClient(http_app(client=studio)) as client:
        assert client.post("/mcp", json={}).status_code == 401
        assert client.post("/mcp", json={}, headers={"Authorization": "Bearer invalid"}).status_code == 401
        monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_AUTH_MODE", "development")
        assert client.post("/mcp", json={}, headers=headers(api_data, "student")).status_code == 401


def test_tool_connection_map_through_rest_and_mcp(api_data):
    from mcp import Client
    from mn_molecular_design_studio.api_client import StudioClient
    from mn_molecular_design_studio.mcp_server import create_server
    from mn_molecular_design_studio.pipeline_graph import connection_compatibility

    data = api_data
    path = "/api/v1/workflows/connections"
    client = data["client"]
    assert client.get(path, headers={"Authorization": ""}).status_code == 401
    assert client.get(path, params={"limit": 101}, headers=headers(data, "student")).status_code == 422

    def get(**params):
        response = client.get(path, params=params, headers=headers(data, "student"))
        assert response.status_code == 200, response.text
        return response.json()

    tiles = get(view="tiles", limit=100)
    by_id = {tile["id"]: tile for tile in tiles["items"]}
    assert set(by_id["optimized-image-boltz2"]["pipeline_operations"]) == {"boltz2"}
    assert set(by_id["group-boltzgen_fold"]["pipeline_operations"]) == {"boltzgen_fold_batch", "prediction_boltzgen_fold"}
    assert by_id["refolding.af2_initial_guess"]["modes"][0]["workflow"] == "refolding.af2_initial_guess_input"
    assert by_id["studio.validation-acceptance"]["modes"][0]["workflow"] == "analysis.validation_acceptance"
    assert get(view="tiles", offset=1, limit=1)["items"] == tiles["items"][1:2]

    workflows = get(view="workflows", query="library.prepare")["items"]
    assert {item["name"] for item in workflows[0]["parameters"] if item["required"]} >= {"source_files", "library_name"}
    operation = get(view="operations", query="md_analysis")["items"][0]
    assert operation["inputs"][0]["artifact_types"] == ["md_trajectory"]
    edges = get(view="connections", operation="md_analysis", direction="incoming", artifact_type="md_trajectory")
    assert edges["items"]
    for edge in edges["items"]:
        compatibility = connection_compatibility(edge["source"], edge["target"], source_handle=edge["output"], target_handle=edge["input"])
        assert edge["target"] == "md_analysis" and compatibility["artifact_type"] == edge["artifact_type"]

    studio = StudioClient("http://studio.test", token=data["student_token"], transport=httpx.ASGITransport(app=data["app"]))

    async def check():
        async with Client(create_server(studio)) as mcp:
            tools = (await mcp.list_tools()).tools
            assert "tool_connections" in {tool.name for tool in tools}
            response = await mcp.call_tool("tool_connections", {"view": "handoffs", "query": "ligand.docking.vina"})
            assert not response.is_error, response.content
            handoffs = response.structured_content["items"]
            assert {row["target"] for row in handoffs if row["source"] == "ligand.docking.vina"} == {
                "analysis.native_ligand_interactions", "analysis.plip", "analysis.pandamap", "analysis.posebusters",
            }
            assert any(row["source"] == "ligand.generation.lddm" and row["target"] == "ligand.docking.vina"
                       for row in handoffs)
            assert all(row["transfer"] == "source_job" for row in handoffs)
            sockets = await mcp.call_tool("tool_connections", {
                "view": "connections", "operation": "filter", "direction": "outgoing", "port": "fail",
            })
            assert not sockets.is_error, sockets.content
            assert sockets.structured_content["items"]
            assert all(row["source"] == "filter" and row["output"] == "fail"
                       for row in sockets.structured_content["items"])

    asyncio.run(check())


def test_mcp_operates_projects_pipelines_datasets_and_saved_results(api_data):
    from mcp import Client
    from mn_molecular_design_studio.api_client import StudioClient
    from mn_molecular_design_studio.mcp_server import create_server
    from mn_molecular_design_studio.molecular_database import MolecularDatabase

    data = api_data
    studio = StudioClient("http://studio.test", token=data["student_token"], transport=httpx.ASGITransport(app=data["app"]))

    async def check():
        async with Client(create_server(studio)) as client:
            async def tool(name, **arguments):
                result = await client.call_tool(name, arguments)
                assert not result.is_error, result.content
                return result.structured_content

            project = await tool("studio_api_request", method="POST", path="/api/v1/projects", body={"name": "MCP practice project"})
            project_id = project["id"]
            assert project_id in {row["id"] for row in (await tool("list_projects"))["items"]}
            operations = await tool("studio_operations")
            assert any(item["path"] == "/api/v1/molecules/browser/results" for item in operations["operations"])
            contract = await tool("studio_operation_schema", method="GET", path="/api/v1/molecules/browser/results")
            assert any(item["name"] == "folder" for item in contract["operation"]["parameters"])

            starter = await tool("create_pipeline_draft", template="md_analysis", project_id=project_id)
            draft = await tool("save_pipeline_draft", project_id=project_id, graph=starter["graph"], version=True)
            history = await tool("list_saved_pipelines", project_id=project_id, pipeline_id=draft["graph"]["id"])
            assert history["items"][0]["revision"] == draft["revision"]

            imported = await tool("studio_api_request", method="POST", path="/api/v1/molecules/imports", body={
                "name": "MCP inputs", "method": "smiles", "text": "CCO ethanol", "project_references": [project_id], "reusable_dataset": True})
            inputs = await tool("list_datasets", project_id=project_id)
            assert [row["id"] for row in inputs["items"]] == [imported["group_id"]]
            one = await tool("dataset_members", group_id=imported["group_id"])
            assert one["items"][0]["name"] == "ethanol" and not one["items"][0]["observations"]

            database = MolecularDatabase()
            saved = database.import_records(data["student"], [{"name": "Saved fixture", "components": [
                {"type": "small_molecule", "chain": "L", "smiles": "CCO"}]}], project_id=project_id,
                name="Saved docking fixture", origin={"collection_kind": "curated", "workflow_role": "docking"})
            structure = database.attach_structure(saved["group_id"], saved["molecule_ids"][0],
                {"name": "pose.pdb", "content": b"HETATM    1  C1  LIG L   1       1.000   2.000   3.000  1.00  0.00           C\nEND\n"},
                context={"run_id": "fixture-run", "candidate_id": "fixture-pose"})
            database.add_observation(saved["group_id"], saved["molecule_ids"][0], "vina", "fixture-run", {"energy": -5.},
                context={"workflow_role": "docking", "task_group": "ligand-docking", "calculation_recorded": True,
                         "representation_id": structure, "candidate_id": "fixture-pose"})
            roots = await tool("browse_saved_results")
            assert "projects" in {row["id"] for row in roots["folders"]}
            results = await tool("browse_saved_results", folder="project:" + project_id, limit=1)
            assert results["items"][0]["row"]["Result · energy"] == -5.
            assert results["items"][0]["molecule_id"] == saved["molecule_ids"][0]
            exported = await tool("studio_api_request", method="GET", path="/api/v1/molecules/browser/results",
                                  query={"folder": "project:" + project_id, "export": "csv"})
            assert exported["requires_authentication"] and "csv" in exported["media_type"]
            await tool("studio_api_request", method="PATCH", path=f"/api/v1/projects/{project_id}/archive", body={"archived": True})
            archived = await tool("list_saved_pipelines", project_id=project_id)
            assert archived["items"][0]["revision"] == draft["revision"]

            unexpected = await client.call_tool("studio_api_request", {"method": "POST",
                "path": f"/api/v1/admin/projects/{project_id}/cleanup/reviews", "body": {"mode": "delete"}})
            assert unexpected.is_error
            denied = data["client"].get("/api/v1/molecules/browser/results", params={"folder": "project:" + project_id},
                                        headers=headers(data, "other"))
            assert denied.status_code == 404
    asyncio.run(check())


def test_project_cleanup_review_and_status_use_existing_guarded_service(api_data):
    client = api_data["client"]
    project = api_data["project"]["id"]
    response = client.post(f"/api/v1/admin/projects/{project}/cleanup/reviews", json={"mode": "delete"})
    assert response.status_code == 201, response.text
    plan = response.json()
    assert plan["status"] == "review" and plan["project_id"] == project
    assert client.get(f"/api/v1/admin/project-cleanup/{plan['id']}").json() == plan
    assert client.get(f"/api/v1/projects/{project}").status_code == 200


def test_mcp_http_forwards_each_users_identity_without_shared_credentials(api_data):
    pytest.importorskip("mcp")
    from mn_molecular_design_studio.api_client import StudioClient
    from mn_molecular_design_studio.mcp_server import http_app
    studio = StudioClient("http://studio.test", token=api_data["admin_token"], transport=httpx.ASGITransport(app=api_data["app"]))
    with TestClient(http_app(client=studio), base_url="http://localhost") as client:
        for role in ("student", "admin"):
            request_headers = {**headers(api_data, role), "Accept": "application/json, text/event-stream", "MCP-Protocol-Version": "2025-11-25"}
            response = client.post("/mcp", headers=request_headers, json={"jsonrpc": "2.0", "id": 1,
                "method": "tools/call", "params": {"name": "studio_identity", "arguments": {}}})
            assert response.status_code == 200, response.text
            messages = [json.loads(line[6:]) for line in response.text.splitlines() if line.startswith("data: ")]
            assert messages, response.text
            result = messages[-1]["result"]
            assert not result.get("isError"), result
            assert result["structuredContent"]["role"] == ("user" if role == "student" else "admin")


def test_workflow_json_schemas_are_valid_and_never_evaluate_type_text():
    from jsonschema import Draft202012Validator
    from mn_molecular_design_studio.workflow_json_schema import type_schema, parameters_schema
    assert type_schema("list[pathlib.Path] | None") == {"anyOf": [{"type": "array", "items": {"type": "string"}}, {"type": "null"}]}
    assert type_schema("__import__('os').system('must-not-run')") == {}
    assert type_schema("one of 'fast', 'slow'") == {"enum": ["fast", "slow"]}
    schema = parameters_schema([{"name": "counts", "type": "tuple[int, int]", "required": True}])
    Draft202012Validator.check_schema(schema)
    Draft202012Validator(schema).validate({"counts": [1, 2]})
