"""Exercise real API/queue/worker bookkeeping without launching prediction containers."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from fastapi.testclient import TestClient

from mn_molecular_design_studio import job_catalog, workspace
from mn_molecular_design_studio.engine.workflow_engine.core import (
    central_scheduler, jobs, local_worker, workflow_queue,
)
from mn_molecular_design_studio.engine.workflow_engine.workflows import optimized_prediction as prediction


PDB = "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C  \nEND\n"


# Queue tests still simulate execution, but their result must honor the requested sequence.
from Bio.SeqUtils import seq3
PREDICTION_PDB = "".join(
    f"ATOM  {i:5d}  CA  {seq3(aa).upper()} A{i:4d}    {i * 3.8:8.3f}{0.:8.3f}{0.:8.3f}  1.00 90.00           C  \n"
    for i, aa in enumerate("ACDEFGHIK", 1)
) + "END\n"


@pytest.fixture
def api_queue(monkeypatch, tmp_path):
    # Keep API initialization, catalog writes, and queued runs inside this test.
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    monkeypatch.setattr(workspace, "ensure_studio_workspace", lambda: tmp_path)
    from mn_molecular_design_studio import api

    root = tmp_path / "runs"
    root.mkdir()
    dispatched = []
    monkeypatch.setattr(api, "ensure_studio_workspace", lambda: tmp_path)
    monkeypatch.setattr(api, "validate_project_id", lambda value: value)
    monkeypatch.setattr(api, "require_project_access", lambda *args: None)
    monkeypatch.setattr(api, "accessible_project_ids", lambda user: [])
    monkeypatch.setattr(api, "_request_directory", lambda *args: tmp_path)
    monkeypatch.setattr(api, "mark_studio_job", lambda result, **kwargs: result)
    monkeypatch.setattr(jobs, "runs_root", lambda: root)
    monkeypatch.setattr(job_catalog, "_results_root", lambda: root)
    monkeypatch.setattr(job_catalog, "catalog_path", lambda: tmp_path / "catalog.sqlite3")
    monkeypatch.setattr(job_catalog, "_ROOT_SCAN_INTERVAL_SECONDS", 0)
    monkeypatch.setattr(prediction, "reference_root", lambda: tmp_path / "reference")
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_CACHE_DIR", str(tmp_path / "cache"))

    def dispatch(run_dir):
        # Capture identity at dispatch time, before a worker can repair it.
        dispatched.append(jobs.read_json(run_dir / "metadata.json"))

    monkeypatch.setattr(workflow_queue, "spawn_worker_for_run", dispatch)
    api.app.dependency_overrides[api._current_user] = lambda: {"id": "test-user", "role": "admin"}
    try:
        with TestClient(api.app) as client:
            yield SimpleNamespace(client=client, root=root, dispatched=dispatched)
    finally:
        api.app.dependency_overrides.pop(api._current_user, None)


def request_for(engine="protenix_v2"):
    profile = prediction.ENGINE_PROFILES[engine]
    task = {"protein_sequences": "ACDEFGHIK"}
    if profile["input_kind"] == "boltzgen":
        task = {"target_pdb_text": PDB, "binder_sequence": "ACDEFGHIK"}
    return {
        "schema_version": 1,
        "workflow": "refolding.optimized_prediction",
        "project_id": "test-project",
        "parameters": {
            "engine": engine, "optimization_mode": profile["modes"][0],
            "job_name": f"test-{engine}", "seed": 101, "samples": 1,
            "gpu_device": "0", "accept_alphafold3_terms": bool(profile.get("terms")),
            "tasks": [task],
        },
    }


def test_mpnn_api_and_worker_reuse_one_named_run(api_queue, monkeypatch):
    from mn_molecular_design_studio import mpnn_input
    from mn_molecular_design_studio.engine.workflow_engine.workflows import sequence_design
    structure = api_queue.root.parent / "input.pdb"
    structure.write_text(PDB)
    weights = api_queue.root.parent / "hyper.pt"
    weights.write_bytes(b"test checkpoint")
    monkeypatch.setattr(mpnn_input, "checkpoint", lambda *a, **k: weights)
    payload = {"schema_version": 1, "workflow": "sequence.mpnn", "project_id": "test-project", "parameters": {
        "structure_path": str(structure), "design_chains": "A", "model_type": "hyper_mpnn", "job_name": "custom-hyper-design"}}
    response = api_queue.client.post("/api/v1/workflows", json=payload)
    assert response.status_code == 202, response.text
    run = Path(response.json()["job"]["run_dir"])
    assert response.json()["job"]["campaign_name"] == "custom-hyper-design"
    assert response.json()["job"]["tool"] == "ligandmpnn"
    def fake_native(*args, existing_job, **kwargs):
        assert existing_job.run_dir == run
        jobs.write_json(run / "result.json", {"success": True})
        jobs.update_status(run, "completed")
        return run
    monkeypatch.setattr(sequence_design, "run_ligandmpnn_sequence_design", fake_native)
    assert local_worker.run_worker_job(run) == 0
    assert list(api_queue.root.glob("*/*/metadata.json")) == [run / "metadata.json"]


@pytest.mark.parametrize("engine", prediction.ENGINE_PROFILES)
def test_api_creates_one_engine_specific_named_job_before_dispatch(api_queue, engine):
    response = api_queue.client.post("/api/v1/workflows", json=request_for(engine))
    assert response.status_code == 202, response.text
    run = Path(response.json()["job"]["run_dir"])
    metadata = jobs.read_json(run / "metadata.json")
    inputs = jobs.read_json(run / "input.json")
    expected = f"optimized_prediction_{engine}"
    assert response.json()["job"]["tool"] == expected
    assert response.json()["job"]["campaign_name"] == f"test-{engine}"
    assert metadata["tool"] == inputs["tool"] == expected
    assert metadata["campaign_name"] == f"test-{engine}"
    assert metadata["status"] == "queued"
    assert api_queue.dispatched[0]["tool"] == expected
    assert api_queue.dispatched[0]["campaign_name"] == f"test-{engine}"
    assert list(api_queue.root.glob("*/*/metadata.json")) == [run / "metadata.json"]
    # The API's SQLite listing and engine filter must see the same queued job.
    listed = api_queue.client.get("/api/v1/jobs", params={"status": "queued", "tool": expected}).json()
    assert len(listed["items"]) == 1
    assert listed["items"][0]["campaign_name"] == f"test-{engine}"


@pytest.mark.parametrize("engine", prediction.ENGINE_PROFILES)
def test_api_selected_prediction_inputs_automatically_share_one_native_job(api_queue, engine):
    payload = request_for(engine)
    parameters = payload["parameters"]
    payload["parameters"] = [parameters, {**parameters, "job_name": "second selected input"}]
    response = api_queue.client.post("/api/v1/workflows", json=payload)
    assert response.status_code == 202, response.text
    run = Path(response.json()["job"]["run_dir"])
    invocation = jobs.read_json(run / "worker_request.json")
    assert invocation["function"] == "run_optimized_prediction"
    assert len(invocation["kwargs"]["tasks"]) == 2
    assert len(api_queue.dispatched) == 1
    assert len(list(api_queue.root.glob("*/*/metadata.json"))) == 1
    assert jobs.read_json(run / "metadata.json")["batch"]["item_count"] == 2


def test_invalid_selected_input_does_not_partially_submit_batch(api_queue):
    payload = request_for("boltz2")
    parameters = payload["parameters"]
    payload["parameters"] = [parameters, {**parameters, "tasks": [{}]}]
    response = api_queue.client.post("/api/v1/workflows", json=payload)
    assert response.status_code == 422
    assert api_queue.dispatched == []
    assert not list(api_queue.root.glob("*/*/metadata.json"))


@pytest.mark.parametrize("engine", ["boltz2", "alphafold3", "af3_torch", "protenix_v2"])
def test_selected_inputs_reach_container_as_one_native_batch(api_queue, monkeypatch, engine):
    reference = api_queue.root.parent / "reference"
    cache = reference / "boltz_models"
    cache.mkdir(parents=True)
    for name in ("boltz2_conf.ckpt", "boltz2_aff.ckpt", "ccd.pkl", "mols.tar"):
        (cache / name).write_bytes(b"test reference")
    payload = request_for(engine)
    parameters = payload["parameters"]
    payload["parameters"] = [parameters, {**parameters, "job_name": "second selected input"}]
    response = api_queue.client.post("/api/v1/workflows", json=payload)
    assert response.status_code == 202, response.text
    run = Path(response.json()["job"]["run_dir"])
    metadata = jobs.read_json(run / "metadata.json")
    jobs.write_json(run / "metadata.json", {**metadata, "scheduler_backend": central_scheduler.BACKEND})
    invocations = []
    def stage(run_dir, **kwargs):
        path = run_dir / "artifacts" / (kwargs["item_id"] + ".json")
        path.write_text("{}")
        return path, {}
    def predict(command, *, cwd, **kwargs):
        invocations.append(command)
        raw = cwd / "artifacts/raw/optimized_prediction"
        manifest = json.loads((raw / "mn_predict_manifest.json").read_text())
        assert len(manifest["items"]) == 2
        rows = []
        for item in manifest["items"]:
            structure = raw / "output" / (item["id"] + ".pdb")
            structure.write_text(PREDICTION_PDB)
            rows.append({"id": item["id"], "status": "completed", "outputs": [structure.relative_to(cwd).as_posix()]})
        (raw / "output/mn-batch-result.json").write_text(json.dumps({"items": rows}))
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(prediction, "_stage_native_input", stage)
    monkeypatch.setattr(prediction.subprocess, "run", predict)
    workflow_queue.run_workflow_call(run, jobs.read_json(run / "worker_request.json"))
    assert len(invocations) == 1
    assert jobs.read_json(run / "metadata.json")["status"] == "completed"
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    assert len(read_candidates(run)) == 2


@pytest.mark.parametrize("changes", [
    {"tasks": []},
    {"tasks": [{}]},
    {"tasks": [{"name": "empty", "protein_sequences": ""}]},
    {"tasks": ["not an object"]},
    {"samples": 0},
    {"seed": -1},
    {"cpu_cores": 0},
    {"num_recycles": 0},
    {"engine_options": {"not_an_engine_option": True}},
    {"engine": "rf3_foundry", "optimization_mode": "fast"},
    {"engine": "alphafold3", "accept_alphafold3_terms": False},
    {"engine": "esmfold", "tasks": [{"ligand_smiles": "C"}]},
])
def test_invalid_requests_never_create_or_dispatch_a_job(api_queue, changes):
    payload = request_for()
    payload["parameters"].update(changes)
    for endpoint in ["/api/v1/workflows/validate", "/api/v1/workflows"]:
        response = api_queue.client.post(endpoint, json=payload)
        assert response.status_code == 422, response.text
    assert not list(api_queue.root.glob("*/*/metadata.json"))
    assert api_queue.dispatched == []


def test_preflight_does_not_generate_random_seeds_or_stage_inputs(api_queue, monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail("Preflight must not generate seeds, stage inputs, or run predictions")

    monkeypatch.setattr(prediction.secrets, "randbelow", forbidden)
    monkeypatch.setattr(prediction, "_stage_native_input", forbidden)
    payload = request_for()
    payload["parameters"]["seed"] = 0
    response = api_queue.client.post("/api/v1/workflows/validate", json=payload)
    assert response.status_code == 200, response.text
    assert response.json()["parameters"]["seed"] == 0
    assert api_queue.dispatched == []


@pytest.mark.parametrize("source_project,passed,selected_project,error", [
    ("test-project", True, None, "Select a project"),
    ("other-project", True, "test-project", "selected project"),
    ("test-project", False, "test-project", "filter-passing"),
])
def test_design_validation_rejects_missing_project_foreign_campaign_or_rejected_design(
    api_queue, design_source_factory, source_project, passed, selected_project, error,
):
    run, candidate = design_source_factory(api_queue.root / "design/source", project_id=source_project, passed=passed)
    payload = request_for("boltz2")
    payload["project_id"] = selected_project
    payload["parameters"].pop("tasks")
    payload["parameters"].update(validation_sources=[{"source_run_dir": str(run), "candidate_ids": [candidate["candidate_id"]]}],
                                   validation_mode="complex")
    for endpoint in ("/api/v1/workflows/validate", "/api/v1/workflows"):
        response = api_queue.client.post(endpoint, json=payload)
        assert response.status_code == 422, response.text
        assert error in response.text
    assert api_queue.dispatched == []
    assert list(api_queue.root.glob("refolding-validation/*/metadata.json")) == []


@pytest.mark.parametrize("mode", ["complex", "binder"])
def test_selected_passing_design_is_queued_once_in_each_mode(api_queue, design_source_factory, mode):
    run, candidate = design_source_factory(api_queue.root / "design/source")
    payload = request_for("boltz2")
    payload["parameters"].pop("tasks")
    payload["parameters"].update(validation_sources=[{"source_run_dir": str(run), "candidate_ids": [candidate["candidate_id"]]}],
                                   validation_mode=mode)
    response = api_queue.client.post("/api/v1/workflows/validate", json=payload)
    assert response.status_code == 200, response.text
    assert api_queue.dispatched == []
    response = api_queue.client.post("/api/v1/workflows", json=payload)
    assert response.status_code == 202, response.text
    assert response.json()["job"]["tool"] == "optimized_prediction_boltz2"
    assert len(api_queue.dispatched) == 1
    assert len(list(api_queue.root.glob("refolding-validation/*/metadata.json"))) == 1


@pytest.mark.parametrize("fail_staging", [False, True])
def test_worker_adopts_same_job_and_never_returns_to_queued(api_queue, monkeypatch, fail_staging):
    response = api_queue.client.post("/api/v1/workflows", json=request_for())
    assert response.status_code == 202, response.text
    run = Path(response.json()["job"]["run_dir"])
    original = jobs.read_json(run / "metadata.json")
    phases = []
    original.update({"scheduler_backend": central_scheduler.BACKEND, "scheduler_job_id": "test-lease"})
    jobs.write_json(run / "metadata.json", original)
    jobs.update_status(run, "running")

    def stage(run_dir, **kwargs):
        phases.append(jobs.read_json(run_dir / "metadata.json")["status"])
        if fail_staging:
            raise ValueError("Simulated preparation failure")
        native = run_dir / "artifacts" / "native.json"
        native.write_text("{}")
        return native, {}

    def predict(command, *, cwd, **kwargs):
        phases.append(jobs.read_json(cwd / "metadata.json")["status"])
        raw = cwd / "artifacts/raw/optimized_prediction"
        manifest = json.loads((raw / "mn_predict_manifest.json").read_text())
        structure = raw / "output/model.pdb"
        structure.write_text(PREDICTION_PDB)
        (raw / "output/mn-batch-result.json").write_text(json.dumps({
            "items": [{"id": manifest["items"][0]["id"], "status": "completed",
                       "outputs": [structure.relative_to(cwd).as_posix()]}],
        }))
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(prediction, "_stage_native_input", stage)
    monkeypatch.setattr(prediction.subprocess, "run", predict)
    result = workflow_queue.run_workflow_call(run, jobs.read_json(run / "worker_request.json"))
    metadata = jobs.read_json(run / "metadata.json")
    assert result == run
    assert metadata["job_code"] == original["job_code"]
    assert metadata["scheduler_job_id"] == "test-lease"
    assert metadata["tool"] == "optimized_prediction_protenix_v2"
    assert metadata["status"] == ("failed" if fail_staging else "completed")
    assert jobs.read_json(run / "result.json")["success"] is (not fail_staging)
    assert phases == (["preparing"] if fail_staging else ["preparing", "running"])
    assert list(api_queue.root.glob("*/*/metadata.json")) == [run / "metadata.json"]


def test_successful_finish_persists_cleanup_of_previous_worker_error(api_queue):
    run = api_queue.root / "design" / "recovered-generator"
    run.mkdir(parents=True)
    jobs.write_json(run / "input.json", {"job_type": "design_generation", "tool": "pxdesign_generator"})
    jobs.write_json(run / "metadata.json", {
        "status": "failed", "worker_error": "previous parser failure",
        "worker_exception_type": "RuntimeError", "warning": "Degraded run: previous parser failure.",
        "warning_level": "degraded",
    })

    jobs.finish_job(run, True, {"metrics": {"candidate_count": 1}})

    metadata = jobs.read_json(run / "metadata.json")
    assert metadata["status"] == "completed"
    assert not {"worker_error", "worker_exception_type", "warning", "warning_level"}.intersection(metadata)
    assert jobs.read_json(run / "result.json")["success"] is True


@pytest.mark.parametrize("status", ["queued", "running", "preparing", "holding", "completed", "failed", "cancelled"])
def test_legacy_generic_rows_keep_status_and_resolve_engine_and_name(api_queue, status):
    run = api_queue.root / "refolding-validation" / "legacy-job"
    run.mkdir(parents=True)
    (run / "metadata.json").write_text(json.dumps({"tool": "optimized_prediction", "status": status}))
    (run / "input.json").write_text(json.dumps({
        "tool": "optimized_prediction", "params": {"engine": "esmfold", "job_name": "my-old-job"},
    }))
    if status in {"preparing", "holding", "completed", "failed", "cancelled"}:
        (run / "result.json").write_text(json.dumps({"success": status == "completed"}))
    for row in [job_catalog._job_payload(run), jobs.collect_jobs(include_hidden=True)[0]]:
        assert row["tool"] == "optimized_prediction_esmfold"
        assert row["campaign_name"] == "my-old-job"
        assert row["status"] == status
    response = api_queue.client.get("/api/v1/jobs", params={"tool": "optimized_prediction_esmfold"})
    assert response.status_code == 200
    assert len(response.json()["items"]) == 1


@pytest.mark.parametrize("backend", ["central", "local"])
def test_dispatch_backends_keep_one_correctly_labelled_run(api_queue, monkeypatch, backend):
    dispatched = []
    monkeypatch.setattr(workflow_queue, "spawn_worker_for_run", local_worker.spawn_worker_for_run)
    monkeypatch.setattr(central_scheduler, "scheduler_mode", lambda: backend)
    monkeypatch.setattr(central_scheduler, "submit_run", lambda run: dispatched.append(run))
    monkeypatch.setattr(local_worker, "ensure_worker_service", lambda: dispatched.append("local-service"))
    response = api_queue.client.post("/api/v1/workflows", json=request_for("boltz2"))
    assert response.status_code == 202, response.text
    run = Path(response.json()["job"]["run_dir"])
    assert len(dispatched) == 1
    assert jobs.read_json(run / "metadata.json")["tool"] == "optimized_prediction_boltz2"
    assert len(list(api_queue.root.glob("*/*/metadata.json"))) == 1


def test_repeated_editor_snapshot_does_not_repeat_submission(monkeypatch):
    from mn_molecular_design_studio import cofolding_input

    model = cofolding_input.default_model()
    model["cards"][0]["value"] = "ACDEFGHIK"
    model["submit_token"] = "first-click"
    monkeypatch.setattr(cofolding_input.st, "session_state", {})
    monkeypatch.setattr(cofolding_input, "_get_editor", lambda: lambda **kwargs: model)
    assert cofolding_input.render_molecule_editor("test", ())[2] is True
    assert cofolding_input.render_molecule_editor("test", ())[2] is False
    model["submit_token"] = "second-click"
    assert cofolding_input.render_molecule_editor("test", ())[2] is True


@pytest.mark.parametrize("fail_replace", [False, True])
def test_metadata_writes_preserve_valid_previous_file_until_atomic_replace(monkeypatch, tmp_path, fail_replace):
    path = tmp_path / "metadata.json"
    path.write_text('{"status": "queued"}')
    original_replace = Path.replace

    def replace(temporary, target):
        assert json.loads(path.read_text()) == {"status": "queued"}
        assert json.loads(temporary.read_text()) == {"status": "running"}
        if fail_replace:
            raise OSError("Simulated interrupted write")
        return original_replace(temporary, target)

    monkeypatch.setattr(Path, "replace", replace)
    if fail_replace:
        with pytest.raises(OSError, match="Simulated interrupted write"):
            jobs.write_json(path, {"status": "running"})
    else:
        jobs.write_json(path, {"status": "running"})
    assert json.loads(path.read_text())["status"] == ("queued" if fail_replace else "running")
    assert not list(tmp_path.glob("*.tmp"))


def test_rfdiffusion3_api_has_specific_identity_and_preflight(api_queue):
    payload = {"schema_version": 1, "workflow": "design.rfdiffusion3", "project_id": "test-project",
               "parameters": {"task": "json", "job_name": "rfd3-api-contract", "num_designs": 1,
                              "json_config": {"unconditional": {"length": "64-64"}}}}
    assert api_queue.client.post("/api/v1/workflows/validate", json=payload).status_code == 200
    response = api_queue.client.post("/api/v1/workflows", json=payload)
    assert response.status_code == 202, response.text
    run = Path(response.json()["job"]["run_dir"])
    metadata = jobs.read_json(run / "metadata.json")
    assert metadata["tool"] == "rfdiffusion3"
    assert metadata["job_type"] == "rfdiffusion3_design"
    assert metadata["campaign_name"] == "rfd3-api-contract"
    assert len(list(api_queue.root.glob("*/*/metadata.json"))) == 1
    payload["parameters"] = {"task": "protein_binder"}
    assert api_queue.client.post("/api/v1/workflows", json=payload).status_code == 422
    assert len(list(api_queue.root.glob("*/*/metadata.json"))) == 1
