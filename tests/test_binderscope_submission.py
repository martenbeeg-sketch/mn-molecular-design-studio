from __future__ import annotations

import copy
import json
import os
from pathlib import Path

import pytest

from test_binderscope import saved_subject
from test_database_browser import USER, OTHER, database


def cohort(database, tmp_path, count=1):
    imported, saved, selector = saved_subject(database, tmp_path)
    sources = []
    with database.connection() as connection:
        columns = [row[1] for row in connection.execute("PRAGMA table_info(molecule_structures)") if row[1] != "id"]
        selected = ",".join(columns)
        for position in range(count):
            identifier = "cohort-structure-" + str(position)
            connection.execute("INSERT INTO molecule_structures(id," + selected + ") SELECT ?," + selected +
                               " FROM molecule_structures WHERE id=?", (identifier, selector["structure_id"]))
            sources.append({**selector, "structure_id": identifier, "assessment": {
                "component_format": "linear protein", "b_factor_interpretation": "plddt_0_100"}})
    return {"schema_version": 1, "workflow": "analysis.binderscope", "project_id": None,
        "parameters": {"use_recorded_settings": True, "cpu_cores": 32}, "resources": {"gpu_device": "none"},
        "binderscope_collection": {"schema": "studio.binderscope_selection.v1", "sources": sources}}


def queued(database, tmp_path, monkeypatch, count=1, payload=None):
    from mn_molecular_design_studio.workspace import save_settings, ensure_studio_workspace

    save_settings({"results_directory": str(tmp_path / "runs"), "reference_directory": str(tmp_path / "reference")},
                  home=tmp_path / "workspace")
    from mn_molecular_design_studio import api, binderscope_submission
    from mn_molecular_design_studio.engine.workflow_engine.core import workflow_queue

    ensure_studio_workspace()
    payload = payload or cohort(database, tmp_path, count)
    monkeypatch.setattr(binderscope_submission, "ReadOnlyLibrary", lambda: database)
    dispatched = []
    monkeypatch.setattr(workflow_queue, "spawn_worker_for_run", lambda run: dispatched.append(run))
    result = api.submit_workflow(payload, base_dir=None, user=USER, credentials=None, idempotency_key="one-cohort")
    assert len(dispatched) == 1
    return payload, result, Path(result["job"]["run_dir"])


def test_505_designs_queue_one_32_cpu_job_without_staging_and_retry_reuses_it(database, tmp_path, monkeypatch):
    from mn_molecular_design_studio import api, binderscope_submission, complex_analysis_sources
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json

    payload = cohort(database, tmp_path, count=505)
    monkeypatch.setattr(complex_analysis_sources, "build_saved_structure_collection", lambda *args, **kwargs: pytest.fail("HTTP must not prepare structures"))
    monkeypatch.setattr(database, "list_molecules", lambda *args, **kwargs: pytest.fail("HTTP must not hydrate the catalog"))
    payload, result, run = queued(database, tmp_path, monkeypatch, payload=payload)
    assert api.validate_workflow(payload, base_dir=None, user=USER, credentials=None)["designs"] == 505
    assert api.workflow_placement_options(payload, base_dir=None, user=USER, credentials=None)["placement"] == "local"
    assert api.submit_workflow(payload, base_dir=None, user=USER, credentials=None, idempotency_key="one-cohort") == result
    request = read_json(run / "worker_request.json")
    assert request["function"] == "run_saved_binderscope"
    assert request["kwargs"]["cpu_cores"] == request["resources"]["cpu_cores"] == 32
    assert read_json(run / "metadata.json")["resource_request"]["cpu_cores"] == 32
    assert request["resources"]["gpu_device"] == "none"
    assert read_json(run / "metadata.json")["database_collection_policy"] == "automatic"
    assert run.is_relative_to(tmp_path / "runs")
    selection = Path(request["kwargs"]["selection_file"]["$path"])
    assert len(read_json(selection)["sources"]) == 505
    assert len(list(run.parent.iterdir())) == 1
    changed = copy.deepcopy(payload)
    changed["parameters"]["cpu_cores"] = 4
    with pytest.raises(Exception, match="different request"):
        api.submit_workflow(changed, base_dir=None, user=USER, credentials=None, idempotency_key="one-cohort")


def test_worker_prepares_then_computes_in_the_same_job(database, tmp_path, monkeypatch):
    from mn_molecular_design_studio import auth, complex_analysis_sources
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import create_job, read_json
    from mn_molecular_design_studio.engine.workflow_engine.core.workflow_queue import run_workflow_call
    from mn_molecular_design_studio.engine.workflow_engine.workflows import binderscope

    payload, result, run = queued(database, tmp_path, monkeypatch, count=3)
    monkeypatch.setattr(auth, "list_users", lambda: [{**USER, "disabled": False, "must_change_password": False}])
    from mn_molecular_design_studio import data_browser
    monkeypatch.setattr(data_browser, "ReadOnlyLibrary", lambda: database)
    calls = []

    def compute(source, candidates, **settings):
        calls.append(settings)
        assert len(candidates.read_text().splitlines()) == 3
        return create_job("analysis", "analysis", "binderscope", {"source": str(source)}, settings).run_dir

    monkeypatch.setattr(binderscope, "run_binderscope", compute)
    assert run_workflow_call(run, read_json(run / "worker_request.json")) == run
    assert calls[0]["cpu_cores"] == 32
    assert len(list(run.parent.iterdir())) == 1
    assert read_json(run / "metadata.json")["studio_owner_id"] == USER["id"]


@pytest.mark.parametrize("problem", ["tampered", "disabled", "access", "cancelled"])
def test_worker_rechecks_selection_and_current_authority_before_computing(database, tmp_path, monkeypatch, problem):
    from mn_molecular_design_studio import auth, data_browser
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
    from mn_molecular_design_studio.engine.workflow_engine.core.workflow_queue import run_workflow_call
    from mn_molecular_design_studio.engine.workflow_engine.workflows import binderscope

    payload, result, run = queued(database, tmp_path, monkeypatch)
    monkeypatch.setattr(auth, "list_users", lambda: [{**USER, "disabled": problem == "disabled", "must_change_password": False}])
    monkeypatch.setattr(data_browser, "ReadOnlyLibrary", lambda: database)
    monkeypatch.setattr(binderscope, "run_binderscope", lambda *args, **kwargs: pytest.fail("Invalid source must not compute"))
    request = read_json(run / "worker_request.json")
    if problem == "tampered":
        Path(request["kwargs"]["selection_file"]["$path"]).write_text("{}")
    if problem == "access":
        with database.connection() as connection:
            connection.execute("UPDATE molecule_groups SET owner_id=?", (OTHER["id"],))
    if problem == "cancelled":
        (run / "cancellation_request.json").write_text("{}")
    with pytest.raises((PermissionError, ValueError, RuntimeError)):
        run_workflow_call(run, request)


def test_http_authorization_rejects_cross_owner_missing_and_mixed_project_inputs(database, tmp_path):
    from mn_molecular_design_studio.binderscope_submission import validate_submission

    payload = cohort(database, tmp_path)
    with pytest.raises(PermissionError):
        validate_submission(payload, OTHER, None, database=database)
    with pytest.raises(PermissionError):
        validate_submission(payload, USER, "another-project", database=database)
    payload["binderscope_collection"]["sources"][0]["structure_id"] = "missing"
    with pytest.raises(PermissionError):
        validate_submission(payload, USER, None, database=database)


@pytest.mark.skipif(os.environ.get("MDS_TEST_BINDERSCOPE_IMAGE") != "1", reason="Explicit offline queued-cohort integration")
def test_one_queued_job_prepares_and_runs_the_real_offline_engine(database, tmp_path, monkeypatch):
    from mn_molecular_design_studio import auth, data_browser
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
    from mn_molecular_design_studio.engine.workflow_engine.core.workflow_queue import run_workflow_call
    from mn_molecular_design_studio.molecular_database_ingest import ingest_run

    payload = cohort(database, tmp_path)
    payload["parameters"]["cpu_cores"] = 1
    payload, result, run = queued(database, tmp_path, monkeypatch, payload=payload)
    monkeypatch.setattr(auth, "list_users", lambda: [{**USER, "disabled": False, "must_change_password": False}])
    monkeypatch.setattr(data_browser, "ReadOnlyLibrary", lambda: database)
    assert run_workflow_call(run, read_json(run / "worker_request.json")) == run
    assert read_json(run / "result.json")["success"]
    assert len(list(run.parent.iterdir())) == 1
    assert ingest_run(run, database=database, force=True) == 1
    with database.connection() as connection:
        assert connection.execute("SELECT count(*) FROM molecule_observations WHERE namespace='binderscope'").fetchone()[0] == 1
