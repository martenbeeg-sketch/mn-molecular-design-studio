import copy
import hashlib
import json
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from test_post_phase2_release import USER, database_subject, prepared_assessment


def test_actual_native_ingestion_cannot_change_admission_identity(database_subject, monkeypatch):
    from mn_molecular_design_studio import api
    from mn_molecular_design_studio.engine.workflow_engine.core import local_worker, workflow_queue
    from mn_molecular_design_studio.scientific_admission import _subject_identity

    body = {"schema_version": 1, "workflow": "analysis.sequence_liability", "database_source": database_subject[-1],
            "parameters": {}, "resources": {"cpu_cores": 1, "gpu_device": "none"}}
    first_source = api._prepare_database_analysis_request(body, USER, None)
    original = json.loads((Path(first_source["parameters"]["source_run_dir"]) / "complex_analysis_source.json").read_text())
    launched = []
    def execute(run):
        launched.append(Path(run))
        assert local_worker.run_worker_job(run) == 0
    monkeypatch.setattr(workflow_queue, "spawn_worker_for_run", execute)
    api.app.dependency_overrides[api._current_user] = lambda: USER
    try:
        with TestClient(api.app) as client:
            first = client.post("/api/v1/workflows", json=body)
            assert first.status_code == 202, first.text
        next_source = api._prepare_database_analysis_request(body, USER, None)
        updated = json.loads((Path(next_source["parameters"]["source_run_dir"]) / "complex_analysis_source.json").read_text())
        assert original["database_identity_components"] != updated["database_identity_components"]
        assert _subject_identity(original) == _subject_identity(updated)
        with TestClient(api.app) as client:
            second = client.post("/api/v1/workflows", json=body)
            assert second.status_code == 202, second.text
        assert first.json()["job"]["run_id"] == second.json()["job"]["run_id"]
        assert len(launched) == 1
    finally:
        api.app.dependency_overrides.clear()


def test_only_descriptive_component_roles_are_ignored():
    from mn_molecular_design_studio.scientific_admission import _subject_identity

    subject = {"database_identity_components": [{"chain": "A", "sequence": "ACDE", "role": ""}],
               "binder_chains": ["A"], "target_chains": ["B"], "reference_context": {"role": "target"}}
    enriched = copy.deepcopy(subject)
    enriched["database_identity_components"][0]["role"] = "binder"
    assert _subject_identity(subject) == _subject_identity(enriched)
    enriched["binder_chains"], enriched["target_chains"] = ["B"], ["A"]
    assert _subject_identity(subject) != _subject_identity(enriched)
    enriched = copy.deepcopy(subject)
    enriched["database_identity_components"][0]["sequence"] = "ACDF"
    assert _subject_identity(subject) != _subject_identity(enriched)
    enriched = copy.deepcopy(subject)
    enriched["reference_context"]["role"] = "binder"
    assert _subject_identity(subject) != _subject_identity(enriched)


@pytest.mark.parametrize("outcome", ["completed", "uncertain", "missing-inputs"])
def test_legacy_admission_is_reconciled_without_losing_uncertainty(database_subject, outcome):
    from mn_molecular_design_studio.api_idempotency import Conflict, submit_once
    from mn_molecular_design_studio.complex_analysis_batch import assessment_protocol
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import workflow_schema
    from mn_molecular_design_studio.runtime import app_home
    from mn_molecular_design_studio.scientific_admission import _freeze_inputs, _subject_identity, submit_database_assessment
    from mn_molecular_design_studio.scientific_settings import scientific_parameters

    payload = prepared_assessment(database_subject)
    method = payload["workflow"]
    manifest = json.loads((Path(payload["parameters"]["source_run_dir"]) / "complex_analysis_source.json").read_text())
    defaults = {field["name"]: field["default"] for field in workflow_schema(method)["parameters"] if "default" in field}
    canonical = {"workflow": method, "source": _subject_identity(manifest, retain_component_roles=True),
                 "settings": scientific_parameters({**defaults, **payload["parameters"]}, method), "protocol": assessment_protocol(method)}
    key = hashlib.sha256(json.dumps(canonical, sort_keys=True, allow_nan=False).encode()).hexdigest()
    scope = {"id": "user:" + USER["id"]}
    root = Path(app_home()) / "scientific-admission"
    root.mkdir()
    def legacy_operation():
        if outcome != "missing-inputs":
            _freeze_inputs(payload, root, hashlib.sha256((scope["id"] + key).encode()).hexdigest())
        if outcome != "completed":
            raise RuntimeError("Legacy scheduler response lost")
        return {"job": {"run_id": "legacy-job"}}
    if outcome == "completed":
        submit_once(scope, "scientific", key, canonical, legacy_operation)
    else:
        with pytest.raises(RuntimeError):
            submit_once(scope, "scientific", key, canonical, legacy_operation)
    current = prepared_assessment(database_subject)
    path = Path(current["parameters"]["source_run_dir"]) / "complex_analysis_source.json"
    enriched = json.loads(path.read_text())
    for component in enriched["database_identity_components"]:
        component["role"] = "binder"
    path.write_text(json.dumps(enriched))
    if outcome == "completed":
        result = submit_database_assessment(USER, None, current, lambda: pytest.fail("Legacy scientific work was duplicated"))
        assert result["job"]["run_id"] == "legacy-job"
    else:
        with pytest.raises(Conflict, match="Legacy|legacy"):
            submit_database_assessment(USER, None, current, lambda: pytest.fail("Uncertain legacy work was repeated"))
