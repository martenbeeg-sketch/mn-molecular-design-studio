import copy
import csv
import hashlib
import io
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Barrier

import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient

from mn_molecular_design_studio import job_catalog, molecule_imports
from mn_molecular_design_studio.complex_analysis_batch import assessment_protocol, evaluation_coverage, scientific_parameters
from mn_molecular_design_studio.complex_analysis_design_review import design_review_page, remember_review_selection, review_measurements
from mn_molecular_design_studio.complex_analysis_shortlist_export import build_shortlist_record, build_shortlist_manifest, manifest_csv
from mn_molecular_design_studio.complex_analysis_sources import build_database_analysis_source
from mn_molecular_design_studio.complex_analysis_ui import _matching_database_analysis_jobs
from mn_molecular_design_studio.comparisons import compatibility_warnings, compare
from mn_molecular_design_studio.evaluations import evaluation_descriptor
from mn_molecular_design_studio.molecular_database import MolecularDatabase, ReviewConflictError
from test_complex_analysis import _two_protein_chain_pdb


USER = {"id": "a" * 32, "role": "user", "_cluster_project_ids": []}


def prepared_assessment(database_subject):
    from mn_molecular_design_studio import api
    selection = database_subject[-1]
    return api._prepare_database_analysis_request({"schema_version": 1, "workflow": "analysis.biotite_interface",
        "database_source": selection, "parameters": {"contact_distance": 5.0}}, USER, None)


@pytest.mark.parametrize("contact_distance", [None, 5, 5.0])
def test_admission_reuses_effective_defaults_across_request_keys(database_subject, contact_distance):
    from mn_molecular_design_studio.scientific_admission import submit_database_assessment
    first = prepared_assessment(database_subject)
    second = copy.deepcopy(first)
    second["parameters"].pop("contact_distance")
    if contact_distance is not None:
        second["parameters"]["contact_distance"] = contact_distance
    calls = []
    def operation():
        calls.append(True)
        return {"job": {"run_id": "one-job"}}
    assert submit_database_assessment(USER, None, first, operation) == submit_database_assessment(USER, None, second, operation)
    assert len(calls) == 1


def test_admitted_inputs_survive_later_repreparation(database_subject):
    from mn_molecular_design_studio.scientific_admission import submit_database_assessment
    payload = prepared_assessment(database_subject)
    original = Path(payload["parameters"]["source_run_dir"])
    original_manifest = (original / "complex_analysis_source.json").read_bytes()
    submit_database_assessment(USER, None, payload, lambda: {"job": {"run_id": "saved"}})
    frozen = Path(payload["parameters"]["source_run_dir"])
    assert frozen != original
    saved_manifest = (frozen / "complex_analysis_source.json").read_bytes()
    (original / "complex_analysis_source.json").write_text("{}")
    assert (frozen / "complex_analysis_source.json").read_bytes() == saved_manifest
    candidate = json.loads(Path(payload["parameters"]["candidates_jsonl"]).read_text())
    assert Path(candidate["complex_pdb"]).is_relative_to(frozen)
    assert Path(candidate["complex_pdb"]).read_bytes() == Path(payload["parameters"]["selected_structure_path"]).read_bytes()
    (original / "complex_analysis_source.json").write_bytes(original_manifest)


def test_invalid_resources_cannot_poison_scientific_admission(database_subject, monkeypatch):
    from mn_molecular_design_studio import api
    from mn_molecular_design_studio.engine.workflow_engine.core import workflow_queue
    monkeypatch.setattr(workflow_queue, "spawn_worker_for_run", lambda *args, **kwargs: None)
    api.app.dependency_overrides[api._current_user] = lambda: USER
    payload = {"schema_version": 1, "workflow": "analysis.biotite_interface",
               "database_source": database_subject[-1], "parameters": {"contact_distance": 5.0}}
    try:
        with TestClient(api.app) as client:
            invalid = client.post("/api/v1/workflows", json={**payload, "resources": {"cpu_cores": 0}})
            assert invalid.status_code == 422, invalid.text
            valid = client.post("/api/v1/workflows", json=payload)
            assert valid.status_code == 202, valid.text
    finally:
        api.app.dependency_overrides.clear()


def test_admission_uncertainty_survives_a_new_caller(database_subject):
    from mn_molecular_design_studio.api_idempotency import Conflict
    from mn_molecular_design_studio.scientific_admission import submit_database_assessment
    payload = prepared_assessment(database_subject)
    def accepted_then_lost_reply():
        raise RuntimeError("accepted by scheduler, reply lost")
    with pytest.raises(RuntimeError):
        submit_database_assessment(USER, None, payload, accepted_then_lost_reply)
    with pytest.raises(Conflict, match="uncertain"):
        submit_database_assessment(USER, None, copy.deepcopy(payload), lambda: pytest.fail("Duplicate admitted after uncertainty"))


def test_admission_subject_identity_ignores_private_staging_paths(database_subject):
    from mn_molecular_design_studio.scientific_admission import submit_database_assessment
    payload = prepared_assessment(database_subject)
    path = Path(payload["parameters"]["source_run_dir"]) / "complex_analysis_source.json"
    manifest = json.loads(path.read_text())
    manifest["reference_context"] = {"reference_structure_path": "/user-one/reference.pdb", "sha256": "c" * 64}
    path.write_text(json.dumps(manifest))
    expected = {"job": {"run_id": "one-job"}}
    assert submit_database_assessment(USER, None, payload, lambda: expected) == expected
    manifest["reference_context"]["reference_structure_path"] = "/user-two/reference.pdb"
    path.write_text(json.dumps(manifest))
    assert submit_database_assessment(USER, None, payload, lambda: pytest.fail("Private paths split scientific identity")) == expected


def test_admission_deduplicates_project_sessions_without_merging_private_owners(database_subject):
    from mn_molecular_design_studio.scientific_admission import submit_database_assessment
    payload = prepared_assessment(database_subject)
    peer_user = {**USER, "id": "b" * 32}
    first = submit_database_assessment(USER, "project", payload, lambda: {"job": {"run_id": "shared-job"}})
    assert submit_database_assessment(peer_user, "project", payload, lambda: pytest.fail("Duplicate project work")) == first
    assert submit_database_assessment(peer_user, None, payload, lambda: {"job": {"run_id": "private-job"}}) != first


def test_admission_requires_exact_previous_job_for_explicit_reassessment(database_subject):
    from mn_molecular_design_studio.api_idempotency import Conflict
    from mn_molecular_design_studio.scientific_admission import submit_database_assessment
    payload = prepared_assessment(database_subject)
    submit_database_assessment(USER, None, payload, lambda: {"job": {"run_id": "original"}})
    retry = {**payload, "execution": {"reassessment_of": "original"}}
    expected = {"job": {"run_id": "retry"}}
    assert submit_database_assessment(USER, None, retry, lambda: expected) == expected
    assert submit_database_assessment(USER, None, retry, lambda: pytest.fail("Retry duplicated")) == expected
    retry["execution"]["reassessment_of"] = "foreign-job"
    with pytest.raises(Conflict):
        submit_database_assessment(USER, None, retry, lambda: pytest.fail("Unauthorized retry"))


def test_replica_assessment_admission_is_fenced_to_its_source_owner(database_subject, monkeypatch):
    from mn_molecular_design_studio import complex_analysis_sources
    database, _imported, _saved, _structure, selection = database_subject
    adapter = build_database_analysis_source(USER, selection, workflow_id="analysis.biotite_interface", project_id=None, database=database)
    adapter["source"]["structure_context"]["_replica_node"] = "another-workstation"
    monkeypatch.setattr(complex_analysis_sources, "build_database_analysis_source", lambda *args, **kwargs: adapter)
    with pytest.raises(ValueError, match="source-owning workstation"):
        complex_analysis_sources.apply_database_analysis_source({"workflow": "analysis.biotite_interface",
            "database_source": selection, "parameters": {}}, USER, None)


@pytest.mark.parametrize("terminal", ["failed", "completed"])
def test_real_queue_reassessment_requires_a_terminal_exact_job(database_subject, monkeypatch, terminal):
    from mn_molecular_design_studio import api
    from mn_molecular_design_studio.engine.workflow_engine.core import workflow_queue
    from mn_molecular_design_studio.runtime import runs_root
    monkeypatch.setattr(workflow_queue, "spawn_worker_for_run", lambda *args, **kwargs: None)
    api.app.dependency_overrides[api._current_user] = lambda: USER
    payload = {"schema_version": 1, "workflow": "analysis.biotite_interface",
               "database_source": database_subject[-1], "parameters": {"contact_distance": 5.0}}
    try:
        with TestClient(api.app) as client:
            first = client.post("/api/v1/workflows", json=payload)
            assert first.status_code == 202, first.text
            run_id = first.json()["job"]["run_id"]
            retry = {**payload, "execution": {"placement": "local", "reassessment_of": run_id}}
            assert client.post("/api/v1/workflows", json=retry).status_code == 409
            run = Path(runs_root()) / "analysis" / run_id
            metadata = json.loads((run / "metadata.json").read_text())
            (run / "metadata.json").write_text(json.dumps({**metadata, "status": terminal}))
            job_catalog._upsert(run.parent.parent, run, refresh_artifacts=False)
            response = client.post("/api/v1/workflows", json=retry)
            assert response.status_code == 202, response.text
            assert response.json()["job"]["run_id"] != run_id
            assert client.post("/api/v1/workflows", json=retry).json()["job"]["run_id"] == response.json()["job"]["run_id"]
    finally:
        api.app.dependency_overrides.clear()


def test_cross_page_selection_is_reauthorized_and_does_not_broaden_campaign(database_subject):
    from mn_molecular_design_studio.complex_analysis_design_review import selected_review_models
    database, imported, saved, structure, _selection = database_subject
    other_group = database.create_group(USER, "Different source")
    database.add_to_group(USER, other_group, [saved["id"]])
    extra = database.attach_structure(other_group, saved["id"], {"name": "other.pdb", "content": b"END\n"}, context={"model_id": "other"})
    selected = selected_review_models(database, USER, [imported["group_id"]], [structure["id"], extra])
    assert [sample["structure_id"] for sample in selected] == [structure["id"]]
    assert selected_review_models(database, {**USER, "id": "b" * 32}, [], [structure["id"]]) == []


def test_different_assessments_cannot_overwrite_each_others_source_manifest(database_subject):
    database, _imported, _saved, _structure, selection = database_subject
    first = build_database_analysis_source(USER, selection, workflow_id="analysis.biotite_interface", project_id=None, database=database)
    path = Path(first["source_run_dir"]) / "complex_analysis_source.json"
    content = path.read_bytes()
    second = build_database_analysis_source(USER, selection, workflow_id="analysis.binder_interface_chemistry", project_id=None, database=database)
    assert first["source_run_dir"] != second["source_run_dir"]
    assert path.read_bytes() == content


def test_real_concurrent_preparation_and_admission_create_one_job(database_subject, monkeypatch):
    from mn_molecular_design_studio import api
    from mn_molecular_design_studio.engine.workflow_engine.core import workflow_queue
    monkeypatch.setattr(workflow_queue, "spawn_worker_for_run", lambda *args, **kwargs: None)
    original_prepare = api._prepare_database_analysis_request
    barrier = Barrier(2)
    def synchronized_prepare(*args, **kwargs):
        prepared = original_prepare(*args, **kwargs)
        barrier.wait(timeout=20)
        return prepared
    monkeypatch.setattr(api, "_prepare_database_analysis_request", synchronized_prepare)
    api.app.dependency_overrides[api._current_user] = lambda: USER
    payload = {"schema_version": 1, "workflow": "analysis.biotite_interface",
               "database_source": database_subject[-1], "parameters": {"contact_distance": 5.0}}
    try:
        with TestClient(api.app) as client:
            with ThreadPoolExecutor(2) as pool:
                responses = list(pool.map(lambda _index: client.post("/api/v1/workflows", json=payload), range(2)))
        assert [response.status_code for response in responses] == [202, 202]
        assert len({response.json()["job"]["run_id"] for response in responses}) == 1
    finally:
        api.app.dependency_overrides.clear()


@pytest.fixture
def database_subject(tmp_path, monkeypatch):
    workspace = tmp_path / "studio"
    workspace.mkdir()
    (workspace / "settings.json").write_text(json.dumps({"results_directory": str(workspace / "runs"), "reference_directory": str(workspace / "reference")}))
    for suffix, path in {"HOME": workspace, "APP_HOME": workspace, "RUN_DIR": workspace / "runs",
                         "REFERENCE_DIR": workspace / "reference"}.items():
        monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_" + suffix, str(path))
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(workspace / "cluster"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_AUTH_MODE", "production")
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_SCHEDULER_MODE", "local")
    database = MolecularDatabase(workspace / "database")
    record = molecule_imports.parse_structure(_two_protein_chain_pdb("ACDEFG", "KLMNPQ").encode(), "source.pdb")[0]
    record["files"][0].update(model_id="0", candidate_id="design-1", binder_chains=["B"], target_chains=["A"], run_id="source")
    imported = database.import_records(USER, [record], name="Audit designs")
    saved = database.list_molecules(USER, molecule_ids=imported["molecule_ids"])["items"][0]
    structure = saved["representations"][0]
    selection = {"molecule_id": saved["id"], "group_id": imported["group_id"], "structure_id": structure["id"],
                 "model_id": "0", "binder_chains": ["B"], "target_chains": ["A"]}
    return database, imported, saved, structure, selection


def saved_evidence(*, protocol=None, metrics=None, **settings):
    protocol = protocol or assessment_protocol("analysis.biotite_interface")
    observation = {"id": "observation", "namespace": "biotite_interface", "run_id": "analysis",
                   "metrics": metrics or {"qc_interface_contact_pairs": 2}, "context": {
                       "source_structure_id": "s", "source_structure_digest": "a" * 64, "representation_id": "s",
                       "model_id": "0", "binder_chains": ["B"], "target_chains": ["A"],
                       "settings": {"contact_distance": 5.0, "assessment_protocol": protocol, **settings},
                       "status": "completed"}}
    descriptor = evaluation_descriptor(observation)
    return observation, descriptor


def matching_row(source, *, protocol=None, parameters=None, tool="biotite_interface"):
    return {"run_id": "existing", "task_group": "analysis", "status": "completed", "tool": tool,
            "owner_id": USER["id"], "run_dir": Path("/unused"), "database_analysis": {
                "source_key": job_catalog.analysis_source_key(source), "source": source,
                "candidate_id": "candidate", "model_id": "0", "binder_chains": ["B"], "target_chains": ["A"],
                "parameters": parameters or {"contact_distance": 5.0}, "protocol": protocol or {"method_version": "obsolete"}}}


def test_inspector_reuse_must_reject_obsolete_protocol():
    source = {"molecule_id": "m", "structure_id": "s", "group_id": "g", "model_id": "0", "structure_digest": "a" * 64}
    arguments = dict(workflow_id="analysis.biotite_interface", source=source, candidate_id="candidate", model_id="0",
                     binder_chains=["B"], target_chains=["A"], parameters={"contact_distance": 5.0},
                     indexed_rows=[matching_row(source)], load_evidence=False)
    assert _matching_database_analysis_jobs(USER, **arguments, expected_protocol={"method_version": "current"}) == []
    actual = _matching_database_analysis_jobs(USER, **arguments)
    print("Inspector default matcher accepted obsolete protocol:", len(actual))
    assert actual == []


def test_known_different_qc_image_digest_must_not_be_compatible():
    observation, descriptor = saved_evidence(image_digest="sha256:old-runtime")
    sample = {"structure_id": "s", "digest": "a" * 64, "model_id": "0", "associated": [(observation, descriptor)]}
    result = evaluation_coverage(sample, "analysis.biotite_interface", parameters={"contact_distance": 5.0},
                                 binder_chains=["B"], target_chains=["A"])
    print("Registered expected protocol:", assessment_protocol("analysis.biotite_interface"), "coverage:", result)
    assert result[0] != "Already evaluated · compatible settings"


def test_subject_paths_must_not_split_scientific_protocols():
    records = []
    for index in range(2):
        observation, descriptor = saved_evidence(selected_structure_path=f"/input/design-{index}.pdb",
                                                selected_model_id=f"model-{index}", candidate_ids=[f"candidate-{index}"])
        records.append({"structure_id": f"structure-{index}", "associated": [(observation, descriptor)]})
    signatures = [sample["associated"][0][1]["parameters_digest"] for sample in records]
    print("Equal scientific settings, subject-specific parameter digests:", signatures)
    assert scientific_parameters(records[0]["associated"][0][0]["context"]["settings"]) == scientific_parameters(records[1]["associated"][0][0]["context"]["settings"])
    series, values = review_measurements(records)
    print("Measurement columns:", len(series), "per-model column IDs:", values)
    assert signatures[0] == signatures[1]


def test_comparison_must_reject_distinct_target_revisions():
    rows = [{"molecule_id": f"m{index}", "name": f"model {index}", "scores": {"binding": index + 0.5},
             "evaluation": {"id": f"e{index}", "status": "completed", "target_id": "same-target",
                            "target_revision": str(index) * 64, "parameters_digest": "same-protocol"}}
            for index in range(2)]
    warnings = compatibility_warnings(rows, ["binding"])
    ranked = compare(rows, fields=["binding"], directions={"binding": "higher"})
    print("Target revision warnings:", warnings, "ranks:", ranked["Comparison rank"].tolist())
    assert warnings


def test_curated_membership_must_retain_reviewable_models(database_subject):
    database, imported, saved, structure, selection = database_subject
    group = database.create_group(USER, "Saved shortlist")
    database.add_to_group(USER, group, imported["molecule_ids"])
    source = design_review_page(database, USER, [imported["group_id"]])
    curated = design_review_page(database, USER, [group])
    print("Source models:", source["total"], "curated memberships:", len(database.list_molecules(USER, group_id=group)["items"]), "curated models:", curated["total"])
    assert curated["total"] == source["total"]


def test_page_changes_must_keep_authorized_selected_models():
    state = {"design-review-selected-cohort": ["page-one"], "complex-analysis-workflow-samples": ["hydrated-one"]}
    selected = remember_review_selection(state, "cohort", ["page-one"], {"page-two"})
    print("Selection after paging:", selected, "hydrated cache:", state.get("complex-analysis-workflow-samples"))
    assert selected == ["page-one"]


def test_metadata_review_export_must_retain_known_sequences(database_subject):
    database, imported, saved, structure, selection = database_subject
    sample = design_review_page(database, USER, [imported["group_id"]])["items"][0]
    exported = build_shortlist_record(sample)
    print("Stored components:", [(item["chain"], item.get("sequence")) for item in sample["components"]])
    print("Exported binder:", exported["binder"], "target:", exported["target"])
    assert exported["binder"]["sequences_by_chain"] == {"B": "KLMNPQ"}
    assert exported["target"]["sequences_by_chain"] == {"A": "ACDEFG"}


@pytest.mark.parametrize("note", ["=1+1", "+1+1", "-1+1", "@SUM(1,1)", "\t=1+1", "\r=1+1"])
def test_csv_must_neutralize_formula_prefixes(database_subject, note):
    database, imported, saved, structure, selection = database_subject
    sample = design_review_page(database, USER, [imported["group_id"]])["items"][0]
    sample["review_annotation"] = {"status": "Reviewed", "note": note}
    manifest = build_shortlist_manifest([sample])
    serialized = manifest_csv(manifest)
    content = serialized.decode() if isinstance(serialized, bytes) else serialized
    cells = list(csv.DictReader(io.StringIO(content)))
    print("Unescaped human text in CSV:", repr(note), [key for key, value in cells[0].items() if value == note])
    assert note not in cells[0].values()


def test_numeric_zero_model_must_survive_review_projection(database_subject):
    database, imported, saved, structure, selection = database_subject
    with database.connection() as connection:
        context = json.loads(connection.execute("SELECT context FROM molecule_structures WHERE id=?", (structure["id"],)).fetchone()[0])
        context["model_id"] = 0
        connection.execute("UPDATE molecule_structures SET context=? WHERE id=?", (json.dumps(context), structure["id"]))
    sample = design_review_page(database, USER, [imported["group_id"]])["items"][0]
    print("Retained numeric model 0 projected as:", repr(sample["model_id"]))
    assert sample["model_id"] == "0"


def test_failed_qc_item_must_not_ingest_as_compatible_success(database_subject, tmp_path):
    from mn_molecular_design_studio.molecular_database_ingest import ingest_run
    database, imported, saved, structure, selection = database_subject
    adapter = build_database_analysis_source(USER, selection, workflow_id="analysis.biotite_interface", project_id=None, database=database)
    run = tmp_path / "studio/runs/analysis/invalid-qc"
    (run / "artifacts/normalized_candidates").mkdir(parents=True)
    candidate = json.loads(Path(adapter["candidates_jsonl"]).read_text())
    candidate.update(candidate_id=candidate["candidate_id"] + "_biotite_interface", source_tool="biotite_interface", tool="biotite_interface", stage="analysis",
                     metrics={"qc_status": "completed", "qc_error": "Required measurements unavailable", "analysis_status": "error"})
    candidate["raw_metadata"]["measurement_status"] = {"status": "missing", "error": "Required measurements unavailable"}
    payloads = {
        "metadata.json": {"status": "completed", "tool": "biotite_interface", "database_collection_policy": "automatic", "upstream_source_kind": "database", "database_source": adapter["source_manifest"]},
        "molecular_design_studio.json": {"application": "mn-molecular-design-studio", "studio_owner_id": USER["id"]},
        "input.json": {"params": {"contact_distance": 5.0, "assessment_protocol": assessment_protocol("analysis.biotite_interface")}},
        "result.json": {"success": True},
    }
    for name, value in payloads.items():
        (run / name).write_text(json.dumps(value))
    (run / "artifacts/normalized_candidates/candidates.jsonl").write_text(json.dumps(candidate) + "\n")
    assert ingest_run(run, database=database, force=True) == 1
    observation = next(item for item in database.list_molecules(USER, molecule_ids=imported["molecule_ids"])["items"][0]["observations"] if item["namespace"] == "biotite_interface")
    descriptor = evaluation_descriptor(observation)
    sample = {"structure_id": structure["id"], "digest": adapter["source"]["structure_digest"], "model_id": "0", "associated": [(observation, descriptor)]}
    coverage = evaluation_coverage(sample, "analysis.biotite_interface", parameters={"contact_distance": 5.0}, binder_chains=["B"], target_chains=["A"])
    print("Failed scientific item persisted status/error:", descriptor["status"], observation["context"]["error"], "coverage:", coverage)
    assert coverage[0] != "Already evaluated · compatible settings"


def test_local_review_cas_survives_real_delete_and_restart(database_subject):
    database, imported, saved, structure, selection = database_subject
    identity = {"source_kind": "database", "subject_id": saved["id"], "structure_ref": structure["id"]}
    absent = database.get_complex_review(USER, **identity)
    stored = database.save_complex_review(USER, **identity, note="new rationale", status="Reviewed", expected_revision=absent["revision"])
    with database.connection() as connection:
        connection.execute("DELETE FROM complex_review_annotations WHERE subject_id=?", (saved["id"],))
    restarted = MolecularDatabase(database.root)
    for stale in (absent, stored):
        with pytest.raises(ReviewConflictError):
            restarted.save_complex_review(USER, **identity, note="stale", expected_revision=stale["revision"])


def test_integrity_snapshot_and_changed_source(database_subject):
    database, imported, saved, structure, selection = database_subject
    snapshot = database.structure_file(USER, structure["id"])
    original = snapshot.read_bytes()
    source = Path(structure["path"])
    source.write_bytes(original + b"REMARK fixture changed\n")
    assert snapshot.read_bytes() == original
    with pytest.raises(ValueError, match="digest mismatch"):
        database.structure_file(USER, structure["id"])


def test_same_scientific_submission_cannot_create_two_queue_jobs(database_subject, monkeypatch):
    from mn_molecular_design_studio import api
    from mn_molecular_design_studio.engine.workflow_engine.core import workflow_queue
    database, imported, saved, structure, selection = database_subject
    monkeypatch.setattr(workflow_queue, "spawn_worker_for_run", lambda *args, **kwargs: None)
    api.app.dependency_overrides[api._current_user] = lambda: USER
    try:
        with TestClient(api.app) as client:
            payload = {"schema_version": 1, "workflow": "analysis.biotite_interface", "database_source": selection,
                       "parameters": {"contact_distance": 5.0}}
            first = client.post("/api/v1/workflows", json=payload)
            second = client.post("/api/v1/workflows", json=payload)
            print("Identical API submissions:", first.status_code, first.text[:350], second.status_code, second.text[:350])
            assert first.status_code == second.status_code == 202
            assert first.json()["job"]["run_id"] == second.json()["job"]["run_id"]
    finally:
        api.app.dependency_overrides.clear()


def test_real_native_chemistry_result_must_be_reusable(database_subject):
    from mn_molecular_design_studio.engine.workflow_engine.workflows import analysis
    from mn_molecular_design_studio.molecular_database_ingest import ingest_run
    database, imported, saved, structure, selection = database_subject
    adapter = build_database_analysis_source(USER, selection, workflow_id="analysis.binder_interface_chemistry", project_id=None, database=database)
    run = analysis.run_binder_interface_chemistry_contract(
        Path(adapter["source_run_dir"]), Path(adapter["candidates_jsonl"]),
        candidate_ids=[adapter["candidate_id"]], selected_structure_path=Path(adapter["selected_structure_path"]),
        selected_model_id="0", binder_chains=["B"], target_chains=["A"])
    (run / "molecular_design_studio.json").write_text(json.dumps({"application": "mn-molecular-design-studio", "studio_owner_id": USER["id"]}))
    ingested = ingest_run(run, database=database, force=True)
    print("Real chemistry ingestion count:", ingested)
    sample = design_review_page(database, USER, [imported["group_id"]])["items"][0]
    result = evaluation_coverage(sample, "analysis.binder_interface_chemistry", parameters={}, binder_chains=["B"], target_chains=["A"])
    observation = next(item for item in sample["associated"] if item[1]["evaluator"] == "binder_interface_chemistry_v1")
    print("Real chemistry job result:", json.loads((run / "result.json").read_text())["success"])
    print("Requested chemistry settings:", {}, "persisted scientific settings:", scientific_parameters(observation[0]["context"]["settings"]), "coverage:", result)
    job_catalog._upsert(run.parent.parent, run, refresh_artifacts=False)
    matches = _matching_database_analysis_jobs(USER, workflow_id="analysis.binder_interface_chemistry", source=adapter["source_manifest"],
        candidate_id=adapter["candidate_id"], model_id="0", binder_chains=["B"], target_chains=["A"], parameters={}, load_evidence=False)
    print("Completed real chemistry job index match count:", len(matches))
    assert result[0] == "Already evaluated · compatible settings"
    assert len(matches) == 1


def test_same_request_idempotency_key_reuses_one_queue_job(database_subject, monkeypatch):
    from mn_molecular_design_studio import api
    from mn_molecular_design_studio.engine.workflow_engine.core import workflow_queue
    database, imported, saved, structure, selection = database_subject
    monkeypatch.setattr(workflow_queue, "spawn_worker_for_run", lambda *args, **kwargs: None)
    api.app.dependency_overrides[api._current_user] = lambda: USER
    try:
        with TestClient(api.app) as client:
            payload = {"schema_version": 1, "workflow": "analysis.biotite_interface", "database_source": selection,
                       "parameters": {"contact_distance": 5.0}}
            responses = [client.post("/api/v1/workflows", json=payload, headers={"Idempotency-Key": "stable-request"}) for _ in range(2)]
            assert all(response.status_code == 202 for response in responses)
            assert responses[0].json()["job"]["run_id"] == responses[1].json()["job"]["run_id"]
    finally:
        api.app.dependency_overrides.clear()


def test_concurrent_api_sessions_must_not_duplicate_scientific_work(database_subject, monkeypatch):
    from mn_molecular_design_studio import api
    from mn_molecular_design_studio.engine.workflow_engine.core import workflow_queue
    database, imported, saved, structure, selection = database_subject
    monkeypatch.setattr(workflow_queue, "spawn_worker_for_run", lambda *args, **kwargs: None)
    original_prepare = api._prepare_database_analysis_request
    payload = {"schema_version": 1, "workflow": "analysis.biotite_interface", "database_source": selection,
               "parameters": {"contact_distance": 5.0}}
    prepared = original_prepare(copy.deepcopy(payload), USER, None)
    monkeypatch.setattr(api, "_prepare_database_analysis_request", lambda *args, **kwargs: copy.deepcopy(prepared))
    barrier = Barrier(2)
    submit_once = api._submit_once
    def synchronized_admission(*args, **kwargs):
        barrier.wait(timeout=10)
        return submit_once(*args, **kwargs)
    monkeypatch.setattr(api, "_submit_once", synchronized_admission)
    api.app.dependency_overrides[api._current_user] = lambda: USER
    try:
        with TestClient(api.app) as client:
            def request(index):
                return client.post("/api/v1/workflows", json=payload, headers={"Idempotency-Key": f"independent-session-{index}"})
            with ThreadPoolExecutor(2) as pool:
                responses = list(pool.map(request, range(2)))
        print("Concurrent independent-session statuses:", [response.status_code for response in responses])
        assert all(response.status_code == 202 for response in responses)
        identifiers = [response.json()["job"]["run_id"] for response in responses]
        print("Concurrent independent-session job IDs:", identifiers)
        assert len(set(identifiers)) == 1
    finally:
        api.app.dependency_overrides.clear()


def test_real_catalog_upgrade_blocks_then_reconciles(database_subject, monkeypatch):
    database, imported, saved, structure, selection = database_subject
    adapter = build_database_analysis_source(USER, selection, workflow_id="analysis.biotite_interface", project_id=None, database=database)
    source = adapter["source_manifest"]
    root = database.root.parent / "runs"
    run = root / "analysis/legacy-analysis"
    run.mkdir(parents=True)
    for name, value in {"metadata.json": {"status": "collecting", "tool": "biotite_interface", "studio_owner_id": USER["id"], "database_source": source},
                        "input.json": {"params": {"contact_distance": 5.0, "selected_model_id": "0", "binder_chains": ["B"], "target_chains": ["A"], "candidate_ids": [adapter["candidate_id"]]}}}.items():
        (run / name).write_text(json.dumps(value))
    job_catalog._upsert(root, run, refresh_artifacts=False)
    with job_catalog._connection() as connection:
        connection.execute("UPDATE jobs SET run_identity_json='{}'")
        connection.execute("UPDATE catalog_projection SET value='old-projection'")
    monkeypatch.setattr(job_catalog, "_SCHEMA_READY_PATH", None)
    with pytest.raises(RuntimeError, match="upgrade"):
        job_catalog.database_analysis_jobs(USER, [source])
    job_catalog._reconcile_locked(root)
    rows = job_catalog.database_analysis_jobs(USER, [source])
    assert len(rows) == 1 and rows[0]["status"] == "collecting"


def test_real_sqlite_lookup_error_must_render_an_unavailable_state(monkeypatch):
    import sqlite3
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import molecular_database
    from test_assessment_submission_audit import SCRIPT
    class Database:
        def group(self, *args, **kwargs):
            return {"project_id": "", "origin": {}}
    monkeypatch.setattr(molecular_database, "MolecularDatabase", Database)
    def locked(*args, **kwargs):
        raise sqlite3.OperationalError("database is locked")
    monkeypatch.setattr(job_catalog, "database_analysis_jobs", locked)
    app = AppTest.from_string(SCRIPT).run()
    print("Real SQLite exception in planner:", [exception.message for exception in app.exception])
    assert not app.exception


def test_results_download_must_reject_recorded_digest_mismatch(database_subject):
    from mn_molecular_design_studio import api
    from mn_molecular_design_studio.engine.workflow_engine.core.artifact_collection import validate_run_artifacts
    database, imported, saved, structure, selection = database_subject
    run = database.root.parent / "runs/analysis/retained-artifact"
    run.mkdir(parents=True)
    original = Path(structure["path"]).read_bytes()
    source = run / "model.pdb"
    source.write_bytes(original)
    for name, value in {"metadata.json": {"status": "completed", "tool": "biotite_interface", "studio_owner_id": USER["id"]},
                        "input.json": {"params": {}}, "result.json": {"success": True},
                        "molecular_design_studio.json": {"application": "mn-molecular-design-studio", "studio_owner_id": USER["id"]},
                        "artifacts.json": {"schema_version": 1, "run_id": run.name, "artifacts": [{"path": "model.pdb", "sha256": hashlib.sha256(original).hexdigest(), "required": True}]}}.items():
        (run / name).write_text(json.dumps(value))
    job_catalog._upsert(run.parent.parent, run, refresh_artifacts=True)
    assert validate_run_artifacts(run)["status"] == "complete"
    changed = original.replace(b"ALA", b"GLY", 1)
    source.write_bytes(changed)
    assert validate_run_artifacts(run)["status"] != "complete"
    api.app.dependency_overrides[api._current_user] = lambda: USER
    try:
        with TestClient(api.app) as client:
            response = client.get("/api/v1/jobs/retained-artifact/artifact", params={"path": "model.pdb"})
            print("Digest-invalid retained Results artifact HTTP:", response.status_code, "changed bytes served:", response.content == changed)
            assert response.status_code in {404, 409, 422}
    finally:
        api.app.dependency_overrides.clear()


def test_bound_receptor_viewer_must_enforce_pose_coordinate_frame(database_subject, tmp_path, monkeypatch):
    from rdkit import Chem
    from mn_molecular_design_studio import presentation_api, bound_docking
    database, imported, saved, structure, selection = database_subject
    original = Path(structure["path"]).read_bytes()
    receptor = tmp_path / "legacy-receptor.pdb"
    receptor.write_bytes(original)
    imported_pose = database.import_records(USER, molecule_imports.parse_smiles("CCO pose"), name="Poses")
    molecule = imported_pose["molecule_ids"][0]
    pose = database.attach_structure(imported_pose["group_id"], molecule,
        {"name": "pose.sdf", "content": (Chem.MolToMolBlock(Chem.MolFromSmiles("CCO")) + "\n$$$$\n").encode()},
        context={"database_structure_policy": "best_docking_pose", "target_path": str(receptor), "target_structure_digest": hashlib.sha256(original).hexdigest()})
    changed = original.replace(b"ALA", b"GLY", 1)
    receptor.write_bytes(changed)
    selected = next(item for item in database.list_molecules(USER, molecule_ids=[molecule])["items"][0]["representations"] if item["id"] == pose)
    with pytest.raises(ValueError, match="coordinate-frame"):
        bound_docking._resolve_pose_receptor(database, USER, selected)
    monkeypatch.setattr("mn_molecular_design_studio.molecular_database.MolecularDatabase", lambda: database)
    app = FastAPI()
    presentation_api.install(app, lambda: USER, resolve_job=lambda *args: None, call=lambda function, *args, **kwargs: function(*args, **kwargs))
    response = TestClient(app).get(f"/api/v1/molecules/{molecule}/viewer/receptor", params={"structure_id": pose})
    print("Coordinate-frame-invalid receptor viewer HTTP:", response.status_code, "changed bytes served:", response.content == changed)
    assert response.status_code in {404, 409, 422}
