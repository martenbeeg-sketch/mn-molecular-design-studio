from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import pytest

from mn_molecular_design_studio.engine.workflow_engine.core import jobs
from mn_molecular_design_studio.engine.workflow_engine.core.artifact_collection import validate_run_artifacts
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates, write_candidates
from mn_molecular_design_studio.engine.workflow_engine.core import portable_paths


def _candidate(path: str, *, candidate_id: str = "design-1") -> dict:
    return {
        "candidate_id": candidate_id,
        "source_tool": "fixture-generator",
        "stage": "generation.backbone_sequence",
        "complex_pdb": path,
        "target_chains": ["A"],
        "binder_chains": ["B"],
        "binder_sequence": "ACDEFG",
        "raw_metadata": {"source_job": "fixture-source"},
    }


def _finish(run: Path, monkeypatch, success=True):
    monkeypatch.setattr("mn_molecular_design_studio.molecular_database_ingest.record_completed_run", lambda _path: None)
    jobs.write_json(run / "input.json", {"job_type": "fixture", "tool": "fixture-generator", "inputs": {}, "params": {}})
    jobs.write_json(run / "metadata.json", {"run_id": run.name, "status": "running", "task_group": "design"})
    jobs.finish_job(run, success, {"metrics": {"fixture": True}})


def test_container_work_paths_are_mapped_only_during_results_collection(tmp_path, monkeypatch):
    run = tmp_path / "runs" / "design" / "run-work"
    output = run / "artifacts" / "raw" / "fixture" / "output" / "model.cif"
    output.parent.mkdir(parents=True)
    output.write_text("data_fixture\n#\n")
    command = ["docker", "run", "-v", f"{run}:/work", "fixture:image", "--output", "/work/artifacts/raw/fixture/output"]
    jobs.write_json(run / "command.json", {"command": command})

    write_candidates(run, "fixture-generator", [_candidate("/work/artifacts/raw/fixture/output/model.cif")])
    persisted = json.loads((run / "artifacts/normalized_candidates/candidates.jsonl").read_text())
    assert persisted["complex_pdb"] == "artifacts/raw/fixture/output/model.cif"
    assert persisted["raw_metadata"]["container_artifact_paths"]["complex_pdb"] == "/work/artifacts/raw/fixture/output/model.cif"
    # The execution command still contains the container-native path verbatim.
    assert jobs.read_json(run / "command.json")["command"] == command
    assert read_candidates(run)[0]["complex_pdb"] == str(output)

    _finish(run, monkeypatch)
    metadata = jobs.read_json(run / "metadata.json")
    assert metadata["status"] == "completed"
    assert metadata["artifact_collection_complete"] is True
    assert metadata["artifact_validation"]["required_artifact_count"] == 1


def test_other_recorded_container_mount_is_materialized_as_portable_results(tmp_path):
    run = tmp_path / "runs" / "design" / "run-workspace"
    workspace = tmp_path / "mounted-workspace"
    structure = workspace / "outputs" / "model.cif"
    structure.parent.mkdir(parents=True)
    structure.write_text("data_workspace\n#\n")
    jobs.write_json(run / "command.json", {"command": ["docker", "run", "-v", f"{workspace}:/workspace:ro", "fixture:image"]})

    write_candidates(run, "fixture-generator", [_candidate("/workspace/outputs/model.cif")])
    persisted = json.loads((run / "artifacts/normalized_candidates/candidates.jsonl").read_text())
    stored = run / persisted["complex_pdb"]
    assert stored.is_file() and stored.read_text() == structure.read_text()
    assert persisted["raw_metadata"]["container_artifact_paths"]["complex_pdb"] == "/workspace/outputs/model.cif"
    assert persisted["raw_metadata"]["materialized_source_paths"]["complex_pdb"] == "/workspace/outputs/model.cif"


def test_boltz_context_output_uses_recorded_mount_mapping(tmp_path):
    run = tmp_path / "runs" / "refolding-validation" / "boltz-context"
    context = tmp_path / "boltz-context-mount"
    structure = context / "processed" / "prediction.cif"
    structure.parent.mkdir(parents=True)
    structure.write_text("data_prediction\n#\n")
    jobs.write_json(run / "command.json", {"command": [
        "docker", "run", "--mount", f"type=bind,source={context},target=/boltz-context", "fixture:image"
    ]})

    write_candidates(run, "fixture-predictor", [_candidate("/boltz-context/processed/prediction.cif")])
    persisted = json.loads((run / "artifacts/normalized_candidates/candidates.jsonl").read_text())
    materialized = run / persisted["complex_pdb"]
    assert materialized.is_file() and materialized.read_text() == structure.read_text()
    assert persisted["raw_metadata"]["container_artifact_paths"]["complex_pdb"] == "/boltz-context/processed/prediction.cif"


def test_container_absolute_symlink_is_materialized_and_kept_as_provenance(tmp_path):
    run = tmp_path / "runs" / "design" / "container-link"
    canonical = run / "artifacts" / "raw" / "fixture" / "output" / "model.cif"
    canonical.parent.mkdir(parents=True)
    canonical.write_text("data_link_target\n#\n")
    alias = canonical.parent / "alias.cif"
    alias.symlink_to("/work/artifacts/raw/fixture/output/model.cif")

    write_candidates(run, "fixture-generator", [_candidate(str(alias))])
    persisted = json.loads((run / "artifacts/normalized_candidates/candidates.jsonl").read_text())
    assert persisted["complex_pdb"] == "artifacts/raw/fixture/output/model.cif"
    assert persisted["raw_metadata"]["container_artifact_provenance"]["complex_pdb"]["absolute_container_symlink_targets"] == [
        "/work/artifacts/raw/fixture/output/model.cif"
    ]


def test_managed_source_structure_uses_existing_cross_run_reference(tmp_path, monkeypatch):
    runs = tmp_path / "managed-runs"
    source = runs / "design" / "source"
    run = runs / "design" / "consumer"
    structure = source / "model.cif"
    structure.parent.mkdir(parents=True)
    structure.write_text("data_shared\n#\n")
    monkeypatch.setattr(portable_paths, "runs_root", lambda: runs)

    write_candidates(run, "fixture-predictor", [_candidate(str(structure))])
    persisted = json.loads((run / "artifacts/normalized_candidates/candidates.jsonl").read_text())
    assert persisted["complex_pdb"] == "runs:///design/source/model.cif"
    assert not (run / "artifacts/materialized_candidates").exists()
    assert read_candidates(run)[0]["complex_pdb"] == str(structure)


def test_missing_required_candidate_output_prevents_successful_finalization(tmp_path, monkeypatch):
    run = tmp_path / "runs" / "design" / "missing-required"
    write_candidates(run, "fixture-generator", [_candidate("/work/artifacts/raw/fixture/output/missing.cif")])

    _finish(run, monkeypatch)
    result = jobs.read_json(run / "result.json")
    metadata = jobs.read_json(run / "metadata.json")
    assert result["success"] is False
    assert "Required scientific output" in result["error"]
    assert metadata["status"] == "failed"
    assert metadata["artifact_collection_complete"] is False
    assert metadata["artifact_validation"]["required_missing"]


def test_direct_completed_status_is_also_gated_by_artifact_collection(tmp_path):
    run = tmp_path / "runs" / "design" / "direct-status"
    write_candidates(run, "fixture-generator", [_candidate("/work/artifacts/raw/fixture/output/missing.cif")])
    jobs.write_json(run / "metadata.json", {"run_id": run.name, "status": "running", "task_group": "design"})
    jobs.write_json(run / "result.json", {"success": True})

    jobs.update_status(run, "completed")
    metadata = jobs.read_json(run / "metadata.json")
    result = jobs.read_json(run / "result.json")
    assert metadata["status"] == "failed"
    assert metadata["artifact_collection_complete"] is False
    assert result["success"] is False


def test_manifest_distinguishes_missing_required_and_optional_execution_artifacts(tmp_path):
    run = tmp_path / "runs" / "structure-jobs" / "manifest-case"
    run.mkdir(parents=True)
    (run / "artifacts.json").write_text(json.dumps({
        "kind": "artifact_manifest", "schema_version": 1, "run_id": run.name,
        "artifacts": [
            {"run_id": run.name, "artifact_type": "predicted_complex", "path": "artifacts/required.cif"},
            {"run_id": run.name, "artifact_type": "temporary_cache", "path": "artifacts/cache.tmp", "required": False},
            {"run_id": run.name, "artifact_type": "temporary_cache", "path": "artifacts/optional-link", "required": False},
        ],
    }))
    dangling = run / "artifacts" / "optional-link"
    dangling.parent.mkdir()
    dangling.symlink_to("/workspace/.cache/unavailable")

    validation = validate_run_artifacts(run)
    assert validation["status"] == "incomplete"
    assert [item["path"] for item in validation["required_missing"]] == ["artifacts/required.cif"]
    assert validation["dangling_symlinks"][0]["classification"] == "optional_execution_artifact"


def test_undeclared_dangling_structure_output_blocks_collection(tmp_path):
    run = tmp_path / "runs" / "design" / "unclassified-output"
    output = run / "artifacts" / "raw" / "engine" / "output" / "model.pdb"
    output.parent.mkdir(parents=True)
    output.symlink_to("/work/artifacts/raw/engine/output/model.pdb")

    validation = validate_run_artifacts(run)
    assert validation["status"] == "incomplete"
    assert validation["dangling_symlinks"][0]["classification"] == "unresolved_scientific_output"


def test_backup_style_copy_preserves_historical_dangling_link_and_checks_canonical_output(tmp_path, monkeypatch):
    run = tmp_path / "runs" / "design" / "historical-shaped"
    canonical = run / "artifacts" / "normalized_candidates" / "model.cif"
    canonical.parent.mkdir(parents=True)
    canonical.write_text("data_canonical\n#\n")
    write_candidates(run, "fixture-generator", [_candidate(str(canonical))])
    dangling = run / "artifacts" / "raw" / "old-engine" / "output" / "legacy.cif"
    dangling.parent.mkdir(parents=True)
    dangling.symlink_to("/workspace/generated/removed/legacy.cif")

    # Backup/restore code must copy the link entry itself, never dereference it.
    backup = tmp_path / "isolated-restore" / "run"
    backup.parent.mkdir(parents=True)
    shutil.copytree(run, backup, symlinks=True)
    restored_link = backup / "artifacts/raw/old-engine/output/legacy.cif"
    assert restored_link.is_symlink()
    assert restored_link.readlink() == Path("/workspace/generated/removed/legacy.cif")

    validation = validate_run_artifacts(backup, strict_dangling_links=False)
    finalization_validation = validate_run_artifacts(backup, strict_dangling_links=True)
    restored_canonical = backup / "artifacts/normalized_candidates/model.cif"
    assert validation["status"] == "complete"
    assert finalization_validation["status"] == "incomplete"
    assert validation["symlinks_followed"] is False
    assert validation["dangling_symlinks"] == [{
        "path": "artifacts/raw/old-engine/output/legacy.cif",
        "target": "/workspace/generated/removed/legacy.cif",
        "classification": "unresolved_scientific_output",
    }]
    assert hashlib.sha256(restored_canonical.read_bytes()).hexdigest() == hashlib.sha256(canonical.read_bytes()).hexdigest()
    assert read_candidates(backup)[0]["candidate_id"] == "design-1"


@pytest.mark.parametrize("suffix", [".cif", ".pdb", ".sdf"])
@pytest.mark.parametrize("container_target", [
    "/tmp/removed/output", "/workspace/processed", "/boltz-context/processed",
    "/work/tmp/output", "/workspace/outputs",
])
def test_scientific_output_precedes_temporary_origin_heuristics(tmp_path, monkeypatch, suffix, container_target):
    run = tmp_path / "runs" / "design" / "scientific-output"
    output = run / "artifacts/raw/engine/output" / ("model" + suffix)
    output.parent.mkdir(parents=True)
    output.symlink_to(container_target + "/model" + suffix)
    validation = validate_run_artifacts(run)
    assert validation["status"] == "incomplete"
    assert validation["dangling_symlinks"][0]["classification"] == "unresolved_scientific_output"
    _finish(run, monkeypatch)
    assert jobs.read_json(run / "metadata.json")["artifact_collection_complete"] is False
    assert jobs.read_json(run / "metadata.json")["status"] == "failed"
    # The other normal finalization entry point must enforce the same gate.
    jobs.update_status(run, "completed")
    assert jobs.read_json(run / "metadata.json")["status"] == "failed"
    assert output.is_symlink()  # Validation never repairs/deletes historical entries.


def test_output_inside_a_cache_named_directory_is_still_scientific(tmp_path):
    run = tmp_path / "design/cache-contained-output"
    output = run / "artifacts/raw/engine/cache/output/model.cif"
    output.parent.mkdir(parents=True)
    output.symlink_to("/workspace/cache/model.cif")
    assert validate_run_artifacts(run)["status"] == "incomplete"


@pytest.mark.parametrize("relative", ["artifacts/processed/prediction.cif", "artifacts/cache/output/confidence.npz"])
def test_scientific_payload_cannot_be_hidden_by_its_own_cache_directory(tmp_path, relative):
    run = tmp_path / "design/cache-payload"
    output = run / relative
    output.parent.mkdir(parents=True)
    output.symlink_to("/workspace/cache/payload")
    assert validate_run_artifacts(run)["status"] == "incomplete"


def test_target_cache_path_alone_cannot_declare_an_unknown_link_optional(tmp_path):
    run = tmp_path / "design/unknown-link"
    run.mkdir(parents=True)
    (run / "opaque-artifact").symlink_to("/tmp/.cache/unknown-payload")
    validation = validate_run_artifacts(run)
    assert validation["status"] == "incomplete"
    assert validation["dangling_symlinks"][0]["classification"] == "unclassified_dangling_link"


def test_explicit_optional_alias_preserves_the_canonical_scientific_output(tmp_path):
    run = tmp_path / "design/declared-alias"
    canonical = run / "artifacts/model.cif"
    alias = run / "artifacts/raw/engine/output/old-model.cif"
    alias.parent.mkdir(parents=True)
    canonical.write_text("data_retained\n")
    alias.symlink_to("/tmp/no-longer-needed/model.cif")
    jobs.write_json(run / "artifacts.json", {
        "schema_version": 1, "run_id": run.name, "artifacts": [
            {"path": "artifacts/model.cif", "artifact_type": "predicted_complex",
             "sha256": hashlib.sha256(canonical.read_bytes()).hexdigest()},
            {"path": "artifacts/raw/engine/output/old-model.cif", "required": False},
        ],
    })
    validation = validate_run_artifacts(run)
    assert validation["status"] == "complete"
    assert validation["required_artifact_count"] == 1
    assert validation["dangling_symlinks"][0]["classification"] == "optional_execution_artifact"
    canonical.unlink()
    assert validate_run_artifacts(run)["status"] == "incomplete"


def test_container_alias_is_optional_only_when_its_exact_canonical_output_is_retained(tmp_path, monkeypatch):
    run = tmp_path / "design/retained-alias"
    canonical = run / "artifacts/raw/engine/output/model.cif"
    canonical.parent.mkdir(parents=True)
    canonical.write_text("data_retained\n")
    alias = canonical.with_name("alias.cif")
    alias.symlink_to("/work/artifacts/raw/engine/output/model.cif")
    write_candidates(run, "fixture-generator", [_candidate(str(alias))])
    validation = validate_run_artifacts(run)
    assert validation["status"] == "complete"
    assert validation["required_artifact_count"] == 1
    assert validation["dangling_symlinks"][0]["classification"] == "retained_scientific_alias"
    _finish(run, monkeypatch)
    assert jobs.read_json(run / "metadata.json")["status"] == "completed"
    canonical.unlink()
    assert validate_run_artifacts(run)["status"] == "incomplete"
