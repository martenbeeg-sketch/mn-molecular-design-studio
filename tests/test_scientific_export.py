import hashlib
from io import BytesIO
import json
from zipfile import ZipFile

import pytest

from mn_molecular_design_studio.engine.workflow_engine.core.artifacts import build_run_zip


@pytest.fixture
def linked_runs(tmp_path, monkeypatch):
    from mn_molecular_design_studio import auth, job_catalog, runtime
    root = tmp_path / "runs"
    parent = root / "campaign/parent"
    child = root / "analysis/child"
    parent.mkdir(parents=True)
    child.mkdir(parents=True)
    structure = child / "artifacts/coordinates.pdb"
    structure.parent.mkdir()
    structure.write_bytes(b"REMARK original coordinates\n")
    digest = hashlib.sha256(structure.read_bytes()).hexdigest()
    (child / "artifacts.json").write_text(json.dumps({"schema_version": 1, "run_id": child.name,
        "artifacts": [{"path": "artifacts/coordinates.pdb", "sha256": digest, "size_bytes": structure.stat().st_size}]}))
    collection = parent / "artifacts/normalized_candidates/candidates.jsonl"
    collection.parent.mkdir(parents=True)
    collection.write_text("".join(json.dumps({"candidate_id": name, "complex_pdb": "runs:///analysis/child/artifacts/coordinates.pdb",
        "metrics": {"evidence/analysis/sequence_qc/qc_sequence_length": 128}}) + "\n" for name in ("one", "two")))
    (parent / "graph_execution.json").write_text('{"status":"completed"}')
    access = {("campaign", "parent"): {"owner_id": "owner"}, ("analysis", "child"): {"owner_id": "owner"}}
    registered = [{"relative_path": "analysis/child/artifacts/coordinates.pdb"}]
    monkeypatch.setattr(runtime, "runs_root", lambda: root)
    from mn_molecular_design_studio.engine.workflow_engine.core import portable_paths
    monkeypatch.setattr(portable_paths, "runs_root", lambda: root)
    monkeypatch.setattr(job_catalog, "get_job_access_record", lambda group, run_id: access.get((group, run_id)))
    monkeypatch.setattr(job_catalog, "job_artifacts", lambda group, run_id: registered)
    monkeypatch.setattr(auth, "can_access_job", lambda user, job: user["id"] == job.get("owner_id"))
    return parent, child, structure, collection, access, registered


def test_pipeline_bundle_retains_authorized_coordinates_and_measurements(linked_runs):
    parent, _, structure, collection, _, _ = linked_runs
    original = collection.read_bytes()
    with ZipFile(BytesIO(build_run_zip(parent, user={"id": "owner"}))) as archive:
        manifest = json.loads(archive.read("scientific_bundle/manifest.json"))
        rows = [json.loads(line) for line in archive.read("scientific_bundle/candidates.jsonl").splitlines()]
        assert manifest["candidate_count"] == 2
        assert len(manifest["structures"]) == 1
        assert rows[0]["complex_pdb"] == rows[1]["complex_pdb"]
        assert archive.read(rows[0]["complex_pdb"]) == structure.read_bytes()
        assert rows[0]["metrics"]["evidence/analysis/sequence_qc/qc_sequence_length"] == 128
        assert "graph_execution.json" in archive.namelist()
    assert collection.read_bytes() == original


@pytest.mark.parametrize("which", [("campaign", "parent"), ("analysis", "child")])
def test_both_parent_and_structure_producer_require_access(linked_runs, which):
    parent, _, _, _, access, _ = linked_runs
    access[which]["owner_id"] = "another-owner"
    with pytest.raises(PermissionError):
        build_run_zip(parent, user={"id": "owner"})


def test_unregistered_linked_structure_is_not_exported(linked_runs):
    parent, _, _, _, _, registered = linked_runs
    registered.clear()
    with pytest.raises(ValueError, match="registered artifact"):
        build_run_zip(parent, user={"id": "owner"})


def test_changed_linked_coordinates_fail_integrity_check(linked_runs):
    parent, _, structure, _, _, _ = linked_runs
    structure.write_bytes(b"REMARK modified coordinates\n")
    with pytest.raises(ValueError, match="digest mismatch"):
        build_run_zip(parent, user={"id": "owner"})


def test_missing_linked_coordinates_cannot_produce_incomplete_successful_bundle(linked_runs):
    parent, _, structure, _, _, _ = linked_runs
    structure.unlink()
    with pytest.raises(ValueError, match="unavailable"):
        build_run_zip(parent, user={"id": "owner"})


def test_run_archive_does_not_follow_outbound_or_dangling_symlinks(tmp_path):
    run = tmp_path / "run"
    artifacts = run / "artifacts"
    artifacts.mkdir(parents=True)
    secret = tmp_path / "secret"
    secret.write_bytes(b"not a job artifact")
    (artifacts / "outbound").symlink_to(secret)
    (artifacts / "directory").symlink_to(tmp_path, target_is_directory=True)
    (artifacts / "dangling").symlink_to("/work/missing")
    (artifacts / "safe.txt").write_text("retained")
    with ZipFile(BytesIO(build_run_zip(run))) as archive:
        assert set(archive.namelist()) == {"artifacts/safe.txt", "export_manifest.json"}
        assert len(json.loads(archive.read("export_manifest.json"))["skipped_symlinks"]) == 3


def test_structure_reference_cannot_escape_managed_results(linked_runs, tmp_path):
    parent, _, _, collection, _, _ = linked_runs
    secret = tmp_path / "secret.pdb"
    secret.write_bytes(b"not a managed artifact")
    collection.write_text(json.dumps({"candidate_id": "malicious", "complex_pdb": str(secret)}) + "\n")
    with pytest.raises(ValueError, match="managed Results"):
        build_run_zip(parent, user={"id": "owner"})


def test_top_level_artifacts_directory_alias_is_not_followed(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    outside = tmp_path / "private-directory"
    outside.mkdir()
    (outside / "private.txt").write_text("not an authorized job artifact")
    (run / "artifacts").symlink_to(outside, target_is_directory=True)
    with ZipFile(BytesIO(build_run_zip(run))) as archive:
        assert archive.namelist() == ["export_manifest.json"]
        assert json.loads(archive.read("export_manifest.json"))["skipped_symlinks"] == ["artifacts"]


def test_dangling_candidate_collection_is_not_reported_as_scientific_bundle(linked_runs):
    parent, _, _, collection, _, _ = linked_runs
    collection.unlink()
    collection.symlink_to(collection.parent / "missing.jsonl")
    with pytest.raises(ValueError, match="collection is unavailable"):
        build_run_zip(parent, user={"id": "owner"})


def test_no_collection_is_reported_as_local_archive_only(linked_runs):
    parent, _, _, collection, _, _ = linked_runs
    collection.unlink()
    with ZipFile(BytesIO(build_run_zip(parent, user={"id": "owner"}))) as archive:
        assert json.loads(archive.read("export_manifest.json"))["scientific_bundle"] is False
