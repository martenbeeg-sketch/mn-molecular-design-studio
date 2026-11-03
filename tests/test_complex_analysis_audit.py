"""Independent regressions using native artifact contracts and negative controls."""
import csv
import json


def test_peppatch_hydrophobic_evidence_excludes_electrostatics_and_keeps_insertions(tmp_path):
    from mn_molecular_design_studio.complex_analysis_liabilities import _pep_patch_findings

    artifact = tmp_path / "residue_patches.csv"
    with artifact.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["method", "chain_id", "residue_id", "patch_id", "value"])
        writer.writeheader()
        writer.writerows([
            {"method": "electrostatic", "chain_id": "A", "residue_id": "10", "patch_id": 1, "value": 3.2},
            {"method": "hydrophobic", "chain_id": "A", "residue_id": "10B", "patch_id": 2, "value": .9},
        ])
    track = {"A": [
        {"chain_id": "A", "author_number": 10, "insertion_code": code,
         "selector": "A:10^" + (code or "_"), "sequence_index": i}
        for i, code in enumerate(["", "B"])
    ]}
    findings, _ = _pep_patch_findings(
        [{"evaluator": "pep_patch", "evaluation_id": "native-eval", "artifacts": [
            {"name": artifact.name, "path": str(artifact)}]}],
        ["A"], track, interface_indices={},
    )
    assert len(findings) == 1
    assert findings[0]["selector"] == "A:10^B"
    assert findings[0]["evidence"]["patch_id"] == 2


def test_pae_lookup_does_not_assign_another_same_length_models_sidecar(tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.workflows.analysis import _find_pae_path
    selected = tmp_path / "model_1.pdb"
    selected.write_text("END\n")
    (tmp_path / "model_10.pdb").write_text("END\n")
    sidecar = tmp_path / "model_10_pae.json"
    sidecar.write_text(json.dumps({"pae": [[0, 2], [3, 0]]}))
    assert _find_pae_path(tmp_path, {}, selected) is None
    assert _find_pae_path(tmp_path, {"raw_metadata": {"pae_path": str(sidecar)}}, selected) == sidecar
    # Missing explicit evidence must remain missing rather than falling back.
    assert _find_pae_path(tmp_path, {"raw_metadata": {"pae_path": "missing.json"}}, selected) is None


def test_database_pae_sidecar_requires_authorized_exact_source_structure(tmp_path, monkeypatch):
    from mn_molecular_design_studio import auth, job_catalog, runtime
    from mn_molecular_design_studio.complex_analysis_sources import database_pae_evidence
    from mn_molecular_design_studio.evaluations import file_digest
    run = tmp_path / "design" / "source-run"
    manifest = run / "artifacts/normalized_candidates/candidates.jsonl"
    manifest.parent.mkdir(parents=True)
    structure = run / "exact-model.cif"
    structure.write_text("data_exact\n")
    sidecar = run / "exact-model.pae.json"
    sidecar.write_text(json.dumps({"pae": [[0, 2], [3, 0]]}))
    manifest.write_text(json.dumps({"candidate_id":"original", "complex_pdb":str(structure),
        "source_tool":"bindcraft2", "raw_metadata":{"pae_path":str(sidecar)}}) + "\n")
    monkeypatch.setattr(runtime, "runs_root", lambda: tmp_path)
    monkeypatch.setattr(job_catalog, "get_job_access_record", lambda *a: {"owner_id":"owner"})
    monkeypatch.setattr(auth, "can_access_job", lambda *a: True)
    source = {"structure_digest":file_digest(structure), "structure_context":{
        "task_group":"design", "run_id":"source-run", "candidate_id":"original"}}
    assert database_pae_evidence({"id":"owner"}, source)["sha256"] == file_digest(sidecar)
    assert database_pae_evidence({"id":"owner"}, {**source,"structure_digest":"different-model"}) is None
    monkeypatch.setattr(auth, "can_access_job", lambda *a: False)
    assert database_pae_evidence({"id":"other"}, source) is None
