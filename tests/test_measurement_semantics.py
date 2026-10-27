"""Execution completion, usable evidence and acceptance are separate decisions."""
import copy
import json
from types import SimpleNamespace

import pytest

from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
from mn_molecular_design_studio.engine.workflow_engine.workflows import analysis, quality


@pytest.mark.parametrize("values, expected", [([.6, .8], "valid"), ([.6, None], "partial"), ([None, None], "failed")])
def test_analysis_batch_reports_measurement_outcome(tmp_path, monkeypatch, values, expected):
    captured = []
    monkeypatch.setattr(analysis, "finish_job", lambda run, success, result: captured.append((success, result)))
    rows = [{"candidate_id": str(i), "metrics": {"ipsae": value, **({"ipsae_error": "original PAE failure"} if value is None else {})}}
            for i, value in enumerate(values)]
    analysis._finish_metric_analysis(SimpleNamespace(run_dir=tmp_path), rows, tool="ipsae", table_name="metrics.csv", metric_names=["ipsae"])
    success, result = captured[0]
    assert result["execution_status"] == "completed"
    assert result["scientific_status"] == expected
    assert success is (expected != "failed")  # usable partial batches may continue through the scheduler
    outputs = read_candidates(tmp_path)
    assert len(outputs) == len(rows)
    assert outputs[0]["raw_metadata"]["source_candidate"]["metrics"] == rows[0]["metrics"]
    if values[-1] is None:
        assert outputs[-1]["metrics"]["ipsae_error"] == "original PAE failure"


def test_analysis_requires_requested_measurement_not_an_auxiliary_score(tmp_path, monkeypatch):
    captured = []
    monkeypatch.setattr(analysis, "finish_job", lambda run, success, result: captured.append((success, result)))
    analysis._finish_metric_analysis(SimpleNamespace(run_dir=tmp_path),
        [{"candidate_id": "no-ipsae", "metrics": {"ipsae": float("nan"), "ipae": 3.}}],
        tool="ipsae", table_name="metrics.csv", metric_names=["ipsae", "ipae"])
    assert captured[0][0] is False


def test_missing_structure_is_retained_for_per_candidate_failure(design_source_factory):
    source, row = design_source_factory()
    path = source / "artifacts/normalized_candidates/candidates.jsonl"
    path.write_text(json.dumps(row) + "\n" + json.dumps({**row, "candidate_id": "missing", "complex_pdb": None}) + "\n")
    assert len(analysis._analysis_source_candidates(source, path)) == 2


@pytest.mark.parametrize("metrics", [{}, {"qc_sequence_length": float("inf")},
    {"qc_sequence_length": 5, "qc_sequence_pi": float("nan"), "qc_sequence_gravy": .1}])
def test_completed_qc_without_valid_required_metrics_is_not_scientific_success(tmp_path, monkeypatch, design_source_factory, metrics):
    source, _ = design_source_factory()
    run = tmp_path / "qc"; run.mkdir()
    monkeypatch.setattr(quality, "_analysis_job", lambda *a, **k: SimpleNamespace(run_dir=run))
    monkeypatch.setattr(quality, "apply_docker_cpu_limits_to_steps", lambda run, steps: steps)
    monkeypatch.setattr(quality, "docker_user_args", lambda: [])
    captured = []
    monkeypatch.setattr(quality, "finish_job", lambda run, success, result: captured.append((success, result)))
    def execute(command, **kwargs):
        if command[:3] == ["docker", "image", "inspect"]:
            return SimpleNamespace(returncode=0, stdout="sha256:test")
        (run / "artifacts/analysis/batch_status.json").write_text(json.dumps({"items": [
            {"id": "item_000000", "status": "completed", "metrics": metrics}]}))
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(quality.subprocess, "run", execute)
    quality.run_sequence_qc(source, source / "artifacts/normalized_candidates/candidates.jsonl")
    assert captured[0][0] is False
    assert captured[0][1]["execution_status"] == "completed"
    assert captured[0][1]["scientific_status"] == "failed"
    assert read_candidates(run)[0]["metrics"]["qc_error"]


def test_rerun_failure_cannot_use_stale_scores(tmp_path, monkeypatch):
    captured = []
    monkeypatch.setattr(analysis, "finish_job", lambda run, success, result: captured.append((success, result)))
    original = {"candidate_id": "old", "metrics": {"ipsae": .9, "ipsae_error": "new scorer failed"}}
    analysis._finish_metric_analysis(SimpleNamespace(run_dir=tmp_path), [copy.deepcopy(original)],
        tool="ipsae", table_name="metrics.csv", metric_names=["ipsae"])
    assert captured[0][0] is False
    assert read_candidates(tmp_path)[0]["metrics"]["analysis_status"] == "error"
