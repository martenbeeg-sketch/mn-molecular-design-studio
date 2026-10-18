from __future__ import annotations

import json
from pathlib import Path

from mn_molecular_design_studio.engine.workflow_engine.core import local_worker


def test_sequence_pipeline_forwards_validation_options_without_duplicate_enqueue_args(monkeypatch, tmp_path):
    source_run = tmp_path / "source"
    source_candidates = source_run / "artifacts" / "normalized_candidates" / "candidates.jsonl"
    source_candidates.parent.mkdir(parents=True)
    source_candidates.write_text(json.dumps({"candidate_id": "backbone-1", "stage": "generation.backbone"}) + "\n")

    pipeline_run = tmp_path / "pipeline"
    pipeline_candidates = pipeline_run / "artifacts" / "normalized_candidates" / "candidates.jsonl"
    pipeline_candidates.parent.mkdir(parents=True)
    pipeline_candidates.write_text(json.dumps({"candidate_id": "sequence-1", "stage": "sequence_design"}) + "\n")
    (pipeline_run / "metadata.json").write_text("{}")
    (pipeline_run / "input.json").write_text(json.dumps({"job_type": "sequence_design", "tool": "ligandmpnn"}))
    (pipeline_run / "result.json").write_text(json.dumps({"success": True, "metrics": {"return_code": 0, "candidate_count": 1}}))
    (pipeline_run / "molecular_design_studio.json").write_text(json.dumps({
        "project_id": "test-project", "studio_owner_id": "test-owner",
    }))
    (pipeline_run / "worker_request.json").write_text(json.dumps({
        "kind": "sequence_design_pipeline",
        "kwargs": {
            "backend": "shared_ligandmpnn",
            "source_run_dir": str(source_run),
            "candidates_jsonl": str(source_candidates),
            "design_kwargs": {},
            "validation_kwargs": {
                "evaluation_name": "custom validation",
                "job_type": "sequence_design_validation",
                "tool_name": "sequence_design_validation",
                "selected_candidate_ids": ["sequence-1"],
                "gpu_device": "0",
                "run_boltz2_initial_guess": True,
            },
            "resume": True,
        },
    }))

    from mn_molecular_design_studio.engine.workflow_engine.workflows import benchmark

    captured = {}
    validation_run = tmp_path / "validation"

    def enqueue(**kwargs):
        captured.update(kwargs)
        return validation_run

    spawned = []
    marked = []
    from mn_molecular_design_studio import workspace

    monkeypatch.setattr(benchmark, "enqueue_candidate_refolding_evaluation", enqueue)
    monkeypatch.setattr(local_worker, "spawn_worker_for_run", lambda path: spawned.append(Path(path)))
    monkeypatch.setattr(workspace, "mark_studio_job", lambda result, **kwargs: marked.append((result, kwargs)))

    assert local_worker.run_worker_job(pipeline_run) == 0
    assert captured["evaluation_name"] == "custom validation"
    assert captured["job_type"] == "sequence_design_validation"
    assert captured["tool_name"] == "sequence_design_validation"
    assert captured["selected_candidate_ids"] == ["sequence-1"]
    assert captured["run_boltz2_initial_guess"] is True
    assert captured["resume"] is True
    assert spawned == [validation_run]
    metadata = json.loads((pipeline_run / "metadata.json").read_text())
    assert metadata["validation_run_dir"] == str(validation_run)
    assert marked == [(
        {"job": {"run_dir": str(validation_run)}},
        {"project_id": "test-project", "owner_id": "test-owner"},
    )]
    result = json.loads((pipeline_run / "result.json").read_text())
    assert result["success"] is True
    metadata = json.loads((pipeline_run / "metadata.json").read_text())
    assert metadata["status"] == "completed"
