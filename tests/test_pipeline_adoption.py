import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from fastapi import HTTPException

from mn_molecular_design_studio import api
from mn_molecular_design_studio.engine.workflow_engine.workflows import pipeline_adoption


PROJECT = "project-a"
OWNER = "user-a"
BINDER = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMGWFRQAPGKGRELVAAISGSGGKTYYPDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAAGGSPHYDGYIPVEDVNFDYWGQGTQVTVSS"


def _write(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload), encoding="utf-8")


def _make_step(root, group, run_id, *, tool, job_type, params=None, inputs=None, candidate=None, stamp="2026-10-23T10:00:00+00:00"):
    run_dir = root / group / run_id
    _write(run_dir / "metadata.json", {
        "run_id": run_id, "task_group": group, "tool": tool, "job_type": job_type,
        "status": "completed", "created_at": stamp, "worker_started_at": stamp,
        "completed_at": stamp, "updated_at": stamp, "job_code": run_id[:5].upper(),
        "studio_owner_id": OWNER, "project_id": PROJECT, "campaign_name": "Original job title",
        "preserve_this": "yes",
    })
    _write(run_dir / "input.json", {"tool": tool, "job_type": job_type, "params": params or {}, "inputs": inputs or {}})
    _write(run_dir / "result.json", {"success": True, "metrics": {"candidate_count": 1}})
    if candidate:
        candidate_path = run_dir / "artifacts" / "normalized_candidates" / "candidates.jsonl"
        candidate_path.parent.mkdir(parents=True, exist_ok=True)
        candidate_path.write_text(json.dumps(candidate) + "\n", encoding="utf-8")
    if tool == "optimized_prediction_boltz2":
        _write(run_dir / "artifacts/raw/optimized_prediction/input_mapping.json", {
            "items": [{"raw_yaml": False, "input_details": {"template_count": 0}}],
        })
    return run_dir


def _adoption_fixture(tmp_path, monkeypatch):
    runs = tmp_path / "runs"
    gen = _make_step(runs, "design", "gen-1", tool="boltzgen_generator", job_type="design_generation", candidate={"candidate_id": "boltzgen_00001"})
    ab = _make_step(
        runs, "design", "ab-1", tool="ligandmpnn", job_type="sequence_design",
        params={"model_type": "ab_mpnn", "sabr_mode": "sabr", "source_run_dir": str(gen), "selected_candidate_ids": ["boltzgen_00001"], "num_seq_per_target": 1, "seed": 271828},
        candidate={"candidate_id": "abmpnn_00001", "binder_sequence": BINDER},
    )
    refold = _make_step(
        runs, "refolding-validation", "fold-1", tool="optimized_prediction_boltz2", job_type="complex_refolding",
        params={"engine": "boltz2", "samples": 1, "num_recycles": 3, "engine_options": {"use_msa": True, "sampling_steps": 200, "seed_count": 1}},
        inputs={"tasks": [{"has_target_structure": False, "restraints": [], "protein_chains": [{"id": "A", "role": "target", "sequence": "ACDE"}, {"id": "B", "role": "binder", "sequence": BINDER}]}]},
    )
    monkeypatch.setattr(pipeline_adoption, "runs_root", lambda: runs)
    monkeypatch.setattr("mn_molecular_design_studio.workspace.mark_studio_job", lambda *a, **k: None)
    serial = {"next": 1}

    def create_job(task_group, job_type, tool, inputs, params):
        run_id = f"summary-{serial['next']}"
        serial["next"] += 1
        run_dir = runs / task_group / run_id
        _write(run_dir / "metadata.json", {"run_id": run_id, "task_group": task_group, "job_type": job_type, "tool": tool, "status": "queued", "job_code": "SUMRY", "created_at": "2026-10-23T11:00:00+00:00"})
        _write(run_dir / "input.json", {"job_type": job_type, "tool": tool, "inputs": inputs, "params": params})
        (run_dir / "command.json").write_text("{}", encoding="utf-8")
        return SimpleNamespace(run_dir=run_dir, run_id=run_id)

    def finish_job(run_dir, success, result):
        _write(run_dir / "result.json", {"success": success, **result})
        metadata_path = run_dir / "metadata.json"
        metadata = json.loads(metadata_path.read_text())
        metadata.update(status="completed", completed_at="2026-10-23T11:00:00+00:00", updated_at="2026-10-23T11:00:00+00:00")
        _write(metadata_path, metadata)

    monkeypatch.setattr(pipeline_adoption, "create_job", create_job)
    monkeypatch.setattr(pipeline_adoption, "finish_job", finish_job)
    monkeypatch.setattr("mn_molecular_design_studio.job_catalog.reconcile", lambda: None)
    steps = [
        {"task_group": "design", "run_id": "gen-1", "stage_label": "Vanilla BoltzGen generation", "selected_candidate_ids": ["boltzgen_00001"]},
        {"task_group": "design", "run_id": "ab-1", "stage_label": "SAbR + AbMPNN sequence design"},
        {"task_group": "refolding-validation", "run_id": "fold-1", "stage_label": "Boltz-2 target-MSA refold"},
    ]
    rows = {f"{s['task_group']}/{s['run_id']}": {"job_code": s["run_id"][:5].upper()} for s in steps}
    group = {"name": "Controlled comparison · Vanilla BoltzGen → SAbR/AbMPNN → Boltz-2", "steps": steps}
    return runs, [gen, ab, refold], [group], rows


def test_adoption_creates_three_completed_display_steps_without_execution(tmp_path, monkeypatch):
    runs, child_dirs, groups, rows = _adoption_fixture(tmp_path, monkeypatch)
    adopted = pipeline_adoption.adopt_completed_pipelines(groups, owner_id=OWNER, project_id=PROJECT, resolved_rows=rows)
    root = Path(adopted[0]["pipeline_run_dir"])
    execution = json.loads((root / "graph_execution.json").read_text())
    graph = json.loads((root / "graph_snapshot.json").read_text())
    assert adopted[0]["step_run_refs"] == ["design/gen-1", "design/ab-1", "refolding-validation/fold-1"]
    assert len(execution["nodes"]) == 3
    assert all(node["status"] == "completed" for node in execution["nodes"].values())
    assert json.loads((root / "result.json").read_text())["metrics"]["completed_steps"] == 3
    boltz_params = graph["nodes"][2]["data"]["params"]
    assert boltz_params["samples"] == boltz_params["diffusion_samples"] == 1
    assert boltz_params["num_recycles"] == boltz_params["recycling_steps"] == 3
    assert boltz_params["sampling_steps"] == 200
    assert boltz_params["use_target_template"] is False
    assert not (root / "worker_request.json").exists()
    for child in child_dirs:
        metadata = json.loads((child / "metadata.json").read_text())
        assert metadata["campaign_id"] == root.name
        assert metadata["preserve_this"] == "yes"
        assert metadata["status"] == "completed"
    root_meta = json.loads((root / "metadata.json").read_text())
    assert root_meta["worker_started_at"] == "2026-10-23T10:00:00+00:00"
    assert root_meta["completed_at"] == "2026-10-23T10:00:00+00:00"
    assert root_meta["adopted_at"]
    assert root_meta["adopted_at"] != root_meta["completed_at"]


def test_adoption_is_idempotent_and_does_not_steal_membership(tmp_path, monkeypatch):
    runs, _, groups, rows = _adoption_fixture(tmp_path, monkeypatch)
    first = pipeline_adoption.adopt_completed_pipelines(groups, owner_id=OWNER, project_id=PROJECT, resolved_rows=rows)
    second = pipeline_adoption.adopt_completed_pipelines(groups, owner_id=OWNER, project_id=PROJECT, resolved_rows=rows)
    assert second[0]["pipeline_run_id"] == first[0]["pipeline_run_id"]
    assert len(list((runs / "campaign").iterdir())) == 1
    different_name = [{**groups[0], "name": "A different display group"}]
    with pytest.raises(ValueError, match="already belongs"):
        pipeline_adoption.adopt_completed_pipelines(different_name, owner_id=OWNER, project_id=PROJECT, resolved_rows=rows)
    assert len(list((runs / "campaign").iterdir())) == 1


def test_undo_restores_original_child_metadata_and_removes_summary(tmp_path, monkeypatch):
    runs, child_dirs, groups, rows = _adoption_fixture(tmp_path, monkeypatch)
    original = [json.loads((child / "metadata.json").read_text()) for child in child_dirs]
    adopted = pipeline_adoption.adopt_completed_pipelines(groups, owner_id=OWNER, project_id=PROJECT, resolved_rows=rows)
    root = Path(adopted[0]["pipeline_run_dir"])
    result = pipeline_adoption.undo_pipeline_adoption(root, owner_id=OWNER, project_id=PROJECT)
    assert result["status"] == "undone"
    assert not root.exists()
    for child, expected in zip(child_dirs, original):
        assert json.loads((child / "metadata.json").read_text()) == expected


@pytest.mark.parametrize("mutation", ["source", "sequence", "project"])
def test_api_adoption_rejects_wrong_source_sequence_or_project(tmp_path, monkeypatch, mutation):
    runs = tmp_path / "runs"
    gen = _make_step(runs, "design", "gen-1", tool="boltzgen_generator", job_type="design_generation", candidate={"candidate_id": "boltzgen_00001"})
    wrong_source = _make_step(runs, "design", "other", tool="boltzgen_generator", job_type="design_generation", candidate={"candidate_id": "boltzgen_00001"})
    ab = _make_step(runs, "design", "ab-1", tool="ligandmpnn", job_type="sequence_design",
                    params={"model_type": "ab_mpnn", "sabr_mode": "sabr", "source_run_dir": str(gen), "selected_candidate_ids": ["boltzgen_00001"]},
                    candidate={"candidate_id": "abmpnn_00001", "binder_sequence": BINDER})
    if mutation == "source":
        ab_params = json.loads((ab / "input.json").read_text())["params"]
        ab_params["source_run_dir"] = str(wrong_source)
        _write(ab / "input.json", {"tool": "ligandmpnn", "job_type": "sequence_design", "params": ab_params, "inputs": {}})
    different = "AAAA" if mutation == "sequence" else BINDER
    refold = _make_step(runs, "refolding-validation", "fold-1", tool="optimized_prediction_boltz2", job_type="complex_refolding",
                        params={"engine": "boltz2", "engine_options": {"use_msa": True}},
                        inputs={"tasks": [{"has_target_structure": False, "restraints": [], "protein_chains": [{"id": "A", "role": "target", "sequence": "ACDE"}, {"id": "B", "role": "binder", "sequence": different}]}]})
    ref = {
        "design/gen-1": {"run_dir": str(gen), "status": "completed", "success": True, "project_id": PROJECT},
        "design/ab-1": {"run_dir": str(ab), "status": "completed", "success": True, "project_id": PROJECT},
        "refolding-validation/fold-1": {"run_dir": str(refold), "status": "completed", "success": True, "project_id": "other-project" if mutation == "project" else PROJECT},
    }
    monkeypatch.setattr(api, "_request_payload_and_project", lambda payload, user: (payload, PROJECT))
    monkeypatch.setattr(api, "_resolve_api_job", lambda reference, user: ref[reference])
    monkeypatch.setattr(
        "mn_molecular_design_studio.job_catalog.get_job_access_record",
        lambda group, run: {"owner_id": OWNER, "project_id": "other-project" if mutation == "project" and run == "fold-1" else PROJECT},
    )
    monkeypatch.setattr("mn_molecular_design_studio.engine.workflow_engine.workflows.pipeline_adoption.adopt_completed_pipelines", lambda *a, **k: [])
    payload = {"project_id": PROJECT, "pipelines": [{"name": "Controlled test", "steps": [
        {"task_group": "design", "run_id": "gen-1", "stage_label": "Vanilla BoltzGen"},
        {"task_group": "design", "run_id": "ab-1", "stage_label": "SAbR + AbMPNN"},
        {"task_group": "refolding-validation", "run_id": "fold-1", "stage_label": "Boltz-2"},
    ]}]}
    with pytest.raises(HTTPException) as error:
        api.adopt_completed_pipeline_groups(payload, {"id": OWNER, "role": "user"})
    assert error.value.status_code == 422


def test_api_accepts_completed_actual_input_shape_and_rows_without_project_field(tmp_path, monkeypatch):
    runs = tmp_path / "runs"
    gen = _make_step(runs, "design", "gen-1", tool="rfantibody", job_type="design_generation", candidate={"candidate_id": "rfantibody_00001"})
    ab = _make_step(runs, "design", "ab-1", tool="ligandmpnn", job_type="sequence_design",
                    params={"model_type": "ab_mpnn", "sabr_mode": "sabr", "source_run_dir": str(gen), "selected_candidate_ids": ["rfantibody_00001"]},
                    candidate={"candidate_id": "abmpnn_00001", "binder_sequence": BINDER})
    fold = _make_step(runs, "refolding-validation", "fold-1", tool="optimized_prediction_boltz2", job_type="complex_refolding",
                      params={"engine": "boltz2", "engine_options": {"use_msa": True}},
                      inputs={"tasks": [{"has_target_structure": False, "restraints": [], "protein_chains": [{"id": "A", "role": "target", "sequence": "ACDE"}, {"id": "B", "role": "binder", "sequence": BINDER}]}]})
    rows = {
        "design/gen-1": {"run_dir": str(gen), "status": "completed", "success": True},
        "design/ab-1": {"run_dir": str(ab), "status": "completed", "success": True},
        "refolding-validation/fold-1": {"run_dir": str(fold), "status": "completed", "success": True},
    }
    monkeypatch.setattr(api, "_request_payload_and_project", lambda payload, user: (payload, PROJECT))
    monkeypatch.setattr(api, "_resolve_api_job", lambda reference, user: rows[reference])
    monkeypatch.setattr("mn_molecular_design_studio.job_catalog.get_job_access_record", lambda group, run: {"owner_id": OWNER, "project_id": PROJECT})
    captured = {}
    monkeypatch.setattr("mn_molecular_design_studio.engine.workflow_engine.workflows.pipeline_adoption.adopt_completed_pipelines", lambda groups, **kwargs: captured.update(groups=groups, kwargs=kwargs) or [{"name": groups[0]["name"]}])
    payload = {"project_id": PROJECT, "pipelines": [{"name": "Controlled test", "steps": [
        {"task_group": "design", "run_id": "gen-1", "stage_label": "MIMOSA RFantibody"},
        {"task_group": "design", "run_id": "ab-1", "stage_label": "SAbR + AbMPNN"},
        {"task_group": "refolding-validation", "run_id": "fold-1", "stage_label": "Boltz-2"},
    ]}]}
    result = api.adopt_completed_pipeline_groups(payload, {"id": OWNER, "role": "user"})
    assert result == [{"name": "Controlled test"}]
    assert captured["groups"][0]["steps"][0]["selected_candidate_ids"] == ["rfantibody_00001"]
