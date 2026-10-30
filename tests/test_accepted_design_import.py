from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
import zipfile

from mn_molecular_design_studio.accepted_design_import import (
    _read_bundle,
    build_accepted_design_bundle,
    promote_accepted_design_bundle,
)


def _source_row(run_dir: Path, candidate_id: str, structure: Path):
    return {
        "candidate_id": candidate_id,
        "source_tool": "rfd4_proteina",
        "stage": "complex_refolding",
        "complex_pdb": str(structure),
        "binder_sequence": "ACDEFG",
        "target_chains": ["A"],
        "binder_chains": ["B"],
        "metrics": {"native_pass_filters": False, "generation_score": 0.7},
        "raw_metadata": {"source_run_dir": str(run_dir), "native_filter_reason": "optional PPI metric unavailable"},
    }


def _validation_run(tmp_path: Path):
    runs = tmp_path / "runs"
    source_dir = runs / "design" / "source-run"
    validation_dir = runs / "refolding-validation" / "validation-run"
    source_dir.mkdir(parents=True)
    validation_dir.mkdir(parents=True)
    accepted_design = source_dir / "accepted-design.cif"
    accepted_design.write_text("accepted native design")
    accepted_prediction = validation_dir / "accepted-refold.cif"
    accepted_prediction.write_text("accepted independent prediction")
    rejected_design = source_dir / "rejected-design.cif"
    rejected_design.write_text("rejected native design")
    rejected_prediction = validation_dir / "rejected-refold.cif"
    rejected_prediction.write_text("rejected independent prediction")
    accepted_source = _source_row(source_dir, "design-accepted", accepted_design)
    rejected_source = _source_row(source_dir, "design-rejected", rejected_design)
    rows = [
        {
            "candidate_id": "validation-accepted",
            "source_tool": "optimized_prediction_boltz2",
            "stage": "complex_refolding",
            "complex_pdb": str(accepted_prediction),
            "binder_sequence": "ACDEFG",
            "metrics": {"engine": "boltz2", "iptm": 0.82, "complex_plddt": 91.2},
            "raw_metadata": {
                "source_candidate": accepted_source,
                "source_run_dir": str(source_dir),
                "validation_acceptance": {"status": "accepted", "reasons": [], "ranking_value": 0.82},
            },
        },
        {
            "candidate_id": "validation-rejected",
            "source_tool": "optimized_prediction_boltz2",
            "stage": "complex_refolding",
            "complex_pdb": str(rejected_prediction),
            "binder_sequence": "ACDEFG",
            "metrics": {"engine": "boltz2", "iptm": 0.2},
            "raw_metadata": {
                "source_candidate": rejected_source,
                "source_run_dir": str(source_dir),
                "validation_acceptance": {"status": "rejected", "reasons": ["below threshold"]},
            },
        },
    ]
    artifact_dir = validation_dir / "artifacts/normalized_candidates"
    artifact_dir.mkdir(parents=True)
    (artifact_dir / "candidates.jsonl").write_text("".join(json.dumps(row) + "\n" for row in rows))
    return runs, validation_dir


def test_bundle_contains_only_accepted_candidate_and_its_design_and_validation_files(tmp_path, monkeypatch):
    runs, validation_dir = _validation_run(tmp_path)
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(runs))

    bundle, manifest = build_accepted_design_bundle(
        validation_dir, node_id="worker-1", validation_job="cluster:validation-job",
    )
    loaded_manifest, candidates, files = _read_bundle(bundle)

    assert manifest == loaded_manifest
    assert manifest["candidate_count"] == 1
    assert [row["candidate_id"] for row in candidates] == ["design-accepted"]
    assert candidates[0]["metrics"]["validation_iptm"] == 0.82
    assert candidates[0]["raw_metadata"]["validation_acceptance"]["status"] == "accepted"
    assert candidates[0]["raw_metadata"]["validation_provenance"]["validation_engine"] == "boltz2"
    assert len(files) == 2
    assert all(b"rejected" not in content for content in files.values())
    assert all(entry["path"].startswith("artifacts/accepted/") for entry in manifest["artifacts"])
    with zipfile.ZipFile(__import__("io").BytesIO(bundle)) as archive:
        assert all("rejected" not in name for name in archive.namelist())


def test_promotion_is_explicit_and_idempotent(tmp_path, monkeypatch):
    from mn_molecular_design_studio.collection_policy import set_run_policy
    from mn_molecular_design_studio.engine.workflow_engine.core import candidates, jobs
    from mn_molecular_design_studio import workspace

    runs, validation_dir = _validation_run(tmp_path)
    app_home = tmp_path / "app"
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(runs))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(app_home))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(tmp_path / "workspace"))
    bundle, _ = build_accepted_design_bundle(validation_dir, node_id="worker-1", validation_job="cluster:validation-job")

    created = []
    completed = []

    def fake_create_job(task_group, job_type, tool, inputs, params):
        run_id = f"import-{len(created) + 1}"
        run_dir = runs / task_group / run_id
        run_dir.mkdir(parents=True)
        (run_dir / "artifacts").mkdir()
        jobs.write_json(run_dir / "metadata.json", {"run_id": run_id, "task_group": task_group, "status": "queued"})
        created.append(run_dir)
        return SimpleNamespace(run_dir=run_dir, run_id=run_id)

    def fake_finish(run_dir, success, result):
        metadata = jobs.read_json(run_dir / "metadata.json")
        metadata["status"] = "completed" if success else "failed"
        jobs.write_json(run_dir / "metadata.json", metadata)
        completed.append(result)

    monkeypatch.setattr(jobs, "create_job", fake_create_job)
    monkeypatch.setattr(jobs, "finish_job", fake_finish)
    monkeypatch.setattr(workspace, "mark_studio_job", lambda *args, **kwargs: {})
    monkeypatch.setattr("mn_molecular_design_studio.collection_policy.set_run_policy", lambda *args, **kwargs: "automatic")

    first = promote_accepted_design_bundle(bundle, user={"id": "researcher"}, project_id="project-1")
    second = promote_accepted_design_bundle(bundle, user={"id": "researcher"}, project_id="project-1")
    imported = candidates.read_candidates(Path(first["run_dir"]))

    assert first["candidate_count"] == 1 and first["already_imported"] is False
    assert second["already_imported"] is True
    assert len(created) == 1 and len(completed) == 1
    assert imported[0]["candidate_id"] == "design-accepted"
    assert imported[0]["metrics"]["validation_iptm"] == 0.82
    assert Path(imported[0]["complex_pdb"]).read_text() == "accepted native design"
    assert Path(imported[0]["raw_metadata"]["validation_complex_pdb"]).read_text() == "accepted independent prediction"


def test_local_results_ui_requires_researcher_click_before_accepted_database_promotion(tmp_path, monkeypatch):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import accepted_design_import

    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    run_dir = tmp_path / "runs/refolding-validation/accepted-run"
    run_dir.mkdir(parents=True)
    (run_dir / "metadata.json").write_text(json.dumps({"database_collection_policy": "manual"}))
    (run_dir / "result.json").write_text(json.dumps({
        "success": True, "metrics": {"validation_acceptance": {"accepted": 2, "rejected": 3}},
    }))
    (run_dir / "molecular_design_studio.json").write_text(json.dumps({"project_id": "project-1"}))
    calls = []
    monkeypatch.setattr(
        accepted_design_import, "build_accepted_design_bundle",
        lambda *args, **kwargs: calls.append((args, kwargs)) or (b"accepted-only", {"candidate_count": 2}),
    )
    monkeypatch.setattr(
        accepted_design_import, "promote_accepted_design_bundle",
        lambda *args, **kwargs: {"candidate_count": 2, "already_imported": False},
    )
    app = AppTest.from_string(
        "from pathlib import Path\n"
        "from mn_molecular_design_studio.molecular_database_results_ui import _render_accepted_design_promotion\n"
        f"_render_accepted_design_promotion(Path({str(run_dir)!r}), {{'id':'researcher'}})\n"
    ).run()

    assert not app.exception
    assert len(app.button) == 1
    assert calls == []
    app.button[0].click().run()
    assert not app.exception
    assert len(calls) == 1
    assert any("Saved 2 accepted design(s)" in item.value for item in app.success)
