import json
from pathlib import Path

from streamlit.testing.v1 import AppTest


def test_partial_pipeline_results_defers_details_and_ignores_empty_evidence(tmp_path, monkeypatch):
    from mn_molecular_design_studio import logical_runs, molecular_database_results_ui
    from mn_molecular_design_studio.auth import development_user

    home = tmp_path / "studio"
    home.mkdir()
    runs = tmp_path / "runs"
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(home))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(home))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(runs))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_AUTH_MODE", "development")
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(home / "cluster"))
    (home / "settings.json").write_text(json.dumps({
        "results_directory": str(runs), "reference_directory": str(home / "reference_files"),
    }))
    run = runs / "campaign" / "partial-pipeline"
    output = run / "artifacts" / "normalized_candidates"
    output.mkdir(parents=True)
    (output / "final_evidence.csv").write_text("\n")
    (run / "stderr.log").write_text("Fixture worker log")
    (run / "artifacts" / "model.pdb").write_text("END\n")
    (run / "metadata.json").write_text(json.dumps({
        "task_group": "campaign", "run_id": run.name, "tool": "visual_pipeline", "status": "failed",
        "job_type": "visual_pipeline", "campaign_name": "Fixture pipeline",
    }))
    (run / "molecular_design_studio.json").write_text(json.dumps({
        "project_id": "fixture", "studio_owner_id": development_user()["id"],
    }))
    (run / "input.json").write_text(json.dumps({"params": {"batch_size": 50}}))
    (run / "result.json").write_text(json.dumps({
        "success": False, "outputs": {"raw_trace": "x" * 600_000},
    }))
    child = runs / "analysis" / "native-step"
    child.mkdir(parents=True)
    (child / "metadata.json").write_text(json.dumps({
        "task_group": "analysis", "run_id": child.name, "tool": "fixture_tool", "status": "completed",
        "campaign_name": "Fixture tool output", "hidden": True,
    }))
    (child / "molecular_design_studio.json").write_text((run / "molecular_design_studio.json").read_text())
    (run / "artifacts/collected.jsonl").write_text(
        json.dumps({"candidate_id": "saved-1", "source_tool": "fixture_tool", "metrics": {"plddt": 77}}) + "\n"
    )
    (run / "graph_snapshot.json").write_text(json.dumps({
        "nodes": [{"id": node, "data": {"label": label}} for node, label in
                  [("refold", "Batch step"), ("failed", "Failed check"), ("tool", "Local tool")]],
        "edges": [{"source": "tool", "target": "refold"}, {"source": "refold", "target": "failed"}],
    }))
    (run / "graph_execution.json").write_text(json.dumps({
        "nodes": {
            "failed": {"status": "failed", "error": "Worker configuration is unavailable",
                       "batches": [{"status": "submitting"}, {"status": "cancelled"}]},
            "refold": {"status": "completed", "candidate_count": 2,
                       "output": {"candidates_jsonl": "artifacts/collected.jsonl"},
                       "batches": [{"status": "completed", "job_run_dir": str(child), "candidate_count": 2},
                                   {"status": "completed", "cluster_job_id": "remote-batch"}]},
            "tool": {"status": "completed", "job_run_dir": str(child)},
        },
    }))

    def reject_eager_work(*args, **kwargs):
        raise AssertionError("The Results view loaded optional content before it was requested")

    monkeypatch.setattr(logical_runs, "catalog", reject_eager_work)
    monkeypatch.setattr(molecular_database_results_ui, "render_run_database_actions", reject_eager_work)
    app_path = Path(__file__).parents[1] / "mn_molecular_design_studio" / "run_app.py"
    app = AppTest.from_file(str(app_path), default_timeout=20)
    app.session_state["page"] = "Job details"
    app.session_state["selected-job-reference"] = f"campaign/{run.name}"
    app.run()

    assert not app.exception
    assert any("no published final candidates" in item.value for item in app.info)
    assert not any("No columns to parse" in item.value for item in app.warning)
    assert not any("Choose a file" in item.label or item.label in {"Log file", "Structure file"} for item in app.selectbox)
    assert not app.json
    assert not any(item.label == "Run result" for item in app.expander)
    assert any(item.proto.label == "Download run record (JSON)" for item in app.get("download_button"))
    assert not any(item.value == "Job files" for item in app.subheader)

    app.session_state[f"job-detail-tabs-{run.name}"] = "Files"
    app.run()
    assert not app.exception
    assert any(item.value == "Job files" for item in app.subheader)
    assert any("Choose a file" in item.label for item in app.selectbox)
    assert not any("no published final candidates" in item.value for item in app.info)

    app.session_state[f"job-detail-tabs-{run.name}"] = "Job settings"
    app.run()
    assert not app.exception
    assert not app.json
    settings = next(table.value for table in app.dataframe if "Field" in table.value.columns)
    assert settings.set_index("Field").loc["Input › Params › Batch size", "Value"] == "50"

    app.session_state[f"job-detail-tabs-{run.name}"] = "Results"
    step_key = f"pipeline-run-step-campaign/{run.name}"
    steps = next(table.value for table in app.dataframe if "Step" in table.value.columns)
    assert list(steps["Step"]) == ["Local tool", "Batch step", "Failed check"]
    assert list(steps["Batches"])[-1] == "2 not submitted"

    app.session_state[step_key] = "failed"
    app.run()
    assert not app.exception
    assert any("Worker configuration is unavailable" in item.value for item in app.error)
    batch_table = next(table.value for table in app.dataframe if "Batch" in table.value.columns)
    assert list(batch_table["Results"]) == ["", ""]

    app.session_state[step_key] = "refold"
    app.session_state[step_key + "-refold-collected"] = True
    app.run()
    assert not app.exception
    batch_table = next(table.value for table in app.dataframe if "Batch" in table.value.columns)
    assert list(batch_table["Results"]) == ["Open", "Open"]
    collected = next(table.value for table in app.dataframe if "Candidate" in table.value.columns)
    assert collected.iloc[0]["pLDDT"] == 77

    app.session_state[step_key] = "tool"
    app.session_state[step_key + "-open"] = "tool"
    app.run()
    assert not app.exception
    assert any(title.value == "Fixture tool output" for title in app.title)
    assert not app.json
    app.button(key="job-details-back").click().run()
    assert not app.exception
    assert app.session_state["selected-job-reference"] == f"campaign/{run.name}"
    assert any(item.label == "Open tool results" for item in app.button)


def test_result_details_preserve_values_and_bound_large_payloads():
    app = AppTest.from_string("""
import streamlit as st
from mn_molecular_design_studio.result_details_ui import render_result_details

render_result_details(st, {
    "metrics": {"plddt": 87.2, "passed": False, "failed_count": 0},
    "diagnostic": "oversized" * 100_000,
    "matrix": [[1, 2, 3]] * 600,
    "records": [{"sample_id": f"sample-{index}", "rank": index} for index in range(500)],
}, table_fields=("records",))
""").run()
    assert not app.exception
    assert not app.json
    fields = app.dataframe[0].value.set_index("Field")["Value"]
    assert fields["Metrics › pLDDT"] == "87.2"
    assert fields["Metrics › Passed"] == "No"
    assert fields["Metrics › Failed count"] == "0"
    assert len(fields["Diagnostic"]) <= 240
    assert fields["Matrix"] == "600 entries"
    assert len(app.dataframe[1].value) == 100
    assert any("100 of 500 entries" in item.value for item in app.caption)
