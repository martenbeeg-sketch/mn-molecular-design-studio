import ast
import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio.cluster import ui, config


@pytest.fixture
def isolated(tmp_path, monkeypatch):
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(tmp_path / "studio"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path / "studio"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(tmp_path / "runs"))
    return tmp_path


def campaign_rows():
    root = dict(cluster_job_id="parent", node_id="organizer", run_id="campaign-run", task_group="campaign",
                tool="visual_pipeline", campaign_name="Screening campaign", status="completed", project_id="p",
                owner_id="scientist", created_at="2026-10-25T00:00:00Z", local_reference="campaign/campaign-run",
                metrics={"round_count": 8, "unique_docked_compounds": 1590790})
    rows = [root]
    for iteration in range(1, 9):
        for batch in range(20):
            node = "organizer" if batch < 6 else "pc004032" if batch < 13 else "pc004033"
            rows.append(dict(root, cluster_job_id=f"batch-{iteration}-{batch}", run_id=f"batch-{iteration}-{batch}",
                tool="unidock_pro", campaign_id=root["run_id"], pipeline_parent_node_id="organizer",
                pipeline_iteration=iteration, node_id=node, node_name=node, metrics={},
                campaign_name=f"Iteration {iteration} docking", local_reference=f"ligand-docking/batch-{iteration}-{batch}"))
        if iteration < 8:
            rows.append(dict(root, cluster_job_id=f"model-{iteration}", run_id=f"model-{iteration}",
                tool="active_learning_model", campaign_id=root["run_id"], pipeline_parent_node_id="organizer",
                pipeline_iteration=iteration, metrics={}, campaign_name=f"Iteration {iteration} neural network"))
    return rows


def test_distributed_children_group_by_organizer_and_iteration():
    rows = campaign_rows()
    groups, mapped = ui._visual_pipeline_groups(rows)
    assert len(groups) == 1 and len(groups["parent"]) == 168
    assert len(mapped) == 168 and set(mapped.values()) == {"parent"}
    iterations = ui._pipeline_iterations(groups["parent"][1:])
    assert list(iterations) == list(range(1, 9))
    assert [len(jobs) for jobs in iterations.values()] == [21] * 7 + [20]
    assert ui._pipeline_step_count_label(rows[0], rows[1:]) == "8 iterations"
    active = dict(rows[0], status="running")
    assert ui._pipeline_step_count_label(active, rows[1:]) == "Iteration 8"
    assert ui._pipeline_iteration_status(active, 8, iterations[8], 8) == "In progress"


def test_grouping_does_not_cross_projects_owners_or_ambiguous_legacy_parents():
    rows = campaign_rows()[:2]
    root, child = rows
    for bad in (dict(child, project_id="other"), dict(child, owner_id="other"),
                dict(child, pipeline_parent_node_id="unknown")):
        assert ui._visual_pipeline_groups([root, bad])[0]["parent"] == [root]
    legacy = dict(child, node_id="third", pipeline_parent_node_id="")
    assert len(ui._visual_pipeline_groups([root, legacy])[0]["parent"]) == 2
    duplicate = dict(root, cluster_job_id="other-parent", node_id="other")
    assert all(len(group) == 1 for group in ui._visual_pipeline_groups([root, duplicate, legacy])[0].values())


def test_job_list_shows_campaign_then_iterations_then_opt_in_batches(isolated, monkeypatch):
    rows = campaign_rows()
    monkeypatch.setattr(ui, "catalog", lambda *a, **kw: ({"items": rows, "total": len(rows)}, False))
    from mn_molecular_design_studio import workspace
    monkeypatch.setattr(workspace, "load_organization", lambda: {"projects": [{"id": "p", "name": "Project"}]})
    app = AppTest.from_string('from mn_molecular_design_studio.cluster.ui import render_runs\n'
        'render_runs({"id":"scientist"}, "token", lambda **kw: {"items":[], "total":0})').run()
    assert not app.exception
    assert len(app.dataframe) == 1 and len(app.dataframe[0].value) == 1
    assert app.dataframe[0].value.iloc[0]["Steps"] == "8 iterations"
    assert app.dataframe[0].value.iloc[0]["Results"] == "Final results"
    app.session_state["cluster-results-pipeline-details-request"] = "parent"
    app.run()
    assert not app.exception
    assert len(app.dataframe) == 2 and len(app.dataframe[1].value) == 8
    assert list(app.dataframe[1].value["Iteration"]) == [f"Iteration {i}" for i in range(1, 9)]
    assert all(app.dataframe[1].value["Batches created"] == 20)
    app.session_state["cluster-results-pipeline-details-request"] = "parent"
    app.session_state["cluster-results-pipeline-iteration-request"] = 2
    app.run()
    assert not app.exception
    assert len(app.dataframe) == 3 and len(app.dataframe[2].value) == 21


def test_final_results_button_uses_parent_native_results_and_distinct_route(monkeypatch):
    root = campaign_rows()[0]
    monkeypatch.setattr(config, "identity", lambda: "organizer")
    assert ui._job_details_reference(root) == "campaign/campaign-run"
    calls = []
    monkeypatch.setattr(ui.st, "session_state", {"clicked": {"row": 0}})
    monkeypatch.setattr(ui, "_navigate_to_job_details", lambda ref, **kw: calls.append((ref, kw)))
    ui._request_results_pipeline_results("clicked", [ui._job_details_reference(root)])
    assert calls == [("campaign/campaign-run", {"final_results": True})]


def test_parent_metadata_survives_job_catalog_indexing(isolated):
    from mn_molecular_design_studio.pipeline_jobs import link_pipeline_child
    from mn_molecular_design_studio import job_catalog
    run = isolated / "runs/ligand-docking/batch"
    run.mkdir(parents=True)
    (run / "metadata.json").write_text(json.dumps({"tool": "unidock_pro", "status": "completed"}))
    link_pipeline_child(run, parent_run_dir="/different/computer/campaign/root", parent_node_id="organizer",
        pipeline_name="A looping campaign", step="Iteration 3 · docking", step_index=4, round_number=3)
    job_catalog._upsert(isolated / "runs", run, refresh_artifacts=False)
    with job_catalog._connection() as c:
        row = c.execute("SELECT * FROM jobs WHERE run_id='batch'").fetchone()
        assert (row["campaign_id"], row["pipeline_parent_node_id"], row["pipeline_iteration"]) == ("root", "organizer", 3)
        assert row["pipeline_name"] == "A looping campaign"


def test_job_inspection_returns_bounded_combined_results_without_import(isolated, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.services import local_automation
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    run = isolated / "runs/campaign/root"
    data = run / "artifacts/active_learning"
    data.mkdir(parents=True)
    rows = [dict(compound_id=f"ID-{i}", smiles="CCO", score=float(i)) for i in range(2000)]
    pq.write_table(pa.Table.from_pylist(rows), data / "docking_scores.parquet")
    monkeypatch.setattr(local_automation, "inspect_job", lambda *a, **kw: {
        "job": {"run_dir": str(run)}, "result": {}, "candidate_count": 0, "campaign_result_path": None})
    result = local_automation.job_results("campaign/root")
    assert result["screening_results"]["compound_count"] == 2000
    assert len(result["screening_results"]["rows"]) == 100
    with MolecularDatabase().connection() as c:
        assert c.execute("SELECT count(*) FROM molecules").fetchone()[0] == 0


def test_combined_results_page_displays_history_and_scores(isolated):
    app = AppTest.from_string('''
from mn_molecular_design_studio.screening_results_ui import render_summary
render_summary({"compound_count": 2000,"rows":[{"compound_id":"X","score":-8,"smiles":"CCO"}],
 "scores_path":"/external/scores.parquet","rounds":[{"round":1,"successful_count":20},{"round":2,"successful_count":20}]},
 {"round_count":2,"best_score":-8,"stop_reason":"stable_top_scores"},name="Campaign",status="completed")
''').run()
    assert not app.exception
    assert app.title[0].value == "Final screening results"
    assert len(app.dataframe[0].value) == 2
    assert app.dataframe[1].value.iloc[0]["compound_id"] == "X"
    # The application must also parse the native route that invokes this view.
    source = Path(__file__).parents[1] / "mn_molecular_design_studio/run_app.py"
    ast.parse(source.read_text())


def test_saturation_uses_complete_top_vectors_and_keeps_zero_change():
    from mn_molecular_design_studio.screening_results_ui import saturation_data
    rows = saturation_data([
        {"round": 3, "best_score": -10, "top_n": 2, "top_scores": [{"score": -10}, {"score": -8}],
         "maximum_rank_score_change": 0, "stable_round_count": 2},
        {"round": 1, "best_score": "-9", "top_n": 2, "top_scores": [{"score": -9}, {"score": -5}]},
        {"round": 2, "best_score": float("nan"), "top_scores": [{"score": -9}, {"score": float("inf")}]},
        {"round": "invalid"}, {"round": 0}, {"round": 1.5},
    ])
    assert [row["iteration"] for row in rows] == [1, 2, 3]
    assert rows[0]["top_mean"] == -7
    assert rows[1]["best_score"] is None and rows[1]["top_mean"] is None
    assert rows[2]["top_mean"] == -9 and rows[2]["maximum_change"] == 0
    assert rows[2]["stable_comparisons"] == 2
    assert saturation_data([]) == []


def test_convergence_charts_use_engine_units_and_actual_stopping_threshold(isolated):
    app = AppTest.from_string('''
from mn_molecular_design_studio.screening_results_ui import render_summary
render_summary({"compound_count": 2,"rows":[],"scores_path":"/external/scores.parquet", "score_units":"REU",
 "score_direction":"Lower is better", "stopping_criteria":{"score_tolerance":0.05,"stable_rounds":3},
 "rounds":[{"round":1,"best_score":-10,"top_n":2,"top_scores":[{"score":-10},{"score":-8}]},
    {"round":2,"best_score":-10,"top_n":2,"top_scores":[{"score":-10},{"score":-8}],
 "maximum_rank_score_change":0,"stable_round_count":1}]},
 {"round_count":2,"best_score":-10},name="Campaign",status="running")
''').run()
    assert not app.exception
    charts = app.get("vega_lite_chart")
    assert len(charts) == 2
    assert "REU" in charts[0].proto.spec and "kcal/mol" not in charts[0].proto.spec
    assert '"Tolerance"' in charts[1].proto.spec
    assert any("Stable comparisons: 1/3 required" in caption.value for caption in app.caption)


def test_top_100_1000_10000_are_shown_together(isolated):
    app = AppTest.from_string('''
from mn_molecular_design_studio.screening_results_ui import render_saturation
render_saturation({"rounds":[{"round":1,"best_score":-17,"score_summaries":[
 {"top_n":100,"mean_score":-16,"count":100},
 {"top_n":1000,"mean_score":-15,"count":1000},
 {"top_n":10000,"mean_score":-14,"count":10000}]}]}, {})
''').run()
    assert not app.exception
    spec = json.loads(app.get("vega_lite_chart")[0].proto.spec)
    assert spec["encoding"]["color"]["scale"]["domain"] == ["Best score", "Top 100", "Top 1,000", "Top 10,000"]


def test_inspection_publishes_declared_units_and_stopping_criteria(isolated):
    from mn_molecular_design_studio.bulk_screening_results import screening_results_payload
    run = isolated / "runs/campaign/root"
    data = run / "artifacts/active_learning"
    data.mkdir(parents=True)
    pq.write_table(pa.Table.from_pylist([{"compound_id": "A", "score": -8}]), data / "docking_scores.parquet")
    (data / "round_history.json").write_text(json.dumps({"score_units": "REU", "rounds": []}))
    (run / "graph_snapshot.json").write_text(json.dumps({"nodes": [{"data": {"kind": "active_decision",
        "params": {"top_n": 100, "score_tolerance": 0.05, "stable_rounds": 3, "max_iterations": 10}}}]}))
    payload = screening_results_payload(run)
    assert payload["score_units"] == "REU"
    assert payload["stopping_criteria"] == {"top_n": 100, "score_tolerance": 0.05, "stable_rounds": 3, "max_iterations": 10}
