import json

from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio.data_browser import foldseek_membership
from mn_molecular_design_studio.foldseek_results_ui import inspect_clustering_pipeline, saved_members


def test_foldseek_results_open_on_shortlist_with_exact_saved_structure(tmp_path):
    candidates = [
        {"candidate_id": "design-1", "metrics": {"diversity_cluster": 1, "diversity_cluster_size": 2},
         "raw_metadata": {"name": "Design one", "studio_molecule_id": "mol-1",
                          "metadata": {"source_group_id": "saved-group", "structure_id": "refolded-1"}}},
        {"candidate_id": "design-2", "metrics": {"diversity_cluster": 1, "diversity_cluster_size": 2},
         "raw_metadata": {"name": "Design two", "studio_molecule_id": "mol-2",
                          "metadata": {"source_group_id": "saved-group", "structure_id": "refolded-2"}}},
    ]
    (tmp_path / "assignments.jsonl").write_text(
        "".join(json.dumps(row) + "\n" for row in candidates))
    (tmp_path / "shortlist.jsonl").write_text(json.dumps(candidates[0]) + "\n")
    (tmp_path / "graph_snapshot.json").write_text(json.dumps({
        "nodes": [{"id": "analysis", "data": {"kind": "foldseek"}},
                  {"id": "select", "data": {"kind": "filter", "params": {"diversity_field": "diversity_cluster"}}}],
        "edges": [{"source": "analysis", "target": "select"}],
    }))
    (tmp_path / "graph_execution.json").write_text(json.dumps({
        "nodes": {"analysis": {"output": {"candidates_jsonl": "assignments.jsonl"}},
                  "select": {"output": {"candidates_jsonl": "shortlist.jsonl"}}},
    }))
    job = {"tool": "visual_pipeline", "run_dir": str(tmp_path)}
    summary = inspect_clustering_pipeline(job)
    assert len(summary["rows"]) == 2
    selected = [row for row in summary["rows"] if row["selected"]]
    expected = {"molecule_id": "mol-1", "group_id": "saved-group", "structure_id": "refolded-1"}
    assert saved_members(selected) == [expected]

    app = AppTest.from_string(
        "from mn_molecular_design_studio.foldseek_results_ui import inspect_clustering_pipeline, render_clustering_results\n"
        f"render_clustering_results(inspect_clustering_pipeline({job!r}), 'project', 'test-run')"
    ).run()
    assert not app.exception
    assert len(app.dataframe[0].value) == 1
    app.get("button_group")[0].set_value("All designs").run()
    assert len(app.dataframe[0].value) == 2
    app.button[0].click().run()
    assert app.session_state["data-browser-open"]["foldseek_members"] == [expected]


def test_foldseek_membership_checks_original_structure_and_collection():
    locator = {"molecule_id": "mol-1", "group_id": "saved-group", "structure_id": "refolded-1"}
    outcome = {"originals": [
        {"id": "source", "molecule_id": "mol-1", "group_id": "saved-group", "structure_ids": ["refolded-1"]},
        {"id": "generated", "molecule_id": "mol-1", "group_id": "saved-group", "structure_ids": ["generated-1"]},
        {"id": "other-group", "molecule_id": "mol-1", "group_id": "other", "structure_ids": ["refolded-1"]},
    ]}
    assert foldseek_membership([outcome], [locator]) == ({"source"}, {"mol-1": "refolded-1"})
    assert foldseek_membership([outcome], [{**locator, "structure_id": "missing"}]) == (set(), {})
