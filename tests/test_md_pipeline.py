import copy
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio import md_pipeline, pipeline_graph
from mn_molecular_design_studio.engine.workflow_engine.tools.md_analysis import CoordinateFluctuations, radius_of_gyration
from mn_molecular_design_studio.engine.workflow_engine.workflows import graph_pipeline as runner


def _json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload))
    return path


def test_md_pipeline_declares_real_trajectory_ports_and_executable_tool():
    graph = pipeline_graph.new_graph("md_analysis")
    plan = pipeline_graph.validate_graph(graph)
    assert plan["order"] == ["md_source", "md_analysis", "results"]
    assert pipeline_graph.NODE_CATALOG["md_analysis"]["workflow"] == "ligand.md.analysis"
    assert pipeline_graph.connection_compatibility("md_source", "md_analysis")["compatible"]
    assert pipeline_graph.connection_compatibility("md_analysis", "results")["compatible"]
    assert not pipeline_graph.connection_compatibility("boltz2", "md_analysis")["compatible"]
    assert not pipeline_graph.connection_compatibility("md_source", "ipsae")["compatible"]
    with pytest.raises(ValueError, match="completed MD job"):
        pipeline_graph.validate_graph(graph, require_inputs=True)


@pytest.mark.parametrize("field,value", [("stride", 1.5), ("stride", True), ("start_ns", float("nan")),
    ("end_ns", "0"), ("end_ns", "invalid"), ("cpu_cores", 0), ("dcd_frame_interval_ps", "10"),
    ("ligand_selection", "")])
def test_md_pipeline_rejects_invalid_analysis_configuration(field, value):
    graph = pipeline_graph.new_graph("md_analysis")
    graph["nodes"][1]["data"]["params"][field] = value
    with pytest.raises(ValueError):
        pipeline_graph.validate_graph(graph)


def test_md_node_request_preserves_exact_lineage_without_joining_segments():
    graph = pipeline_graph.new_graph("md_analysis")
    params = graph["nodes"][1]["data"]["params"]
    params.update(end_ns="25", dcd_frame_interval_ps="10", dcd_first_time_ps="1000")
    source = {"md_source": {"source_run_dir": "/runs/ligand-md/source", "topology_path": "/runs/ligand-md/source/system.pdb",
        "trajectory_path": "/runs/ligand-md/source/trajectory.dcd",
        "source_segments": [{"study_id": "study-a", "run_id": "segment-2", "previous_duration_ns": 1., "segment_duration_ns": 99.}]}}
    request = runner.node_request("md_analysis", params, source, "test")
    assert request["workflow"] == "ligand.md.analysis"
    assert request["parameters"]["end_ns"] == 25.
    assert request["parameters"]["source_segments"] == source["md_source"]["source_segments"]
    assert request["parameters"]["dcd_first_time_ps"] == 1000.
    with pytest.raises(ValueError, match="Saved MD run"):
        runner.node_request("md_analysis", params, {}, "test")


@pytest.fixture
def md_source(tmp_path, monkeypatch):
    from mn_molecular_design_studio import auth, job_catalog, md_analysis_ui
    from mn_molecular_design_studio.engine.workflow_engine import runtime
    from mn_molecular_design_studio.engine.workflow_engine.workflows import md_analysis
    root = tmp_path / "runs"
    source = root / "ligand-md" / "source"
    _json(source / "metadata.json", {"status": "completed"})
    topology = source / "system.pdb"
    trajectory = source / "trajectory.xtc"
    topology.write_bytes(b"topology")
    trajectory.write_bytes(b"trajectory")
    monkeypatch.setattr(runtime, "runs_root", lambda: root)
    monkeypatch.setattr(md_analysis, "runs_root", lambda: root)
    access = {"owner_id": "owner", "project_id": "project"}
    monkeypatch.setattr(job_catalog, "get_job_access_record", lambda *args: access)
    monkeypatch.setattr(auth, "list_users", lambda: [{"id": "owner", "role": "user"}, {"id": "other", "role": "user"}])
    monkeypatch.setattr(auth, "require_job_access", lambda user, record: None if user["id"] == "owner" else (_ for _ in ()).throw(PermissionError("not accessible")))
    monkeypatch.setattr(md_analysis_ui, "_source_options_for_job", lambda row: [{"option_id": "segment-2", "label": "Ligand A · replica 2",
        "trajectory_path": trajectory, "topology_paths": [topology], "study": {"study_id": "study-a"},
        "segment": {"run_id": "production-2", "replica_index": 2, "segment_duration_ns": 99.}}])
    return {"source_run_dir": str(source), "option_id": "segment-2", "topology_relative_path": "system.pdb"}, access


def test_pipeline_source_requires_access_exact_topology_and_project(md_source, tmp_path):
    params, access = md_source
    original = (Path(params["source_run_dir"]) / "trajectory.xtc").read_bytes()
    output = md_pipeline.source_output(params, owner_id="owner", project_id="project")
    assert output["md_source"]["source_segments"][0]["run_id"] == "production-2"
    assert output["md_source"]["expected_source_sha256"]["trajectory"] == hashlib.sha256(original).hexdigest()
    assert runner._cache_valid(output)
    assert (Path(params["source_run_dir"]) / "trajectory.xtc").read_bytes() == original
    with pytest.raises(PermissionError):
        md_pipeline.source_output(params, owner_id="other", project_id="project")
    with pytest.raises(ValueError, match="pipeline project"):
        md_pipeline.source_output(params, owner_id="owner", project_id="other")
    with pytest.raises(ValueError, match="no longer available"):
        md_pipeline.source_output({**params, "option_id": "other-replica"}, owner_id="owner", project_id="project")
    other = tmp_path / "other.pdb"
    other.write_bytes(b"different topology")
    with pytest.raises(ValueError, match="does not belong"):
        md_pipeline.source_output({**params, "topology_relative_path": str(other)}, owner_id="owner", project_id="project")
    (Path(params["source_run_dir"]) / "trajectory.xtc").write_bytes(b"changed")
    assert not runner._cache_valid(output)


def _analysis_child(tmp_path, *, success=True):
    child = tmp_path / "md-analysis" / "child"
    _json(child / "metadata.json", {"status": "completed" if success else "failed", "task_group": "md-analysis"})
    _json(child / "result.json", {"success": success, "outputs": {}})
    summary = {"success": success, "schema_version": 1, "source_study_id": "study-a",
        "source": {"trajectory": {"sha256": "a" * 64}}, "source_segments": [{"run_id": "production-2"}],
        "metrics": {"selected_frame_count": 5}, "artifacts": {"timeseries_csv": "artifacts/md_analysis/timeseries.csv"}}
    _json(child / "artifacts/md_analysis/summary.json", summary)
    (child / "artifacts/md_analysis/timeseries.csv").write_text("time_ns,ligand_rmsd_angstrom\n0,0\n")
    return child


def test_md_analysis_cache_requires_all_numerical_outputs_and_exact_hashes(tmp_path):
    child = _analysis_child(tmp_path)
    output = md_pipeline.analysis_output(child)
    assert output["analysis_run"]["source_study_id"] == "study-a"
    assert runner._cache_valid(output)
    (child / "artifacts/md_analysis/timeseries.csv").write_text("changed")
    assert not runner._cache_valid(output)


@pytest.mark.parametrize("success", [True, False])
def test_md_pipeline_executes_an_analysis_child_and_collects_it_not_designs(tmp_path, monkeypatch, md_source, success):
    from mn_molecular_design_studio import cache_provenance, workspace
    params, access = md_source
    graph = pipeline_graph.new_graph("md_analysis")
    graph["nodes"][0]["data"]["params"] = params
    monkeypatch.setattr(pipeline_graph, "app_home", lambda: tmp_path / "home")
    monkeypatch.setattr(workspace, "mark_studio_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(cache_provenance, "execution_identity", lambda *args: {"runtime": "test"})
    parent = tmp_path / "campaign" / "parent"
    _json(parent / "metadata.json", {"status": "queued", "job_type": "visual_pipeline"})
    _json(parent / "input.json", {})
    monkeypatch.setattr(runner, "create_job", lambda *args, **kwargs: SimpleNamespace(run_dir=parent))
    child = _analysis_child(tmp_path, success=success)
    requests = []
    monkeypatch.setattr(runner, "_submit", lambda request, *args, **kwargs: (requests.append(copy.deepcopy(request)) or child))
    run = runner.run_graph_pipeline(graph, "MD test", "owner", "project", reuse_completed=False)
    result = json.loads((run / "result.json").read_text())
    state = json.loads((run / "graph_execution.json").read_text())
    assert len(requests) == 1
    assert requests[0]["workflow"] == "ligand.md.analysis"
    assert requests[0]["parameters"]["source_segments"][0]["run_id"] == "production-2"
    assert result["success"] is success
    assert state["nodes"]["results"]["status"] == ("completed" if success else "blocked")
    assert result["metrics"]["candidate_count"] == 0
    if success:
        assert result["outputs"]["surface_analyses"][0]["kind"] == "md_analysis"
        assert state["nodes"]["md_analysis"]["output"]["content_manifest"]


def test_streaming_fluctuations_are_population_rmsf_about_window_mean():
    accumulator = CoordinateFluctuations(2)
    frames = np.asarray([[[0, 0, 0], [1, 0, 0]], [[2, 0, 0], [1, 0, 0]], [[4, 0, 0], [1, 0, 0]]], dtype=float)
    for frame in frames:
        accumulator.update(frame)
    expected = np.sqrt(np.mean(np.sum((frames - frames.mean(axis=0)) ** 2, axis=2), axis=0))
    assert accumulator.rmsf() == pytest.approx(expected)
    with pytest.raises(ValueError, match="finite"):
        accumulator.update(np.full((2, 3), float("nan")))
    assert radius_of_gyration([[0, 0, 0], [2, 0, 0]], [1, 1]) == pytest.approx(1.)
    assert radius_of_gyration([[0, 0, 0], [2, 0, 0]], [1, 3]) == pytest.approx(np.sqrt(.75))


def test_md_pipeline_submission_requests_cpu_and_keeps_child_lineage(tmp_path, monkeypatch):
    from mn_molecular_design_studio import workspace
    from mn_molecular_design_studio.cluster import config
    from mn_molecular_design_studio.engine.workflow_engine.services import workflow_automation
    parent = tmp_path / "parent"
    _json(parent / "metadata.json", {})
    _json(parent / "input.json", {})
    child = _analysis_child(tmp_path)
    captured = {}

    def submit(path):
        captured.update(json.loads(Path(path).read_text()))
        return {"job": {"run_dir": str(child)}}

    monkeypatch.setattr(workflow_automation, "submit_workflow_request", submit)
    monkeypatch.setattr(workspace, "mark_studio_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(config, "identity", lambda: "local-node")
    runner._submit({"workflow": "ligand.md.analysis", "parameters": {"cpu_cores": 3}}, parent,
        {"id": "md_analysis"}, owner_id="owner", project_id="project", gpu_device="all")
    assert captured["resources"] == {"gpu_device": "none", "cpu_cores": 3, "min_vram_gb": 0}
    assert captured["database_collection_policy"] == "none"
    assert json.loads((child / "metadata.json").read_text())["pipeline_parent_node_id"] == "local-node"


def test_md_source_picker_uses_readable_recorded_segments(md_source, monkeypatch):
    from mn_molecular_design_studio import md_analysis_ui
    params, access = md_source
    source = Path(params["source_run_dir"])
    monkeypatch.setattr(md_analysis_ui, "_cached_source_options_for_job", lambda row: [{"option_id": "segment-2",
        "label": "Ligand A · replica 2 · 99 ns", "topology_paths": [source / "system.pdb"]}])
    app = AppTest.from_string(f'''
import streamlit as st
from mn_molecular_design_studio.pipeline_graph import new_graph
from mn_molecular_design_studio.pipeline_graph_ui import _render_node_input
graph = new_graph("md_analysis")
values = _render_node_input(graph["nodes"][0], "owner", "project", lambda **kwargs: {{"items": [{{"run_dir": {str(source)!r}, "campaign_name": "Ligand campaign"}}]}}, "test", graph)
st.json(values)
''').run()
    assert not app.exception
    assert app.selectbox[0].value == str(source)
    assert app.selectbox[1].options == ["Ligand A · replica 2 · 99 ns"]
    assert app.selectbox[2].value == "system.pdb"


def test_rmsf_and_rg_results_are_native_plots_not_just_report_download(tmp_path):
    output = tmp_path / "run/artifacts/md_analysis"
    output.mkdir(parents=True)
    (output / "timeseries.csv").write_text("time_ns,protein_aligned_backbone_rmsd_angstrom,ligand_rmsd_angstrom,protein_rg_angstrom\n0,0,0,10\n1,1,2,11\n")
    (output / "rmsf.csv").write_text("selection,atom_index,atom,residue,rmsf_angstrom\nprotein_ca,1,CA,A:1:ALA,0.5\nligand_heavy,5,C1,A:501:LIG,1.5\n")
    (output / "residue_contacts.csv").write_text("residue,contact_occupancy_fraction\nA:1:ALA,0.5\n")
    summary = {"success": True, "schema_version": 1, "artifacts": {"rmsf_csv": "artifacts/md_analysis/rmsf.csv"}}
    app = AppTest.from_string(f'''
import streamlit as st
from pathlib import Path
from mn_molecular_design_studio.md_results_ui import _render_md_analysis_summary
_render_md_analysis_summary(st, {summary!r}, run_root=Path({str(output.parents[1])!r}))
''').run()
    assert not app.exception
    assert len(app.get("arrow_vega_lite_chart")) + len(app.get("vega_lite_chart")) == 4
