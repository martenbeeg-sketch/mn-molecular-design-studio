from __future__ import annotations

from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest
from typer.testing import CliRunner

from mn_molecular_design_studio.cluster import cli, config, transport, ui


@pytest.fixture
def isolated(monkeypatch, tmp_path):
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(tmp_path / "studio"))
    return tmp_path


def test_settings_render_without_enabling_sharing(isolated):
    at = AppTest.from_string('from mn_molecular_design_studio.cluster.ui import settings\nsettings({"id":"admin", "role":"admin"}, "token")').run()
    assert not at.exception
    assert at.selectbox[0].value == "disabled"
    assert config.load()["mode"] == "disabled"


def test_workstation_resources_are_visible_to_an_ordinary_user(isolated, monkeypatch):
    config.save({"mode": "worker", "coordinator_url": "http://127.0.0.1:18766", "default_placement": "auto"})
    node = {"id": "remote", "name": "4032", "availability": "online", "snapshot": {"resources": {
        "free_gpus": [{"gpu_id": 0}, {"gpu_id": 1}], "cpu_free": 16, "ram_available_gb": 64}}}
    monkeypatch.setattr(ui, "_request", lambda *a, **k: {"nodes": [node]})
    at = AppTest.from_string('from mn_molecular_design_studio.cluster.ui import execution_selector\nexecution_selector({"id":"scientist", "role":"user"}, "token")').run()
    assert not at.exception
    assert any("2 free GPUs" in c.value for c in at.caption)
    assert at.selectbox[0].value == "auto"


def test_shared_results_can_page_beyond_first_hundred_jobs(isolated, monkeypatch):
    rows = [{"cluster_job_id": f"job-{i}", "node_id": "remote", "node_name": "4032", "project_id": "p", "campaign_name": f"Design {i}",
             "status": "completed", "availability": "offline", "created_at": f"{i:05d}"} for i in range(121)]
    monkeypatch.setattr(ui, "catalog", lambda *a: ({"items": rows, "total": len(rows)}, False))
    from mn_molecular_design_studio import workspace
    monkeypatch.setattr(workspace, "load_organization", lambda: {"projects": [{"id": "p", "name": "Project"}]})
    at = AppTest.from_string('from mn_molecular_design_studio.cluster.ui import render_runs\nrender_runs({"id":"scientist"}, "token", lambda **kw: {"items":[], "total":0})').run()
    assert not at.exception
    assert len([b for b in at.button if b.label == "View"]) == 50
    at.number_input[0].set_value(3).run()
    assert not at.exception
    assert len([b for b in at.button if b.label == "View"]) == 21


def test_remote_candidate_results_include_sequence_and_scores():
    rows = ui._candidate_result_rows([{
        "candidate_id": "rf2_model_1",
        "binder_sequence": "QVQLVESGGGLVQPGGSLRLSCAAS",
        "binder_length": 127,
        "metrics": {"pred_lddt": 0.9, "interaction_pae": 9.57, "hotspot_atom_contact_fraction": 1.0},
    }])
    assert rows == [{
        "Rank": 1,
        "Design": "rf2_model_1",
        "Sequence": "QVQLVESGGGLVQPGGSLRLSCAAS",
        "Length": 127,
        "pLDDT": 0.9,
        "Interaction PAE": 9.57,
        "Hotspot atom contacts": 1.0,
    }]


def test_candidate_structure_uses_confidence_colors_only_with_per_atom_scores():
    candidate = {"metrics": {"pred_lddt": 0.9}}
    atom_lines = "\n".join(("ATOM  " + " " * 54 + f"{score:6.2f}") for score in (55.0, 90.0))
    assert ui._candidate_structure_color(candidate) == "chain-id"
    assert ui._candidate_structure_color(candidate, atom_lines) == "plddt"
    assert ui._candidate_structure_color({"metrics": {"binder_length": 127}}, atom_lines) == "chain-id"


def test_remote_candidate_structure_path_is_run_relative_and_confined():
    data = {"job": {"run_dir": "/data/results/runs/job-1"}}
    assert ui._candidate_structure_path({
        "complex_pdb": "/data/results/runs/job-1/artifacts/design/model.pdb",
    }, data) == ("artifacts/design/model.pdb", "complex_pdb")
    assert ui._candidate_structure_path({
        "complex_pdb": "/data/results/runs/other/model.pdb",
        "binder_pdb": "../../outside.pdb",
    }, data) is None


def test_remote_job_displays_sequence_and_loads_a_representative_structure(isolated, monkeypatch):
    import importlib

    candidate = {
        "candidate_id": "rf2_model_1",
        "binder_sequence": "QVQLVESGGGLVQPGGSLRLSCAAS",
        "complex_pdb": "/data/results/runs/job-1/artifacts/design/model.pdb",
        "binder_length": 127,
        "metrics": {"pred_lddt": 0.9, "interaction_pae": 9.57},
    }
    row = {
        "cluster_job_id": "job-1", "campaign_name": "Remote antibody", "tool": "rfantibody",
        "status": "completed", "availability": "online", "node_name": "pc004032",
        "local_reference": "design/run-1", "project_id": "project-1", "artifacts": [],
    }
    structure_requests = []
    viewer_calls = []

    def request(path, token, payload=None, method="GET", **kwargs):
        if path == "/jobs/job-1":
            return row
        if path == "/jobs/job-1/actions" and payload["action"] == "inspect":
            return {"command_id": "inspect-command"}
        if path == "/commands/inspect-command":
            return {"status": "completed", "response": {
                "job": {"run_dir": "/data/results/runs/job-1"},
                "results": {"candidates": [candidate]}, "input": {}, "result": {},
            }}
        if path == "/jobs/job-1/actions" and payload["action"] == "artifact":
            structure_requests.append(payload["parameters"]["path"])
            return {"command_id": "structure-command"}
        if path == "/commands/structure-command":
            return {"status": "completed", "response": {"name": "model.pdb", "size_bytes": 12}}
        if path == "/commands/structure-command/file":
            return b"ATOM sample\n"
        raise AssertionError(f"Unexpected cluster request: {path} {payload}")

    monkeypatch.setattr(ui, "_request", request)
    viewer = importlib.import_module("mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer")
    monkeypatch.setattr(viewer, "molstar_custom_component", lambda *args, **kwargs: viewer_calls.append((args, kwargs)))
    at = AppTest.from_string(
        'from mn_molecular_design_studio.cluster.ui import render_job\n'
        'render_job("cluster:job-1", {"id":"scientist"}, "token", lambda *a, **k: {})'
    ).run()
    assert not at.exception
    next(button for button in at.button if button.label == "Load full results from workstation").click().run()
    assert not at.exception
    assert structure_requests == ["artifacts/design/model.pdb"]
    assert viewer_calls
    visualization = viewer_calls[0][0][0][0]
    assert visualization.color == "chain-id"
    assert visualization.representation_type == "cartoon"
    assert any("QVQLVESGGGLVQPGGSLRLSCAAS" in str(table.value) for table in at.dataframe)
    assert any("How to interpret these scores" in block.value for block in at.markdown)


def test_services_use_the_app_conda_and_private_configuration_home(isolated, monkeypatch):
    monkeypatch.setattr(Path, "home", lambda: isolated)
    runner = CliRunner()
    result = runner.invoke(cli.app, ["init", "--name", "workstation"])
    assert result.exit_code == 0, result.output
    result = runner.invoke(cli.app, ["install-services"])
    assert result.exit_code == 0, result.output
    unit = (isolated / ".config/systemd/user/studio-agent.service").read_text()
    import sys
    assert sys.executable in unit and str(config.home()) in unit
    secret = config.load()["node_token"]
    assert secret not in unit and secret not in runner.invoke(cli.app, ["status"]).output


def test_local_placement_refuses_unavailable_results_disk(isolated, monkeypatch):
    from mn_molecular_design_studio.cluster import service
    from mn_molecular_design_studio.cluster.agent import Agent
    config.save({"mode": "coordinator", "coordinator_url": "http://127.0.0.1:8766"})
    monkeypatch.setattr(Agent, "storage_online", lambda self: False)
    with pytest.raises(transport.ClusterUnavailable, match="Results storage"):
        service.submit_from_api({"execution": {"placement": "local"}}, {"id": "user"}, request_dir=isolated, token="", authorize_path=lambda p: None)
