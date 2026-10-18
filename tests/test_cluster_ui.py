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
