from __future__ import annotations

from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio.cluster import ui


def test_shared_results_view_uses_same_tab_job_route(monkeypatch, tmp_path):
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(tmp_path / "studio"))
    rows = [{
        "cluster_job_id": "job-123",
        "node_id": "remote-node",
        "node_name": "pc004032",
        "project_id": "project-1",
        "campaign_name": "Docking job",
        "tool": "vina",
        "status": "running",
        "availability": "online",
        "created_at": "2026-10-20T10:00:00Z",
    }]
    monkeypatch.setattr(
        ui,
        "catalog",
        lambda *args, **kwargs: ({"items": rows, "total": len(rows)}, False),
    )

    from mn_molecular_design_studio import workspace

    monkeypatch.setattr(workspace, "load_organization", lambda: {"projects": [{"id": "project-1", "name": "Project"}]})
    app = AppTest.from_string(
        'from mn_molecular_design_studio.cluster.ui import render_runs\n'
        'render_runs({"id":"scientist"}, "token", lambda **kw: {"items":[], "total":0})'
    ).run()

    assert not app.exception
    table = app.dataframe[0].value
    assert table.loc["job-123", "View"] == "View"
    calls = []
    monkeypatch.setattr(ui.st, "session_state", {"clicked": {"row": 0}})
    monkeypatch.setattr(ui, "_navigate_to_job_details", lambda ref: calls.append(ref))
    ui._request_results_job_view("clicked", [ui._job_details_reference(rows[0])])
    assert calls == ["cluster:job-123"]
