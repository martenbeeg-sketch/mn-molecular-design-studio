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
    view_links = [block.value for block in app.markdown if "job_ref=" in block.value]
    assert len(view_links) == 1
    assert 'href="?job_ref=cluster%3Ajob-123"' in view_links[0]
    assert 'target="_self"' in view_links[0]
