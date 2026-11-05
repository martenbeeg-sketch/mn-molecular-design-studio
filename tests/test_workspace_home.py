from __future__ import annotations

import json

from mn_molecular_design_studio import workspace


def test_app_home_is_used_when_workspace_home_is_not_set(tmp_path, monkeypatch):
    app_home = tmp_path / "isolated-studio"
    results = app_home / "saved-results"
    references = app_home / "references"
    app_home.mkdir()
    (app_home / "settings.json").write_text(json.dumps({
        "results_directory": str(results),
        "reference_directory": str(references),
    }))
    monkeypatch.delenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", raising=False)
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(app_home))

    assert workspace._studio_home() == app_home.resolve()
    assert workspace.ensure_studio_workspace() == app_home.resolve()
    assert workspace._studio_home() == app_home.resolve()
    assert results.is_dir()
    assert references.is_dir()


def test_workspace_home_override_has_precedence_over_app_home(tmp_path, monkeypatch):
    workspace_home = tmp_path / "workspace"
    app_home = tmp_path / "runtime"
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(workspace_home))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(app_home))

    assert workspace._studio_home() == workspace_home.resolve()
