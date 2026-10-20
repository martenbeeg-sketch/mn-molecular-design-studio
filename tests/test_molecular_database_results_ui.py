from __future__ import annotations

import hashlib
from pathlib import Path

from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio.molecular_database import MolecularDatabase


USER = {"id": "a" * 32, "role": "user", "_cluster_project_ids": []}


def _isolate_workspace(tmp_path, monkeypatch):
    workspace = tmp_path / "workspace"
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(workspace))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(workspace))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(workspace / "runs"))


def _result_rows(group_id: str):
    group = {
        "id": group_id,
        "name": "Source compounds",
        "project_id": "",
        "owner_id": USER["id"],
        "import_method": "user_uploads",
        "origin": {"type": "user"},
    }
    return [
        {"id": "molecule-a", "name": "Ethanol", "type": "small_molecule", "components": [{"chain": "A", "smiles": "CCO"}], "groups": [group], "structures": [], "observations": []},
        {"id": "molecule-b", "name": "Acetone", "type": "small_molecule", "components": [{"chain": "A", "smiles": "CC(=O)C"}], "groups": [group], "structures": [], "observations": []},
    ]


def test_results_action_saves_exact_selected_subset(tmp_path, monkeypatch):
    _isolate_workspace(tmp_path, monkeypatch)
    from mn_molecular_design_studio import molecular_database_results, molecular_database_results_ui

    database = MolecularDatabase(tmp_path / "database")
    created = database.import_records(
        USER,
        [{"name": "Starting compound", "components": [{"chain": "A", "type": "small_molecule", "smiles": "CC"}]}],
        name="Existing curated set",
    )
    source_rows = _result_rows(created["group_id"])
    calls = []
    monkeypatch.setattr(molecular_database_results_ui, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(molecular_database_results, "list_run_molecules", lambda *args, **kwargs: source_rows)

    def save_selection(run_dir, user, molecule_ids, **kwargs):
        calls.append((list(molecule_ids), kwargs))
        return {"group_id": kwargs["group_id"], "molecule_ids": list(molecule_ids), "saved_count": len(molecule_ids)}

    monkeypatch.setattr(molecular_database_results, "save_run_selection", save_selection)
    run_dir = (tmp_path / "runs" / "vina" / "run-1").resolve()
    run_key = hashlib.sha1(str(run_dir).encode()).hexdigest()[:12]
    script = f'''\
from pathlib import Path
from mn_molecular_design_studio.molecular_database_results_ui import render_run_database_actions
render_run_database_actions(Path({str(run_dir)!r}), {USER!r}, [])
'''
    app = AppTest.from_string(script, default_timeout=10).run()
    assert not app.exception
    app.multiselect(key=f"run-database-selection-{run_key}").set_value(["molecule-b"]).run()
    app.button(key=f"run-database-save-{run_key}").click().run()

    assert not app.exception
    assert calls[0][0] == ["molecule-b"]
    assert calls[0][1]["group_id"] == created["group_id"]
    assert any("Saved 1 selected molecule" in item.value for item in app.success)


def test_open_in_database_preselects_source_group_and_all_indexed_mode(tmp_path, monkeypatch):
    _isolate_workspace(tmp_path, monkeypatch)
    from mn_molecular_design_studio import molecular_database_results, molecular_database_results_ui

    source_group = {
        "id": "automatic-output-group",
        "name": "Docking inputs",
        "project_id": "",
        "owner_id": USER["id"],
        "import_method": "tool_inputs",
        "origin": {"type": "job"},
    }
    record = {"id": "molecule-a", "name": "Ethanol", "type": "small_molecule", "components": [{"chain": "A", "smiles": "CCO"}], "groups": [source_group], "source_group_id": source_group["id"], "structures": [], "observations": []}
    database = MolecularDatabase(tmp_path / "database")
    database.create_group(
        USER,
        source_group["name"],
        import_method="tool_inputs",
        origin={"type": "job"},
        group_id=source_group["id"],
    )
    monkeypatch.setattr(molecular_database_results_ui, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(molecular_database_results, "list_run_molecules", lambda *args, **kwargs: [record])
    run_dir = (tmp_path / "runs" / "vina" / "run-2").resolve()
    run_key = hashlib.sha1(str(run_dir).encode()).hexdigest()[:12]
    script = f'''\
from pathlib import Path
from mn_molecular_design_studio.molecular_database_results_ui import render_run_database_actions
render_run_database_actions(Path({str(run_dir)!r}), {USER!r}, [])
'''
    app = AppTest.from_string(script, default_timeout=10).run()
    assert not app.exception
    app.button(key=f"run-database-open-{run_key}").click().run()

    assert not app.exception
    assert app.session_state["page"] == "Database"
    assert app.session_state["molecular-preselect-group"] == source_group["id"]
    assert app.session_state["molecular-preselect-collection"] == "all"
    assert app.session_state["molecular-preselect-molecule"] == "molecule-a"


def test_database_defaults_to_library_and_passes_collection_to_queries(tmp_path, monkeypatch):
    _isolate_workspace(tmp_path, monkeypatch)
    from mn_molecular_design_studio import database_ui

    calls = []
    database = MolecularDatabase(tmp_path / "database")

    class DatabaseSpy:
        def groups(self, user, **kwargs):
            calls.append(("groups", kwargs.get("collection")))
            return database.groups(user, **kwargs)

        def metric_fields(self, user, **kwargs):
            calls.append(("metric_fields", kwargs.get("collection")))
            return database.metric_fields(user, **kwargs)

        def list_molecules(self, user, **kwargs):
            calls.append(("list_molecules", kwargs.get("collection")))
            return database.list_molecules(user, **kwargs)

    monkeypatch.setattr(database_ui, "MolecularDatabase", DatabaseSpy)
    monkeypatch.setattr(database_ui, "sync_workspace", lambda *args, **kwargs: [])
    script = f'''\
from mn_molecular_design_studio.database_ui import render_database
render_database({USER!r}, [])
'''
    app = AppTest.from_string(script, default_timeout=10).run()
    assert not app.exception
    assert app.session_state["molecular-collection-view"] == "Library"
    assert {kind for _, kind in calls} == {"library"}

    app.session_state["molecular-collection-view"] = "All indexed"
    app.run()
    assert not app.exception
    assert {kind for _, kind in calls if kind is not None} == {"all", "library"}
    assert calls[-1] == ("list_molecules", "all")
