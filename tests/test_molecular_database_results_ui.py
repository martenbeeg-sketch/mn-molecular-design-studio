from __future__ import annotations

import hashlib
from pathlib import Path

from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio.molecular_database import MolecularDatabase


USER = {"id": "a" * 32, "role": "user", "_cluster_project_ids": []}


def test_database_table_exposes_compound_id_name_and_smiles():
    from mn_molecular_design_studio.database_ui import _table_row

    row = _table_row({
        "id": "mol-1", "name": "Bepristat 1a", "type": "small_molecule",
        "components": [{"chain": "L", "smiles": "CCN"}],
        "metadata": {"compound_id": "70701262"}, "scores": {}, "structures": [],
        "created_at": "2026-10-21T12:00:00Z", "origin": {},
    })
    assert row["values"]["Compound ID"] == "70701262"
    assert row["values"]["Name"] == "Bepristat 1a"
    assert row["values"]["SMILES"] == "CCN"
    assert row["values"]["L"] == "CCN"


def test_passing_collection_target_label_uses_selected_results_provenance():
    from mn_molecular_design_studio.molecular_database_results_ui import _collection_target_labels

    passing = [{"candidate_id": "bindcraft2_00001", "target_identity": ""}]
    records = [{"id": "result-1", "metadata": {
        "target_identity": "Target sequence (chain A · 115 aa)",
        "target_sequences_by_chain": {"A": "A" * 115},
    }}]
    assert _collection_target_labels(passing, records, ["result-1"]) == ["Target sequence (chain A · 115 aa)"]


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


def test_open_in_database_prefers_generated_group_over_tool_inputs(tmp_path, monkeypatch):
    _isolate_workspace(tmp_path, monkeypatch)
    from mn_molecular_design_studio import molecular_database_results, molecular_database_results_ui

    confirmed_user = {"id": "pc002975", "role": "admin", "_cluster_project_ids": []}

    input_group = {
        "id": "tool-input-group",
        "name": "Tool inputs",
        "project_id": "",
        "owner_id": USER["id"],
        "import_method": "tool_inputs",
        "origin": {"type": "job", "collection_kind": "inputs"},
    }
    output_group = {
        "id": "automatic-output-group",
        "name": "RFD4 smoke test generation",
        "project_id": "",
        "owner_id": confirmed_user["id"],
        "import_method": "tool_output",
        "origin": {"type": "job", "collection_kind": "generated", "is_test": True,
                   "run_id": "20261016-234215-9257efd9"},
    }
    input_record = {
        "id": "molecule-input", "name": "1AAY", "type": "protein",
        "components": [{"chain": "A", "sequence": "M"}], "groups": [input_group],
        "source_group_id": input_group["id"], "run_group_ids": [input_group["id"]],
        "run_groups": [{**input_group, "id": input_group["id"]}], "structures": [], "observations": [],
    }
    output_record = {
        "id": "molecule-generated", "name": "rfd4_000001_dna_binder", "type": "complex",
        "components": [{"chain": "A", "type": "dna"}, {"chain": "B", "type": "dna"}, {"chain": "C", "type": "protein"}],
        "groups": [output_group], "source_group_id": output_group["id"],
        "run_group_ids": [output_group["id"]], "run_groups": [{**output_group, "id": output_group["id"]}],
        "structures": [], "observations": [],
    }
    database = MolecularDatabase(tmp_path / "database")
    database.create_group(
        USER,
        input_group["name"],
        import_method="tool_inputs",
        origin=input_group["origin"],
        group_id=input_group["id"],
    )
    database.create_group(
        confirmed_user,
        output_group["name"],
        import_method="tool_output",
        origin=output_group["origin"],
        group_id=output_group["id"],
    )
    monkeypatch.setattr(molecular_database_results_ui, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(molecular_database_results, "list_run_molecules", lambda *args, **kwargs: [input_record, output_record])
    run_dir = (tmp_path / "runs" / "vina" / "run-2").resolve()
    run_key = hashlib.sha1(str(run_dir).encode()).hexdigest()[:12]
    script = f'''\
from pathlib import Path
from mn_molecular_design_studio.molecular_database_results_ui import render_run_database_actions
render_run_database_actions(Path({str(run_dir)!r}), {confirmed_user!r}, [])
'''
    app = AppTest.from_string(script, default_timeout=10).run()
    assert not app.exception
    app.session_state["molecular-hide-smoke-test-groups"] = True
    app.button(key=f"run-database-open-{run_key}").click().run()

    assert not app.exception
    assert app.session_state["page"] == "Database"
    assert app.session_state["molecular-preselect-group"] == output_group["id"]
    assert app.session_state["molecular-preselect-collection"] == "library"
    assert app.session_state["molecular-preselect-molecule"] == "molecule-generated"
    assert app.session_state["molecular-hide-smoke-test-groups"] is False


def test_results_only_candidates_are_not_opened_as_database_records(tmp_path, monkeypatch):
    _isolate_workspace(tmp_path, monkeypatch)
    from mn_molecular_design_studio import molecular_database_results, molecular_database_results_ui

    database = MolecularDatabase(tmp_path / "database")
    row = {"id": "result-1", "molecule_id": "mol-1", "name": "Unselected candidate",
           "type": "protein", "components": [{"type": "protein", "chain": "A", "sequence": "ACDE"}],
           "groups": [], "run_groups": [], "run_group_ids": [], "structures": [], "representations": [],
           "result_only": True, "acceptance_status": "unknown"}
    monkeypatch.setattr(molecular_database_results_ui, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(molecular_database_results, "list_run_molecules", lambda *args, **kwargs: [row])
    run_dir = (tmp_path / "runs" / "design" / "run-1").resolve()
    script = f'''\
from pathlib import Path
from mn_molecular_design_studio.molecular_database_results_ui import render_run_database_actions
render_run_database_actions(Path({str(run_dir)!r}), {USER!r}, [])
'''
    app = AppTest.from_string(script, default_timeout=10).run()
    assert not app.exception
    assert not any(button.label == "Open in Database" for button in app.button)
    assert any("remain in Results" in item.value for item in app.caption)


def test_database_defaults_to_library_and_passes_collection_to_queries(tmp_path, monkeypatch):
    _isolate_workspace(tmp_path, monkeypatch)
    from mn_molecular_design_studio import database_ui

    calls = []
    database = MolecularDatabase(tmp_path / "database")

    class DatabaseSpy:
        path = database.path

        def evaluation_options(self, user, **kwargs):
            return database.evaluation_options(user, **kwargs)

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


def test_saved_collection_navigation_survives_an_empty_automatic_output_list(tmp_path, monkeypatch):
    _isolate_workspace(tmp_path, monkeypatch)
    from mn_molecular_design_studio import molecular_database_results
    monkeypatch.setattr(molecular_database_results, 'list_run_molecules', lambda *args, **kwargs: [])
    monkeypatch.setattr(molecular_database_results, 'list_saved_run_collections', lambda *args, **kwargs: [
        {'id':'accepted-run', 'name':'Accepted binders', 'molecule_ids':['binder-1']}])
    app = AppTest.from_string(f'''
from pathlib import Path
from mn_molecular_design_studio.molecular_database_results_ui import render_run_database_actions
render_run_database_actions(Path({str(tmp_path / 'runs/design/original')!r}), {USER!r}, [])
''').run()
    assert not app.exception
    app.button(key='run-saved-collection-accepted-run').click().run()
    assert not app.exception
    assert app.session_state['complex-analysis-selected-collections'] == ['accepted-run']
    assert app.session_state['complex-analysis-collections-open'] is True
    assert app.session_state['page'] == 'Complex Analysis'


def test_results_collection_requires_explicit_native_stage_and_opens_collection(tmp_path, monkeypatch):
    _isolate_workspace(tmp_path, monkeypatch)
    from mn_molecular_design_studio import molecular_database_results, molecular_database_results_ui

    run_dir = (tmp_path / "runs" / "design" / "native-only").resolve()
    run_key = hashlib.sha1(str(run_dir).encode()).hexdigest()[:12]
    inventory = {
        "candidate_count": 2, "model_count": 2, "run_id": "native-only", "project_id": "",
        "source_tool": "BindCraft2",
        "stages": [{"stage_id": "native", "label": "Native generator acceptance",
                    "passed": [{"candidate_id": "pass-id", "model_id": "model-1",
                                "source_tool": "bindcraft2", "generator_engine": "BindCraft2",
                                "target_identity": "PD-L1",
                                "decision": {"stage_id": "native", "cohort_status": "passed",
                                             "policy": {"native": "bindcraft"}}}],
                    "failed": [{"candidate_id": "fail-id", "model_id": "model-2"}],
                    "pending": [], "unrecorded": [], "passed_design_count": 1,
                    "failed_model_count": 1, "pending_model_count": 0, "unrecorded_model_count": 0}],
    }
    records = [{"id": "result-pass", "result_only": True, "candidate_id": "pass-id", "model_id": "model-1"},
               {"id": "result-fail", "result_only": True, "candidate_id": "fail-id", "model_id": "model-2"}]
    saved = []
    monkeypatch.setattr(molecular_database_results, "run_acceptance_inventory", lambda *_args: inventory)
    monkeypatch.setattr(molecular_database_results, "save_run_selection",
                        lambda _run, _user, ids, **kwargs: saved.append((list(ids), kwargs))
                        or {"group_id": "design-collection-test", "saved_count": 1})
    monkeypatch.setattr(molecular_database_results_ui, "MolecularDatabase", lambda: object())
    app = AppTest.from_string(f'''\
from pathlib import Path
from mn_molecular_design_studio.molecular_database_results_ui import _render_run_passing_collection
_render_run_passing_collection(Path({str(run_dir)!r}), {USER!r}, object(), {records!r}, {run_key!r})
''').run()
    assert not app.exception
    assert app.selectbox(key=f"run-collect-acceptance-stage-{run_key}").value == ""
    assert any("Native generator acceptance: 1 passing design(s)" in item.value for item in app.caption)
    assert any("No final decision is recorded for this run" in item.value for item in app.info)
    assert not saved

    app.selectbox(key=f"run-collect-acceptance-stage-{run_key}").set_value("native").run()
    assert not app.exception
    app.button(key=f"run-collect-passing-{run_key}-native").click().run()
    assert not app.exception
    assert saved[0][0] == ["result-pass"]
    assert saved[0][1]["candidate_models"][0]["candidate_id"] == "pass-id"
    assert saved[0][1]["collection_metadata"]["acceptance_stage_id"] == "native"
    assert saved[0][1]["collection_metadata"]["generator"] == "BindCraft2"

    app.button(key=f"run-collect-open-{run_key}").click().run()
    assert not app.exception
    assert app.session_state["complex-analysis-selected-collections"] == ["design-collection-test"]
    assert app.session_state["page"] == "Complex Analysis"


def test_external_screening_results_wait_for_exact_user_selection(tmp_path, monkeypatch):
    _isolate_workspace(tmp_path, monkeypatch)
    from mn_molecular_design_studio import bulk_screening_results, molecular_database_results_ui

    run_dir = (tmp_path / "runs" / "campaign" / "screening-run").resolve()
    results_path = run_dir / "artifacts" / "active_learning" / "docking_scores.parquet"
    results_path.parent.mkdir(parents=True)
    results_path.write_bytes(b"fixture marker")
    result = {
        "compound_count": 1_590_790,
        "rows": [
            {"compound_id": f"compound-{index}", "compound_name": f"Compound {index}",
             "score": float(-10 + index), "best_round": 3}
            for index in range(1, 4)
        ],
        "scores_path": str(results_path), "files": {},
        "selection": {"group_id": None, "imported_count": 0},
    }
    imported = []
    monkeypatch.setattr(bulk_screening_results, "external_results", lambda *_args: result)
    monkeypatch.setattr(
        bulk_screening_results,
        "import_final_results",
        lambda *args, **kwargs: (
            imported.append((args[2], kwargs.get("selected_compound_ids"))) or
            {"imported_count": len(kwargs.get("selected_compound_ids") or [])}
        ),
    )
    monkeypatch.setattr(molecular_database_results_ui, "MolecularDatabase", lambda: None)
    run_key = hashlib.sha1(str(run_dir).encode()).hexdigest()[:12]
    app = AppTest.from_string(f'''\
from pathlib import Path
from mn_molecular_design_studio.molecular_database_results_ui import render_run_database_actions
render_run_database_actions(Path({str(run_dir)!r}), {USER!r}, [])
''').run()
    assert not app.exception
    assert app.multiselect(key=f"external-screening-{run_key}-selected-compounds").value == []
    assert app.number_input(key=f"external-screening-{run_key}-count").value == 0
    assert app.button(key=f"external-screening-{run_key}-import-selected").disabled
    assert app.button(key=f"external-screening-{run_key}-import").disabled
    assert imported == []

    app.multiselect(key=f"external-screening-{run_key}-selected-compounds").set_value(
        ["compound-1", "compound-3"],
    ).run()
    app.button(key=f"external-screening-{run_key}-import-selected").click().run()
    assert not app.exception
    assert imported == [(0, ["compound-1", "compound-3"])]


def test_historical_spiro_shortlist_result_opens_in_all_indexed(tmp_path, monkeypatch):
    _isolate_workspace(tmp_path, monkeypatch)
    from mn_molecular_design_studio import bulk_screening_results

    run_dir = (tmp_path / "runs" / "campaign" / "historical-screening").resolve()
    results_path = run_dir / "artifacts" / "active_learning" / "docking_scores.parquet"
    results_path.parent.mkdir(parents=True)
    results_path.write_bytes(b"fixture marker")
    group_id = "screening-final-a5ba730c962dc449fb40bbab7702a25b"
    monkeypatch.setattr(bulk_screening_results, "external_results", lambda *_args: {
        "compound_count": 1_590_790, "rows": [], "scores_path": str(results_path), "files": {},
        "selection": {"group_id": group_id, "imported_count": 1000},
    })
    run_key = hashlib.sha1(str(run_dir).encode()).hexdigest()[:12]
    app = AppTest.from_string(f'''\
from pathlib import Path
from mn_molecular_design_studio.molecular_database_results_ui import render_run_database_actions
render_run_database_actions(Path({str(run_dir)!r}), {USER!r}, [])
''').run()
    assert not app.exception
    assert any("original, unreviewed computational shortlist" in item.value for item in app.caption)
    app.button(key=f"external-screening-{run_key}-open").click().run()
    assert not app.exception
    assert app.session_state["page"] == "Database"
    assert app.session_state["molecular-preselect-group"] == group_id
    assert app.session_state["molecular-preselect-collection"] == "all"
