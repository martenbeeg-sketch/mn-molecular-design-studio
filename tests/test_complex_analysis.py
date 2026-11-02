from __future__ import annotations

import json
import hashlib
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio import molecule_imports
from mn_molecular_design_studio.complex_analysis_ui import component_inventory
from mn_molecular_design_studio.molecular_database import MolecularDatabase


USER = {"id": "a" * 32, "role": "user", "_cluster_project_ids": []}
OTHER_USER = {"id": "b" * 32, "role": "user", "_cluster_project_ids": []}


def _atom(serial: int, record: str, atom: str, residue: str, chain: str, number: int, x: float) -> str:
    element = atom[0]
    return (
        f"{record:<6}{serial:5d} {atom:^4s} {residue:>3s} {chain}{number:4d}    "
        f"{x:8.3f}{1.0:8.3f}{1.0:8.3f}{1.0:6.2f}{20.0:6.2f}          {element:>2s}"
    )


def _complex_pdb(second_type: str, *, x: float = 2.0) -> bytes:
    lines = [_atom(1, "ATOM", "CA", "ALA", "A", 1, 1.0)]
    if second_type == "protein":
        lines.append(_atom(2, "ATOM", "CA", "GLY", "B", 1, x))
    elif second_type == "small_molecule":
        lines.append(_atom(2, "HETATM", "C1", "LIG", "L", 1, x))
    elif second_type == "dna":
        lines.append(_atom(2, "ATOM", "P", "DA", "D", 1, x))
    elif second_type == "rna":
        lines.append(_atom(2, "ATOM", "P", "A", "R", 1, x))
    else:
        raise ValueError(second_type)
    return ("\n".join([*lines, "END", ""])).encode()


def _isolated_database(tmp_path, monkeypatch) -> MolecularDatabase:
    workspace = tmp_path / "workspace"
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(workspace))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(workspace))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(workspace / "runs"))
    return MolecularDatabase(workspace / "database")


@pytest.mark.parametrize(
    ("second_type", "expected"),
    [
        ("protein", "protein"),
        ("small_molecule", "small_molecule"),
        ("dna", "dna"),
        ("rna", "rna"),
    ],
)
def test_component_inventory_identifies_supported_complex_components(second_type, expected):
    components = component_inventory(_complex_pdb(second_type), "complex.pdb")

    assert {item["type"] for item in components} == {"protein", expected}
    assert {item["source_chain"] for item in components} == {
        "A", {"protein": "B", "small_molecule": "L", "dna": "D", "rna": "R"}[second_type]
    }


def test_structure_display_reuses_molstar_with_all_chains(tmp_path, monkeypatch):
    path = tmp_path / "complex.pdb"
    path.write_bytes(_complex_pdb("rna"))
    captured = []
    from mn_molecular_design_studio.engine.workflow_engine.app.components import molstar_viewer
    from mn_molecular_design_studio.complex_analysis_ui import _show_structure

    monkeypatch.setattr(molstar_viewer, "molstar_custom_component", lambda structures, **kwargs: captured.append((structures, kwargs)))
    script = f'''\
from pathlib import Path
from mn_molecular_design_studio.complex_analysis_ui import _show_structure
_show_structure(Path({str(path)!r}), {{}}, [], key="complex-test")
'''
    app = AppTest.from_string(script, default_timeout=10).run()

    assert not app.exception
    assert len(captured) == 1
    visualization = captured[0][0][0]
    assert " ALA A" in visualization.pdb and " A R" in visualization.pdb
    assert visualization.color == "chain-id"
    assert visualization.representation_type == "cartoon+ball-and-stick"


def test_database_workspace_selects_exact_model_evaluation_and_shows_imported_structure(tmp_path, monkeypatch):
    database = _isolated_database(tmp_path, monkeypatch)
    model_a = _complex_pdb("dna", x=2.0)
    model_b = _complex_pdb("dna", x=9.0)
    imported_record = molecule_imports.parse_structure(model_a, "model-a.pdb")[0]
    imported_record["name"] = "Imported protein–DNA complex"
    imported_record["files"][0].update({"model_index": 1, "run_id": "prediction-1", "candidate_id": "complex-1"})
    imported = database.import_records(USER, [imported_record], name="Imported structures")
    molecule_id = imported["molecule_ids"][0]
    initial = database.list_molecules(USER, molecule_ids=[molecule_id])["items"][0]["representations"][0]
    model_b_id = database.attach_structure(
        imported["group_id"], molecule_id,
        {"name": "model-b.pdb", "content": model_b, "label": "Model 2", "model_index": 2,
         "run_id": "prediction-1", "candidate_id": "complex-1"},
    )
    observation_id = database.add_observation(
        imported["group_id"], molecule_id, "boltz2", "prediction-1", {"iptm": 0.87, "pae": 3.2},
        context={"representation_id": model_b_id, "target_id": "dna-target", "candidate_id": "complex-1",
                 "candidate": {"raw_metadata": {"confidence_outputs": {"ptm": 0.93}}}},
    )
    from mn_molecular_design_studio.evaluations import evaluation_descriptor

    record = database.list_molecules(USER, molecule_ids=[molecule_id])["items"][0]
    observation = next(row for row in record["observations"] if row["id"] == observation_id)
    evaluation_id = evaluation_descriptor(observation)["id"]
    shown = []

    from mn_molecular_design_studio import complex_analysis_ui, molecular_database

    monkeypatch.setattr(molecular_database, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(complex_analysis_ui, "_show_structure", lambda path, context, components, **kwargs: shown.append(path.read_bytes()))
    script = f'''\
import streamlit as st
from mn_molecular_design_studio.complex_analysis_ui import render_complex_analysis_workspace
st.session_state["complex-analysis-source-kind"] = "Database"
st.session_state["complex-analysis-open"] = {{"kind": "database", "molecule_id": {molecule_id!r}, "evaluation_id": {evaluation_id!r}, "structure_id": {model_b_id!r}}}
render_complex_analysis_workspace({USER!r})
'''
    app = AppTest.from_string(script, default_timeout=10).run()

    assert not app.exception
    assert len(shown) == 1
    assert shown[0] == model_b
    assert any("No saved evaluation is associated" in item.value for item in app.info) is False
    saved_json = [item.value for item in app.json]
    assert any("0.87" in json.dumps(value) for value in saved_json)
    assert any("0.93" in json.dumps(value) for value in saved_json)
    assert initial["id"] != model_b_id

    app.multiselect(key=f"complex-analysis-group-members-{model_b_id}-group-1").set_value(["dna:D"]).run()
    assert not app.exception
    assert any("Group 1: dna chain D" in item.value for item in app.caption)


def test_results_workspace_displays_generic_predicted_complex_and_saved_confidence(tmp_path, monkeypatch):
    database = _isolated_database(tmp_path, monkeypatch)
    structure_path = tmp_path / "predicted-complex.pdb"
    structure_path.write_bytes(_complex_pdb("small_molecule"))
    candidate = {
        "candidate_id": "prediction-1",
        "model_index": 1,
        "stage": "complex_refolding",
        "source_tool": "boltz2",
        "metrics": {"iptm": 0.76},
        "raw_metadata": {"confidence_outputs": {"ptm": 0.81, "plddt": 88.4}},
    }
    job = {"task_group": "boltz2", "run_id": "predicted-1", "status": "completed", "tool": "boltz2"}
    option = {
        "path": structure_path,
        "relative_path": "predicted-complex.pdb",
        "candidate": candidate,
        "candidate_id": candidate["candidate_id"],
        "label": structure_path.name,
    }
    shown = []

    from mn_molecular_design_studio import complex_analysis_ui, molecular_database

    monkeypatch.setattr(molecular_database, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(complex_analysis_ui, "_accessible_result_jobs", lambda user: [job])
    monkeypatch.setattr(complex_analysis_ui, "_job_structure_options", lambda user, selected_job: [option])
    monkeypatch.setattr(complex_analysis_ui, "_show_structure", lambda path, context, components, **kwargs: shown.append(path))
    script = f'''\
import streamlit as st
from mn_molecular_design_studio.complex_analysis_ui import render_complex_analysis_workspace
st.session_state["complex-analysis-source-kind"] = "Results"
render_complex_analysis_workspace({USER!r})
'''
    app = AppTest.from_string(script, default_timeout=10).run()

    assert not app.exception
    assert shown == [structure_path]
    assert any("0.76" in json.dumps(item.value) and "0.81" in json.dumps(item.value) for item in app.json)


def test_results_workspace_submits_exact_biotite_model_and_component_groups(tmp_path, monkeypatch):
    database = _isolated_database(tmp_path, monkeypatch)
    from mn_molecular_design_studio import complex_analysis_ui, molecular_database

    source_run = tmp_path / "workspace" / "runs" / "prediction" / "source-job"
    source_run.mkdir(parents=True)
    structure_path = source_run / "model-2.pdb"
    structure_path.write_bytes(_complex_pdb("protein"))
    candidate = {
        "candidate_id": "generic-complex-42", "model_index": 2,
        "stage": "complex_refolding", "source_tool": "generic-predictor",
        "complex_pdb": "model-2.pdb", "binder_chains": ["B"], "target_chains": ["A"],
        "metrics": {"confidence": 0.71}, "raw_metadata": {"model_index": 2},
    }
    candidate_file = source_run / "artifacts" / "normalized_candidates" / "candidates.jsonl"
    candidate_file.parent.mkdir(parents=True)
    candidate_file.write_text(json.dumps(candidate) + "\n")
    (source_run / "metadata.json").write_text(json.dumps({"tool": "generic-predictor", "status": "completed"}))
    (source_run / "input.json").write_text(json.dumps({"tool": "generic-predictor", "params": {}}))
    (source_run / "result.json").write_text(json.dumps({"success": True, "metrics": {}}))
    (source_run / "molecular_design_studio.json").write_text(json.dumps({"studio_owner_id": USER["id"], "project_id": "test-project"}))
    job = {"task_group": "prediction", "run_id": "source-job", "status": "completed",
           "tool": "generic-predictor", "run_dir": str(source_run), "project_id": "test-project"}
    option = {"path": structure_path, "relative_path": "model-2.pdb", "candidate": candidate,
              "candidate_id": candidate["candidate_id"], "label": "Model 2"}

    monkeypatch.setattr(molecular_database, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(complex_analysis_ui, "_accessible_result_jobs", lambda user: [job])
    monkeypatch.setattr(complex_analysis_ui, "_job_structure_options", lambda user, selected_job: [option])
    monkeypatch.setattr(complex_analysis_ui, "_show_structure", lambda *args, **kwargs: None)
    monkeypatch.setattr(complex_analysis_ui, "_show_human_review", lambda *args, **kwargs: None)
    monkeypatch.setattr(complex_analysis_ui, "_matching_biotite_jobs", lambda *args, **kwargs: [])
    submitted = []

    def submit(schema, workflow, parameters, project, settings, title):
        submitted.append((schema, workflow, dict(parameters), project, settings, title))
        import streamlit as st
        st.session_state["last_submitted_job"] = {"run_id": "analysis-1"}

    key = hashlib.sha256(b"prediction/source-job:model-2.pdb").hexdigest()[:14]
    script = f'''\
import streamlit as st
from mn_molecular_design_studio.complex_analysis_ui import render_complex_analysis_workspace
st.session_state["complex-analysis-source-kind"] = "Results"
render_complex_analysis_workspace(
    {USER!r},
    workflow_schema=lambda workflow_id: {{"workflow": workflow_id, "schema_version": 1, "parameters": [{{"name": "contact_distance", "default": 5.0}}]}},
    project_selector=lambda prefix, **kwargs: "test-project",
    submit={submit!r},
)
'''
    # AppTest executes the submitted callback in its own script namespace.
    script = script.replace(f"submit={submit!r}", "submit=submit")
    script = script.replace("render_complex_analysis_workspace(\n", "def submit(schema, workflow, parameters, project, settings, title):\n    st.session_state['submitted_parameters'] = parameters\n    st.session_state['submitted_workflow'] = workflow\n    st.session_state['last_submitted_job'] = {'run_id': 'analysis-1'}\nrender_complex_analysis_workspace(\n")
    app = AppTest.from_string(script, default_timeout=10).run()

    assert not app.exception
    run_button = app.button(key=f"complex-analysis-run-biotite-{key}")
    assert run_button is not None and not run_button.disabled
    run_button.click().run()
    assert not app.exception
    assert app.session_state["submitted_workflow"] == "analysis.biotite_interface"
    parameters = app.session_state["submitted_parameters"]
    assert parameters == {
        "candidate_ids": ["generic-complex-42"],
        "selected_structure_path": str(structure_path.resolve()),
        "selected_model_id": "2",
        "binder_chains": ["B"], "target_chains": ["A"], "contact_distance": 5.0,
        "source_job": "prediction/source-job",
    }


def test_biotite_matching_reuses_only_accessible_exact_analysis_jobs(monkeypatch, tmp_path):
    from mn_molecular_design_studio import complex_analysis_ui
    from mn_molecular_design_studio import auth, job_catalog

    source = {"task_group": "prediction", "run_id": "source", "candidate_id": "candidate"}
    accessible_run = tmp_path / "analysis" / "accessible"
    private_run = tmp_path / "analysis" / "private"
    for run_dir, owner in ((accessible_run, USER["id"]), (private_run, OTHER_USER["id"])):
        (run_dir / "artifacts" / "analysis").mkdir(parents=True)
        (run_dir / "metadata.json").write_text(json.dumps({
            "upstream_task_group": "prediction", "upstream_run_id": "source", "status": "completed",
        }))
        manifest = {
            "tool": "biotite_interface", "settings": {"contact_distance": 5.0},
            "items": [{"id": "item_000000", "candidate_id": "candidate", "structure_digest": "sha256-structure",
                       "binder_chains": ["B"], "target_chains": ["A"]}],
        }
        (run_dir / "artifacts" / "analysis" / "batch_manifest.json").write_text(json.dumps(manifest))
        (run_dir / "artifacts" / "analysis" / "batch_status.json").write_text(json.dumps({
            "items": [{"id": "item_000000", "status": "failed", "error": "Selected chains were absent."}],
        }))
    rows = [
        {"task_group": "analysis", "run_id": "accessible", "run_dir": str(accessible_run), "status": "completed", "owner_id": USER["id"]},
        {"task_group": "analysis", "run_id": "private", "run_dir": str(private_run), "status": "completed", "owner_id": OTHER_USER["id"]},
    ]
    monkeypatch.setattr(job_catalog, "search_jobs", lambda **kwargs: {"items": rows})
    monkeypatch.setattr(job_catalog, "get_job_access_record", lambda group, run: {"owner_id": USER["id"] if run == "accessible" else OTHER_USER["id"]})
    monkeypatch.setattr(auth, "can_access_job", lambda user, record: record.get("owner_id") == user["id"])
    monkeypatch.setattr(auth, "accessible_project_ids", lambda user: [])

    matches = complex_analysis_ui._matching_biotite_jobs(
        USER, source_task_group=source["task_group"], source_run_id=source["run_id"],
        candidate_id=source["candidate_id"], structure_digest="sha256-structure",
        binder_chains=["B"], target_chains=["A"], contact_distance=5.0,
    )

    assert [item["run_id"] for item in matches] == ["accessible"]
    assert matches[0]["manifest"]["tool"] == "biotite_interface"
    assert matches[0]["status"] == "failed"
    assert matches[0]["error"] == "Selected chains were absent."


def test_database_workspace_handles_missing_metadata_and_denies_private_structure(tmp_path, monkeypatch):
    database = _isolated_database(tmp_path, monkeypatch)
    imported_record = molecule_imports.parse_structure(_complex_pdb("protein"), "private-complex.pdb")[0]
    imported = database.import_records(OTHER_USER, [imported_record], name="Private PPI")
    row = database.list_molecules(OTHER_USER, molecule_ids=imported["molecule_ids"])["items"][0]
    structure_id = row["representations"][0]["id"]
    from mn_molecular_design_studio import complex_analysis_ui, molecular_database

    monkeypatch.setattr(molecular_database, "MolecularDatabase", lambda: database)
    shown = []
    monkeypatch.setattr(complex_analysis_ui, "_show_structure", lambda *args, **kwargs: shown.append(args[0]))
    script = f'''\
import streamlit as st
from mn_molecular_design_studio.complex_analysis_ui import render_complex_analysis_workspace
st.session_state["complex-analysis-source-kind"] = "Database"
st.session_state["complex-analysis-open"] = {{"kind": "database", "molecule_id": {row['id']!r}}}
render_complex_analysis_workspace({USER!r})
'''
    app = AppTest.from_string(script, default_timeout=10).run()

    assert not app.exception
    assert not shown
    assert any("No accessible Database records" in item.value for item in app.info)
    owner_script = f'''\
import streamlit as st
from mn_molecular_design_studio.complex_analysis_ui import render_complex_analysis_workspace
st.session_state["complex-analysis-source-kind"] = "Database"
st.session_state["complex-analysis-open"] = {{"kind": "database", "molecule_id": {row['id']!r}}}
render_complex_analysis_workspace({OTHER_USER!r})
'''
    owner_app = AppTest.from_string(owner_script, default_timeout=10).run()
    assert not owner_app.exception
    assert any("No saved evaluation is associated" in item.value for item in owner_app.info)
    with pytest.raises(PermissionError):
        database.structure_file(USER, structure_id)


def test_completed_results_structures_follow_existing_job_permissions(tmp_path, monkeypatch):
    database = _isolated_database(tmp_path, monkeypatch)
    run_dir = tmp_path / "workspace" / "runs" / "boltz2" / "private-result"
    run_dir.mkdir(parents=True)
    (run_dir / "molecular_design_studio.json").write_text(json.dumps({"studio_owner_id": OTHER_USER["id"], "project_id": ""}))
    (run_dir / "metadata.json").write_text(json.dumps({"tool": "boltz2", "status": "completed", "campaign_name": "Private prediction"}))
    (run_dir / "input.json").write_text(json.dumps({"tool": "boltz2", "params": {}}))
    (run_dir / "result.json").write_text(json.dumps({"success": True, "metrics": {}}))
    structure_path = run_dir / "complex.pdb"
    structure_path.write_bytes(_complex_pdb("protein"))
    from mn_molecular_design_studio import complex_analysis_ui
    from mn_molecular_design_studio.job_catalog import reconcile

    reconcile()
    assert complex_analysis_ui._accessible_result_jobs(USER) == []
    visible = complex_analysis_ui._accessible_result_jobs(OTHER_USER)
    assert len(visible) == 1
    assert complex_analysis_ui._job_structure_options(USER, visible[0]) == []
    options = complex_analysis_ui._job_structure_options(OTHER_USER, visible[0])
    assert [item["path"] for item in options] == [structure_path]
    saved = database.save_complex_review(
        OTHER_USER, source_kind="results", subject_id="boltz2/private-result",
        structure_ref="complex.pdb", status="Reviewed", note="Manual review note.",
    )
    assert saved["note"] == "Manual review note."
    with pytest.raises(PermissionError):
        database.get_complex_review(
            USER, source_kind="results", subject_id="boltz2/private-result", structure_ref="complex.pdb",
        )


def test_database_review_annotations_persist_separately_and_follow_record_access(tmp_path, monkeypatch):
    database = _isolated_database(tmp_path, monkeypatch)
    imported_record = molecule_imports.parse_structure(_complex_pdb("dna"), "reviewed-complex.pdb")[0]
    imported = database.import_records(USER, [imported_record], name="Reviewed complex")
    record = database.list_molecules(USER, molecule_ids=imported["molecule_ids"])["items"][0]
    structure_id = record["representations"][0]["id"]
    before_observations = list(record["observations"])

    saved = database.save_complex_review(
        USER, source_kind="database", subject_id=record["id"], structure_ref=structure_id,
        status="Needs review", note="Check source numbering before interpreting this model.",
    )
    reloaded = MolecularDatabase(database.root).get_complex_review(
        USER, source_kind="database", subject_id=record["id"], structure_ref=structure_id,
    )

    assert saved["status"] == reloaded["status"] == "Needs review"
    assert saved["note"] == reloaded["note"]
    assert reloaded["updated_at"]
    assert database.list_molecules(USER, molecule_ids=[record["id"]])["items"][0]["observations"] == before_observations
    with pytest.raises(PermissionError):
        database.get_complex_review(
            OTHER_USER, source_kind="database", subject_id=record["id"], structure_ref=structure_id,
        )
    with pytest.raises(ValueError):
        database.save_complex_review(
            USER, source_kind="database", subject_id=record["id"], structure_ref=structure_id,
            status="Predicted binder", note="Not a supported review status.",
        )


def test_workspace_review_form_saves_persistent_annotation(tmp_path, monkeypatch):
    database = _isolated_database(tmp_path, monkeypatch)
    imported_record = molecule_imports.parse_structure(_complex_pdb("protein"), "review-form.pdb")[0]
    imported = database.import_records(USER, [imported_record], name="Review form complex")
    record = database.list_molecules(USER, molecule_ids=imported["molecule_ids"])["items"][0]
    structure_id = record["representations"][0]["id"]
    token = hashlib.sha256("\0".join(("database", record["id"], structure_id, "")).encode()).hexdigest()[:16]
    from mn_molecular_design_studio import complex_analysis_ui, molecular_database

    monkeypatch.setattr(molecular_database, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(complex_analysis_ui, "_show_structure", lambda *args, **kwargs: None)
    script = f'''\
import streamlit as st
from mn_molecular_design_studio.complex_analysis_ui import render_complex_analysis_workspace
st.session_state["complex-analysis-source-kind"] = "Database"
st.session_state["complex-analysis-open"] = {{"kind": "database", "molecule_id": {record['id']!r}, "structure_id": {structure_id!r}}}
render_complex_analysis_workspace({USER!r})
'''
    app = AppTest.from_string(script, default_timeout=10).run()
    app.selectbox(key=f"complex-analysis-review-status-{token}").set_value("Reviewed").run()
    app.text_area(key=f"complex-analysis-review-note-{token}").set_value("Human review completed.").run()
    app.button(key=f"complex-analysis-review-save-{token}").click().run()

    assert not app.exception
    saved = database.get_complex_review(
        USER, source_kind="database", subject_id=record["id"], structure_ref=structure_id,
    )
    assert saved["status"] == "Reviewed"
    assert saved["note"] == "Human review completed."
