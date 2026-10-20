from __future__ import annotations

import uuid

import pytest
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio import docking_preparations as preparations


OWNER = "a" * 32
WORKFLOW = "ligand.docking.vina"


@pytest.fixture(autouse=True)
def isolated_preparations(tmp_path, monkeypatch):
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path))


def test_queued_prediction_and_preparation_survive_a_new_session():
    preparation_id = str(uuid.uuid4())
    state = {
        "ligand-vina-job-name": "Pocket example", "vina-target-path": "/saved/target.pdb",
        "vina-compound-input": {"path": "/saved/compounds.json"},
        "vina-boltz-selected-chains": ["B"], "vina-box-size-x": 25.0,
        "vina-pending-prediction": {
            "job": {"cluster_job_id": "queued-id", "status": "waiting_for_workstation"},
            "selected_chains": ["B"], "project_id": "project-1", "target_hash": "target",
        },
        "studio-session-token": "must-not-save", "vina-protein-upload-0": object(),
        "ligand-vina-run": True,
    }
    preparations.save(OWNER, "project-1", WORKFLOW, preparation_id, state)
    restored = preparations.list_saved(OWNER, "project-1", WORKFLOW)[0]
    assert restored["state"] == {key: value for key, value in state.items() if key not in {
        "studio-session-token", "vina-protein-upload-0", "ligand-vina-run",
    }}
    assert restored["state"]["vina-pending-prediction"]["job"]["cluster_job_id"] == "queued-id"
    assert preparations.list_saved(OWNER, "project-2", WORKFLOW) == []
    assert preparations.list_saved("b" * 32, "project-1", WORKFLOW) == []
    assert preparations.list_saved(OWNER, "project-1", "ligand.docking.gnina") == []


@pytest.mark.parametrize("project", ["../escape", "/absolute", "", "a/b"])
def test_preparation_paths_cannot_escape_the_project(project):
    with pytest.raises(ValueError, match="valid project"):
        preparations.save(OWNER, project, WORKFLOW, str(uuid.uuid4()), {})


def test_preparation_picker_restores_form_after_navigation_and_new_browser():
    script = f'''
import streamlit as st
from mn_molecular_design_studio.vina_ui import _render_preparation_picker, _save_current_preparation
page = st.radio("Page", ["Docking", "Results"])
if page == "Docking":
    _render_preparation_picker({OWNER!r}, "project-1", {WORKFLOW!r})
    st.text_input("Job", value="default", key="ligand-vina-job-name")
    st.number_input("Width", value=20.0, key="vina-box-size-x")
    _save_current_preparation()
'''
    app = AppTest.from_string(script).run()
    assert not app.exception
    app.text_input(key="ligand-vina-job-name").set_value("Saved receptor").run()
    app.number_input(key="vina-box-size-x").set_value(33.0).run()
    app.radio[0].set_value("Results").run()
    app.radio[0].set_value("Docking").run()
    assert not app.exception
    assert app.text_input(key="ligand-vina-job-name").value == "Saved receptor"
    assert app.number_input(key="vina-box-size-x").value == 33.0
    reopened = AppTest.from_string(script).run()
    assert not reopened.exception
    assert reopened.text_input(key="ligand-vina-job-name").value == "Saved receptor"
    assert reopened.number_input(key="vina-box-size-x").value == 33.0
    reopened.button(key="docking-new-preparation").click().run()
    assert not reopened.exception
    assert reopened.text_input(key="ligand-vina-job-name").value == "default"
    assert len(preparations.list_saved(OWNER, "project-1", WORKFLOW)) == 2


def test_waiting_cluster_job_is_remembered_without_a_run_directory():
    app = AppTest.from_string(f'''
import streamlit as st
from mn_molecular_design_studio.vina_ui import _render_preparation_picker, _remember_boltz_prediction
_render_preparation_picker({OWNER!r}, "project-1", {WORKFLOW!r})
st.session_state["last_submitted_job"] = {{"cluster_job_id": "waiting-id", "status": "waiting_for_workstation"}}
_remember_boltz_prediction({{"project_id": "project-1", "target_hash": "target", "selected_chains": ["B"]}})
''').run()
    assert not app.exception
    saved = preparations.list_saved(OWNER, "project-1", WORKFLOW)[0]
    assert saved["state"]["vina-pending-prediction"]["job"]["cluster_job_id"] == "waiting-id"


def test_live_prediction_status_refreshes_and_applies_completed_contacts(monkeypatch):
    from mn_molecular_design_studio import vina_inputs

    monkeypatch.setattr(vina_inputs, "predicted_contact_residues", lambda *args, **kwargs: {("B", 448), ("B", 450)})
    monkeypatch.setattr(vina_inputs, "predicted_ligand_overlay_pdb", lambda *args, **kwargs: "HETATM PREDICTED\nEND\n")
    monkeypatch.setattr(vina_inputs, "protein_sequence_entries", lambda *args: {"B": [(448, "V"), (450, "E")]})
    script = '''
import streamlit as st
from mn_molecular_design_studio.vina_ui import _render_prediction_progress
def state(job):
    completed = st.session_state.get("test-completed", False)
    return {"job": {**job, "status": "completed" if completed else "running"},
            "workstation": "pc002975", "structures": ["PREDICTED"] if completed else []}
pending = st.session_state.setdefault("vina-pending-prediction", {
    "job": {"cluster_job_id": "job-id", "job_code": "99072"}, "selected_chains": ["B"]})
_render_prediction_progress(pending, "TARGET", 5.0, state)
'''
    app = AppTest.from_string(script).run()
    assert not app.exception
    assert any("Running on pc002975" in item.value for item in app.info)
    app.session_state["test-completed"] = True
    app.run()
    assert not app.exception
    assert any("Completed on pc002975" in item.value for item in app.success)
    assert app.session_state["ligand-vina-selected-residues"] == [("B", 448), ("B", 450)]
    assert app.session_state["vina-prediction-overlay-pdb"].startswith("HETATM")


def test_remote_prediction_structure_paths_are_scoped_to_the_original_run():
    job = {"node_id": "node-a", "task_group": "refolding", "run_id": "run-1", "location": "/remote/runs/refolding/run-1"}
    relative = "artifacts/raw/prediction.pdb"
    assert preparations.artifact_relative_path("runs:///refolding/run-1/" + relative, job) == relative
    assert preparations.artifact_relative_path("cluster://node-a/refolding/run-1/" + relative, job) == relative
    assert preparations.artifact_relative_path(job["location"] + "/" + relative, job) == relative
    for invalid in ["../secret.pdb", "cluster://node-b/refolding/run-1/" + relative,
                    "runs:///refolding/other/" + relative, "/remote/other/file.pdb", "metadata.json"]:
        with pytest.raises(ValueError):
            preparations.artifact_relative_path(invalid, job)


def test_full_vina_page_resumes_waiting_prediction_and_uses_its_finished_pocket():
    script = f'''
import hashlib
import streamlit as st
from mn_molecular_design_studio import vina_inputs
from mn_molecular_design_studio.vina_ui import render_docking_setup
if not st.session_state.get("test-initialized"):
    target, compounds = vina_inputs.stage_example({OWNER!r}, project_id="project-1")
    st.session_state.update({{
        "test-initialized": True, "vina-target-path": str(target), "vina-compound-input": compounds,
        "vina-default-target-signature": hashlib.sha1(target.read_bytes()).hexdigest()[:12],
        "ligand-vina-pocket-mode": "Predict (Boltz)", "vina-boltz-selected-chains": ["A"],
        "ligand-vina-job-name": "autodock-vina-example",
    }})
def submit(schema, workflow, parameters, project, settings, title):
    st.session_state["last_submitted_job"] = {{
        "cluster_job_id": "queued-before-assignment", "job_code": "TEST1", "project_id": project,
        "tool": "optimized_prediction_boltz2", "status": "waiting_for_workstation",
        "campaign_name": parameters["job_name"],
    }}
def status(job):
    from pathlib import Path
    complete = st.session_state.get("test-complete", False)
    return {{"job": {{**job, "status": "completed" if complete else "running"}},
             "workstation": "pc002975", "structures": [Path(st.session_state["vina-target-path"]).read_text()] if complete else []}}
render_docking_setup({WORKFLOW!r}, {{}}, owner_id={OWNER!r}, compound_sets=[],
    project_selector=lambda key: "project-1", submit=submit, title="AutoDock Vina", prediction_state=status)
'''
    app = AppTest.from_string(script, default_timeout=20).run()
    assert not app.exception
    app.button(key="vina-predict-pocket").click().run()
    assert not app.exception
    assert any("Running on pc002975" in item.value for item in app.info)
    restored = AppTest.from_string(script, default_timeout=20).run()
    assert not restored.exception
    assert any("Job TEST1" in item.value for item in restored.info)
    assert restored.multiselect(key="vina-boltz-selected-chains").value == ["A"]
    restored.session_state["test-complete"] = True
    restored.run()
    assert not restored.exception
    assert any("Completed on pc002975" in item.value for item in restored.success)
    assert not restored.button(key="ligand-vina-run").disabled
    assert restored.session_state["vina-applied-prediction"]
