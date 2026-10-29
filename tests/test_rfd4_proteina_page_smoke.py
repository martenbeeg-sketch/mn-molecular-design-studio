"""Streamlit page smoke coverage for the RFD4-Proteina setup form."""
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio import rfd4_proteina_input as inputs


REPO = Path(__file__).resolve().parents[1]
TARGET = REPO / "mn_molecular_design_studio/assets/proteina_complexa/PD-L1.pdb"


def page(*, folding=False, target=False):
    source = f'''from pathlib import Path
import streamlit as st
from mn_molecular_design_studio import rfd4_proteina_ui as ui
ui.user_input_root=lambda owner:Path('/tmp/rfd4-proteina-page-smoke')/owner
'''
    if target:
        source += f"st.session_state['rfd4-target-fetched-target']={str(TARGET)!r}\n"
    source += f'''ui.render_setup({{}},owner_id='page-smoke-user',
    project_selector=lambda _: 'project-smoke',submit=lambda *args:None,folding={folding!r})
'''
    return AppTest.from_string(source).run()


def pills(app, key: str):
    """Find a Streamlit pill group through AppTest's button-group query."""
    for element in app.get("button_group"):
        if element.key == key:
            return element
    raise AssertionError(f"No pills widget found with key {key!r}")


@pytest.mark.parametrize("task", [*inputs.TASKS, "idr_peptide_binder", "cyclic_peptide_binder"])
def test_each_design_task_page_renders(task):
    app = page()
    assert not app.exception, f"initial render: {app.exception}"
    pills(app, "rfd4-task-pills").set_value(task).run()
    assert not app.exception, f"{task}: {app.exception}"
    assert pills(app, "rfd4-task-pills").value == task
    validate = app.button(key="rfd4-validate-prepare")
    submit = app.button(key="rfd4-submit")
    if task in {"idr_peptide_binder", "cyclic_peptide_binder"}:
        assert validate.disabled and submit.disabled
    else:
        # Tasks that need an uploaded structure, sequence, or motif should
        # still render their controls while keeping actions disabled.
        assert isinstance(validate.disabled, bool)


def test_idr_sequence_only_form_submits_native_condition_spec(tmp_path):
    source = f'''from pathlib import Path
import streamlit as st
from mn_molecular_design_studio import rfd4_proteina_ui as ui
ui.user_input_root=lambda owner:Path({str(tmp_path)!r})/owner/'inputs'
st.session_state['submitted']=[]
def capture(*args): st.session_state['submitted'].append(args)
ui.render_setup({{}},owner_id='idr-page-user',project_selector=lambda _: 'project-idr',submit=capture)
'''
    app = AppTest.from_string(source).run()
    pills(app, "rfd4-task-pills").set_value("idr_peptide_binder").run()
    app.text_area(key="rfd4-idr-preview-fasta").set_value(">target_idr\nACDEFGHIKLMNPQRSTVWY").run()
    app.text_input(key="rfd4-idr-preview-hotspots").set_value("A5, A8-12").run()
    assert not app.exception
    assert not app.button(key="rfd4-validate-prepare").disabled
    assert not app.button(key="rfd4-submit").disabled
    app.button(key="rfd4-submit").click().run()
    assert not app.exception
    _schema, workflow, parameters, project, _options, _title = app.session_state["submitted"][0]
    assert workflow == "design.rfd4_proteina"
    assert project == "project-idr"
    assert parameters["task"] == "protein_binder"
    assert parameters["input_mode"] == "condition_spec"
    assert parameters["target_chains"] == ["A"]
    assert parameters["target_structure"] == parameters["input_files"][0]
    assert "C_SEQ" in parameters["condition_spec"]
    assert "C_HOT" in parameters["condition_spec"]
    assert "C_CRD:" not in parameters["condition_spec"]
    assert "C_CRDGRP:" not in parameters["condition_spec"]


def test_idr_condition_spec_builder_keeps_sequence_and_hotspots_without_coordinate_conditions():
    import yaml

    from mn_molecular_design_studio import rfd4_proteina_input as idr_inputs

    sequence = idr_inputs.parse_idr_fasta(">idR target\nACDEFGHIKLMNPQRSTVWY")
    spec = yaml.safe_load(idr_inputs.idr_binder_condition_spec(
        sequence, binder_length_min=70, binder_length_max=110, hotspots="A5,A8-12"
    ))["idr_peptide_binder"]
    assert spec["input"] == "target.cif"
    assert spec["contig"] == "70-110, /0, A1-20"
    assert spec["conditions"]["C_SEQ"][True] == [{"select": "sel('A')"}]
    assert spec["conditions"]["C_HOT"][True] == [
        {"select": "sel('A/*/5')"}, {"select": "sel('A/*/[8-12]')"}
    ]
    assert not ({"C_CRD", "C_CRDGRP"} & set(spec["conditions"]))


def test_idr_sequence_carrier_uses_extended_backbone_geometry(tmp_path):
    import numpy as np

    from mn_molecular_design_studio import rfd4_proteina_input as idr_inputs
    from mn_molecular_design_studio.design_validation import protein_chains, read_structure, sequence

    path = idr_inputs.write_idr_sequence_carrier("ACDEFGHIKL", tmp_path / "target.cif")
    structure = read_structure(path)
    chains = protein_chains(structure)
    assert sequence(chains["A"]) == "ACDEFGHIKL"
    residues = list(chains["A"])
    assert len(residues) == 10
    for left, right in zip(residues, residues[1:]):
        assert np.linalg.norm(left["C"].coord - right["N"].coord) == pytest.approx(1.329, abs=0.01)


@pytest.mark.parametrize("operation", list(inputs.OPERATIONS))
def test_existing_structure_operations_render(operation):
    app = page(target=True)
    pills(app, "rfd4-task-pills").set_value("existing_structure").run()
    pills(app, "rfd4-operation").set_value(operation).run()
    assert not app.exception, f"{operation}: {app.exception}"
    assert pills(app, "rfd4-operation").value == operation
    assert app.button(key="rfd4-submit").disabled is False


@pytest.mark.parametrize("framework", inputs.FRAMEWORKS)
def test_nanobody_framework_choices_render(framework):
    app = page(target=True)
    pills(app, "rfd4-task-pills").set_value("nanobody").run()
    pills(app, "rfd4-framework").set_value(framework).run()
    assert not app.exception, f"{framework}: {app.exception}"
    assert pills(app, "rfd4-framework").value == framework
    if framework == "Your own nanobody":
        app.multiselect(key="rfd4-design-chains-nanobody").set_value(["A"]).run()
        assert app.text_input(key="rfd4-cdr-nanobody").value == ""
    else:
        assert any("built-in framework" in item.value for item in app.caption)


def test_folding_page_switches_between_sequence_and_structure():
    app = page(folding=True, target=True)
    assert not app.exception
    assert app.radio[0].value == "Protein sequence"
    app.text_area[0].set_value("ACDEFGHIK").run()
    assert not app.exception
    assert not app.button(key="rfd4-submit").disabled
    app.radio[0].set_value("Structure").run()
    assert not app.exception
    assert any(widget.label == "Enter PDB Code" for widget in app.text_input)
    assert not app.button(key="rfd4-submit").disabled


def test_load_example_button_populates_page_inputs(tmp_path):
    source = f'''from pathlib import Path
from mn_molecular_design_studio import rfd4_proteina_ui as ui
ui.user_input_root=lambda owner:Path({str(tmp_path)!r})/owner/'inputs'
ui.render_setup({{}},owner_id='page-smoke-user',project_selector=lambda _: 'project-smoke',submit=lambda *args:None)
'''
    app = AppTest.from_string(source).run()
    app.button(key="rfd4-load-example").click().run()
    assert not app.exception
    assert pills(app, "rfd4-task-pills").value == "protein_binder"
    assert app.multiselect(key="rfd4-target-chains-protein_binder").value == ["A"]
    assert app.text_input(key="rfd4-hotspots-protein_binder").value == "A40,A99,A107"
    assert app.number_input(key="rfd4-min-protein_binder").value == 60
    assert app.number_input(key="rfd4-max-protein_binder").value == 120
    assert not app.button(key="rfd4-submit").disabled


def test_validate_prepare_button_calls_prepare_only_validator(monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.workflows import rfd4_proteina_prepare

    captured = []
    monkeypatch.setattr(rfd4_proteina_prepare, "validate_inputs", lambda payload: captured.append(payload) or {
        "count": 1, "items": [{"name": "prepared input", "atoms": 18,
                                 "designed_chains": [], "target_chains": [], "warnings": []}],
    })
    app = page()
    pills(app, "rfd4-task-pills").set_value("unconditional").run()
    app.button(key="rfd4-validate-prepare").click().run()
    assert not app.exception
    assert len(captured) == 1
    assert captured[0]["task"] == "unconditional"
    assert any("input preparation passed" in item.value for item in app.success)


def test_submit_button_sends_design_payload_and_project():
    app = AppTest.from_string('''import streamlit as st
from mn_molecular_design_studio import rfd4_proteina_ui as ui
st.session_state['submitted']=[]
def capture(*args): st.session_state['submitted'].append(args)
ui.render_setup({},owner_id='page-smoke-user',project_selector=lambda _: 'project-smoke',submit=capture)
''').run()
    pills(app, "rfd4-task-pills").set_value("unconditional").run()
    app.button(key="rfd4-submit").click().run()
    assert not app.exception
    assert len(app.session_state["submitted"]) == 1
    schema, workflow, payload, project, options, title = app.session_state["submitted"][0]
    assert workflow == "design.rfd4_proteina"
    assert payload["task"] == "unconditional"
    assert payload["length_min"] == 100 and payload["length_max"] == 100
    assert project == "project-smoke"
    assert title == "RFD4-Proteina"


@pytest.mark.parametrize("use_structure", [False, True])
def test_folding_submission_payload_tracks_selected_input_mode(use_structure):
    source = f'''from pathlib import Path
import streamlit as st
from mn_molecular_design_studio import rfd4_proteina_ui as ui
ui.user_input_root=lambda owner:Path('/tmp/rfd4-proteina-page-smoke')/owner
st.session_state['rfd4-target-fetched-target']={str(TARGET)!r}
st.session_state['submitted']=[]
def capture(*args): st.session_state['submitted'].append(args)
ui.render_setup({{}},owner_id='page-smoke-user',project_selector=lambda _: 'project-smoke',
                submit=capture,folding=True)
'''
    app = AppTest.from_string(source).run()
    if use_structure:
        app.radio[0].set_value("Structure").run()
    else:
        app.text_area[0].set_value("ACDEFGHIK").run()
    assert not app.button(key="rfd4-submit").disabled
    app.button(key="rfd4-submit").click().run()
    assert not app.exception
    schema, workflow, payload, project, options, title = app.session_state["submitted"][0]
    assert workflow == "refolding.rfd4_proteina"
    assert project == "project-smoke"
    assert payload["options"]["operation"] == "structure_prediction"
    if use_structure:
        assert payload["target_structure"] == str(TARGET)
        assert "protein_sequences" not in payload
    else:
        assert payload["protein_sequences"] == "ACDEFGHIK"
        assert "target_structure" not in payload


def test_results_table_separates_designed_and_target_residue_counts():
    app = AppTest.from_string('''from pathlib import Path
from mn_molecular_design_studio import rfd4_proteina_ui as ui
ui.render_results('run-1',Path('/tmp/nonexistent-rfd4-run'),[{
    'candidate_id':'design-1','binder_sequence':'ACDE','binder_chains':['B'],
    'target_chains':['A'],'metrics':{'generated_residue_count':24},
}],{}, {})
''').run()
    assert not app.exception
    frame = app.dataframe[0].value
    assert list(frame.columns) == ["Design", "Sequence", "Designed Chains", "Designed Residues", "Target Chains", "Total Residues"]
    assert frame.iloc[0]["Designed Chains"] == "B"
    assert frame.iloc[0]["Designed Residues"] == 4
    assert frame.iloc[0]["Target Chains"] == "A"
    assert frame.iloc[0]["Total Residues"] == 24
