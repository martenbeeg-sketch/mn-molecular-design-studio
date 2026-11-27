import json
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio import dataset_modelling_ui
from mn_molecular_design_studio.engine.workflow_engine.workflows.ligand_tools import inspect_complex_gaps
from mn_molecular_design_studio.protein_modelling import loop_residues, observed_chains, parse_options


PUBLIC_PDB = Path(dataset_modelling_ui.__file__).parent / "assets/inverse_folding/1CRN.pdb"


@pytest.fixture
def gapped_pdb(tmp_path):
    source = tmp_path / "missing-residue.pdb"
    source.write_text("\n".join(line for line in PUBLIC_PDB.read_text().splitlines()
                               if not (line.startswith("ATOM  ") and line[21:26] == "A  21")) + "\n")
    return source


def test_modelling_requires_confirmation_and_a_connected_loop(gapped_pdb):
    gap, = [row for row in inspect_complex_gaps(gapped_pdb) if row["author_start"] == 21]
    assert gap["suggestion"] == "T"
    selection = {"chain": "A", "left_resseq": 20, "right_resseq": 22, "sequence": "T"}
    with pytest.raises(ValueError, match="Confirm"):
        parse_options(json.dumps({"gaps": [selection]}))
    assert parse_options(json.dumps({"gaps": [{**selection, "confirmed": True}]}))["gaps"][0]["sequence"] == "T"
    loop = {"chain": "A", "start": "20", "end": "22"}
    assert loop_residues(observed_chains(PUBLIC_PDB.read_text()), loop) == ["20", "21", "22"]
    with pytest.raises(ValueError, match="backbone break"):
        loop_residues(observed_chains(gapped_pdb.read_text()), loop)


def test_dataset_gap_editor_queues_only_confirmed_input(gapped_pdb, monkeypatch):
    expander = dataset_modelling_ui.st.expander
    monkeypatch.setattr(dataset_modelling_ui.st, "expander", lambda label, **kwargs: expander(label, expanded=True, **kwargs))
    app = AppTest.from_string(
        "import streamlit as st\nfrom pathlib import Path\n"
        "from mn_molecular_design_studio.dataset_modelling_ui import render_modelling\n"
        f"source = Path({str(gapped_pdb)!r})\n"
        "def capture(schema, workflow, parameters, project, settings, label):\n"
        "    st.session_state['submitted'] = (workflow, parameters, project)\n"
        "    return {'job': {'run_id': 'test-modelling-job'}}\n"
        "render_modelling(source, source.read_text(), 'test', [{'id':'p1','name':'Project'}], ['p1'], capture)",
        default_timeout=30,
    ).run()
    assert not app.exception
    assert app.button(key="test-run-modeller").disabled
    gap = next(value for value in app.multiselect[0].options if "G20–P22" in value)
    app.multiselect[0].set_value([gap]).run()
    assert app.text_input[0].value == "T"
    assert app.button(key="test-run-modeller").disabled
    app.checkbox[0].set_value(True).run()
    app.button(key="test-run-modeller").click().run()
    assert not app.exception
    workflow, parameters, project = app.session_state["submitted"]
    assert workflow == "ligand.complex.prepare" and project == "p1"
    assert json.loads(parameters["modelling_json"]) == {"gaps": [
        {"chain": "A", "left_resseq": 20, "right_resseq": 22, "sequence": "T", "confirmed": True}
    ]}
