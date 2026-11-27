import json
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio import dataset_modelling_ui
from mn_molecular_design_studio.engine.workflow_engine.workflows.ligand_tools import inspect_complex_gaps
from mn_molecular_design_studio.protein_modelling import loop_residues, observed_chains, parse_options, substitution_sites


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
        parse_options({"gaps": [selection]})
    assert parse_options({"gaps": [{**selection, "confirmed": True}]})["gaps"][0]["sequence"] == "T"
    loop = {"chain": "A", "start": "20", "end": "22"}
    assert loop_residues(observed_chains(PUBLIC_PDB.read_text()), loop) == ["20", "21", "22"]
    with pytest.raises(ValueError, match="backbone break"):
        loop_residues(observed_chains(gapped_pdb.read_text()), loop)


def test_dataset_gap_editor_queues_only_confirmed_input(gapped_pdb, monkeypatch):
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(gapped_pdb.parent / "workspace"))
    app = AppTest.from_string(
        "import streamlit as st\nfrom pathlib import Path\n"
        "from mn_molecular_design_studio.dataset_modelling_ui import render_protein_tools\n"
        "from mn_molecular_design_studio.preparation_history import extend_history, sha256\n"
        f"source = Path({str(gapped_pdb)!r})\n"
        "def capture(schema, workflow, parameters, project, settings, label):\n"
        "    st.session_state['submitted'] = (workflow, parameters, project)\n"
        "    return {'job': {'run_id': 'test-modelling-job'}}\n"
        "content = source.read_bytes()\n"
        "history = extend_history(None, sha256(content), sha256(content), [{'label':'Trim structure'}])\n"
        "render_protein_tools({'id':'a' * 32}, source, content, history, 'test', [{'id':'p1','name':'Project'}], ['p1'], capture, viewer_key='viewer')",
        default_timeout=30,
    ).run()
    assert not app.exception
    app.selectbox(key="test-operation").set_value("Rebuild missing residues or a loop").run()
    assert app.button(key="test-run-preparation").disabled
    gap = next(value for value in app.multiselect[0].options if "G20–P22" in value)
    app.session_state["viewer"] = _pick(1, 20, "GLY")
    app.run()
    assert app.multiselect[0].value == [gap]
    assert app.text_input[0].value == "T"
    assert app.button(key="test-run-preparation").disabled
    app.checkbox[0].set_value(True).run()
    app.button(key="test-run-preparation").click().run()
    assert not app.exception
    workflow, parameters, project = app.session_state["submitted"]
    assert workflow == "ligand.complex.prepare" and project == "p1"
    assert parameters["modelling_options"] == {"gaps": [
        {"chain": "A", "left_resseq": 20, "right_resseq": 22, "sequence": "T", "confirmed": True}
    ]}


def _pick(event, number, name, insertion=""):
    return json.dumps({"residuePick": {
        "eventId": event, "chainId": "A", "residueNumber": number,
        "insertionCode": insertion, "residueName": name,
    }})


def test_viewer_picks_preserve_author_identity_and_require_confirmation_to_trim():
    app = AppTest.from_string(
        "import streamlit as st\nfrom pathlib import Path\n"
        "from mn_molecular_design_studio.dataset_modelling_ui import render_protein_tools\n"
        f"source = Path({str(PUBLIC_PDB)!r})\n"
        "st.slider('Trim', 1, 46, (1, 46), key='test-A')\n"
        "render_protein_tools({}, source, source.read_bytes(), {}, 'test', [{'id':'p1','name':'Project'}],"
        " ['p1'], None, viewer_key='viewer')",
    ).run()
    app.selectbox(key="test-operation").set_value("Substitute amino acids").run()
    app.session_state["viewer"] = _pick(1, 9, "ALA", insertion="A")
    app.run()
    assert app.multiselect(key="test-A-substitutions").value == []
    app.session_state["viewer"] = _pick(2, 9, "ALA")
    app.run()
    assert app.multiselect(key="test-A-substitutions").value == ["9"]
    app.run()
    assert app.multiselect(key="test-A-substitutions").value == ["9"]
    app.session_state["viewer"] = _pick(3, 9, "ALA")
    app.run()
    assert app.multiselect(key="test-A-substitutions").value == []

    app.selectbox(key="test-operation").set_value("Remodel an existing loop").run()
    app.session_state["viewer"] = _pick(4, 9, "ALA")
    app.run()
    app.session_state["viewer"] = _pick(5, 7, "ILE")
    app.run()
    assert app.selectbox(key="test-loop-A-start").value == "7"
    assert app.selectbox(key="test-loop-A-end").value == "9"

    app.selectbox(key="test-operation").set_value("Extend the C-terminus").run()
    app.text_input(key="test-extension-sequence").set_value("G").run()
    app.session_state["viewer"] = _pick(6, 9, "ALA")
    app.run()
    assert app.slider(key="test-A").value == (1, 46)
    assert app.button(key="test-run-preparation").disabled
    app.button(key="test-apply-attachment").click().run()
    assert not app.exception
    assert app.slider(key="test-A").value == (1, 9)


def test_amino_acid_substitution_uses_an_exact_observed_residue(gapped_pdb):
    change = {"chain": "A", "residue": "21", "amino_acid": "a"}
    options = parse_options({"substitutions": [change]})
    site, = substitution_sites(observed_chains(PUBLIC_PDB.read_text()), options["substitutions"])
    assert (site["resseq"], site["resname"], site["target"]) == ("21", "THR", "ALA")
    with pytest.raises(ValueError, match="observed backbone coordinates"):
        substitution_sites(observed_chains(gapped_pdb.read_text()), options["substitutions"])
    with pytest.raises(ValueError, match="single substitution"):
        parse_options({"substitutions": [change, change]})
