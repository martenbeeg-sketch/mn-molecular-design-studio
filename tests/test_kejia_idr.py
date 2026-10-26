import json
from pathlib import Path

import pytest

from mn_molecular_design_studio import kejia_idr_input as inputs
from mn_molecular_design_studio.design_validation import protein_chains, read_structure, sequence
from mn_molecular_design_studio.engine.workflow_engine.workflows.kejia_idr import normalize_output, stage_template
from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import node_request
from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import workflow_schema
from mn_molecular_design_studio.pipeline_graph import new_graph, validate_graph


@pytest.fixture
def template(tmp_path):
    path = tmp_path / "two_chains.pdb"
    lines, serial = [], 0
    for chain in ("X", "Y"):
        for residue in range(11, 15):
            for atom, element in (("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")):
                serial += 1
                lines.append(f"ATOM  {serial:5d} {atom:^4s} GLY {chain}{residue:4d}    {float(serial):8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{20.:6.2f}          {element:>2s}")
        lines.append("TER")
    path.write_text("\n".join(lines) + "\nEND\n")
    return path


def test_single_fasta_and_whitespace_are_normalized():
    assert inputs.normalize_sequence(">fragment\n fncr ewcwn\n") == inputs.EXAMPLE_SEQUENCE


@pytest.mark.parametrize("value", ["", "ACD", "A" * 41, "ACDX", ">one\nAAAA\n>two\nGGGG", "ACD*"])
def test_invalid_target_fragments_are_rejected(value):
    with pytest.raises(ValueError):
        inputs.normalize_sequence(value)


def test_template_roles_and_cpu_settings_are_checked(template):
    request = dict(target_sequence="GGGG", template_structure=template, binder_chain="X", peptide_chain="Y")
    assert inputs.prepare_request(request)["templates"][0]["binder_chain"] == "X"
    with pytest.raises(ValueError, match="distinct"):
        inputs.prepare_request(request | {"peptide_chain": "X"})
    with pytest.raises(ValueError, match="finite"):
        inputs.prepare_request(request | {"max_farep": float("nan")})
    with pytest.raises(ValueError, match="integer"):
        inputs.prepare_request(request | {"flex_docking_repeats": True})
    with pytest.raises(ValueError, match="not both"):
        inputs.prepare_request(request | {"template_library": "walle"})


def test_staging_preserves_chains_and_fragment_metadata(template, tmp_path):
    output = tmp_path / "template_00001"
    output.mkdir()
    staged = output / "candidate.pdb"
    stage_template(template, staged, "X", "Y")
    chains = protein_chains(read_structure(staged))
    assert list(chains) == ["A", "B"]
    assert sequence(chains["B"]) == "GGGG"
    assert [residue.id[1] for residue in chains["B"]] == [1, 2, 3, 4]
    record = {"tag": "candidate", "pdb": staged.name, "target_sequence": "GGGG", "scores": {
        "target_peptide_start": "2", "target_peptide_end": "5", "delta_fa_rep": 1.2,
    }}
    (output / "threading.json").write_text(json.dumps({"records": [record]}))
    rows = normalize_output(output, "AGGGGA", {"path": str(template)})
    assert rows[0]["raw_metadata"]["target_fragment_range"] == [2, 5]
    assert rows[0]["binder_sequence"] is None
    assert set(protein_chains(read_structure(Path(rows[0]["target_pdb"])))) == {"B"}
    with pytest.raises(ValueError, match="does not match"):
        normalize_output(output, "AAAAAA", {})
    record["pdb"] = "../outside.pdb"
    (output / "threading.json").write_text(json.dumps({"records": [record]}))
    with pytest.raises(ValueError, match="invalid structure path"):
        normalize_output(output, "AGGGGA", {})


def test_pipeline_starter_and_registered_workflow(template):
    graph = new_graph("kejia_idr")
    assert [node["data"]["kind"] for node in graph["nodes"]] == ["kejia_idr", "ligandmpnn", "af2", "results"]
    validate_graph(graph)
    with pytest.raises(ValueError, match="4–40"):
        validate_graph(graph, require_inputs=True)
    graph["nodes"][0]["data"]["params"].update(target_sequence="GGGG", template_structure=str(template), binder_chain="X", peptide_chain="Y")
    validate_graph(graph, require_inputs=True)
    payload = node_request("kejia_idr", graph["nodes"][0]["data"]["params"], {}, "example")
    assert payload["workflow"] == "design.kejia_idr"
    assert payload["parameters"]["job_name"] == "example"
    schema = workflow_schema("design.kejia_idr")
    assert next(p for p in schema["parameters"] if p["name"] == "target_sequence")["required"] is True
    assert "gpu_device" not in {p["name"] for p in schema["parameters"]}
    af2 = graph["nodes"][2]["data"]["params"]
    assert af2["model_count"] == 5
    payload = node_request("af2", af2, {"run_dir": "/example", "candidates_jsonl": "/example/candidates.jsonl"}, "af2-example")
    assert payload["parameters"]["use_initial_guess"] is True
    assert new_graph("kejia_idr_generated")["nodes"][0]["data"]["kind"] == "backbone_source"


def test_cluster_detects_cpu_image_and_template_requirements():
    from mn_molecular_design_studio.cluster.capabilities import requirements
    result = requirements({"workflow": "design.kejia_idr", "parameters": {"template_library": "walle"}}, installed_images={})
    assert inputs.IMAGE in result["images"]
    assert "kejia-idr/source/templates/walle" in result["references"]
    assert result["gpu_count"] == 0


def test_form_example_submission_and_pipeline_handoff(monkeypatch, tmp_path, template):
    import importlib
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import kejia_idr_ui
    example = tmp_path / "conf_2.pdb"
    stage_template(template, example, "X", "Y")
    monkeypatch.setattr(inputs, "example_template", lambda: example)
    monkeypatch.setattr(kejia_idr_ui, "user_input_root", lambda _owner: tmp_path / "inputs")
    viewer = importlib.import_module("mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer")
    viewer_calls = []
    monkeypatch.setattr(viewer, "molstar_custom_component", lambda *args, **kwargs: viewer_calls.append((args, kwargs)))
    app = AppTest.from_string('''
import streamlit as st
from mn_molecular_design_studio.kejia_idr_ui import render_setup
def submit(schema, workflow, parameters, project, settings, label):
    st.session_state["test-submission"] = parameters
render_setup({"workflow":"design.kejia_idr"}, owner_id="test", project_selector=lambda key: "project", submit=submit)
''').run()
    assert not app.exception
    assert not viewer_calls
    next(button for button in app.button if button.label == "Load Tamarind example inputs").click().run()
    assert not app.exception
    assert app.text_area[0].value == inputs.EXAMPLE_SEQUENCE
    visualization = viewer_calls[-1][0][0][0]
    assert visualization.pdb == example.read_text()
    assert visualization.representation_type is None
    assert [(chain.chain_id, chain.representation_type) for chain in visualization.chains] == [
        ("A", "cartoon"), ("B", "cartoon+ball-and-stick"),
    ]
    assert visualization.chains[0].color_params != visualization.chains[1].color_params
    assert viewer_calls[-1][1]["show_controls"] is True
    assert any("Template peptide: chain B · 4 residues" in caption.value for caption in app.caption)
    next(button for button in app.button if button.label == "Run template threading").click().run()
    assert app.session_state["test-submission"]["target_sequence"] == inputs.EXAMPLE_SEQUENCE
    next(button for button in app.button if button.label == "Open full pipeline").click().run()
    assert not app.exception
    pending = app.session_state["pipeline-pending-draft"]
    assert pending["project_id"] == "project"
    validate_graph(pending["graph"], require_inputs=True)


def test_template_preview_tracks_chain_roles_for_mmcif(monkeypatch, tmp_path, template):
    import importlib
    from Bio.PDB import MMCIFIO
    from streamlit.testing.v1 import AppTest

    path = tmp_path / "template.cif"
    writer = MMCIFIO()
    writer.set_structure(read_structure(template))
    writer.save(str(path))
    viewer = importlib.import_module("mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer")
    calls = []
    monkeypatch.setattr(viewer, "molstar_custom_component", lambda *args, **kwargs: calls.append((args, kwargs)))
    app = AppTest.from_string('''
import streamlit as st
from pathlib import Path
from mn_molecular_design_studio.kejia_idr_ui import _render_template_preview
_render_template_preview(Path(st.session_state["template"]), "Y", "X")
''')
    app.session_state["template"] = str(path)
    app.run()
    assert not app.exception
    visualization = calls[-1][0][0][0]
    assert "loop_" in visualization.pdb
    assert [chain.chain_id for chain in visualization.chains] == ["Y", "X"]
    assert visualization.chains[1].label == "Template peptide (target position)"
