import hashlib
import inspect
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import pytest

from mn_molecular_design_studio import peptide_trap_input as inputs
from mn_molecular_design_studio.engine.workflow_engine.tools.peptide_trap_runtime import build_protocol_xml
from mn_molecular_design_studio.engine.workflow_engine.workflows import peptide_trap as runner
from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import node_request
from mn_molecular_design_studio.pipeline_graph import new_graph, validate_graph


NATIVE_XML = '''<ROSETTASCRIPTS><SCOREFXNS/><RESIDUE_SELECTORS/><TASKOPERATIONS>
<OperateOnResidueSubset name="repack" selector="chainA"><RestrictToRepackingRLT/></OperateOnResidueSubset>
</TASKOPERATIONS><MOVERS><AtomTree name="docking_foldtree"/>
<AddCompositionConstraintMover name="5trp"/><AddCompositionConstraintMover name="alacomp"/>
<FastDesign name="FastDesign_ppi" repeats="2"/><FastRelax name="FastRelax" repeats="2"/>
</MOVERS><FILTERS><SSPrediction cmd="/software/psipred"/></FILTERS><PROTOCOLS/></ROSETTASCRIPTS>'''


def _structure(path, peptide_length):
    lines, serial = [], 0
    for chain, count in (("A", peptide_length), ("B", 4)):
        for residue in range(1, count + 1):
            for atom, element in (("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")):
                serial += 1
                lines.append(f"ATOM  {serial:5d} {atom:^4s} GLY {chain}{residue:4d}    {float(serial):8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{20.:6.2f}          {element:>2s}")
        lines.append("TER")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\nEND\n")
    return path


@pytest.fixture
def library(monkeypatch, tmp_path):
    source = tmp_path / "source"
    for key, info in inputs.SCAFFOLDS.items():
        _structure(source / "pdb_models" / info["file"], len(info["pattern"]))
    (source / "design_after_seqthread.xml").write_text(NATIVE_XML)
    manifest = {"archive_sha256": inputs.ARCHIVE_SHA256, "files": {
        path.relative_to(source).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in source.rglob("*") if path.is_file()
    }}
    (tmp_path / "source_manifest.json").write_text(json.dumps(manifest))
    monkeypatch.setattr(inputs, "library_root", lambda: tmp_path)
    return tmp_path


def test_author_masks_registers_and_proline_exclusion():
    assert inputs.match_registers(inputs.EXAMPLE_SEQUENCE, ["C37"]) == [{
        "scaffold": "C37", "start": 1, "end": 8, "sequence": "DVRFNFRE", "burial_pattern": "10101011"}]
    assert inputs.match_registers("DDDDDDDD", ["C37"]) == []
    assert inputs.match_registers("DVRPNFRE", ["C37"]) == []
    assert inputs.match_registers("PVRFNFRP", ["C37"])
    assert inputs.match_registers("GQRIRVRITG", ["C104.1"])
    rows = inputs.match_registers("A" * 12, ["C37", "C104.1"])
    assert len(rows) == 8
    assert rows[-1]["start"] == 3 and rows[-1]["end"] == 12


@pytest.mark.parametrize("value", ["", "AAAAAAA", "A" * 5001, "AAAAAAAX", ">a\nAAAAAAAA\n>b\nAAAAAAAA"])
def test_invalid_sequences(value):
    with pytest.raises(ValueError):
        inputs.normalize_sequence(value)


def test_requests_limit_registers_and_verify_assets(library):
    p = inputs.prepare_request({"target_sequence": ">target\n" + "A" * 12, "max_registers": 2})
    assert len(p["all_matches"]) == 8 and len(p["registers"]) == 2
    for overrides in ({"scaffold_ids": []}, {"scaffold_ids": ["C37", "C37"]},
                      {"scaffold_ids": ["not-a-scaffold"]}, {"max_backbone_rmsd": float("nan")},
                      {"design_repeats": True}, {"designs_per_register": 0}):
        with pytest.raises(ValueError):
            inputs.prepare_request({"target_sequence": inputs.EXAMPLE_SEQUENCE} | overrides)
    inputs.template_path("C37").write_text("changed")
    with pytest.raises(ValueError, match="modified"):
        inputs.validate_library(["C37"])


def test_protocol_retains_native_movers_and_fixed_target_identity():
    xml = ET.fromstring(build_protocol_xml(NATIVE_XML, 1, 3))
    assert xml.find("FILTERS") is None
    assert xml.find(".//RestrictToRepackingRLT") is not None
    assert xml.find(".//FastDesign").get("repeats") == "1"
    assert xml.find(".//FastRelax").get("repeats") == "3"
    assert [step.get("mover_name") for step in xml.find("PROTOCOLS")] == [
        "docking_foldtree", "5trp", "alacomp", "FastDesign_ppi", "FastRelax"]


def test_output_sequence_integrity_and_geometry_gate(library, tmp_path):
    output = tmp_path / "results"
    output.mkdir()
    path = _structure(output / "design.pdb", 8)
    register = inputs.match_registers("GGGGGGGG", ["C37"])[0]
    record = {"pdb": path.name, "target_sequence": "GGGGGGGG", "binder_sequence": "GGGG",
              "accepted": True, "metrics": {"passes_template_geometry": True}}
    report = output / "designs.json"
    report.write_text(json.dumps({"records": [record]}))
    rows = runner.normalize_output(output, register, "GGGGGGGG")
    assert rows[0]["target_chains"] == ["A"] and rows[0]["binder_chains"] == ["B"]
    assert rows[0]["stage"] == "generation.backbone_sequence"
    record["accepted"] = False
    report.write_text(json.dumps({"records": [record]}))
    assert runner.normalize_output(output, register, "GGGGGGGG") == []
    record.update(accepted=True, target_sequence="AAAAAAAA")
    report.write_text(json.dumps({"records": [record]}))
    with pytest.raises(ValueError, match="preserve"):
        runner.normalize_output(output, register, "GGGGGGGG")
    record["pdb"] = "../outside.pdb"
    report.write_text(json.dumps({"records": [record]}))
    with pytest.raises(ValueError, match="invalid structure path"):
        runner.normalize_output(output, register, "GGGGGGGG")


def test_pipeline_schema_cpu_capabilities_and_continuation(library):
    from mn_molecular_design_studio.cluster.capabilities import requirements
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import WORKFLOWS, _get_function
    graph = new_graph("peptide_trap")
    assert [node["data"]["kind"] for node in graph["nodes"]] == ["peptide_trap", "af2", "rosetta_metrics", "results"]
    validate_graph(graph)
    graph["nodes"][0]["data"]["params"]["target_sequence"] = inputs.EXAMPLE_SEQUENCE
    validate_graph(graph, require_inputs=True)
    payload = node_request("peptide_trap", graph["nodes"][0]["data"]["params"], {}, "test")
    inspect.signature(_get_function(WORKFLOWS[payload["workflow"]])).bind(**payload["parameters"])
    capability = requirements(payload, installed_images={})
    assert capability["gpu_count"] == 0 and inputs.IMAGE in capability["images"]
    assert any(path.endswith("source_manifest.json") for path in capability["references"])
    assert new_graph("peptide_trap_generated")["nodes"][0]["data"]["kind"] == "generation_source"


def test_example_form_submission_and_pipeline(monkeypatch, library):
    import importlib
    from streamlit.testing.v1 import AppTest
    viewer = importlib.import_module("mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer")
    calls = []
    monkeypatch.setattr(viewer, "molstar_custom_component", lambda *a, **kw: calls.append((a, kw)))
    app = AppTest.from_string('''
import streamlit as st
from mn_molecular_design_studio.peptide_trap_ui import render_setup
def submit(schema, workflow, parameters, project, settings, label):
    st.session_state["submitted"] = parameters
render_setup({"workflow":"design.peptide_trap"},owner_id="owner",project_selector=lambda k:"project",submit=submit)
''').run()
    assert not app.exception
    next(button for button in app.button if button.label == "Load C37 paper example").click().run()
    assert not app.exception and app.text_area[0].value == inputs.EXAMPLE_SEQUENCE
    assert calls[-1][0][0][0].chains[1].chain_id == "A"
    next(button for button in app.button if button.label == "Generate peptide traps").click().run()
    assert app.session_state["submitted"]["scaffold_ids"] == ["C37"]
    next(button for button in app.button if button.label == "Open full pipeline").click().run()
    assert not app.exception
    validate_graph(app.session_state["pipeline-pending-draft"]["graph"], require_inputs=True)


def test_canvas_scaffold_selector_preserves_list_type():
    from streamlit.testing.v1 import AppTest
    app = AppTest.from_string('''
import streamlit as st
from mn_molecular_design_studio.pipeline_node_ui import _render_field
st.session_state["selected"] = _render_field({"name":"scaffold_ids", "type":"json", "label":"Scaffolds", "default":["C37"]}, {"scaffold_ids":["C37"]}, "test")
''').run()
    assert not app.exception
    assert app.session_state["selected"] == ["C37"]
    app.multiselect[0].set_value(["C104.1"]).run()
    assert app.session_state["selected"] == ["C104.1"]
