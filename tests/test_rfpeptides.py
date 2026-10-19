"""RFpeptides form, workflow and target-preparation contracts."""
import json
from pathlib import Path

import pytest

from mn_molecular_design_studio import rfpeptides_input as inputs
from mn_molecular_design_studio.design_validation import protein_chains, read_structure
from mn_molecular_design_studio.engine.workflow_engine.core.manifests import load_manifest
from mn_molecular_design_studio.engine.workflow_engine.services import workflow_automation as automation

EXAMPLE_TARGET = Path(__file__).resolve().parents[1] / "mn_molecular_design_studio" / "resources" / "rfpeptides" / "7zkr_GABARAP.pdb"


def request_parameters(**changes):
    return {
        "target_pdb": str(EXAMPLE_TARGET),
        "target_chains": ["A"],
        "binder_length": "12-18",
        "hotspots": "A51,A52,A50",
        "number_of_designs": 1,
        "diffusion_steps": 50,
    } | changes


def test_rfpeptides_form_defaults_match_tamarind_fields():
    normalized = inputs.validate_request(request_parameters())
    assert normalized["binder_length"] == "12-18"
    assert normalized["number_of_designs"] == 1
    assert normalized["diffusion_steps"] == 50
    assert normalized["hotspots"] == ["A51", "A52", "A50"]


def test_contig_format_supports_one_or_multiple_target_chains():
    assert inputs.build_contig("12-18", ["B"], {"B": 120}) == "12-18 B1-120/0"
    assert inputs.build_contig("16", ["B", "C"], {"B": 120, "C": 80}) == "16 B1-120/0 C1-80/0"


def test_rfpeptides_input_page_renders_without_submitting():
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_string('''
from mn_molecular_design_studio.rfpeptides_ui import render_setup
from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import workflow_schema
render_setup(workflow_schema("design.rfpeptides"), owner_id="pytest-rfpeptides",
    project_selector=lambda _: None, submit=lambda *args: None)
''').run(timeout=20)
    assert not app.exception
    assert any(widget.label == "PDB File" for widget in app.button_group)
    assert {widget.label for widget in app.text_input} >= {"Job Name", "Binder Length", "Target PDB code"}


@pytest.mark.parametrize("length", ["", "12-8", "2", "201", "12-x"])
def test_invalid_peptide_lengths_are_rejected(length):
    with pytest.raises(ValueError):
        inputs.validate_request(request_parameters(binder_length=length))


@pytest.mark.parametrize("field,value", [
    ("number_of_designs", 0),
    ("number_of_designs", 100001),
    ("diffusion_steps", 0),
    ("diffusion_steps", 201),
    ("hotspots", "B51"),
    ("hotspots", "A9999"),
    ("target_chains", []),
])
def test_invalid_rfpeptides_requests_are_rejected(field, value):
    with pytest.raises(ValueError):
        inputs.validate_request(request_parameters(**{field: value}))


def test_target_preparation_remaps_chains_and_hotspots_for_cyclic_binder(tmp_path):
    output = tmp_path / "target.pdb"
    mapping = inputs.prepare_target(EXAMPLE_TARGET, ["A"], ["A51", "A52"], output)
    chains = protein_chains(read_structure(output))
    assert mapping["chain_map"] == {"A": "B"}
    assert mapping["target_chains"] == ["B"]
    assert mapping["hotspots"] == ["B49", "B50"]
    assert list(chains) == ["B"]
    assert len(chains["B"]) == mapping["chain_lengths"]["B"]


def test_workflow_registration_manifest_and_api_validation(tmp_path):
    schema = automation.workflow_schema("design.rfpeptides")
    assert schema["workflow"] == "design.rfpeptides"
    assert {parameter["name"] for parameter in schema["parameters"]} >= {
        "binder_length", "hotspots", "number_of_designs", "diffusion_steps",
    }
    manifest = load_manifest("rfpeptides")
    assert manifest["image"] == inputs.IMAGE
    assert manifest["gpu"] is True
    request = {
        "schema_version": 1,
        "workflow": "design.rfpeptides",
        "parameters": request_parameters(),
        "resources": {"cpu_cores": 4, "gpu_device": "0"},
    }
    path = tmp_path / "request.json"
    path.write_text(json.dumps(request))
    assert automation.validate_workflow_request(path)["valid"] is True
    request["parameters"]["binder_length"] = "18-12"
    path.write_text(json.dumps(request))
    with pytest.raises(automation.AutomationError):
        automation.validate_workflow_request(path)
