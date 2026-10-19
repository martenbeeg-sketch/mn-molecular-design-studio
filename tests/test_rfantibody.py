"""RFantibody input and workflow contracts, without starting a GPU job."""
import json
from pathlib import Path

import pytest

from mn_molecular_design_studio import rfantibody_input as inputs
from mn_molecular_design_studio.engine.workflow_engine.services import workflow_automation as automation
from mn_molecular_design_studio.engine.workflow_engine.workflows import rfantibody as workflow
from mn_molecular_design_studio.design_validation import protein_chains, read_structure


EXAMPLE_TARGET = inputs.EXAMPLE_INPUTS / "flu_HA.pdb"


def nanobody_parameters(**changes):
    return {
        "task": "nanobody",
        "target_pdb": str(EXAMPLE_TARGET),
        "target_chains": ["B"],
        "hotspots": "B146,B170,B177",
        "framework": "h-NbBCII10",
        "regions": ["H1", "H2", "H3"],
        "number_of_designs": 2,
    } | changes


def test_loop_selections_and_lengths_follow_rfantibody_syntax():
    assert inputs.build_design_loops("nanobody", ["H1", "H3"], {"H1": "10", "H3": "10-15"}) == "H1:10,H3:10-15"
    assert inputs.build_design_loops("antibody", ["L1", "H3"], {"L1": "auto", "H3": "12"}) == "H3:12,L1:"


@pytest.mark.parametrize("regions", [[], ["L1"], ["X1"]])
def test_invalid_nanobody_regions_are_rejected(regions):
    with pytest.raises(ValueError):
        inputs.normalize_regions("nanobody", regions)


def test_hotspots_global_and_per_residue_bias_are_normalized():
    assert inputs.normalize_hotspots("B146, B170", ["B"]) == ["B146", "B170"]
    assert inputs.normalize_global_bias("W:3.0,P:-1") == "W:3,P:-1"
    assert json.loads(inputs.normalize_per_residue_bias('{"H":{"31":{"W":2}}}')) == {"H": {"31": {"W": 2.0}}}


@pytest.mark.parametrize("field,value", [
    ("hotspots", "A146"),
    ("hotspots", "B9999"),
    ("number_of_designs", 0),
    ("mpnn_temperature", 3.0),
    ("bias_aas", "W:nope"),
    ("bias_by_res_json", "not json"),
])
def test_invalid_rfantibody_inputs_are_rejected_before_job(field, value):
    with pytest.raises(ValueError):
        inputs.validate_request(nanobody_parameters(**{field: value}))


def test_tamarind_supported_design_setting_ranges_are_accepted():
    normalized = inputs.validate_request(nanobody_parameters(
        number_of_designs=100000,
        number_of_sequences=100,
        rf2_repeats=20,
        num_recycles=0,
        num_connections=100,
        mpnn_temperature=0.01,
    ))
    assert normalized["number_of_designs"] == 100000
    assert normalized["number_of_sequences"] == 100
    assert normalized["rf2_repeats"] == 20
    assert normalized["num_recycles"] == 0


def test_target_is_prepared_with_valid_chain_mapping_and_hotspots(tmp_path):
    prepared = workflow._prepare_target(EXAMPLE_TARGET, ["B"], ["B146", "B170"], tmp_path / "target.pdb")
    assert prepared["chain_map"] == {"B": "B"}
    assert prepared["hotspots"] == ["B146", "B170"]
    assert list(protein_chains(read_structure(tmp_path / "target.pdb"))) == ["B"]


def test_queued_workflow_job_parameters_serialize_resolved_structure_paths(monkeypatch, tmp_path):
    references = tmp_path / "references"
    antibody_weights = references / "rfantibody"
    antibody_weights.mkdir(parents=True)
    for name in ("RFdiffusion_Ab.pt", "ProteinMPNN_v48_noise_0.2.pt", "RF2_ab.pt"):
        (antibody_weights / name).touch()
    monkeypatch.setattr(workflow, "reference_root", lambda: references)
    captured = {}

    class StopAfterCapture(Exception):
        pass

    def capture_job(_task_group, _job_type, _tool, _inputs, parameters):
        captured.update(parameters)
        json.dumps(parameters)
        raise StopAfterCapture

    monkeypatch.setattr(workflow, "create_job", capture_job)
    with pytest.raises(StopAfterCapture):
        workflow.run_rfantibody_design(
            task="nanobody",
            target_pdb=EXAMPLE_TARGET,
            target_chains=["B"],
            hotspots="B146,B170,B177",
            framework="h-NbBCII10",
            regions=["H1", "H2", "H3"],
        )

    assert captured["target_pdb"] == str(EXAMPLE_TARGET)
    assert captured["framework_pdb"] is None


def test_workflow_registration_and_schema_validation(tmp_path):
    schema = automation.workflow_schema("design.rfantibody")
    assert schema["workflow"] == "design.rfantibody"
    assert any(parameter["name"] == "rf2_repeats" for parameter in schema["parameters"])
    request = {
        "schema_version": 1,
        "workflow": "design.rfantibody",
        "parameters": nanobody_parameters(number_of_designs=1),
        "resources": {"cpu_cores": 2, "gpu_device": "0"},
    }
    path = tmp_path / "request.json"
    path.write_text(json.dumps(request))
    assert automation.validate_workflow_request(path)["valid"] is True
    request["parameters"]["regions"] = ["L1"]
    path.write_text(json.dumps(request))
    with pytest.raises(automation.AutomationError):
        automation.validate_workflow_request(path)
