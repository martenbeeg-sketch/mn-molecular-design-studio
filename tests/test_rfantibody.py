"""RFantibody input and workflow contracts, without starting a GPU job."""
import json
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

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


def test_proteinmpnn_command_matches_pinned_image_cli():
    command = workflow._proteinmpnn_command(
        input_quiver="/work/input.qv",
        output_quiver="/work/output.qv",
        loops="H1,H2,H3",
        sequences_per_structure=1,
        temperature=0.1,
        augment_eps=0.05,
        num_connections=48,
        omit_aas="CX",
        checkpoint="/ref/ProteinMPNN.pt",
    )

    assert command[:2] == ["python", "scripts/proteinmpnn_interface_design.py"]
    assert command[command.index("-num_connections") + 1] == "48"
    assert "--num-connections" not in command
    assert "--bias-aas" not in command
    assert "--bias-by-res-json" not in command


def test_rf2_repeats_have_distinct_candidate_and_binder_file_ids(monkeypatch, tmp_path):
    root = tmp_path / "artifacts" / "raw" / "rfantibody"
    (root / "structures").mkdir(parents=True)
    for repeat in (1, 2):
        model = root / "rf2_models" / f"model_{repeat:02d}"
        model.mkdir(parents=True)
        (model / "samples_design_0_dldesign_0_best.pdb").write_text("fixture")
    monkeypatch.setattr(workflow, "read_structure", lambda _: object())
    monkeypatch.setattr(workflow, "protein_chains", lambda _: {"H": [object()], "A": [object()]})
    monkeypatch.setattr(workflow, "sequence", lambda _: "ACD")
    monkeypatch.setattr(workflow, "_write_binder_structure", lambda _src, dst, _chains: dst.write_text("binder"))
    monkeypatch.setattr(workflow, "_score_from_pdb", lambda _: {"pred_lddt": 0.9, "interaction_pae": 1.0})

    candidates = workflow._collect_candidates(tmp_path, tmp_path / "target.pdb", ["A"], [],
                                               calculate_epitope_distance=False)

    assert len(candidates) == 2
    assert [candidate["metrics"]["rf2_repeat"] for candidate in candidates] == [1, 2]
    assert len({candidate["candidate_id"] for candidate in candidates}) == 2
    assert len({candidate["binder_pdb"] for candidate in candidates}) == 2
    assert all(Path(candidate["binder_pdb"]).is_file() for candidate in candidates)


def test_unsupported_rfantibody_bias_settings_are_rejected_before_submission():
    with pytest.raises(ValueError, match="does not support global or per-residue"):
        inputs.validate_request(nanobody_parameters(bias_aas="W:2"))
    with pytest.raises(ValueError, match="does not support global or per-residue"):
        inputs.validate_request(nanobody_parameters(bias_by_res_json='{"H":{"31":{"W":2}}}'))


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


def test_rfantibody_results_end_with_score_interpretation(monkeypatch, tmp_path):
    from mn_molecular_design_studio import rfantibody_ui

    structure = tmp_path / "complex.pdb"
    structure.write_text("ATOM sample\n")
    monkeypatch.setattr(rfantibody_ui, "molstar_custom_component", lambda *args, **kwargs: None)
    candidate = {
        "candidate_id": "test-design",
        "binder_sequence": "ACDE",
        "complex_pdb": str(structure),
        "binder_length": 4,
        "metrics": {"pred_lddt": 0.9, "interaction_pae": 8.0},
    }
    at = AppTest.from_string(
        "from pathlib import Path\n"
        "from mn_molecular_design_studio.rfantibody_ui import render_results\n"
        f"render_results('test-run', Path({str(tmp_path)!r}), "
        f"[{candidate!r}], {{}}, {{}})"
    ).run()
    assert not at.exception
    assert any("How to interpret these scores" in block.value for block in at.markdown)
    assert any("<10 Å" in block.value and "0.90 corresponds to 90/100" in block.value for block in at.markdown)


def test_rfantibody_mimosa_results_separate_stages_and_default_to_best_rf2(monkeypatch, tmp_path):
    from mn_molecular_design_studio import rfantibody_ui

    structure = tmp_path / "complex.pdb"
    structure.write_text("ATOM sample\n")
    monkeypatch.setattr(rfantibody_ui, "molstar_custom_component", lambda *args, **kwargs: None)
    candidates = [
        {"candidate_id": "generated", "stage": "generation.backbone", "binder_chains": ["H"], "complex_pdb": str(structure), "metrics": {"generated_motif_rmsd_angstrom": 1.2}},
        {"candidate_id": "rf2-1", "stage": "complex_refolding", "complex_pdb": str(structure), "metrics": {"rf2_repeat": 1, "pred_lddt": 0.90, "interaction_pae": 5.0}},
        {"candidate_id": "rf2-2", "stage": "complex_refolding", "complex_pdb": str(structure), "metrics": {"rf2_repeat": 2, "pred_lddt": 0.92, "interaction_pae": 7.0}},
    ]
    at = AppTest.from_string(
        "from pathlib import Path\n"
        "from mn_molecular_design_studio.rfantibody_ui import render_results\n"
        f"render_results('mimosa-run', Path({str(tmp_path)!r}), "
        f"{candidates!r}, {{}}, {{}})"
    ).run()
    assert not at.exception
    assert any("1 generated backbone · 2 RF2 predictions" in block.value for block in at.caption)
    assert at.selectbox[0].value == 2
    table = at.dataframe[0].value
    assert table.loc[0, "Stage"] == "Generated backbone"
    assert table.loc[1, "Stage"] == "RF2 prediction"
    assert table.loc[0, "Rank"] == "—"
    assert table.loc[0, "Sequence"] == "—"


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
