"""Contract, chain-safety and shared-form checks; no inference is launched."""
import copy
import inspect
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from mn_molecular_design_studio import mpnn_input as inputs
from mn_molecular_design_studio.engine.workflow_engine.core.gpu import docker_gpu_args
from mn_molecular_design_studio.engine.workflow_engine.workflows import mpnn, sequence_design as native


def test_model_checkpoints_are_noise_specific_and_custom_models_use_protein_architecture(tmp_path):
    assert inputs.checkpoint("protein_mpnn", .02, root=tmp_path).name == "proteinmpnn_v_48_002.pt"
    assert inputs.checkpoint("ligand_mpnn", .05, root=tmp_path).name == "ligandmpnn_v_32_005_25.pt"
    assert inputs.checkpoint("hyper_mpnn", .02, root=tmp_path).name == "v48_002_epoch240_hyper.pt"
    assert inputs.checkpoint("ab_mpnn", root=tmp_path).name == "abmpnn.pt"
    assert inputs.MODELS["hyper_mpnn"]["native"] == inputs.MODELS["ab_mpnn"]["native"] == "protein_mpnn"
    with pytest.raises(ValueError, match="no checkpoint"):
        inputs.checkpoint("ligand_mpnn", .02)


@pytest.mark.parametrize("device", ["0,1", "device=0,1"])
def test_docker_gpu_device_lists_quote_csv_selector(device):
    assert docker_gpu_args(device) == ["--gpus", '"device=0,1"']


@pytest.mark.parametrize("settings", [
    {"sampling_temp": 0}, {"sampling_temp": float("nan")}, {"num_seq_per_target": 0},
    {"seed": -1}, {"seed": True}, {"omit_aas": inputs.AA},
    {"bias_aa": "A:nan"}, {"bias_aa_per_residue": {"A1": {"Q": float("inf")}}},
    {"omit_aa_per_residue": {"A1": inputs.AA}}, {"backend": "foundry", "model_type": "hyper_mpnn"},
    {"use_side_chain_context": True},
])
def test_invalid_constraints_are_rejected_before_submission(settings):
    with pytest.raises((ValueError, TypeError)):
        inputs.prepare(settings)


def test_author_residue_ranges_insertion_codes_and_blank_omission():
    info = {"residues": {k: {"kind": "protein"} for k in ["B26", "B27", "B27A"]}}
    assert inputs.residues("B26-27 B27A", info, chains=["B"]) == ["B26", "B27", "B27A"]
    with pytest.raises(ValueError, match="outside"):
        inputs.residues("A26", info, chains=["B"])
    with pytest.raises(ValueError, match="absent"):
        inputs.residues("B26-29", info)
    p = inputs.prepare({"omit_aas": "empty", "bias_aa": "W:3,P:-1", "bias_aa_per_residue": '{"A1":{"G":2}}'})
    assert p["omit_aas"] == "" and p["bias_aa_per_residue"]["A1"]["G"] == 2
    assert native._normalize_omit_aas("empty") == native._normalize_omit_aas("") == ""


def test_auto_interface_lock_requires_a_valid_sequence_matching_the_structure(design_source_factory):
    source, candidate = design_source_factory()
    assert native._candidate_sequence_matches_structure(candidate, source / "designed.cif", ["B"])
    candidate["binder_sequence"] = "KLMN"
    assert not native._candidate_sequence_matches_structure(candidate, source / "designed.cif", ["B"])
    candidate["binder_sequence"] = "KLMN?"
    assert not native._candidate_sequence_matches_structure(candidate, source / "designed.cif", ["B"])


def test_auto_interface_lock_is_applied_to_shared_mpnn(capture_runner, monkeypatch):
    source, _, collection, _, commands, _ = capture_runner
    monkeypatch.setattr(native, "_candidate_sequence_matches_structure", lambda *a: True)
    monkeypatch.setattr(native, "_interface_contacts_for_structure", lambda *a, **k: ["Z1", "Z2"])
    mpnn.run_mpnn_sequence_design(source, collection, model_type="hyper_mpnn", fix_interface_residues=True)
    fixed = commands[0][commands[0].index("--fixed_residues") + 1].split()
    assert {"Z1", "Z2"} <= set(fixed)


def test_auto_interface_lock_is_applied_to_foundry_mpnn(capture_runner, monkeypatch):
    source, candidate, collection, _, commands, _ = capture_runner
    candidate["complex_pdb"] = str(source / "designed.cif")
    collection.write_text(json.dumps(candidate) + "\n")
    monkeypatch.setattr(native, "_candidate_sequence_matches_structure", lambda *a: True)
    monkeypatch.setattr(native, "_interface_contacts_for_structure", lambda *a, **k: ["B1"])
    mpnn.run_mpnn_sequence_design(source, collection, fix_interface_residues=True)
    designed = json.loads(commands[0][commands[0].index("--designed_residues") + 1])
    assert "B1" not in designed and set(designed) == {"B2", "B3", "B4", "B5"}


@pytest.fixture
def capture_runner(monkeypatch, tmp_path, design_source_factory):
    source, candidate = design_source_factory()
    original = source / "designed.cif"
    pdb = source / "designed.pdb"
    mpnn._pdb_copy(original, pdb)
    candidate.update(stage="generation.backbone", complex_pdb=str(pdb))
    collection = source / "artifacts/normalized_candidates/candidates.jsonl"
    collection.write_text(json.dumps(candidate) + "\n")
    run = tmp_path / "mpnn-run"
    run.mkdir()
    monkeypatch.setattr(mpnn, "create_job", lambda *a, **k: SimpleNamespace(run_dir=run))
    commands = []
    monkeypatch.setattr(native, "_run_shell_steps", lambda directory, steps: commands.extend(s["command"] for s in steps) or 1)
    monkeypatch.setattr(native, "finish_job", lambda *a, **k: None)
    weights = tmp_path / "weights.pt"
    weights.write_bytes(b"checkpoint-placeholder")
    monkeypatch.setattr(inputs, "checkpoint", lambda *a, **k: weights)
    return source, candidate, collection, run, commands, weights


@pytest.mark.parametrize("model", list(inputs.MODELS))
def test_all_models_select_their_checkpoint_and_keep_target_fixed(capture_runner, model):
    source, candidate, collection, run, commands, weights = capture_runner
    result = mpnn.run_mpnn_sequence_design(source, collection, model_type=model, designed_residues="B2-3",
        fixed_residues="B1", bias_aa_per_residue={"B2": {"G": 3}}, omit_aa_per_residue={"B3": "C"},
        omit_aas="empty", seed=123)
    command = commands[0]
    from mn_molecular_design_studio.engine.workflow_engine.workflows.chain_roles import BINDER_CHAIN_ORDER
    binder = BINDER_CHAIN_ORDER[0]
    assert result == run
    assert command[command.index("--model_type") + 1] == inputs.MODELS[model]["native"]
    assert command[command.index("--chains_to_design") + 1] == binder
    assert str(weights) + ":/mpnn-checkpoint.pt:ro" in command
    assert command[command.index("--omit_AA") + 1] == ""
    fixed = command[command.index("--fixed_residues") + 1].split()
    assert set(fixed) == {binder + suffix for suffix in ("1", "4", "5")}
    bias = json.loads((run / "artifacts/raw/ligandmpnn/inputs" / candidate["candidate_id"] / "bias_aa_per_residue.json").read_text())
    assert bias == {binder + "2": {"G": 3}}
    assert json.loads(collection.read_text())["binder_chains"] == ["B"]


def test_native_foundry_receives_the_same_constraints(capture_runner):
    source, candidate, collection, run, commands, _ = capture_runner
    candidate["complex_pdb"] = str(source / "designed.cif")
    collection.write_text(json.dumps(candidate) + "\n")
    mpnn.run_mpnn_sequence_design(source, collection, designed_residues="B2-3", sampling_temp=.25,
                                 bias_aa="W:2", omit_aa_per_residue={"B2": "C"}, seed=99)
    command = commands[0]
    assert "mpnn" in command and "--designed_chains" not in command
    assert json.loads(command[command.index("--designed_residues") + 1]) == ["B2", "B3"]
    assert command[command.index("--temperature") + 1] == "0.25"
    assert json.loads(command[command.index("--bias") + 1]) == {"TRP": 2}
    assert json.loads(command[command.index("--omit_per_residue") + 1]) == {"B2": ["CYS"]}
    assert command[command.index("--seed") + 1] == "99"


def test_foundry_output_roles_follow_explicit_chain_selection(capture_runner, monkeypatch):
    import shutil
    source, candidate, collection, run, _, _ = capture_runner
    candidate["complex_pdb"] = str(source / "designed.cif")
    collection.write_text(json.dumps(candidate) + "\n")
    def fake_output(directory, steps):
        output = directory / "artifacts/raw/foundry_mpnn/output" / candidate["candidate_id"] / "designed.cif"
        output.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source / "designed.cif", output)
        return 0
    monkeypatch.setattr(native, "_run_shell_steps", fake_output)
    mpnn.run_mpnn_sequence_design(source, collection, design_chains="A")
    row = json.loads((run / "artifacts/normalized_candidates/candidates.jsonl").read_text())
    assert row["binder_chains"] == ["A"] and row["target_chains"] == ["B"]
    assert row["binder_sequence"] == "ACDEFG"
    assert row["target_pdb"].endswith("inputs/passing-design.cif")


def test_foundry_output_maps_remapped_cif_chain_ids_and_keeps_dna_as_context(tmp_path, monkeypatch):
    input_cif = tmp_path / "input.cif"
    input_cif.write_text("input")
    run_dir = tmp_path / "run"
    output_cif = run_dir / "artifacts/raw/foundry_mpnn/output/input/designed.cif"
    output_cif.parent.mkdir(parents=True)
    output_cif.write_text("output")
    input_rows = [
        {"auth_asym_id": "A", "label_asym_id": "A"},
        {"auth_asym_id": "B", "label_asym_id": "C"},
        {"auth_asym_id": "C", "label_asym_id": "F"},
    ]
    output_rows = [
        {"auth_asym_id": "A", "label_asym_id": "A", "label_atom_id": "P", "auth_seq_id": "1", "pdbx_PDB_ins_code": "?", "label_comp_id": "DA"},
        {"auth_asym_id": "C", "label_asym_id": "C", "label_atom_id": "P", "auth_seq_id": "1", "pdbx_PDB_ins_code": "?", "label_comp_id": "DT"},
        {"auth_asym_id": "F", "label_asym_id": "F", "label_atom_id": "CA", "auth_seq_id": "1", "pdbx_PDB_ins_code": "?", "label_comp_id": "ALA"},
        {"auth_asym_id": "F", "label_asym_id": "F", "label_atom_id": "CA", "auth_seq_id": "2", "pdbx_PDB_ins_code": "?", "label_comp_id": "GLY"},
    ]
    monkeypatch.setattr(native, "_cif_atom_rows", lambda path: input_rows if Path(path) == input_cif else output_rows)
    source = [{"candidate_id": "input", "binder_chains": ["C"], "target_chains": ["A", "B"], "source_tool": "structure_input"}]
    rows = native._normalize_foundry_mpnn_candidates(run_dir, source, {"input": input_cif}, {})
    assert len(rows) == 1
    assert rows[0]["binder_chains"] == ["F"]
    assert rows[0]["target_chains"] == ["A", "C"]
    assert rows[0]["binder_sequence"] == "AG"


def test_structure_input_is_staged_and_not_rewritten(capture_runner):
    source, _, _, run, commands, _ = capture_runner
    path = source / "designed.pdb"
    before = path.read_bytes()
    mpnn.run_mpnn_sequence_design(structure_path=path, design_chains="B", model_type="ab_mpnn")
    assert path.read_bytes() == before
    assert commands[0][commands[0].index("--model_type") + 1] == "protein_mpnn"
    assert (run / "artifacts/mpnn_source/candidates.jsonl").is_file()


def test_symmetry_never_ties_a_fixed_target(capture_runner):
    source, _, _, _, _, _ = capture_runner
    with pytest.raises(ValueError, match="equally sized"):
        inputs.prepare({"structure_path": str(source / "designed.pdb"), "design_chains": "A,B", "homo_oligomer": 2})


def test_batch_constraints_are_not_overwritten_by_the_next_candidate(capture_runner):
    source, candidate, collection, run, commands, _ = capture_runner
    second = dict(candidate, candidate_id="second-design")
    collection.write_text(json.dumps(candidate) + "\n" + json.dumps(second) + "\n")
    mpnn.run_mpnn_sequence_design(source, collection, bias_aa_per_residue={"B2": {"W": 2}})
    paths = [c[c.index("--bias_AA_per_residue") + 1] for c in commands]
    assert len(set(paths)) == 2
    assert all((run / p.removeprefix("/work/")).is_file() for p in paths)


def test_cif_long_chain_ids_and_author_residue_constraints_are_remapped(capture_runner):
    from Bio.PDB import MMCIFIO
    from mn_molecular_design_studio.design_validation import read_structure
    source, candidate, collection, run, commands, _ = capture_runner
    structure = read_structure(source / "designed.cif")
    structure[0]["B"].id = "AA"
    path = source / "long.cif"
    writer = MMCIFIO(); writer.set_structure(structure); writer.save(str(path))
    candidate.update(complex_pdb=str(path), binder_chains=["AA"])
    collection.write_text(json.dumps(candidate) + "\n")
    mpnn.run_mpnn_sequence_design(source, collection, backend="shared_ligandmpnn", model_type="hyper_mpnn",
                                  design_chains="AA", designed_residues="AA2-3", bias_aa_per_residue={"AA2": {"W": 2}})
    command = commands[0]
    assert set(command[command.index("--fixed_residues") + 1].split()) == {"Z1", "Z4", "Z5"}


def test_ca_only_unknown_generator_backbone_is_validated_after_staging(capture_runner):
    source, candidate, collection, run, commands, _ = capture_runner
    pdb_path = Path(candidate["complex_pdb"])
    lines = []
    for line in pdb_path.read_text().splitlines():
        if line.startswith("ATOM  ") and line[21] == "B":
            line = line[:17] + "UNK" + line[20:]
        lines.append(line)
    pdb_path.write_text("\n".join(lines) + "\n")

    mpnn.run_mpnn_sequence_design(
        source,
        collection,
        model_type="ligand_mpnn",
        fixed_residues="B2",
        fix_interface_residues=True,
    )

    command = commands[0]
    assert command[command.index("--chains_to_design") + 1] == "Z"
    assert command[command.index("--fixed_residues") + 1] == "Z2"
    staged = run / "artifacts/raw/ligandmpnn/inputs" / f"{candidate['candidate_id']}.pdb"
    assert staged.is_file()
    assert inputs.inventory(staged)["chains"]["Z"]["kind"] == "protein"


def test_ligandmpnn_preserves_nucleotide_context(tmp_path):
    path = inputs.ASSETS / "1BC8.pdb"
    staged = tmp_path / "staged.pdb"
    native._stage_ligandmpnn_role_normalized_input(source_pdb=path, staged_pdb=staged,
        candidate={"binder_chains": ["C"], "target_chains": ["A", "B"]}, requested_design_chains="C")
    original = [line[17:20].strip() for line in path.read_text().splitlines() if line.startswith("ATOM") and line[21] in "AB"]
    output = [line[17:20].strip() for line in staged.read_text().splitlines() if line.startswith("ATOM") and line[21] in "AB"]
    assert output and set(original) == set(output)
    assert "GLY" not in set(output)


def test_every_mpnn_template_uses_the_common_request_and_form():
    from mn_molecular_design_studio import pipeline_graph as graph
    from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import node_request
    for template in ["foundry", "complexa", "pxdesign", "rfdiffusion", "genie3", "protpardelle"]:
        g = graph.new_graph(template)
        node = next(n for n in g["nodes"] if n["data"]["kind"] in {"ligandmpnn", "foundry_mpnn"})
        p = copy.deepcopy(node["data"]["params"])
        p.update(model_type="hyper_mpnn", noise_level="0.1", bias_aa_per_residue='{"B5":{"W":2}}')
        request = node_request(node["data"]["kind"], p, {"run_dir": "/source", "candidates_jsonl": "/source/candidates.jsonl"}, "test")
        assert request["workflow"] == "sequence.mpnn"
        assert request["parameters"]["noise_level"] == .1
        assert request["parameters"]["bias_aa_per_residue"] == {"B5": {"W": 2}}
        inspect.signature(mpnn.run_mpnn_sequence_design).bind(**request["parameters"])


def test_common_tile_and_pipeline_forms_render_all_models_without_jobs():
    from streamlit.testing.v1 import AppTest
    app = AppTest.from_string('''
from mn_molecular_design_studio.mpnn_ui import render_setup
from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import workflow_schema
render_setup(workflow_schema("sequence.mpnn"), owner_id="test", project_selector=lambda _:None,
    search_jobs=lambda **k:{"items":[]}, submit=lambda *a:(_ for _ in ()).throw(AssertionError("No job expected")))
''').run()
    assert not app.exception
    for model in inputs.MODELS:
        app.button_group(key="mpnn-field-model").set_value(model).run()
        assert not app.exception
        assert app.button(key="mpnn-submit").disabled
    node = AppTest.from_string('''
from mn_molecular_design_studio.pipeline_graph import new_graph
from mn_molecular_design_studio.pipeline_node_ui import render_node_form
g=new_graph("foundry")
render_node_form(g["nodes"][2],g,key="common",owner_id="test",render_input=lambda:{})
''').run()
    assert not node.exception
    assert node.button_group(key="common-mpnn-field-model")


def test_api_model_contract_validates_before_any_job(tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import validate_workflow_request, AutomationError
    path = tmp_path / "request.json"
    path.write_text(json.dumps({"schema_version": 1, "workflow": "sequence.mpnn", "parameters": {"model_type": "ab_mpnn"}}))
    with pytest.raises(AutomationError, match="Upload a structure"):
        validate_workflow_request(path)


def test_source_collection_cannot_be_taken_from_a_different_run(tmp_path):
    with pytest.raises(ValueError, match="belong"):
        inputs.prepare({"source_run_dir": str(tmp_path / "one"), "candidates_jsonl": str(tmp_path / "two/candidates.jsonl")})


def test_verification_keeps_multiple_designed_chains_separate(design_source_factory):
    run, candidate = design_source_factory()
    candidate["binder_chains"] = ["A", "B"]
    (run / "artifacts/normalized_candidates/candidates.jsonl").write_text(json.dumps(candidate) + "\n")
    request = mpnn.verification_request(run)
    assert request["workflow"] == "refolding.optimized_prediction"
    assert request["parameters"]["engine"] == "esmfold"
    assert [t["protein_sequences"] for t in request["parameters"]["tasks"]] == ["ACDEFG", "KLMNP"]


def test_api_source_is_scoped_to_selected_project(monkeypatch, tmp_path):
    from fastapi import HTTPException
    from mn_molecular_design_studio import api, job_catalog
    run = tmp_path / "design/run"
    run.mkdir(parents=True)
    monkeypatch.setattr(api, "_resolve_api_job", lambda *a: {"task_group": "design"})
    monkeypatch.setattr(job_catalog, "get_job_access_record", lambda *a: {"project_id": "other"})
    with pytest.raises(HTTPException, match="selected project"):
        api._guard_source_job({"workflow": "sequence.mpnn", "parameters": {"source_run_dir": str(run)}}, {"id": "test"}, project_id="chosen")
