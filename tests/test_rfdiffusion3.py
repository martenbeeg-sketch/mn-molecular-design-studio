"""RFdiffusion3 task, staging, upload and result contracts, without a GPU."""
import gzip
import io
import json
import shutil
import zipfile
from pathlib import Path
from types import SimpleNamespace

import pytest

from mn_molecular_design_studio import rfdiffusion3_input as inputs
from mn_molecular_design_studio.engine.workflow_engine.workflows import rfdiffusion3 as workflow
from mn_molecular_design_studio.engine.workflow_engine.services import workflow_automation as automation
from mn_molecular_design_studio.engine.workflow_engine.core import jobs


def task_parameters(task, mode="atom"):
    p = {"task": task, "num_designs": 1, "length_min": 100, "length_max": 100}
    names = {"protein_binder": "5o45_cropped.pdb", "enzyme": "M0255_1mg5.pdb", "nucleic_acid_binder": "1bna.pdb", "ligand_binder": "IAI.pdb"}
    if task == "json":
        return p | {"json_config": {"unconditional": {"length": "64-64", "is_non_loopy": True}}}
    p["target_structure"] = str(inputs.ASSETS / names[task])
    if task == "protein_binder":
        p.update(target_chains=["A"], hotspots="A56,A115,A123")
    elif task == "nucleic_acid_binder":
        p.update(target_chains=["A", "B"], ori_token=[24, 20, 10])
    elif task == "ligand_binder":
        p.update(ligands=["IAI"], atom_conditions=[{"condition": "Fixed", "ligand": "IAI", "atoms": ""}, {"condition": "Buried", "ligand": "IAI", "atoms": "ALL"}])
    else:
        p.update(ligands=["NAI", "ACT"], scaffolding_mode=mode, fixed_atoms={"A108": "ND2,CG", "A139": "OG,CB,CA", "A152": "OH,CZ", "A156": "NZ,CE,CD"}, unindex="all", motif_contig="20,A108,20,A139,20,A152,20,A156,20")
    return p


@pytest.mark.parametrize("task,mode", [(t,"atom") for t in inputs.TASKS] + [("enzyme","residue")])
def test_all_task_forms_compile(task, mode):
    prepared = inputs.prepare_request(task_parameters(task, mode))
    spec = prepared["specifications"]["spec_001"]
    assert any(spec.get(k) for k in ("length", "contig", "unindex"))
    if task == "enzyme" and mode == "atom":
        assert spec["unindex"] == "A108,A139,A152,A156"
        assert spec["select_fixed_atoms"]["A108"] == "ND2,CG"
    if task == "nucleic_acid_binder":
        assert "A1-12,/0,B13-24" in spec["contig"]
    if task == "ligand_binder":
        assert spec["select_fixed_atoms"] == {"IAI": ""}


@pytest.mark.parametrize("changes", [
    {"num_designs": 0}, {"diffusion_batch_size": 0}, {"length_min": 200, "length_max": 100},
    {"step_scale": float("nan")}, {"target_chains": ["Z"]}, {"hotspots": "A9999"},
    {"hotspots": "A56ZZ"}, {"ori_token": [1,2]}, {"ligands": ["MISSING"]},
])
def test_bad_form_rejected_before_job(changes):
    with pytest.raises(ValueError):
        inputs.prepare_request(task_parameters("protein_binder") | changes)


@pytest.mark.parametrize("spec", [
    {"length": "64", "cyclic_chains": ["A"]}, {"contig": "A1-10"}, {"length": "64", "ligand": "ATP"},
    {"length": "64", "input": "/etc/passwd"}, {"length": "64", "dialect": 1}, {},
])
def test_custom_json_cannot_bypass_input_validation(spec):
    with pytest.raises(ValueError):
        inputs.prepare_request({"task": "json", "json_config": {"test": spec}})


def test_atom_selection_and_duplicates_checked():
    p = task_parameters("enzyme")
    with pytest.raises(ValueError, match="has no atom"):
        inputs.prepare_request(p | {"fixed_atoms": {"A108": "DOES_NOT_EXIST"}})
    with pytest.raises(ValueError, match="at least one"):
        inputs.prepare_request(p | {"fixed_atoms": {"A108": ""}})
    p = task_parameters("ligand_binder")
    with pytest.raises(ValueError, match="Combine atoms"):
        inputs.prepare_request(p | {"atom_conditions": p["atom_conditions"] * 2})


def test_nonprotein_upload_opt_in(monkeypatch, tmp_path):
    from mn_molecular_design_studio import target_editor
    monkeypatch.setattr(target_editor, "user_input_root", lambda _: tmp_path)
    pdb = (inputs.ASSETS / "1bna.pdb").read_bytes()
    path = target_editor._target_storage_path(pdb, "dna.pdb", owner_id="test", allow_non_protein=True)
    assert inputs.inventory(path)["chains"]["A"]["kind"] == "dna"
    with pytest.raises(ValueError):
        target_editor._target_storage_path(b"invalid", "invalid.pdb", owner_id="test", allow_non_protein=True)


def test_zip_paths_are_not_extracted(monkeypatch, tmp_path):
    from mn_molecular_design_studio import rfdiffusion3_ui as ui
    monkeypatch.setattr(ui, "user_input_root", lambda _: tmp_path)
    data = io.BytesIO()
    with zipfile.ZipFile(data, "w") as z:
        z.writestr("../../dna.pdb", (inputs.ASSETS / "1bna.pdb").read_bytes())
    upload = SimpleNamespace(name="input.zip", getvalue=data.getvalue)
    files = ui.store_json_structures([upload], "test")
    assert len(files) == 1 and files[0].name == "dna.pdb"
    files[0].resolve().relative_to(tmp_path)
    with pytest.raises(ValueError, match="unique filenames"):
        ui.store_json_structures([upload, upload], "test")


def test_new_workflow_registered_and_preflight_checks_inputs():
    schema = automation.workflow_schema("design.rfdiffusion3")
    assert schema["workflow"] == "design.rfdiffusion3"
    spec = automation._workflow_spec("design.rfdiffusion3") if hasattr(automation, "_workflow_spec") else next(s for s in automation._SPECS if s.workflow_id == "design.rfdiffusion3")
    automation._validate_workflow_inputs(spec, task_parameters("json"))
    with pytest.raises(automation.AutomationError):
        automation._validate_workflow_inputs(spec, {"task": "protein_binder"})


def test_runner_stages_inputs_uses_existing_image_and_limits_final_batch(monkeypatch, tmp_path):
    run = tmp_path / "run"; run.mkdir()
    refs = tmp_path / "reference/foundry"; refs.mkdir(parents=True)
    (refs / "rfd3_latest.ckpt").touch()
    monkeypatch.setattr(workflow, "create_job", lambda *a: SimpleNamespace(run_dir=run, run_id="test-12345678"))
    monkeypatch.setattr(workflow, "update_status", lambda *a, **k: None)
    monkeypatch.setattr(workflow, "reference_root", lambda: refs.parent)
    finished = []
    monkeypatch.setattr(workflow, "finish_job", lambda *a: finished.append(a))
    def execute(command, **kwargs):
        assert inputs.IMAGE in command and command[command.index("--entrypoint") + 1] == "rfd3"
        assert "diffusion_batch_size=2" in command and "n_batches=2" in command
        assert "low_memory_mode=False" in command
        prepared = json.loads((run / "artifacts/raw/rfdiffusion3/inputs.json").read_text())
        assert prepared["spec_001"]["input"].startswith("/work/")
        assert (run / "artifacts/raw/rfdiffusion3/inputs/5o45_cropped.pdb").is_file()
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(workflow.subprocess, "run", execute)
    requests = []
    monkeypatch.setattr(workflow, "_normalize", lambda _, prepared, count: requests.append(count) or [])
    workflow.run_rfdiffusion3_design(**(task_parameters("protein_binder") | {"num_designs": 3, "diffusion_batch_size": 2}))
    assert requests == [3]
    assert finished and finished[0][1] is False  # No normalized outputs: cannot claim success.


@pytest.mark.parametrize("task", inputs.TASKS)
def test_each_loaded_task_form_renders(monkeypatch, tmp_path, task):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import rfdiffusion3_ui, target_editor
    monkeypatch.setattr(rfdiffusion3_ui, "user_input_root", lambda _: tmp_path)
    monkeypatch.setattr(target_editor, "user_input_root", lambda _: tmp_path)
    app = AppTest.from_string('''
from mn_molecular_design_studio.rfdiffusion3_ui import render_setup
render_setup({}, owner_id="test", header=lambda *a: None,
             project_selector=lambda *a: "test-project", submit=lambda *a: None)
''')
    app.session_state["rfd3-task-pills"] = task
    if task != "json":
        source = tmp_path / rfdiffusion3_ui.EXAMPLES[task]
        shutil.copy2(inputs.ASSETS / source.name, source)
        app.session_state[f"rfd3-{task}-fetched-target"] = str(source)
    app.run(timeout=15)
    assert not app.exception, [str(e.value) for e in app.exception]
    assert not app.error, [e.value for e in app.error]
    assert any(button.label == "Submit design" and not button.disabled for button in app.button)


def test_json_references_use_staged_filenames_but_uploaded_paths_remain_guarded(monkeypatch, tmp_path):
    from mn_molecular_design_studio import api
    from fastapi import HTTPException
    monkeypatch.setattr(api, "user_input_root", lambda _: tmp_path / "owned")
    monkeypatch.setattr(api, "reference_root", lambda: tmp_path / "reference")
    monkeypatch.setattr(api, "path_is_within_accessible_job", lambda *a, **k: False)
    monkeypatch.setattr(api, "accessible_project_ids", lambda _: [])
    user = {"id": "test", "role": "user"}
    request = {"json_config": {"dna": {"input": "1bna.pdb", "contig": "A1-12,/0,100"}}}
    api._validate_user_path_values(request, user, tmp_path)
    with pytest.raises(HTTPException) as error:
        api._validate_user_path_values(request | {"input_files": [str(inputs.ASSETS / "1bna.pdb")]}, user, tmp_path)
    assert error.value.status_code == 403


def test_native_normalizer_reads_gzip_metrics_and_keeps_motif_map(design_source_factory, tmp_path):
    source, _ = design_source_factory()
    run = tmp_path / "rfd3"
    output = run / "artifacts/raw/rfdiffusion3/output"; output.mkdir(parents=True)
    native = output / "studio_spec_001_0_model_0.cif.gz"
    native.write_bytes(gzip.compress((source / "designed.cif").read_bytes()))
    native.with_name(native.name.removesuffix(".cif.gz") + ".json").write_text(json.dumps({"seed": 42, "metrics": {"radius_of_gyration": 12.5}, "diffused_index_map": {"A108": "B1"}}))
    prepared = inputs.prepare_request(task_parameters("enzyme"))
    candidates = workflow._normalize(run, prepared, 1)
    assert len(candidates) == 1
    assert candidates[0]["metrics"]["radius_of_gyration"] == 12.5
    assert candidates[0]["raw_metadata"]["fixed_residues"] == ["B1"]
    assert candidates[0]["stage"] == "generation.backbone"
    assert candidates[0]["binder_sequence"] is None
    assert "native_pass_filters" not in candidates[0]["metrics"]
    assert Path(candidates[0]["complex_pdb"]).is_file()


def test_foundry_rf3_summary_is_discovered_in_nested_output(tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.workflows import design as design_workflow

    output = tmp_path / "rf3" / "candidate" / "input" / "input"
    output.mkdir(parents=True)
    summary = output / "input_summary_confidences.json"
    summary.write_text(json.dumps({"ranking_score": 0.82, "iptm": 0.71, "ptm": 0.88, "has_clash": False}))

    discovered = design_workflow._find_rf3_summary(tmp_path / "rf3" / "candidate")

    assert discovered == summary
    assert design_workflow._rf3_summary_metrics(discovered) == {
        "ranking_score": 0.82,
        "ptm": 0.88,
        "iptm": 0.71,
        "has_clash": False,
    }


def test_native_normalizer_retains_dna_and_maps_protein_only_as_designed(design_source_factory, tmp_path):
    from Bio.PDB import MMCIFIO
    from mn_molecular_design_studio.design_validation import read_structure
    source, _ = design_source_factory()
    dna = read_structure(inputs.ASSETS / "1bna.pdb")
    protein = next(read_structure(source / "designed.cif").get_models())["A"].copy()
    protein.id = "C"
    next(dna.get_models()).add(protein)
    run = tmp_path / "rfd3-dna"
    output = run / "artifacts/raw/rfdiffusion3/output"; output.mkdir(parents=True)
    writer = MMCIFIO(); writer.set_structure(dna)
    writer.save(str(output / "studio_spec_001_0_model_0.cif"))
    candidates = workflow._normalize(run, inputs.prepare_request(task_parameters("nucleic_acid_binder")), 1)
    candidate = candidates[0]
    assert candidate["binder_chains"] == ["C"]
    info = inputs.inventory(Path(candidate["complex_pdb"]))
    assert info["chains"]["A"]["kind"] == info["chains"]["B"]["kind"] == "dna"
    assert candidate["metrics"]["generated_residue_count"] == 6
    assert candidate["raw_metadata"]["input_target_chains"] == ["A", "B"]


def test_generation_target_alignment_uses_original_input_chain_ids(design_source_factory):
    from Bio.PDB import Chain, Model, PDBIO, Structure
    from mn_molecular_design_studio.design_validation import read_structure
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import _with_generation_structure_metrics

    source, _ = design_source_factory()
    output = read_structure(source / "designed.cif")
    target = Structure.Structure("target")
    model = Model.Model(0); target.add(model)
    target_chain = output[0]["B"].copy(); target_chain.id = "A"; model.add(target_chain)
    target_path = source / "target.pdb"
    writer = PDBIO(); writer.set_structure(target); writer.save(str(target_path))
    candidate = {
        "candidate_id": "rf3-alignment", "source_tool": "rfdiffusion3_foundry",
        "complex_pdb": str(source / "designed.cif"), "target_pdb": str(target_path),
        "binder_chains": ["A"], "target_chains": ["B"], "metrics": {},
        "raw_metadata": {"input_target_chains": ["A"]},
    }

    result = _with_generation_structure_metrics(source, candidate)

    assert result["metrics"]["generation_target_alignment_status"] == "ok"
    assert result["metrics"]["generation_target_alignment_ca_count"] == 5
    assert result["metrics"]["generation_target_alignment_rmsd"] < 1e-6


def test_sequence_design_keeps_rf3_motif_residues_after_chain_role_remapping():
    from mn_molecular_design_studio.engine.workflow_engine.workflows.sequence_design import _fixed_residues_for_candidate
    roles = {"chain_roles": {"chain_map": {
        "binder": [{"original_chain": "B", "engine_chain": "A"}],
        "targets": [{"original_chain": "C", "engine_chain": "B"}],
    }}}
    rf3_candidate = {"source_tool": "rfdiffusion3_foundry", "raw_metadata": {"fixed_residues": ["B4", "B5"]}}
    assert _fixed_residues_for_candidate(rf3_candidate, ["C9", "B4"], roles) == ["A4", "A5", "B9"]
    regular_candidate = {"source_tool": "boltz2", "raw_metadata": {"fixed_residues": ["B4"]}}
    assert _fixed_residues_for_candidate(regular_candidate, ["C9"], roles) == ["A4", "B9"]
