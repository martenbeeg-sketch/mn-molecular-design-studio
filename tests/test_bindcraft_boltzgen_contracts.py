from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from mn_molecular_design_studio.engine.workflow_engine.core import jobs
from mn_molecular_design_studio.engine.workflow_engine.core import portable_paths, workflow_queue
from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import workflow_schema
from mn_molecular_design_studio.engine.workflow_engine.workflows import bindcraft2, design
from mn_molecular_design_studio.job_catalog import _job_payload


def test_bindcraft2_settings_keep_job_name_modality_hotspots_and_lengths():
    settings = bindcraft2.build_bindcraft2_settings(
        target_chains=["A"],
        binder_length="65-90",
        hotspots="A45,A52",
        campaign_name="bindcraft2-contract",
        modality="VHH",
        design_properties=["forced_targeting"],
        campaign_seed=0,
        number_of_final_designs=3,
        max_trajectories=12,
    )

    assert settings["campaign_name"] == "bindcraft2-contract"
    assert settings["modality"] == "VHH"
    assert settings["targets"][0]["chains"] == "A"
    assert settings["targets"][0]["hotspots"] == "A45,A52"
    assert settings["binder_lengths"] == [65, 90]
    assert settings["forced_targeting"] is True
    assert settings["campaign_seed"] == 0
    assert settings["number_of_final_designs"] == 3
    assert settings["max_trajectories"] == 12


def test_bindcraft2_bounded_attempt_replaces_local_targets_and_records_its_seed():
    settings = bindcraft2.build_bindcraft2_settings(target_chains=["A"], binder_length="60")
    attempt = bindcraft2.bounded_bindcraft2_settings(settings, 947)
    assert attempt["number_of_final_designs"] == 1
    assert attempt["max_trajectories"] == 1
    assert attempt["campaign_seed"] == 947
    assert attempt["save_failed_trajectories"] is True
    assert settings["number_of_final_designs"] == 1
    with pytest.raises(ValueError, match="nonnegative"):
        bindcraft2.bounded_bindcraft2_settings(settings, -1)


def test_bindcraft2_rejects_invalid_hotspot_and_incompatible_modality():
    with pytest.raises(ValueError, match="selected target chains"):
        bindcraft2.build_bindcraft2_settings(
            target_chains=["A"], binder_length="80", hotspots="B45",
        )

    with pytest.raises(ValueError, match="fold_switch objective"):
        bindcraft2.build_bindcraft2_settings(
            target_chains=["A"], binder_length="80", modality=["VHH", "fold_switch"],
        )


def test_bindcraft2_normalizes_a_single_disordered_target_fasta():
    assert bindcraft2.parse_disordered_target_fasta(">idR target\nacde\nFGH\n") == (
        "idR target", "ACDEFGH",
    )
    with pytest.raises(ValueError, match="single >header"):
        bindcraft2.parse_disordered_target_fasta(">one\nACD\n>two\nEFG\n")


def test_bindcraft2_campaign_state_is_the_authoritative_attempt_count(tmp_path):
    run_dir = tmp_path / "run"
    output = run_dir / "artifacts/raw/bindcraft2/output"
    output.mkdir(parents=True)
    (output / ".campaign_state.json").write_text(
        json.dumps({"trajectories": 3, "accepted": 1, "rejections": {"terminated": {"screen": 2}}})
    )
    assert bindcraft2.read_bindcraft2_campaign_state(run_dir) == {
        "trajectories": 3,
        "accepted": 1,
        "rejections": {"terminated": {"screen": 2}},
    }
    with pytest.raises(RuntimeError, match="invalid trajectories"):
        (output / ".campaign_state.json").write_text(
            json.dumps({"trajectories": 0.5, "accepted": 0, "rejections": {}})
        )
        bindcraft2.read_bindcraft2_campaign_state(run_dir)


def test_bindcraft2_distributed_tile_creates_a_cpu_controller_parent(monkeypatch, tmp_path):
    target = tmp_path / "target.pdb"
    target.write_text(
        "ATOM      1  CA  ALA A   1      10.000  11.000  12.000  1.00 20.00           C\nEND\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(jobs, "runs_root", lambda: tmp_path / "runs")
    monkeypatch.setattr(bindcraft2, "missing_bindcraft2_parameters", lambda *_args, **_kwargs: [])
    monkeypatch.setattr(bindcraft2, "bindcraft2_image_available", lambda *_args, **_kwargs: True)
    dispatched = []
    monkeypatch.setattr(bindcraft2, "spawn_worker_for_run", lambda path: dispatched.append(path))

    run_dir = bindcraft2.enqueue_bindcraft2_design(
        target_pdb=target,
        target_chains=["A"],
        binder_length="60",
        distributed_campaign={
            "enabled": True,
            "accepted_target": 2,
            "max_attempts": 4,
            "batch_size": 2,
            "max_concurrent_batches": 2,
        },
        gpu_device="auto",
    )
    params = jobs.read_json(run_dir / "input.json")["params"]
    assert jobs.read_json(run_dir / "metadata.json")["job_type"] == "bindcraft2_campaign"
    assert params["gpu_device"] == "none"
    assert params["distributed_campaign"]["gpu_selection"] == "auto"
    assert jobs.read_json(run_dir / "worker_request.json")["kind"] == "bindcraft2_campaign_controller"
    assert dispatched == [run_dir]


def test_bindcraft_stages_job_name_and_settings_without_running_docker(monkeypatch, tmp_path: Path):
    runs_root = tmp_path / "runs"
    target_pdb = tmp_path / "target.pdb"
    target_pdb.write_text(
        "ATOM      1  CA  ALA A   1      10.000  11.000  12.000  1.00 20.00           C\nEND\n",
        encoding="utf-8",
    )
    captured: dict[str, object] = {}

    def fake_container_runner(run_dir: Path, steps: list[dict]) -> int:
        captured["steps"] = steps
        return 1  # Verify staged files without starting the BindCraft container.

    monkeypatch.setattr(jobs, "runs_root", lambda: runs_root)
    monkeypatch.setattr(design, "load_manifest", lambda _name: {"image": "bindcraft-test:unit"})
    monkeypatch.setattr(design, "_run_shell_steps", fake_container_runner)

    run_dir = design.run_bindcraft(
        target_pdb=target_pdb,
        target_chains=["A"],
        binder_length="60-80",
        hotspots="A1",
        campaign_name="bindcraft-contract",
        number_of_final_designs=2,
        max_trajectories=7,
        num_seqs_override=8,
        max_mpnn_sequences_override=10,
        enable_mpnn=False,
    )

    raw_dir = run_dir / "artifacts" / "raw" / "bindcraft"
    bindcraft_input = json.loads((raw_dir / "input.json").read_text(encoding="utf-8"))
    advanced = json.loads((raw_dir / "settings_advanced.json").read_text(encoding="utf-8"))
    assert bindcraft_input["binder_name"] == "bindcraft-contract"
    assert bindcraft_input["chains"] == "A"
    assert bindcraft_input["target_hotspot_residues"] == "A1"
    assert bindcraft_input["lengths"] == [60, 80]
    assert advanced["max_trajectories"] == 7
    assert advanced["num_seqs"] == 8
    assert advanced["max_mpnn_sequences"] == 10
    assert advanced["enable_mpnn"] is False
    assert _job_payload(run_dir)["campaign_name"] == "bindcraft-contract"
    assert "bindcraft-test:unit" in captured["steps"][0]["command"]


def test_boltzgen_schema_exposes_task_ligand_and_pipeline_controls():
    schema = workflow_schema("design.boltzgen")
    parameters = {row["name"]: row for row in schema["parameters"]}

    for name in ("task", "ligand_type", "ligand_value", "run_vanilla_pipeline", "budget"):
        assert name in parameters
    assert parameters["task"]["type"] == "str"
    assert parameters["task"]["default"] == "protein"


def test_boltzgen_stages_small_molecule_smiles_before_container_execution(monkeypatch, tmp_path: Path):
    runs_root = tmp_path / "runs"
    captured: dict[str, object] = {}

    def fake_container_runner(run_dir: Path, steps: list[dict]) -> int:
        captured["run_dir"] = run_dir
        captured["steps"] = steps
        return 1  # Stop after staging; this test never starts Docker.

    monkeypatch.setattr(jobs, "runs_root", lambda: runs_root)
    monkeypatch.setattr(design, "load_manifest", lambda _name: {"image": "boltzgen-test:unit"})
    monkeypatch.setattr(design, "_run_shell_steps", fake_container_runner)

    run_dir = design.run_boltzgen(
        task="small_molecule",
        ligand_type="smiles",
        ligand_value="NCC(=O)O",
        binder_length="60-80",
        campaign_name="boltzgen-small-molecule-test",
        num_designs=3,
        budget=2,
        run_vanilla_pipeline=True,
    )

    input_yaml = run_dir / "artifacts" / "raw" / "boltzgen" / "input" / "target_binder.yaml"
    payload = yaml.safe_load(input_yaml.read_text(encoding="utf-8"))
    assert payload["entities"][0]["protein"]["sequence"] == "60..80"
    assert payload["entities"][1]["ligand"] == {"id": "B", "smiles": "NCC(=O)O"}

    input_params = json.loads((run_dir / "input.json").read_text(encoding="utf-8"))["params"]
    assert input_params["campaign_name"] == "boltzgen-small-molecule-test"
    assert input_params["protocol"] == "protein-small_molecule"
    assert input_params["budget"] == 2
    assert jobs.read_json(run_dir / "input.json")["inputs"]["target_pdb"] == ""
    command = captured["steps"][0]["command"]
    assert "--steps design inverse_folding folding design_folding affinity analysis filtering" in command[-1]


def test_boltzgen_nanobody_stages_optional_target_hotspots_as_binding_guidance(monkeypatch, tmp_path: Path):
    target_pdb = tmp_path / "target.pdb"
    target_pdb.write_text(
        "ATOM      1  CA  ALA X  41      10.000  11.000  12.000  1.00 20.00           C\n"
        "ATOM      2  CA  GLY X  42      13.000  14.000  15.000  1.00 20.00           C\n"
        "END\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(jobs, "runs_root", lambda: tmp_path / "runs")
    monkeypatch.setattr(design, "load_manifest", lambda _name: {"image": "boltzgen-test:unit"})
    monkeypatch.setattr(design, "_run_shell_steps", lambda *_args: 1)

    run_dir = design.run_boltzgen(
        task="nanobody",
        target_pdb=target_pdb,
        target_chains=["X"],
        hotspots="X42",
        scaffold_name="7eow",
        num_designs=1,
        budget=1,
        run_vanilla_pipeline=False,
    )

    input_yaml = run_dir / "artifacts" / "raw" / "boltzgen" / "input" / "target_binder.yaml"
    payload = yaml.safe_load(input_yaml.read_text(encoding="utf-8"))
    target_file = payload["entities"][0]["file"]
    assert target_file["binding_types"] == [{"chain": {"id": "X", "binding": 2}}]


def test_boltzgen_requires_hotspots_when_site_adherence_filter_is_enabled():
    with pytest.raises(ValueError, match="Select at least one binding-site residue"):
        design.run_boltzgen(
            task="small_molecule",
            ligand_type="ccd",
            ligand_value="ATP",
            filter_binding_site_adherence=True,
        )


def test_boltzgen_stages_queued_inline_yaml_without_path_prefix(monkeypatch, tmp_path):
    runs_root = tmp_path / "runs"
    request_dir = runs_root / "design" / "request"
    target = tmp_path / "target.pdb"
    target.write_text(
        "ATOM      1  CA  ALA B   1      10.000  11.000  12.000  1.00 20.00           C\nEND\n"
    )
    template = (
        "entities:\n  - protein:\n      id: {{DESIGN_CHAIN}}\n      sequence: 60..60\n"
        "  - file:\n      path: {{TARGET_FILE}}\n      include:\n"
        "        - chain:\n            id: {{TARGET_CHAIN}}\n"
    )
    monkeypatch.setattr(jobs, "runs_root", lambda: runs_root)
    monkeypatch.setattr(portable_paths, "runs_root", lambda: runs_root)
    monkeypatch.setattr(design, "load_manifest", lambda _: {"image": "boltzgen-test:unit"})
    monkeypatch.setattr(design, "_run_shell_steps", lambda *args: 1)
    request_path = request_dir / "worker_request.json"
    kwargs = dict(target_pdb=target, target_chains=["B"], task="yaml", design_yaml=template)
    jobs.write_json(request_path, {"kwargs": workflow_queue._invocation_value(kwargs)})
    restored = workflow_queue._restore_invocation_value(jobs.read_json(request_path)["kwargs"])

    run_dir = design.run_boltzgen(**restored)
    staged_yaml = run_dir / "artifacts/raw/boltzgen/input/target_binder.yaml"
    payload = yaml.safe_load(staged_yaml.read_text())
    assert list(payload) == ["entities"]
    assert payload["entities"][0]["protein"] == {"id": "A", "sequence": "60..60"}
    assert payload["entities"][1]["file"]["path"] == "target.pdb"
    assert jobs.read_json(run_dir / "input.json")["params"]["design_yaml"] == template


@pytest.mark.parametrize("target_crop", ["DEFG", "ACGH"])
@pytest.mark.parametrize("inline_sequence", [True, False])
def test_bindcraft2_sequence_result_keeps_cropped_target_and_binder_roles(monkeypatch, tmp_path, target_crop, inline_sequence):
    ranked = tmp_path / "artifacts/raw/bindcraft2/output/3_Ranked"
    ranked.mkdir(parents=True)
    (ranked / "!_Ranked.csv").write_text("design,Binder_Sequence\naccepted,LLLL\n")
    (ranked / "accepted.cif").touch()
    (tmp_path / "artifacts/raw/bindcraft2/target.fasta").write_text(">target\nACDEFGH\n")
    monkeypatch.setattr(bindcraft2, "_sequences_by_chain", lambda _: {"A": target_crop, "B": "LLLL"})
    monkeypatch.setattr(bindcraft2, "write_candidates", lambda _run, _tool, candidates: candidates)

    settings = {"target_type": "sequence", "target_chains": []}
    if inline_sequence:
        settings["target_sequence"] = "ACDEFGH"
    candidates = bindcraft2.normalize_bindcraft2_candidates(tmp_path, settings)
    assert candidates[0]["target_chains"] == ["A"]
    assert candidates[0]["binder_chains"] == ["B"]
    assert candidates[0]["binder_sequence"] == "LLLL"
    assert candidates[0]["metrics"]["target_chain_inference"] == "fasta_crop_or_binder_sequence_match"


def test_bindcraft2_ranked_candidates_use_native_refolded_pass_status(monkeypatch, tmp_path):
    ranked = tmp_path / "artifacts/raw/bindcraft2/output/3_Ranked"
    refolded = tmp_path / "artifacts/raw/bindcraft2/output/2_Refolded"
    ranked.mkdir(parents=True)
    refolded.mkdir(parents=True)
    (ranked / "!_Ranked.csv").write_text(
        "rank,design,Binder_Sequence,hash\n1,passed,LLLL,shared-hash\n2,rejected,MMMM,shared-hash\n"
    )
    (refolded / "!_Refolded.csv").write_text(
        "design,Binder_Sequence,hash,outcome,failed_filters\npassed,LLLL,shared-hash,passed,\nrejected,MMMM,shared-hash,rejected,i_pTM\n"
    )
    (ranked / "passed.cif").touch()
    (ranked / "rejected.cif").touch()
    (tmp_path / "artifacts/raw/bindcraft2/target.fasta").write_text(">target\nACDEFGH\n")
    monkeypatch.setattr(
        bindcraft2,
        "_sequences_by_chain",
        lambda path: {"A": "ACDEFGH", "B": "LLLL" if path.stem == "passed" else "MMMM"},
    )
    monkeypatch.setattr(bindcraft2, "write_candidates", lambda _run, _tool, candidates: candidates)

    candidates = bindcraft2.normalize_bindcraft2_candidates(
        tmp_path, {"target_type": "sequence", "target_sequence": "ACDEFGH"}
    )

    assert [candidate["metrics"]["bindcraft2_rank"] for candidate in candidates] == [1, 2]
    assert [candidate["metrics"]["native_pass_filters"] for candidate in candidates] == [True, False]
    assert candidates[1]["raw_metadata"]["bindcraft2_refolded_outcome"] == "rejected"
    assert candidates[1]["raw_metadata"]["bindcraft2_failed_filters"] == "i_pTM"
