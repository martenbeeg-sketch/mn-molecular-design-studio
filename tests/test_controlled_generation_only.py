from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import workflow_schema
from mn_molecular_design_studio.engine.workflow_engine.workflows import design, rfantibody


def _write_structure(path: Path, chains: tuple[str, ...]) -> None:
    rows = []
    serial = 1
    for chain_index, chain in enumerate(chains):
        for residue_id, residue_name in enumerate(("GLY", "GLU", "THR", "GLY"), 1):
            x = 3.0 * residue_id
            for atom, dx, element in (("N", -1.0, "N"), ("CA", 0.0, "C"), ("C", 1.0, "C")):
                rows.append(
                    f"ATOM  {serial:5d} {atom:<4s} {residue_name} {chain}{residue_id:4d}    "
                    f"{x + dx + chain_index * 20.0:8.3f}{0.0:8.3f}{0.0:8.3f}"
                    f"{1.0:6.2f}{30.0:6.2f}          {element:>2s}\n"
                )
                serial += 1
        rows.append("TER\n")
    rows.append("END\n")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(rows))


def test_generator_only_parameters_are_public_in_workflow_schemas() -> None:
    bg = workflow_schema("design.boltzgen_mimosa")
    rfa = workflow_schema("design.rfantibody_mimosa")
    bg_params = {item["name"] for item in bg["parameters"]}
    rfa_params = {item["name"] for item in rfa["parameters"]}
    assert {"run_vanilla_pipeline", "nanobody_cdr_insertions"} <= bg_params
    assert "generation_only" in rfa_params


def test_boltzgen_fixed_nanobody_insertions_materialize_in_job_yaml(tmp_path, monkeypatch) -> None:
    target = tmp_path / "target.pdb"
    source = tmp_path / "mimosa-source.pdb"
    _write_structure(target, ("A",))
    _write_structure(source, ("P", "X"))
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    captured = {}

    monkeypatch.setattr(design, "create_job", lambda *args: SimpleNamespace(run_dir=run_dir))
    monkeypatch.setattr(design, "update_status", lambda *args, **kwargs: None)
    def capture_steps(path, steps):
        captured["steps"] = steps
        return 0

    monkeypatch.setattr(design, "_run_shell_steps", capture_steps)
    monkeypatch.setattr(design, "_finish_design_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(design.subprocess, "run", lambda *args, **kwargs: SimpleNamespace(returncode=0))

    design.run_boltzgen(
        target_pdb=target,
        target_chains=["A"],
        hotspots="A1",
        num_designs=1,
        budget=1,
        sampling_steps=20,
        run_vanilla_pipeline=False,
        task="nanobody",
        scaffold_name="7eow",
        nanobody_cdr_insertions={"CDR1": 3, "CDR2": 3, "CDR3": 7},
        mimosa_condition="sequence_fixed",
        mimosa_complex=source,
        mimosa_binder_chain="P",
        mimosa_complex_target_chain="X",
        mimosa_target_chain="A",
        mimosa_cdr="CDR3",
        mimosa_design_loop="CDR3",
        mimosa_source_positions=["2"],
    )

    params = json.loads((run_dir / "input.json").read_text())["params"]
    assert params["generator_only_tile"] is True
    assert params["mimosa_condition"] == "sequence_fixed"
    assert params["resolved_nanobody_cdr_insertions"] == {"CDR1": 3, "CDR2": 3, "CDR3": 7}
    scaffold = run_dir / "artifacts/raw/boltzgen/input/7eow.yaml"
    insertion_lines = [line.strip() for line in scaffold.read_text().splitlines() if "num_residues:" in line]
    assert insertion_lines == ["num_residues: 3", "num_residues: 3", "num_residues: 7"]
    command = captured["steps"][0]["command"][-1]
    assert "--steps design" in command
    assert "run-generation-only" in command
    assert "--steps inverse_folding" not in command
    source_scaffold = Path(design.__file__).parent.parent / "assets/boltzgen_nanobody_scaffolds/7eow.yaml"
    assert "num_residues: 1..5" in source_scaffold.read_text()


def test_rfantibody_mimosa_generation_only_skips_mpnn_and_rf2(tmp_path, monkeypatch) -> None:
    target = tmp_path / "target.pdb"
    source = tmp_path / "mimosa-source.pdb"
    _write_structure(target, ("A",))
    _write_structure(source, ("P", "X"))
    ref = tmp_path / "reference/rfantibody"
    ref.mkdir(parents=True)
    (ref / "RFdiffusion_Ab.pt").write_bytes(b"fixture")
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    commands = []
    candidate = {
        "candidate_id": "rfdiffusion_design",
        "source_tool": "rfantibody",
        "stage": "generation.backbone",
        "complex_pdb": str(target),
        "binder_chains": ["H"],
        "target_chains": ["A"],
        "metrics": {"result_kind": "generated"},
        "raw_metadata": {"generation_tag": "design", "mimosa_diagnostics": {"fixed_identity_count": 4}},
    }

    monkeypatch.setattr(rfantibody, "reference_root", lambda: tmp_path / "reference")
    monkeypatch.setattr(rfantibody, "create_job", lambda *args: SimpleNamespace(run_dir=run_dir, run_id="test-run-id"))
    monkeypatch.setattr(rfantibody, "update_status", lambda *args, **kwargs: None)
    monkeypatch.setattr(rfantibody, "finish_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(rfantibody, "_collect_generation_candidates", lambda *args, **kwargs: [candidate])
    monkeypatch.setattr(rfantibody, "_collect_candidates", lambda *args, **kwargs: pytest.fail("generation-only must not collect RF2 outputs"))
    monkeypatch.setattr(
        rfantibody.subprocess,
        "run",
        lambda command, *args, **kwargs: commands.append(command) or SimpleNamespace(returncode=0),
    )

    result = rfantibody.run_rfantibody_design(
        task="nanobody",
        target_pdb=target,
        target_chains=["A"],
        framework="h-NbBCII10",
        regions=["H1", "H2", "H3"],
        loop_lengths={"H1": "7", "H2": "6", "H3": "10"},
        number_of_designs=1,
        generation_only=True,
        mimosa_enabled=True,
        mimosa_condition="sequence_fixed",
        mimosa_complex=source,
        mimosa_binder_chain="P",
        mimosa_complex_target_chain="X",
        mimosa_target_chain="A",
        mimosa_cdr="CDR3",
        mimosa_design_loop="H3",
        mimosa_source_positions=["2"],
    )

    request = json.loads((run_dir / "input.json").read_text())
    assert request["params"]["generation_only"] is True
    assert request["params"]["use_abmpnn"] is False
    script = json.loads((run_dir / "command.json").read_text())["command"][-1]
    assert "--mimosa-packet" in script and "--mimosa-condition sequence_fixed" in script
    assert "proteinmpnn_interface_design.py" not in script
    assert "rf2 --input-quiver" not in script
    saved = json.loads((result / "artifacts/normalized_candidates/candidates.jsonl").read_text().splitlines()[0])
    assert saved["raw_metadata"]["generation_tag"] == "design"
    assert saved["raw_metadata"]["mimosa_diagnostics"]["fixed_identity_count"] == 4
    assert len(commands) == 2


def test_boltzgen_generation_only_motif_identities_reach_common_abmpnn(tmp_path) -> None:
    from Bio.SeqUtils import seq1

    from mn_molecular_design_studio.design_validation import read_structure
    from mn_molecular_design_studio.engine.workflow_engine.workflows.sequence_design import (
        _verify_mimosa_fixed_identity_sequence,
    )
    from mn_molecular_design_studio.mimosa_motif import (
        apply_fixed_motif_identities,
        hydrate_candidate_fixed_motif,
    )

    run_dir = tmp_path / "run"
    structure = run_dir / "artifacts/raw/boltzgen/run-generation-only/intermediate_designs/design.pdb"
    metadata = structure.with_suffix(".npz")
    _write_structure(structure, ("H", "X"))
    np.savez(
        metadata,
        mimosa_condition=np.asarray("sequence_fixed"),
        # BoltzGen amino-acid token IDs are canonical-order + 2; token 3 is R.
        mimosa_fixed_aa=np.asarray([-1, 3, -1, -1, -1, -1, -1, -1]),
        source_context=np.asarray(json.dumps({"chains": [
            {"source_chain": "B", "output_chain": "H", "indices": [0, 1, 2, 3]},
            {"source_chain": "X", "output_chain": "X", "indices": [0, 1, 2, 3]},
        ]})),
    )
    candidate = hydrate_candidate_fixed_motif({
        "candidate_id": "boltzgen_00001",
        "source_tool": "boltzgen_generator",
        "complex_pdb": str(structure.relative_to(run_dir)),
        "binder_chains": ["H"],
        "target_chains": ["X"],
        "raw_metadata": {"metadata_path": str(metadata.relative_to(run_dir))},
    }, run_dir)
    assert candidate["raw_metadata"]["motif_mapping_verified"] is True
    assert candidate["raw_metadata"]["mimosa_fixed_residue_identities"] == {"H2": "R"}

    staged = tmp_path / "abmpnn-input.pdb"
    staged.write_text(structure.read_text())
    assert apply_fixed_motif_identities(candidate, staged, {}) == {"H2": "R"}
    model = next(read_structure(staged).get_models())
    sequence = "".join(seq1(residue.resname, undef_code="X") for residue in model["H"] if "CA" in residue)
    assert sequence == "GRTG"
    result = _verify_mimosa_fixed_identity_sequence(candidate, sequence, ["H"], staged, {})
    assert result == {"status": "verified", "count": 1, "identities": {"H2": "R"}}
    with pytest.raises(ValueError, match="changed fixed MIMOSA residue"):
        _verify_mimosa_fixed_identity_sequence(candidate, "GETG", ["H"], staged, {})
