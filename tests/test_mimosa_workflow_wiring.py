from __future__ import annotations

import json
from pathlib import Path

import pytest

from mn_molecular_design_studio import rfantibody_input
from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import (
    WORKFLOWS,
    workflow_schema,
)
from mn_molecular_design_studio.engine.workflow_engine.workflows import rfantibody
from mn_molecular_design_studio.output_vocabulary import candidate_output


def _write_test_complex(path: Path) -> None:
    rows = []
    serial = 1
    for chain_id, offset in (("X", 0.0), ("H", 10.0), ("L", 10.0)):
        for residue_id in range(1, 4):
            x = offset + residue_id * 3.0
            for atom, dx, element in (("N", -1.0, "N"), ("CA", 0.0, "C"), ("C", 1.0, "C")):
                rows.append(
                    f"ATOM  {serial:5d} {atom:<4s} GLY {chain_id}{residue_id:4d}    "
                    f"{x + dx:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{30.0:6.2f}          {element:>2s}\n"
                )
                serial += 1
        rows.append("TER\n")
    rows.append("END\n")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(rows))


@pytest.mark.parametrize(
    "workflow_id",
    ["design.rfantibody_mimosa", "design.boltzgen_mimosa"],
)
def test_mimosa_workflows_are_registered_with_selection_fields(workflow_id: str) -> None:
    assert workflow_id in WORKFLOWS
    schema = workflow_schema(workflow_id)
    parameter_names = {item["name"] for item in schema["parameters"]}
    assert {"mimosa_condition", "mimosa_source_positions"}.issubset(parameter_names)


def test_rfantibody_mimosa_input_keeps_selected_source_positions() -> None:
    normalized = rfantibody_input.validate_request(
        {
            "task": "nanobody",
            "framework": "h-NbBCII10",
            "regions": ["H3"],
            "target_chains": ["A"],
            "mimosa_enabled": True,
            "mimosa_condition": "sequence_fixed",
            "mimosa_complex": "/tmp/source-complex.pdb",
            "mimosa_binder_chain": "B",
            "mimosa_complex_target_chain": "A",
            "mimosa_target_chain": "A",
            "mimosa_cdr": "CDR3",
            "mimosa_design_loop": "H3",
            "mimosa_source_positions": ["21", "22A"],
        },
        require_files=False,
    )

    assert normalized["mimosa_condition"] == "sequence_fixed"
    assert normalized["mimosa_source_positions"] == ["21", "22A"]


def test_rfantibody_mimosa_accepts_scattered_motif_across_all_heavy_cdrs() -> None:
    normalized = rfantibody_input.validate_request(
        {
            "task": "nanobody",
            "framework": "h-NbBCII10",
            "regions": ["H1", "H2", "H3"],
            "target_chains": ["E"],
            "mimosa_enabled": True,
            "mimosa_condition": "sequence_fixed",
            "mimosa_complex": "/tmp/source-complex.cif",
            "mimosa_binder_chain": "A,B",
            "mimosa_complex_target_chain": "E",
            "mimosa_target_chain": "E",
            "mimosa_cdr": "CDR3",
            "mimosa_design_loop": "all_cdrs",
            "mimosa_source_positions": ["A:33", "B:94"],
            "diffusion_deterministic": True,
        },
        require_files=False,
    )

    assert normalized["mimosa_design_loop"] == "ALL_CDRS"
    assert normalized["mimosa_binder_chain"] == "A,B"
    assert normalized["mimosa_source_positions"] == ["A:33", "B:94"]
    assert normalized["diffusion_deterministic"] is True


def test_rfantibody_mimosa_requires_at_least_one_source_residue() -> None:
    with pytest.raises(ValueError, match="Select at least one contacting amino acid"):
        rfantibody_input.validate_request(
            {
                "task": "nanobody",
                "framework": "h-NbBCII10",
                "regions": ["H3"],
                "target_chains": ["A"],
                "mimosa_enabled": True,
                "mimosa_condition": "geometry_only",
                "mimosa_complex": "/tmp/source-complex.pdb",
                "mimosa_binder_chain": "B",
                "mimosa_complex_target_chain": "A",
                "mimosa_target_chain": "A",
                "mimosa_cdr": "CDR3",
                "mimosa_design_loop": "H3",
                "mimosa_source_positions": [],
            },
            require_files=False,
        )


def test_rfantibody_generation_candidate_keeps_metrics_and_links_refold(tmp_path: Path) -> None:
    raw = tmp_path / "artifacts/raw/rfantibody"
    stages = raw / "stages"
    qv = stages / "1_rfdiffusion.qv"
    pdb = tmp_path / "design.pdb"
    _write_test_complex(pdb)
    pdb_text = pdb.read_text()
    qv.parent.mkdir(parents=True)
    qv.write_text("QV_TAG accepted\nQV_SCORE accepted test=1\n" + pdb_text + "\n")
    rejected = stages / "rejected_backbones/rejected.pdb"
    rejected.parent.mkdir()
    rejected.write_text(pdb_text)
    report = {
        "candidates": [
            {
                "tag": "accepted",
                "backbone_integrity_status": "passed",
                "backbone_bad_bond_count": 0,
                "backbone_checked_bond_count": 2,
                "generated_motif_rmsd_angstrom": 1.25,
                "generated_motif_potential_half_sum_angstrom_squared": 3.125,
                "fixed_residues": [{"generated_pose_index": 2}],
                "sampler_diagnostics": {"num_particles": 4, "effective_forward_batch_sizes": [4]},
            },
            {
                "tag": "rejected",
                "backbone_integrity_status": "rejected",
                "backbone_bad_bond_count": 1,
                "backbone_checked_bond_count": 2,
                "generated_motif_rmsd_angstrom": 1.75,
                "generated_motif_potential_half_sum_angstrom_squared": 6.125,
                "fixed_residues": [],
                "sampler_diagnostics": {"num_particles": 4, "effective_forward_batch_sizes": [4]},
            },
        ]
    }
    (stages / "mimosa_integrity_results.json").write_text(json.dumps(report))

    generated = rfantibody._collect_generation_candidates(
        tmp_path, pdb, ["X"], ["X1"]
    )
    by_tag = {row["raw_metadata"]["generation_tag"]: row for row in generated}
    assert set(by_tag) == {"accepted", "rejected"}
    assert by_tag["accepted"]["stage"] == "generation.backbone"
    assert by_tag["accepted"]["metrics"]["generated_motif_rmsd_angstrom"] == 1.25
    assert by_tag["accepted"]["metrics"]["fixed_motif_residue_count"] == 1
    assert by_tag["rejected"]["metrics"]["result_kind"] == "structural_rejection"
    assert candidate_output(by_tag["accepted"])["metrics"]["mimosa.generated_motif_rmsd"] == 1.25

    refolded_path = raw / "rf2_models/model_01/accepted_dldesign_0_best.pdb"
    _write_test_complex(refolded_path)
    refolds = rfantibody._collect_candidates(
        tmp_path, pdb, ["X"], ["X1"], calculate_epitope_distance=False,
        generation_candidates=generated,
    )
    assert refolds[0]["parents"] == [by_tag["accepted"]["candidate_id"]]
    assert "mimosa.generated_motif_rmsd" not in candidate_output(refolds[0])["metrics"]


def test_refolded_motif_metric_is_recomputed_in_target_aligned_frame(tmp_path: Path) -> None:
    structure = tmp_path / "refold.pdb"
    _write_test_complex(structure)
    stages = tmp_path / "artifacts/raw/rfantibody/stages"
    stages.mkdir(parents=True)
    # RFdiffusion's real sidecar convention includes the generated stage
    # prefix: 1_rfdiffusion_<QV tag>.mimosa.json.
    (stages / "1_rfdiffusion_accepted.mimosa.json").write_text(json.dumps({
        "source": {"target_sequence": "GGG"},
        "assignment_method": "sliding_window",
        "guided_motif_residues": 1,
        "target_alignment": [
            {"target_offset": i, "ca": [3.0 * (i + 1), 0.0, 0.0]} for i in range(3)
        ],
        "guided_loop_positions": [
            {"chain_id": "H", "author_residue_id": str(i), "insertion_code": ""}
            for i in range(10, 13)
        ],
        "motif_assignments": [{
            "generated_pose_index": 0,
            "source_ca": [13.0, 0.0, 0.0],
            # Sampler author numbering is retained for lineage, but output
            # PDB/QV residues are renumbered 1..N.
            "generated_position": {"chain_id": "H", "author_residue_id": "10", "insertion_code": ""},
        }],
    }))
    qv_lines = structure.read_text().splitlines(keepends=True)
    qv_lines = [line for chain_id in ("H", "X") for line in qv_lines
                if line.startswith("ATOM") and line[21] == chain_id]
    (stages / "1_rfdiffusion.qv").write_text("QV_TAG accepted\n" + "".join(qv_lines))

    assert rfantibody._refolded_motif_rmsd(tmp_path, "accepted", structure) == pytest.approx(0.0)


def test_refolded_motif_metric_maps_global_pose_index_into_light_chain(tmp_path: Path) -> None:
    structure = tmp_path / "refold.pdb"
    _write_test_complex(structure)
    stages = tmp_path / "artifacts/raw/rfantibody/stages"
    stages.mkdir(parents=True)
    (stages / "1_rfdiffusion_light.mimosa.json").write_text(json.dumps({
        "source": {"target_sequence": "GGG"},
        "assignment_method": "sliding_window",
        "guided_motif_residues": 1,
        "target_alignment": [
            {"target_offset": i, "ca": [3.0 * (i + 1), 0.0, 0.0]} for i in range(3)
        ],
        "guided_loop_positions": [
            {"chain_id": "L", "author_residue_id": str(i), "insertion_code": "",
             "zero_based_pose_index": i - 17}
            for i in range(20, 23)
        ],
        "motif_assignments": [{
            "generated_pose_index": 3,
            "source_ca": [13.0, 0.0, 0.0],
            "generated_position": {"chain_id": "L", "author_residue_id": "20", "insertion_code": ""},
        }],
    }))
    qv_lines = structure.read_text().splitlines(keepends=True)
    qv_lines = [line for chain_id in ("H", "L", "X") for line in qv_lines
                if line.startswith("ATOM") and line[21] == chain_id]
    (stages / "1_rfdiffusion.qv").write_text("QV_TAG light\n" + "".join(qv_lines))

    assert rfantibody._refolded_motif_rmsd(tmp_path, "light", structure) == pytest.approx(0.0)


def test_geometry_only_sidecar_reports_zero_fixed_identities() -> None:
    metrics = rfantibody._mimosa_generation_metrics({
        "condition": "geometry_only",
        "fixed_residues": [{"zero_based_pose_index": 4}],
        "generated_motif_rmsd_angstrom": 1.1,
    })
    assert metrics["fixed_motif_residue_count"] == 0


def test_rfantibody_baseline_ignores_dormant_mimosa_condition() -> None:
    normalized = rfantibody_input.validate_request(
        {
            "task": "nanobody",
            "framework": "h-NbBCII10",
            "regions": ["H3"],
            "target_chains": ["A"],
            "mimosa_enabled": False,
            "mimosa_condition": "sequence_fixed",
            "mimosa_source_positions": [],
        },
        require_files=False,
    )

    assert normalized["mimosa_enabled"] is False
    assert normalized["mimosa_condition"] == "unconstrained"
    assert normalized["mimosa_source_positions"] == []
