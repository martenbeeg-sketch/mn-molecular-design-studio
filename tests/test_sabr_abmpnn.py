from pathlib import Path
import json

import pytest

from mn_molecular_design_studio import mpnn_input, sabr_runtime
from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import (
    validate_workflow_request,
    workflow_schema,
)
from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import node_request
from mn_molecular_design_studio.pipeline_graph import new_graph, validate_graph


def test_sabr_detects_cdrs_and_maps_fixed_framework_residues():
    structure = mpnn_input.ASSETS / "9F6O.pdb"
    annotation = sabr_runtime.annotate_file(structure, ["B"])
    cdrs = annotation["chains"]["B"]["cdr_residues"]

    assert all(cdrs[name] for name in sabr_runtime.CDR_NAMES)
    fixed, selected = sabr_runtime.fixed_residues_from_annotation(
        structure, ["B"], ["CDR3"], cdrs
    )
    assert selected["CDR3"] == cdrs["CDR3"]
    assert not set(selected["CDR3"]) & set(fixed)
    assert set(cdrs["CDR1"] + cdrs["CDR2"]) <= set(fixed)


def test_abmpnn_and_sabr_workflows_are_exposed():
    abmpnn = workflow_schema("sequence.abmpnn")
    sabr = workflow_schema("analysis.sabr")

    assert abmpnn["workflow"] == "sequence.abmpnn"
    assert {"cdr_regions", "num_seq_per_target", "structure_path"} <= {
        field["name"] for field in abmpnn["parameters"]
    }
    assert sabr["workflow"] == "analysis.sabr"
    assert {"structure_path", "scheme", "chain_ids"} <= {
        field["name"] for field in sabr["parameters"]
    }


def test_standalone_sabr_request_validates_against_vendored_workflow(tmp_path: Path):
    structure = mpnn_input.ASSETS / "9F6O.pdb"
    request = tmp_path / "sabr-request.json"
    request.write_text(json.dumps({
        "schema_version": 1,
        "workflow": "analysis.sabr",
        "parameters": {"structure_path": str(structure), "chain_ids": ["B"]},
    }))

    result = validate_workflow_request(request)

    assert result["valid"] is True
    assert result["workflow"] == "analysis.sabr"


def test_abmpnn_pipeline_template_routes_designs_through_sabr_and_analysis():
    plan = validate_graph(new_graph("abmpnn_cdr_redesign"))
    kinds = [node["data"]["kind"] for node in plan["graph"]["nodes"]]
    edges = {(edge["source"], edge["target"]) for edge in plan["graph"]["edges"]}

    assert kinds == ["project_design_source", "sabr", "abmpnn", "boltz2", "sequence_qc", "results"]
    assert ("project_design_source-0", "sabr-1") in edges
    assert ("sabr-1", "abmpnn-2") in edges
    assert ("abmpnn-2", "boltz2-3") in edges
    assert ("boltz2-3", "sequence_qc-4") in edges
    assert ("sequence_qc-4", "results-5") in edges


def test_pipeline_project_source_preserves_hydrated_mimosa_fixed_identities():
    from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import _selected_source_rows

    unhydrated = {"candidate_id": "rf-design-1", "raw_metadata": {"generation_tag": "samples_design_0"}}
    hydrated = {
        "candidate_id": "rf-design-1",
        "raw_metadata": {
            "generation_tag": "samples_design_0",
            "fixed_motif_required": True,
            "motif_mapping_verified": True,
            "mimosa_fixed_residue_identities": {"H100": "E", "H101": "T", "H102": "G"},
        },
    }

    selected = _selected_source_rows("project_design_source", [hydrated], [unhydrated])

    assert selected == [hydrated]
    assert selected[0]["raw_metadata"]["mimosa_fixed_residue_identities"] == {
        "H100": "E", "H101": "T", "H102": "G",
    }


def test_project_design_source_accepts_completed_app_generator_jobs(tmp_path: Path):
    from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import source_project_designs

    run = tmp_path / "design-run"
    normalized = run / "artifacts/normalized_candidates"
    normalized.mkdir(parents=True)
    (run / "molecular_design_studio.json").write_text(json.dumps({"project_id": "test-project"}))
    (run / "metadata.json").write_text(json.dumps({"task_group": "design", "status": "completed"}))
    (run / "result.json").write_text(json.dumps({"success": True}))
    structure = mpnn_input.ASSETS / "9F6O.pdb"
    (normalized / "candidates.jsonl").write_text(json.dumps({
        "candidate_id": "candidate-1", "stage": "generation.backbone_sequence",
        "complex_pdb": str(structure), "binder_chains": ["B"],
    }) + "\n")

    source, rows = source_project_designs(run, "test-project", ["candidate-1"])

    assert source == run.resolve()
    assert [row["candidate_id"] for row in rows] == ["candidate-1"]
    assert Path(rows[0]["complex_pdb"]).is_file()


def test_abmpnn_hydrates_and_applies_rfantibody_fixed_mimosa_identities(tmp_path: Path):
    from mn_molecular_design_studio.mimosa_motif import (
        apply_fixed_motif_identities,
        hydrate_candidate_fixed_motif,
    )
    from mn_molecular_design_studio.engine.workflow_engine.workflows import sequence_design
    from mn_molecular_design_studio.design_validation import read_structure
    from Bio.SeqUtils import seq1

    run = tmp_path / "source-run"
    stages = run / "artifacts/raw/rfantibody/stages"
    stages.mkdir(parents=True)
    structure = run / "artifacts/raw/rfantibody/structures/generation/design.pdb"
    structure.parent.mkdir(parents=True)
    # The sampler's author numbering (H114) differs from the QV candidate
    # writer's chain-local numbering (H3). Pose index 2 identifies the actual
    # fixed residue in the emitted QV order.
    pdb_lines = []
    qv_ca = []
    serial = 1
    for chain_id, offset in (("H", 0.0), ("X", 40.0)):
        for residue_index in range(1, 5):
            x = offset + residue_index * 3.0
            for atom, dx, element in (("N", -1.0, "N"), ("CA", 0.0, "C"), ("C", 1.0, "C")):
                line = (
                    f"ATOM  {serial:5d} {atom:<4s} GLY {chain_id}{residue_index:4d}    "
                    f"{x + dx:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{30.0:6.2f}          {element:>2s}\n"
                )
                pdb_lines.append(line)
                if atom == "CA":
                    qv_ca.append(line)
                serial += 1
        pdb_lines.append("TER\n")
    structure.write_text("".join(pdb_lines) + "END\n")
    tag = "design"
    (stages / "1_rfdiffusion.qv").write_text(
        f"QV_TAG {tag}\n" + "".join(qv_ca) + "QV_TAG next\n"
    )
    (stages / f"1_rfdiffusion_{tag}.mimosa.json").write_text(json.dumps({
        "condition": "sequence_fixed",
        "fixed_residues": [{
            "zero_based_pose_index": 2,
            "generated_pose_index": 2,
            "amino_acid": "E",
            "generated_position": {"chain_id": "H", "author_residue_id": "114", "insertion_code": ""},
        }],
    }))
    candidate = hydrate_candidate_fixed_motif({
        "candidate_id": "design-1",
        "source_tool": "rfantibody",
        "complex_pdb": str(structure.relative_to(run)),
        "binder_chains": ["H"],
        "target_chains": ["X"],
        "raw_metadata": {"generation_tag": tag},
    }, run)

    raw = candidate["raw_metadata"]
    assert raw["motif_mapping_verified"] is True
    assert raw["mimosa_fixed_residue_identities"] == {"H3": "E"}
    assert raw["mimosa_fixed_residue_input_mismatches"] == [
        {"residue": "H3", "expected": "E", "observed": "G"}
    ]

    staged = tmp_path / "staged.pdb"
    staged.write_text(structure.read_text())
    patches = apply_fixed_motif_identities(candidate, staged, {})
    assert patches == {"H3": "E"}
    assert any(line.startswith("ATOM") and line[21] == "H" and line[22:26].strip() == "3" and line[17:20] == "GLU"
               for line in staged.read_text().splitlines())
    model = next(read_structure(staged).get_models())
    binder_sequence = "".join(seq1(residue.resname, undef_code="X") for residue in model["H"] if "CA" in residue)
    result = sequence_design._verify_mimosa_fixed_identity_sequence(
        candidate, binder_sequence, ["H"], staged, {"source_to_engine_residues": {"H3": "H3"}}
    )
    assert result == {"status": "verified", "count": 1, "identities": {"H3": "E"}}
    with pytest.raises(ValueError, match="changed fixed MIMOSA residue"):
        sequence_design._verify_mimosa_fixed_identity_sequence(
            candidate, "GGGG", ["H"], staged, {"source_to_engine_residues": {"H3": "H3"}}
        )


def test_abmpnn_pipeline_adapter_uses_selected_cdrs_and_dedicated_workflow(tmp_path: Path):
    source_run = tmp_path / "source"
    source_run.mkdir()
    params = {
        "sequences": 2,
        "redesign_cdr1": False,
        "redesign_cdr2": True,
        "redesign_cdr3": True,
        "sampling_temp": 0.15,
        "omit_aas": "C",
        "seed": 27,
    }

    request = node_request(
        "abmpnn", params,
        {"run_dir": str(source_run), "candidates_jsonl": str(source_run / "candidates.jsonl")},
        "Test AbMPNN",
    )

    assert request["workflow"] == "sequence.abmpnn"
    assert request["parameters"]["cdr_regions"] == ["CDR2", "CDR3"]
    assert request["parameters"]["num_seq_per_target"] == 2
    assert request["parameters"]["sampling_temp"] == 0.15
    assert "model_type" not in request["parameters"]
