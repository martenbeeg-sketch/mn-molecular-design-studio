import json
import math

import pytest

from mn_molecular_design_studio.output_vocabulary import (
    PREFIX, candidate_output, engine_name, fasta_metrics, normalize_metrics,
)
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import normalize_candidate, write_candidates, read_candidates
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio import molecule_imports


@pytest.mark.parametrize("tool", ["pxdesign_generator", "genie3_generator", "rfdiffusion_generator", "rfdiffusion_classic", "rfdiffusion_beta_pairing", "rfdiffusion3_foundry", "protpardelle_generator", "protpardelle_1c", "boltzgen", "boltzgen_generator", "BindCraft2", "bindcraft2_generator", "proteina_complexa", "rfantibody", "rfpeptides", "esmfold2_binder_design"])
def test_generators_share_only_explicitly_equivalent_geometry(tool):
    native = {"binder_ca_radius_of_gyration": 12.5, "radius_of_gyration": 99,
              "binder_radius_of_gyration": 13, "binder_geometry_status": "ok", "binder_helix_fraction": .6,
              "helix_fraction": .9, "binder_ca_trace_helix_like_fraction": .4,
              "complex_plddt": .92, "binder_plddt": .88, "confidence_score": .7}
    output = normalize_metrics(native, tool)
    assert output["metrics"]["geometry.binder.ca_rg"] == 12.5
    assert output["metrics"]["geometry.binder.heavy_atom_rg"] == 13
    assert output["metrics"]["secondary_structure.binder.dssp.helix_fraction"] == .6
    assert output["metrics"]["secondary_structure.binder.ca_trace.helix_like_fraction"] == .4
    assert 99 not in output["metrics"].values()
    assert not any("plddt" in name or "confidence" in name for name in output["metrics"])
    assert native["helix_fraction"] == .9


def test_conflicts_invalid_values_and_model_scopes_are_not_guessed():
    result = normalize_metrics({"binder_end_to_end_distance": 11, "binder_ca_end_to_end_distance": 12,
                                "binder_helix_fraction": 65, "binder_sheet_fraction": True,
                                "binder_coil_fraction": float("nan")})
    assert result["metrics"] == {}
    assert "geometry.binder.ca_end_to_end" in result["conflicts"]
    assert len(result["issues"]) == 3
    assert normalize_metrics({"score": 1.2, "seq_rec": .4}, "unknown")["metrics"] == {}


@pytest.mark.parametrize("tool", ["protein_mpnn", "ligand_mpnn", "soluble_mpnn", "hyper_mpnn", "ab_mpnn"])
def test_mpnn_headers_recovered_and_models_kept_separate(tool):
    header = "input, id=1, T=0.1, overall_confidence=3.0e-1, ligand_confidence=.5, seq_rec=0.2500"
    candidate = {"stage": "sequence_design", "source_tool": tool, "binder_sequence": "ACD:EF",
                 "metrics": {}, "raw_metadata": {"fasta_header": header}}
    result = candidate_output(candidate)
    model = engine_name(tool)
    assert result["metrics"][f"sequence.{model}.redesigned_nll"] == pytest.approx(-math.log(.3))
    assert result["metrics"]["sequence.recovery_designed"] == .25
    assert result["metrics"]["design.sequence_length"] == 5
    assert result["provenance"][f"sequence.{model}.redesigned_nll"]["aliases"][0]["source_field"] == "fasta_header.overall_confidence"
    assert candidate["metrics"] == {}


def test_fasta_parser_uses_exact_keys_and_preserves_native_scale(tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.workflows.sequence_design import _read_fasta_sequences
    header = "global_score=2.4, score=1.2e+0, seq_recovery=.7, overall_confidence=0.3010"
    assert fasta_metrics(header)["score"] == 1.2  # not the suffix of global_score
    fasta = tmp_path / "seq.fasta"
    fasta.write_text(">" + header + "\nACDE\n")
    result = _read_fasta_sequences(fasta)[0]
    assert result["metrics"]["global_score"] == 2.4
    common = normalize_metrics(result["metrics"], "protein_mpnn")
    assert common["metrics"]["sequence.proteinmpnn.global_nll"] == 2.4
    assert common["metrics"]["sequence.proteinmpnn.redesigned_nll"] == 1.2


def test_generator_confidence_has_declared_units_scope_and_stage():
    native = {"pLDDT": .88, "Target_pLDDT": .9, "Unbound_Binder_pLDDT": .8, "i_pAE": .2, "i_pTM": .7, "pTM": .75}
    trajectory = candidate_output({"source_tool": "bindcraft2_generator", "stage": "generation.backbone_sequence", "metrics": native})
    validation = candidate_output({"source_tool": "bindcraft2", "stage": "complex_refolding", "metrics": native})
    assert trajectory["engine"] == validation["engine"] == "bindcraft2"
    assert trajectory["metrics"]["confidence.generator.target.plddt"] == 90
    assert trajectory["metrics"]["confidence.generator.bound_binder.plddt"] == 88
    assert trajectory["metrics"]["confidence.generator.unbound_binder.plddt"] == 80
    assert trajectory["metrics"]["confidence.generator.interface.mean_pae"] == pytest.approx(6.2)
    assert validation["metrics"]["confidence.validation.complex.iptm"] == .7
    assert "confidence.generator.complex.iptm" not in validation["metrics"]
    assert not any("binder.plddt" == key for key in trajectory["metrics"])
    # The same spelling from an unaudited tool is not assigned BindCraft's scope.
    assert candidate_output({"source_tool": "other", "stage": "generation.backbone_sequence", "metrics": native})["metrics"] == {}
    assert native["i_pAE"] == .2
    bad_scale = candidate_output({"source_tool": "bindcraft2_generator", "stage": "generation.backbone_sequence", "metrics": {"Target_pLDDT": 90}})
    assert bad_scale["metrics"] == {}


def test_rf_boltzgen_and_bindcraft_validation_tm_scores_share_family_not_method():
    result = [candidate_output({"source_tool": tool, "stage": "complex_refolding", "metrics": {key: .8}})
              for tool, key in [("rfdiffusion3_foundry", "iptm"), ("boltzgen", "iptm"), ("bindcraft2", "i_pTM")]]
    assert all(row["metrics"]["confidence.validation.complex.iptm"] == .8 for row in result)
    assert len({row["provenance"]["confidence.validation.complex.iptm"]["method"] for row in result}) == 3


def test_legacy_ambiguous_radius_is_not_renamed_without_shared_method_marker():
    assert normalize_metrics({"binder_radius_of_gyration": 12}, "rfdiffusion_classic")["metrics"] == {}


def test_out_of_range_numbers_cannot_break_database_projection():
    assert normalize_metrics({"binder_ca_radius_of_gyration": 10 ** 400})["metrics"] == {}
    output = candidate_output({"source_tool": "bindcraft2_generator", "stage": "generation.backbone_sequence", "metrics": {"i_pAE": 1e308}})
    assert output["metrics"] == {}
    assert output["issues"]


def test_output_is_idempotent_and_backbones_do_not_become_sequences(tmp_path):
    candidate = {"candidate_id": "design", "source_tool": "rfdiffusion_generator", "stage": "generation.backbone",
                 "binder_length": "80-120", "metrics": {"binder_residue_count": 90}, "raw_metadata": {"native": "retained"}}
    row = normalize_candidate(candidate)
    assert row == normalize_candidate(row)
    assert row["normalized_output"]["kind"] == "backbone_proposal"
    assert "design.sequence_length" not in row["normalized_output"]["metrics"]
    write_candidates(tmp_path, "rfdiffusion_generator", [row])
    assert read_candidates(tmp_path)[0] == row
    assert row["metrics"] == candidate["metrics"]


def test_database_common_fields_filter_sort_and_sources_from_old_observations(tmp_path):
    user = {"id": "owner", "role": "user"}
    db = MolecularDatabase(tmp_path)
    imported = db.import_records(user, molecule_imports.parse_sequences(">one\nACD\n>two\nACE"))
    one, two = imported["molecule_ids"]
    group = imported["group_id"]
    db.add_observation(group, one, "genie3_generator", "old", {"binder_ca_radius_of_gyration": 12}, created_at="2026-01-01")
    db.add_observation(group, two, "rfdiffusion3_foundry", "old", {"binder_ca_radius_of_gyration": 9}, created_at="2026-01-01")
    db.add_observation(group, one, "protein_mpnn", "new", {}, context={"candidate": {
        "source_tool": "protein_mpnn", "stage": "sequence_design", "binder_sequence": "ACD",
        "raw_metadata": {"fasta_header": "id=1, overall_confidence=0.5, seq_rec=0.25"}}}, created_at="2026-02-01")
    field = PREFIX + "geometry.binder.ca_rg"
    assert field in db.metric_fields(user)
    assert [row["name"] for row in db.list_molecules(user, sort=field, descending=False)["items"]] == ["two", "one"]
    # A newer sequence evaluation must not borrow geometry from an older model.
    assert db.list_molecules(user, metric_filters=[{"field": field, "operator": ">", "value": 10}])["items"] == []
    rows = db.list_molecules(user, evaluation_filters={"run_id": "old"}, metric_filters=[{"field": field, "operator": ">", "value": 10}])["items"]
    assert len(rows) == 1
    assert rows[0]["scores"][field] == 12
    assert rows[0]["score_sources"][field]["run_id"] == "old"
    assert "common.sequence.recovery_designed" not in rows[0]["scores"]
    latest = db.list_molecules(user, molecule_ids=[one])["items"][0]
    assert latest["scores"]["common.sequence.recovery_designed"] == .25
    assert rows[0]["scores"]["genie3_generator.binder_ca_radius_of_gyration"] == 12
    assert rows[0]["observations"][0]["metrics"] == {}
    # A private observation about the same identity must not leak via common fields.
    private = db.import_records({"id": "other", "role": "user"}, molecule_imports.parse_sequences("ACD"))
    db.add_observation(private["group_id"], one, "secret", "private", {"binder_radius_of_gyration": 99})
    assert "common.geometry.binder.heavy_atom_rg" not in db.metric_fields(user)
    assert "common.geometry.binder.heavy_atom_rg" not in db.list_molecules(user)["items"][0]["scores"]


def test_pipeline_filter_accepts_same_common_field_for_different_tools():
    from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import filter_candidates
    candidates = [{"candidate_id": "a", "source_tool": "genie3_generator", "metrics": {"binder_ca_radius_of_gyration": 10}},
                  {"candidate_id": "b", "source_tool": "rfdiffusion3_foundry", "metrics": {"binder_ca_radius_of_gyration": 20}}]
    result = filter_candidates(candidates, {"metric": "common.geometry.binder.ca_rg", "operator": "<", "threshold": 15,
                                           "passing_only": False, "higher_is_better": False, "limit": 0})
    assert [c["candidate_id"] for c in result] == ["a"]


def test_every_design_tile_is_covered_by_shared_contract():
    import ast
    from pathlib import Path
    from mn_molecular_design_studio.output_vocabulary import GENERATOR_WORKFLOWS, vocabulary
    tree = ast.parse((Path(__file__).parents[1] / "mn_molecular_design_studio/run_app.py").read_text())
    tiles = next(ast.literal_eval(node.value) for node in tree.body if isinstance(node, ast.Assign)
                 and any(isinstance(target, ast.Name) and target.id == "TOOL_CATALOG_WORKFLOW_IDS" for target in node.targets))
    assert {workflow for workflow in tiles if workflow.startswith("design.")} == set(GENERATOR_WORKFLOWS)
    assert set(vocabulary()["generator_workflows"]) == set(GENERATOR_WORKFLOWS)


def test_output_dictionary_api_and_accessible_metric_fields(tmp_path, monkeypatch):
    from fastapi import FastAPI
    from fastapi.testclient import TestClient
    from mn_molecular_design_studio import molecular_database_api
    user = {"id": "owner", "role": "user"}
    db = MolecularDatabase(tmp_path)
    imported = db.import_records(user, molecule_imports.parse_sequences("ACD"))
    db.add_observation(imported["group_id"], imported["molecule_ids"][0], "rfdiffusion3_foundry", "run", {"binder_ca_radius_of_gyration": 11})
    monkeypatch.setattr(molecular_database_api, "MolecularDatabase", lambda: db)
    monkeypatch.setattr(molecular_database_api, "sync_workspace", lambda user: [])
    app = FastAPI(); molecular_database_api.install(app, lambda: user)
    client = TestClient(app)
    response = client.get("/api/v1/molecules/output-vocabulary")
    assert response.status_code == 200
    assert response.json()["version"] == "studio.output.v1"
    assert "common.geometry.binder.ca_rg" in client.get("/api/v1/molecules/metrics").json()["items"]
