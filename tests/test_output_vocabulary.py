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


def test_mol1_native_measurements_keep_complex_and_binder_scopes_distinct():
    result = candidate_output({
        "source_tool": "optimized_prediction_mol1", "stage": "complex_refolding",
        "target_chains": ["A"], "binder_chains": ["B"],
        "metrics": {"ptm": .81, "iptm": .72, "plddt_mean": 79.0, "complex_plddt": 79.0,
                    "binder_plddt": 84.0, "mol1_model_plddt_mean": 80.0},
    })
    assert result["metrics"]["confidence.validation.complex.ptm"] == .81
    assert result["metrics"]["confidence.validation.complex.iptm"] == .72
    assert result["metrics"]["confidence.validation.structure.mean_ca_plddt"] == 79.0
    assert result["metrics"]["confidence.validation.bound_binder.plddt"] == 84.0
    assert result["metrics"]["confidence.validation.model.mol1_plddt_mean"] == 80.0


def test_mol1_monomer_does_not_report_complex_or_interface_measurements():
    result = candidate_output({
        "candidate_id": "minisog-monomer", "source_tool": "optimized_prediction_mol1",
        "stage": "complex_refolding", "metrics": {
            "ptm": .83376, "iptm": 0.0, "plddt_mean": 88.314, "complex_plddt": 88.314,
        },
        "raw_metadata": {
            "api_engine": "mol1", "protein_chains": [{"id": "A", "sequence": "ACDE"}],
            "confidence_outputs": {"plddt": 84.976, "plddt_mean": 84.976, "iptm": 0.0},
        },
    })
    assert result["stage"] == "monomer_refolding"
    assert result["kind"] == "refolded_monomer"
    assert result["metrics"]["confidence.validation.structure.mean_ca_plddt"] == 88.314
    assert result["metrics"]["confidence.validation.model.mol1_plddt_mean"] == 84.976
    assert result["metrics"]["confidence.validation.monomer.ptm"] == .83376
    assert "confidence.validation.complex.iptm" not in result["metrics"]
    assert "confidence.validation.complex.ptm" not in result["metrics"]
    assert "confidence.validation.bound_binder.plddt" not in result["metrics"]


def test_bindcraft2_scaffold_and_induced_fit_metrics_keep_their_native_meaning():
    result = candidate_output({
        "source_tool": "BindCraft2", "stage": "complex_refolding",
        "metrics": {
            "Scaffold_Sequence_Retained_Fraction": 1.0,
            "Scaffold_Framework_RMSD": 1.8,
            "Off_Paratope_Contact_Fraction": 0.1,
            "Framework_Packing_Fraction": 0.4,
            "Induced_Fit_Interface_RMSD": 5.2,
            "Induced_Fit_TM": 0.55,
            "Interface_Residues": 11,
            "Backbone_Clashes": 0,
            "i_pAE": 0.33,
        },
    })
    metrics = result["metrics"]
    assert metrics["sequence.native.bindcraft2.scaffold_sequence_retained_fraction"] == 1.0
    assert metrics["geometry.native.bindcraft2.scaffold_framework_rmsd"] == 1.8
    assert metrics["interface.native.bindcraft2.off_paratope_contact_fraction"] == 0.1
    assert metrics["structure.native.bindcraft2.framework_packing_fraction"] == 0.4
    assert metrics["geometry.native.bindcraft2.induced_fit_interface_rmsd"] == 5.2
    assert metrics["confidence.native.bindcraft2.induced_fit_tm"] == 0.55
    assert metrics["interface.native.bindcraft2.binder_residue_count_4a"] == 11
    assert metrics["geometry.native.bindcraft2.backbone_ca_clashes_2_5a"] == 0
    assert metrics["confidence.native.bindcraft2.interface.ipae_normalized"] == 0.33


def test_cofold_metrics_keep_measured_values_units_and_separate_analysis_sources():
    from mn_molecular_design_studio.output_vocabulary import metric_value, vocabulary

    prediction = candidate_output({
        "candidate_id": "pdl1-fixture", "source_tool": "optimized_prediction_boltz2",
        "stage": "complex_refolding", "binder_sequence": "ACDE",
        "metrics": {"iptm": .67, "ipae": 4.8, "binder_plddt": 81.5,
                    "target_aligned_binder_rmsd": 2.1},
        "raw_metadata": {"prediction_protocol": {"engine": "boltz2", "samples": 1}},
    })
    assert prediction["metrics"]["confidence.validation.complex.iptm"] == .67
    assert prediction["metrics"]["confidence.validation.interface.mean_pae"] == 4.8
    assert prediction["metrics"]["confidence.validation.bound_binder.plddt"] == 81.5
    assert prediction["metrics"]["geometry.validation.binder.target_aligned_rmsd"] == 2.1
    assert prediction["provenance"]["confidence.validation.interface.mean_pae"]["unit"] == "Å"
    assert prediction["provenance"]["geometry.validation.binder.target_aligned_rmsd"]["method"].startswith("Studio matched target Cα")
    assert metric_value({"source_tool": "optimized_prediction_boltz2", "stage": "complex_refolding",
                         "metrics": {"iptm": .67, "ipae": 4.8}}, "common.confidence.validation.interface.mean_pae") == 4.8
    assert metric_value({"source_tool": "optimized_prediction_boltz2", "stage": "complex_refolding",
                         "metrics": {"iptm": .67, "ipae": 4.8}}, "common.confidence.validation.complex.iptm") == .67

    ipsae = candidate_output({"source_tool": "ipsae", "stage": "analysis", "metrics": {"ipsae_min": .42}})
    contacts = candidate_output({"source_tool": "biotite_interface", "stage": "analysis",
                                 "metrics": {"qc_interface_contact_pairs": 57}})
    assert ipsae["metrics"]["interface.ipsae.minimum"] == .42
    assert contacts["metrics"]["interface.contacts"] == 57
    definitions = {row["key"]: row for row in vocabulary()["measures"]}
    assert definitions["confidence.validation.interface.mean_pae"]["range_semantics"].startswith("Physical/value-domain")
    assert definitions["interface.ipsae.minimum"]["direction"] == "higher"

    native_routes = [
        ("pxdesign", "geometry.native.pxdesign.bound_unbound_rmsd", True),
        ("genie3", "geometry.native.genie3.binder_self_consistency_rmsd", True),
        ("protpardelle_1c", "geometry.native.protpardelle.binder_scaffold_rmsd", False),
        ("proteina_complexa", "geometry.native.proteina_complexa.binder_scaffold_rmsd", True),
    ]
    for tool, rmsd_field, has_iptm in native_routes:
        metrics = {"binder_plddt": 82.0, "ipae": 5.2, "binder_rmsd": 2.3}
        if has_iptm:
            metrics["iptm"] = .72
        if tool == "proteina_complexa":
            metrics.update(ipsae_min=.41, ipsae_max=.57, ipsae_avg=.49)
        native = candidate_output({"source_tool": tool, "stage": "complex_refolding", "metrics": metrics})
        assert native["metrics"]["confidence.validation.bound_binder.plddt"] == 82.0
        assert native["metrics"]["confidence.validation.interface.mean_pae"] == 5.2
        assert native["metrics"][rmsd_field] == 2.3
        assert native["provenance"]["confidence.validation.interface.mean_pae"]["method"]
        if has_iptm:
            assert native["metrics"]["confidence.validation.complex.iptm"] == .72
        if tool == "proteina_complexa":
            assert native["metrics"]["interface.ipsae.minimum"] == .41


def test_tool_scoped_measurement_value_wins_when_equivalent_aliases_are_present():
    row = candidate_output({
        "source_tool": "pxdesign", "stage": "complex_refolding",
        # PXDesign's engine-scoped interface PAE takes precedence over the
        # generic binder_pae alias when the source values agree within rounding.
        "metrics": {"binder_pae": 5.2, "ipae": 5.200000005},
    })
    field = "confidence.validation.interface.mean_pae"
    assert row["metrics"][field] == 5.200000005
    assert row["provenance"][field]["method"] == "PXDesign native AF2-IG"


def test_rfd4_oracle_aliases_preserve_native_scale_and_require_rfd4_generation():
    native = {"oracle.overall_plddt": 83.5, "oracle.overall_pae": 3.4178,
              "oracle.overall_ptm": 45.0, "oracle.plddt_coi": 83.5,
              "oracle.pae_coi": 3.4178, "oracle.ptm_coi": 45.0,
              "binder_ca_radius_of_gyration": 12.5}
    candidate = {"source_tool": "rfd4_proteina", "stage": "generation.backbone_sequence",
                 "metrics": native}
    result = candidate_output(candidate)
    assert result["metrics"]["confidence.generator.complex.plddt"] == 83.5
    assert result["metrics"]["confidence.generator.complex.pae"] == 3.4178
    assert result["metrics"]["confidence.generator.complex.ptm"] == .45
    assert result["metrics"]["confidence.generator.chain_of_interest.plddt"] == 83.5
    assert result["metrics"]["confidence.generator.chain_of_interest.pae"] == 3.4178
    assert result["metrics"]["confidence.generator.chain_of_interest.ptm"] == .45
    assert result["metrics"]["geometry.binder.ca_rg"] == 12.5
    assert result["provenance"]["confidence.generator.complex.ptm"]["aliases"][0] == {
        "source_field": "oracle.overall_ptm", "raw_value": 45.0, "value": .45,
        "conversion": "percent_to_fraction", "method": "RFD4 AF3 oracle"}
    assert native["oracle.overall_ptm"] == 45.0
    for out_of_scope in (candidate_output({**candidate, "source_tool": "other"}),
                         candidate_output({**candidate, "stage": "complex_refolding"})):
        assert not any(key.startswith("confidence.generator.") for key in out_of_scope["metrics"])


def test_database_projection_keeps_shared_geometry_guard_from_candidate_context():
    native = {"binder_radius_of_gyration": 13.6}
    candidate = {"source_tool": "rfd4_proteina", "stage": "generation.backbone_sequence",
                 "metrics": {**native, "binder_geometry_status": "ok"}}
    observation_context = {"candidate": {"source_tool": "rfd4_proteina",
                                          "stage": "generation.backbone_sequence",
                                          "binder_geometry_status": "ok"}}
    assert candidate_output(candidate)["metrics"]["geometry.binder.heavy_atom_rg"] == 13.6
    assert normalize_metrics(native, "rfd4_proteina", observation_context)["metrics"]["geometry.binder.heavy_atom_rg"] == 13.6


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
