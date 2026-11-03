from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from mn_molecular_design_studio import molecule_imports
from mn_molecular_design_studio.molecular_database import MolecularDatabase


USER = {"id": "e" * 32, "role": "user", "_cluster_project_ids": []}
OTHER_USER = {"id": "f" * 32, "role": "user", "_cluster_project_ids": []}
AA3 = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY",
       "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN",
       "P": "PRO", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR"}


def _write_structure(path: Path, sequence: str = "ACNGMFWY", *, missing=(), sg_coordinates=None,
                     add_ssbond=False, partner_sequence="", partner_y_offset=11.0,
                     insertion_codes=None, b_factor=45.0):
    insertion_codes = insertion_codes or {}
    sg_coordinates = sg_coordinates or {}
    lines = ["HEADER    CONTROLLED SEQUENCE LIABILITY FIXTURE"]
    for chain_id, seq in (("A", sequence), ("B", partner_sequence)):
        if not seq:
            continue
        lines.append(f"SEQRES   1 {chain_id} {len(seq):4d}  " + " ".join(AA3[letter] for letter in seq))
    if add_ssbond:
        row = [" "] * 80
        row[0:6] = "SSBOND"
        row[7:10] = f"{1:3d}"
        row[11:14] = "CYS"
        row[15] = "A"
        row[17:21] = f"{1:4d}"
        row[21] = " "
        row[25:28] = "CYS"
        row[29] = "A"
        row[31:35] = f"{2:4d}"
        row[35] = " "
        lines.append("".join(row))
    serial = 1
    for chain_id, seq, y_offset in (("A", sequence, 0.0), ("B", partner_sequence, partner_y_offset)):
        if not seq:
            continue
        for index, letter in enumerate(seq, start=1):
            if (chain_id, index) in missing:
                continue
            number, insertion = insertion_codes.get((chain_id, index), (index, " "))
            x = (index - 1) * 3.8
            atoms = [
                ("N", (x - 1.2, y_offset, 0.0), "N"),
                ("CA", (x, y_offset, 0.0), "C"),
                ("C", (x + 1.3, y_offset, 0.1), "C"),
                ("O", (x + 1.7, y_offset + 0.8, 0.1), "O"),
                ("CB", (x, y_offset, 1.5), "C"),
            ]
            if letter == "C" and not (chain_id, index, "no_sg") in sg_coordinates:
                default_sg = (x, y_offset, 3.0)
                atoms.append(("SG", sg_coordinates.get((chain_id, index), default_sg), "S"))
            for atom_name, coord, element in atoms:
                lines.append(
                    f"ATOM  {serial:5d} {atom_name:>4s} {AA3[letter]:>3s} {chain_id}{number:4d}{insertion:1s}   "
                    f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}{1.0:6.2f}{b_factor:6.2f}          {element:>2s}  "
                )
                serial += 1
        lines.append(f"TER   {serial:5d}      {AA3[seq[-1]]} {chain_id}{len(seq):4d}")
    path.write_text("\n".join([*lines, "END", ""]), encoding="utf-8")
    return path


def _identity(path: Path, model_id="model-1"):
    from mn_molecular_design_studio.evaluations import file_digest
    return {"molecule_id": "molecule-1", "structure_id": "structure-1", "group_id": "group-1",
            "structure_digest": file_digest(path), "model_id": model_id, "evaluation_id": ""}


def _calculate(path: Path, *, chains=("A",), partners=(), **kwargs):
    from mn_molecular_design_studio.complex_analysis_liabilities import calculate_sequence_liabilities
    return calculate_sequence_liabilities(
        path, component_chains=list(chains), partner_chains=list(partners), model_id="model-1",
        source_identity=_identity(path), **kwargs,
    )


def test_sequence_motifs_cover_terminal_dipeptides_and_unknown_residues_break_alerts():
    from mn_molecular_design_studio.complex_analysis_liabilities import _residue_findings

    rows = [{"chain_id": "A", "author_number": index + 1, "insertion_code": "", "author_residue": str(index + 1),
             "residue_name": "ALA", "amino_acid": aa, "selector": f"A:{index + 1}", "coordinates_present": False}
            for index, aa in enumerate("NGDSNPS")]
    findings = _residue_findings("NGDSNPS", rows, {}, set(), exposure_status="unavailable")
    kinds = {(row["kind"], row["motif"]) for row in findings}
    assert ("deamidation_ng", "NG") in kinds
    assert ("isomerization_ds", "DS") in kinds
    assert not any(row["kind"] == "n_glycosylation_sequon" for row in findings)  # N-P-S is excluded

    unknown_rows = [dict(row, amino_acid=aa) for row, aa in zip(rows, "NXSGMFW")]
    unknown = _residue_findings("NXSGMFW", unknown_rows, {}, set(), exposure_status="unavailable")
    assert not any(row["kind"] == "n_glycosylation_sequon" for row in unknown)
    assert not any(row["kind"] == "deamidation_ng" for row in unknown)


def test_saved_sequence_properties_and_confidence_are_kept_method_specific(tmp_path):
    path = _write_structure(tmp_path / "properties.pdb", "ACDEFGH")
    from mn_molecular_design_studio.complex_analysis_liabilities import calculate_sequence_liabilities
    from mn_molecular_design_studio.engine.workflow_engine.workflows.quality_runner import sequence_metrics

    saved_qc = {"evaluator": "sequence_qc", "evaluation_id": "qc-1", "run_id": "run-qc",
                "sequence": "ACDEFGH", "source_digest": _identity(path)["structure_digest"],
                "chain_id": "A", "metrics": sequence_metrics("ACDEFGH", 7.4)}
    confidence = {"evaluation_id": "pred-1", "evaluator": "boltz2", "run_id": "run-pred",
                  "structure_digest": _identity(path)["structure_digest"], "model_id": "model-1",
                  "native_scalar_metrics": {"iptm": 0.62}, "scope": "global scalar"}
    result = calculate_sequence_liabilities(
        path, component_chains=["A"], partner_chains=[], model_id="model-1", source_identity=_identity(path),
        reused_evidence=[saved_qc], saved_confidence_context=confidence,
    )
    assert result["sequence_properties"]["A"]["status"] == "reused_saved_evaluation"
    assert result["sequence_properties"]["A"]["metrics"]["qc_sequence_pi"] == pytest.approx(saved_qc["metrics"]["qc_sequence_pi"])
    wrong_ph = {**saved_qc, "metrics": {**saved_qc["metrics"], "qc_sequence_ph": 8.0}}
    from mn_molecular_design_studio.complex_analysis_liabilities import _safe_protein_metrics
    assert _safe_protein_metrics("ACDEFGH", wrong_ph)["status"] == "measured"
    saved_conf = next(row for row in result["confidence_evidence"] if row.get("status") == "evaluation_context_only")
    assert saved_conf["evaluation_id"] == "pred-1"
    assert saved_conf["metrics"] == {"iptm": 0.62}
    assert "not residue-local" in saved_conf["scope"]
    assert not any(row.get("status") == "residue_local_values" for row in result["confidence_evidence"])
    assert any(row["data"] == "residue-local structural confidence" for row in result["missing_data"])

    mismatched = dict(confidence, structure_digest="stale")
    result = calculate_sequence_liabilities(
        path, component_chains=["A"], partner_chains=[], model_id="model-1", source_identity=_identity(path),
        saved_confidence_context=mismatched,
    )
    assert any(row.get("status") == "not_reused_stale_or_mismatched" for row in result["confidence_evidence"])


def test_cyclic_format_and_unresolved_cysteine_are_not_called_free_termini_or_unpaired(tmp_path):
    path = _write_structure(tmp_path / "incomplete.pdb", "CNGM", missing={("A", 1)})
    result = _calculate(path, component_format="cyclic_construct")
    assert result["sequence_coverage"]["A"]["complete_coordinate_coverage"] is False
    assert result["sequence_coverage"]["A"]["terminal_assessment"]["status"] == "not_applicable"
    assert any(item["kind"] == "cysteine_sg_coordinates_missing" for item in result["findings"])
    assert not any(item["kind"] == "unpaired_cysteine" for item in result["findings"])
    assert any(row["data"] == "free termini" and row["status"] == "not_applicable" for row in result["missing_data"])


def test_ambiguous_disulfide_matching_is_not_forced_and_recorded_pairs_take_precedence(tmp_path):
    from mn_molecular_design_studio.complex_analysis_liabilities import (
        _disulfide_assessment, _optimal_matching, _selected_residues,
    )
    from mn_molecular_design_studio.complex_analysis_tag_fusion import _read_structure

    nodes = [("A", index, "") for index in range(1, 5)]
    edges = {(nodes[i], nodes[j]): 3.2 for i in range(4) for j in range(i + 1, 4)}
    _score, count, _pairs = _optimal_matching(nodes, edges)
    assert count > 1

    ambiguous_path = _write_structure(tmp_path / "ambiguous.pdb", "CCCC", sg_coordinates={
        ("A", 1): (0, 0, 0), ("A", 2): (0, 0, 0), ("A", 3): (0, 0, 0), ("A", 4): (0, 0, 0),
    })
    _structure, model = _read_structure(ambiguous_path)
    residue_index = _selected_residues(model, ["A"])
    inferred = _disulfide_assessment(ambiguous_path, model, ["A"], residue_index)
    assert inferred["status"] == "inferred_geometry"
    assert not inferred["pairs"]
    assert all(item["kind"] == "ambiguous_disulfide_pairing" for item in inferred["findings"])

    recorded_path = _write_structure(tmp_path / "recorded.pdb", "CC", add_ssbond=True)
    _structure, recorded_model = _read_structure(recorded_path)
    recorded_index = _selected_residues(recorded_model, ["A"])
    recorded = _disulfide_assessment(recorded_path, recorded_model, ["A"], recorded_index)
    assert recorded["basis"] == "explicit_structure_record"
    assert len(recorded["pairs"]) == 1
    assert recorded["pairs"][0]["basis"] == "recorded"


def test_interchain_disulfide_is_recognized_and_unresolved_cys_stays_unknown_with_records(tmp_path):
    from mn_molecular_design_studio.complex_analysis_liabilities import _disulfide_assessment, _selected_residues
    from mn_molecular_design_studio.complex_analysis_tag_fusion import _read_structure

    interchain_path = _write_structure(
        tmp_path / "interchain.pdb", "C", partner_sequence="C", partner_y_offset=0.0,
        sg_coordinates={("A", 1): (0.0, 0.0, 3.0), ("B", 1): (2.05, 0.0, 3.0)},
    )
    _structure, interchain_model = _read_structure(interchain_path)
    selected_index = _selected_residues(interchain_model, ["A"])
    interchain = _disulfide_assessment(interchain_path, interchain_model, ["A"], selected_index)
    assert interchain["status"] == "inferred_geometry"
    assert len(interchain["pairs"]) == 1
    assert {interchain["pairs"][0]["first"][0], interchain["pairs"][0]["second"][0]} == {"A", "B"}
    assert not any(row["kind"] == "unpaired_cysteine" for row in interchain["findings"])

    recorded_missing_path = _write_structure(
        tmp_path / "recorded-with-missing.pdb", "CCC", add_ssbond=True,
        sg_coordinates={("A", 3, "no_sg"): (0.0, 0.0, 0.0)},
    )
    _structure, recorded_missing_model = _read_structure(recorded_missing_path)
    recorded_missing_index = _selected_residues(recorded_missing_model, ["A"])
    recorded_missing = _disulfide_assessment(
        recorded_missing_path, recorded_missing_model, ["A"], recorded_missing_index,
    )
    third_cys = [row for row in recorded_missing["findings"] if row["author_number"] == 3]
    assert len(third_cys) == 1
    assert third_cys[0]["kind"] == "cysteine_sg_coordinates_missing"


def test_missing_dssp_and_sasa_are_reported_unknown_not_favorable(tmp_path, monkeypatch):
    from mn_molecular_design_studio import complex_analysis_tag_fusion
    path = _write_structure(tmp_path / "missing-sasa.pdb", "ACMFWY")
    monkeypatch.setattr(complex_analysis_tag_fusion, "_dssp_data", lambda *_args: ({}, {"status": "unavailable", "reason": "fixture"}))
    monkeypatch.setattr(complex_analysis_tag_fusion, "_sasa_data", lambda *_args: ({}, {"status": "unavailable", "reason": "fixture"}))
    result = _calculate(path)
    oxidation = [item for item in result["findings"] if item["kind"] in {
        "methionine_oxidation_susceptibility", "tryptophan_oxidation_susceptibility"}]
    assert oxidation
    assert all(item["basis"] == "sequence_only" for item in oxidation)
    assert all(item["evidence"]["relative_sasa"] is None for item in oxidation)
    assert {row["data"] for row in result["missing_data"]} >= {
        "DSSP relative accessibility / secondary structure", "absolute residue SASA"}


def test_ligand_interface_evidence_requires_exact_saved_pose_identity(tmp_path):
    from mn_molecular_design_studio.complex_analysis_liabilities import _matching_saved_evidence

    artifact = tmp_path / "plip.csv"
    artifact.write_text("protein_chain,protein_residue_number,interaction_type\nA,10,hydrogen bond\n")
    ligand = {"chain": "L", "residue_key": "L:101", "ccd": "LIG", "smiles": "C[C@H](O)N"}
    saved = [{"evaluator": "plip", "evaluation_id": "eval-1", "run_id": "run-1",
              "ligand_selection": dict(ligand), "artifacts": [{"name": "interactions.csv", "path": str(artifact)}]}]
    residues, provenance, _confidence = _matching_saved_evidence(saved, ["A"], [], interface_ligand_selection=ligand)
    assert residues == {("A", 10, "")}
    assert provenance["ligand_selection"] == ligand
    residues, provenance, _confidence = _matching_saved_evidence(
        saved, ["A"], [], interface_ligand_selection={**ligand, "residue_key": "L:102"},
    )
    assert not residues
    assert provenance["status"] == "unavailable"


def test_database_source_acl_and_exact_saved_structure_round_trip(tmp_path, monkeypatch):
    from mn_molecular_design_studio.complex_analysis_sources import apply_database_analysis_source, build_database_analysis_source
    from mn_molecular_design_studio.engine.workflow_engine.workflows.analysis import run_sequence_liability_contract
    from mn_molecular_design_studio.evaluations import evaluation_descriptor
    from mn_molecular_design_studio.molecular_database_ingest import ingest_run
    from mn_molecular_design_studio.workspace import mark_studio_job

    workspace = tmp_path / "real-1iep-workspace"
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(workspace))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(workspace))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(workspace / "workdir" / "runs"))
    database = MolecularDatabase(workspace / "database")
    real_structure = Path(__file__).parents[1] / "mn_molecular_design_studio" / "assets" / "vina" / "1iep.pdb"
    record = molecule_imports.parse_structure(real_structure.read_bytes(), real_structure.name)[0]
    imported = database.import_records(USER, [record], name="Saved experimental 1IEP structure")
    saved_record = database.list_molecules(USER, molecule_ids=imported["molecule_ids"], collection="all")["items"][0]
    structure = saved_record["representations"][0]
    model_id = str((structure.get("context") or {}).get("model_id") or structure["id"])
    source = {"molecule_id": saved_record["id"], "structure_id": structure["id"],
              "group_id": imported["group_id"], "model_id": model_id,
              "binder_chains": ["A"], "target_chains": []}
    adapter = build_database_analysis_source(USER, source, workflow_id="analysis.sequence_liability",
                                             project_id=None, database=database)
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import (
        WORKFLOWS, validate_workflow_request, workflow_schema,
    )
    assert "analysis.sequence_liability" in WORKFLOWS
    schema_names = {row["name"] for row in workflow_schema("analysis.sequence_liability")["parameters"]}
    assert {"component_format", "b_factor_interpretation", "interface_ligand_selection"} <= schema_names
    request = apply_database_analysis_source({
        "schema_version": 1, "workflow": "analysis.sequence_liability", "database_source": source,
        "parameters": {"component_format": "linear_protein", "b_factor_interpretation": "uninterpreted"},
    }, USER, None)
    request_path = tmp_path / "sequence-liability-workflow-request.json"
    request_path.write_text(json.dumps({key: request[key] for key in ("schema_version", "workflow", "parameters")}))
    validated = validate_workflow_request(request_path)
    assert validated["valid"] is True, validated
    assert validated["parameters"]["interface_ligand_selection"] == {}
    with pytest.raises(PermissionError):
        build_database_analysis_source(OTHER_USER, source, workflow_id="analysis.sequence_liability",
                                       project_id=None, database=database)

    run_dir = run_sequence_liability_contract(
        Path(adapter["source_run_dir"]), Path(adapter["candidates_jsonl"]),
        candidate_ids=[adapter["candidate_id"]],
        selected_structure_path=Path(adapter["selected_structure_path"]),
        selected_model_id=adapter["selected_model_id"], binder_chains=["A"], target_chains=[],
    )
    payload_path = run_dir / "artifacts" / "analysis" / "sequence_liability" / "sequence_liability.json"
    result = json.loads(payload_path.read_text(encoding="utf-8"))
    assert result["source"]["molecule_id"] == saved_record["id"]
    assert result["source"]["structure_id"] == structure["id"]
    assert result["source"]["structure_digest"] == adapter["source"]["structure_digest"]
    assert result["sequence_properties"]["A"]["metrics"]["qc_sequence_length"] > 100
    assert result["findings"]
    assert all(item["basis"] in {"sequence_only", "structure_supported", "sequence_and_structure"}
               for item in result["findings"])

    mark_studio_job({"run_dir": str(run_dir)}, project_id=None, owner_id=USER["id"], allow_private_source=True)
    assert ingest_run(run_dir, database=database, force=True) == 1
    reopened = database.list_molecules(USER, molecule_ids=[saved_record["id"]], collection="all")["items"][0]
    observation = next(item for item in reopened["observations"] if item.get("namespace") == "sequence_liability")
    assert evaluation_descriptor(observation)["evaluator"] == "sequence_liability"
    assert observation["context"]["representation_id"] == structure["id"]
    saved_result = observation["context"]["candidate"]["raw_metadata"]["sequence_liability"]
    assert saved_result["measurement_counts"] == result["measurement_counts"]

    from mn_molecular_design_studio import complex_analysis_ui, molecular_database
    monkeypatch.setattr(molecular_database, "MolecularDatabase", lambda: database)
    inspector = complex_analysis_ui._database_inspector_data(
        USER, reopened, structure, database.structure_file(USER, structure["id"]),
        saved_record["components"], [(observation, evaluation_descriptor(observation))],
    )
    shown = next(item for item in inspector["evaluations"] if item.get("evaluator") == "sequence_liability")
    assert shown["sequence_liability_result"]["measurement_counts"] == result["measurement_counts"]
    saved_artifact = next(row for row in shown["artifacts"] if row.get("name") == "sequence_liability.json")
    assert Path(saved_artifact["path"]).resolve() == payload_path.resolve()
    highlights = [row for row in inspector["evidence"] if row.get("evaluator") == "sequence_liability"]
    assert highlights
    assert any(row.get("selectors") for row in highlights)

    # A second real saved protein structure is a negative control for motif
    # and free-Cys alerts. 1CRN carries explicit PDB SSBOND records.
    negative_path = Path(__file__).parents[1] / "mn_molecular_design_studio" / "assets" / "inverse_folding" / "1CRN.pdb"
    negative_record = molecule_imports.parse_structure(negative_path.read_bytes(), negative_path.name)[0]
    negative_import = database.import_records(USER, [negative_record], name="Saved experimental 1CRN control")
    negative_saved = database.list_molecules(USER, molecule_ids=negative_import["molecule_ids"], collection="all")["items"][0]
    negative_structure = negative_saved["representations"][0]
    negative_model = str((negative_structure.get("context") or {}).get("model_id") or negative_structure["id"])
    negative_adapter = build_database_analysis_source(USER, {
        "molecule_id": negative_saved["id"], "structure_id": negative_structure["id"],
        "group_id": negative_import["group_id"], "model_id": negative_model,
        "binder_chains": ["A"], "target_chains": [],
    }, workflow_id="analysis.sequence_liability", project_id=None, database=database)
    negative_run = run_sequence_liability_contract(
        Path(negative_adapter["source_run_dir"]), Path(negative_adapter["candidates_jsonl"]),
        candidate_ids=[negative_adapter["candidate_id"]],
        selected_structure_path=Path(negative_adapter["selected_structure_path"]),
        selected_model_id=negative_adapter["selected_model_id"], binder_chains=["A"], target_chains=[],
    )
    negative_result = json.loads((negative_run / "artifacts" / "analysis" / "sequence_liability" / "sequence_liability.json").read_text())
    assert negative_result["disulfides"]["basis"] == "explicit_structure_record"
    assert not any(item["kind"] == "unpaired_cysteine" for item in negative_result["findings"])
    assert not any(item["kind"] in {"deamidation_ng", "deamidation_ns", "isomerization_dg", "isomerization_ds",
                                     "n_glycosylation_sequon", "hydrophobic_stretch"}
                   for item in negative_result["findings"])


def test_pep_patch_native_residue_outputs_are_flattened_for_reuse():
    from mn_molecular_design_studio.engine.workflow_engine.workflows.quality_runner import flatten_pep_patch_residues

    native = [{"method": "hydrophobic", "patch_config": {"scale": "Ja"},
               "per_residue_patches": [
                   {"chain_id": "A", "residue_id": 12, "patch_id": 1, "value": 0.8},
                   {"chain_id": "A", "residue_id": 13, "patch_id": 1, "value": 0.9},
               ]}]
    rows = flatten_pep_patch_residues(native)
    assert len(rows) == 2
    assert rows[0]["method"] == "hydrophobic"
    assert rows[0]["patch_config"] == {"scale": "Ja"}
    assert rows[1]["residue_id"] == 13


def test_saved_pep_patch_annotations_link_author_residues_without_recalculating(tmp_path):
    from mn_molecular_design_studio.complex_analysis_liabilities import _pep_patch_findings

    artifact = tmp_path / "residue_patches.csv"
    artifact.write_text(
        'method,patch_config,chain_id,residue_id,patch_id,patch_labels,value\n'
        'hydrophobic,"{""scale"":""Ja""}",A,10,1,"{""1"":42.0}",0.81\n'
    )
    chain_rows = {"A": [{"chain_id": "A", "author_number": 10, "insertion_code": "B",
                         "author_residue": "10B", "residue_name": "TYR", "amino_acid": "Y",
                         "selector": "A:10^B", "sequence_index": 0}]}
    findings, summaries = _pep_patch_findings(
        [{"evaluator": "pep_patch", "evaluation_id": "pep-eval", "run_id": "pep-run",
          "method_version": "registered-v1",
          "artifacts": [{"name": "residue_patches.csv", "path": str(artifact)}]}],
        ["A"], chain_rows, interface_indices={"A": {0}},
    )
    assert not summaries
    assert len(findings) == 1
    assert findings[0]["selector"] == "A:10^B"
    assert findings[0]["interface_status"] == "in_interface"
    assert findings[0]["evidence"]["residue_patch_area_angstrom2"] == pytest.approx(42.0)
    assert findings[0]["evidence"]["saved_patch_metrics_evaluation_id"] == "pep-eval"


def test_inspector_renders_liability_provenance_and_highlights_exact_residue_selector():
    from streamlit.testing.v1 import AppTest

    item = {
        "evaluation": {"id": "evaluation-liability", "run_id": "run-liability"},
        "sequence_liability_result": {
            "assessment_status": "partial", "method": "method", "method_version": "v1",
            "source": {"structure_id": "structure-1", "structure_digest": "digest"},
            "selected_component_chains": ["A"], "selected_partner_chains": ["B"],
            "parameters": {"ph": 7.4}, "sequence_properties": {}, "sequence_coverage": {},
            "findings": [{"finding_id": "liability-00001", "severity": "warn", "basis": "sequence_only",
                          "kind": "deamidation_ng", "residue_ref": "A:10^A", "motif": "NG",
                          "interface_status": "unknown", "message": "Sequence susceptibility only."}],
            "missing_data": [{"data": "DSSP", "status": "unavailable", "reason": "Not installed."}],
            "confidence_evidence": [], "interpretation_limitations": ["No experimental modification measured."],
            "attribution": "BinderScope MIT adaptation",
        },
    }
    evidence = [{"id": "sequence_liability:liability-00001:evaluation-liability",
                 "evaluator": "sequence_liability", "label": "WARN · deamidation_ng · A:10^A",
                 "selectors": ["A:10^A"]}]
    script = f"""
from mn_molecular_design_studio.complex_analysis_ui import _render_inspector_liabilities
_render_inspector_liabilities([{item!r}], "fixture", {evidence!r})
    """
    app = AppTest.from_string(script, default_timeout=15).run()
    assert not app.exception
    assert any("deamidation_ng" in str(frame.value.to_dict()) for frame in app.dataframe)
    button = next(button for button in app.button if button.label == "Highlight this liability in Mol*")
    button.click().run()
    assert not app.exception
    assert app.session_state["fixture-inspector-highlights"] == ["A:10^A"]
