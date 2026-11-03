from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest


USER = {"id": "a" * 32, "role": "user", "_cluster_project_ids": []}
OTHER_USER = {"id": "b" * 32, "role": "user", "_cluster_project_ids": []}
_AA3 = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY", "H": "HIS"}


def _write_complex(path: Path, *, missing_n: bool = False, far_partner: bool = False,
                   n_terminal_insertion: bool = False, low_confidence: bool = False) -> Path:
    sequence = "ACDEFGH"
    lines = ["HEADER    STUDIO SYNTHETIC GEOMETRY FIXTURE"]
    for chain in ("A", "B"):
        lines.append(f"SEQRES   1 {chain}    7  " + " ".join(_AA3[aa] for aa in sequence))
    serial = 1
    for chain in ("A", "B"):
        y_offset = 12.0 if chain == "A" else (42.0 if far_partner else 16.0)
        for index, amino_acid in enumerate(sequence):
            if chain == "A" and missing_n and index == 0:
                continue
            x = index * 3.8
            z = (0.0, 0.8, 0.0, -0.8, 0.0, 0.8, 0.0)[index]
            residue_atoms = (
                ("N", (x - 1.2, y_offset + 0.2, z - 0.2), "N"),
                ("CA", (x, y_offset, z), "C"),
                ("C", (x + 1.3, y_offset + 0.3, z + 0.2), "C"),
                ("O", (x + 1.7, y_offset + 1.1, z + 0.3), "O"),
                ("CB", (x, y_offset, z + 1.5), "C"),
            )
            b_factor = 52.0 if low_confidence and chain == "A" and index < 5 else 78.0
            for atom_name, (ax, ay, az), element in residue_atoms:
                lines.append(
                    f"ATOM  {serial:5d} {atom_name:>4s} {_AA3[amino_acid]:>3s} {chain}{index + 1:4d}    "
                    f"{ax:8.3f}{ay:8.3f}{az:8.3f}{1.0:6.2f}{b_factor:6.2f}          {element:>2s}  "
                )
                serial += 1
        lines.append(f"TER   {serial:5d}      {_AA3[sequence[-1]]} {chain}{len(sequence):4d}")
    if n_terminal_insertion:
        for index, line in enumerate(lines):
            if line.startswith("ATOM  ") and line[21:22] == "A" and line[22:26].strip() == "1":
                lines[index] = line[:22] + f"{10:4d}" + "B" + line[27:]
    path.write_text("\n".join([*lines, "END", ""]), encoding="utf-8")
    return path


def _identity(path: Path, *, evaluation_id: str = "") -> dict:
    from mn_molecular_design_studio.evaluations import file_digest

    return {
        "molecule_id": "molecule-1", "structure_id": "structure-1", "group_id": "group-1",
        "structure_digest": file_digest(path), "model_id": "model-1", "evaluation_id": evaluation_id,
    }


def _calculate(path: Path, **kwargs):
    from mn_molecular_design_studio.complex_analysis_tag_fusion import calculate_tag_fusion_suitability

    params = {
        "component_chain": "A", "interface_chains": ["B"], "model_id": "model-1",
        "source_identity": _identity(path), "termini_confirmed": True,
        "b_factor_interpretation": "plddt_0_100",
    }
    params.update(kwargs)
    return calculate_tag_fusion_suitability(path, **params)


def test_terminal_evidence_preserves_author_identity_confidence_scope_and_clearance(tmp_path, monkeypatch):
    from mn_molecular_design_studio import complex_analysis_tag_fusion as evaluator

    path = _write_complex(tmp_path / "two-chain-complex.pdb", n_terminal_insertion=True)
    # Exercise full structure parsing, contact detection, SASA, terminal
    # geometry, and confidence handling while keeping the grid-path solver
    # isolated in its own focused test below.
    monkeypatch.setattr(evaluator, "_tag_clearances", lambda *_args: {
        "his6": {"label": "His6", "kind": "flexible_affinity_tag",
                 "status": "clearance_path_without_added_linker", "minimum_linker_residues_tested": 0},
    })
    saved_confidence = {
        "evaluation_id": "evaluation-1", "evaluator": "boltz2", "run_id": "run-1",
        "structure_digest": hashlib.sha256(path.read_bytes()).hexdigest(), "model_id": "model-1",
        "native_scalar_metrics": {"iptm": 0.6},
    }
    result = _calculate(path, source_identity=_identity(path, evaluation_id="evaluation-1"),
                        saved_confidence_context=saved_confidence)

    assert result["assessment_status"] == "measured"
    assert result["interface"]["observed_interface"] is True
    assert result["interface"]["component_residue_count_contacting_partner"] > 0
    assert result["terminal_geometry"]["N"]["terminal_residue"]["selector"] == "A:10^B"
    assert result["terminal_geometry"]["C"]["terminal_residue"]["selector"] == "A:7^_"
    assert result["terminal_geometry"]["N"]["terminal_window_mean_plddt"] == pytest.approx(78.0)
    assert result["terminal_geometry"]["N"]["tag_clearance"]["his6"]["minimum_linker_residues_tested"] == 0
    assert result["saved_confidence_context"]["used_for_terminal_grading"] is False
    assert "not localized to terminal residues" in result["saved_confidence_context"]["scope"]
    assert result["source_structure_unchanged"] is True
    assert any("not experimental proof" in item for item in result["interpretation_limitations"])


def test_folded_fusion_clearance_pass_remains_distinct_from_interface_attachment_advisory(tmp_path, monkeypatch):
    from mn_molecular_design_studio import complex_analysis_tag_fusion as evaluator

    path = _write_complex(tmp_path / "terminal-interface.pdb")
    passing = {
        tag["key"]: {"label": tag["label"], "kind": tag["kind"],
                     "status": "clearance_cone_without_added_linker", "minimum_linker_residues_tested": 0,
                     "minimum_clearance_angstrom": 5.0}
        for tag in evaluator.TAGS
    }
    monkeypatch.setattr(evaluator, "_tag_clearances", lambda *_args: {key: dict(value) for key, value in passing.items()})
    result = _calculate(path)
    n_terminal = result["terminal_geometry"]["N"]
    assert n_terminal["contact_at_terminal_residue"] is True
    for key in ("gfp_snap", "fc_mbp"):
        clearance = n_terminal["tag_clearance"][key]
        assert clearance["status"] == "clearance_cone_without_added_linker"
        assert clearance["suitability_advisory"]["status"] == "not_recommended_terminal_interface_contact"
        assert "sampled geometric envelope may clear" in clearance["suitability_advisory"]["reason"]
        assert result["tag_comparison"][key]["preferred_end_advisory"] == "neither_end_suitable_under_method"


def test_declared_incomplete_terminus_is_not_scored_and_user_confirmation_does_not_override_it(tmp_path, monkeypatch):
    from mn_molecular_design_studio import complex_analysis_tag_fusion as evaluator

    path = _write_complex(tmp_path / "truncated-terminal.pdb", missing_n=True)
    monkeypatch.setattr(evaluator, "_tag_clearances", lambda *_args: {"his6": {"status": "clear"}})
    result = _calculate(path, termini_confirmed=True)

    assert result["terminal_completeness"]["declared_sequence"]["first_position_resolved"] is False
    assert result["terminal_geometry"]["N"]["assessment_status"] == "not_assessed_incomplete_terminal_structure"
    assert result["terminal_geometry"]["N"]["tag_clearance"] == {}
    assert result["terminal_geometry"]["C"]["tag_clearance"]
    assert result["assessment_status"] == "partial_terminal_evidence"


@pytest.mark.parametrize("component_format,reason", [
    ("cyclic_peptide", "Cyclic peptides"),
    ("multichain_construct", "multichain construct"),
])
def test_unsupported_component_formats_are_explicit_not_invalid(tmp_path, component_format, reason):
    path = _write_complex(tmp_path / f"{component_format}.pdb")
    result = _calculate(path, component_format=component_format)

    assert result["assessment_status"] == "not_applicable"
    assert reason in result["not_applicable_reason"]
    assert result["checks_not_performed"]


def test_unconfirmed_observed_ends_are_provisional_and_do_not_receive_tag_recommendations(tmp_path, monkeypatch):
    from mn_molecular_design_studio import complex_analysis_tag_fusion as evaluator

    path = _write_complex(tmp_path / "unconfirmed-construct.pdb")
    monkeypatch.setattr(evaluator, "_tag_clearances", lambda *_args: {"his6": {"status": "clear"}})
    result = _calculate(path, termini_confirmed=False)

    assert result["assessment_status"] == "provisional_construct_termini_unconfirmed"
    assert result["terminal_geometry"]["N"]["assessment_status"] == "provisional_observed_chain_endpoint"
    assert result["terminal_geometry"]["N"]["tag_clearance"] == {}
    assert result["tag_comparison"]["his6"]["preferred_end_advisory"] == "not_assessed"


def test_clearance_models_are_distinct_and_method_parameters_are_saved():
    import numpy as np
    from mn_molecular_design_studio.complex_analysis_tag_fusion import _tag_clearances

    result = _tag_clearances(np.zeros(3), np.asarray([1.0, 0.0, 0.0]), np.empty((0, 3)))
    assert result["his6"]["kind"] == "flexible_affinity_tag"
    assert result["twin_strep"]["tag_length_residues"] > result["his6"]["tag_length_residues"]
    assert result["gfp_snap"]["kind"] == "globular_fusion"
    assert result["fc_mbp"]["effective_spherical_radius_angstrom"] > result["gfp_snap"]["effective_spherical_radius_angstrom"]
    assert result["his6"]["status"] == "clearance_path_without_added_linker"
    assert result["gfp_snap"]["status"] == "clearance_cone_without_added_linker"


def test_inspector_terminal_measurement_rows_keep_numeric_and_categorical_values_arrow_safe():
    from mn_molecular_design_studio.complex_analysis_ui import (
        _display_dssp_assignment, _tag_fusion_terminal_metric_rows,
    )

    rows = _tag_fusion_terminal_metric_rows({
        "distance_to_selected_interface_angstrom": 6.666,
        "secondary_structure": "-",
        "terminal_window_secondary_structure": ["-", "H"],
        "flexibility_evidence": "not_available",
        "terminal_window_mean_plddt": None,
    })

    assert rows == [
        {"Evidence": "Distance to selected interface (Å)", "Stored value": "6.666"},
        {"Evidence": "Secondary structure", "Stored value": "-"},
        {"Evidence": "Flexibility evidence", "Stored value": "not_available"},
    ]
    assert all(isinstance(row["Stored value"], str) for row in rows)
    path_rows = _tag_fusion_terminal_metric_rows({"secondary_structure": "/runs/analysis/old-run/-"})
    assert path_rows == [{"Evidence": "Secondary structure", "Stored value": "-"}]
    assert _display_dssp_assignment("/runs/analysis/old-run/-") == "-"
    assert _display_dssp_assignment("/runs/analysis/old-run/not-a-DSSP-code") == "unavailable (invalid stored assignment)"


def test_geometry_and_clearance_are_rigid_transform_covariant_and_sasa_sampling_is_stable(tmp_path):
    from mn_molecular_design_studio import complex_analysis_tag_fusion as evaluator

    original = _write_complex(tmp_path / "rigid-original.pdb")
    rotation = np.asarray([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    translation = np.asarray([31.25, -17.5, 8.125])
    moved_lines = []
    for line in original.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            xyz = np.asarray([float(line[start:start + 8]) for start in (30, 38, 46)])
            moved = rotation @ xyz + translation
            line = line[:30] + "".join(f"{value:8.3f}" for value in moved) + line[54:]
        moved_lines.append(line)
    moved = tmp_path / "rigid-moved.pdb"
    moved.write_text("\n".join([*moved_lines, ""]), encoding="utf-8")

    baseline = _calculate(original)
    transformed = _calculate(moved, source_identity=_identity(moved))
    assert baseline["method_version"] == transformed["method_version"]
    for end in ("N", "C"):
        first = baseline["terminal_geometry"][end]
        second = transformed["terminal_geometry"][end]
        assert second["distance_to_selected_interface_angstrom"] == pytest.approx(
            first["distance_to_selected_interface_angstrom"], abs=0.01)
        assert second["aim_toward_partner_degrees"] == pytest.approx(first["aim_toward_partner_degrees"], abs=0.05)
        assert second["terminal_window_mean_plddt"] == first["terminal_window_mean_plddt"]
        assert second["terminal_window_secondary_structure"] == first["terminal_window_secondary_structure"]
        assert second["tag_clearance"] == first["tag_clearance"]
        # Shrake-Rupley is a finite point approximation, so its residual
        # rotational sensitivity is bounded rather than claimed to be zero.
        assert second["terminal_residue_sasa_angstrom2"] == pytest.approx(
            first["terminal_residue_sasa_angstrom2"], abs=1.0)
    assert baseline["sasa"]["sphere_points"] == evaluator.SASA_SPHERE_POINTS == 2400


def test_low_confidence_is_reported_as_possible_flexibility_and_missing_confidence_stays_missing(tmp_path, monkeypatch):
    from mn_molecular_design_studio import complex_analysis_tag_fusion as evaluator

    path = _write_complex(tmp_path / "low-confidence.pdb", low_confidence=True)
    monkeypatch.setattr(evaluator, "_tag_clearances", lambda *_args: {"his6": {"status": "clear"}})
    monkeypatch.setattr(evaluator, "_dssp_data", lambda _path, model: (
        {(str(chain.id), residue.id): {"secondary_structure": "-", "relative_accessibility": 0.5}
         for chain in model for residue in chain}, {"status": "measured", "method": "test fixture"}))
    result = _calculate(path)
    n_terminal = result["terminal_geometry"]["N"]
    assert n_terminal["terminal_window_mean_plddt"] == pytest.approx(52.0)
    assert n_terminal["flexibility_evidence"] == "loop_rich_low_declared_plddt_possible_flexibility"

    untrusted_b = _calculate(path, b_factor_interpretation="uninterpreted")
    assert untrusted_b["terminal_geometry"]["N"]["terminal_window_mean_plddt"] is None
    assert untrusted_b["terminal_geometry"]["N"]["terminal_window_mean_b_factor"] == pytest.approx(52.0)
    assert "unavailable unless separately supplied" in untrusted_b["confidence_interpretation"]


def test_registered_workflows_and_database_adapter_preserve_existing_structure_and_ingest_evaluation(tmp_path, monkeypatch):
    from mn_molecular_design_studio import job_catalog, molecule_imports
    from mn_molecular_design_studio.complex_analysis_sources import (
        apply_database_analysis_source, build_database_analysis_source,
    )
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import WORKFLOWS, workflow_schema
    from mn_molecular_design_studio.engine.workflow_engine.workflows.analysis import run_tag_fusion_suitability_contract
    from mn_molecular_design_studio.evaluations import evaluation_descriptor, file_digest
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    from mn_molecular_design_studio.molecular_database_ingest import ingest_run
    from mn_molecular_design_studio.workspace import mark_studio_job

    workspace = tmp_path / "workspace"
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(workspace))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(workspace))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(workspace / "runs"))
    database = MolecularDatabase(workspace / "database")
    source_path = _write_complex(tmp_path / "database-complex.pdb")
    source_bytes = source_path.read_bytes()
    source_record = molecule_imports.parse_structure(source_bytes, source_path.name)[0]
    imported = database.import_records(USER, [source_record], name="Terminal geometry fixture")
    source = database.list_molecules(USER, molecule_ids=imported["molecule_ids"], collection="all")["items"][0]
    structure = source["representations"][0]
    model_id = str((structure.get("context") or {}).get("model_id") or structure["id"])
    selection = {
        "molecule_id": source["id"], "structure_id": structure["id"],
        "group_id": imported["group_id"], "model_id": model_id,
        "binder_chains": ["A"], "target_chains": ["B"],
    }
    assert {"analysis.tag_fusion_suitability", "analysis.tag_fusion_context"} <= set(WORKFLOWS)
    schema_names = {row["name"] for row in workflow_schema("analysis.tag_fusion_suitability")["parameters"]}
    assert {"component_format", "termini_confirmed", "b_factor_interpretation", "interface_contact_cutoff"} <= schema_names

    adapter = build_database_analysis_source(
        USER, selection, workflow_id="analysis.tag_fusion_suitability", project_id=None, database=database,
    )
    contextual_selection = {**selection, "source_component_chain": "A"}
    with pytest.raises(ValueError, match="explicit, current, unique exact sequence association"):
        build_database_analysis_source(
            USER, contextual_selection, workflow_id="analysis.tag_fusion_context", project_id=None, database=database,
        )
    with pytest.raises(ValueError, match="caller-supplied paths or identities"):
        apply_database_analysis_source({
            "schema_version": 1, "workflow": "analysis.tag_fusion_context",
            "database_source": contextual_selection, "parameters": {"reference_context": {"status": "validated"}},
        }, USER, None)
    with pytest.raises(PermissionError):
        build_database_analysis_source(
            OTHER_USER, selection, workflow_id="analysis.tag_fusion_suitability", project_id=None, database=database,
        )
    from mn_molecular_design_studio.complex_analysis_sources import apply_database_analysis_source
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import validate_workflow_request
    request = apply_database_analysis_source({
        "schema_version": 1, "workflow": "analysis.tag_fusion_suitability",
        "database_source": selection,
        "parameters": {"interface_contact_cutoff": 4.5, "termini_confirmed": True,
                       "component_format": "single_chain_linear_protein",
                       "b_factor_interpretation": "uninterpreted"},
    }, USER, None)
    request_path = tmp_path / "tag-fusion-workflow-request.json"
    request_path.write_text(json.dumps({key: request[key] for key in ("schema_version", "workflow", "parameters")}))
    validated = validate_workflow_request(request_path)
    assert validated["valid"] is True, validated
    original_digest = file_digest(database.structure_file(USER, structure["id"]))
    run_dir = run_tag_fusion_suitability_contract(
        Path(adapter["source_run_dir"]), Path(adapter["candidates_jsonl"]),
        candidate_ids=[adapter["candidate_id"]],
        selected_structure_path=Path(adapter["selected_structure_path"]),
        selected_model_id=adapter["selected_model_id"], binder_chains=["A"], target_chains=["B"],
        termini_confirmed=True,
    )
    payload = json.loads((run_dir / "artifacts/analysis/tag_fusion_suitability.json").read_text())
    assert payload["source"]["structure_id"] == structure["id"]
    assert payload["source"]["structure_digest"] == original_digest
    assert file_digest(database.structure_file(USER, structure["id"])) == original_digest
    assert payload["assessment_status"] == "measured"
    mark_studio_job({"run_dir": str(run_dir)}, project_id=None, owner_id=USER["id"], allow_private_source=True)
    assert ingest_run(run_dir, database=database, force=True) == 1
    # Artifact lookup is authorized through the catalog ACL record and must not
    # depend on the separate full-text Results index returning this Evaluation.
    monkeypatch.setattr(job_catalog, "search_jobs", lambda **_kwargs: {"items": [], "total": 0})
    reopened = database.list_molecules(USER, molecule_ids=[source["id"]], collection="all")["items"][0]
    evaluations = [row for row in reopened["observations"] if row.get("namespace") == "tag_fusion_suitability"]
    assert len(evaluations) == 1
    descriptor = evaluation_descriptor(evaluations[0])
    assert descriptor["evaluator"] == "tag_fusion_suitability"
    raw = evaluations[0]["context"]["candidate"]["raw_metadata"]
    assert raw["tag_fusion_suitability"]["source"]["structure_id"] == structure["id"]
    assert raw["measurement_provenance"]["tag_fusion_suitability"]["method_version"]
    terminal = raw["tag_fusion_suitability"]["terminal_geometry"]["N"]
    assert terminal["secondary_structure"] in {"-", "H", "B", "E", "G", "I", "T", "S"}
    assert all(code in {"-", "H", "B", "E", "G", "I", "T", "S"}
               for code in terminal["terminal_window_secondary_structure"])
    assert len(reopened["representations"]) == len(source["representations"])
    from mn_molecular_design_studio.complex_analysis_ui import _database_inspector_data, _evaluation_options
    from mn_molecular_design_studio.complex_analysis_ui import component_inventory
    reopened_structure = next(item for item in reopened["representations"] if item["id"] == structure["id"])
    reopened_path = database.structure_file(USER, structure["id"])
    reopened_components = component_inventory(reopened_path.read_bytes(), reopened_path.name)
    _options, evaluations_by_id = _evaluation_options(reopened)
    associated = list(evaluations_by_id.values())
    # Older candidate-metadata round-trips could turn short DSSP labels into
    # run-relative paths. The Inspector should reopen the authorized result
    # artifact and retain its original scientific values.
    for observation, evaluation in associated:
        if evaluation.get("evaluator") != "tag_fusion_suitability":
            continue
        metadata = observation["context"]["candidate"]["raw_metadata"]
        metadata["tag_fusion_suitability"]["terminal_geometry"]["N"]["secondary_structure"] = "/old/run/-"
        metadata["tag_fusion_suitability"]["terminal_geometry"]["N"]["terminal_window_secondary_structure"] = [
            "/old/run/-", "/old/run/H", "/old/run/E",
        ]
    inspector = _database_inspector_data(
        USER, reopened, reopened_structure, reopened_path, reopened_components, associated,
    )
    saved_tag_item = next(item for item in inspector["evaluations"] if item["evaluator"] == "tag_fusion_suitability")
    assert saved_tag_item["tag_fusion_result"]["source"]["structure_digest"] == original_digest
    assert saved_tag_item["tag_fusion_result"]["terminal_geometry"]["N"]["secondary_structure"] == "-"
    assert saved_tag_item["tag_fusion_result"]["terminal_geometry"]["N"]["terminal_window_secondary_structure"]
    assert all(not str(code).startswith("/")
               for code in saved_tag_item["tag_fusion_result"]["terminal_geometry"]["N"]["terminal_window_secondary_structure"])
    assert any(row["label"].startswith("N-terminus geometry") and row["selectors"] for row in inspector["evidence"])


def test_reference_context_requires_and_reuses_an_ingested_exact_placement_evaluation(tmp_path, monkeypatch):
    from mn_molecular_design_studio import molecule_imports
    from mn_molecular_design_studio.complex_analysis_sources import build_database_analysis_source
    from mn_molecular_design_studio.engine.workflow_engine.workflows.analysis import (
        run_reference_placement_contract, run_tag_fusion_context_contract,
    )
    from mn_molecular_design_studio.evaluations import file_digest
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    from mn_molecular_design_studio.molecular_database_ingest import ingest_run
    from mn_molecular_design_studio.workspace import mark_studio_job

    workspace = tmp_path / "workspace-context"
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(workspace))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(workspace))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(workspace / "runs"))
    database = MolecularDatabase(workspace / "database")
    source_path = _write_complex(tmp_path / "associated-source.pdb")
    ref_lines = ["HEADER    EXPLICITLY ASSOCIATED REFERENCE FIXTURE",
                 "SEQRES   1 R    7  ALA CYS ASP GLU PHE GLY HIS"]
    serial = 1
    rotation = np.asarray([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    translation = np.asarray([11.0, -7.0, 3.0])
    for line in source_path.read_text(encoding="utf-8").splitlines():
        if line.startswith("ATOM  ") and line[21:22] == "B":
            xyz = np.asarray([float(line[start:start + 8]) for start in (30, 38, 46)])
            moved = rotation @ xyz + translation
            ref_lines.append(line[:6] + f"{serial:5d}" + line[11:21] + "R" + line[22:30]
                             + "".join(f"{value:8.3f}" for value in moved) + line[54:])
            serial += 1
    ref_lines.extend(["TER", "END", ""])
    reference_path = tmp_path / "associated-reference.pdb"
    reference_path.write_text("\n".join(ref_lines), encoding="utf-8")
    source_record = molecule_imports.parse_structure(source_path.read_bytes(), source_path.name)[0]
    reference_record = molecule_imports.parse_structure(reference_path.read_bytes(), reference_path.name)[0]
    source_import = database.import_records(USER, [source_record], name="Context source")
    reference_import = database.import_records(USER, [reference_record], name="Explicit reference")
    source = database.list_molecules(USER, molecule_ids=source_import["molecule_ids"], collection="all")["items"][0]
    reference = database.list_molecules(USER, molecule_ids=reference_import["molecule_ids"], collection="all")["items"][0]
    source_structure, reference_structure = source["representations"][0], reference["representations"][0]
    model_id = str((source_structure.get("context") or {}).get("model_id") or source_structure["id"])
    association = database.save_complex_reference_context(
        USER, source_molecule_id=source["id"], source_structure_id=source_structure["id"],
        source_chain="B", reference_molecule_id=reference["id"],
        reference_structure_id=reference_structure["id"], reference_chain="R",
    )
    assert association["status"] == "validated"
    placement_selection = {
        "molecule_id": source["id"], "structure_id": source_structure["id"],
        "group_id": source_import["group_id"], "model_id": model_id,
        "source_component_chain": "B", "binder_chains": [], "target_chains": [],
    }
    placement_adapter = build_database_analysis_source(
        USER, placement_selection, workflow_id="analysis.reference_placement",
        project_id=None, database=database,
    )
    placement_run = run_reference_placement_contract(
        Path(placement_adapter["source_run_dir"]), Path(placement_adapter["candidates_jsonl"]),
        candidate_ids=[placement_adapter["candidate_id"]],
        selected_structure_path=Path(placement_adapter["selected_structure_path"]),
        selected_model_id=placement_adapter["selected_model_id"], source_component_chain="B",
    )
    mark_studio_job({"run_dir": str(placement_run)}, project_id=None,
                    owner_id=USER["id"], allow_private_source=True)
    assert ingest_run(placement_run, database=database, force=True) == 1

    contextual_selection = {
        **placement_selection, "source_component_chain": "B",
        "binder_chains": ["A"], "target_chains": ["B"],
    }
    contextual_adapter = build_database_analysis_source(
        USER, contextual_selection, workflow_id="analysis.tag_fusion_context",
        project_id=None, database=database,
    )
    staged_context = contextual_adapter["source_manifest"]["tag_fusion_context"]
    assert staged_context["placement_evaluation_id"]
    assert staged_context["placement_run_id"] == placement_run.name
    assert staged_context["reference_structure_id"] == reference_structure["id"]
    assert staged_context["reference_structure_digest"] == file_digest(database.structure_file(USER, reference_structure["id"]))
    placement_matrix = np.asarray(staged_context["placement_result"]["matrix_4x4"], dtype=float)
    assert not np.allclose(placement_matrix, np.eye(4))
    contextual_run = run_tag_fusion_context_contract(
        Path(contextual_adapter["source_run_dir"]), Path(contextual_adapter["candidates_jsonl"]),
        candidate_ids=[contextual_adapter["candidate_id"]],
        selected_structure_path=Path(contextual_adapter["selected_structure_path"]),
        selected_model_id=model_id, binder_chains=["A"], target_chains=["B"],
        source_component_chain="B", termini_confirmed=True,
    )
    result = json.loads((contextual_run / "artifacts/analysis/tag_fusion_suitability.json").read_text())
    assert result["assessment_scope"] == "validated_reference_context"
    assert result["reference_context"]["placement_evaluation_id"] == staged_context["placement_evaluation_id"]
    assert result["reference_context"]["reference_structure_id"] == reference_structure["id"]
    assert result["reference_context"]["reference_component_chain"] == "R"
    assert any("unmodeled glycans" in item and "biological assembly" in item
               for item in result["checks_not_performed"])
    assert result["reference_context"]["explicit_branch_entity_annotations"] == []
    ca_source = []
    for line in Path(contextual_adapter["selected_structure_path"]).read_text(encoding="utf-8").splitlines():
        if line.startswith("ATOM  ") and line[21:22] == "A" and line[12:16].strip() == "CA":
            ca_source.append(np.asarray([float(line[start:start + 8]) for start in (30, 38, 46)]))
    source_direction = ca_source[0] - np.mean(ca_source[1:7], axis=0)
    source_direction /= np.linalg.norm(source_direction)
    expected_reference_direction = placement_matrix[:3, :3] @ source_direction
    assert result["terminal_geometry"]["N"]["direction"] == pytest.approx(expected_reference_direction, abs=1e-4)
    with database.connection() as connection:
        connection.execute("UPDATE molecule_groups SET owner_id=? WHERE id=?",
                           (OTHER_USER["id"], reference_structure["group_id"]))
    with pytest.raises(ValueError, match="explicit, current, unique exact sequence association"):
        build_database_analysis_source(
            USER, contextual_selection, workflow_id="analysis.tag_fusion_context",
            project_id=None, database=database,
        )


def test_inspector_terminal_buttons_select_exact_molstar_residue(tmp_path):
    from streamlit.testing.v1 import AppTest

    script = '''\
from mn_molecular_design_studio.complex_analysis_ui import _render_inspector_tag_fusion
items = [{"evaluator": "tag_fusion_suitability", "evaluation": {"id": "eval-tag", "run_id": "run-tag"},
          "tag_fusion_terminal_selectors": {"N": "A:10^B", "C": "A:20^_"},
          "tag_fusion_result": {"method": "test", "method_version": "v1", "assessment_scope": "intrinsic_complex",
            "assessment_status": "measured", "terminal_geometry": {
              "N": {"assessment_status": "measured", "terminal_residue": {"chain_id": "A", "author_number": 10, "insertion_code": "B", "residue_name": "SER", "selector": "A:10^B"}, "tag_clearance": {}},
              "C": {"assessment_status": "measured", "terminal_residue": {"chain_id": "A", "author_number": 20, "insertion_code": "", "residue_name": "GLY", "selector": "A:20^_"}, "tag_clearance": {}}},
            "tag_comparison": {}, "source": {}, "parameters": {}, "method_parameters": {}, "interpretation_limitations": []}}]
_render_inspector_tag_fusion(items, "fixture")
'''
    app = AppTest.from_string(script, default_timeout=20).run()
    assert not app.exception
    app.button(key="fixture-tag-fusion-highlight-N").click().run()
    assert app.session_state["fixture-inspector-highlights"] == ["A:10^B"]
    app.button(key="fixture-tag-fusion-highlight-C").click().run()
    assert app.session_state["fixture-inspector-highlights"] == ["A:20^_"]


def test_database_workspace_launches_intrinsic_tag_fusion_with_selected_chain_provenance(tmp_path, monkeypatch):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import complex_analysis_ui, molecule_imports, molecular_database
    from mn_molecular_design_studio.molecular_database import MolecularDatabase

    workspace = tmp_path / "workspace-ui"
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(workspace))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(workspace))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(workspace / "runs"))
    database = MolecularDatabase(workspace / "database")
    path = _write_complex(tmp_path / "ui-complex.pdb")
    parsed = molecule_imports.parse_structure(path.read_bytes(), path.name)[0]
    imported = database.import_records(USER, [parsed], name="Tag/fusion UI fixture")
    record = database.list_molecules(USER, molecule_ids=imported["molecule_ids"], collection="all")["items"][0]
    structure = record["representations"][0]
    structure_path = database.structure_file(USER, structure["id"])
    components = complex_analysis_ui.component_inventory(structure_path.read_bytes(), structure_path.name)
    monkeypatch.setattr(molecular_database, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(complex_analysis_ui, "_show_structure", lambda *_args, **_kwargs: None)
    source_key = hashlib.sha256(f"database:{record['id']}:{structure['id']}".encode()).hexdigest()[:14]
    script = f'''\
import streamlit as st
from mn_molecular_design_studio.complex_analysis_ui import _render_database_analysis
USER = {USER!r}
RECORD = {record!r}
STRUCTURE = {structure!r}
PATH = {str(structure_path)!r}
COMPONENTS = {components!r}
def schema_provider(workflow_id):
    return {{"workflow": workflow_id, "schema_version": 1, "parameters": []}}
def submit(schema, workflow, parameters, project, settings, title):
    st.session_state["submitted_workflow"] = workflow
    st.session_state["submitted_parameters"] = parameters
_render_database_analysis(USER, RECORD, STRUCTURE, __import__("pathlib").Path(PATH), COMPONENTS, "",
    workflow_schema=schema_provider, submit=submit)
'''
    app = AppTest.from_string(script, default_timeout=20).run()
    assert not app.exception
    run_key = f"complex-analysis-db-run-analysis-tag_fusion_suitability-{source_key}-tag-fusion"
    app.button(key=run_key).click().run()
    assert app.session_state["submitted_workflow"] == "analysis.tag_fusion_suitability"
    parameters = app.session_state["submitted_parameters"]
    assert parameters["database_source"]["molecule_id"] == record["id"]
    assert parameters["database_source"]["structure_id"] == structure["id"]
    assert parameters["database_source"]["binder_chains"] == ["A"]
    assert parameters["database_source"]["target_chains"] == ["B"]
    assert parameters["termini_confirmed"] is False
