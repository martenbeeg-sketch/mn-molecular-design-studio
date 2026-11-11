from __future__ import annotations

from pathlib import Path

import pytest


def _write_simple_complex(path: Path, *, two_models: bool = False, include_hydrogens: bool = True,
                          binder_offset: float = 0.0, insertion_codes: dict | None = None,
                          preceding_binder_backbone: bool = False, add_ligand: bool = False) -> None:
    rows = [
        ("N", "SER", "A", 1, -2.0, 0.0, 0.0, "N"),
        ("CA", "SER", "A", 1, -1.0, 0.0, 0.0, "C"),
        ("C", "SER", "A", 1, -1.0, 1.0, 0.0, "C"),
        ("O", "SER", "A", 1, -1.0, 2.0, 0.0, "O"),
        ("CB", "SER", "A", 1, -0.5, 0.5, 0.0, "C"),
        ("OG", "SER", "A", 1, 0.0, 0.0, 0.0, "O"),
        ("H", "SER", "A", 1, 0.0, 0.0, 1.0, "H"),
        ("N", "ASP", "B", 1, 5.0, 0.0, 0.0, "N"),
        ("CA", "ASP", "B", 1, 4.0, 0.0, 0.0, "C"),
        ("C", "ASP", "B", 1, 4.0, 1.0, 0.0, "C"),
        ("O", "ASP", "B", 1, 4.0, 2.0, 0.0, "O"),
        ("CB", "ASP", "B", 1, 3.0, 0.0, 1.5, "C"),
        ("CG", "ASP", "B", 1, 1.0, 0.0, 2.8, "C"),
        ("OD1", "ASP", "B", 1, 0.0, 0.0, 2.8, "O"),
        ("OD2", "ASP", "B", 1, 1.0, 1.0, 3.5, "O"),
    ]
    lines = []
    models = (1, 2) if two_models else (1,)
    serial = 1
    insertion_codes = insertion_codes or {}
    for model_number in models:
        if two_models:
            lines.append(f"MODEL     {model_number:4d}\n")
        model_rows = list(rows)
        if preceding_binder_backbone:
            model_rows.extend([
                ("C", "GLY", "B", 0, 3.5, 0.0, 0.0, "C"),
                ("O", "GLY", "B", 0, 4.5, 0.0, 0.0, "O"),
            ])
        model_rows = sorted(model_rows, key=lambda row: (row[2], row[3]))
        for name, resname, chain, number, x, y, z, element in model_rows:
            if element == "H" and not include_hydrogens:
                continue
            if preceding_binder_backbone and chain == "B" and number == 1 and name == "N":
                x = 2.0
            if chain == "B":
                x += binder_offset
            insertion_code = (str(insertion_codes.get((chain, number), ""))[:1] or " ")
            lines.append(
                f"ATOM  {serial:5d} {name:>4s} {resname:>3s} {chain}{number:4d}{insertion_code:1s}   "
                f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{80.00:6.2f}          {element:>2s}\n"
            )
            serial += 1
        if add_ligand:
            lines.append(f"HETATM{serial:5d} {'C1':>4s} LIG C{1:4d}    "
                         f"{10.0:8.3f}{10.0:8.3f}{10.0:8.3f}{1.00:6.2f}{20.00:6.2f}          {'C':>2s}\n")
            serial += 1
        if two_models:
            lines.append("ENDMDL\n")
    lines.append("END\n")
    path.write_text("".join(lines), encoding="ascii")


def test_interface_chemistry_records_surface_hbonds_and_author_residue_evidence(tmp_path):
    from mn_molecular_design_studio.complex_analysis_inspector import structure_digest
    from mn_molecular_design_studio.complex_analysis_interface_chemistry import assess_protein_interface

    structure = tmp_path / "interface.pdb"
    _write_simple_complex(structure)
    digest = structure_digest(structure)
    result = assess_protein_interface(structure, binder_chains=["B"], target_chains=["A"], source_identity={
        "molecule_id": "mol-stage", "structure_id": "structure-stage", "group_id": "group-stage",
        "structure_digest": digest, "model_id": "model-stage", "source_candidate_id": "candidate-stage",
    })

    metrics = result["metrics"]
    assert metrics["interface_total_delta_sasa_angstrom2"] > 0
    assert metrics["interface_buried_surface_area_half_delta_sasa_angstrom2"] == pytest.approx(
        metrics["interface_total_delta_sasa_angstrom2"] / 2, abs=0.001)
    assert metrics["interface_total_delta_sasa_angstrom2"] == pytest.approx(
        metrics["binder_delta_sasa_angstrom2"] + metrics["target_delta_sasa_angstrom2"], abs=0.002)
    assert metrics["hydrophobic_cs_buried_fraction"] is not None
    # The explicit OG→OD1 bond and a separate geometry-only backbone N→OD1
    # contact both satisfy the declared method; keep the evidence bases distinct.
    assert metrics["interfacial_hydrogen_bond_count"] == 2
    assert metrics["interfacial_hydrogen_bond_evidence_basis_counts"]["explicit"] == 1
    assert any(row["donor"]["atom"] == "OG" and row["acceptor"]["atom"] == "OD1"
               and row["evidence_basis"] == "explicit" for row in result["hydrogen_bonds"])
    assert metrics["contacting_residue_pair_count"] > 0
    assert {row["chain_id"] for row in result["residue_annotations"]} == {"A", "B"}
    assert all(row["selector"] for row in result["residue_annotations"])
    assert result["source"]["structure_digest"] == digest
    assert result["method_version"] == "studio.ppi-interface-chemistry.v1"
    assert "not energies" in result["methodological_assumptions"][3]
    assert "distinct from BinderScope's Interface_Hydrophobicity" in result["methodological_assumptions"][1]


def test_interface_chemistry_rejects_digest_mismatch_and_multi_model_input(tmp_path):
    from mn_molecular_design_studio.complex_analysis_interface_chemistry import assess_protein_interface

    structure = tmp_path / "multi.pdb"
    _write_simple_complex(structure, two_models=True)
    with pytest.raises(ValueError, match="digest"):
        assess_protein_interface(structure, binder_chains=["B"], target_chains=["A"], source_identity={
            "structure_digest": "unrelated", "model_id": "model-stage",
        })
    from mn_molecular_design_studio.complex_analysis_inspector import structure_digest
    with pytest.raises(ValueError, match="single-model"):
        assess_protein_interface(structure, binder_chains=["B"], target_chains=["A"], source_identity={
            "structure_digest": structure_digest(structure), "model_id": "model-stage",
        })


def test_interface_chemistry_zero_contacts_are_not_reported_as_zero_buried_fraction(tmp_path):
    from mn_molecular_design_studio.complex_analysis_inspector import structure_digest
    from mn_molecular_design_studio.complex_analysis_interface_chemistry import assess_protein_interface

    structure = tmp_path / "separated.pdb"
    _write_simple_complex(structure, binder_offset=100.0)
    result = assess_protein_interface(structure, binder_chains=["B"], target_chains=["A"], source_identity={
        "structure_digest": structure_digest(structure), "model_id": "separated-model",
    })
    assert result["metrics"]["contacting_residue_pair_count"] == 0
    assert result["metrics"]["target_interface_residue_count"] == 0
    assert result["metrics"]["buried_unpaired_polar_atom_count"] == 0
    assert result["metrics"]["hydrophobic_cs_buried_fraction"] is None
    assert result["contacts"] == []


def test_interface_chemistry_labels_missing_hydrogen_evidence_and_preserves_insertion_codes(tmp_path):
    from mn_molecular_design_studio.complex_analysis_inspector import structure_digest
    from mn_molecular_design_studio.complex_analysis_interface_chemistry import assess_protein_interface

    structure = tmp_path / "no-hydrogen-insertion.pdb"
    _write_simple_complex(structure, include_hydrogens=False, preceding_binder_backbone=True,
                          insertion_codes={("B", 1): "A"})
    result = assess_protein_interface(structure, binder_chains=["B"], target_chains=["A"], source_identity={
        "structure_digest": structure_digest(structure), "model_id": "insertion-model",
    })
    basis = result["metrics"]["interfacial_hydrogen_bond_evidence_basis_counts"]
    assert basis["explicit"] == 0
    assert basis["computed_backbone_amide"] > 0
    assert basis["heavy_atom_angle_proxy"] > 0
    assert any(row["chain_id"] == "B" and row["author_number"] == 1
               and row["insertion_code"] == "A" and row["selector"]
               for row in result["residue_annotations"])


def test_interface_chemistry_rejects_ligand_only_selected_partner(tmp_path):
    from mn_molecular_design_studio.complex_analysis_inspector import structure_digest
    from mn_molecular_design_studio.complex_analysis_interface_chemistry import assess_protein_interface

    structure = tmp_path / "protein-ligand.pdb"
    _write_simple_complex(structure, add_ligand=True)
    with pytest.raises(ValueError, match="no supported amino-acid residues"):
        assess_protein_interface(structure, binder_chains=["B"], target_chains=["C"], source_identity={
            "structure_digest": structure_digest(structure), "model_id": "ligand-model",
        })


def test_interface_chemistry_is_in_comprehensive_profile_and_reuses_its_versioned_evaluation():
    from mn_molecular_design_studio.complex_analysis_interface_chemistry import METHOD_VERSION
    from mn_molecular_design_studio.complex_analysis_batch import (
        ASSESSMENT_METHODS, COMPREHENSIVE_PROFILE, assessment_applicability, evaluation_coverage,
    )

    method = "analysis.binder_interface_chemistry"
    assert method in COMPREHENSIVE_PROFILE
    sample = {"kind": "protein-protein", "structure_id": "s1", "digest": "d1", "model_id": "m1",
              "components": [{"type": "protein", "source_chain": "A", "sequence": "A" * 40},
                             {"type": "protein", "source_chain": "B", "sequence": "G" * 40}]}
    assert assessment_applicability(sample, method, ["B"], ["A"]) == (True, "")
    assert evaluation_coverage(sample, method, parameters={}, binder_chains=["B"], target_chains=["A"])[0] == "Missing evaluation"
    observation = {"namespace": "binder_interface_chemistry_v1", "context": {
        "representation_id": "s1", "source_structure_digest": "d1", "settings": {"method_version": METHOD_VERSION},
        "candidate": {"candidate_id": "db-candidate", "model_id": "m1", "binder_chains": ["B"], "target_chains": ["A"]},
    }}
    descriptor = {"id": "evaluation-stage", "evaluator": "binder_interface_chemistry_v1", "representation_id": "s1"}
    sample["associated"] = [(observation, descriptor)]
    assert evaluation_coverage(sample, method, parameters={}, binder_chains=["B"], target_chains=["A"])[0] == "Already evaluated · compatible settings"
    assert ASSESSMENT_METHODS[method]["evaluators"] == {"binder_interface_chemistry_v1"}


def test_registered_interface_chemistry_contract_writes_result_candidate_and_exact_provenance(tmp_path, monkeypatch):
    import json
    from types import SimpleNamespace

    from mn_molecular_design_studio.complex_analysis_inspector import structure_digest
    from mn_molecular_design_studio.engine.workflow_engine.workflows import analysis

    source = tmp_path / "authorized-source"
    source.mkdir()
    structure = source / "selected.pdb"
    _write_simple_complex(structure)
    digest = structure_digest(structure)
    identity = {"molecule_id": "mol-stage", "structure_id": "structure-stage", "group_id": "group-stage",
                "project_id": "project-stage", "model_id": "model-stage", "structure_digest": digest,
                "candidate_id": "database-stage"}
    manifest = {"schema": "studio.complex_analysis_database_source.v1", **identity,
                "binder_chains": ["B"], "target_chains": ["A"]}
    (source / "complex_analysis_source.json").write_text(json.dumps(manifest), encoding="utf-8")
    (source / "metadata.json").write_text(json.dumps({"source_kind": "database", "run_id": "source-stage"}), encoding="utf-8")
    candidate = {"candidate_id": "database-stage", "source_molecule_id": "mol-stage", "source_tool": "molecular_database",
                 "stage": "database_structure", "complex_pdb": str(structure), "model_index": "model-stage",
                 "binder_chains": ["B"], "target_chains": ["A"], "metrics": {},
                 "raw_metadata": {"complex_analysis_selection": {"structure_path": str(structure),
                     "structure_digest": digest, "model_id": "model-stage", "candidate_id": "database-stage",
                     "binder_chains": ["B"], "target_chains": ["A"]}}}
    candidates = source / "analysis_subjects.jsonl"
    candidates.write_text(json.dumps(candidate) + "\n", encoding="utf-8")
    run = tmp_path / "analysis-run"
    run.mkdir()
    monkeypatch.setattr(analysis, "_analysis_job", lambda *args, **kwargs: SimpleNamespace(run_dir=run))
    monkeypatch.setattr(analysis, "finish_job", lambda *args, **kwargs: None)

    result_dir = analysis.run_binder_interface_chemistry_contract(
        source, candidates, candidate_ids=["database-stage"], selected_structure_path=structure,
        selected_model_id="model-stage", binder_chains=["B"], target_chains=["A"],
    )

    result = json.loads((result_dir / "artifacts/analysis/binder_interface_chemistry.json").read_text())
    output = json.loads((result_dir / "artifacts/normalized_candidates/candidates.jsonl").read_text().splitlines()[0])
    assert result["source"]["molecule_id"] == identity["molecule_id"]
    assert result["source"]["structure_digest"] == digest
    assert result["source"]["model_id"] == "model-stage"
    assert output["metrics"]["binder_interface_chemistry_total_delta_sasa_angstrom2"] == result["metrics"]["interface_total_delta_sasa_angstrom2"]
    assert output["parents"] == ["database-stage"]
    assert output["raw_metadata"]["binder_interface_chemistry"]["residue_annotations"]
    assert output["raw_metadata"]["measurement_provenance"]["binder_interface_chemistry"]["method_version"] == "studio.ppi-interface-chemistry.v1"


def test_inspector_lists_the_persisted_interface_chemistry_artifact_once(tmp_path, monkeypatch):
    import json

    from mn_molecular_design_studio import complex_analysis_ui

    run = tmp_path / "analysis-run"
    artifact_dir = run / "artifacts/analysis"
    artifact_dir.mkdir(parents=True)
    (artifact_dir / "binder_interface_chemistry.json").write_text(json.dumps({
        "residue_annotations": [{"chain_id": "B", "author_number": 7}],
    }), encoding="utf-8")
    monkeypatch.setattr(complex_analysis_ui, "_authorized_evaluation_run_dir", lambda _user, _observation: run)
    observation = {"context": {"representation_id": "structure-1", "source_structure_digest": "digest-1",
                                "candidate": {"candidate_id": "candidate-1"}}}

    artifacts = complex_analysis_ui._saved_evaluation_artifacts(
        {"id": "user-1"}, observation, {"id": "structure-1"}, "digest-1",
    )

    selected = [item for item in artifacts if item.get("name") == "binder_interface_chemistry.json"]
    assert len(selected) == 1
    assert selected[0]["rows"] == [{"chain_id": "B", "author_number": 7}]
