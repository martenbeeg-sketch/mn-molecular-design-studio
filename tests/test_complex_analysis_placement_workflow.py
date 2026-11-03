from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from mn_molecular_design_studio import molecule_imports
from mn_molecular_design_studio.molecular_database import MolecularDatabase


USER = {"id": "c" * 32, "role": "user", "_cluster_project_ids": []}
OTHER_USER = {"id": "d" * 32, "role": "user", "_cluster_project_ids": []}
_BASE = np.asarray([[0, 0, 0], [3, 0, 0], [3, 2, 0], [3, 2, 3], [5, 3, 2]], dtype=float)
_ROTATION = np.asarray([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
_TRANSLATION = np.asarray([10.0, -4.0, 2.0])


def _write_chain(path: Path, chains: list[tuple[str, str, np.ndarray]]) -> Path:
    three = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
             "G": "GLY", "H": "HIS", "S": "SER", "T": "THR"}
    lines, serial = [], 1
    for chain, sequence, coordinates in chains:
        for index, (amino_acid, (x, y, z)) in enumerate(zip(sequence, coordinates)):
            lines.append(
                f"ATOM  {serial:5d} {'CA':^4s} {three[amino_acid]:>3s} {chain}{index + 1:4d} "
                f"  {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{20.0:6.2f}          C  "
            )
            serial += 1
        lines.append("TER")
    path.write_text("\n".join([*lines, "END", ""]))
    return path


def _database_pair(tmp_path, monkeypatch):
    workspace = tmp_path / "workspace"
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(workspace))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(workspace))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(workspace / "runs"))
    database = MolecularDatabase(workspace / "database")
    source_path = _write_chain(
        tmp_path / "source-complex.pdb",
        [("A", "CDEFG", _BASE), ("B", "AG", np.asarray([[20, 0, 0], [21, 1, 1]], dtype=float))],
    )
    reference_coords = np.asarray([[20, 1, -2], *(_BASE @ _ROTATION + _TRANSLATION), [30, 2, 4]], dtype=float)
    reference_path = _write_chain(tmp_path / "reference.pdb", [("R", "ACDEFGH", reference_coords)])
    source_record = molecule_imports.parse_structure(source_path.read_bytes(), source_path.name)[0]
    reference_record = molecule_imports.parse_structure(reference_path.read_bytes(), reference_path.name)[0]
    source_import = database.import_records(USER, [source_record], name="Explicitly associated complex")
    reference_import = database.import_records(USER, [reference_record], name="Reference construct")
    source = database.list_molecules(USER, molecule_ids=source_import["molecule_ids"], collection="all")["items"][0]
    reference = database.list_molecules(USER, molecule_ids=reference_import["molecule_ids"], collection="all")["items"][0]
    source_structure = source["representations"][0]
    reference_structure = reference["representations"][0]
    saved = database.save_complex_reference_context(
        USER, source_molecule_id=source["id"], source_structure_id=source_structure["id"],
        source_chain="A", reference_molecule_id=reference["id"],
        reference_structure_id=reference_structure["id"], reference_chain="R",
    )
    selection = {
        "molecule_id": source["id"], "structure_id": source_structure["id"],
        "group_id": source_import["group_id"],
        "model_id": (source_structure.get("context") or {}).get("model_id") or source_structure["id"],
        "binder_chains": [], "target_chains": [], "source_component_chain": "A",
    }
    return database, source, source_structure, reference, reference_structure, saved, selection


def test_database_adapter_requires_authorized_current_explicit_association(tmp_path, monkeypatch):
    from mn_molecular_design_studio.complex_analysis_sources import (
        apply_database_analysis_source, build_database_analysis_source,
    )

    database, source, source_structure, reference, reference_structure, saved, selection = _database_pair(tmp_path, monkeypatch)
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import WORKFLOWS, workflow_schema
    assert "analysis.reference_placement" in WORKFLOWS
    schema_names = {item["name"] for item in workflow_schema("analysis.reference_placement")["parameters"]}
    assert {"source_component_chain", "reference_context_id", "reference_structure_digest", "method_version"} <= schema_names
    assert saved["status"] == "validated"
    adapter = build_database_analysis_source(
        USER, selection, workflow_id="analysis.reference_placement", project_id=None, database=database,
    )
    reference_context = adapter["source_manifest"]["reference_context"]
    assert reference_context["id"] == saved["id"]
    assert reference_context["mapping"]["mapping_status"] == "validated"
    assert reference_context["reference_molecule_id"] == reference["id"]
    assert Path(reference_context["structure_path"]).is_file()
    assert reference_context["reference_structure_digest"] == saved["reference_structure_digest"]
    assert adapter["source"]["structure_id"] == source_structure["id"]
    assert reference_structure["id"] == saved["reference_structure_id"]

    with pytest.raises(ValueError, match="caller-supplied paths or identities"):
        apply_database_analysis_source({
            "schema_version": 1, "workflow": "analysis.reference_placement",
            "database_source": selection,
            "parameters": {"source_component_chain": "spoof", "reference_structure_path": "/tmp/attacker.pdb"},
        }, USER, None)
    request = apply_database_analysis_source({
        "schema_version": 1, "workflow": "analysis.reference_placement",
        "database_source": selection, "parameters": {"source_component_chain": "spoof"},
    }, USER, None)
    assert request["parameters"]["source_component_chain"] == "A"
    assert request["parameters"]["reference_context_id"] == saved["id"]
    assert request["parameters"]["reference_structure_id"] == saved["reference_structure_id"]
    assert request["parameters"]["reference_structure_digest"] == saved["reference_structure_digest"]
    assert request["parameters"]["source_run_dir"] == adapter["source_run_dir"]
    assert request["parameters"]["source_component_chain"] == "A"
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import validate_workflow_request
    request_path = tmp_path / "registered-placement-request.json"
    request_path.write_text(json.dumps({key: request[key] for key in ("schema_version", "workflow", "parameters")}))
    validated_request = validate_workflow_request(request_path)
    assert validated_request["valid"] is True, validated_request
    assert validated_request["workflow"] == "analysis.reference_placement"
    assert validated_request["parameters"]["reference_context_id"] == saved["id"]

    with pytest.raises(PermissionError):
        build_database_analysis_source(
            OTHER_USER, selection, workflow_id="analysis.reference_placement", project_id=None, database=database,
        )

    with database.connection() as connection:
        connection.execute("UPDATE molecule_groups SET owner_id=? WHERE id=?", (OTHER_USER["id"], reference_structure["group_id"]))
    hidden_reference = database.get_complex_reference_context(
        USER, source_molecule_id=source["id"], source_structure_id=source_structure["id"], source_chain="A",
    )
    assert hidden_reference["status"] == "reference_unavailable"
    with pytest.raises(ValueError, match="explicit, current, unique exact"):
        build_database_analysis_source(
            USER, selection, workflow_id="analysis.reference_placement", project_id=None, database=database,
        )


def test_stale_or_unvalidated_association_cannot_be_staged_for_placement(tmp_path, monkeypatch):
    from mn_molecular_design_studio.complex_analysis_sources import build_database_analysis_source

    database, source, source_structure, _reference, reference_structure, saved, selection = _database_pair(tmp_path, monkeypatch)
    reference_file = database.structure_file(USER, reference_structure["id"])
    reference_file.write_bytes(reference_file.read_bytes() + b"REMARK changed after association\n")
    stale = database.get_complex_reference_context(
        USER, source_molecule_id=source["id"], source_structure_id=source_structure["id"], source_chain="A",
    )
    assert stale["status"] == "stale"
    with pytest.raises(ValueError, match="explicit, current, unique exact"):
        build_database_analysis_source(
            USER, selection, workflow_id="analysis.reference_placement", project_id=None, database=database,
        )


def test_registered_placement_contract_ingests_evaluation_on_existing_database_structure(tmp_path, monkeypatch):
    from mn_molecular_design_studio.complex_analysis_sources import build_database_analysis_source
    from mn_molecular_design_studio.engine.workflow_engine.workflows.analysis import run_reference_placement_contract
    from mn_molecular_design_studio.evaluations import evaluation_descriptor
    from mn_molecular_design_studio.molecular_database_ingest import ingest_run

    database, source, source_structure, reference, _reference_structure, saved, selection = _database_pair(tmp_path, monkeypatch)
    adapter = build_database_analysis_source(
        USER, selection, workflow_id="analysis.reference_placement", project_id=None, database=database,
    )
    run_dir = run_reference_placement_contract(
        Path(adapter["source_run_dir"]), Path(adapter["candidates_jsonl"]),
        candidate_ids=[adapter["candidate_id"]],
        selected_structure_path=Path(adapter["selected_structure_path"]),
        selected_model_id=adapter["selected_model_id"], source_component_chain="A",
    )
    result = json.loads((run_dir / "artifacts/analysis/reference_placement/reference_placement.json").read_text())
    assert result["aligned_pair_count"] == 5
    assert result["rmsd_angstrom"] == pytest.approx(0.0, abs=1e-5)
    assert result["source"]["structure_id"] == source_structure["id"]
    assert result["reference"]["structure_id"] == saved["reference_structure_id"]
    assert result["reference"]["structure_digest"] == saved["reference_structure_digest"]
    with pytest.raises(ValueError, match="does not match the authorized saved association"):
        run_reference_placement_contract(
            Path(adapter["source_run_dir"]), Path(adapter["candidates_jsonl"]),
            candidate_ids=[adapter["candidate_id"]],
            selected_structure_path=Path(adapter["selected_structure_path"]),
            selected_model_id=adapter["selected_model_id"], source_component_chain="A",
            reference_model_id="incorrect-reference-model",
        )
    from mn_molecular_design_studio.workspace import mark_studio_job
    mark_studio_job({"run_dir": str(run_dir)}, project_id=None, owner_id=USER["id"], allow_private_source=True)
    assert ingest_run(run_dir, database=database, force=True) == 1

    after = database.list_molecules(USER, molecule_ids=[source["id"]], collection="all")["items"][0]
    observations = [item for item in after["observations"] if item.get("namespace") == "reference_placement"]
    assert len(observations) == 1
    assert evaluation_descriptor(observations[0])["evaluator"] == "reference_placement"
    assert observations[0]["context"]["representation_id"] == source_structure["id"]
    assert observations[0]["context"]["source_structure_digest"] == result["source"]["structure_digest"]
    observation_candidate = observations[0]["context"]["candidate"]
    observation_raw = observation_candidate["raw_metadata"]
    assert observation_raw["reference_placement"]["matrix_4x4"] == result["matrix_4x4"]
    assert observation_raw["measurement_provenance"]["reference_placement"]["reference"]["structure_id"] == saved["reference_structure_id"]
    from mn_molecular_design_studio.evaluations import file_digest
    assert file_digest(database.structure_file(USER, source_structure["id"])) == result["source"]["structure_digest"]
    assert len(after["representations"]) == len(source["representations"])

    from mn_molecular_design_studio.complex_analysis_ui import _reference_placement_match
    reopened = database.get_complex_reference_context(
        USER, source_molecule_id=source["id"], source_structure_id=source_structure["id"], source_chain="A",
    )
    model_id = selection["model_id"]
    matched = _reference_placement_match(USER, source, source_structure, reopened, model_id)
    assert matched is not None
    assert matched["result"]["reference"]["molecule_id"] == reference["id"]
    assert matched["result"]["matrix_4x4"] == result["matrix_4x4"]


def test_placement_overlay_passes_reference_first_and_validated_transform_to_molstar(tmp_path, monkeypatch):
    from mn_molecular_design_studio import complex_analysis_ui
    from mn_molecular_design_studio.engine.workflow_engine.app.components import molstar_viewer

    reference = _write_chain(tmp_path / "reference.pdb", [("R", "ACDEFGH",
                             np.asarray([[index, index % 2, index % 3] for index in range(7)], dtype=float))])
    source = _write_chain(tmp_path / "source.pdb", [("A", "CDEFG", _BASE)])
    matrix = [[0, -1, 0, 10], [1, 0, 0, -4], [0, 0, 1, 2], [0, 0, 0, 1]]
    captured = []
    monkeypatch.setattr(molstar_viewer, "molstar_custom_component",
                        lambda structures, **kwargs: captured.append((structures, kwargs)))

    complex_analysis_ui._show_reference_placement_overlay(
        reference, source, {}, {"matrix_4x4": matrix}, key="overlay",
        highlighted_selections=["A:4"],
    )

    assert len(captured) == 1
    structures, options = captured[0]
    assert len(structures) == 2
    assert " R" in structures[0].pdb
    assert " A" in structures[1].pdb
    assert structures[0].transform_matrix is None
    assert structures[1].transform_matrix == matrix
    assert structures[1].highlighted_selections == ["A:4"]
    assert options["key"] == "overlay"


def test_database_workspace_launches_placement_from_saved_association(tmp_path, monkeypatch):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import complex_analysis_ui, molecular_database

    database, source, source_structure, reference, _reference_structure, saved, _selection = _database_pair(tmp_path, monkeypatch)
    monkeypatch.setattr(molecular_database, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(complex_analysis_ui, "_show_structure", lambda *args, **kwargs: None)
    script = f'''\
import streamlit as st
from mn_molecular_design_studio.complex_analysis_ui import render_complex_analysis_workspace
def schema_provider(workflow_id):
    return {{"workflow": workflow_id, "schema_version": 1, "parameters": []}}
def submit(schema, workflow, parameters, project, settings, title):
    st.session_state["submitted_workflow"] = workflow
    st.session_state["submitted_parameters"] = parameters
st.session_state["complex-analysis-source-kind"] = "Database"
st.session_state["complex-analysis-open"] = {{"kind":"database", "molecule_id":{source['id']!r}, "structure_id":{source_structure['id']!r}}}
render_complex_analysis_workspace({USER!r}, workflow_schema=schema_provider, submit=submit)
'''
    app = AppTest.from_string(script, default_timeout=20).run()
    assert not app.exception
    key = __import__("hashlib").sha256(f"database:{source['id']}:{source_structure['id']}".encode()).hexdigest()[:14]
    run_key = f"complex-analysis-db-run-analysis-reference_placement-{key}-reference-placement"
    app.button(key=run_key).click().run()

    assert not app.exception
    assert app.session_state["submitted_workflow"] == "analysis.reference_placement"
    request = app.session_state["submitted_parameters"]
    assert request["database_source"]["molecule_id"] == source["id"]
    assert request["database_source"]["structure_id"] == source_structure["id"]
    assert request["database_source"]["source_component_chain"] == "A"
    assert request["reference_context_id"] == saved["id"]
    assert request["reference_molecule_id"] == reference["id"]
