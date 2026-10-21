from __future__ import annotations

import json
from pathlib import Path

import pytest

from mn_molecular_design_studio import molecule_imports
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio.molecular_database_ingest import ingest_run


A = {"id": "a" * 32, "role": "user", "_cluster_project_ids": ["project-1"]}
B = {"id": "b" * 32, "role": "user", "_cluster_project_ids": []}


@pytest.fixture
def database(tmp_path, monkeypatch):
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(tmp_path))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(tmp_path / "runs"))
    return MolecularDatabase(tmp_path / "database")


def _make_run(database, task, run_id, *, owner=A["id"], project="project-1", tool="boltz2",
              workflow="refolding.optimized_prediction", status="completed", success=True, candidates=(), metrics=None):
    from mn_molecular_design_studio.workspace import STUDIO_SOURCE_ID

    run_dir = database.root.parent / "runs" / task / run_id
    (run_dir / "artifacts/normalized_candidates").mkdir(parents=True, exist_ok=True)
    (run_dir / "molecular_design_studio.json").write_text(json.dumps({
        "application": STUDIO_SOURCE_ID, "studio_owner_id": owner, "project_id": project,
    }))
    (run_dir / "metadata.json").write_text(json.dumps({
        "studio_owner_id": owner, "project_id": project, "task_group": task,
        "tool": tool, "workflow": workflow, "status": status, "campaign_name": run_id,
    }))
    (run_dir / "input.json").write_text(json.dumps({"tool": tool, "workflow": workflow, "params": {}}))
    (run_dir / "result.json").write_text(json.dumps({"success": success, "metrics": metrics or {}}))
    if candidates:
        (run_dir / "artifacts/normalized_candidates/candidates.jsonl").write_text(
            "".join(json.dumps(candidate) + "\n" for candidate in candidates))
    from mn_molecular_design_studio.job_catalog import reconcile
    reconcile()
    return run_dir


def _pdb_acd():
    lines = []
    for serial, (residue, x) in enumerate((("ALA", 1.0), ("CYS", 2.0), ("ASP", 3.0)), start=1):
        lines.append(f"ATOM  {serial:5d}  CA  {residue} A{serial:4d}    {x:8.3f}{1.0:8.3f}{1.0:8.3f}{1.00:6.2f}{20.00:6.2f}           C")
    return "\n".join([*lines, "TER", "END", ""])


def _pdb_acd_with_generic_ligand():
    ligand = [
        "HETATM   10  C1  LIG B   1       8.000   1.000   1.000  1.00 20.00           C",
        "HETATM   11  C2  LIG B   1       9.000   1.000   1.000  1.00 20.00           C",
    ]
    return _pdb_acd().replace("TER\nEND", "TER\n" + "\n".join(ligand) + "\nEND")


def test_library_collection_is_curated_by_origin_and_accumulates_accessible_scores(database):
    imported = database.import_records(A, molecule_imports.parse_smiles("CCO ethanol"), name="Imported")
    molecule_id = imported["molecule_ids"][0]
    old_auto = database.import_records(A, molecule_imports.parse_smiles("CCO ethanol"), name="Old auto",
                                       project_id="project-1", import_method="tool_output", origin={"type": "job", "run_id": "old"})
    database.add_observation(old_auto["group_id"], molecule_id, "vina", "old-run", {"energy": -7.1})
    generated = database.import_records(A, molecule_imports.parse_smiles("CCN generated"), name="Generated",
                                        project_id="project-1", import_method="tool_output",
                                        origin={"type": "job", "collection_kind": "generated"})

    library = database.groups(A, collection="library")
    all_groups = database.groups(A, collection="all")
    assert {row["id"] for row in library} == {imported["group_id"], generated["group_id"]}
    assert {row["id"] for row in all_groups} >= {imported["group_id"], old_auto["group_id"], generated["group_id"]}
    row = database.list_molecules(A, group_id=imported["group_id"], collection="library")["items"][0]
    assert row["scores"] == {"vina.energy": -7.1}
    assert database.metric_fields(A, group_id=imported["group_id"], collection="library") == ["vina.energy"]


def test_failed_or_empty_runs_do_not_create_groups_and_backbones_stay_all_indexed(database):
    empty = _make_run(database, "empty", "no-output", workflow="design_generation", tool="rfdiffusion3", candidates=())
    assert ingest_run(empty, database=database) == 0
    failed = _make_run(database, "failed", "candidate", workflow="design_generation", tool="bindcraft2", status="failed", success=False,
                       candidates=[{"candidate_id": "failed-1", "stage": "sequence_design", "binder_sequence": "ACDE"}])
    assert ingest_run(failed, database=database) == 0

    backbone = _make_run(database, "design_generation", "backbone", workflow="design_generation", tool="rfdiffusion3",
                         candidates=[{"candidate_id": "bb-1", "stage": "generation.backbone", "binder_pdb": "backbone.pdb"}])
    (backbone / "backbone.pdb").write_text(_pdb_acd())
    assert ingest_run(backbone, database=database) == 1
    assert database.groups(A, collection="library") == []
    all_groups = database.groups(A, collection="all")
    assert len(all_groups) == 1
    assert all_groups[0]["origin"]["collection_kind"] == "intermediate"


def test_generated_complex_refolding_is_library_output(database):
    run = _make_run(database, "design", "boltzgen-final", workflow="design_generation", tool="boltzgen",
                    candidates=[{"candidate_id": "design-1", "stage": "complex_refolding", "binder_sequence": "ACDE",
                                 "metrics": {"confidence": 0.91}}])
    assert ingest_run(run, database=database) == 1
    library = database.groups(A, project_id="project-1", collection="library")
    assert len(library) == 1
    assert library[0]["origin"]["collection_kind"] == "generated"
    assert database.list_molecules(A, group_id=library[0]["id"], collection="library")["total"] == 1


def test_generated_complex_reindex_keeps_binder_identity_and_only_passing_designs(database):
    candidate = {
        "candidate_id": "bindcraft2_00001", "stage": "complex_refolding",
        "binder_sequence": "EFG", "binder_chains": ["B"],
        "metrics": {"native_pass_filters": True},
    }
    run = _make_run(database, "design", "bindcraft2-reindex", workflow="design_generation",
                    tool="bindcraft2", candidates=[candidate])

    assert ingest_run(run, database=database) == 1
    group = next(group for group in database.groups(A, collection="all") if group["origin"].get("run_id") == "bindcraft2-reindex")
    original_id = database.list_molecules(A, group_id=group["id"])["items"][0]["id"]

    # Recreate the earlier duplicate full-complex row produced when structure
    # parsing changed a candidate's identity from binder-only to target+binder.
    complex_text = "\n".join([
        "ATOM      1  CA  ALA A   1       1.000   1.000   1.000  1.00 20.00           C",
        "ATOM      2  CA  CYS A   2       2.000   1.000   1.000  1.00 20.00           C",
        "ATOM      3  CA  ASP A   3       3.000   1.000   1.000  1.00 20.00           C",
        "ATOM      4  CA  GLU B   1       8.000   1.000   1.000  1.00 20.00           C",
        "ATOM      5  CA  PHE B   2       9.000   1.000   1.000  1.00 20.00           C",
        "ATOM      6  CA  GLY B   3      10.000   1.000   1.000  1.00 20.00           C",
        "END", "",
    ])
    complex_path = run / "complex.pdb"
    complex_path.write_text(complex_text)
    full_components = molecule_imports.parse_structure(complex_text.encode(), "complex.pdb")[0]["components"]
    with database.connection() as connection:
        duplicate_id = database._put(connection, group["id"], {
            "name": "bindcraft2_00001", "components": full_components, "files": [],
        })
    database.attach_structure(group["id"], duplicate_id, {"path": str(complex_path), "label": "Complex Pdb"})
    assert duplicate_id != original_id

    passing = {**candidate, "complex_pdb": "complex.pdb"}
    rejected = {
        "candidate_id": "bindcraft2_00002", "stage": "complex_refolding",
        "binder_sequence": "VVV", "binder_chains": ["B"],
        "metrics": {"native_pass_filters": False},
    }
    (run / "artifacts/normalized_candidates/candidates.jsonl").write_text(
        "".join(json.dumps(row) + "\n" for row in (passing, rejected))
    )

    assert ingest_run(run, database=database, force=True) == 1
    rows = database.list_molecules(A, group_id=group["id"])["items"]
    assert len(rows) == 1
    assert rows[0]["id"] == original_id
    assert [(component["chain"], component["sequence"]) for component in rows[0]["components"]] == [("B", "EFG")]
    assert len(rows[0]["structures"]) == 1
    assert rows[0]["structures"][0]["context"]["portable"] is True


@pytest.mark.parametrize("complex_identity", [False, True])
def test_protein_qc_reuses_the_exact_frozen_database_identity(database, complex_identity):
    binder = {"type": "protein", "chain": "B", "role": "binder", "sequence": "ACD"}
    target = {"type": "protein", "chain": "A", "role": "target", "sequence": "ACD"}
    components = [target, binder] if complex_identity else [binder]
    imported = database.import_records(A, [{"name": "Source design", "components": components}], project_id="project-1")
    molecule = imported["molecule_ids"][0]
    run = _make_run(database, "analysis", "qc-identity", tool="sequence_qc", workflow="analysis.sequence_qc", candidates=[{
        "candidate_id": "qc-source", "stage": "analysis", "binder_sequence": "ACD", "binder_chains": ["B"], "target_chains": ["A"],
        "complex_pdb": "complex.pdb", "metrics": {"qc_status": "completed", "qc_sequence_length": 3},
        "raw_metadata": {"analysis_subject": "protein", "studio_molecule_id": molecule, "components": components},
    }])
    protein_lines = [line for line in _pdb_acd().splitlines() if line.startswith("ATOM")]
    (run / "complex.pdb").write_text("\n".join([*protein_lines, "TER", *[line[:21] + "B" + line[22:] for line in protein_lines], "TER", "END", ""]))
    assert ingest_run(run, database=database) == 1
    rows = database.list_molecules(A, collection="all")["items"]
    assert len(rows) == 1
    assert rows[0]["id"] == molecule
    assert rows[0]["scores"]["sequence_qc.qc_sequence_length"] == 3
    assert rows[0]["structures"]


def test_prediction_indexes_only_the_best_structure_per_input(database):
    run = _make_run(database, "optimized-prediction", "best-model", tool="optimized_prediction_boltz2",
                    workflow="refolding.optimized_prediction", candidates=[
        {"candidate_id": "model-a", "stage": "complex_refolding", "complex_pdb": "model-a.pdb",
         "raw_metadata": {"batch_item_id": "input-1", "compound_name": "Probe", "ligand_smiles": ["CCO"]},
         "metrics": {"ranking_score": .72}},
        {"candidate_id": "model-b", "stage": "complex_refolding", "complex_pdb": "model-b.pdb",
         "raw_metadata": {"batch_item_id": "input-1", "compound_name": "Probe", "ligand_smiles": ["CCO"]},
         "metrics": {"ranking_score": .91}},
    ])
    structure = _pdb_acd_with_generic_ligand()
    (run / "model-a.pdb").write_text(structure)
    (run / "model-b.pdb").write_text(structure)

    assert ingest_run(run, database=database) == 1
    row = database.list_molecules(A, collection="all")["items"][0]
    assert len(row["structures"]) == 1
    assert row["structures"][0]["context"]["candidate_id"] == "model-b"
    assert row["structures"][0]["context"]["database_structure_policy"] == "best_prediction"
    assert row["structures"][0]["context"]["portable"] is True
    assert Path(row["structures"][0]["path"]).is_relative_to(database.root / "molecular_files")
    assert len(row["observations"]) == 2

    # Reindexing after a changed ranking replaces the old attached representative.
    candidate_file = run / "artifacts/normalized_candidates/candidates.jsonl"
    candidates = [json.loads(line) for line in candidate_file.read_text().splitlines()]
    candidates[0]["metrics"]["ranking_score"] = .99
    candidates[1]["metrics"]["ranking_score"] = .4
    candidate_file.write_text("".join(json.dumps(candidate) + "\n" for candidate in candidates))
    assert ingest_run(run, database=database, force=True) == 1
    row = database.list_molecules(A, collection="all")["items"][0]
    assert len(row["structures"]) == 1
    assert row["structures"][0]["context"]["candidate_id"] == "model-a"


def test_docking_indexes_only_the_top_pose_for_each_compound(database):
    from mn_molecular_design_studio import ligand_compounds

    run = _make_run(database, "vina", "best-pose", tool="vina", workflow="ligand.docking")
    ligand_compounds.write_compound_set(
        run / "artifacts/input/compound_set.json", name="Inputs",
        compounds=[{"compound_id": "CMP1", "compound_name": "Ethanol", "smiles": "CCO"}],
        owner_id=A["id"], project_id="project-1", created_at="now",
    )
    poses = run / "artifacts/poses"
    poses.mkdir(parents=True)
    (poses / "pose-1.pdbqt").write_text("REMARK top-ranked pose\n")
    (poses / "pose-2.pdbqt").write_text("REMARK lower-ranked pose\n")
    (run / "artifacts/docking_scores.csv").write_text(
        "compound_id\tstatus\tpose_rank\tbest_score_kcal_mol\tpose_path\n"
        "CMP1\tcompleted\t1\t-7.1\tartifacts/poses/pose-1.pdbqt\n"
        "CMP1\tcompleted\t2\t-8.0\tartifacts/poses/pose-2.pdbqt\n"
    )

    assert ingest_run(run, database=database) == 1
    row = next(row for row in database.list_molecules(A, collection="all")["items"] if row["name"] == "Ethanol")
    assert len(row["structures"]) == 1
    assert row["structures"][0]["context"]["database_structure_policy"] == "best_docking_pose"
    assert row["structures"][0]["context"]["portable"] is True
    assert Path(row["structures"][0]["path"]).suffix == ".pdbqt"
    assert "top-ranked pose" in Path(row["structures"][0]["path"]).read_text()
    assert len(row["observations"]) == 2


def test_generator_contracts_ingest_common_fields_and_preserve_output_identity(database):
    definitions = [
        ("bindcraft2_generator", "generation.backbone_sequence", "ACD", {"pTM": .8, "binder_ca_radius_of_gyration": 10}),
        ("boltzgen_generator", "generation.backbone_sequence", "ACE", {"binder_ca_radius_of_gyration": 11}),
        ("rfdiffusion3_foundry", "generation.backbone", None, {"binder_ca_radius_of_gyration": 12}),
    ]
    for index, (tool, stage, sequence, metrics) in enumerate(definitions):
        candidate = {"candidate_id": f"design-{index}", "source_tool": tool, "stage": stage,
                     "binder_sequence": sequence, "binder_chains": ["A"], "metrics": metrics}
        if not sequence:
            candidate["binder_pdb"] = "backbone.pdb"
        run = _make_run(database, "design", f"generator-{index}", tool=tool, workflow="design_generation", candidates=[candidate])
        if not sequence:
            (run / "backbone.pdb").write_text(_pdb_acd())
        assert ingest_run(run, database=database) == 1
        assert ingest_run(run, database=database) == 0
    page = database.list_molecules(A, sort="common.geometry.binder.ca_rg", descending=False)
    assert page["total"] == 3
    assert [row["scores"]["common.geometry.binder.ca_rg"] for row in page["items"]] == [10, 11, 12]
    first = page["items"][0]
    assert first["scores"]["common.confidence.generator.complex.ptm"] == .8
    assert first["metadata"]["output_kind"] == "designed_structure"
    assert len(first["observations"]) == 1
    assert first["components"][0]["role"] == "binder"
    assert len(database.list_molecules(A, collection="library")["items"]) == 2
    assert page["items"][2]["metadata"]["output_kind"] == "backbone_proposal"
    assert "common.design.sequence_length" not in page["items"][2]["scores"]


def test_sequence_only_design_is_not_identified_as_target_structure(database):
    from mn_molecular_design_studio.molecular_database_ingest import _record
    target = database.root / "target.pdb"
    target.write_text(_pdb_acd())
    record = _record({"candidate_id": "new-design", "stage": "sequence_design", "binder_sequence": "KKK",
                      "binder_chains": ["Z"], "target_chains": ["A"], "target_pdb": str(target),
                      "source_tool": "esm_if1", "raw_metadata": {"source_tool": "rfdiffusion_generator"}}, database.root)
    assert record["components"][0]["sequence"] == "KKK"
    assert record["components"][0]["chain"] == "Z"
    assert record["components"][0]["role"] == "binder"
    assert record["metadata"]["source_tool"] == "esm_if1"
    assert record["files"] == []


def test_repeated_evaluations_reuse_original_group_and_accumulate_once(database):
    components = [{"chain": "A", "type": "protein", "sequence": "ACDE"}]
    imported = database.import_records(A, [{"name": "input", "components": components}], name="Original", project_id="project-1")
    molecule_id = imported["molecule_ids"][0]
    candidate = {"candidate_id": molecule_id, "stage": "complex_refolding", "binder_sequence": "ACDE",
                 "metrics": {"confidence": 0.8}, "raw_metadata": {"studio_molecule_id": molecule_id, "components": components}}
    run_one = _make_run(database, "optimized-prediction", "eval-1", candidates=[candidate])
    run_two = _make_run(database, "optimized-prediction", "eval-2", candidates=[candidate])

    assert ingest_run(run_one, database=database) == 1
    assert ingest_run(run_one, database=database) == 0
    assert ingest_run(run_two, database=database) == 1
    rows = database.list_molecules(A, group_id=imported["group_id"], collection="library")["items"]
    assert len(rows) == 1
    assert len(rows[0]["observations"]) == 2
    assert {row["run_id"] for row in rows[0]["observations"]} == {"eval-1", "eval-2"}
    assert {row["id"] for row in database.groups(A, collection="all")} == {imported["group_id"]}


@pytest.mark.parametrize("mapping_path", [
    "artifacts/raw/optimized_prediction/input_mapping.json",
    "artifacts/raw/af3_jax_cofolding/input_mapping.json",
])
def test_batched_prediction_inputs_restore_compound_identity_for_database(database, mapping_path):
    candidates = [
        {"candidate_id": "batch-1-model", "source_tool": "optimized_prediction_boltz2",
         "stage": "complex_refolding", "complex_pdb": "compound-1.pdb",
         "metrics": {"ranking_score": .8}, "raw_metadata": {"batch_item_id": "task-1"}},
        {"candidate_id": "batch-2-model", "source_tool": "optimized_prediction_boltz2",
         "stage": "complex_refolding", "complex_pdb": "compound-2.pdb",
         "metrics": {"ranking_score": .9}, "raw_metadata": {"batch_item_id": "task-2"}},
    ]
    run = _make_run(database, "refolding-validation", "compound-batch", candidates=candidates)
    for filename in ("compound-1.pdb", "compound-2.pdb"):
        (run / filename).write_text(_pdb_acd_with_generic_ligand())
    (run / "input.json").write_text(json.dumps({
        "tool": "optimized_prediction_boltz2",
        "workflow": "refolding.optimized_prediction",
        "inputs": {"tasks": [
            {"name": "Ethanol", "compound_id": "CHEM-101", "compound_name": "Ethanol",
             "ligand_smiles": [{"id": "B", "smiles": "CCO"}]},
            {"name": "Cyclopropane", "compound_id": "CHEM-102", "compound_name": "Cyclopropane",
             "ligand_smiles": [{"id": "B", "smiles": "C1CC1"}]},
        ]},
        "params": {"engine": "boltz2"},
    }))
    mapping_file = run / mapping_path
    mapping_file.parent.mkdir(parents=True, exist_ok=True)
    mapping_file.write_text(json.dumps({"items": [
        {"id": "task-1", "name": "Ethanol"},
        {"id": "task-2", "name": "Cyclopropane"},
    ]}))

    assert ingest_run(run, database=database) == 2
    assert ingest_run(run, database=database) == 0
    records = database.list_molecules(A, collection="all")["items"]
    assert {record["name"] for record in records} == {"Ethanol", "Cyclopropane"}
    ligand_smiles = {
        record["name"]: next(component["smiles"] for component in record["components"] if component["type"] == "small_molecule")
        for record in records
    }
    assert ligand_smiles == {"Ethanol": "CCO", "Cyclopropane": "C1CC1"}
    assert {record["metadata"]["compound_id"] for record in records} == {"CHEM-101", "CHEM-102"}


def test_save_run_selection_is_exact_repeatable_and_copies_run_context(database):
    from mn_molecular_design_studio.molecular_database_results import list_run_molecules, save_run_selection

    run = _make_run(database, "proteinmpnn", "results", tool="proteinmpnn", workflow="sequence_design",
                    candidates=[
                        {"candidate_id": "chosen", "stage": "sequence_design", "binder_sequence": "ACD",
                         "binder_pdb": "chosen.pdb", "metrics": {"rank": 1}},
                        {"candidate_id": "excluded", "stage": "sequence_design", "binder_sequence": "ACDE",
                         "metrics": {"rank": 2}},
                    ])
    (run / "chosen.pdb").write_text(_pdb_acd())
    assert ingest_run(run, database=database) == 2
    run_records = list_run_molecules(run, A, database=database)
    assert len(run_records) == 2
    chosen = next(row for row in run_records if row["name"] == "chosen")
    assert chosen["run_group_ids"] and chosen["run_groups"][0]["collection_kind"] == "generated"

    first = save_run_selection(run, A, [chosen["id"]], database=database, name="Saved", project_id="")
    again = save_run_selection(run, A, [chosen["id"]], database=database, group_id=first["group_id"], project_id="")
    second = save_run_selection(run, A, [chosen["id"]], database=database, name="Saved again", project_id="")
    assert first["molecule_ids"] == again["molecule_ids"] == second["molecule_ids"] == [chosen["id"]]
    assert first["saved_count"] == again["saved_count"] == second["saved_count"] == 1
    assert database.group(first["group_id"], A)["project_id"] == ""
    assert database.list_molecules(A, group_id=first["group_id"])["total"] == 1
    assert database.list_molecules(A, group_id=second["group_id"])["total"] == 1

    with database.connection() as connection:
        for destination in (first["group_id"], second["group_id"]):
            assert connection.execute("SELECT count(*) FROM molecule_observations WHERE group_id=?", (destination,)).fetchone()[0] == 1
            assert connection.execute("SELECT count(*) FROM molecule_structures WHERE group_id=?", (destination,)).fetchone()[0] == 1
            copied_sources = connection.execute("SELECT count(*) FROM molecule_sources WHERE group_id=? AND json_extract(provenance,'$.saved_selection')=1", (destination,)).fetchone()[0]
            assert copied_sources == 1
        members = connection.execute("SELECT molecule_id FROM molecule_members WHERE group_id=?", (first["group_id"],)).fetchall()
    assert [row[0] for row in members] == [chosen["id"]]


def test_run_result_acl_destination_schema_and_project_are_enforced(database):
    from mn_molecular_design_studio.molecular_database_results import list_run_molecules, save_run_selection

    run = _make_run(database, "proteinmpnn", "protected", tool="proteinmpnn", workflow="sequence_design",
                    candidates=[{"candidate_id": "selected", "stage": "sequence_design", "binder_sequence": "ACD",
                                 "metrics": {"rank": 1}}])
    row = list_run_molecules(run, A, database=database)[0]
    with pytest.raises(PermissionError):
        list_run_molecules(run, B, database=database)
    with pytest.raises(PermissionError):
        save_run_selection(run, B, [row["id"]], database=database, project_id="")

    constrained = database.create_group(A, "Constrained", project_id="project-1", import_method="user_uploads",
                                        origin={"type": "user", "collection_kind": "curated"})
    database.set_schema(A, constrained, {"fields": [{"name": "required", "type": "text", "required": True}]})
    with pytest.raises(ValueError, match="Required property is missing"):
        save_run_selection(run, A, [row["id"]], database=database, group_id=constrained, project_id="project-1")
    assert database.list_molecules(A, group_id=constrained)["total"] == 0
    with pytest.raises(PermissionError):
        save_run_selection(run, A, [row["id"]], database=database, project_id="project-2")


def test_run_results_api_lists_and_promotes_only_selected_ids(database, monkeypatch):
    from fastapi.testclient import TestClient
    from mn_molecular_design_studio.api import app, _current_user
    from mn_molecular_design_studio import molecular_database_api
    from mn_molecular_design_studio import runtime

    monkeypatch.setattr(molecular_database_api, "sync_workspace", lambda *args, **kwargs: [])
    # The application module may already have loaded a site-level runtime path
    # before this isolated test fixture set its environment variables.
    monkeypatch.setattr(runtime, "runs_root", lambda: database.root.parent / "runs")
    run = _make_run(database, "proteinmpnn", "api-results", tool="proteinmpnn", workflow="sequence_design",
                    candidates=[{"candidate_id": "first", "stage": "sequence_design", "binder_sequence": "ACD"},
                                {"candidate_id": "second", "stage": "sequence_design", "binder_sequence": "ACDE"}])
    from mn_molecular_design_studio.molecular_database_results import list_run_molecules
    assert len(list_run_molecules(run, A, database=database)) == 2
    app.dependency_overrides[_current_user] = lambda: A
    try:
        with TestClient(app) as client:
            response = client.get("/api/v1/molecules/jobs/proteinmpnn/api-results")
            assert response.status_code == 200, response.text
            selected = response.json()["items"][0]["id"]
            saved = client.post("/api/v1/molecules/jobs/proteinmpnn/api-results/selections", json={
                "molecule_ids": [selected], "name": "Only selected", "project_id": "",
            })
            assert saved.status_code == 200, saved.text
            assert saved.json()["molecule_ids"] == [selected]
            assert database.list_molecules(A, group_id=saved.json()["group_id"])["total"] == 1
            app.dependency_overrides[_current_user] = lambda: B
            denied = client.get("/api/v1/molecules/jobs/proteinmpnn/api-results")
            assert denied.status_code == 404
    finally:
        app.dependency_overrides.pop(_current_user, None)


def test_sync_throttle_preserves_last_indexing_warning(database, monkeypatch):
    from mn_molecular_design_studio import molecular_database_ingest as ingest

    monkeypatch.setattr(ingest, "migrate_legacy_groups", lambda *args: (_ for _ in ()).throw(RuntimeError("migration issue")))
    monkeypatch.setattr("mn_molecular_design_studio.auth.accessible_project_ids", lambda user: [])
    monkeypatch.setattr("mn_molecular_design_studio.job_catalog.search_jobs",
                        lambda **kwargs: {"items": [], "next_cursor": None})

    first = ingest.sync_workspace(A, force=True, database=database)
    assert first == [{"run_id": "legacy-migration", "error": "migration issue"}]
    # A normal UI rerun is throttled, but must continue showing the warning.
    assert ingest.sync_workspace(A, database=database) == first

    monkeypatch.setattr(ingest, "migrate_legacy_groups", lambda *args: None)
    assert ingest.sync_workspace(A, force=True, database=database) == []
    assert ingest.sync_workspace(A, database=database) == []
