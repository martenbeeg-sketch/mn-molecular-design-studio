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
