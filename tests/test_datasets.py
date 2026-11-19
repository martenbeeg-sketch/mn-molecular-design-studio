import pytest

from mn_molecular_design_studio import datasets, workspace
from mn_molecular_design_studio.data_browser import ReadOnlyLibrary
from mn_molecular_design_studio.database_browser import build_browser_index, molecular_subjects, subject_results
from mn_molecular_design_studio.database_outcomes import scientific_outcomes
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from test_database_browser import OTHER, USER, database


ADMIN = {**USER, "role": "admin"}


def compounds():
    return [{"name": "Ethanol", "components": [{"chain": "L", "type": "small_molecule", "smiles": "CCO"}]}]


def projects():
    return [workspace.create_project(name, owner_id=USER["id"])["id"] for name in ("Docking screen", "Cofolding screen")]


def test_dataset_is_independent_and_referenced_by_multiple_projects_without_copying(database):
    first, second = projects()
    source = database.import_records(ADMIN, datasets.input_records(compounds()), name="Reusable library",
        origin=datasets.input_origin(ADMIN, [first, second]))
    group = database.group(source["group_id"], ADMIN)
    assert group["project_id"] == ""
    assert datasets.is_input_dataset(group)
    assert datasets.project_references(group) == sorted([first, second])
    with database.connection() as connection:
        assert connection.execute("SELECT count(*) FROM molecules").fetchone()[0] == 1
        assert [row[0] for row in connection.execute("SELECT project_id FROM dataset_project_links WHERE group_id=? ORDER BY project_id",
            (group["id"],))] == sorted([first, second])
    assert datasets.catalog(database, ADMIN, project_id=first)[0]["id"] == source["group_id"]
    assert datasets.catalog(database, ADMIN, project_id=second)[0]["id"] == source["group_id"]
    assert database.groups(OTHER) == []


def test_dataset_links_do_not_grant_project_members_read_or_edit_access(database, monkeypatch):
    project = projects()[0]
    saved = database.import_records(ADMIN, compounds(), name="Private reusable inputs", origin=datasets.input_origin(ADMIN, [project]))
    monkeypatch.setattr("mn_molecular_design_studio.auth.accessible_project_ids", lambda user: [project])
    assert datasets.catalog(database, OTHER, project_id=project) == []
    with pytest.raises(PermissionError):
        database.group(saved["group_id"], OTHER)
    with pytest.raises(PermissionError):
        datasets.save_project_references(database, OTHER, saved["group_id"], [], expected_revision="unused")


def test_reference_edits_are_revision_guarded_and_do_not_change_identity_or_members(database):
    first, second = projects()
    saved = database.import_records(ADMIN, compounds(), name="Inputs", origin=datasets.input_origin(ADMIN))
    group = database.group(saved["group_id"], ADMIN)
    with database.connection() as connection:
        before = [tuple(row) for row in connection.execute("SELECT * FROM molecule_members")]
    updated = datasets.save_project_references(database, ADMIN, group["id"], [first, second], expected_revision=datasets.source_revision(group))
    assert datasets.project_references(updated) == sorted([first, second])
    with pytest.raises(datasets.DatasetConflictError):
        datasets.save_project_references(database, ADMIN, group["id"], [], expected_revision=datasets.source_revision(group))
    with database.connection() as connection:
        assert before == [tuple(row) for row in connection.execute("SELECT * FROM molecule_members")]
        assert connection.execute("SELECT count(*) FROM molecules").fetchone()[0] == 1
    current = database.group(group["id"], ADMIN)
    datasets.save_project_references(database, ADMIN, group["id"], [], expected_revision=datasets.source_revision(current))
    with database.connection() as connection:
        assert not connection.execute("SELECT * FROM dataset_project_links").fetchall()


def test_reference_update_rejects_a_stale_older_writer_and_deleted_source(database):
    saved = datasets.create_dataset(database, ADMIN, "Independent")
    with database.connection() as connection:
        connection.execute("UPDATE molecule_groups SET origin=json_set(origin,'$.legacy_note','old application edit') WHERE id=?", (saved["id"],))
    with pytest.raises(datasets.DatasetConflictError):
        datasets.save_project_references(database, ADMIN, saved["id"], [], expected_revision=datasets.source_revision(saved))
    with database.connection() as connection:
        connection.execute("DELETE FROM molecule_groups WHERE id=?", (saved["id"],))
    with pytest.raises(PermissionError):
        datasets.save_project_references(database, ADMIN, saved["id"], [], expected_revision=datasets.source_revision(saved))


def test_private_import_needs_no_project_and_raw_design_metadata_is_not_admission(database):
    records = compounds()
    records[0]["metadata"] = {"source_tool": "rfd4_proteina", "workflow_role": "design_generation", "acceptance_status": "passed"}
    saved = database.import_records(ADMIN, datasets.input_records(records), name="Unreviewed supplied sequences",
        origin=datasets.input_origin(ADMIN))
    index = build_browser_index(database, ADMIN, include_related=True)
    results = subject_results(database, index, molecular_subjects(index)[0])
    assert scientific_outcomes(index, results) == []
    assert database.group(saved["group_id"], ADMIN)["project_id"] == ""
    assert datasets.catalog(database, ADMIN)[0]["molecule_count"] == 1


def test_input_members_added_to_legacy_dataset_do_not_become_design_outcomes(database):
    legacy = database.create_group(ADMIN, "Original reference library")
    records = compounds()
    records[0]["metadata"] = {"source_tool": "rfd4_proteina", "workflow_role": "design_generation", "acceptance_status": "passed"}
    database.import_records(ADMIN, datasets.input_records(records), group_id=legacy)
    index = build_browser_index(database, ADMIN)
    assert scientific_outcomes(index, subject_results(database, index, molecular_subjects(index)[0])) == []


def test_project_references_validate_access_existence_cleanup_and_replica_authority(database, monkeypatch):
    project = projects()[0]
    saved = datasets.create_dataset(database, ADMIN, "Library")
    with pytest.raises(ValueError, match="no longer exists"):
        datasets.save_project_references(database, ADMIN, saved["id"], ["not-a-project"], expected_revision=datasets.source_revision(saved))
    monkeypatch.setattr("mn_molecular_design_studio.auth.accessible_project_ids", lambda user: [])
    with pytest.raises(PermissionError):
        datasets.save_project_references(database, USER, saved["id"], [project], expected_revision=datasets.source_revision(saved))
    with database.connection() as connection:
        connection.execute("UPDATE molecule_groups SET origin=json_set(origin,'$._replica_node','other-node') WHERE id=?", (saved["id"],))
    with pytest.raises(PermissionError, match="source workstation"):
        datasets.save_project_references(database, ADMIN, saved["id"], [], expected_revision="unused")


def test_dataset_projection_survives_replication_restart_and_group_deletion(database, monkeypatch, tmp_path):
    from mn_molecular_design_studio.molecular_database_peer import apply_record
    project = projects()[0]
    source = datasets.create_dataset(database, ADMIN, "Shared identity", project_ids=[project])
    with database.connection() as connection:
        raw = dict(connection.execute("SELECT * FROM molecule_groups WHERE id=?", (source["id"],)).fetchone())
    replica = MolecularDatabase(tmp_path / "replica")
    key = "molecule_group:source-node:" + source["id"]
    apply_record("molecule_group", key, {"node_id": "source-node", "group": raw}, "source-node", database=replica)
    reopened = MolecularDatabase(replica.root)
    mirrored = reopened.groups(ADMIN)[0]
    assert datasets.project_references(mirrored) == [project]
    with reopened.connection() as connection:
        assert connection.execute("SELECT project_id FROM dataset_project_links WHERE group_id=?", (mirrored["id"],)).fetchone()[0] == project
    apply_record("molecule_group", key, None, "source-node", database=reopened)
    with reopened.connection() as connection:
        assert connection.execute("SELECT count(*) FROM dataset_project_links").fetchone()[0] == 0


def test_legacy_reader_and_schema_upgrade_preserve_members_and_labels(database):
    saved = database.import_records(ADMIN, compounds(), name="Existing input library")
    with database.connection() as connection:
        before = [tuple(row) for row in connection.execute("SELECT * FROM molecule_members")]
    reader = ReadOnlyLibrary(database.root)
    assert datasets.catalog(reader, ADMIN)[0]["id"] == saved["group_id"]
    with database.connection() as connection:
        assert before == [tuple(row) for row in connection.execute("SELECT * FROM molecule_members")]


def test_old_schema_reads_without_mutation_then_backfills_project_reference_index(database):
    import sqlite3
    project = projects()[0]
    saved = datasets.create_dataset(database, ADMIN, "Portable library", project_ids=[project])
    with sqlite3.connect(database.path) as connection:
        for action in ("insert", "update", "delete"):
            connection.execute("DROP TRIGGER dataset_projects_" + action)
        connection.execute("DROP TABLE dataset_project_links")
    reader = ReadOnlyLibrary(database.root)
    assert datasets.catalog(reader, ADMIN, project_id=project)[0]["id"] == saved["id"]
    with sqlite3.connect(database.path) as connection:
        assert not connection.execute("SELECT 1 FROM sqlite_master WHERE name='dataset_project_links'").fetchone()
    with database.connection() as connection:
        assert connection.execute("SELECT project_id FROM dataset_project_links WHERE group_id=?", (saved["id"],)).fetchone()[0] == project
    assert database.group(saved["id"], ADMIN)["origin"] == saved["origin"]


def test_cleanup_guard_protects_a_project_bound_dataset_referenced_elsewhere(database):
    first, second = projects()
    saved = database.import_records(ADMIN, compounds(), name="Old project-shared library", project_id=first)
    group = database.group(saved["group_id"], ADMIN)
    updated = datasets.save_project_references(database, ADMIN, group["id"], [second], expected_revision=datasets.source_revision(group))
    assert updated["project_id"] == first
    assert datasets.cleanup_references(updated, first) == [second]


def test_source_project_cleanup_review_is_blocked_for_referenced_dataset(database, tmp_path):
    from mn_molecular_design_studio.project_lifecycle import ProjectNode
    first, second = projects()
    saved = database.import_records(ADMIN, compounds(), name="Shared reference library", project_id=first)
    group = database.group(saved["group_id"], ADMIN)
    datasets.save_project_references(database, ADMIN, group["id"], [second], expected_revision=datasets.source_revision(group))
    node = ProjectNode(home=workspace._studio_home(), runs=tmp_path / "runs", database=database,
        node_id="local", job_reader=lambda: [])
    plan = node.inspect(ADMIN, first, mode="delete", include_collections=True)
    assert any("Reusable dataset" in reason for reason in plan["blocked"])
    assert database.groups(ADMIN)[0]["molecule_count"] == 1


def test_dataset_view_does_not_export_another_groups_scores_or_coordinates(database, tmp_path):
    from mn_molecular_design_studio.molecule_exports import export_molecules
    saved = database.import_records(ADMIN, compounds(), name="Pure inputs", origin=datasets.input_origin(ADMIN))
    result = database.import_records(ADMIN, compounds(), name="Docking results")
    molecule_id = saved["molecule_ids"][0]
    database.add_observation(result["group_id"], molecule_id, "vina", "calculation", {"affinity": -10})
    view = datasets.DatasetView(database, saved["group_id"])
    row = view.list_molecules(ADMIN)["items"][0]
    assert row["scores"] == {} and row["observations"] == []
    content, suffix, mime = export_molecules(view, ADMIN, [molecule_id], "CSV")
    assert b"affinity" not in content and suffix == "csv"
    assert [group["id"] for group in row["groups"]] == [saved["group_id"]]


def test_dataset_entry_labels_and_aliases_are_loaded_and_searched_from_canonical_sql(database):
    from mn_molecular_design_studio import entity_properties
    saved = database.import_records(ADMIN, compounds(), name="Reference entries", origin=datasets.input_origin(ADMIN))
    molecule_id = saved["molecule_ids"][0]
    snapshot = entity_properties.entry_snapshot(database, ADMIN, saved["group_id"], molecule_id)
    entity_properties.save_group_properties(database, ADMIN, snapshot, "Screening reference ethanol", {"notes": "Use as reference"},
        ["EtOH standard"], expected_revision=snapshot["revision"])
    view = datasets.DatasetView(database, saved["group_id"])
    record = view.list_molecules(ADMIN, query="EtOH standard")["items"][0]
    assert record["name"] == "Screening reference ethanol"
    assert record["metadata"]["notes"] == "Use as reference"
    assert record["id"] == molecule_id


def test_dataset_preview_does_not_hydrate_scientific_observations(database, monkeypatch):
    from contextlib import contextmanager
    saved = database.import_records(ADMIN, compounds(), name="Inputs", origin=datasets.input_origin(ADMIN))
    observed = []
    original = database.connection

    @contextmanager
    def traced():
        with original() as connection:
            connection.set_trace_callback(lambda statement: observed.append(statement))
            yield connection

    monkeypatch.setattr(database, "connection", traced)
    view = datasets.DatasetView(database, saved["group_id"])
    assert view.list_molecules(ADMIN)["total"] == 1
    assert not any(statement.upper().startswith("SELECT") and any(table in statement.lower() for table in
        ("molecule_observations", "molecule_evaluations")) for statement in observed)


def test_connected_results_follow_exact_chemistry_across_distinct_pose_identities(database):
    from mn_molecular_design_studio.data_browser import dataset_results_scope
    from mn_molecular_design_studio.scientific_browser_tree import scientific_scopes
    library = database.import_records(ADMIN, datasets.input_records(compounds()), name="Input compounds", origin=datasets.input_origin(ADMIN))
    pose = {**compounds()[0], "identity_context": {"pose": "vina-pose"}}
    prediction = database.import_records(ADMIN, [pose], name="Docking outcomes", import_method="tool_output",
        origin={"collection_kind": "generated", "task_group": "ligand-docking", "tool": "vina", "run_id": "vina-run"})
    unrelated = database.import_records(ADMIN, [{**pose, "components": [{"chain": "L", "type": "small_molecule", "smiles": "CCN"}]}],
        name="Same name, different chemistry", import_method="tool_output", origin={"collection_kind": "generated", "task_group": "ligand-docking", "run_id": "other-run"})
    assert library["molecule_ids"] != prediction["molecule_ids"]
    index = build_browser_index(database, ADMIN, collection="library", include_related=True)
    scopes = scientific_scopes(index, [])
    identifier = dataset_results_scope(scopes, index, library["group_id"], library["molecule_ids"][0])
    assert identifier is not None
    results = subject_results(database, index, scopes[identifier]["subject"])
    result_ids = {result["molecule_id"] for result in results}
    assert prediction["molecule_ids"][0] in result_ids
    assert unrelated["molecule_ids"][0] not in result_ids
    whole_dataset = dataset_results_scope(scopes, index, library["group_id"])
    results = subject_results(database, index, scopes[whole_dataset]["subject"])
    assert prediction["molecule_ids"][0] in {result["molecule_id"] for result in results}
    assert unrelated["molecule_ids"][0] not in {result["molecule_id"] for result in results}
