from __future__ import annotations

from mn_molecular_design_studio.complex_analysis_catalog import authorized_structure_catalog


def test_source_context_uses_native_job_evidence_and_source_qualified_index(tmp_path, monkeypatch):
    import pytest
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    from mn_molecular_design_studio.complex_analysis_catalog import collection_source_contexts
    from mn_molecular_design_studio.scientific_lineage_ingest import snapshot_from_candidates
    from mn_molecular_design_studio.scientific_lineage import put_snapshot

    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path / "app"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(tmp_path / "app"))
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    database = MolecularDatabase(tmp_path / "db")
    user = {"id": "owner", "role": "admin"}
    origin = {"run_id": "same-run-id", "task_group": "design", "tool": "bindcraft2",
              "source_node": "remote-a", "workstation": "PC004032"}
    imported = database.import_records(user, [{"name": "saved binder", "components": [
        {"type": "protein", "sequence": "ACDE", "chain": "B"}]}], origin=origin)
    group_id, molecule_id = imported["group_id"], imported["molecule_ids"][0]
    database.add_observation(group_id, molecule_id, "bindcraft2", "same-run-id", {},
        context={"target_id": "original-target", "settings": {"target_name": "Tau R3"}})
    database.add_observation(group_id, molecule_id, "biotite_interface", "analysis-run", {},
        context={"target_id": "unrelated-target", "settings": {"target_name": "Wrong target"}})
    groups = [group for group in database.groups(user, collection="all") if group["id"] == group_id]
    context = collection_source_contexts(database, user, groups)[group_id]
    assert [target["target_id"] for target in context["targets"]] == ["original-target"]
    for node, value in (("remote-a", 3665), ("remote-b", 999)):
        s = snapshot_from_candidates({"run_id": "same-run-id", "task_group": "design",
            "owner_id": "owner", "project_id": "", "tool": "bindcraft2", "workstation": "PC004032"}, [], node_id=node)
        s["counts"] = [{"name": "trajectory_attempt_count", "value": value, "unit": "trajectories"}]
        put_snapshot(database, s, author=node)
    context = collection_source_contexts(database, user, groups)[group_id]
    assert context["job"]["node_id"] == "remote-a"
    assert context["counts"][0]["value"] == 3665
    assert context["campaign_ref"] == ""  # One-job projection isn't a recorded shared parent.
    with pytest.raises(PermissionError):
        collection_source_contexts(database, {"id": "stranger", "role": "user"}, groups)


class FakeDatabase:
    def __init__(self):
        self.rows = [{"id": f"molecule-{index}", "name": f"Record {index}"} for index in range(5)]
        self.visible = {row["id"] for row in self.rows}
        self.calls = []

    def list_structure_records(self, user, *, collection, query, limit, offset=0, molecule_ids=None,
                               project_id="", group_id=""):
        self.calls.append({"user": user["id"], "query": query, "limit": limit,
                           "offset": offset, "molecule_ids": molecule_ids})
        rows = [row for row in self.rows if query.casefold() in row["name"].casefold()]
        if molecule_ids is not None:
            rows = [row for row in rows if row["id"] in set(molecule_ids) and row["id"] in self.visible]
        else:
            rows = [row for row in rows if row["id"] in self.visible]
        items = rows[offset:offset + limit]
        return {"items": items,
                "next_offset": offset + len(items) if offset + len(items) < len(rows) else None}


def test_catalog_pages_broad_discovery_then_acl_revalidates_cached_ids():
    database = FakeDatabase()
    state = {}
    user = {"id": "user-1"}
    first = authorized_structure_catalog(database, user, session_state=state, cache_key="catalog", limit=2)
    assert [row["id"] for row in first["records"]] == ["molecule-0", "molecule-1"]
    assert first["next_offset"] == 2
    assert database.calls[-1]["molecule_ids"] is None
    assert database.calls[-1]["limit"] == 2

    # Revoke one record between Streamlit reruns. Revalidation uses the current
    # canonical ACL-filtered Database query, not cached display metadata.
    database.visible.remove("molecule-1")
    repeated = authorized_structure_catalog(database, user, session_state=state, cache_key="catalog", limit=2)
    assert [row["id"] for row in repeated["records"]] == ["molecule-0"]
    assert database.calls[-1]["molecule_ids"] == ["molecule-0", "molecule-1"]
    assert database.calls[-1]["limit"] == 2
    assert sum(call["molecule_ids"] is None for call in database.calls) == 1

    next_page = authorized_structure_catalog(database, user, session_state=state, cache_key="catalog", limit=2, advance=True)
    assert [row["id"] for row in next_page["records"]] == ["molecule-0", "molecule-2", "molecule-3", "molecule-4"]
    assert database.calls[-2]["offset"] == 0
    assert database.calls[-2]["limit"] == 4
    assert database.calls[-1]["molecule_ids"] == ["molecule-0", "molecule-2", "molecule-3", "molecule-4"]


def test_catalog_refresh_query_and_user_changes_start_independent_pages():
    database = FakeDatabase()
    state = {}
    first = authorized_structure_catalog(database, {"id": "user-1"}, session_state=state, cache_key="catalog", limit=2)
    assert first["discovered"]
    changed_query = authorized_structure_catalog(database, {"id": "user-1"}, session_state=state,
                                                 cache_key="catalog", query="Record 4", limit=2)
    assert changed_query["discovered"]
    assert [item["id"] for item in changed_query["records"]] == ["molecule-4"]
    changed_user = authorized_structure_catalog(database, {"id": "user-2"}, session_state=state,
                                                 cache_key="catalog", query="Record 4", limit=2)
    assert changed_user["discovered"]
    assert database.calls[-1]["user"] == "user-2"


def test_catalog_acl_revocation_of_every_loaded_row_does_not_skip_from_stale_offset():
    database = FakeDatabase()
    state = {}
    user = {"id": "user-1"}
    authorized_structure_catalog(database, user, session_state=state, cache_key="catalog", limit=2)
    database.visible.clear()
    revoked = authorized_structure_catalog(database, user, session_state=state, cache_key="catalog", limit=2)
    assert revoked["records"] == []
    next_page = authorized_structure_catalog(database, user, session_state=state,
                                             cache_key="catalog", limit=2, advance=True)
    assert next_page["records"] == []
    assert next_page["next_offset"] is None
    assert database.calls[-1]["offset"] == 0


def test_eligible_catalog_scans_past_unrelated_records_before_showing_empty_page():
    database = FakeDatabase()
    database.rows = [
        {"id": f"molecule-{index}", "name": f"Record {index}", "complex": index in {2, 4}}
        for index in range(6)
    ]
    database.visible = {row["id"] for row in database.rows}
    state = {}
    eligible = lambda row: row.get("complex") is True

    first = authorized_structure_catalog(
        database, {"id": "user-1"}, session_state=state, cache_key="eligible",
        limit=2, eligible=eligible,
    )

    assert [row["id"] for row in first["records"]] == ["molecule-2", "molecule-4"]
    assert first["authorized_count"] == 2
    # The page is full before the raw catalog is exhausted. The next click
    # continues at the first unscanned row, where the last molecule lives.
    assert first["next_offset"] == 5
    assert [call["offset"] for call in database.calls if call["molecule_ids"] is None] == [0, 2, 4]


def test_eligible_catalog_paginates_by_eligible_count_and_revalidates_acl():
    database = FakeDatabase()
    database.rows = [
        {"id": f"molecule-{index}", "name": f"Record {index}", "complex": index >= 1}
        for index in range(6)
    ]
    database.visible = {row["id"] for row in database.rows}
    state = {}
    eligible = lambda row: row.get("complex") is True
    user = {"id": "user-1"}

    first = authorized_structure_catalog(database, user, session_state=state, cache_key="eligible",
                                         limit=2, eligible=eligible)
    assert [row["id"] for row in first["records"]] == ["molecule-1", "molecule-2"]
    assert first["next_offset"] == 3

    repeated = authorized_structure_catalog(database, user, session_state=state, cache_key="eligible",
                                            limit=2, eligible=eligible)
    assert [row["id"] for row in repeated["records"]] == ["molecule-1", "molecule-2"]

    database.visible.remove("molecule-1")
    next_page = authorized_structure_catalog(database, user, session_state=state, cache_key="eligible",
                                             limit=2, eligible=eligible, advance=True)
    assert [row["id"] for row in next_page["records"]] == [
        "molecule-2", "molecule-3", "molecule-4", "molecule-5",
    ]
    assert next_page["next_offset"] is None
