from __future__ import annotations

from mn_molecular_design_studio.complex_analysis_catalog import authorized_structure_catalog


class FakeDatabase:
    def __init__(self):
        self.rows = [{"id": f"molecule-{index}", "name": f"Record {index}"} for index in range(5)]
        self.visible = {row["id"] for row in self.rows}
        self.calls = []

    def list_structure_records(self, user, *, collection, query, limit, offset=0, molecule_ids=None):
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
