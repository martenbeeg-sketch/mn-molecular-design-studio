from __future__ import annotations

import os

import pytest

from mn_molecular_design_studio import database_browser_cache
from mn_molecular_design_studio.database_browser import component_subjects
from mn_molecular_design_studio.database_browser_cache import BrowserCache
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from test_database_browser import USER, _import, _record, database


def test_unchanged_metadata_and_results_are_reused_without_hydration(database, monkeypatch):
    _import(database)
    calls = []
    original = database_browser_cache.build_browser_index

    def build(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(database_browser_cache, "build_browser_index", build)
    monkeypatch.setattr(database, "list_molecules", lambda *args, **kwargs: pytest.fail("Metadata cache must not hydrate results"))
    cache = BrowserCache(database.path)
    first = cache.get_index(database, USER, scope="all", include_related=True)
    assert cache.get_index(database, USER, scope="all", include_related=True) is first
    subject = component_subjects(first)[0]
    rows = cache.get_results(database, first, subject)
    assert cache.get_results(database, first, subject) is rows
    assert len(calls) == 1
    cache.close()


def test_peer_edits_invalidate_cache_without_a_local_capture_revision(database):
    imported = _import(database)
    observation = database.add_observation(imported["group_id"], imported["molecule_ids"][0], "tool", "run-1", {"score": 1})
    with database.connection() as connection:
        connection.execute("UPDATE molecule_groups SET origin=json_set(origin,'$._replica_node','peer') WHERE id=?", (imported["group_id"],))
    cache = BrowserCache(database.path)
    first = cache.get_index(database, USER, scope="all", include_related=True)
    with database.connection() as connection:
        revision = connection.execute("SELECT local_revision FROM molecular_database_sync_state").fetchone()[0]
        connection.execute("UPDATE molecule_observations SET metrics=? WHERE id=?", ('{"score":7}', observation))
    with database.connection() as connection:
        assert connection.execute("SELECT local_revision FROM molecular_database_sync_state").fetchone()[0] == revision
    assert cache.get_index(database, USER, scope="all", include_related=True) is not first
    cache.close()


def test_deletion_invalidates_cached_outputs(database):
    imported = _import(database)
    cache = BrowserCache(database.path)
    first = cache.get_index(database, USER, scope="all", include_related=True)
    subject = component_subjects(first)[0]
    assert cache.get_results(database, first, subject)
    with database.connection() as connection:
        connection.execute("DELETE FROM molecule_members WHERE group_id=?", (imported["group_id"],))
    second = cache.get_index(database, USER, scope="all", include_related=True)
    assert not second["molecules"]
    with pytest.raises(ValueError, match="no longer current"):
        cache.get_results(database, first, subject)
    cache.close()


def test_project_revocation_invalidates_cache_even_without_a_library_commit(database):
    member = {**USER, "_cluster_project_ids": ["project"]}
    _import(database, user=member, project_id="project")
    cache = BrowserCache(database.path)
    first = cache.get_index(database, member, scope="all", include_related=True)
    assert first["molecules"]
    revoked = cache.get_index(database, USER, scope="all", include_related=True)
    assert not revoked["molecules"] and revoked is not first
    cache.close()


def test_replaced_database_and_reopened_watch_do_not_reuse_old_snapshots(database, tmp_path):
    _import(database)
    cache = BrowserCache(database.path)
    first = cache.get_index(database, USER, scope="all", include_related=True)
    cache.close()
    reopened = cache.get_index(database, USER, scope="all", include_related=True)
    assert reopened is not first
    replacement = MolecularDatabase(tmp_path / "replacement")
    imported = _import(replacement, records=[_record("KLMN")])
    os.replace(replacement.path, database.path)
    changed = cache.get_index(database, USER, scope="all", include_related=True)
    assert set(changed["molecules"]) == set(imported["molecule_ids"])
    assert changed is not reopened
    cache.close()
