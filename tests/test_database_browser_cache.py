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
    subjects = cache.get_subjects(first, role="targets")
    assert cache.get_subjects(first, role="targets") is subjects
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
    with pytest.raises(ValueError, match="no longer current"):
        cache.get_subjects(first, role="targets")
    assert not cache.get_subjects(second, role="targets")
    cache.close()


def test_project_revocation_invalidates_cache_even_without_a_library_commit(database):
    member = {**USER, "_cluster_project_ids": ["project"]}
    _import(database, user=member, project_id="project")
    cache = BrowserCache(database.path)
    first = cache.get_index(database, member, scope="all", include_related=True)
    assert first["molecules"]
    assert cache.get_subjects(first)
    revoked = cache.get_index(database, USER, scope="all", include_related=True)
    assert not revoked["molecules"] and revoked is not first
    assert not cache.get_subjects(revoked)
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


def test_saved_folders_are_reused_after_other_folders_and_raw_cache_eviction(database, monkeypatch):
    from mn_molecular_design_studio import data_browser

    for sequence in ["ACD", "EFG", "HIK", "LMN", "PQR", "STV"]:
        _import(database, run_id=sequence, kind="curated", method="user_uploads", records=[_record(sequence)])
    cache = BrowserCache(database.path)
    index = cache.get_index(database, USER, scope="all", include_related=True)
    subjects = component_subjects(index)
    calls = []
    original = data_browser.saved_results

    def project(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(data_browser, "saved_results", project)
    first = cache.get_saved_results(database, index, subjects[0])
    assert first
    for subject in subjects[1:]:
        cache.get_saved_results(database, index, subject)
    assert len(cache.results) == 4 and len(calls) == len(subjects)
    assert cache.get_saved_results(database, index, subjects[0]) is first
    assert len(calls) == len(subjects)
    cache.close()


def test_saved_scope_cache_matches_full_subject_not_only_its_display_identifier(database):
    for sequence in ["ACD", "EFG"]:
        _import(database, run_id=sequence, kind="curated", method="user_uploads", records=[_record(sequence)])
    cache = BrowserCache(database.path)
    index = cache.get_index(database, USER, scope="all", include_related=True)
    first, second = [{**subject, "id": "same-navigation-label"} for subject in component_subjects(index)]
    first_rows = cache.get_saved_results(database, index, first)
    second_rows = cache.get_saved_results(database, index, second)
    assert {row["molecule_id"] for row in first_rows} == first["molecule_ids"]
    assert {row["molecule_id"] for row in second_rows} == second["molecule_ids"]
    assert cache.get_saved_results(database, index, first) is first_rows
    cache.close()


def test_saved_cache_is_bounded_by_folder_count_and_rows(database, monkeypatch):
    monkeypatch.setattr(database_browser_cache, "MAX_SAVED_SCOPES", 2)
    monkeypatch.setattr(database_browser_cache, "MAX_SAVED_ROWS", 2)
    for sequence in ["ACD", "EFG", "HIK"]:
        _import(database, run_id=sequence, kind="curated", method="user_uploads", records=[_record(sequence)])
    cache = BrowserCache(database.path)
    index = cache.get_index(database, USER, scope="all", include_related=True)
    for subject in component_subjects(index):
        cache.get_saved_results(database, index, subject)
    assert len(cache.saved) == 2 and sum(entry[1] for entry in cache.saved.values()) == 2
    oversized = {"id": "all", "molecule_ids": set(index["molecules"]), "group_ids": set(index["groups"])}
    assert len(cache.get_saved_results(database, index, oversized)) == 3
    assert len(cache.saved) == 2
    cache.close()


def test_saved_reference_fallback_is_a_separate_cache_contract(database):
    from test_scientific_browser_tree import target

    target(database)
    cache = BrowserCache(database.path)
    index = cache.get_index(database, USER, scope="all", include_related=True)
    subject = component_subjects(index)[0]
    assert not cache.get_saved_results(database, index, subject)
    reference = cache.get_saved_results(database, index, subject, fallback_references=True)
    assert reference and reference[0]["outcome_kind"] == "Reference molecule"
    assert cache.get_saved_results(database, index, subject, fallback_references=True) is reference
    assert not cache.get_saved_results(database, index, subject)
    cache.close()


@pytest.mark.parametrize("change", ["peer_edit", "deletion", "revocation", "reopen", "replacement"])
def test_saved_cache_does_not_survive_database_or_permission_changes(database, change, tmp_path):
    user = {**USER, "_cluster_project_ids": ["shared"]}
    imported = _import(database, user=user, project_id="shared", kind="curated", method="user_uploads")
    cache = BrowserCache(database.path)
    first = cache.get_index(database, user, scope="all", include_related=True)
    subject = component_subjects(first)[0]
    assert cache.get_saved_results(database, first, subject) and cache.saved
    if change == "revocation":
        user = {**USER, "_cluster_project_ids": []}
    elif change == "reopen":
        cache.close()
    elif change == "replacement":
        replacement = MolecularDatabase(tmp_path / "replacement")
        _import(replacement, user=user, project_id="shared", records=[_record("KLMN")])
        os.replace(replacement.path, database.path)
    else:
        with database.connection() as connection:
            if change == "peer_edit":
                connection.execute("UPDATE molecule_groups SET origin=json_set(origin,'$._replica_node','peer'),name='Changed' WHERE id=?", (imported["group_id"],))
            else:
                connection.execute("DELETE FROM molecule_members WHERE group_id=?", (imported["group_id"],))
    second = cache.get_index(database, user, scope="all", include_related=True)
    assert second is not first and not cache.saved
    with pytest.raises(ValueError, match="no longer current"):
        cache.get_saved_results(database, first, subject)
    if change in {"revocation", "deletion"}:
        assert not second["molecules"]
    cache.close()


def test_commit_during_projection_is_not_retained_as_a_stable_saved_snapshot(database, monkeypatch):
    from mn_molecular_design_studio import data_browser

    imported = _import(database, kind="curated", method="user_uploads")
    cache = BrowserCache(database.path)
    index = cache.get_index(database, USER, scope="all", include_related=True)
    subject = component_subjects(index)[0]
    original = data_browser.saved_results

    def project(*args, **kwargs):
        inventory = original(*args, **kwargs)
        with database.connection() as connection:
            connection.execute("UPDATE molecule_members SET name='Changed during read' WHERE group_id=?", (imported["group_id"],))
        return inventory

    monkeypatch.setattr(data_browser, "saved_results", project)
    assert cache.get_saved_results(database, index, subject)
    assert not cache.saved
    assert cache.get_index(database, USER, scope="all", include_related=True) is not index
    cache.close()
