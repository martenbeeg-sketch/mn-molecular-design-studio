"""AUD-005: local forms must not erase newer local or replicated rationale."""
from concurrent.futures import ThreadPoolExecutor
from threading import Barrier

import pytest
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio import molecular_database as db_module

USER = {"id": "a" * 32, "role": "user", "_cluster_project_ids": []}


@pytest.fixture
def reviews(tmp_path, monkeypatch):
    root = tmp_path / "studio"
    for key, value in {"HOME": root, "APP_HOME": root, "RUN_DIR": root / "runs"}.items():
        monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_" + key, str(value))
    database = db_module.MolecularDatabase(root / "database")
    identities, samples = [], []
    for index, sequence in enumerate(("ACDEFG", "ACDEFA")):
        imported = database.import_records(USER, [{"name": f"design {index}", "components": [
            {"type": "protein", "chain": "B", "sequence": sequence}]}], name=f"cohort {index}")
        molecule = imported["molecule_ids"][0]
        database.attach_structure(imported["group_id"], molecule, {"name": "model.cif", "content": b"data_fixture\n#\n"},
            context={"candidate_id": f"candidate-{index}", "model_id": "model-0"})
        record = database.list_molecules(USER, molecule_ids=[molecule])["items"][0]
        structure = record["representations"][0]
        identity = {"source_kind": "database", "subject_id": molecule, "structure_ref": structure["id"]}
        identities.append(identity)
        samples.append({"source_kind": "Database", "molecule_id": molecule, "structure_id": structure["id"],
                        "label": f"design {index}", "record": record, "structure": structure, "digest": "a" * 64})
    return database, identities, samples


def save(database, identity, note, status="Reviewed"):
    current = database.get_complex_review(USER, **identity)
    # This helper also exercises the pre-fix implementation for the UI reproductions.
    version = {"expected_revision": current["revision"]} if "revision" in current else {}
    return database.save_complex_review(USER, **identity, status=status, note=note, **version)


def test_stale_local_editor_cannot_overwrite_newer_rationale(reviews):
    database, (identity, _), _ = reviews
    original = save(database, identity, "Original")
    newer = save(database, identity, "Newer rationale")
    with pytest.raises(ValueError, match="changed|conflict"):
        database.save_complex_review(USER, **identity, note="Stale rationale", status="Reviewed",
                                     expected_revision=original["revision"])
    assert database.get_complex_review(USER, **identity) == newer


def test_two_simultaneous_local_editors_compare_and_set_atomically(reviews):
    database, (identity, _), _ = reviews
    base = save(database, identity, "Original")
    barrier = Barrier(2)
    def edit(note):
        db = db_module.MolecularDatabase(database.root)
        barrier.wait(timeout=5)
        try:
            return db.save_complex_review(USER, **identity, note=note, status="Reviewed", expected_revision=base["revision"])
        except ValueError:
            return None
    with ThreadPoolExecutor(2) as pool:
        outcomes = list(pool.map(edit, ["First new rationale", "Second new rationale"]))
    assert sum(value is not None for value in outcomes) == 1
    assert database.get_complex_review(USER, **identity)["note"] in {"First new rationale", "Second new rationale"}


@pytest.mark.parametrize("clear", ["status", "delete", "delete-recreate-delete"])
def test_edit_versus_clear_and_aba_stays_protected_after_restart(reviews, clear):
    database, (identity, _), _ = reviews
    absent = database.get_complex_review(USER, **identity)
    original = save(database, identity, "Original")
    if clear == "status":
        save(database, identity, "", "Unreviewed")
    else:
        with database.connection() as connection:
            connection.execute("DELETE FROM complex_review_annotations WHERE subject_id=?", (identity["subject_id"],))
        if clear == "delete-recreate-delete":
            save(database, identity, "Intervening edit")
            with database.connection() as connection:
                connection.execute("DELETE FROM complex_review_annotations WHERE subject_id=?", (identity["subject_id"],))
    restarted = db_module.MolecularDatabase(database.root)
    for stale in (original, absent):
        with pytest.raises(ValueError, match="changed|conflict"):
            restarted.save_complex_review(USER, **identity, status="Reviewed", note="Resurrected stale decision",
                                          expected_revision=stale["revision"])
    assert not restarted.get_complex_review(USER, **identity)["note"]


def test_unchanged_saves_and_blind_updates(reviews):
    database, (identity, _), _ = reviews
    saved = save(database, identity, "Stable")
    assert save(database, identity, "Stable") == saved
    with pytest.raises(ValueError, match="revision|changed|conflict"):
        database.save_complex_review(USER, **identity, note="Blind overwrite", status="Reviewed")
    assert database.get_complex_review(USER, **identity) == saved


def test_multi_design_save_is_atomic_when_one_member_is_stale(reviews):
    database, identities, _ = reviews
    bases = [save(database, identity, "Original") for identity in identities]
    save(database, identities[1], "Newer second design")
    with pytest.raises(ValueError, match="changed|conflict"):
        database.save_complex_reviews(USER, [{**identity, "status": "Reviewed", "note": "Batch stale",
            "expected_revision": base["revision"]} for identity, base in zip(identities, bases)])
    assert database.get_complex_review(USER, **identities[0])["note"] == "Original"
    assert database.get_complex_review(USER, **identities[1])["note"] == "Newer second design"


def test_existing_legacy_writer_invalidates_editor_without_changing_wire_contract(reviews):
    database, (identity, _), _ = reviews
    base = save(database, identity, "Original")
    with database.connection() as connection:
        # Old code / peer apply writes the existing annotation columns only.
        connection.execute("UPDATE complex_review_annotations SET note=? WHERE subject_id=?",
                           ("Legacy/peer update with unchanged timestamp", identity["subject_id"]))
    with pytest.raises(ValueError, match="changed|conflict"):
        database.save_complex_review(USER, **identity, status="Reviewed", note="Stale",
                                     expected_revision=base["revision"])


def app_for(reviews, monkeypatch):
    database, _, samples = reviews
    monkeypatch.setattr(db_module, "MolecularDatabase", lambda *a, **k: database)
    return AppTest.from_string(f'''
import streamlit as st
from mn_molecular_design_studio.complex_analysis_shortlist_export import render_shortlist_review_export
if "audit_samples" not in st.session_state:
    st.session_state.audit_samples = {samples[:1]!r}
render_shortlist_review_export({USER!r}, st.session_state.audit_samples, key_prefix="audit")
''', default_timeout=10).run()


def test_untouched_form_refreshes_and_exports_latest_persisted_review(reviews, monkeypatch):
    database, (identity, _), _ = reviews
    save(database, identity, "Original rationale")
    app = app_for(reviews, monkeypatch)
    assert not app.exception
    save(database, identity, "Newer rationale saved elsewhere")
    app.run()
    assert app.text_area[0].value == "Newer rationale saved elsewhere"
    assert app.session_state["audit_samples"][0]["review_annotation"]["note"] == "Newer rationale saved elsewhere"
    app.button(key="audit-save-review").click().run()
    assert not app.exception
    assert database.get_complex_review(USER, **identity)["note"] == "Newer rationale saved elsewhere"


def test_dirty_form_preserves_typing_reports_conflict_and_requires_explicit_reload(reviews, monkeypatch):
    database, (identity, _), _ = reviews
    save(database, identity, "Original rationale")
    app = app_for(reviews, monkeypatch)
    app.text_area[0].set_value("Unsaved local typing").run()
    save(database, identity, "Newer remote rationale")
    app.button(key="audit-save-review").click().run()
    assert not app.exception
    assert app.text_area[0].value == "Unsaved local typing"
    assert database.get_complex_review(USER, **identity)["note"] == "Newer remote rationale"
    assert any("changed" in message.value.casefold() for message in [*app.error, *app.warning])
    reload_button = next(b for b in app.button if b.label == "Load latest saved review")
    reload_button.click().run()
    assert app.text_area[0].value == "Newer remote rationale"


def test_revoked_access_blocks_save_and_does_not_export_cached_rationale(reviews, monkeypatch):
    database, (identity, _), _ = reviews
    save(database, identity, "Original rationale")
    app = app_for(reviews, monkeypatch)
    def denied(*args, **kwargs):
        raise PermissionError("Access revoked")
    monkeypatch.setattr(database, "get_complex_review", denied)
    app.run()
    assert not app.exception
    assert any("Access revoked" in warning.value for warning in app.warning)
    assert app.button(key="audit-save-review").disabled
    assert not app.get("download_button")


def test_user_switch_does_not_reuse_another_accounts_unsaved_editor(reviews, monkeypatch):
    database, (identity, _), samples = reviews
    save(database, identity, "First user's saved note")
    monkeypatch.setattr(db_module, "MolecularDatabase", lambda *a, **k: database)
    app = AppTest.from_string(f'''
import streamlit as st
from mn_molecular_design_studio.complex_analysis_shortlist_export import render_shortlist_review_export
user = st.session_state.get("active_user", {USER!r})
render_shortlist_review_export(user, {samples[:1]!r}, key_prefix="accounts")
''', default_timeout=10).run()
    app.text_area[0].set_value("First user's unsaved draft").run()
    app.session_state["active_user"] = {"id": "b" * 32, "role": "admin"}
    app.run()
    assert not app.exception
    assert app.text_area[0].value == ""
    assert database.get_complex_review(USER, **identity)["note"] == "First user's saved note"


def test_batch_permissions_are_rechecked_without_partial_changes(reviews):
    database, identities, _ = reviews
    bases = [save(database, identity, "Original") for identity in identities]
    with database.connection() as connection:
        connection.execute("UPDATE molecule_groups SET owner_id=? WHERE id IN "
                           "(SELECT group_id FROM molecule_structures WHERE id=?)",
                           ("b" * 32, identities[1]["structure_ref"]))
    with pytest.raises(PermissionError):
        database.save_complex_reviews(USER, [{**identity, "status": "Reviewed", "note": "Batch draft",
            "expected_revision": base["revision"]} for identity, base in zip(identities, bases)])
    assert database.get_complex_review(USER, **identities[0])["note"] == "Original"


def test_duplicate_batch_identity_does_not_write_twice(reviews):
    database, (identity, _), _ = reviews
    base = save(database, identity, "Original")
    decision = {**identity, "status": "Reviewed", "note": "New", "expected_revision": base["revision"]}
    with pytest.raises(ValueError, match="once"):
        database.save_complex_reviews(USER, [decision, decision])
    assert database.get_complex_review(USER, **identity) == base


def test_multi_design_ui_conflict_preserves_every_draft_and_saves_nothing(reviews, monkeypatch):
    database, identities, samples = reviews
    for identity in identities:
        save(database, identity, "Original")
    monkeypatch.setattr(db_module, "MolecularDatabase", lambda *a, **k: database)
    app = AppTest.from_string(f'''
from mn_molecular_design_studio.complex_analysis_shortlist_export import render_shortlist_review_export
render_shortlist_review_export({USER!r}, {samples!r}, key_prefix="multi")
''', default_timeout=10).run()
    app.text_area[0].set_value("First local draft").run()
    app.text_area[1].set_value("Second local draft").run()
    save(database, identities[1], "Second remote rationale")
    app.button(key="multi-save-review").click().run()
    assert not app.exception
    assert any("changed" in error.value for error in app.error)
    assert [field.value for field in app.text_area] == ["First local draft", "Second local draft"]
    assert [database.get_complex_review(USER, **identity)["note"] for identity in identities] == ["Original", "Second remote rationale"]


def test_inspector_uses_the_same_conflict_protected_editor(reviews, monkeypatch):
    database, (identity, _), _ = reviews
    save(database, identity, "Original")
    monkeypatch.setattr(db_module, "MolecularDatabase", lambda *a, **k: database)
    app = AppTest.from_string(f'''
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio.complex_analysis_ui import _show_human_review
_show_human_review(MolecularDatabase(), {USER!r}, **{identity!r})
''', default_timeout=10).run()
    app.text_area[0].set_value("Inspector local draft").run()
    save(database, identity, "Newer peer rationale")
    next(button for button in app.button if button.label == "Save review").click().run()
    assert not app.exception
    assert app.text_area[0].value == "Inspector local draft"
    assert any("changed" in error.value for error in app.error)
    assert database.get_complex_review(USER, **identity)["note"] == "Newer peer rationale"


def test_existing_review_survives_additive_version_metadata_without_rewrite(reviews):
    database, (identity, _), _ = reviews
    original = save(database, identity, "Historical human rationale")
    with database.connection() as connection:
        for action in ("insert", "update", "delete"):
            connection.execute(f"DROP TRIGGER complex_review_version_{action}")
        connection.execute("DROP TABLE complex_review_versions")
    reopened = db_module.MolecularDatabase(database.root)
    restored = reopened.get_complex_review(USER, **identity)
    for field in ("status", "note", "created_at", "updated_at"):
        assert restored[field] == original[field]
    saved = reopened.save_complex_review(USER, **identity, status="Reviewed", note="Explicit later review",
                                         expected_revision=restored["revision"])
    assert saved["created_at"] == original["created_at"]
    with reopened.connection() as connection:
        assert connection.execute("SELECT count(*) FROM complex_review_annotations").fetchone()[0] == 1


def test_export_serializer_receives_latest_persisted_rationale(reviews, monkeypatch):
    from mn_molecular_design_studio import complex_analysis_shortlist_export as export
    database, (identity, _), _ = reviews
    save(database, identity, "Original")
    manifests = []
    original = export.manifest_json
    def capture(manifest):
        manifests.append(manifest)
        return original(manifest)
    monkeypatch.setattr(export, "manifest_json", capture)
    app = app_for(reviews, monkeypatch)
    save(database, identity, "Newer persisted rationale")
    app.run()
    assert not app.exception
    assert manifests[-1]["records"][0]["human_review"]["rationale"] == "Newer persisted rationale"
