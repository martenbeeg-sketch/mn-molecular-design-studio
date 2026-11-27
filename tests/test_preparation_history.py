from collections import defaultdict
from io import BytesIO
import json
from pathlib import Path
from zipfile import ZipFile

import pytest

from mn_molecular_design_studio import molecule_imports
from mn_molecular_design_studio.legacy_preparation_history import recover_histories
from mn_molecular_design_studio.molecular_preparation import cropped_record
from mn_molecular_design_studio.preparation_history import MANIFEST, extend_history, portable_bundle, sha256
from test_database_browser import USER, database
from test_molecular_database_policy import _pdb_efg_acd_complex


def test_bundle_preserves_history_after_import_and_rejects_a_different_structure(database):
    content = _pdb_efg_acd_complex().encode()
    saved = extend_history(None, None, sha256(content), [{"label": "MODELLER gap reconstruction"}])
    bundle = portable_bundle("model.pdb", content, saved)
    record, = molecule_imports.parse_archive(bundle)
    assert record["files"][0]["preparation_history"] == {**saved, "imported_from_bundle": True}
    derived = cropped_record(record, ["A"])["files"][0]
    assert [step["label"] for step in derived["preparation_history"]["steps"]] == [
        "MODELLER gap reconstruction", "Select chains"]
    assert derived["preparation_history"]["source_sha256"] == sha256(content)

    raw, = molecule_imports.parse_structure(content, "model.pdb")
    group_id = database.import_records(USER, [raw], name="Protein inputs")["group_id"]
    database.import_records(USER, [record], group_id=group_id)
    with database.connection() as connection:
        stored, = connection.execute("SELECT context FROM molecule_structures WHERE group_id=?", (group_id,)).fetchall()
    assert json.loads(stored["context"])["preparation_history"] == {**saved, "imported_from_bundle": True}

    altered = BytesIO()
    with ZipFile(altered, "w") as archive:
        archive.writestr("model.pdb", content + b"\n")
        archive.writestr(MANIFEST, json.dumps(saved))
    with pytest.raises(ValueError, match="checksum"):
        molecule_imports.parse_archive(altered.getvalue())


def test_archive_recovery_updates_only_the_exact_saved_structure(database, tmp_path):
    original = _pdb_efg_acd_complex().encode()
    database.import_records(USER, molecule_imports.parse_structure(original, "model.pdb"), name="Test dataset")
    with database.connection() as connection:
        before = connection.execute("SELECT * FROM molecule_structures").fetchone()
    runs = tmp_path / "legacy"
    job = runs / "terminal-repair" / "old-run"
    job.mkdir(parents=True)
    (job / "metadata.json").write_text(json.dumps({"status": "completed", "job_type": "terminal_repair",
                                                   "chain": "A", "extension_sequence": "GG"}))
    (job / "prepared.pdb").write_bytes(original)
    (job / "same-name-different-data.pdb").write_bytes(original.replace(b"  1.00", b"  0.99"))
    targets = {len(original): {sha256(original)}}
    history = recover_histories(targets, (runs,))[sha256(original)]
    assert history["steps"][0]["label"] == "MODELLER repair"
    assert [item["path"] for item in history["artifacts"]] == ["terminal-repair/old-run/prepared.pdb"]
    assert database.record_preparation_history(USER, before["id"], history)
    assert not database.record_preparation_history(USER, before["id"], history)
    with database.connection() as connection:
        after = connection.execute("SELECT * FROM molecule_structures").fetchone()
    assert [before[key] for key in ("id", "path", "molecule_id", "group_id")] == [
        after[key] for key in ("id", "path", "molecule_id", "group_id")]
    assert json.loads(after["context"])["preparation_history"] == history
    assert database.structure_file(USER, before["id"]).read_bytes() == original
