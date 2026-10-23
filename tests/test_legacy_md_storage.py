from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import pytest

from mn_molecular_design_studio.legacy_md_storage import (
    LegacyMDStorageError,
    legacy_md_run_inventory,
    stage_legacy_md_run,
)


def _stage(source: Path, destination_run: Path, *, workstation: str = "remote-md-01"):
    return stage_legacy_md_run(
        source,
        destination_run,
        source_run_id=source.name,
        source_workstation=workstation,
    )


def test_stage_copies_all_files_and_writes_verified_provenance_manifest(tmp_path):
    source = tmp_path / "source-run"
    (source / "artifacts" / "md").mkdir(parents=True)
    (source / "metadata.json").write_text('{"run": "source"}')
    (source / "artifacts" / "md" / "large.xtc").write_bytes(b"trajectory-bytes")
    # Small and extension-agnostic files must be included too.
    (source / "artifacts" / "md" / "tiny.bin").write_bytes(b"x")
    destination_run = tmp_path / "destination-run"

    inventory = legacy_md_run_inventory(source)
    assert inventory["total_size_bytes"] == sum(
        path.stat().st_size for path in source.rglob("*") if path.is_file()
    )
    result = _stage(source, destination_run)
    copied_root = destination_run / "artifacts" / "legacy" / source.name

    assert result["skipped"] is False
    assert (copied_root / "artifacts/md/large.xtc").read_bytes() == b"trajectory-bytes"
    assert (copied_root / "artifacts/md/tiny.bin").read_bytes() == b"x"
    source_stat = (source / "artifacts/md/large.xtc").stat()
    copied_stat = (copied_root / "artifacts/md/large.xtc").stat()
    assert source_stat.st_ino != copied_stat.st_ino
    assert not (copied_root / "artifacts/md/large.xtc").is_symlink()

    manifest = json.loads((copied_root / "legacy_copy_manifest.json").read_text())
    assert manifest["source_run_id"] == source.name
    assert manifest["source_run_path"] == str(source.resolve())
    assert manifest["source_workstation"] == "remote-md-01"
    assert manifest["files"] == [
        {
            "relative_path": "artifacts/md/large.xtc",
            "size_bytes": len(b"trajectory-bytes"),
            "sha256": hashlib.sha256(b"trajectory-bytes").hexdigest(),
        },
        {
            "relative_path": "artifacts/md/tiny.bin",
            "size_bytes": 1,
            "sha256": hashlib.sha256(b"x").hexdigest(),
        },
        {
            "relative_path": "metadata.json",
            "size_bytes": len('{"run": "source"}'),
            "sha256": hashlib.sha256(b'{"run": "source"}').hexdigest(),
        },
    ]


def test_identical_copy_is_idempotently_validated_and_existing_bytes_are_checked(tmp_path):
    source = tmp_path / "source-run"
    source.mkdir()
    (source / "output.dcd").write_bytes(b"coordinates")
    destination_run = tmp_path / "destination-run"

    first = _stage(source, destination_run)
    second = _stage(source, destination_run)
    assert first["skipped"] is False
    assert second["skipped"] is True

    destination = destination_run / "artifacts/legacy/source-run/output.dcd"
    destination.write_bytes(b"tampered")
    with pytest.raises(LegacyMDStorageError, match="failed byte verification"):
        _stage(source, destination_run)


def test_source_file_symlink_within_run_is_copied_as_bytes(tmp_path):
    source = tmp_path / "source-run"
    source.mkdir()
    (source / "actual.dat").write_bytes(b"safe")
    (source / "linked.dat").symlink_to("actual.dat")
    destination_run = tmp_path / "destination-run"

    _stage(source, destination_run)
    copied_root = destination_run / "artifacts/legacy/source-run"
    assert (copied_root / "linked.dat").read_bytes() == b"safe"
    assert not (copied_root / "linked.dat").is_symlink()


def test_external_symlink_fails_before_final_destination_is_created(tmp_path):
    source = tmp_path / "source-run"
    source.mkdir()
    external = tmp_path / "secret.dat"
    external.write_bytes(b"outside")
    (source / "escape.dat").symlink_to(external)
    destination_run = tmp_path / "destination-run"

    with pytest.raises(LegacyMDStorageError, match="Symlink escapes"):
        _stage(source, destination_run)
    assert not (destination_run / "artifacts/legacy/source-run").exists()


def test_broken_links_special_files_and_reserved_manifest_paths_fail(tmp_path):
    source = tmp_path / "source-run"
    source.mkdir()
    (source / "broken.dat").symlink_to("missing.dat")
    with pytest.raises((LegacyMDStorageError, OSError)):
        legacy_md_run_inventory(source)

    (source / "broken.dat").unlink()
    (source / "legacy_copy_manifest.json").write_text("reserved")
    with pytest.raises(LegacyMDStorageError, match="reserved manifest"):
        legacy_md_run_inventory(source)

    if hasattr(os, "mkfifo"):
        (source / "pipe").unlink(missing_ok=True)
        (source / "legacy_copy_manifest.json").unlink()
        os.mkfifo(source / "pipe")
        with pytest.raises(LegacyMDStorageError, match="non-regular file"):
            legacy_md_run_inventory(source)


def test_source_and_destination_overlap_and_unsafe_run_ids_are_rejected(tmp_path):
    source = tmp_path / "source-run"
    source.mkdir()
    (source / "data").write_text("x")
    with pytest.raises(LegacyMDStorageError, match="must not overlap"):
        stage_legacy_md_run(source, source / "nested", source_run_id="nested")
    with pytest.raises(LegacyMDStorageError, match="single safe path component"):
        stage_legacy_md_run(source, tmp_path / "dest", source_run_id="../escape")


def test_existing_conflicting_copy_is_not_overwritten(tmp_path):
    source = tmp_path / "source-run"
    source.mkdir()
    (source / "data").write_text("new")
    destination_run = tmp_path / "destination-run"
    conflicting = destination_run / "artifacts/legacy/source-run"
    conflicting.mkdir(parents=True)
    (conflicting / "data").write_text("old")

    with pytest.raises(LegacyMDStorageError, match="no valid manifest"):
        _stage(source, destination_run)
    assert (conflicting / "data").read_text() == "old"


def test_source_change_during_copy_fails_and_discards_staging(tmp_path, monkeypatch):
    import mn_molecular_design_studio.legacy_md_storage as storage

    source = tmp_path / "source-run"
    source.mkdir()
    source_file = source / "data"
    source_file.write_bytes(b"before")
    destination_run = tmp_path / "destination-run"
    original_sha = storage._sha256_file
    changed = False

    def mutate_after_hash(path):
        nonlocal changed
        digest = original_sha(path)
        if path == source_file and not changed:
            changed = True
            source_file.write_bytes(b"changed-after-hash")
        return digest

    monkeypatch.setattr(storage, "_sha256_file", mutate_after_hash)
    with pytest.raises(LegacyMDStorageError, match="changed while hashing"):
        _stage(source, destination_run)
    assert not (destination_run / "artifacts/legacy/source-run").exists()


def test_persisted_copy_corruption_is_detected_before_final_rename(tmp_path, monkeypatch):
    import mn_molecular_design_studio.legacy_md_storage as storage

    source = tmp_path / "source-run"
    source.mkdir()
    (source / "data").write_bytes(b"original")
    destination_run = tmp_path / "destination-run"
    original_sha = storage._sha256_file

    def corrupt_staged_copy(path):
        if ".staging-" in str(path.parent) and path.name == "data":
            path.write_bytes(b"corrupted-on-disk")
        return original_sha(path)

    monkeypatch.setattr(storage, "_sha256_file", corrupt_staged_copy)
    with pytest.raises(LegacyMDStorageError, match="readback verification"):
        _stage(source, destination_run)
    assert not (destination_run / "artifacts/legacy/source-run").exists()
