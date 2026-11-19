import pytest

from mn_molecular_design_studio.cluster import replication
from test_molecular_database_peer import roots


def test_large_molecular_snapshots_use_bounded_pages_without_skipping_records(roots):
    replica = replication.Replica()
    author = replication.config.identity()
    with replica.db() as connection:
        for number in range(12):
            replica._write(connection, f"molecule_entry:{author}:{number}", "molecule_entry", {"payload": "x" * 1_500_000})
    first = replica.export()
    assert first["more"] and len(first["records"]) < 12
    assert len(replication.encoded(first["records"]).encode()) <= replication.MAX_METADATA_PAGE_BYTES
    second = replica.export(first["cursor"])
    assert not second["more"]
    assert len({row["key"] for row in first["records"] + second["records"]}) == 12


def test_large_molecular_record_is_accepted_but_other_metadata_stays_bounded(roots, monkeypatch):
    monkeypatch.setattr(replication, "capture_identity", lambda: None)
    monkeypatch.setattr(replication.Replica, "_apply", lambda *args, **kwargs: None)
    replica = replication.Replica()
    author = replication.config.identity()
    record = {"kind": "molecule_entry", "key": f"molecule_entry:{author}:test", "author": author,
        "counter": 1, "value": {"payload": "x" * (3 * 1024 * 1024)}}
    assert replica.receive([record]) == 1
    record.update(kind="project", key="project:test", counter=2)
    with pytest.raises(ValueError, match="size limit"):
        replica.receive([record])


def test_oversized_metadata_pages_are_rejected_before_application(roots, monkeypatch):
    monkeypatch.setattr(replication, "capture_identity", lambda: pytest.fail("Oversized pages must fail before writes."))
    replica = replication.Replica()
    with pytest.raises(ValueError, match="page exceeds"):
        replica.receive([{"value": "x" * replication.MAX_METADATA_PAGE_BYTES}])
