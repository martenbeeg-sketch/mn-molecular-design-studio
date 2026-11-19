"""Invalidate browser snapshots on real data edits, not unrelated SQL bookkeeping."""
from __future__ import annotations

TABLES = (
    "molecules", "molecule_groups", "molecule_members", "molecule_sources", "molecule_structures",
    "molecule_observations", "complex_review_annotations", "complex_reference_contexts",
    "entity_names", "entity_name_details",
)


def install_schema(connection):
    connection.execute("CREATE TABLE IF NOT EXISTS molecular_browser_revision (singleton INTEGER PRIMARY KEY CHECK(singleton=1), revision INTEGER NOT NULL)")
    connection.execute("INSERT OR IGNORE INTO molecular_browser_revision VALUES (1,0)")
    present = {row[0] for row in connection.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    for table in TABLES:
        if table not in present:
            continue
        columns = [row[1] for row in connection.execute(f'PRAGMA table_info("{table}")')]
        changed = " OR ".join(f'OLD."{column}" IS NOT NEW."{column}"' for column in columns)
        for action in ("insert", "update", "delete"):
            condition = " WHEN " + changed if action == "update" else ""
            connection.execute(f'CREATE TRIGGER IF NOT EXISTS browser_revision_v1_{table}_{action} '
                f'AFTER {action.upper()} ON "{table}"{condition} BEGIN '
                "UPDATE molecular_browser_revision SET revision=revision+1 WHERE singleton=1; END")


def revision(connection):
    if not connection.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='molecular_browser_revision'").fetchone():
        return None
    row = connection.execute("SELECT revision FROM molecular_browser_revision WHERE singleton=1").fetchone()
    return row[0] if row else None
