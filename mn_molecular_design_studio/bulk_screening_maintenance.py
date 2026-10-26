"""Narrow, set-based removal of accidental automatic bulk-run indexing.

Callers must back up the Database and preserve external results/pose archives
before applying this maintenance operation. This module never deletes files.
"""
from __future__ import annotations

import json

from mn_molecular_design_studio.molecular_database import MolecularDatabase, _now


def remove_automatic_bulk_records(references, user, *, owner_id, project_id="", database=None, dry_run=True):
    from mn_molecular_design_studio.auth import require_job_access
    require_job_access(user, {"owner_id": owner_id, "project_id": project_id})
    references = sorted(set(tuple(str(part) for part in reference) for reference in references))
    if not references or any(len(r) != 2 or any(not p or p in {".", ".."} or "/" in p or "\\" in p for p in r) for r in references):
        raise ValueError("Provide explicit job references for this campaign.")
    database = database or MolecularDatabase()
    with database.connection() as c:
        c.execute("CREATE TEMP TABLE cleanup_runs (task TEXT,run TEXT,PRIMARY KEY(task,run)) WITHOUT ROWID")
        c.executemany("INSERT INTO cleanup_runs VALUES (?,?)", references)
        c.execute("CREATE TEMP TABLE cleanup_sources AS SELECT s.* FROM molecule_sources s "
                  "JOIN molecule_groups g ON g.id=s.group_id JOIN cleanup_runs r "
                  "ON r.task=json_extract(s.provenance,'$.task_group') AND r.run=json_extract(s.provenance,'$.run_id') "
                  "WHERE json_extract(s.provenance,'$.saved_selection') IS NOT 1 "
                  "AND json_extract(g.origin,'$._replica_node') IS NULL AND g.owner_id=? AND g.project_id=?",
                  (owner_id, project_id))
        groups = [r[0] for r in c.execute("SELECT DISTINCT group_id FROM cleanup_sources")]
        for group in groups:
            database.group(group, user, write=True)
        c.execute("CREATE TEMP TABLE cleanup_pairs (group_id TEXT,molecule_id TEXT,PRIMARY KEY(group_id,molecule_id)) WITHOUT ROWID")
        c.execute("INSERT INTO cleanup_pairs SELECT DISTINCT group_id,molecule_id FROM cleanup_sources")
        result = {"run_count": len(references), "automatic_sources": c.execute("SELECT count(*) FROM cleanup_sources").fetchone()[0],
                  "affected_entries": c.execute("SELECT count(*) FROM cleanup_pairs").fetchone()[0],
                  "dry_run": dry_run, "preserved_saved_selections": True}
        if dry_run:
            return result
        c.execute("DELETE FROM molecule_sources WHERE source_key IN (SELECT source_key FROM cleanup_sources)")
        for table in ("molecule_observations", "molecule_structures"):
            removed = c.execute(f"DELETE FROM {table} WHERE (group_id,molecule_id) IN (SELECT group_id,molecule_id FROM cleanup_pairs) "
                "AND json_extract(context,'$.saved_selection') IS NOT 1 AND EXISTS "
                "(SELECT 1 FROM cleanup_runs r WHERE r.task=json_extract(context,'$.task_group') AND r.run=json_extract(context,'$.run_id'))").rowcount
            result["removed_" + table.removeprefix("molecule_")] = removed
        c.execute("CREATE TEMP TABLE cleanup_members AS SELECT m.group_id,m.molecule_id FROM molecule_members m "
                  "JOIN cleanup_pairs p ON p.group_id=m.group_id AND p.molecule_id=m.molecule_id "
                  "JOIN molecule_groups g ON g.id=m.group_id WHERE json_extract(g.origin,'$.type')='job' "
                  "AND (json_extract(g.origin,'$.collection_kind')='inputs' OR EXISTS "
                  "(SELECT 1 FROM cleanup_runs r WHERE r.task=json_extract(g.origin,'$.task_group') AND r.run=json_extract(g.origin,'$.run_id'))) "
                  "AND json_extract(g.origin,'$.collection_kind') IS NOT 'curated' "
                  "AND NOT EXISTS(SELECT 1 FROM molecule_sources s WHERE s.group_id=m.group_id AND s.molecule_id=m.molecule_id) "
                  "AND NOT EXISTS(SELECT 1 FROM molecule_structures s WHERE s.group_id=m.group_id AND s.molecule_id=m.molecule_id) "
                  "AND NOT EXISTS(SELECT 1 FROM molecule_observations o WHERE o.group_id=m.group_id AND o.molecule_id=m.molecule_id)")
        result["removed_members"] = c.execute("DELETE FROM molecule_members WHERE (group_id,molecule_id) IN (SELECT group_id,molecule_id FROM cleanup_members)").rowcount
        # Chemical identities shared with any curated group, import report or
        # unrelated job survive. Only affected, fully unreferenced identities go.
        result["removed_identities"] = c.execute("DELETE FROM molecules WHERE id IN (SELECT molecule_id FROM cleanup_pairs) "
            "AND NOT EXISTS(SELECT 1 FROM molecule_members WHERE molecule_id=molecules.id) "
            "AND NOT EXISTS(SELECT 1 FROM molecule_sources WHERE molecule_id=molecules.id) "
            "AND NOT EXISTS(SELECT 1 FROM molecule_structures WHERE molecule_id=molecules.id) "
            "AND NOT EXISTS(SELECT 1 FROM molecule_observations WHERE molecule_id=molecules.id) "
            "AND NOT EXISTS(SELECT 1 FROM molecule_import_rows WHERE molecule_id=molecules.id)").rowcount
        for group in groups:
            c.execute("UPDATE molecule_groups SET updated_at=? WHERE id=?", (_now(), group))
        return result


def capture_bulk_deletions(replica, *, database=None, node_id):
    """Publish existing replica entry tombstones with one clock allocation.

    Avoid recapturing hundreds of thousands of removed molecules. Never create
    records for entries that have not previously entered peer replication.
    """
    database = database or MolecularDatabase()
    with replica.db() as peer, database.connection() as c:
        removed = []
        for row in peer.execute("SELECT key,value FROM records WHERE kind='molecule_entry' AND author=? AND value!='null'", (node_id,)):
            value = json.loads(row["value"])
            if value and not c.execute("SELECT 1 FROM molecule_members WHERE group_id=? AND molecule_id=?",
                                      (value["group_id"], value["molecule"]["id"])).fetchone():
                removed.append(row["key"])
        maximum = peer.execute("SELECT coalesce(max(counter),0),coalesce(max(seq),0) FROM records").fetchone()
        clock = replica.state("clock", 0)
        counter, seq = max(maximum[0], int(clock)), maximum[1]
        peer.executemany("INSERT OR REPLACE INTO records VALUES (?, 'molecule_entry', ?, ?, 'null', ?)",
                         [(key, node_id, counter + i, seq + i) for i, key in enumerate(removed, 1)])
        peer.execute("INSERT OR REPLACE INTO state VALUES ('clock',?)", (str(counter + len(removed)),))
    return len(removed)
