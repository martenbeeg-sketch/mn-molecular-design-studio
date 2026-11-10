"""Versioned, rebuildable scientific lineage projection in the existing Database.

Job ownership is authoritative. Snapshots contain metadata, never file bytes or
new canonical molecular identities. Existing job replication carries snapshots.
"""
from __future__ import annotations

import hashlib
import json
import re

VERSION = "studio.scientific-lineage.v1"
KINDS = {"design", "folding", "ligand_cofolding", "docking", "dynamics", "analysis", "unknown"}
MAX_BYTES = 1_500_000  # Existing replication record envelope, not a campaign limit.
PAGE_BYTES = 700_000   # Keep a generous margin below the existing 2 MB peer-record limit.


def encoded(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _publication_revision(manifest, page_hashes):
    # Count/acceptance/construct changes must publish even when retained links
    # themselves did not change. Legacy v1 manifests used only page hashes.
    if "source_summary" not in manifest:
        value = page_hashes
    else:
        value = {"page_sha256": page_hashes, "source_summary": manifest.get("source_summary"),
                 "target": manifest.get("target"),
                 "counts": manifest.get("counts")}
    return hashlib.sha256(encoded(value).encode()).hexdigest()


def acceptance_counts(outputs):
    # Acceptance is meaningful for candidate/model/pose occurrences, not MD
    # replica/segment records or Evaluation rows.
    rows = [o for o in outputs if o.get("kind") not in {
        "evaluation", "retained_representation", "md_replica", "md_segment", "md_study",
    }]
    stages = sorted({d["stage_id"] for o in rows for d in o.get("acceptance", []) if d.get("stage_id")})
    result = {}
    for stage in stages:
        counts = dict(passed=0, rejected=0, unresolved=0, unrecorded=0)
        for row in rows:
            states = {d.get("cohort_status") for d in row.get("acceptance", []) if d.get("stage_id") == stage}
            status = "unrecorded" if not states else next(iter(states)) if len(states) == 1 else "unresolved"
            counts[status if status in counts else "unresolved"] += 1
        result[stage] = dict(counts, unit="candidate/model/pose occurrences", scope="source-authored saved outputs")
    return result


def _source_summary(snapshot):
    """Compact source-authored counts; detailed Results remain owner-local."""
    outputs = [o for o in snapshot.get("outputs", [])
               if o.get("kind") not in {"evaluation", "retained_representation"}]
    by_kind, by_stage = {}, {}
    for row in outputs:
        kind = str(row.get("kind") or "unknown")
        stage = str(row.get("stage") or "unrecorded")
        by_kind[kind] = by_kind.get(kind, 0) + 1
        by_stage[stage] = by_stage.get(stage, 0) + 1
    candidate_kinds = {"design", "backbone", "sequence_variant", "prediction_model", "docking_pose"}
    retained_occurrences = [o for o in outputs if o.get("kind") in candidate_kinds and o.get("retained")]
    retained_structure_refs = sum(len(o.get("retained") or []) for o in retained_occurrences)
    # Preserve producer-reported count definitions and the normalized row count
    # independently; neither is relabeled as another acceptance stage.
    return {
        "source_output_occurrences": len(outputs),
        "candidate_model_pose_occurrences": sum(by_kind.get(kind, 0) for kind in candidate_kinds),
        "generated_design_occurrences": by_kind.get("design", 0) + by_kind.get("backbone", 0),
        "output_kind_counts": by_kind,
        "stage_occurrence_counts": by_stage,
        "acceptance_counts": acceptance_counts(outputs),
        "producer_counts": list(snapshot.get("counts") or []),
        "database_retained_candidate_model_occurrences": len(retained_occurrences),
        "database_retained_structure_references": retained_structure_refs,
        "mapping_warning_counts": _warning_counts(snapshot.get("warnings") or []),
        "source_summary_complete": True,
    }


def _warning_counts(warnings):
    categories = {}
    for warning in warnings:
        text = str(warning)
        if text.startswith("No exact retained structure mapping"):
            key = "candidate_model_without_exact_retained_structure"
        elif text.startswith("Repeated candidate/model key"):
            key = "repeated_candidate_model_occurrence"
        elif text.startswith("Multiple target constructs"):
            key = "multiple_target_constructs"
        elif text.startswith("Parent run recorded"):
            key = "parent_link_missing_source_workstation"
        elif text.startswith("Retained representation"):
            key = "retained_representation_source_mapping_unresolved"
        elif text.startswith("Recorded target contexts differ"):
            key = "target_contexts_differ"
        else:
            key = "other_recorded_mapping_warning"
        categories[key] = categories.get(key, 0) + 1
    return categories


def _retained_link(output):
    """Publish a relational pointer, never a second copy of a Database design."""
    retained = []
    for row in output.get("retained") or []:
        if not all(row.get(key) for key in ("molecule_id", "structure_id", "source_node", "owner_structure_id")):
            continue
        retained.append({key: row[key] for key in (
            "molecule_id", "structure_id", "sha256", "source_node", "owner_structure_id",
        ) if key in row})
    if not retained:
        return None
    return {
        "ref": output["ref"], "job_ref": output.get("job_ref", ""),
        "kind": output.get("kind", "unknown"),
        "candidate_id": str(output.get("candidate_id") or ""),
        "native_id": str(output.get("native_id") or output.get("candidate_id") or ""),
        "model_id": str(output.get("model_id") or ""),
        "stage": str(output.get("stage") or ""),
        "method": str(output.get("method") or ""),
        "acceptance": list(output.get("acceptance") or []),
        "target_construct_key": _target_construct_key(output["target_construct"])
            if output.get("target_construct") else "",
        "retained": retained,
        "data_location": recorded_data_location(output),
    }


def _job_relationships(snapshot):
    """Keep only compact job/controller relationships, not output file/score edges."""
    result = []
    job_ref = snapshot.get("job", {}).get("ref", "")
    for edge in snapshot.get("dependencies", []):
        if edge.get("role") != "execution_child_of":
            continue
        source, target = str(edge.get("source_ref") or ""), str(edge.get("target_ref") or "")
        if not source.startswith(("cluster://", "cluster-job:", "node-run:")) or not target.startswith(("cluster://", "cluster-job:", "node-run:")):
            continue
        result.append({"source_ref": source, "target_ref": target, "role": "execution_child_of",
                       "status": str(edge.get("status") or "recorded"),
                       "evidence": {key: value for key, value in (edge.get("evidence") or {}).items()
                                    if key in {"source", "attempt_id"}}})
    # Parent links may be held in job metadata; snapshot construction already
    # emits their explicit source-qualified or unresolved execution edge.
    return sorted(result, key=lambda row: (row["source_ref"], row["target_ref"]))


def identifier(prefix, *values):
    return prefix + hashlib.sha256(encoded(values).encode()).hexdigest()[:32]


def job_reference(node, task_group, run_id):
    if any(not isinstance(p, str) or not re.fullmatch(r"[A-Za-z0-9_.-]+", p) or p in {".", ".."} for p in (node, task_group, run_id)):
        raise ValueError("A source-qualified job identity is required.")
    return f"cluster://{node}/{task_group}/{run_id}"


def install_schema(c):
    c.executescript("""
        CREATE TABLE IF NOT EXISTS scientific_jobs (
            ref TEXT PRIMARY KEY, node_id TEXT NOT NULL, owner_id TEXT NOT NULL,
            project_id TEXT NOT NULL, task_group TEXT NOT NULL, run_id TEXT NOT NULL,
            cluster_id TEXT NOT NULL, snapshot TEXT NOT NULL,
            job_metadata TEXT NOT NULL DEFAULT '{}', facts TEXT NOT NULL DEFAULT '{}',
            publication_revision TEXT NOT NULL DEFAULT '', publication_pages INTEGER NOT NULL DEFAULT 0,
            publication_state TEXT NOT NULL DEFAULT 'legacy',
            UNIQUE(node_id,task_group,run_id));
        CREATE INDEX IF NOT EXISTS scientific_job_author ON scientific_jobs(node_id);
        CREATE INDEX IF NOT EXISTS scientific_job_cluster ON scientific_jobs(cluster_id);
        CREATE TABLE IF NOT EXISTS scientific_outputs (
            ref TEXT PRIMARY KEY, job_ref TEXT NOT NULL REFERENCES scientific_jobs(ref) ON DELETE CASCADE,
            kind TEXT NOT NULL, candidate_id TEXT NOT NULL, model_id TEXT NOT NULL,
            stage TEXT NOT NULL, payload TEXT NOT NULL);
        CREATE INDEX IF NOT EXISTS scientific_output_job ON scientific_outputs(job_ref,ref);
        CREATE TABLE IF NOT EXISTS scientific_dependencies (
            job_ref TEXT NOT NULL REFERENCES scientific_jobs(ref) ON DELETE CASCADE,
            source_ref TEXT NOT NULL, target_ref TEXT NOT NULL, role TEXT NOT NULL,
            status TEXT NOT NULL, evidence TEXT NOT NULL,
            PRIMARY KEY(job_ref,source_ref,target_ref,role));
        CREATE INDEX IF NOT EXISTS scientific_dependency_target ON scientific_dependencies(target_ref);
        CREATE INDEX IF NOT EXISTS scientific_dependency_source ON scientific_dependencies(source_ref);
        CREATE TABLE IF NOT EXISTS scientific_lineage_pages (
            job_ref TEXT NOT NULL REFERENCES scientific_jobs(ref) ON DELETE CASCADE,
            page_index INTEGER NOT NULL, publication_revision TEXT NOT NULL,
            page_sha256 TEXT NOT NULL, payload TEXT NOT NULL,
            PRIMARY KEY(job_ref,page_index));
        CREATE TABLE IF NOT EXISTS scientific_job_aliases (
            alias TEXT NOT NULL, job_ref TEXT NOT NULL REFERENCES scientific_jobs(ref) ON DELETE CASCADE,
            PRIMARY KEY(alias,job_ref));
        CREATE INDEX IF NOT EXISTS scientific_job_alias_lookup ON scientific_job_aliases(alias,job_ref);
        CREATE TABLE IF NOT EXISTS scientific_execution_edges (
            child_ref TEXT NOT NULL REFERENCES scientific_jobs(ref) ON DELETE CASCADE,
            parent_ref TEXT NOT NULL REFERENCES scientific_jobs(ref) ON DELETE CASCADE,
            PRIMARY KEY(child_ref,parent_ref));
        CREATE INDEX IF NOT EXISTS scientific_execution_parent ON scientific_execution_edges(parent_ref,child_ref);
        CREATE TABLE IF NOT EXISTS scientific_campaigns (
            ref TEXT PRIMARY KEY, owner_id TEXT NOT NULL, project_id TEXT NOT NULL,
            kind TEXT NOT NULL, name TEXT NOT NULL, target_key TEXT NOT NULL,
            updated_at TEXT NOT NULL, summary TEXT NOT NULL);
        CREATE INDEX IF NOT EXISTS scientific_campaign_access
            ON scientific_campaigns(project_id,owner_id,kind,target_key,updated_at DESC,ref);
        CREATE TABLE IF NOT EXISTS scientific_campaign_jobs (
            campaign_ref TEXT NOT NULL REFERENCES scientific_campaigns(ref) ON DELETE CASCADE,
            job_ref TEXT NOT NULL REFERENCES scientific_jobs(ref) ON DELETE CASCADE,
            PRIMARY KEY(campaign_ref,job_ref));
        CREATE INDEX IF NOT EXISTS scientific_campaign_job ON scientific_campaign_jobs(job_ref);
        CREATE TABLE IF NOT EXISTS scientific_target_constructs (
            target_key TEXT PRIMARY KEY, biological_identity_key TEXT NOT NULL DEFAULT '',
            description TEXT NOT NULL DEFAULT '', target_json TEXT NOT NULL);
        CREATE INDEX IF NOT EXISTS scientific_target_biological
            ON scientific_target_constructs(biological_identity_key,target_key);
        CREATE TABLE IF NOT EXISTS scientific_campaign_targets (
            campaign_ref TEXT NOT NULL REFERENCES scientific_campaigns(ref) ON DELETE CASCADE,
            target_key TEXT NOT NULL REFERENCES scientific_target_constructs(target_key),
            relationship_json TEXT NOT NULL DEFAULT '{}', PRIMARY KEY(campaign_ref,target_key));
        CREATE INDEX IF NOT EXISTS scientific_campaign_target_lookup
            ON scientific_campaign_targets(target_key,campaign_ref);
    """)
    columns = {r[1] for r in c.execute("PRAGMA table_info(scientific_jobs)")}
    for name in ("job_metadata", "facts"):
        if name not in columns:
            c.execute(f"ALTER TABLE scientific_jobs ADD COLUMN {name} TEXT NOT NULL DEFAULT '{{}}'")
    for name, declaration in (("publication_revision", "TEXT NOT NULL DEFAULT ''"),
                              ("publication_pages", "INTEGER NOT NULL DEFAULT 0"),
                              ("publication_state", "TEXT NOT NULL DEFAULT 'legacy'")):
        if name not in columns:
            c.execute(f"ALTER TABLE scientific_jobs ADD COLUMN {name} {declaration}")


def validate_snapshot(snapshot, author, envelope=None):
    if not isinstance(snapshot, dict) or snapshot.get("schema") != VERSION:
        raise ValueError("Unsupported scientific lineage snapshot.")
    job = snapshot.get("job") or {}
    if not isinstance(job, dict):
        raise ValueError("Scientific job metadata must be an object.")
    if job.get("node_id") != author or job.get("ref") != job_reference(author, job.get("task_group", ""), job.get("run_id", "")):
        raise PermissionError("Only the scientific source workstation may publish its job.")
    if job.get("kind") not in KINDS:
        raise ValueError("Invalid scientific lineage snapshot.")
    if envelope is not None and any(str(job.get(k) or "") != str(envelope.get(k) or "")
                                    for k in ("owner_id", "project_id", "node_id", "run_id", "task_group")):
        raise PermissionError("Scientific lineage ownership differs from its job envelope.")
    # Full source snapshots are split into bounded, replayable pages before
    # they enter replication. Manifests contain no output population.
    if "publication" in snapshot:
        publication = snapshot.get("publication") or {}
        if (publication.get("job_ref") != job.get("ref")
                or type(publication.get("page_count")) is not int
                or type(publication.get("output_count")) is not int
                or type(publication.get("dependency_count")) is not int
                or min(publication.get("page_count", -1), publication.get("output_count", -1),
                       publication.get("dependency_count", -1)) < 0
                or ("source_summary" in snapshot and
                    (not isinstance(snapshot.get("source_summary"), dict)
                     or type(snapshot["source_summary"].get("source_output_occurrences")) is not int
                     or snapshot["source_summary"].get("source_output_occurrences", -1) < 0
                     or not isinstance(snapshot["source_summary"].get("acceptance_counts", {}), dict)))
                or ("source_output_count" in publication and
                    (type(publication["source_output_count"]) is not int or publication["source_output_count"] < 0))
                or ("source_summary" in snapshot and
                    (publication.get("source_output_count") != snapshot["source_summary"].get("source_output_occurrences")
                     or publication.get("output_count") != snapshot["source_summary"].get("database_retained_candidate_model_occurrences")))
                or len(encoded(snapshot).encode()) > MAX_BYTES):
            raise ValueError("Invalid or oversized scientific lineage publication manifest.")
        page_hashes = publication.get("page_sha256")
        if (not isinstance(page_hashes, list) or len(page_hashes) != publication["page_count"]
                or any(not isinstance(value, str) or not re.fullmatch(r"[a-f0-9]{64}", value) for value in page_hashes)
                or _publication_revision(snapshot, page_hashes) != publication.get("revision")):
            raise ValueError("Scientific publication manifest hashes are inconsistent.")
        return job
    outputs = snapshot.get("outputs", [])
    if not isinstance(outputs, list) or any(not isinstance(o, dict) for o in outputs):
        raise ValueError("Scientific output records must be objects.")
    refs = [o.get("ref") for o in outputs]
    if any(not isinstance(r, str) for r in refs) or len(refs) != len(set(refs)) or any(not r.startswith(job["ref"] + "#") for r in refs):
        raise ValueError("Output occurrences must have unique job-qualified identities.")
    dependencies = snapshot.get("dependencies", [])
    if not isinstance(dependencies, list) or any(not isinstance(e, dict) or not all(isinstance(e.get(k), str) and e[k] for k in ("source_ref", "target_ref", "role"))
           for e in dependencies):
        raise ValueError("Typed scientific input references are required.")
    return job


def _page_rows(snapshot):
    rows = ([{"type": "output", "record": row} for row in snapshot.get("outputs", [])]
            + [{"type": "dependency", "record": row} for row in snapshot.get("dependencies", [])])
    rows.sort(key=lambda row: (row["type"], str(row["record"].get("ref") or row["record"].get("source_ref") or ""),
                               str(row["record"].get("target_ref") or "")))
    pages, current, size = [], [], 2
    for row in rows:
        encoded_row = encoded(row)
        if len(encoded_row.encode()) + 2 > PAGE_BYTES:
            raise ValueError("A single scientific lineage occurrence exceeds the safe peer page size.")
        if current and size + len(encoded_row.encode()) + 1 > PAGE_BYTES:
            pages.append(current)
            current, size = [], 2
        current.append(row)
        size += len(encoded_row.encode()) + 1
    if current:
        pages.append(current)
    return pages


def publication_for_snapshot(snapshot):
    """Publish aggregate source facts plus links for Database-retained outputs.

    Detailed rejected/unresolved candidates, unselected models, score tables and
    trajectory records stay in owner-local Results. The optional pages contain
    only relational pointers to already replicated Database structures.
    """
    job = snapshot["job"]
    summary = _source_summary(snapshot)
    links = [_retained_link(row) for row in snapshot.get("outputs", [])]
    links = [row for row in links if row is not None]
    compact = {"outputs": links, "dependencies": _job_relationships(snapshot)}
    pages = _page_rows(compact)
    page_hashes = [hashlib.sha256(encoded(page).encode()).hexdigest() for page in pages]
    manifest = {key: value for key, value in snapshot.items()
                if key not in {"outputs", "dependencies", "warnings"}}
    manifest["warnings"] = []
    manifest["source_summary"] = summary
    revision = _publication_revision(manifest, page_hashes)
    manifest["publication"] = {"job_ref": job["ref"], "revision": revision,
        "page_count": len(pages), "output_count": len(links),
        "dependency_count": len(compact["dependencies"]),
        "source_output_count": summary["source_output_occurrences"],
        "page_sha256": page_hashes}
    envelope_pages = [{"schema": VERSION, "job_ref": job["ref"],
        "page_index": index, "page_sha256": page_hashes[index], "items": page}
        for index, page in enumerate(pages)]
    if len(encoded(manifest).encode()) > MAX_BYTES:
        raise ValueError("Scientific lineage manifest is too large to replicate safely.")
    return manifest, envelope_pages


def put_snapshot(database, snapshot, *, author, envelope=None, rebuild=True):
    if snapshot.get("publication"):
        raise ValueError("Use put_manifest for an already paged scientific publication.")
    manifest, pages = publication_for_snapshot(snapshot)
    put_manifest(database, manifest, author=author, envelope=envelope, rebuild=rebuild)
    for page in pages:
        put_publication_page(database, page, author=author, rebuild=rebuild)
    return manifest["publication"]


def put_manifest(database, manifest, *, author, envelope=None, rebuild=True):
    job = validate_snapshot(manifest, author, envelope)
    publication = manifest["publication"]
    source_summary = manifest.get("source_summary") or {}
    facts = {"target": manifest.get("target") or {}, "counts": manifest.get("counts") or [],
             "warnings": [], "source_summary": source_summary, "indexed_occurrences": 0,
             "expected_occurrences": publication["output_count"],
             "source_output_occurrences": int(source_summary.get("source_output_occurrences", 0)),
             "scientific_output_occurrences": int(source_summary.get("source_output_occurrences", 0)),
             "acceptance_counts": source_summary.get("acceptance_counts", {}),
             "retained_refs": [], "publication_state": "incomplete",
             "publication_revision": publication["revision"], "publication_pages": publication["page_count"]}
    with database.connection() as c:
        previous = c.execute("SELECT owner_id,project_id,publication_revision,publication_state,facts FROM scientific_jobs WHERE ref=?", (job["ref"],)).fetchone()
        if previous and (previous["owner_id"], previous["project_id"]) != (job.get("owner_id", ""), job.get("project_id", "")):
            raise PermissionError("A lineage update cannot transfer project or user ownership.")
        if previous and previous["publication_revision"] == publication["revision"] and previous["publication_state"] == "complete":
            prior_facts = json.loads(previous["facts"])
            facts.update({key: prior_facts[key] for key in ("acceptance_counts", "retained_refs", "indexed_occurrences",
                "scientific_output_occurrences", "source_output_occurrences", "source_summary") if key in prior_facts})
            facts["publication_state"] = "complete"
        c.execute("INSERT INTO scientific_jobs (ref,node_id,owner_id,project_id,task_group,run_id,cluster_id,snapshot,job_metadata,facts,publication_revision,publication_pages,publication_state) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?) ON CONFLICT(ref) DO UPDATE SET "
                  "cluster_id=excluded.cluster_id,snapshot=excluded.snapshot,job_metadata=excluded.job_metadata,facts=excluded.facts,"
                  "publication_revision=excluded.publication_revision,publication_pages=excluded.publication_pages,publication_state=excluded.publication_state",
                  (job["ref"], author, job.get("owner_id", ""), job.get("project_id", ""),
                   job["task_group"], job["run_id"], job.get("cluster_job_id", ""), encoded(manifest), encoded(job), encoded(facts),
                   publication["revision"], publication["page_count"], facts["publication_state"]))
        _index_job_aliases(c, job)
        if publication["page_count"] == 0:
            if publication["output_count"] or publication["dependency_count"]:
                raise ValueError("Empty scientific publication has nonzero expected counts.")
            c.execute("DELETE FROM scientific_outputs WHERE job_ref=?", (job["ref"],))
            c.execute("DELETE FROM scientific_dependencies WHERE job_ref=?", (job["ref"],))
            c.execute("DELETE FROM scientific_lineage_pages WHERE job_ref=?", (job["ref"],))
            facts.update(publication_state="complete", publication_revision=publication["revision"],
                         indexed_occurrences=0, expected_occurrences=0,
                         scientific_output_occurrences=int(source_summary.get("source_output_occurrences", 0)),
                         source_output_occurrences=int(source_summary.get("source_output_occurrences", 0)),
                         acceptance_counts=source_summary.get("acceptance_counts", {}))
            c.execute("UPDATE scientific_jobs SET facts=?,publication_state='complete' WHERE ref=?", (encoded(facts), job["ref"]))
        else:
            _try_finalize_publication(c, job["ref"], rebuild=False)
        if rebuild:
            rebuild_campaigns(c, affected_job_ref=job["ref"])


def _validate_page(page, author):
    if not isinstance(page, dict) or page.get("schema") != VERSION:
        raise ValueError("Unsupported scientific lineage page.")
    ref = page.get("job_ref") or ""
    if not ref.startswith(f"cluster://{author}/") or page.get("page_index", -1) < 0:
        raise PermissionError("Only the source workstation may publish this scientific page.")
    if not isinstance(page.get("items"), list) or len(encoded(page).encode()) > MAX_BYTES:
        raise ValueError("Invalid or oversized scientific lineage page.")
    digest = hashlib.sha256(encoded(page["items"]).encode()).hexdigest()
    if digest != page.get("page_sha256"):
        raise ValueError("Scientific lineage page digest mismatch.")
    for item in page["items"]:
        if not isinstance(item, dict) or item.get("type") not in {"output", "dependency"} or not isinstance(item.get("record"), dict):
            raise ValueError("Invalid scientific lineage page row.")
    return ref


def put_publication_page(database, page, *, author, rebuild=True):
    ref = _validate_page(page, author)
    with database.connection() as c:
        job = c.execute("SELECT owner_id,project_id,node_id,publication_revision,publication_pages,snapshot FROM scientific_jobs WHERE ref=?", (ref,)).fetchone()
        if not job or job["node_id"] != author:
            raise PermissionError("Scientific page has no authorized source job manifest.")
        manifest = json.loads(job["snapshot"])
        publication = manifest.get("publication") or {}
        if "source_summary" in manifest:
            output_fields = {"ref", "job_ref", "kind", "candidate_id", "native_id", "model_id",
                             "stage", "method", "acceptance", "target_construct_key", "retained", "data_location"}
            retained_fields = {"molecule_id", "structure_id", "sha256", "source_node", "owner_structure_id"}
            for item in page.get("items", []):
                record = item.get("record") or {}
                if item["type"] == "output":
                    if set(record) - output_fields or any(set(row) - retained_fields for row in record.get("retained", [])):
                        raise ValueError("Scientific campaign pages may contain only compact retained Database links.")
                elif record.get("role") != "execution_child_of":
                    raise ValueError("Scientific campaign pages may contain only compact job relationships.")
        expected_hashes = publication.get("page_sha256") or []
        if page["page_index"] >= job["publication_pages"] or page["page_sha256"] != expected_hashes[page["page_index"]]:
            raise ValueError("Scientific page does not match the active publication manifest.")
        if job["publication_revision"] == publication.get("revision") and job["publication_pages"] and c.execute(
            "SELECT publication_state FROM scientific_jobs WHERE ref=?", (ref,)).fetchone()[0] == "complete":
            return True
        c.execute("INSERT INTO scientific_lineage_pages VALUES (?,?,?,?,?) ON CONFLICT(job_ref,page_index) DO UPDATE SET "
                  "publication_revision=excluded.publication_revision,page_sha256=excluded.page_sha256,payload=excluded.payload",
                  (ref, page["page_index"], job["publication_revision"], page["page_sha256"], encoded(page)))
        return _try_finalize_publication(c, ref, rebuild=rebuild)


def _try_finalize_publication(c, ref, *, rebuild):
    job = c.execute("SELECT publication_revision,publication_pages,publication_state,snapshot,facts FROM scientific_jobs WHERE ref=?", (ref,)).fetchone()
    if not job:
        return False
    manifest = json.loads(job["snapshot"])
    publication = manifest.get("publication") or {}
    expected_count = int(job["publication_pages"])
    if expected_count == 0:
        return True
    if job["publication_state"] == "complete" and json.loads(job["facts"]).get("publication_revision") == job["publication_revision"]:
        return True
    ready = {r["page_index"]: r for r in c.execute(
        "SELECT page_index,page_sha256,payload FROM scientific_lineage_pages WHERE job_ref=?", (ref,))}
    expected_hashes = publication.get("page_sha256") or []
    if len(expected_hashes) != expected_count or any(i not in ready or ready[i]["page_sha256"] != expected_hashes[i] for i in range(expected_count)):
        return False
    actual_hashes = [ready[i]["page_sha256"] for i in range(expected_count)]
    if _publication_revision(manifest, actual_hashes) != publication.get("revision"):
        raise ValueError("Scientific publication digest mismatch.")
    outputs, dependencies = [], []
    for index in range(expected_count):
        stored = json.loads(ready[index]["payload"])
        if hashlib.sha256(encoded(stored["items"]).encode()).hexdigest() != expected_hashes[index]:
            raise ValueError("Stored scientific page failed digest validation.")
        for item in stored["items"]:
            (outputs if item["type"] == "output" else dependencies).append(item["record"])
    if len(outputs) != publication.get("output_count") or len(dependencies) != publication.get("dependency_count"):
        raise ValueError("Scientific publication is incomplete; expected occurrence counts differ.")
    c.execute("DELETE FROM scientific_outputs WHERE job_ref=?", (ref,))
    c.execute("DELETE FROM scientific_dependencies WHERE job_ref=?", (ref,))
    for output in outputs:
        c.execute("INSERT INTO scientific_outputs VALUES (?,?,?,?,?,?,?)",
                  (output["ref"], ref, output["kind"], str(output.get("candidate_id") or ""),
                   str(output.get("model_id") or ""), str(output.get("stage") or ""), encoded(output)))
    for edge in dependencies:
        c.execute("INSERT OR IGNORE INTO scientific_dependencies VALUES (?,?,?,?,?,?)",
                  (ref, edge["source_ref"], edge["target_ref"], edge["role"], edge.get("status", "unresolved"), encoded(edge.get("evidence") or {})))
    current_job = json.loads(c.execute("SELECT job_metadata FROM scientific_jobs WHERE ref=?", (ref,)).fetchone()[0])
    _index_job_aliases(c, current_job)
    facts = json.loads(job["facts"])
    source_summary = manifest.get("source_summary") or {}
    facts.update(publication_state="complete", publication_revision=job["publication_revision"],
        indexed_occurrences=len(outputs), expected_occurrences=len(outputs),
        source_output_occurrences=int(source_summary.get("source_output_occurrences", len(outputs))),
        scientific_output_occurrences=int(source_summary.get("source_output_occurrences", len(outputs))),
        source_summary=source_summary,
        acceptance_counts=source_summary.get("acceptance_counts", acceptance_counts(outputs)),
        retained_refs=sorted({(r.get("source_node", ""), r.get("owner_structure_id", "")) for o in outputs for r in o.get("retained", [])}))
    c.execute("UPDATE scientific_jobs SET facts=?,publication_state='complete' WHERE ref=?", (encoded(facts), ref))
    c.execute("DELETE FROM scientific_lineage_pages WHERE job_ref=? AND page_index>=?", (ref, expected_count))
    # Keep source pages only on the source workstation for idempotent publish
    # and changed-page retries. Receivers retain normalized rows, not a second
    # serialized copy of the same output population.
    try:
        from mn_molecular_design_studio.cluster.config import identity as local_identity
        local_is_author = local_identity() == c.execute("SELECT node_id FROM scientific_jobs WHERE ref=?", (ref,)).fetchone()[0]
    except Exception:
        local_is_author = False
    if not local_is_author:
        c.execute("DELETE FROM scientific_lineage_pages WHERE job_ref=?", (ref,))
    if rebuild:
        rebuild_campaigns(c, affected_job_ref=ref)
    return True


def _index_job_aliases(c, job):
    ref = job["ref"]
    aliases = {job["ref"], f"node-run:{job['node_id']}:{job['run_id']}"}
    if job.get("cluster_job_id"):
        aliases.add("cluster-job:" + str(job["cluster_job_id"]))
    # These explicit identifiers are also grouping evidence in logical_runs.
    # Scope them by owner/project and hash the tuple so names never become
    # scientific identity or cross-project aliases.
    for field in ("studio_run_id", "campaign_group_id", "distributed_campaign_id"):
        value = str(job.get(field) or "")
        if value:
            aliases.add("logical-group:" + identifier("", job.get("owner_id", ""),
                job.get("project_id", ""), field, value))
    c.execute("DELETE FROM scientific_job_aliases WHERE job_ref=?", (ref,))
    c.executemany("INSERT OR IGNORE INTO scientific_job_aliases VALUES (?,?)", [(alias, ref) for alias in aliases])
    # Resolve only edges touching this job's aliases. A newly arrived worker or
    # controller can complete a previously unresolved relationship.
    placeholders = ",".join("?" for _ in aliases)
    edges = c.execute(f"SELECT job_ref,source_ref,target_ref FROM scientific_dependencies WHERE role='execution_child_of' AND (source_ref IN ({placeholders}) OR target_ref IN ({placeholders}))",
                      [*aliases, *aliases]).fetchall()
    for edge in edges:
        source = c.execute("SELECT job_ref FROM scientific_job_aliases WHERE alias=?", (edge["source_ref"],)).fetchall()
        target = c.execute("SELECT job_ref FROM scientific_job_aliases WHERE alias=?", (edge["target_ref"],)).fetchall()
        if len(source) != 1 or len(target) != 1 or source[0][0] == target[0][0]:
            continue
        child, parent = source[0][0], target[0][0]
        child_job = c.execute("SELECT owner_id,project_id FROM scientific_jobs WHERE ref=?", (child,)).fetchone()
        parent_job = c.execute("SELECT owner_id,project_id FROM scientific_jobs WHERE ref=?", (parent,)).fetchone()
        if child_job and parent_job and tuple(child_job) == tuple(parent_job):
            c.execute("INSERT OR IGNORE INTO scientific_execution_edges VALUES (?,?)", (child, parent))


def rebuild_campaigns(c, *, affected_job_ref=None):
    """Incrementally regroup only the affected execution component at write time."""
    from mn_molecular_design_studio.logical_runs import group_runs
    old_campaigns = set()
    if affected_job_ref:
        old_campaigns = {r[0] for r in c.execute("SELECT campaign_ref FROM scientific_campaign_jobs WHERE job_ref=?", (affected_job_ref,))}
        impacted = {affected_job_ref}
        # Jobs sharing an explicit Studio run/group identity may not yet have
        # a resolved execution edge (for example a controller or analysis job
        # can arrive before its parent). Rebuild that indexed logical group
        # alongside the affected execution component.
        group_aliases = [r[0] for r in c.execute(
            "SELECT alias FROM scientific_job_aliases WHERE job_ref=? AND alias LIKE 'logical-group:%'",
            (affected_job_ref,))]
        for alias in group_aliases:
            impacted.update(r[0] for r in c.execute(
                "SELECT a.job_ref FROM scientific_job_aliases a JOIN scientific_jobs j ON j.ref=a.job_ref "
                "JOIN scientific_jobs source ON source.ref=? "
                "WHERE a.alias=? AND j.owner_id=source.owner_id AND j.project_id=source.project_id",
                (affected_job_ref, alias)))
        if old_campaigns:
            marks = ",".join("?" for _ in old_campaigns)
            impacted.update(r[0] for r in c.execute(f"SELECT job_ref FROM scientific_campaign_jobs WHERE campaign_ref IN ({marks})", sorted(old_campaigns)))
        queue = list(impacted)
        while queue:
            current = queue.pop()
            connected = {r[0] for r in c.execute(
                "SELECT parent_ref FROM scientific_execution_edges WHERE child_ref=? UNION "
                "SELECT child_ref FROM scientific_execution_edges WHERE parent_ref=?", (current, current))}
            for ref in connected - impacted:
                impacted.add(ref)
                queue.append(ref)
        if not impacted:
            return
        marks = ",".join("?" for _ in impacted)
        old_campaigns.update(r[0] for r in c.execute(
            f"SELECT DISTINCT campaign_ref FROM scientific_campaign_jobs WHERE job_ref IN ({marks})", sorted(impacted)))
        marks = ",".join("?" for _ in impacted)
        selected = c.execute(f"SELECT ref,job_metadata,facts,publication_state FROM scientific_jobs WHERE ref IN ({marks}) AND job_metadata!='{{}}'", sorted(impacted)).fetchall()
        c.execute(f"DELETE FROM scientific_campaign_jobs WHERE campaign_ref IN ({','.join('?' for _ in old_campaigns)})" if old_campaigns else "SELECT 1", sorted(old_campaigns)) if old_campaigns else None
        if old_campaigns:
            c.executemany("DELETE FROM scientific_campaign_targets WHERE campaign_ref=?", [(ref,) for ref in old_campaigns])
            c.executemany("DELETE FROM scientific_campaigns WHERE ref=?", [(ref,) for ref in old_campaigns])
    else:
        selected = c.execute("SELECT ref,job_metadata,facts,publication_state FROM scientific_jobs WHERE job_metadata!='{}'").fetchall()
        c.execute("DELETE FROM scientific_campaign_jobs")
        c.execute("DELETE FROM scientific_campaign_targets")
        c.execute("DELETE FROM scientific_campaigns")
    snapshots = {r["ref"]: {"job": json.loads(r["job_metadata"]), **json.loads(r["facts"]),
                             "publication_state": r["publication_state"]}
                 for r in selected}
    rows = {ref: dict(s["job"], campaign_name="", node_name=s["job"].get("workstation", "")) for ref, s in snapshots.items()}
    parents = {}
    marks = ",".join("?" for _ in rows) or "NULL"
    for edge in c.execute(f"SELECT child_ref,parent_ref FROM scientific_execution_edges WHERE child_ref IN ({marks}) OR parent_ref IN ({marks})", [*rows, *rows]):
        if edge["child_ref"] in rows and edge["parent_ref"] in rows:
            parents.setdefault(edge["child_ref"], set()).add(edge["parent_ref"])
    for ref, row in rows.items():
        # Never guess a parent from names or arbitrarily choose conflicting receipts.
        row.pop("parent_run_id", None)
        row.pop("campaign_id", None)
        if len(parents.get(ref, set())) == 1:
            parent = rows[next(iter(parents[ref]))]
            row.update(parent_run_id=parent["run_id"], parent_node_id=parent["node_id"], pipeline_parent_node_id=parent["node_id"])
    for campaign in group_runs(rows.values()):
        jobs = campaign["jobs"]
        root = snapshots[campaign["root"]["ref"]]["job"]
        kind = root["kind"]
        if kind == "unknown":
            activities = {j["kind"] for j in jobs} - {"analysis", "unknown"}
            kind = "design" if "design" in activities else next(iter(activities)) if len(activities) == 1 else "unknown"
        target_records = {}
        for j in jobs:
            target = snapshots[j["ref"]].get("target") or {}
            candidates = target.get("recorded_constructs", []) if target.get("mapping_status") == "multiple_constructs" else ([target] if target else [])
            for item in candidates:
                if item and (item.get("sequences_by_chain") or
                             any(component.get("sequence") for component in item.get("components") or [])):
                    normalized_target = _target_construct_record(item)
                    target_records[_target_construct_key(normalized_target)] = normalized_target
        targets = {encoded(target) for target in target_records.values()}
        counts = [dict(count, job_ref=j["ref"]) for j in jobs for count in snapshots[j["ref"]].get("counts", [])]
        sections = [{"id": out["id"], "name": out["name"], "job_refs": [j["ref"] for j in out["sources"]]}
                    for out in campaign["outputs"]]
        nonempty = [section for section in sections if any(snapshots[r].get("scientific_output_occurrences", 0) for r in section["job_refs"])]
        primary = next((s for s in nonempty if s["id"] == "results"), None)
        selected_section = primary or (nonempty[0] if len(nonempty) == 1 else None)
        sources = set(selected_section["job_refs"]) if selected_section else set()
        decisions = {}
        for ref in sources:
            for stage, stage_counts in snapshots[ref].get("acceptance_counts", {}).items():
                result = decisions.setdefault(stage, dict(passed=0, rejected=0, unresolved=0, unrecorded=0,
                    unit="candidate/model/pose occurrences", scope="source-authored aggregate; details remain on owner workstation"))
                for key in ("passed", "rejected", "unresolved", "unrecorded"):
                    result[key] += stage_counts[key]
        expected_links = sum(int(snapshots[j["ref"]].get("expected_occurrences", snapshots[j["ref"]].get("indexed_occurrences", 0)) or 0) for j in jobs)
        indexed_links = sum(int(snapshots[j["ref"]].get("indexed_occurrences", 0) or 0) for j in jobs)
        source_outputs = sum(int(snapshots[ref].get("source_output_occurrences", snapshots[ref].get("scientific_output_occurrences", 0)) or 0)
                             for ref in sources)
        source_summaries = [snapshots[ref].get("source_summary") or {} for ref in sources]
        all_job_summaries = [snapshots[j["ref"]].get("source_summary") or {} for j in jobs]
        kind_counts, stage_counts = {}, {}
        for source_summary in source_summaries:
            for key, value in (source_summary.get("output_kind_counts") or {}).items():
                kind_counts[key] = kind_counts.get(key, 0) + int(value or 0)
            for key, value in (source_summary.get("stage_occurrence_counts") or {}).items():
                stage_counts[key] = stage_counts.get(key, 0) + int(value or 0)
        complete_jobs = sum(snapshots[j["ref"]].get("publication_state") != "incomplete" for j in jobs)
        summary = {"targets": [json.loads(t) for t in sorted(targets)], "counts": counts,
                   "activity_kinds": sorted({str(j.get("kind") or "unknown") for j in jobs
                                              if j.get("kind") not in {"", "unknown"}}),
                   "jobs": len(jobs), "workstations": sorted({j.get("workstation") or j["node_id"] for j in jobs}),
                   "methods": sorted({j.get("method") for j in jobs if j.get("method")}),
                   "indexed_occurrences": indexed_links,
                   "retained_link_count": indexed_links,
                   "expected_retained_link_count": expected_links,
                   "source_output_occurrences": source_outputs,
                   "authoritative_output_occurrences": source_outputs,
                   "candidate_model_pose_occurrences": sum(int(s.get("candidate_model_pose_occurrences", 0) or 0) for s in source_summaries),
                   "generated_design_occurrences": sum(int(s.get("generated_design_occurrences", 0) or 0) for s in source_summaries),
                   "output_kind_counts": kind_counts, "stage_occurrence_counts": stage_counts,
                   "database_retained_candidate_model_occurrences": sum(int(s.get("database_retained_candidate_model_occurrences", 0) or 0) for s in all_job_summaries),
                   "database_retained_structure_references": sum(int(s.get("database_retained_structure_references", 0) or 0) for s in all_job_summaries),
                   "authoritative_source_retained_candidate_model_occurrences": sum(int(s.get("database_retained_candidate_model_occurrences", 0) or 0) for s in source_summaries),
                   "all_execution_job_retained_link_count": indexed_links,
                   "source_summary_complete": all(s.get("source_summary_complete") is True for s in source_summaries),
                   "source_job_refs": sorted(sources),
                   "publication": {"complete": complete_jobs == len(jobs), "complete_jobs": complete_jobs,
                       "incomplete_jobs": len(jobs) - complete_jobs,
                       "job_count": len(jobs), "expected_occurrences": expected_links,
                       "indexed_occurrences": indexed_links,
                       "content": "Database-retained references only; source Results are owner-local"},
                   "output_sections": sections, "source_selection": selected_section["id"] if selected_section else "unresolved / no scientific output",
                   "retained_representation_refs": sum(int((snapshots[ref].get("source_summary") or {}).get("database_retained_structure_references", 0) or 0) for ref in sources),
                   "acceptance_counts": decisions,
                   "mapping_warning_counts": {kind: sum(int((snapshots[ref].get("source_summary") or {}).get("mapping_warning_counts", {}).get(kind, 0) or 0) for ref in sources)
                                              for kind in sorted({kind for ref in sources for kind in (snapshots[ref].get("source_summary") or {}).get("mapping_warning_counts", {})})},
                   "mapping_warnings": []}
        target_key = next(iter(target_records)) if len(target_records) == 1 else ""
        updated = max((str(j.get("updated_at") or j.get("created_at") or "") for j in jobs), default="")
        c.execute("INSERT INTO scientific_campaigns VALUES (?,?,?,?,?,?,?,?)",
                  (campaign["id"], root.get("owner_id", ""), root.get("project_id", ""), kind,
                   root.get("name") or root["run_id"], target_key, updated, encoded(summary)))
        c.executemany("INSERT INTO scientific_campaign_jobs VALUES (?,?)", [(campaign["id"], j["ref"]) for j in jobs])
        for target_key, target in target_records.items():
            identity = _biological_identity_key(target)
            c.execute("INSERT OR IGNORE INTO scientific_target_constructs (target_key,biological_identity_key,description,target_json) VALUES (?,?,?,?)",
                      (target_key, identity, str(target.get("label") or ""), encoded(target)))
            c.execute("INSERT INTO scientific_campaign_targets VALUES (?,?,?)",
                      (campaign["id"], target_key, encoded({"relationship": "uses_exact_construct"})))


def _target_construct_record(target):
    """Keep only explicit target/construct and preparation facts; names are labels, not identity."""
    fields = ("biological_target_id", "target_molecule_id", "uniprot_id", "accession", "target_id",
              "biological_target", "label", "sequences_by_chain", "residue_boundaries",
              "target_residue_boundaries", "boundaries", "mutations", "chain_roles", "target_chain_roles",
              "target_component_roles", "components", "assembly_id", "structure_digest",
              "target_structure_digest", "model_id", "preparation")
    return {key: target[key] for key in fields if key in target and target[key] is not None}


def _target_construct_key(target):
    # Labels are useful UI context, but never create or merge target identity.
    identity = {key: value for key, value in target.items() if key not in {"label", "biological_target"}}
    return identifier("construct-", identity)


def _biological_identity_key(target):
    identity = {key: target[key] for key in ("biological_target_id", "target_molecule_id", "uniprot_id", "accession") if target.get(key)}
    # target_id is allowed only as an explicitly recorded stable identifier.
    if not identity and target.get("target_id"):
        identity = {"target_id": target["target_id"]}
    return identifier("biological-target-", identity) if identity else ""


def recorded_data_location(output):
    """Describe recorded storage without probing a remote filesystem."""
    states = {str(a.get("availability") or "unresolved") for a in output.get("artifacts", [])}
    if "owner_present" in states:
        owner_results = "present_on_owner_at_index_time"
    elif "owner_local_unverified" in states:
        owner_results = "owner_reference_recorded_live_reachability_not_checked"
    elif "missing" in states:
        owner_results = "missing_at_owner_at_index_time"
    elif states & {"inaccessible", "unresolved"}:
        owner_results = "source_availability_unresolved"
    else:
        owner_results = "no_results_artifact_reference_recorded"
    retained = output.get("retained", [])
    return {"owner_results": owner_results, "source_node_id": output.get("job_ref", "").split("/")[2] if output.get("job_ref", "").startswith("cluster://") else "",
            "database_representation_count": len(retained),
            "database_representation": "retained_reference" if retained else "none_mapped",
            "coordinate_locality": "not_checked_until_structure_is_opened",
            "remote_reachability": "not_checked_during_index_navigation"}


def _access(user, alias="c"):
    if not user.get("id"):
        raise PermissionError("A signed-in user is required.")
    if user.get("role") == "admin":
        return "1", []
    from mn_molecular_design_studio.auth import accessible_project_ids
    projects = list(accessible_project_ids(user))
    return (f"(({alias}.project_id='' AND {alias}.owner_id=?) OR {alias}.project_id IN "
            f"({','.join('?' for _ in projects) or 'NULL'}))", [user["id"], *projects])


def campaigns(database, user, *, kind="", project_id="", target="", query="", limit=25, offset=0):
    limit, offset = max(1, min(int(limit), 100)), max(0, int(offset))
    clause, params = _access(user)
    kinds = ([str(value) for value in kind if value] if isinstance(kind, (list, tuple, set))
             else [str(kind)] if kind else [])
    if kinds:
        # A logical campaign can contain several scientific activities (for
        # example design → independent folding). Keep c.kind as its root
        # identity, while allowing a campaign to be found by any recorded
        # stage in its compact summary.
        marks = ",".join("?" for _ in kinds)
        clause += (f" AND (c.kind IN ({marks}) OR EXISTS (SELECT 1 FROM json_each(c.summary,'$.activity_kinds') a "
                   f"WHERE a.value IN ({marks})))")
        params.extend([*kinds, *kinds])
    if project_id:
        clause += " AND c.project_id=?"
        params.append(project_id)
    if target:
        clause += " AND (c.target_key=? OR EXISTS (SELECT 1 FROM scientific_campaign_targets ct JOIN scientific_target_constructs t USING(target_key) WHERE ct.campaign_ref=c.ref AND (t.target_key=? OR t.biological_identity_key=?)))"
        params.extend([target, target, target])
    if query:
        clause += " AND c.name LIKE ?"
        params.append("%" + query + "%")
    with database.connection() as c:
        total = c.execute("SELECT count(*) FROM scientific_campaigns c WHERE " + clause, params).fetchone()[0]
        rows = c.execute("SELECT * FROM scientific_campaigns c WHERE " + clause +
                         " ORDER BY c.updated_at DESC,c.ref LIMIT ? OFFSET ?", [*params, limit, offset]).fetchall()
    items = [dict(row, summary=json.loads(row["summary"])) for row in rows]
    return {"items": items, "total": total, "next_offset": offset + limit if offset + limit < total else None}


def campaign_detail(database, user, ref, *, limit=100, offset=0, job_limit=100, job_offset=0,
                    dependency_limit=100, dependency_offset=0):
    clause, params = _access(user)
    with database.connection() as c:
        campaign = c.execute("SELECT c.* FROM scientific_campaigns c WHERE c.ref=? AND " + clause, [ref, *params]).fetchone()
        if not campaign:
            raise PermissionError("Scientific campaign not found or inaccessible.")
        job_total = c.execute("SELECT count(*) FROM scientific_jobs j JOIN scientific_campaign_jobs m ON m.job_ref=j.ref WHERE m.campaign_ref=?", (ref,)).fetchone()[0]
        incomplete_job_count = c.execute("SELECT count(*) FROM scientific_jobs j JOIN scientific_campaign_jobs m ON m.job_ref=j.ref WHERE m.campaign_ref=? AND j.publication_state='incomplete'", (ref,)).fetchone()[0]
        jobs = c.execute("SELECT j.job_metadata,j.facts,j.publication_state,j.publication_revision,j.publication_pages FROM scientific_jobs j JOIN scientific_campaign_jobs m ON m.job_ref=j.ref WHERE m.campaign_ref=? ORDER BY j.ref LIMIT ? OFFSET ?",
                         (ref, max(1, min(int(job_limit), 500)), max(0, int(job_offset)))).fetchall()
        total = c.execute("SELECT count(*) FROM scientific_outputs o JOIN scientific_jobs j ON j.ref=o.job_ref JOIN scientific_campaign_jobs m ON m.job_ref=o.job_ref WHERE m.campaign_ref=? AND j.publication_state!='incomplete'", (ref,)).fetchone()[0]
        outputs = c.execute("SELECT o.payload FROM scientific_outputs o JOIN scientific_jobs j ON j.ref=o.job_ref JOIN scientific_campaign_jobs m ON m.job_ref=o.job_ref WHERE m.campaign_ref=? AND j.publication_state!='incomplete' ORDER BY o.ref LIMIT ? OFFSET ?",
                            (ref, max(1, min(int(limit), 500)), max(0, int(offset)))).fetchall()
        dependency_total = c.execute("SELECT count(*) FROM scientific_dependencies d JOIN scientific_campaign_jobs m ON m.job_ref=d.job_ref WHERE m.campaign_ref=?", (ref,)).fetchone()[0]
        edges = c.execute("SELECT d.* FROM scientific_dependencies d JOIN scientific_campaign_jobs m ON m.job_ref=d.job_ref WHERE m.campaign_ref=? ORDER BY d.job_ref,d.source_ref,d.target_ref LIMIT ? OFFSET ?",
                          (ref, max(1, min(int(dependency_limit), 500)), max(0, int(dependency_offset)))).fetchall()
    job_status = [{**json.loads(row["job_metadata"]), "publication_state": row["publication_state"],
                   "publication_revision": row["publication_revision"],
                   "publication_pages": row["publication_pages"],
                   "publication": json.loads(row["facts"]).get("publication_state", "legacy")}
                  for row in jobs]
    return {"campaign": dict(campaign, summary=json.loads(campaign["summary"])),
            "jobs": job_status,
            "job_total": job_total, "incomplete_job_count": incomplete_job_count,
            "next_job_offset": job_offset + len(job_status) if job_offset + len(job_status) < job_total else None,
            "incomplete_jobs": [j["ref"] for j in job_status if j["publication_state"] == "incomplete"],
            "outputs": [json.loads(row[0]) for row in outputs], "output_total": total,
            "dependencies": [dict(row, evidence=json.loads(row["evidence"])) for row in edges],
            "dependency_total": dependency_total,
            "next_dependency_offset": dependency_offset + len(edges) if dependency_offset + len(edges) < dependency_total else None}


def retained_models(database, user, campaign_ref, *, stage="final", limit=100, offset=0, include_execution_refs=False):
    """Exact retained representations, reauthorized and resolved to local peer IDs."""
    from mn_molecular_design_studio.molecular_database import _id
    clause, parameters = _access(user)
    limit, offset = max(1, min(int(limit), 500)), max(0, int(offset))
    with database.connection() as c:
        campaign = c.execute("SELECT summary FROM scientific_campaigns c WHERE c.ref=? AND " + clause, [campaign_ref, *parameters]).fetchone()
        if not campaign:
            raise PermissionError("Scientific campaign not found or inaccessible.")
        source_refs = json.loads(campaign[0])["source_job_refs"]
        where, params = "m.campaign_ref=? AND o.kind!='evaluation' AND j.publication_state!='incomplete'", [campaign_ref]
        if not include_execution_refs:
            where += " AND o.job_ref IN (SELECT value FROM json_each(?))"
            params.append(encoded(source_refs))
        if stage != "all":
            where += (" AND EXISTS(SELECT 1 FROM json_each(o.payload,'$.acceptance') d WHERE json_extract(d.value,'$.stage_id')=? AND json_extract(d.value,'$.cohort_status')='passed')"
                      " AND NOT EXISTS(SELECT 1 FROM json_each(o.payload,'$.acceptance') d WHERE json_extract(d.value,'$.stage_id')=? AND coalesce(json_extract(d.value,'$.cohort_status'),'unresolved')!='passed')")
            params += [stage, stage]
        query = " FROM scientific_outputs o JOIN scientific_jobs j ON j.ref=o.job_ref JOIN scientific_campaign_jobs m ON m.job_ref=o.job_ref WHERE " + where
        total = c.execute("SELECT count(*)" + query, params).fetchone()[0]
        outputs = [json.loads(row[0]) for row in c.execute("SELECT o.payload" + query + " ORDER BY o.ref LIMIT ? OFFSET ?", [*params, limit, offset])]
    models = {}
    with database.connection() as c:
        for output in outputs:
            decisions = [d for d in output.get("acceptance", []) if d.get("stage_id") == stage]
            if stage != "all" and (not decisions or {d.get("cohort_status") for d in decisions} != {"passed"}):
                continue
            for retained in output.get("retained", []):
                source_id = retained["owner_structure_id"]
                peer_id = _id("peer-structure-", retained["source_node"] + ":" + source_id)
                matches = c.execute("SELECT s.*,g.owner_id,g.project_id FROM molecule_structures s JOIN molecule_groups g ON s.group_id=g.id WHERE s.id IN (?,?)", (source_id, peer_id)).fetchall()
                for row in matches:
                    context = json.loads(row["context"])
                    # Local and source IDs may coincide in a staging copy. The
                    # source context and digest still have to agree.
                    if context.get("_replica_node") and context["_replica_node"] != retained["source_node"]:
                        continue
                    if row["molecule_id"] != retained["molecule_id"] or not database._can_read(row, user):
                        continue
                    digest = context.get("portable_sha256") or context.get("content_sha256") or ""
                    if retained.get("sha256") and digest != retained["sha256"]:
                        continue
                    evaluations = c.execute("SELECT descriptor FROM molecule_evaluations WHERE molecule_id=? AND group_id=? AND "
                        "(json_extract(descriptor,'$.representation_id')=? OR "
                        "(?!='' AND json_extract(descriptor,'$.structure_digest')=? AND run_id=? "
                        "AND json_extract(descriptor,'$.candidate_id') IN (?,?) AND json_extract(descriptor,'$.model_id')=?))",
                        (row["molecule_id"], row["group_id"], row["id"], digest, digest,
                         output["ref"].split("#", 1)[0].rsplit("/", 1)[-1], output.get("candidate_id", ""),
                         output.get("native_id", ""), output.get("model_id", ""))).fetchall()
                    location = recorded_data_location(output)
                    location["database_representation"] = "authorized_structure_reference_matches_recorded_identity_and_digest"
                    location["coordinate_bytes"] = "not_loaded_until_explicit_open"
                    models[row["id"]] = {"structure_id": row["id"], "molecule_id": row["molecule_id"],
                        "group_id": row["group_id"], "sha256": digest, "output_ref": output["ref"],
                        "candidate_id": output.get("candidate_id"), "model_id": output.get("model_id"),
                        "acceptance": decisions, "evaluations": [json.loads(e[0]) for e in evaluations],
                        "data_location": location}
    return {"items": list(models.values()), "output_total": total,
            "next_offset": offset + limit if offset + limit < total else None}


def analysis_placement(artifacts, *, local_node):
    """No transfer or dispatch: describe data locality to the existing scheduler."""
    artifacts = [a for a in artifacts if a.get("required", True)]
    owners = {a.get("node_id") for a in artifacts if a.get("storage") == "results"}
    if any(a.get("availability") in {"missing", "unresolved", "inaccessible"} for a in artifacts):
        return {"status": "blocked", "reason": "A required scientific artifact is unavailable."}
    if len(owners) > 1:
        return {"status": "requires_explicit_staging", "reason": "Results inputs have multiple owning workstations."}
    return {"status": "owner_local" if owners else "retained_inputs", "node_id": next(iter(owners), local_node)}
