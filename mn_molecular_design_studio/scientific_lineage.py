"""Versioned, rebuildable scientific lineage projection in the existing Database.

Job ownership is authoritative. Snapshots contain metadata, never file bytes or
new canonical molecular identities. Existing job replication carries snapshots.
"""
from __future__ import annotations

import hashlib
import json
import re

VERSION = "studio.scientific-lineage.v1"
KINDS = {"design", "folding", "docking", "dynamics", "analysis", "unknown"}
MAX_BYTES = 1_500_000  # Leave room within the existing 2 MB replication envelope.


def encoded(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def acceptance_counts(outputs):
    rows = [o for o in outputs if o.get("kind") not in {"evaluation", "retained_representation"}]
    stages = sorted({d["stage_id"] for o in rows for d in o.get("acceptance", []) if d.get("stage_id")})
    result = {}
    for stage in stages:
        counts = dict(passed=0, rejected=0, unresolved=0, unrecorded=0)
        for row in rows:
            states = {d.get("cohort_status") for d in row.get("acceptance", []) if d.get("stage_id") == stage}
            status = "unrecorded" if not states else next(iter(states)) if len(states) == 1 else "unresolved"
            counts[status if status in counts else "unresolved"] += 1
        result[stage] = dict(counts, unit="candidate/model occurrences", scope="authoritative saved outputs")
    return result


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
    """)
    columns = {r[1] for r in c.execute("PRAGMA table_info(scientific_jobs)")}
    for name in ("job_metadata", "facts"):
        if name not in columns:
            c.execute(f"ALTER TABLE scientific_jobs ADD COLUMN {name} TEXT NOT NULL DEFAULT '{{}}'")


def validate_snapshot(snapshot, author, envelope=None):
    if not isinstance(snapshot, dict) or snapshot.get("schema") != VERSION:
        raise ValueError("Unsupported scientific lineage snapshot.")
    job = snapshot.get("job") or {}
    if not isinstance(job, dict):
        raise ValueError("Scientific job metadata must be an object.")
    if job.get("node_id") != author or job.get("ref") != job_reference(author, job.get("task_group", ""), job.get("run_id", "")):
        raise PermissionError("Only the scientific source workstation may publish its job.")
    if job.get("kind") not in KINDS or len(encoded(snapshot).encode()) > MAX_BYTES:
        raise ValueError("Invalid or oversized scientific lineage snapshot.")
    if envelope is not None and any(str(job.get(k) or "") != str(envelope.get(k) or "")
                                    for k in ("owner_id", "project_id", "node_id", "run_id", "task_group")):
        raise PermissionError("Scientific lineage ownership differs from its job envelope.")
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


def put_snapshot(database, snapshot, *, author, envelope=None, rebuild=True):
    job = validate_snapshot(snapshot, author, envelope)
    outputs = snapshot.get("outputs", [])
    facts = {"target": snapshot.get("target") or {}, "counts": snapshot.get("counts") or [],
             "warnings": snapshot.get("warnings") or [], "indexed_occurrences": len(outputs),
             "scientific_output_occurrences": len([o for o in outputs if o.get("kind") not in {"evaluation", "retained_representation"}]),
             "acceptance_counts": acceptance_counts(outputs),
             "retained_refs": sorted({(r.get("source_node", ""), r.get("owner_structure_id", "")) for o in outputs for r in o.get("retained", [])})}
    with database.connection() as c:
        previous = c.execute("SELECT owner_id,project_id FROM scientific_jobs WHERE ref=?", (job["ref"],)).fetchone()
        if previous and tuple(previous) != (job.get("owner_id", ""), job.get("project_id", "")):
            raise PermissionError("A lineage update cannot transfer project or user ownership.")
        c.execute("INSERT INTO scientific_jobs (ref,node_id,owner_id,project_id,task_group,run_id,cluster_id,snapshot,job_metadata,facts) VALUES (?,?,?,?,?,?,?,?,?,?) ON CONFLICT(ref) DO UPDATE SET "
                  "cluster_id=excluded.cluster_id,snapshot=excluded.snapshot,job_metadata=excluded.job_metadata,facts=excluded.facts",
                  (job["ref"], author, job.get("owner_id", ""), job.get("project_id", ""),
                   job["task_group"], job["run_id"], job.get("cluster_job_id", ""), encoded(snapshot), encoded(job), encoded(facts)))
        c.execute("DELETE FROM scientific_outputs WHERE job_ref=?", (job["ref"],))
        c.execute("DELETE FROM scientific_dependencies WHERE job_ref=?", (job["ref"],))
        for row in snapshot.get("outputs", []):
            c.execute("INSERT INTO scientific_outputs VALUES (?,?,?,?,?,?,?)",
                      (row["ref"], job["ref"], row["kind"], str(row.get("candidate_id") or ""),
                       str(row.get("model_id") or ""), str(row.get("stage") or ""), encoded(row)))
        for edge in snapshot.get("dependencies", []):
            c.execute("INSERT OR IGNORE INTO scientific_dependencies VALUES (?,?,?,?,?,?)",
                      (job["ref"], edge["source_ref"], edge["target_ref"], edge["role"],
                       edge.get("status", "unresolved"), encoded(edge.get("evidence") or {})))
        if rebuild:
            rebuild_campaigns(c)


def rebuild_campaigns(c):
    """Metadata-only write-time grouping; never executed by a page query."""
    from mn_molecular_design_studio.logical_runs import group_runs
    snapshots = {r["ref"]: {"job": json.loads(r["job_metadata"]), **json.loads(r["facts"])}
                 for r in c.execute("SELECT ref,job_metadata,facts FROM scientific_jobs") if r["job_metadata"] != "{}"}
    rows = {ref: dict(s["job"], campaign_name="", node_name=s["job"].get("workstation", "")) for ref, s in snapshots.items()}
    alias_candidates = {}
    for ref, row in rows.items():
        for alias in (["cluster-job:" + row["cluster_job_id"]] if row.get("cluster_job_id") else []) + [f"node-run:{row['node_id']}:{row['run_id']}"]:
            alias_candidates.setdefault(alias, set()).add(ref)
    aliases = {alias: next(iter(refs)) for alias, refs in alias_candidates.items() if len(refs) == 1}
    parents = {}
    for edge in c.execute("SELECT * FROM scientific_dependencies WHERE role='execution_child_of' AND status='recorded'"):
        child = aliases.get(edge["source_ref"], edge["source_ref"])
        parent = aliases.get(edge["target_ref"], edge["target_ref"])
        if child in rows and parent in rows and child != parent:
            a, b = rows[child], rows[parent]
            if (a.get("owner_id"), a.get("project_id")) == (b.get("owner_id"), b.get("project_id")):
                parents.setdefault(child, set()).add(parent)
    for ref, row in rows.items():
        # Never guess a parent from names or arbitrarily choose conflicting receipts.
        row.pop("parent_run_id", None)
        row.pop("campaign_id", None)
        if len(parents.get(ref, set())) == 1:
            parent = rows[next(iter(parents[ref]))]
            row.update(parent_run_id=parent["run_id"], parent_node_id=parent["node_id"], pipeline_parent_node_id=parent["node_id"])
    c.execute("DELETE FROM scientific_campaign_jobs")
    c.execute("DELETE FROM scientific_campaigns")
    for campaign in group_runs(rows.values()):
        jobs = campaign["jobs"]
        root = snapshots[campaign["root"]["ref"]]["job"]
        kind = root["kind"]
        if kind == "unknown":
            activities = {j["kind"] for j in jobs} - {"analysis", "unknown"}
            kind = "design" if "design" in activities else next(iter(activities)) if len(activities) == 1 else "unknown"
        targets = {encoded(snapshots[j["ref"]].get("target") or {}) for j in jobs
                   if snapshots[j["ref"]].get("target")}
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
                    unit="candidate/model occurrences", scope="authoritative saved outputs"))
                for key in ("passed", "rejected", "unresolved", "unrecorded"):
                    result[key] += stage_counts[key]
        summary = {"targets": [json.loads(t) for t in sorted(targets)], "counts": counts,
                   "jobs": len(jobs), "workstations": sorted({j.get("workstation") or j["node_id"] for j in jobs}),
                   "methods": sorted({j.get("method") for j in jobs if j.get("method")}),
                   "indexed_occurrences": sum(snapshots[j["ref"]].get("indexed_occurrences", 0) for j in jobs),
                   "authoritative_output_occurrences": sum(snapshots[r].get("scientific_output_occurrences", 0) for r in sources), "source_job_refs": sorted(sources),
                   "output_sections": sections, "source_selection": selected_section["id"] if selected_section else "unresolved / no scientific output",
                   "retained_representation_refs": len({tuple(r) for j in jobs for r in snapshots[j["ref"]].get("retained_refs", [])}),
                   "acceptance_counts": decisions,
                   "mapping_warnings": [dict(job_ref=j["ref"], warning=w) for j in jobs
                                        for w in snapshots[j["ref"]].get("warnings", [])]}
        exact_target = json.loads(next(iter(targets))) if len(targets) == 1 else {}
        target_key = identifier("construct-", exact_target) if exact_target and exact_target.get("mapping_status") != "multiple_constructs" else ""
        updated = max((str(j.get("updated_at") or j.get("created_at") or "") for j in jobs), default="")
        c.execute("INSERT INTO scientific_campaigns VALUES (?,?,?,?,?,?,?,?)",
                  (campaign["id"], root.get("owner_id", ""), root.get("project_id", ""), kind,
                   root.get("name") or root["run_id"], target_key, updated, encoded(summary)))
        c.executemany("INSERT INTO scientific_campaign_jobs VALUES (?,?)", [(campaign["id"], j["ref"]) for j in jobs])


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
    for field, value in (("kind", kind), ("project_id", project_id), ("target_key", target)):
        if value:
            clause += f" AND c.{field}=?"
            params.append(value)
    if query:
        clause += " AND c.name LIKE ?"
        params.append("%" + query + "%")
    with database.connection() as c:
        total = c.execute("SELECT count(*) FROM scientific_campaigns c WHERE " + clause, params).fetchone()[0]
        rows = c.execute("SELECT * FROM scientific_campaigns c WHERE " + clause +
                         " ORDER BY c.updated_at DESC,c.ref LIMIT ? OFFSET ?", [*params, limit, offset]).fetchall()
    items = [dict(row, summary=json.loads(row["summary"])) for row in rows]
    return {"items": items, "total": total, "next_offset": offset + limit if offset + limit < total else None}


def campaign_detail(database, user, ref, *, limit=100, offset=0):
    clause, params = _access(user)
    with database.connection() as c:
        campaign = c.execute("SELECT c.* FROM scientific_campaigns c WHERE c.ref=? AND " + clause, [ref, *params]).fetchone()
        if not campaign:
            raise PermissionError("Scientific campaign not found or inaccessible.")
        jobs = c.execute("SELECT j.job_metadata FROM scientific_jobs j JOIN scientific_campaign_jobs m ON m.job_ref=j.ref WHERE m.campaign_ref=?", (ref,)).fetchall()
        total = c.execute("SELECT count(*) FROM scientific_outputs o JOIN scientific_campaign_jobs m ON m.job_ref=o.job_ref WHERE m.campaign_ref=?", (ref,)).fetchone()[0]
        outputs = c.execute("SELECT o.payload FROM scientific_outputs o JOIN scientific_campaign_jobs m ON m.job_ref=o.job_ref WHERE m.campaign_ref=? ORDER BY o.ref LIMIT ? OFFSET ?",
                            (ref, max(1, min(int(limit), 500)), max(0, int(offset)))).fetchall()
        edges = c.execute("SELECT d.* FROM scientific_dependencies d JOIN scientific_campaign_jobs m ON m.job_ref=d.job_ref WHERE m.campaign_ref=?", (ref,)).fetchall()
    return {"campaign": dict(campaign, summary=json.loads(campaign["summary"])),
            "jobs": [json.loads(row["job_metadata"]) for row in jobs],
            "outputs": [json.loads(row[0]) for row in outputs], "output_total": total,
            "dependencies": [dict(row, evidence=json.loads(row["evidence"])) for row in edges]}


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
        where, params = "m.campaign_ref=? AND o.kind!='evaluation'", [campaign_ref]
        if not include_execution_refs:
            where += " AND o.job_ref IN (SELECT value FROM json_each(?))"
            params.append(encoded(source_refs))
        if stage != "all":
            where += (" AND EXISTS(SELECT 1 FROM json_each(o.payload,'$.acceptance') d WHERE json_extract(d.value,'$.stage_id')=? AND json_extract(d.value,'$.cohort_status')='passed')"
                      " AND NOT EXISTS(SELECT 1 FROM json_each(o.payload,'$.acceptance') d WHERE json_extract(d.value,'$.stage_id')=? AND coalesce(json_extract(d.value,'$.cohort_status'),'unresolved')!='passed')")
            params += [stage, stage]
        query = " FROM scientific_outputs o JOIN scientific_campaign_jobs m ON m.job_ref=o.job_ref WHERE " + where
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
                    models[row["id"]] = {"structure_id": row["id"], "molecule_id": row["molecule_id"],
                        "group_id": row["group_id"], "sha256": digest, "output_ref": output["ref"],
                        "candidate_id": output.get("candidate_id"), "model_id": output.get("model_id"),
                        "acceptance": decisions, "evaluations": [json.loads(e[0]) for e in evaluations]}
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
