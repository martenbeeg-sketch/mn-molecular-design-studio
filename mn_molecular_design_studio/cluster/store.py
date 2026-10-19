"""Coordinator metadata only. Run data is never resolved on this machine."""
from __future__ import annotations

import hashlib
import hmac
import json
import secrets
import sqlite3
import time
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

from . import config

TERMINAL = {"completed", "failed", "cancelled", "succeeded"}


def stamp(epoch=None):
    return datetime.fromtimestamp(time.time() if epoch is None else epoch, timezone.utc).isoformat(timespec="seconds")


def encoded(value):
    return json.dumps(value, allow_nan=False, separators=(",", ":"))


def visible(row, user):
    return (user.get("role") == "admin" or row.get("owner_id") == user.get("id")
            or row.get("project_id") in set(user.get("project_ids") or []))


class Store:
    def __init__(self, root: Path | None = None, *, clock=time.time):
        self.root = Path(root) if root else config.home()
        self.root.mkdir(parents=True, exist_ok=True, mode=0o700)
        self.path = self.root / "coordinator.sqlite3"
        self.clock = clock
        with self.connect() as db:
            db.executescript("""
                PRAGMA journal_mode=WAL;
                CREATE TABLE IF NOT EXISTS nodes (
                    id TEXT PRIMARY KEY, name TEXT NOT NULL, token_hash TEXT NOT NULL,
                    disabled INTEGER NOT NULL DEFAULT 0, last_seen REAL,
                    snapshot TEXT NOT NULL DEFAULT '{}');
                CREATE TABLE IF NOT EXISTS jobs (
                    id TEXT PRIMARY KEY, node_id TEXT, local_ref TEXT,
                    owner_id TEXT NOT NULL, project_id TEXT NOT NULL,
                    payload TEXT NOT NULL, created_at TEXT NOT NULL,
                    UNIQUE(node_id,local_ref));
                CREATE TABLE IF NOT EXISTS submissions (
                    job_id TEXT PRIMARY KEY, owner_id TEXT NOT NULL,
                    idempotency_key TEXT NOT NULL, fingerprint TEXT NOT NULL,
                    request TEXT NOT NULL, requirements TEXT NOT NULL,
                    origin_node TEXT NOT NULL, placement TEXT NOT NULL,
                    affinity_node TEXT, command_id TEXT, cancelled INTEGER NOT NULL DEFAULT 0,
                    UNIQUE(owner_id,idempotency_key));
                CREATE TABLE IF NOT EXISTS commands (
                    id TEXT PRIMARY KEY, node_id TEXT NOT NULL,
                    job_id TEXT NOT NULL, kind TEXT NOT NULL, payload TEXT NOT NULL,
                    owner_id TEXT NOT NULL, project_id TEXT NOT NULL,
                    status TEXT NOT NULL DEFAULT 'queued', lease_token TEXT,
                    lease_until REAL, response TEXT, error TEXT,
                    created_at TEXT NOT NULL);
                CREATE INDEX IF NOT EXISTS cluster_commands_queue ON commands(node_id,status,created_at);
            """)
            if "cancelled" not in {row[1] for row in db.execute("PRAGMA table_info(submissions)")}:
                db.execute("ALTER TABLE submissions ADD COLUMN cancelled INTEGER NOT NULL DEFAULT 0")
        self.path.chmod(0o600)

    @contextmanager
    def connect(self):
        db = sqlite3.connect(self.path, timeout=20)
        db.row_factory = sqlite3.Row
        db.execute("PRAGMA busy_timeout=20000")
        try:
            with db:
                yield db
        finally:
            db.close()

    def enroll(self, name, node_id=None):
        node_id = str(uuid.UUID(node_id)) if node_id else str(uuid.uuid4())
        token = secrets.token_urlsafe(40)
        if not isinstance(name, str):
            raise ValueError("Enter a workstation name.")
        name = name.strip()
        if not name or len(name) > 120:
            raise ValueError("Enter a workstation name of at most 120 characters.")
        with self.connect() as db:
            if db.execute("SELECT 1 FROM nodes WHERE id=?", (node_id,)).fetchone():
                raise ValueError("This workstation is already enrolled; rotate its token instead.")
            db.execute("INSERT INTO nodes(id,name,token_hash) VALUES (?,?,?)",
                       (node_id, name, hashlib.sha256(token.encode()).hexdigest()))
        return {"node_id": node_id, "name": name, "node_token": token}

    def authenticate_node(self, token):
        digest = hashlib.sha256(str(token).encode()).hexdigest()
        with self.connect() as db:
            row = db.execute("SELECT * FROM nodes WHERE token_hash=? AND disabled=0", (digest,)).fetchone()
        if not row or not hmac.compare_digest(row["token_hash"], digest):
            raise PermissionError("Unknown or revoked workstation credential.")
        return dict(row)

    def nodes(self):
        with self.connect() as db:
            rows = db.execute("SELECT * FROM nodes ORDER BY name,id").fetchall()
        result = []
        for row in rows:
            item = dict(row)
            item.pop("token_hash", None)
            snapshot = json.loads(item.pop("snapshot"))
            seen = item.pop("last_seen")
            age = self.clock() - seen if seen is not None else float("inf")
            item.update(snapshot=snapshot, last_seen_at=stamp(seen) if seen is not None else None,
                        availability="disabled" if item["disabled"] else "online" if age < 45 else "unreachable" if age < 120 else "offline")
            result.append(item)
        return result

    def heartbeat(self, node_id, snapshot, jobs, known_references=None):
        if not isinstance(snapshot, dict) or not isinstance(jobs, list) or len(jobs) > 1000:
            raise ValueError("Invalid workstation heartbeat.")
        with self.connect() as db:
            db.execute("BEGIN IMMEDIATE")
            db.execute("UPDATE nodes SET snapshot=?,last_seen=? WHERE id=? AND disabled=0",
                       (encoded(snapshot), self.clock(), node_id))
            for raw in jobs:
                row = dict(raw)
                group, run = str(row.get("task_group") or ""), str(row.get("run_id") or "")
                if not group or not run or "/" in group + run or ".." in {group, run}:
                    raise ValueError("Invalid local job reference.")
                local_ref = group + "/" + run
                existing = db.execute("SELECT * FROM jobs WHERE node_id=? AND local_ref=?", (node_id, local_ref)).fetchone()
                jid = str(row.get("cluster_job_id") or config.job_id(node_id, local_ref))
                superseded_payload = {}
                if existing and existing["id"] != jid:
                    # Older agents could publish the deterministic machine-local
                    # ID before recovering the real cluster submission ID. If the
                    # original submission record is present, promote it and retire
                    # the fallback row so the shared catalog has one identity.
                    if existing["owner_id"] != str(row.get("owner_id") or "") or existing["project_id"] != str(row.get("project_id") or ""):
                        raise PermissionError("A workstation cannot change job ownership or project while reconciling a job identity.")
                    canonical = db.execute("SELECT * FROM jobs WHERE id=?", (jid,)).fetchone()
                    if canonical:
                        if (canonical["node_id"] not in {None, node_id}
                                or canonical["local_ref"] not in {None, local_ref}
                                or canonical["owner_id"] != str(row.get("owner_id") or "")
                                or canonical["project_id"] != str(row.get("project_id") or "")):
                            raise PermissionError("The reported cluster submission ID belongs to another job.")
                        superseded_payload = json.loads(existing["payload"])
                        db.execute("DELETE FROM jobs WHERE id=?", (existing["id"],))
                    else:
                        db.execute("UPDATE jobs SET id=? WHERE id=?", (jid, existing["id"]))
                collision = db.execute("SELECT * FROM jobs WHERE id=?", (jid,)).fetchone()
                if collision and (collision["node_id"] != node_id or collision["owner_id"] != str(row.get("owner_id") or "") or collision["project_id"] != str(row.get("project_id") or "")):
                    raise PermissionError("A workstation cannot change job ownership or adopt another workstation's run.")
                if collision:
                    previous = json.loads(collision["payload"])
                    row = previous | superseded_payload | row
                    if previous.get("cancel_requested") or superseded_payload.get("cancel_requested"):
                        row["cancel_requested"] = True
                row.update(cluster_job_id=jid, node_id=node_id, local_reference=local_ref,
                           run_dir=config.artifact_uri(node_id, group, run), location=row.get("location") or raw.get("run_dir"),
                           data_state="present", wait_reason="")
                db.execute("INSERT INTO jobs VALUES (?,?,?,?,?,?,?) ON CONFLICT(id) DO UPDATE SET local_ref=excluded.local_ref,payload=excluded.payload",
                           (jid, node_id, local_ref, str(row.get("owner_id") or ""), str(row.get("project_id") or ""), encoded(row), row.get("created_at") or stamp(self.clock())))
            if known_references is not None and snapshot.get("storage_available", True):
                known = set(known_references)
                for old in db.execute("SELECT id,local_ref,payload FROM jobs WHERE node_id=? AND local_ref IS NOT NULL", (node_id,)).fetchall():
                    data = json.loads(old["payload"])
                    state = "present" if old["local_ref"] in known else "not_found"
                    if data.get("data_state") != state:
                        data["data_state"] = state
                        db.execute("UPDATE jobs SET payload=? WHERE id=?", (encoded(data), old["id"]))
        self.dispatch()

    def _job(self, db, jid):
        row = db.execute("SELECT * FROM jobs WHERE id=?", (str(jid).removeprefix("cluster:"),)).fetchone()
        if not row:
            raise LookupError("Job not found.")
        return dict(row)

    def job(self, jid, user):
        with self.connect() as db:
            row = self._job(db, jid)
        if not visible(row, user):
            raise LookupError("Job not found.")
        item = json.loads(row["payload"])
        node = next((n for n in self.nodes() if n["id"] == row["node_id"]), None)
        item.update(node_name=node["name"] if node else "Not assigned", availability=node["availability"] if node else "pending",
                    last_seen_at=node["last_seen_at"] if node else None)
        if node and node["availability"] == "online" and not node["snapshot"].get("storage_available", True):
            item["availability"] = "storage_unavailable"
        elif node and node["availability"] == "online" and item.get("data_state") == "not_found":
            item["availability"] = "files_unavailable"
        return item

    def search(self, user, *, node_id=None, project_id=None, query="", statuses=None, tools=None, limit=100, offset=0):
        if type(limit) is not int or not 1 <= limit <= 1000 or type(offset) is not int or offset < 0:
            raise ValueError("Invalid result page.")
        with self.connect() as db:
            rows = db.execute("SELECT * FROM jobs ORDER BY created_at DESC,id DESC").fetchall()
        selected = []
        for raw in rows:
            row = dict(raw)
            if not visible(row, user) or node_id and row["node_id"] != node_id or project_id and row["project_id"] != project_id:
                continue
            item = json.loads(row["payload"])
            if statuses and item.get("status") not in statuses or tools and item.get("tool") not in tools:
                continue
            if query and str(query).casefold() not in " ".join(str(item.get(k) or "") for k in ("campaign_name", "pipeline_name", "run_id", "tool", "job_code")).casefold():
                continue
            selected.append(row["id"])
        items = [self.job(j, user) for j in selected[offset:offset + limit]]
        # Detailed manifests and request graphs are fetched when opening a job.
        return {"items": [{k: v for k, v in item.items() if k not in {"artifacts", "request_summary"}} for item in items], "total": len(selected)}

    def submit(self, request, requirements, user, *, origin_node, placement="auto", affinity_node=None, idempotency_key=None):
        if placement not in {"auto", "local"} and placement not in {n["id"] for n in self.nodes()}:
            raise ValueError("Select a registered workstation.")
        if placement == "local":
            placement = origin_node
        if affinity_node and placement != "auto" and placement != affinity_node:
            raise ValueError("This job must run on the workstation that owns its input data.")
        key = str(idempotency_key or uuid.uuid4())
        fingerprint = hashlib.sha256(json.dumps([request, requirements, origin_node, placement, affinity_node], sort_keys=True, allow_nan=False).encode()).hexdigest()
        with self.connect() as db:
            db.execute("BEGIN IMMEDIATE")
            previous = db.execute("SELECT * FROM submissions WHERE owner_id=? AND idempotency_key=?", (user["id"], key)).fetchone()
            if previous:
                if previous["fingerprint"] != fingerprint:
                    raise ValueError("This submission ID was already used for different inputs.")
                jid = previous["job_id"]
            else:
                jid = str(uuid.uuid4())
                payload = request.get("payload", {})
                params = payload.get("parameters", {})
                workflow = str(payload.get("workflow") or "visual_pipeline")
                graph = payload.get("graph", {})
                name = params.get("job_name") or params.get("campaign_name") or graph.get("name") or params.get("pipeline_name") or workflow or "Pipeline"
                pipeline_name = ""
                if workflow == "visual_pipeline":
                    pipeline_name = str(params.get("pipeline_name") or graph.get("name") or name)
                    tool = "visual_pipeline"
                elif workflow == "refolding.optimized_prediction":
                    engine = str(params.get("engine") or "")
                    tool = "optimized_prediction_" + engine if engine else workflow
                else:
                    try:
                        from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import WORKFLOWS
                        tool = WORKFLOWS[workflow].tool
                    except (ImportError, KeyError):
                        tool = workflow
                preferred = affinity_node or (placement if placement != "auto" else None)
                row = {"cluster_job_id": jid, "node_id": preferred, "owner_id": user["id"], "project_id": request["project_id"],
                       "status": "waiting_for_workstation", "campaign_name": name, "created_at": stamp(self.clock()),
                       "tool": tool, "pipeline_name": pipeline_name, "job_code": jid[:5].upper(),
                       "request_summary": request["payload"], "wait_reason": "Waiting for a compatible workstation."}
                db.execute("INSERT INTO jobs VALUES (?,?,?,?,?,?,?)", (jid, None, None, user["id"], request["project_id"], encoded(row), row["created_at"]))
                db.execute("INSERT INTO submissions(job_id,owner_id,idempotency_key,fingerprint,request,requirements,origin_node,placement,affinity_node,command_id) VALUES (?,?,?,?,?,?,?,?,?,NULL)",
                           (jid, user["id"], key, fingerprint, encoded(request), encoded(requirements), origin_node, placement, affinity_node))
        self.dispatch()
        return {"job": self.job(jid, user)}

    @staticmethod
    def compatible(node, needs):
        snapshot = node["snapshot"]
        if node["availability"] != "online" or not snapshot.get("accept_jobs") or not snapshot.get("scheduler_online"):
            return False
        if needs.get("code_revision") and snapshot.get("code_revision") != needs["code_revision"]:
            return False
        for tag, digest in needs.get("images", {}).items():
            present = snapshot.get("images", {}).get(tag)
            if not present or digest and present != digest:
                return False
        if any(not snapshot.get("references", {}).get(p) for p in needs.get("references", [])):
            return False
        return True

    @staticmethod
    def capacity(node, needs, reserved):
        snapshot = node["snapshot"]
        available = snapshot.get("resources", {})
        if reserved:
            return False  # Await a fresh local-scheduler acknowledgement before another dispatch.
        gpu_count = int(needs.get("gpu_count", 0))
        if needs.get("all_gpus"):
            gpu_count = max(gpu_count, len(available.get("gpus", [])))
        requested = set(str(v) for v in needs.get("gpu_ids", []))
        gpus = [g for g in available.get("free_gpus", []) if float(g.get("free_vram_gb", 0)) >= float(needs.get("min_vram_gb", 0))
                and (not requested or str(g["gpu_id"]) in requested)]
        return (len(gpus) >= gpu_count and float(available.get("cpu_free", 0)) >= float(needs.get("cpu_threads", 1))
                and float(available.get("ram_available_gb", 0)) >= float(needs.get("ram_gb", 0))
                and float(available.get("scratch_free_gb", 0)) >= float(needs.get("scratch_gb", 0)))

    def dispatch(self):
        nodes = self.nodes()
        with self.connect() as db:
            db.execute("BEGIN IMMEDIATE")
            if config.is_peer():
                online = {n["id"] for n in nodes if n["availability"] == "online"}
                queued = db.execute("SELECT s.job_id,s.command_id,c.node_id FROM submissions s JOIN commands c ON c.id=s.command_id WHERE s.placement='auto' AND s.affinity_node IS NULL AND s.cancelled=0 AND c.status='queued'").fetchall()
                for waiting in queued:
                    if waiting["node_id"] in online:
                        continue
                    # Claim and reassignment use the same transaction: once claimed,
                    # a job remains pinned even if its acceptance reply is lost.
                    db.execute("UPDATE commands SET status='cancelled' WHERE id=? AND status='queued'", (waiting["command_id"],))
                    db.execute("UPDATE submissions SET command_id=NULL WHERE job_id=?", (waiting["job_id"],))
                    value = json.loads(self._job(db, waiting["job_id"])["payload"])
                    value.update(node_id=None, status="waiting_for_workstation", wait_reason="The selected workstation disconnected before accepting; selecting an available workstation.")
                    db.execute("UPDATE jobs SET node_id=NULL,payload=? WHERE id=?", (encoded(value), waiting["job_id"]))
            rows = db.execute("SELECT * FROM submissions WHERE command_id IS NULL AND cancelled=0").fetchall()
            for row in rows:
                target = row["affinity_node"] or (row["placement"] if row["placement"] != "auto" else None)
                needs = json.loads(row["requirements"])
                choices = [n for n in nodes if (not target or n["id"] == target) and self.compatible(n, needs)]
                if not target:
                    choices = [n for n in choices if self.capacity(n, needs, db.execute(
                        "SELECT 1 FROM commands WHERE node_id=? AND kind IN ('workflow','pipeline') AND status IN ('queued','claimed') LIMIT 1", (n["id"],)).fetchone())]
                choices.sort(key=lambda n: (n["id"] != row["origin_node"], -float(n["snapshot"].get("resources", {}).get("cpu_free", 0)), n["id"]))
                if not choices:
                    continue
                node_id = choices[0]["id"]
                cmd = str(uuid.uuid4())
                request = json.loads(row["request"])
                db.execute("INSERT INTO commands(id,node_id,job_id,kind,payload,owner_id,project_id,created_at) VALUES (?,?,?,?,?,?,?,?)",
                           (cmd, node_id, row["job_id"], request["kind"], row["request"], row["owner_id"], request["project_id"], stamp(self.clock())))
                db.execute("UPDATE submissions SET command_id=? WHERE job_id=?", (cmd, row["job_id"]))
                job = json.loads(self._job(db, row["job_id"])["payload"])
                job.update(node_id=node_id, status="dispatching", wait_reason="Inputs are being staged on the selected workstation.")
                db.execute("UPDATE jobs SET node_id=?,payload=? WHERE id=?", (node_id, encoded(job), row["job_id"]))

    def action(self, jid, kind, payload, user):
        if kind not in {"inspect", "artifact", "cancel"}:
            raise ValueError("Unsupported workstation action.")
        job = self.job(jid, user)
        if kind == "cancel" and not job.get("local_reference"):
            with self.connect() as db:
                db.execute("BEGIN IMMEDIATE")
                submission = db.execute("SELECT command_id FROM submissions WHERE job_id=?", (job["cluster_job_id"],)).fetchone()
                command = db.execute("SELECT status FROM commands WHERE id=?", (submission["command_id"],)).fetchone() if submission and submission["command_id"] else None
                if not command or command["status"] == "queued":
                    if submission and submission["command_id"]:
                        db.execute("UPDATE commands SET status='cancelled' WHERE id=?", (submission["command_id"],))
                    db.execute("UPDATE submissions SET cancelled=1 WHERE job_id=?", (job["cluster_job_id"],))
                    job.update(status="cancelled", wait_reason="Cancelled before the workstation accepted this job.")
                else:
                    job.update(cancel_requested=True, wait_reason="Cancellation requested; waiting for the workstation to confirm acceptance.")
                db.execute("UPDATE jobs SET payload=? WHERE id=?", (encoded(job), job["cluster_job_id"]))
            return {"status": job["status"], "cancel_requested": job.get("cancel_requested", False)}
        if not job.get("local_reference"):
            raise ValueError("This job has not been accepted by its workstation yet.")
        command_id = str(uuid.uuid4())
        with self.connect() as db:
            db.execute("INSERT INTO commands(id,node_id,job_id,kind,payload,owner_id,project_id,created_at) VALUES (?,?,?,?,?,?,?,?)",
                       (command_id, job["node_id"], job["cluster_job_id"], kind, encoded(payload | {"local_reference": job["local_reference"]}), user["id"], job["project_id"], stamp(self.clock())))
        return {"command_id": command_id, "status": "queued", "node_id": job["node_id"]}

    def claim(self, node_id):
        self.dispatch()
        with self.connect() as db:
            db.execute("BEGIN IMMEDIATE")
            # Expired commands are retried only on the original machine. Its receipt
            # journal prevents repeating a scientific submission after a lost reply.
            row = db.execute("SELECT * FROM commands WHERE node_id=? AND (status='queued' OR (status='claimed' AND lease_until<?)) ORDER BY created_at,id LIMIT 1",
                             (node_id, self.clock())).fetchone()
            if not row:
                return None
            token = secrets.token_urlsafe(24)
            db.execute("UPDATE commands SET status='claimed',lease_token=?,lease_until=? WHERE id=?", (token, self.clock() + 180, row["id"]))
            return {"id": row["id"], "job_id": row["job_id"], "kind": row["kind"], "payload": json.loads(row["payload"]),
                    "owner_id": row["owner_id"], "project_id": row["project_id"], "lease_token": token}

    def finish(self, node_id, command_id, lease_token, response=None, error=None):
        with self.connect() as db:
            db.execute("BEGIN IMMEDIATE")
            row = db.execute("SELECT * FROM commands WHERE id=? AND node_id=?", (command_id, node_id)).fetchone()
            if not row or row["status"] != "claimed" or not hmac.compare_digest(str(row["lease_token"]), str(lease_token)):
                raise PermissionError("This command lease is no longer owned by this workstation.")
            db.execute("UPDATE commands SET status=?,response=?,error=?,lease_token=NULL WHERE id=?",
                       ("failed" if error else "completed", encoded(response or {}), str(error)[:8000] if error else None, command_id))
            if row["kind"] in {"workflow", "pipeline"}:
                job = json.loads(self._job(db, row["job_id"])["payload"])
                if error:
                    job.update(status="submission_uncertain" if "uncertain" in str(error).casefold() else "failed", error=str(error), wait_reason="")
                else:
                    accepted = dict(response["job"])
                    local_ref = accepted["task_group"] + "/" + accepted["run_id"]
                    job.update(accepted, cluster_job_id=row["job_id"], node_id=node_id,
                               local_reference=local_ref, location=accepted.get("run_dir"),
                               run_dir=config.artifact_uri(node_id, accepted["task_group"], accepted["run_id"]), wait_reason="")
                    # Heartbeats may have indexed the accepted run before its reply.
                    existing = db.execute("SELECT id FROM jobs WHERE node_id=? AND local_ref=?", (node_id, local_ref)).fetchone()
                    if existing and existing["id"] != row["job_id"]:
                        db.execute("DELETE FROM jobs WHERE id=?", (existing["id"],))
                    db.execute("UPDATE jobs SET local_ref=? WHERE id=?", (local_ref, row["job_id"]))
                    if job.get("cancel_requested"):
                        db.execute("INSERT INTO commands(id,node_id,job_id,kind,payload,owner_id,project_id,created_at) VALUES (?,?,?,?,?,?,?,?)",
                                   (str(uuid.uuid4()), node_id, row["job_id"], "cancel", encoded({"local_reference": local_ref}),
                                    row["owner_id"], row["project_id"], stamp(self.clock())))
                db.execute("UPDATE jobs SET payload=? WHERE id=?", (encoded(job), row["job_id"]))

    def command(self, command_id, user):
        with self.connect() as db:
            row = db.execute("SELECT * FROM commands WHERE id=?", (command_id,)).fetchone()
        if not row or not visible(dict(row), user):
            raise LookupError("Operation not found.")
        return {"command_id": row["id"], "status": row["status"], "response": json.loads(row["response"] or "{}"), "error": row["error"]}

    def cleanup_transfers(self, *, max_age_seconds=86400):
        """Expire requested download copies, retaining inputs awaiting acceptance."""
        with self.connect() as db:
            pending = db.execute("SELECT s.request FROM submissions s LEFT JOIN commands c ON c.id=s.command_id WHERE s.cancelled=0 AND (s.command_id IS NULL OR c.status IN ('queued','claimed'))").fetchall()
        keep = {json.loads(r["request"]).get("bundle_id") for r in pending}
        removed = 0
        for path in (self.root / "transfer-cache").glob("*"):
            if path.name not in keep and path.is_file() and self.clock() - path.stat().st_mtime > max_age_seconds:
                path.unlink()
                removed += 1
        return removed

    def reference_probes(self):
        with self.connect() as db:
            rows = db.execute("SELECT s.requirements FROM submissions s LEFT JOIN commands c ON c.id=s.command_id WHERE s.cancelled=0 AND (s.command_id IS NULL OR c.status IN ('queued','claimed'))").fetchall()
        return sorted({p for row in rows for p in json.loads(row["requirements"]).get("references", [])})[:10000]
