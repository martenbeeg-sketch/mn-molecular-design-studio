"""Coordinator metadata only. Run data is never resolved on this machine."""
from __future__ import annotations

import hashlib
import hmac
import json
import re
import secrets
import sqlite3
import time
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

from . import config

TERMINAL = {"completed", "failed", "cancelled", "succeeded"}
ACTIVE_RESOURCE_STATUSES = {
    "running", "queued", "preparing", "pending", "starting", "dispatching",
    "cancelling", "cancel_requested", "submitted",
}


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
                CREATE TABLE IF NOT EXISTS library_transfers (
                    job_id TEXT PRIMARY KEY, node_id TEXT NOT NULL, origin_node TEXT NOT NULL,
                    destination TEXT NOT NULL, library_id TEXT, files TEXT NOT NULL DEFAULT '[]',
                    status TEXT NOT NULL DEFAULT 'waiting', total_bytes INTEGER NOT NULL DEFAULT 0,
                    transferred_bytes INTEGER NOT NULL DEFAULT 0, error TEXT, updated_at TEXT NOT NULL);
                CREATE UNIQUE INDEX IF NOT EXISTS library_transfers_destination ON library_transfers(destination);
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
                external_library_record = row.pop("external_library_record", None)
                group, run = str(row.get("task_group") or ""), str(row.get("run_id") or "")
                if not group or not run or "/" in group + run or ".." in {group, run}:
                    raise ValueError("Invalid local job reference.")
                local_ref = group + "/" + run
                existing = db.execute("SELECT * FROM jobs WHERE node_id=? AND local_ref=?", (node_id, local_ref)).fetchone()
                jid = str(row.get("cluster_job_id") or config.job_id(node_id, local_ref))
                purged = db.execute(
                    "SELECT 1 FROM commands WHERE job_id=? AND kind='purge' AND status='completed' "
                    "AND json_extract(response,'$.status')='purged' LIMIT 1",
                    (jid,),
                ).fetchone()
                if purged:
                    # A successful purge is authoritative even if an older
                    # workstation catalog still reports the removed run.
                    continue
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
                if collision:
                    # A pipeline child can be indexed before its project marker
                    # has propagated into the shared catalog. Allow the owning
                    # workstation to backfill those legacy empty fields for the
                    # exact same local run, while preserving any established
                    # ownership and rejecting conflicting identities.
                    if collision["node_id"] not in {None, node_id} or collision["local_ref"] not in {None, local_ref}:
                        raise PermissionError("A workstation cannot change job ownership or adopt another workstation's run.")
                    for field in ("owner_id", "project_id"):
                        previous_value = str(collision[field] or "")
                        incoming_value = str(row.get(field) or "")
                        if previous_value and incoming_value and previous_value != incoming_value:
                            raise PermissionError("A workstation cannot change job ownership or project.")
                        if previous_value and not incoming_value:
                            row[field] = previous_value
                if collision:
                    previous = json.loads(collision["payload"])
                    row = previous | superseded_payload | row
                    if previous.get("hidden") or superseded_payload.get("hidden"):
                        row["hidden"] = True
                    if previous.get("cancel_requested") or superseded_payload.get("cancel_requested"):
                        row["cancel_requested"] = True
                row.update(cluster_job_id=jid, node_id=node_id, local_reference=local_ref,
                           run_dir=config.artifact_uri(node_id, group, run), location=row.get("location") or raw.get("run_dir"),
                           data_state="present", wait_reason="")
                db.execute("INSERT INTO jobs VALUES (?,?,?,?,?,?,?) ON CONFLICT(id) DO UPDATE SET local_ref=excluded.local_ref,owner_id=excluded.owner_id,project_id=excluded.project_id,payload=excluded.payload",
                           (jid, node_id, local_ref, str(row.get("owner_id") or ""), str(row.get("project_id") or ""), encoded(row), row.get("created_at") or stamp(self.clock())))
                if external_library_record is not None:
                    if row.get("tool") not in {"library_preparation", "library_centroids"}:
                        raise PermissionError("Only a library preparation job can publish external library metadata.")
                    from mn_molecular_design_studio.external_libraries import LibraryCatalog
                    LibraryCatalog().import_remote_record(external_library_record, node_id=node_id)
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
            transfer = db.execute(
                "SELECT status,total_bytes,transferred_bytes,error FROM library_transfers WHERE job_id=?",
                (row["id"],),
            ).fetchone()
        if not visible(row, user):
            raise LookupError("Job not found.")
        item = json.loads(row["payload"])
        node = next((n for n in self.nodes() if n["id"] == row["node_id"]), None)
        item.update(node_name=node["name"] if node else "Not assigned", availability=node["availability"] if node else "pending",
                    last_seen_at=node["last_seen_at"] if node else None)
        if transfer:
            item["library_transfer"] = {
                "status": transfer["status"], "total_bytes": transfer["total_bytes"],
                "transferred_bytes": transfer["transferred_bytes"], "error": transfer["error"],
            }
        if node and node["availability"] == "online" and not node["snapshot"].get("storage_available", True):
            item["availability"] = "storage_unavailable"
        elif node and node["availability"] == "online" and item.get("data_state") == "not_found":
            item["availability"] = "files_unavailable"
        return item

    def search(self, user, *, node_id=None, project_id=None, query="", statuses=None, tools=None, limit=100, offset=0, include_hidden=False):
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
            if not include_hidden and item.get("hidden"):
                continue
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
                destination = request.get("organizer_output_directory")
                if destination:
                    if request.get("kind") != "workflow" or request.get("payload", {}).get("workflow") not in {"library.prepare", "library.centroids"}:
                        raise ValueError("Only remote library preparation jobs can register an output transfer.")
                    if db.execute("SELECT 1 FROM library_transfers WHERE destination=?", (str(destination),)).fetchone():
                        raise ValueError("This organizer destination is already assigned to another library job.")
                    db.execute(
                        "INSERT INTO library_transfers(job_id,node_id,origin_node,destination,status,updated_at) VALUES (?,?,?,?,?,?)",
                        (jid, placement, origin_node, str(destination), "waiting", stamp(self.clock())),
                    )
        self.dispatch()
        return {"job": self.job(jid, user)}

    @staticmethod
    def compatible(node, needs):
        snapshot = node["snapshot"]
        if node["availability"] != "online" or not snapshot.get("accept_jobs") or not snapshot.get("scheduler_online"):
            return False
        if (needs.get("code_revision") and snapshot.get("code_revision") != needs["code_revision"]
                and not needs.get("_allow_code_revision_refresh")):
            return False
        for tag, digest in needs.get("images", {}).items():
            present = snapshot.get("images", {}).get(tag)
            if not present or digest and present != digest:
                return False
        if any(not snapshot.get("references", {}).get(p) for p in needs.get("references", [])):
            return False
        if any(snapshot.get('reference_hashes', {}).get(path) != checksum
               for path, checksum in needs.get('reference_hashes', {}).items()):
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
        from .gpu_runtime import compatible as runtime_compatible
        gpus = [g for g in available.get("free_gpus", []) if float(g.get("free_vram_gb", 0)) >= float(needs.get("min_vram_gb", 0))
                and (not requested or str(g["gpu_id"]) in requested) and runtime_compatible(g, needs.get("gpu_runtime_profiles", {}))]
        return (len(gpus) >= gpu_count and float(available.get("cpu_free", 0)) >= float(needs.get("cpu_threads", 1))
                and float(available.get("ram_available_gb", 0)) >= float(needs.get("ram_gb", 0))
                and float(available.get("scratch_free_gb", 0)) >= float(needs.get("scratch_gb", 0)))

    @staticmethod
    def _active_job_resource_usage(db, node_id):
        """Read running allocations from the synchronized results catalog.

        The scheduler resource snapshot is useful for free-memory and CPU
        telemetry, but older or externally-started runs can be absent from its
        reservation list. The synchronized job catalog is the second source of
        truth for resources already in use on a workstation.
        """
        busy_gpu_ids = set()
        all_gpus_busy = False
        cpu_in_use = 0.0
        rows = db.execute("SELECT payload FROM jobs WHERE node_id=?", (node_id,)).fetchall()
        for row in rows:
            try:
                job = json.loads(row["payload"])
            except (TypeError, json.JSONDecodeError):
                continue
            if str(job.get("status") or "").casefold() not in ACTIVE_RESOURCE_STATUSES:
                continue

            request = job.get("scheduler_resource_request") or job.get("resource_request")
            if isinstance(request, dict):
                cpu_in_use += max(0.0, float(request.get("cpu_threads") or request.get("cpu_cores") or 0))
                gpu_ids = request.get("gpu_ids")
                if isinstance(gpu_ids, (list, tuple, set)):
                    busy_gpu_ids.update(str(value) for value in gpu_ids)
                elif request.get("gpu"):
                    selected = str(request.get("gpu_device") or "").strip().lower()
                    if selected in {"", "all", "auto", "any"}:
                        all_gpus_busy = True
                    else:
                        busy_gpu_ids.update(re.findall(r"\d+", selected))

            queue_resource = str(job.get("queue_resource") or "").strip().lower()
            if queue_resource.startswith("gpu"):
                selected = queue_resource.partition(":")[2].strip()
                if not selected or selected in {"all", "any", "auto"}:
                    all_gpus_busy = True
                else:
                    busy_gpu_ids.update(re.findall(r"\d+", selected))

            allocation = str(job.get("resource_allocation") or "")
            match = re.search(r"(\d+(?:\.\d+)?)\s+CPU cores?", allocation, flags=re.IGNORECASE)
            if match and not (isinstance(request, dict) and (request.get("cpu_threads") or request.get("cpu_cores"))):
                cpu_in_use += float(match.group(1))

        return {"busy_gpu_ids": busy_gpu_ids, "all_gpus_busy": all_gpus_busy, "cpu_in_use": cpu_in_use}

    @staticmethod
    def _with_active_job_resources(node, usage):
        effective = dict(node)
        snapshot = dict(node.get("snapshot") or {})
        resources = dict(snapshot.get("resources") or {})
        if usage["all_gpus_busy"]:
            resources["free_gpus"] = []
        elif usage["busy_gpu_ids"]:
            resources["free_gpus"] = [
                gpu for gpu in resources.get("free_gpus", [])
                if str(gpu.get("gpu_id")) not in usage["busy_gpu_ids"]
            ]
        if usage["cpu_in_use"] > 0:
            reported_free = float(resources.get("cpu_free", 0))
            if resources.get("cpu_total") is not None:
                active_limit = max(0.0, float(resources["cpu_total"]) - usage["cpu_in_use"])
                resources["cpu_free"] = min(reported_free, active_limit)
            else:
                resources["cpu_free"] = max(0.0, reported_free - usage["cpu_in_use"])
        snapshot["resources"] = resources
        effective["snapshot"] = snapshot
        return effective

    def placement_options(self, needs, *, origin_node, affinity_node=None):
        """List currently compatible and available workstations without queuing."""
        nodes = self.nodes()
        options = []
        with self.connect() as db:
            for node in nodes:
                if affinity_node and node["id"] != affinity_node:
                    continue
                if not self.compatible(node, needs):
                    continue
                reserved = db.execute(
                    "SELECT 1 FROM commands WHERE node_id=? AND kind IN ('workflow','pipeline') "
                    "AND status IN ('queued','claimed') LIMIT 1", (node["id"],),
                ).fetchone()
                usage = self._active_job_resource_usage(db, node["id"])
                effective = self._with_active_job_resources(node, usage)
                if not self.capacity(effective, needs, reserved):
                    continue
                resources = effective["snapshot"].get("resources", {})
                options.append({
                    "node_id": node["id"], "name": node["name"],
                    "cpu_free": resources.get("cpu_free", 0),
                    "free_gpus": resources.get("free_gpus", []),
                })
        options.sort(key=lambda item: (
            item["node_id"] != origin_node, -float(item.get("cpu_free", 0)), item["node_id"],
        ))
        return options

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
                choices = []
                wait_reason = "Waiting for a compatible workstation."
                for node in nodes:
                    if target and node["id"] != target:
                        continue
                    if node["availability"] != "online":
                        wait_reason = f"Waiting for selected workstation {node['name']} to come online."
                        continue
                    # A job pinned to its submitting workstation should not
                    # become permanently stuck after an unrelated Studio code
                    # update. Keep the exact image and checkpoint checks, but
                    # allow that same workstation to validate the current
                    # request against its own code when accepting the job.
                    refresh_local_code = (
                        row["placement"] != "auto"
                        and row["placement"] == row["origin_node"]
                        and node["id"] == row["origin_node"]
                    )
                    compatible_needs = dict(needs)
                    if refresh_local_code:
                        compatible_needs["_allow_code_revision_refresh"] = True
                    if not self.compatible(node, compatible_needs):
                        snapshot = node["snapshot"]
                        missing_images = [
                            tag for tag, digest in needs.get("images", {}).items()
                            if not snapshot.get("images", {}).get(tag)
                            or digest and snapshot["images"].get(tag) != digest
                        ]
                        missing_references = [
                            path for path in needs.get("references", [])
                            if not snapshot.get("references", {}).get(path)
                        ]
                        if missing_images:
                            wait_reason = f"{node['name']} is missing the required image: {', '.join(missing_images)}."
                        elif missing_references:
                            wait_reason = f"{node['name']} is missing required reference files: {', '.join(missing_references)}."
                        elif any(snapshot.get('reference_hashes', {}).get(path) != checksum for path, checksum in needs.get('reference_hashes', {}).items()):
                            wait_reason = f"{node['name']} has different checkpoint contents or has not yet reported their SHA-256 identities."
                        elif needs.get("code_revision") and snapshot.get("code_revision") != needs["code_revision"]:
                            wait_reason = f"{node['name']} is running a different Studio code revision."
                        else:
                            wait_reason = f"{node['name']} is not accepting jobs or its scheduler is offline."
                        continue
                    reserved = db.execute(
                        "SELECT 1 FROM commands WHERE node_id=? AND kind IN ('workflow','pipeline') "
                        "AND status IN ('queued','claimed') LIMIT 1", (node["id"],),
                    ).fetchone()
                    usage = self._active_job_resource_usage(db, node["id"])
                    effective = self._with_active_job_resources(node, usage)
                    if self.capacity(effective, needs, reserved):
                        choices.append((node, effective, refresh_local_code))
                    elif target:
                        resources = effective["snapshot"].get("resources", {})
                        requested_gpus = {str(value) for value in needs.get("gpu_ids", [])}
                        free_gpus = resources.get("free_gpus", [])
                        free_ids = {str(gpu.get("gpu_id")) for gpu in free_gpus
                                    if float(gpu.get("free_vram_gb", 0)) >= float(needs.get("min_vram_gb", 0))}
                        busy_requested = sorted(requested_gpus.intersection(usage["busy_gpu_ids"]))
                        if busy_requested:
                            label = ", ".join(f"GPU {gpu}" for gpu in busy_requested)
                            wait_reason = f"Waiting for {label} to become available on {node['name']}."
                        elif reserved:
                            wait_reason = f"Waiting for {node['name']} to accept its queued job before staging another."
                        elif requested_gpus and not requested_gpus.intersection(free_ids):
                            label = ", ".join(f"GPU {gpu}" for gpu in sorted(requested_gpus))
                            wait_reason = f"Waiting for {label} with enough free memory on {node['name']}."
                        elif float(resources.get("cpu_free", 0)) < float(needs.get("cpu_threads", 1)):
                            wait_reason = f"Waiting for {needs.get('cpu_threads', 1)} free CPU threads on {node['name']}."
                        else:
                            wait_reason = f"Waiting for the requested GPU and memory to become available on {node['name']}."
                choices.sort(key=lambda pair: (
                    pair[0]["id"] != row["origin_node"],
                    -float(pair[1]["snapshot"].get("resources", {}).get("cpu_free", 0)),
                    pair[0]["id"],
                ))
                if not choices:
                    job = json.loads(self._job(db, row["job_id"])["payload"])
                    job.update(status="waiting_for_workstation", wait_reason=wait_reason)
                    db.execute("UPDATE jobs SET payload=? WHERE id=?", (encoded(job), row["job_id"]))
                    continue
                node_id, refresh_local_code = choices[0][0]["id"], choices[0][2]
                cmd = str(uuid.uuid4())
                request = json.loads(row["request"])
                request["origin_node"] = row["origin_node"]
                if refresh_local_code:
                    request["allow_code_revision_refresh"] = True
                db.execute("INSERT INTO commands(id,node_id,job_id,kind,payload,owner_id,project_id,created_at) VALUES (?,?,?,?,?,?,?,?)",
                           (cmd, node_id, row["job_id"], request["kind"], encoded(request), row["owner_id"], request["project_id"], stamp(self.clock())))
                db.execute("UPDATE submissions SET command_id=? WHERE job_id=?", (cmd, row["job_id"]))
                job = json.loads(self._job(db, row["job_id"])["payload"])
                job.update(node_id=node_id, status="dispatching", wait_reason="Inputs are being staged on the selected workstation.")
                db.execute("UPDATE jobs SET node_id=?,payload=? WHERE id=?", (node_id, encoded(job), row["job_id"]))

    def action(self, jid, kind, payload, user):
        if kind not in {"inspect", "artifact", "cancel", "hide", "restore", "retry", "change_gpu", "purge", "md_source_options", "collect_accepted_designs"}:
            raise ValueError("Unsupported workstation action.")
        job = self.job(jid, user)
        if kind == "collect_accepted_designs":
            if job.get("status") not in {"completed", "succeeded"}:
                raise ValueError("Only completed validation results can be promoted to the shared Database.")
            counts = (job.get("candidate_summary") or {}).get("acceptance_counts") or {}
            if int(counts.get("accepted") or 0) < 1:
                raise ValueError("This validation job has no explicitly accepted designs to promote.")
        if kind == "purge" and job.get("status") not in TERMINAL:
            raise ValueError("Only finished jobs can be permanently removed.")
        if kind == "change_gpu":
            if job.get("status") != "waiting_for_workstation" or job.get("local_reference"):
                raise ValueError("A GPU can only be changed while a job is waiting and has not been accepted.")
            try:
                gpu_id = int(payload["gpu_id"])
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError("Provide the target gpu_id as a non-negative integer.") from exc
            if gpu_id < 0:
                raise ValueError("gpu_id must be non-negative.")
            with self.connect() as db:
                db.execute("BEGIN IMMEDIATE")
                submission = db.execute(
                    "SELECT request,requirements,command_id,placement,origin_node FROM submissions WHERE job_id=?",
                    (job["cluster_job_id"],),
                ).fetchone()
                if (not submission or submission["command_id"]
                        or submission["placement"] == "auto"
                        or submission["placement"] != submission["origin_node"]):
                    raise ValueError("GPU changes are supported only for an unclaimed job pinned to its submitting workstation.")
                node = next((item for item in self.nodes() if item["id"] == submission["placement"]), None)
                if not node or node["availability"] != "online":
                    raise ValueError("The selected workstation is not online.")
                requirements = json.loads(submission["requirements"])
                minimum_vram = float(requirements.get("min_vram_gb", 0))
                available = next((gpu for gpu in node["snapshot"].get("resources", {}).get("free_gpus", [])
                                  if int(gpu.get("gpu_id", -1)) == gpu_id
                                  and float(gpu.get("free_vram_gb", 0)) >= minimum_vram), None)
                if not available:
                    raise ValueError(f"GPU {gpu_id} is not currently free with the required memory on {node['name']}.")
                request = json.loads(submission["request"])
                workflow_payload = request.get("payload") or {}
                parameters = workflow_payload.setdefault("parameters", {})
                parameters["gpu_device"] = str(gpu_id)
                resources = workflow_payload.get("resources")
                if isinstance(resources, dict) and "gpu_device" in resources:
                    resources["gpu_device"] = str(gpu_id)
                requirements["gpu_count"] = max(1, int(requirements.get("gpu_count", 1)))
                requirements["gpu_ids"] = [gpu_id]
                request["requirements"] = requirements
                db.execute("UPDATE submissions SET request=?,requirements=? WHERE job_id=?",
                           (encoded(request), encoded(requirements), job["cluster_job_id"]))
                job["request_summary"] = workflow_payload
                job.update(wait_reason=f"Rechecking placement on GPU {gpu_id}.")
                db.execute("UPDATE jobs SET payload=? WHERE id=?", (encoded(job), job["cluster_job_id"]))
            self.dispatch()
            return self.job(job["cluster_job_id"], user)
        if kind == "retry":
            if job.get("status") != "failed" or job.get("local_reference"):
                raise ValueError("Only failed jobs that were not accepted by a workstation can be retried.")
            if job.get("error") != "The workstation's engine images, code, or checkpoints no longer match this submission.":
                raise ValueError("This job did not fail during workstation compatibility checks and cannot be retried here.")
            with self.connect() as db:
                db.execute("BEGIN IMMEDIATE")
                submission = db.execute("SELECT command_id FROM submissions WHERE job_id=?", (job["cluster_job_id"],)).fetchone()
                if not submission or not submission["command_id"]:
                    raise ValueError("No failed workstation acceptance attempt is available to retry.")
                previous = db.execute("SELECT status FROM commands WHERE id=?", (submission["command_id"],)).fetchone()
                if not previous or previous["status"] != "failed":
                    raise ValueError("The workstation acceptance attempt is not in a failed state.")
                db.execute("UPDATE submissions SET command_id=NULL,cancelled=0 WHERE job_id=?", (job["cluster_job_id"],))
                job.update(status="waiting_for_workstation", error="", wait_reason="Rechecking the selected workstation.")
                db.execute("UPDATE jobs SET payload=? WHERE id=?", (encoded(job), job["cluster_job_id"]))
            self.dispatch()
            return self.job(job["cluster_job_id"], user)
        if kind in {"hide", "restore"}:
            job["hidden"] = kind == "hide"
            with self.connect() as db:
                db.execute("UPDATE jobs SET payload=? WHERE id=?", (encoded(job), job["cluster_job_id"]))
            return {"hidden": job["hidden"]}
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
        if kind == "purge":
            with self.connect() as db:
                existing = db.execute(
                    "SELECT id,status FROM commands WHERE job_id=? AND kind='purge' "
                    "AND status IN ('queued','claimed') ORDER BY created_at DESC LIMIT 1",
                    (job["cluster_job_id"],),
                ).fetchone()
            if existing:
                return {"command_id": existing["id"], "status": existing["status"], "node_id": job["node_id"]}
        command_id = str(uuid.uuid4())
        with self.connect() as db:
            db.execute("INSERT INTO commands(id,node_id,job_id,kind,payload,owner_id,project_id,created_at) VALUES (?,?,?,?,?,?,?,?)",
                       (command_id, job["node_id"], job["cluster_job_id"], kind, encoded(payload | {"local_reference": job["local_reference"]}), user["id"], job["project_id"], stamp(self.clock())))
        return {"command_id": command_id, "status": "queued", "node_id": job["node_id"]}

    def purge_context(self, jid, user):
        """Return the terminal job and staging command needed for local cleanup."""
        job = self.job(jid, user)
        if job.get("status") not in TERMINAL:
            raise ValueError("Only finished jobs can be permanently removed.")
        with self.connect() as db:
            submission = db.execute(
                "SELECT command_id,request FROM submissions WHERE job_id=?",
                (job["cluster_job_id"],),
            ).fetchone()
        request = json.loads(submission["request"]) if submission else {}
        return {
            "job": job,
            "command_id": str(submission["command_id"] or "") if submission else "",
            "bundle_id": str(request.get("bundle_id") or ""),
        }

    def purge_metadata(self, jid, user):
        """Remove coordinator records for an already-cleaned terminal job."""
        job_id = str(jid).removeprefix("cluster:")
        with self.connect() as db:
            db.execute("BEGIN IMMEDIATE")
            row = self._job(db, job_id)
            if not visible(row, user):
                raise LookupError("Job not found.")
            payload = json.loads(row["payload"])
            if payload.get("status") not in TERMINAL:
                raise ValueError("Only finished jobs can be permanently removed.")
            submission = db.execute(
                "SELECT request FROM submissions WHERE job_id=?", (job_id,)
            ).fetchone()
            request = json.loads(submission["request"]) if submission else {}
            bundle_id = str(request.get("bundle_id") or "")
            command_count = db.execute(
                "SELECT count(*) FROM commands WHERE job_id=?", (job_id,)
            ).fetchone()[0]
            db.execute("DELETE FROM commands WHERE job_id=?", (job_id,))
            db.execute("DELETE FROM submissions WHERE job_id=?", (job_id,))
            db.execute("DELETE FROM library_transfers WHERE job_id=?", (job_id,))
            db.execute("DELETE FROM jobs WHERE id=?", (job_id,))

            # Input bundles can be shared across a compound batch. Return a
            # digest for disk cleanup only when no remaining job or command
            # references it.
            still_referenced = False
            if bundle_id:
                for table, columns in (("submissions", ("request",)), ("commands", ("payload", "response"))):
                    for row in db.execute("SELECT " + ",".join(columns) + " FROM " + table):
                        for column in columns:
                            try:
                                value = json.loads(row[column] or "{}")
                            except (TypeError, ValueError):
                                continue
                            if bundle_id in json.dumps(value, separators=(",", ":")):
                                still_referenced = True
                                break
                        if still_referenced:
                            break
                    if still_referenced:
                        break
            return {
                "cluster_job_id": job_id,
                "removed_commands": int(command_count),
                "removed_submission": submission is not None,
                "unreferenced_bundle_id": bundle_id if bundle_id and not still_referenced else "",
            }

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
            return {"id": row["id"], "job_id": row["job_id"], "node_id": row["node_id"],
                    "kind": row["kind"], "payload": json.loads(row["payload"]),
                    "owner_id": row["owner_id"], "project_id": row["project_id"], "lease_token": token}

    def finish(self, node_id, command_id, lease_token, response=None, error=None):
        unreferenced_bundle_id = ""
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
            elif row["kind"] == "purge" and not error:
                # The owning workstation has already removed the run directory
                # and run-derived Database rows. Keep this completed command so
                # the Results page can show the final remote-operation status.
                submission = db.execute(
                    "SELECT request FROM submissions WHERE job_id=?", (row["job_id"],)
                ).fetchone()
                request = json.loads(submission["request"]) if submission else {}
                bundle_id = str(request.get("bundle_id") or "")
                db.execute("DELETE FROM commands WHERE job_id=? AND id<>?", (row["job_id"], command_id))
                db.execute("DELETE FROM submissions WHERE job_id=?", (row["job_id"],))
                db.execute("DELETE FROM library_transfers WHERE job_id=?", (row["job_id"],))
                db.execute("DELETE FROM jobs WHERE id=?", (row["job_id"],))
                still_referenced = False
                if bundle_id:
                    for table, columns in (("submissions", ("request",)), ("commands", ("payload", "response"))):
                        for item in db.execute("SELECT " + ",".join(columns) + " FROM " + table):
                            for column in columns:
                                try:
                                    value = json.loads(item[column] or "{}")
                                except (TypeError, ValueError):
                                    continue
                                if bundle_id in json.dumps(value, separators=(",", ":")):
                                    still_referenced = True
                                    break
                            if still_referenced:
                                break
                        if still_referenced:
                            break
                if bundle_id and not still_referenced:
                    unreferenced_bundle_id = bundle_id

        if re.fullmatch(r"[0-9a-f]{64}", unreferenced_bundle_id):
            blob = self.root / "transfer-cache" / unreferenced_bundle_id
            if blob.is_file() and not blob.is_symlink():
                blob.unlink()

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
            rows = db.execute(
                "SELECT s.requirements,c.status AS command_status,j.payload AS job_payload "
                "FROM submissions s "
                "LEFT JOIN commands c ON c.id=s.command_id "
                "LEFT JOIN jobs j ON j.id=s.job_id "
                "WHERE s.cancelled=0 ORDER BY s.rowid DESC"
            )
            probes = set()
            for row in rows:
                job_status = ""
                try:
                    job_status = str(json.loads(row["job_payload"] or "{}").get("status") or "").casefold()
                except (TypeError, ValueError, AttributeError):
                    pass
                command_pending = row["command_status"] in {None, "queued", "claimed"}
                if not command_pending and job_status not in ACTIVE_RESOURCE_STATUSES:
                    continue
                probes.update(json.loads(row["requirements"]).get("references", []))
                if len(probes) >= 10000:
                    break
        return sorted(probes)[:10000]
