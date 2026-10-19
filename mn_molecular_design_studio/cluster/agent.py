"""Outbound workstation agent; local schedulers continue owning computation."""
from __future__ import annotations

import fcntl
import hashlib
import json
import logging
import signal
import threading
import time
from pathlib import Path
from urllib.parse import urlsplit

from . import bundles, capabilities, config, context, progress as job_progress, transport

log = logging.getLogger(__name__)


class Agent:
    def __init__(self):
        self.settings = config.load()
        self.node_id = config.identity()
        self._inventory_time = 0
        self._inventory = {}
        self._published = {}
        self._artifact_signatures = {}
        self._artifact_manifests = {}
        self._known_references = set()
        self._reference_probes = []
        self._stop = threading.Event()
        self._peer_retry = {}
        self._poll_index = 0

    def snapshot(self):
        if time.monotonic() - self._inventory_time > 120 or not self._inventory:
            self._inventory = {"code_revision": capabilities.code_revision(), "images": capabilities.images()}
            self._inventory_time = time.monotonic()
        online, resources = capabilities.resources()
        storage_online = self.storage_online()
        return self._inventory | {"scheduler_online": online, "resources": resources,
                                  "references": capabilities.reference_inventory(self._reference_probes),
                                  "accept_jobs": config.load()["accept_jobs"] and storage_online,
                                  "storage_available": storage_online}

    def storage_online(self):
        from mn_molecular_design_studio.runtime import runs_root
        root = runs_root()
        state = config.home() / "storage-location.json"
        marker = root / ".studio-workstation"
        previous = json.loads(state.read_text()) if state.exists() else {}
        if previous.get("root") == str(root):
            return marker.is_file() and marker.read_text().strip() == self.node_id
        if not root.is_dir():
            return False
        if marker.exists() and marker.read_text().strip() != self.node_id:
            raise ValueError("This Results directory is owned by another workstation. Select a local Results directory.")
        marker.write_text(self.node_id + "\n")
        config.atomic_json(state, {"root": str(root), "node_id": self.node_id})
        return True

    def catalog(self):
        from mn_molecular_design_studio import job_catalog
        from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
        result, offset = [], 0
        self._known_references = set()
        while True:
            page = job_catalog.search_jobs(access_admin=True, limit=1000, offset=offset)
            for row in page["items"]:
                run = Path(row["run_dir"])
                artifact_signature = (
                    row.get("metadata_mtime_ns"),
                    row.get("input_mtime_ns"),
                    row.get("result_mtime_ns"),
                )
                progress_mtime = job_progress.bindcraft2_progress_mtime_ns(run)
                signature = (*artifact_signature, progress_mtime)
                reference = row["task_group"] + "/" + row["run_id"]
                self._known_references.add(reference)
                if self._published.get(reference) == signature:
                    continue
                metadata = read_json(run / "metadata.json")
                marker = read_json(run / "molecular_design_studio.json")
                metrics = read_json(run / "result.json").get("metrics", {})
                jid = (
                    metadata.get("cluster_job_id")
                    or marker.get("cluster_job_id")
                    or metadata.get("cluster_submission_id")
                    or marker.get("cluster_submission_id")
                    or config.job_id(self.node_id, reference)
                )
                if self._artifact_signatures.get(reference) != artifact_signature:
                    manifest = job_catalog.job_artifacts(row["task_group"], row["run_id"])
                    self._artifact_manifests[reference] = manifest[:1000]
                    self._artifact_signatures[reference] = artifact_signature
                progress = job_progress.bindcraft2_progress(run)
                item = dict(row, cluster_job_id=jid, metrics=metrics, progress=progress,
                            artifacts=self._artifact_manifests.get(reference, []),
                            location=str(run), last_reported_at=time.time(), _signature=list(signature))
                result.append(item)
            offset += len(page["items"])
            if offset >= page["total"] or not page["items"]:
                return result

    def heartbeat(self):
        if config.is_peer():
            return self.peer_heartbeat()
        storage_online = self.storage_online()
        jobs = self.catalog() if storage_online else []
        batches = [jobs[i:i + 1000] for i in range(0, len(jobs), 1000)] or [[]]
        snapshot = self.snapshot()
        for index, batch in enumerate(batches):
            body = {"snapshot": snapshot, "jobs": batch}
            if storage_online and index == len(batches) - 1:
                body["known_references"] = sorted(self._known_references)
            response = transport.node_request("/heartbeat", body, method="POST", timeout=45)
            for row in batch:
                self._published[row["task_group"] + "/" + row["run_id"]] = tuple(row["_signature"])
            if response.get("node_id") != self.node_id:
                raise RuntimeError("This credential belongs to another workstation. Check the enrollment file.")
            if self.settings["mode"] == "worker":
                config.atomic_json(config.home() / "organization-cache.json", response["organization"])
            self._reference_probes = response.get("reference_probes", [])
        config.atomic_json(config.home() / "agent-status.json", {"last_heartbeat": time.time(), "node_id": self.node_id, "online": True})
        if self.settings["mode"] == "coordinator":
            from .store import Store
            Store().cleanup_transfers()

    def peer_heartbeat(self):
        from .store import Store
        from .replication import Replica, capture_identity, sync_peer, locked
        store = Store()
        jobs = self.catalog() if self.storage_online() else []
        probes = set(store.reference_probes())
        for path in config.home().glob("peer-probes-*.json"):
            try:
                probes.update(json.loads(path.read_text()).get("references", []))
            except (OSError, ValueError):
                continue
        self._reference_probes = sorted(probes)
        snapshot = self.snapshot()
        for index in range(0, max(1, len(jobs)), 1000):
            batch = jobs[index:index+1000]
            known = sorted(self._known_references) if index + 1000 >= len(jobs) and self.storage_online() else None
            store.heartbeat(self.node_id, snapshot, batch, known)
            for row in batch:
                self._published[row["task_group"] + "/" + row["run_id"]] = tuple(row["_signature"])
        with locked():
            capture_identity()
            Replica().capture_jobs()
        for peer in config.load().get("peers", []):
            try:
                sync_peer(peer, snapshot)
            except Exception as exc:
                log.warning("Workstation %s synchronization paused: %s", peer.get("name",peer["id"]), exc)
        config.atomic_json(config.home() / "agent-status.json", {"last_heartbeat":time.time(), "node_id":self.node_id, "online":True})
        store.cleanup_transfers()

    def execute(self, command):
        directory = config.home() / "receipts"
        directory.mkdir(parents=True, exist_ok=True, mode=0o700)
        receipt = directory / (command["id"] + ".json")
        if receipt.exists():
            previous = json.loads(receipt.read_text())
            if previous["state"] == "finished":
                return previous["response"], previous.get("error")
            # A crash between local scheduler submission and recording its reply
            # cannot safely be treated as a new submission.
            return {}, "Submission state is uncertain after an agent restart. Reconcile the existing local run before retrying."
        config.atomic_json(receipt, {"state": "accepting", "job_id": command["job_id"], "command_id": command["id"]})
        response, error = {}, None
        try:
            if command["kind"] in {"workflow", "pipeline"}:
                response = self.submit(command)
            elif command["kind"] == "inspect":
                from mn_molecular_design_studio.engine.workflow_engine.services.local_automation import inspect_job, job_results
                ref = command["payload"]["local_reference"]
                response = inspect_job(ref) | {"results": job_results(ref, include_candidates=True)}
            elif command["kind"] == "cancel":
                from mn_molecular_design_studio.engine.workflow_engine.services.local_automation import cancel_job
                response = {"job": cancel_job(command["payload"]["local_reference"], reason="Cancelled through the shared workstation catalog")}
            elif command["kind"] == "artifact":
                response = self.artifact(command)
            else:
                raise ValueError("Unsupported workstation command.")
        except Exception as exc:
            error = str(exc)
            log.exception("Workstation operation failed: %s", command["id"])
        config.atomic_json(receipt, {"state": "finished", "response": response, "error": error})
        return response, error

    def submit(self, command):
        from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import submit_workflow_request, validate_workflow_request
        from mn_molecular_design_studio.workspace import mark_studio_job
        from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, write_json
        payload = command["payload"]
        needs = payload.get("requirements", {})
        snapshot = self.snapshot()
        from .store import Store
        if not Store.compatible({"availability": "online", "snapshot": snapshot}, needs):
            raise ValueError("The workstation's engine images, code, or checkpoints no longer match this submission.")
        staging = config.home() / "accepted-inputs" / command["id"]
        staging.mkdir(parents=True, exist_ok=True, mode=0o700)
        if payload.get("bundle_id"):
            digest = payload["bundle_id"]
            blob = bundles.blob_path(digest)
            if not blob.is_file():
                data = transport.node_request("/blobs/" + digest, raw=True, timeout=90)
                if hashlib.sha256(data).hexdigest() != digest:
                    raise ValueError("Input bundle failed its SHA-256 integrity check.")
                blob.write_bytes(data)
            bundles.unpack(blob, staging, expected_digest=digest)
        prepared = bundles.restore(payload["payload"], staging=staging, node_id=self.node_id)
        if payload["payload"].get("source_job", "").startswith("cluster://"):
            prepared["source_job"] = urlsplit(payload["payload"]["source_job"]).path.lstrip("/")
        with context.submission(command):
            if command["kind"] == "pipeline":
                from mn_molecular_design_studio.pipeline_graph import clean_graph, validate_graph
                from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import enqueue_graph_pipeline
                from mn_molecular_design_studio.engine.workflow_engine.services.local_automation import job_reference
                graph = clean_graph(prepared["graph"])
                validate_graph(graph, require_inputs=True)
                run = enqueue_graph_pipeline(graph, pipeline_name=graph["name"], owner_id=command["owner_id"],
                                             project_id=command["project_id"], reuse_completed=prepared.get("reuse_completed", True), gpu_device=prepared.get("gpu_device", "all"))
                result = {"job": job_reference(run)}
            else:
                request_file = staging / "request.json"
                config.atomic_json(request_file, prepared)
                validate_workflow_request(request_file)
                result = submit_workflow_request(request_file)
            result = mark_studio_job(result, project_id=command["project_id"], owner_id=command["owner_id"])
        run = Path(result["job"]["run_dir"])
        metadata = read_json(run / "metadata.json")
        metadata.update(cluster_job_id=command["job_id"], execution_node_id=self.node_id)
        write_json(run / "metadata.json", metadata)
        result["job"].update(cluster_job_id=command["job_id"], node_id=self.node_id)
        return result

    def artifact(self, command):
        from mn_molecular_design_studio.engine.workflow_engine.services.local_automation import resolve_job
        run = Path(resolve_job(command["payload"]["local_reference"])["run_dir"]).resolve()
        relative = command["payload"].get("path")
        if not isinstance(relative, str) or Path(relative).is_absolute():
            raise ValueError("Select a run-relative artifact path.")
        path = (run / relative).resolve()
        if not path.is_relative_to(run) or not path.is_file() or path.stat().st_size > bundles.MAX_BYTES:
            raise ValueError("This artifact is unavailable or exceeds the 512 MB download limit.")
        response = transport.node_request("/blobs", binary=path.read_bytes(), method="POST", timeout=90)
        if response["blob_id"] != bundles.digest_file(path):
            raise ValueError("The downloaded artifact failed its integrity check.")
        return response | {"name": path.name}

    def _poll_server(self, timeout=15, expected=None):
        reply = transport.node_request("/claim", method="POST", timeout=timeout)
        if expected and reply.get("server_node_id") != expected:
            raise transport.ClusterUnavailable("The configured peer address belongs to another workstation.")
        command = reply.get("command")
        if command:
            response, error = self.execute(command)
            transport.node_request("/commands/" + command["id"] + "/finish",
                                   {"lease_token": command["lease_token"], "response": response, "error": error}, method="POST", timeout=30)
        return bool(command)

    def poll(self):
        if not config.is_peer():
            return self._poll_server()
        settings = config.load()
        servers = [(settings["coordinator_url"], self.node_id), *[(p["url"],p["id"]) for p in settings.get("peers",[])]]
        addresses = [address for address,_ in servers]
        for index in range(len(addresses)):
            position = (self._poll_index + index) % len(addresses)
            address = addresses[position]
            if self._peer_retry.get(address, 0) > time.monotonic():
                continue
            try:
                with transport.command_server(address):
                    active = self._poll_server(timeout=3, expected=servers[position][1])
                self._peer_retry.pop(address, None)
                if active:
                    self._poll_index = (position + 1) % len(addresses)
                    return True
            except transport.ClusterUnavailable:
                self._peer_retry[address] = time.monotonic() + 15
        return False

    def run(self, *, once=False):
        if self.settings["mode"] == "disabled" or not self.settings["node_token"]:
            raise ValueError("Enroll and configure this workstation before starting its agent.")
        with (config.home() / "agent.lock").open("w") as lock:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as exc:
                raise RuntimeError("A workstation agent is already running on this machine.") from exc
            self.heartbeat()
            if once:
                self.poll()
                return
            def heartbeat_loop():
                while not self._stop.wait(15):
                    try:
                        self.heartbeat()
                    except Exception as exc:
                        log.warning("Catalog synchronization paused: %s", exc)
            thread = threading.Thread(target=heartbeat_loop, daemon=True)
            thread.start()
            signal.signal(signal.SIGTERM, lambda *_: self._stop.set())
            signal.signal(signal.SIGINT, lambda *_: self._stop.set())
            try:
                while not self._stop.is_set() and config.enabled():
                    try:
                        active = self.poll()
                    except Exception as exc:
                        log.warning("Coordinator connection paused: %s", exc)
                        active = False
                    self._stop.wait(0.2 if active else 2)
            finally:
                self._stop.set()
                thread.join(timeout=5)
