"""Outbound workstation agent; local schedulers continue owning computation."""
from __future__ import annotations

import fcntl
import hashlib
import json
import logging
import shutil
import signal
import threading
import time
import uuid
from pathlib import Path
from urllib.parse import urlsplit

from . import bundles, capabilities, config, context, progress as job_progress, transport
from mn_molecular_design_studio.native_filter_metrics import native_filter_evidence

log = logging.getLogger(__name__)


def _accepted_candidate_artifacts(run: Path, job: dict, artifacts: list[dict]) -> list[dict]:
    """Expose only campaign control files and files referenced by accepted rows."""
    prefix = f"{job.get('task_group') or ''}/{job.get('run_id') or ''}/".lstrip("/")
    run_path = Path(run).resolve()
    allowed = {
        prefix + name for name in (
            "metadata.json", "input.json",
            "artifacts/distributed_campaign/accepted_candidates.jsonl",
            "artifacts/rfd4_proteina/request.json",
            "artifacts/rfd4_proteina/native_command.json",
            "artifacts/rfd4_proteina/effective_config.yaml",
            "artifacts/rfd4_proteina/prepared_inputs.json",
        )
    }

    def collect(value):
        if isinstance(value, dict):
            for nested in value.values():
                collect(nested)
        elif isinstance(value, list):
            for nested in value:
                collect(nested)
        elif isinstance(value, str):
            if value.startswith(("runs:///", "run:///")):
                reference = value.split("///", 1)[1]
                if reference.startswith(prefix):
                    allowed.add(reference)
                return
            if value.startswith("artifacts/"):
                if value.startswith("artifacts/rfd4_proteina/input/"):
                    return
                allowed.add(prefix + value)
                return
            candidate_path = Path(value)
            if candidate_path.is_absolute():
                try:
                    relative = candidate_path.resolve(strict=False).relative_to(run_path)
                except ValueError:
                    return
                if relative.as_posix().startswith("artifacts/rfd4_proteina/input/"):
                    return
                allowed.add(prefix + relative.as_posix())

    accepted_path = run_path / "artifacts" / "distributed_campaign" / "accepted_candidates.jsonl"
    try:
        with accepted_path.open(encoding="utf-8") as handle:
            for line in handle:
                if line.strip():
                    collect(json.loads(line))
    except (OSError, UnicodeError, json.JSONDecodeError):
        pass
    result = []
    for artifact in artifacts:
        relative = str(artifact.get("relative_path") or artifact.get("path") or "")
        if relative in allowed:
            result.append(artifact)
    return result


def _candidate_summary(run: Path, job: dict, acceptance_policy: dict | None = None) -> dict:
    """Publish campaign counts only; full candidate records stay on this worker."""
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
    from mn_molecular_design_studio.target_identity import target_reference_type

    input_data = read_json(Path(run) / "input.json")
    metadata_data = read_json(Path(run) / "metadata.json")
    reference_type = target_reference_type(input_data, metadata_data, run_dir=run)
    path = Path(run) / "artifacts" / "normalized_candidates" / "candidates.jsonl"
    if not path.is_file():
        return {"schema_version": 1, "status": "unavailable", "candidate_count": None}
    counts = {"passed": 0, "failed": 0, "unavailable": 0}
    acceptance = {"accepted": 0, "rejected": 0, "unresolved": 0, "not_evaluated": 0}
    native_vs_acceptance = ({native: dict.fromkeys(acceptance, 0)
                             for native in ("passed", "failed", "unavailable")}
                            if acceptance_policy else None)
    candidate_count = 0
    eligible_count = 0
    structure_count = 0
    complex_count = 0
    target_systems = {}
    accepted_rows = []
    run_root = Path(run).resolve()

    def structure_file(value):
        from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path

        path_value = resolve_stored_path(value, run_dir=run_root, must_exist=True) if value else None
        if path_value is None:
            return None
        try:
            resolved = path_value.resolve(strict=True)
            resolved.relative_to(run_root)
        except (OSError, ValueError):
            return None
        if resolved.is_file() and resolved.name.casefold().endswith(
                (".pdb", ".ent", ".cif", ".mmcif", ".bcif", ".pdbqt", ".pdb.gz", ".cif.gz", ".mmcif.gz")):
            return resolved
        return None

    try:
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                row = json.loads(line)
                if not isinstance(row, dict) or not str(row.get("candidate_id") or "").strip():
                    raise ValueError("Normalized candidate row has no stable identity.")
                raw = row.get("raw_metadata")
                if not isinstance(raw, dict):
                    raw = {}
                    row["raw_metadata"] = raw
                complex_path = structure_file(row.get("complex_pdb"))
                candidate_structure = any(structure_file(row.get(key) or raw.get("structure_file"))
                                           for key in ("complex_pdb", "binder_pdb", "structure"))
                if candidate_structure:
                    structure_count += 1
                if complex_path is not None:
                    complex_count += 1
                target_chains = [str(value) for value in row.get("target_chains") or raw.get("target_chains") or []]
                sequences_by_chain = raw.get("sequences_by_chain") if isinstance(raw.get("sequences_by_chain"), dict) else {}
                target_sequences = [str(sequences_by_chain.get(chain) or "") for chain in target_chains
                                    if sequences_by_chain.get(chain)]
                if not target_sequences:
                    target_sequences = [str(component.get("sequence") or "") for component in raw.get("components") or []
                                        if isinstance(component, dict) and component.get("sequence")
                                        and str(component.get("role") or "").casefold() == "target"]
                normalized_target = ["".join(sequence.split()).upper() for sequence in target_sequences if sequence]
                if normalized_target:
                    target_digest = hashlib.sha256(":".join(sorted(normalized_target)).encode()).hexdigest()
                    identity = "sequence:" + target_digest
                    target_systems[identity] = {"identity": identity, "identity_kind": "sequence",
                                                "sequence_length": sum(map(len, normalized_target)),
                                                "reference_type": reference_type,
                                                "label": f"Sequence-defined target · {sum(map(len, normalized_target))} aa · {target_digest[:8]}"}
                decision = (raw.get("validation_acceptance")
                            if isinstance(raw.get("validation_acceptance"), dict) else {})
                eligible = (
                    raw.get("candidate_pool_level") != "prefilter_generated"
                    and raw.get("result_kind") not in {"generation_only", "generator_trajectory"}
                    and str(row.get("stage") or "") not in {"generation.backbone", "generation_backbone", "generator_trajectory"}
                )
                if eligible:
                    eligible_count += 1
                    filter_value = native_filter_evidence(row)["value"]
                    counts["passed" if filter_value is True else "failed" if filter_value is False else "unavailable"] += 1
                    if acceptance_policy:
                        from mn_molecular_design_studio.engine.workflow_engine.workflows.distributed_campaign import (
                            _native_acceptance_geometry,
                            _stage_selectors,
                        )

                        _native_acceptance_geometry(run, row, acceptance_policy)
                        decision = _stage_selectors(row, row, {"acceptance_policy": acceptance_policy}).get(
                            "validation_acceptance", {}
                        )
                        raw["validation_acceptance"] = decision
                        status = str(decision.get("status") or "unresolved")
                        status = status if status in acceptance else "not_evaluated"
                        acceptance[status] += 1
                        native_status = ("passed" if filter_value is True else
                                         "failed" if filter_value is False else "unavailable")
                        native_vs_acceptance[native_status][status] += 1
                        if status == "accepted":
                            accepted_rows.append(row)
                        candidate_count += 1
                        continue
                decision_status = str(decision.get("status") or "not_evaluated")
                acceptance[decision_status if decision_status in acceptance else "not_evaluated"] += 1
                candidate_count += 1
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as exc:
        return {"schema_version": 1, "status": "invalid", "candidate_count": None, "error": str(exc)}
    if acceptance_policy:
        accepted_path = Path(run) / "artifacts" / "distributed_campaign" / "accepted_candidates.jsonl"
        accepted_path.parent.mkdir(parents=True, exist_ok=True)
        temporary = accepted_path.with_suffix(".jsonl.tmp")
        temporary.write_text(
            "".join(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n" for row in accepted_rows),
            encoding="utf-8",
        )
        temporary.replace(accepted_path)
    if candidate_count and not target_systems:
        from mn_molecular_design_studio.target_identity import prediction_input_systems

        target_systems.update({system["identity"]: system for system in prediction_input_systems(input_data)})
    summary = {
        "schema_version": 1,
        "status": "complete",
        "candidate_count": candidate_count,
        "eligible_count": eligible_count,
        "structure_count": structure_count,
        "complex_count": complex_count,
        "target_systems": sorted(target_systems.values(), key=lambda item: item["identity"]),
        "target_reference_type": reference_type,
        "native_filter_counts": counts,
        "acceptance_counts": acceptance,
    }
    if acceptance_policy:
        summary["acceptance_evaluation"] = "worker_local_saved_evidence"
        summary["native_filter_vs_acceptance"] = native_vs_acceptance
    return summary


class Agent:
    def __init__(self):
        self.settings = config.load()
        self.node_id = config.identity()
        self._inventory_time = 0
        self._inventory = {}
        self._published = {}
        self._artifact_signatures = {}
        self._artifact_manifests = {}
        self._artifact_manifest_counts = {}
        self._artifact_manifest_complete = {}
        self._reference_hash_cache = {}
        self._known_references = set()
        self._reference_probes = []
        self._stop = threading.Event()
        self._peer_retry = {}
        self._poll_index = 0
        self._processed_purge_receipts = set()
        self._library_transfer_threads = {}
        self._library_transfer_lock = threading.Lock()

    def snapshot(self):
        if time.monotonic() - self._inventory_time > 120 or not self._inventory:
            self._inventory = {"code_revision": capabilities.code_revision(), "images": capabilities.images()}
            self._inventory_time = time.monotonic()
        online, resources = capabilities.resources()
        storage_online = self.storage_online()
        return self._inventory | {"scheduler_online": online, "resources": resources,
                                  "references": capabilities.reference_inventory(self._reference_probes),
                                  "reference_hashes": capabilities.reference_fingerprints(
                                      self._reference_probes, self._reference_hash_cache,
                                  ),
                                  "accept_jobs": config.load()["accept_jobs"] and storage_online,
                                  "storage_available": storage_online}

    def _refresh_reference_probes(self, store=None):
        """Refresh model-path probes from local requests and peer requests."""
        probes = set()
        if store is None and config.is_peer():
            from .store import Store
            store = Store()
        if store is not None:
            try:
                probes.update(store.reference_probes())
            except Exception:
                # A temporary catalog lock must not suppress liveness. Keep the
                # last known probes and retry on the next pulse.
                probes.update(self._reference_probes)
        for path in config.home().glob("peer-probes-*.json"):
            try:
                values = json.loads(path.read_text()).get("references", [])
                if isinstance(values, list):
                    probes.update(value for value in values if isinstance(value, str))
            except (OSError, ValueError):
                continue
        self._reference_probes = sorted(probes)
        return self._reference_probes

    def _restore_published_signatures(self, store):
        """Avoid replaying the whole local Results catalog after an agent restart."""
        if self._published:
            return
        try:
            with store.connect() as db:
                rows = db.execute(
                    "SELECT local_ref,json_extract(payload,'$._signature') AS signature "
                    "FROM jobs WHERE node_id=? AND local_ref IS NOT NULL "
                    "AND json_type(payload,'$._signature')='array'",
                    (self.node_id,),
                ).fetchall()
            for row in rows:
                try:
                    signature = json.loads(row["signature"])
                except (TypeError, ValueError):
                    continue
                if isinstance(signature, list):
                    self._published[str(row["local_ref"])] = tuple(signature)
        except Exception:
            # A catalog lock or an older database schema should delay this
            # optimization, not prevent the workstation agent from starting.
            log.warning("Could not restore the local job sync cursor; will reconcile the catalog normally.")

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
                result_data = read_json(run / "result.json")
                metrics = result_data.get("metrics", {})
                acceptance_policy = metadata.get("distributed_acceptance_policy")
                output_policy = metadata.get("distributed_output_policy")
                compact_candidates = _candidate_summary(
                    run, row, acceptance_policy,
                ) if row.get("status") in {"completed", "succeeded"} else None
                if output_policy == "accepted_candidates_only":
                    metrics = {}
                jid = (
                    metadata.get("cluster_job_id")
                    or marker.get("cluster_job_id")
                    or metadata.get("cluster_submission_id")
                    or marker.get("cluster_submission_id")
                    or config.job_id(self.node_id, reference)
                )
                if self._artifact_signatures.get(reference) != artifact_signature:
                    manifest = job_catalog.job_artifacts(row["task_group"], row["run_id"])
                    if output_policy == "accepted_candidates_only":
                        manifest = _accepted_candidate_artifacts(run, row, manifest)
                    self._artifact_manifest_counts[reference] = len(manifest)
                    self._artifact_manifest_complete[reference] = len(manifest) <= 1000
                    self._artifact_manifests[reference] = manifest[:1000]
                    self._artifact_signatures[reference] = artifact_signature
                progress = job_progress.bindcraft2_progress(run)
                item = dict(row, cluster_job_id=jid, metrics=metrics, progress=progress,
                            artifacts=self._artifact_manifests.get(reference, []),
                            artifact_count=self._artifact_manifest_counts.get(reference, 0),
                            artifact_manifest_complete=self._artifact_manifest_complete.get(reference, False),
                            location=str(run), last_reported_at=time.time(), _signature=list(signature))
                if compact_candidates is not None:
                    item["candidate_summary"] = compact_candidates
                item["worker_metadata"] = {
                    "resource_allocation": metadata.get("resource_allocation") or {},
                    "scheduler_worker_id": metadata.get("scheduler_worker_id") or metadata.get("worker_id"),
                    "worker_started_at": metadata.get("worker_started_at"),
                    "completed_at": metadata.get("completed_at"),
                    "software_environment": metadata.get("software_environment") or {},
                }
                if self.settings.get("mode") == "worker":
                    outputs = result_data.get("outputs")
                    library_id = outputs.get("library_id") if isinstance(outputs, dict) else None
                    if library_id:
                        from mn_molecular_design_studio.external_libraries import LibraryCatalog
                        item["external_library_record"] = LibraryCatalog().get(library_id)
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
        self._schedule_library_transfers()

    def peer_heartbeat(self):
        from .store import Store
        from .replication import Replica, _changed_peer_jobs, capture_identity, sync_peer, locked
        store = Store()
        self._restore_published_signatures(store)
        # Reconcile completed purges before rebuilding the catalog or publishing
        # the peer mirror. Otherwise the stale mirrored row can be sent again.
        self._reconcile_purged_peer_jobs()
        # Required model probes must be published before catalog reconciliation.
        # On a workstation with a large historical Results catalog, catalog()
        # can take minutes on its first pass. Waiting for it here leaves the
        # node's capability snapshot without the very checkpoint paths needed
        # by newly submitted work, even though those files are installed.
        self._refresh_reference_probes(store)
        snapshot = self.snapshot()
        # Publish current capabilities and dispatch compatible submissions
        # before the potentially expensive catalog scan and peer replication.
        # Molecular artifacts remain on their producing workstation; this
        # heartbeat contains only capabilities and the normal compact job
        # metadata used by the existing cluster controller.
        store.heartbeat(self.node_id, snapshot, [])

        jobs = self.catalog() if self.storage_online() else []
        peer_jobs = _changed_peer_jobs(jobs, self._published)
        for index in range(0, max(1, len(jobs)), 1000):
            batch = jobs[index:index+1000]
            known = sorted(self._known_references) if index + 1000 >= len(jobs) and self.storage_online() else None
            store.heartbeat(self.node_id, snapshot, batch, known)
            for row in batch:
                self._published[row["task_group"] + "/" + row["run_id"]] = tuple(row["_signature"])
        with locked():
            capture_identity()
            replica = Replica()
            replica.compact_legacy_job_payloads()
            replica.capture_jobs()
        for peer in config.load().get("peers", []):
            try:
                # Job status/count summaries are an immediate path around the
                # durable replication cursor, but unchanged catalog rows do
                # not need to be resent every 15-second heartbeat.
                sync_peer(peer, snapshot, peer_jobs)
            except Exception as exc:
                log.warning("Workstation %s synchronization paused: %s", peer.get("name",peer["id"]), exc)
        config.atomic_json(config.home() / "agent-status.json", {"last_heartbeat":time.time(), "node_id":self.node_id, "online":True})
        store.cleanup_transfers()
        self._schedule_library_transfers()

    def publish_liveness(self):
        """Refresh local and remote node status without scanning or syncing jobs."""
        settings = config.load()
        store = None
        previous_local = {}
        if settings["mode"] == "peer":
            from .store import Store
            store = Store()
            self._refresh_reference_probes(store)
            previous_local = next(
                (node.get("snapshot", {}) for node in store.nodes() if node["id"] == self.node_id),
                {},
            )
        snapshot = self.snapshot()
        if settings["mode"] == "peer":
            changed = any(
                previous_local.get(key) != snapshot.get(key)
                for key in ("images", "references", "reference_hashes", "code_revision", "accept_jobs", "scheduler_online")
            )
            store.pulse(self.node_id, snapshot)
            for peer in settings.get("peers", []):
                try:
                    response = transport.node_request(
                        "/peer/hello", {
                            "node_id": self.node_id,
                            "snapshot": snapshot,
                            "reference_probes": self._reference_probes,
                            "jobs": [],
                        }, method="POST",
                        base_url=peer["url"], timeout=4,
                    )
                    if response.get("node_id") != peer["id"]:
                        raise transport.ClusterUnavailable(
                            "The configured peer address belongs to another workstation."
                        )
                    previous_peer = next(
                        (node.get("snapshot", {}) for node in store.nodes() if node["id"] == peer["id"]),
                        {},
                    )
                    peer_snapshot = response.get("snapshot") or {}
                    changed = changed or any(
                        previous_peer.get(key) != peer_snapshot.get(key)
                        for key in ("images", "references", "reference_hashes", "code_revision", "accept_jobs", "scheduler_online")
                    )
                    store.pulse(peer["id"], peer_snapshot)
                except Exception as exc:
                    log.warning(
                        "Workstation %s liveness check paused: %s",
                        peer.get("name", peer["id"]), exc,
                    )
            if changed:
                # Capability changes can make a queued model workflow eligible.
                # Avoid a full scheduling pass for routine resource telemetry.
                store.dispatch()
        else:
            response = transport.node_request(
                "/pulse", {"snapshot": snapshot}, method="POST", timeout=10,
            )
            if response.get("node_id") != self.node_id:
                raise RuntimeError(
                    "This credential belongs to another workstation. Check the enrollment file."
                )
        config.atomic_json(config.home() / "agent-status.json", {
            "last_heartbeat": time.time(), "node_id": self.node_id, "online": True,
        })

    def _heartbeat_interval(self):
        try:
            # Keep direct pulses comfortably below the 45-second unreachable
            # threshold even if an older config requests a slower interval.
            return max(5, min(15, int(config.load().get("heartbeat_seconds", 15))))
        except (TypeError, ValueError):
            return 15

    def _liveness_loop(self):
        while not self._stop.is_set():
            try:
                self.publish_liveness()
            except Exception as exc:
                log.warning("Workstation liveness update paused: %s", exc)
            if self._stop.wait(self._heartbeat_interval()):
                break

    def _schedule_library_transfers(self):
        """Retry finished external-library output copies after each successful sync."""
        from mn_molecular_design_studio import job_catalog
        from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json

        try:
            from mn_molecular_design_studio.runtime import runs_root
            page_offset = 0
            while True:
                page = job_catalog.search_jobs(access_admin=True, limit=1000, offset=page_offset)
                for row in page["items"]:
                    if row.get("tool") not in {"library_preparation", "library_centroids"}:
                        continue
                    run = Path(row["run_dir"])
                    metadata = read_json(run / "metadata.json")
                    result_data = read_json(run / "result.json")
                    outputs = result_data.get("outputs") if isinstance(result_data.get("outputs"), dict) else {}
                    manifest = outputs.get("feature_manifest") or outputs.get("centroid_manifest")
                    library_id = outputs.get("library_id")
                    organizer_node = metadata.get("organizer_node_id")
                    destination = metadata.get("organizer_output_directory")
                    job_id = str(metadata.get("cluster_job_id") or row.get("cluster_job_id") or "")
                    if (not manifest or not library_id or not organizer_node or not destination
                            or organizer_node == self.node_id or not job_id
                            or str(row.get("status") or "").casefold() not in {"completed", "succeeded"}):
                        continue
                    receipt_dir = config.home() / "library-transfer-receipts"
                    receipt_dir.mkdir(parents=True, exist_ok=True, mode=0o700)
                    receipt = receipt_dir / f"{job_id}.json"
                    if receipt.exists():
                        continue
                    with self._library_transfer_lock:
                        active = self._library_transfer_threads.get(job_id)
                        if active and active.is_alive():
                            continue
                        thread = threading.Thread(
                            target=self._transfer_library_outputs,
                            args=(job_id, str(library_id), Path(manifest), str(organizer_node), str(destination), receipt),
                            name=f"library-transfer-{job_id[:8]}", daemon=True,
                        )
                        self._library_transfer_threads[job_id] = thread
                        thread.start()
                page_offset += len(page["items"])
                if page_offset >= page["total"] or not page["items"]:
                    break
        except Exception:
            log.exception("Could not schedule completed library output transfers")

    def _transfer_library_outputs(self, job_id, library_id, manifest_path, organizer_node, destination, receipt):
        from mn_molecular_design_studio.external_libraries import FEATURE_SCHEMA, CENTROID_SCHEMA, read_feature_manifest, digest_file
        import json

        try:
            settings = config.load()
            base_url = None
            if settings["mode"] == "peer":
                peer = next((item for item in settings.get("peers", []) if item.get("id") == organizer_node), None)
                if peer is None:
                    raise ValueError("The originating organizer is not configured as a trusted peer.")
                base_url = peer["url"]
            elif settings["mode"] not in {"worker", "coordinator"}:
                raise ValueError("No organizer route is configured for the completed library output.")
            elif organizer_node == self.node_id:
                raise ValueError("A remote output transfer cannot target the compute workstation itself.")

            manifest_path = Path(manifest_path).resolve(strict=True)
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            if manifest.get("schema_version") == FEATURE_SCHEMA:
                manifest_path, manifest = read_feature_manifest(manifest_path)
            elif manifest.get("schema_version") != CENTROID_SCHEMA or manifest.get("status") != "completed":
                raise ValueError("The completed library manifest type is unsupported.")
            if manifest.get("library_id") != library_id:
                raise ValueError("The completed fingerprint manifest belongs to a different library.")
            output_dir = manifest_path.parent.resolve()
            files = []

            def add_file(name, expected_sha=None, expected_size=None):
                path = (output_dir / name).resolve(strict=True)
                if not path.is_relative_to(output_dir) or not path.is_file():
                    raise ValueError("A generated library output points outside its manifest directory.")
                size = path.stat().st_size
                digest = digest_file(path)
                if expected_sha and digest != expected_sha or expected_size is not None and size != expected_size:
                    raise ValueError(f"Generated output {name} changed after fingerprint preparation.")
                files.append({"path": name, "size_bytes": size, "sha256": digest})

            add_file("manifest.json")
            if manifest.get("schema_version") == FEATURE_SCHEMA:
                rejection = manifest.get("rejections", {})
                if rejection.get("path") != "rejected.parquet":
                    raise ValueError("The rejected-row manifest path is unsupported.")
                add_file(rejection["path"], rejection.get("sha256"))
                for shard in manifest.get("shards", []):
                    add_file(shard["path"], shard.get("sha256"), shard.get("size_bytes"))
            else:
                for item in manifest.get("files", []):
                    add_file(item["path"], item.get("sha256"))

            prefix = f"/api/v1/cluster/node/library-transfers/{job_id}"
            response = transport.node_request(prefix + "/start", {"library_id": library_id, "files": files},
                                               method="POST", timeout=60, base_url=base_url)
            if response.get("status") == "waiting_for_space":
                raise ValueError(response.get("error") or "Waiting for organizer disk space.")
            offsets = {item["path"]: item["offset"] for item in response.get("files", [])}
            for item in files:
                path = output_dir / item["path"]
                offset = int(offsets.get(item["path"], 0))
                while offset < item["size_bytes"]:
                    with path.open("rb") as stream:
                        stream.seek(offset)
                        block = stream.read(min(8 * 1024 * 1024, item["size_bytes"] - offset))
                    if not block:
                        raise ValueError(f"Generated output {item['path']} ended before its declared size.")
                    try:
                        chunk_result = transport.node_request(
                            prefix + "/files/" + item["path"] + f"?offset={offset}",
                            binary=block, method="PUT", timeout=180, base_url=base_url,
                        )
                        offset = int(chunk_result["offset"])
                    except ValueError:
                        current = transport.node_request(prefix, base_url=base_url, timeout=30)
                        refreshed = {entry["path"]: entry["offset"] for entry in current.get("files", [])}
                        new_offset = int(refreshed.get(item["path"], 0))
                        if new_offset <= offset:
                            raise
                        offset = new_offset
            completion = transport.node_request(prefix + "/complete", method="POST", timeout=180, base_url=base_url)
            if completion.get("status") != "completed":
                raise ValueError("The organizer did not confirm the completed fingerprint output transfer.")
            config.atomic_json(receipt, {"status": "completed", "job_id": job_id,
                                         "library_id": library_id, "completed_at": time.time()})
        except Exception as exc:
            log.warning("Library output transfer for job %s will retry: %s", job_id, exc)

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
            # Purging is idempotent: if a crash happens after deleting the run
            # but before acknowledging the command, retrying finishes cleanup.
            if command["kind"] != "purge":
                return {}, "Submission state is uncertain after an agent restart. Reconcile the existing local run before retrying."
        config.atomic_json(receipt, {"state": "accepting", "job_id": command["job_id"], "command_id": command["id"]})
        response, error = {}, None
        try:
            if command["kind"] in {"workflow", "pipeline"}:
                response = self.submit(command)
            elif command["kind"] == "inspect":
                from mn_molecular_design_studio.engine.workflow_engine.services.local_automation import inspect_job, job_results
                ref = command["payload"]["local_reference"]
                # The coordinator has already authorized this source job. Run
                # results also include native batches hidden from normal lists.
                response = inspect_job(ref, include_hidden=True) | {"results": job_results(ref, include_candidates=True, include_hidden=True)}
            elif command["kind"] == "md_source_options":
                from mn_molecular_design_studio.auth import list_users
                from mn_molecular_design_studio.engine.workflow_engine.services.local_automation import resolve_job
                from mn_molecular_design_studio.md_analysis_ui import md_source_options_for_job
                actor = next((u for u in list_users() if u["id"] == command["owner_id"] and not u.get("disabled")), None)
                if actor is None:
                    raise PermissionError("The MD source request's account is unavailable.")
                source_job = resolve_job(command["payload"]["local_reference"])
                response = md_source_options_for_job(source_job, actor)
            elif command["kind"] == "cancel":
                from mn_molecular_design_studio.engine.workflow_engine.services.local_automation import cancel_job
                response = {"job": cancel_job(command["payload"]["local_reference"], reason="Cancelled through the shared workstation catalog")}
            elif command["kind"] == "artifact":
                response = self.artifact(command)
            elif command["kind"] == "collect_accepted_designs":
                response = self.collect_accepted_designs(command)
            elif command["kind"] == "purge":
                response = self.purge(command)
            else:
                raise ValueError("Unsupported workstation command.")
        except Exception as exc:
            error = str(exc)
            log.exception("Workstation operation failed: %s", command["id"])
        config.atomic_json(receipt, {
            "state": "finished", "job_id": command["job_id"],
            "command_id": command["id"], "kind": command["kind"],
            "response": response, "error": error,
        })
        if command["kind"] == "purge" and not error:
            self._reconcile_purged_peer_jobs()
        return response, error

    def _reconcile_purged_peer_jobs(self):
        """Publish tombstones for jobs whose remote purge has completed.

        In peer mode the synchronized Store mirror can outlive the run folder.
        If it is left in place after a purge, ``capture_jobs`` republishes it
        and the coordinator recreates the supposedly deleted Results row.
        Receipts also let a restarted agent repair mirrors from older purges.
        """
        if not config.is_peer():
            return
        receipts_dir = config.home() / "receipts"
        purge_receipts = []
        for receipt_path in receipts_dir.glob("*.json"):
            if receipt_path.name in self._processed_purge_receipts:
                continue
            try:
                receipt = json.loads(receipt_path.read_text())
            except (OSError, json.JSONDecodeError):
                continue
            if not isinstance(receipt, dict):
                continue
            response = receipt.get("response")
            if (receipt.get("state") != "finished" or receipt.get("error")
                    or not isinstance(response, dict) or response.get("status") != "purged"):
                continue
            purge_receipts.append((receipt_path.name, str(receipt.get("job_id") or ""),
                                   str(receipt.get("command_id") or receipt_path.stem)))
        if not purge_receipts:
            return

        from .replication import Replica, capture_identity, locked
        from .store import Store

        processed_receipts = set()
        with Store().connect() as db:
            for receipt_name, job_id, command_id in purge_receipts:
                if not job_id:
                    command = db.execute("SELECT job_id FROM commands WHERE id=?", (command_id,)).fetchone()
                    job_id = str(command["job_id"] or "") if command else ""
                if not job_id:
                    # Keep the receipt eligible for a later retry in case its
                    # synchronized command record has not arrived yet.
                    continue
                row = db.execute("SELECT node_id,local_ref FROM jobs WHERE id=?", (job_id,)).fetchone()
                if row and row["node_id"] == self.node_id and row["local_ref"]:
                    db.execute("DELETE FROM commands WHERE job_id=?", (job_id,))
                    db.execute("DELETE FROM submissions WHERE job_id=?", (job_id,))
                    db.execute("DELETE FROM jobs WHERE id=?", (job_id,))
                processed_receipts.add(receipt_name)

        with locked():
            capture_identity()
            Replica().capture_jobs()
        self._processed_purge_receipts.update(processed_receipts)

    def purge(self, command):
        """Remove a terminal run and its run-derived records on this workstation."""
        payload = command.get("payload") or {}
        reference = str(payload.get("local_reference") or "")
        parts = reference.split("/")
        if len(parts) != 2 or any(part in {"", ".", ".."} or "\\" in part for part in parts):
            raise ValueError("The job has an invalid local run reference.")
        task_group, run_id = parts

        from mn_molecular_design_studio.engine.workflow_engine.core.jobs import delete_job_run, runs_root
        root = Path(runs_root()).expanduser().resolve()
        run_dir = (root / task_group / run_id).resolve()
        run_dir.relative_to(root)
        run_files_removed = False
        if run_dir.exists():
            if not run_dir.is_dir():
                raise ValueError("The job run path is not a directory.")
            try:
                metadata = json.loads((run_dir / "metadata.json").read_text())
            except (OSError, ValueError):
                metadata = {}
            status = str(metadata.get("status") or "").casefold()
            if status not in {"completed", "succeeded", "cancelled", "failed"}:
                raise ValueError("Only finished runs can be permanently removed.")

        actor = dict(payload.get("actor") or {})
        if not actor.get("id"):
            raise ValueError("The purge request is missing its authorized user.")
        if "_cluster_project_ids" not in actor and isinstance(actor.get("project_ids"), list):
            actor["_cluster_project_ids"] = actor["project_ids"]
        from mn_molecular_design_studio.molecular_database_results import remove_run_records
        database_result = remove_run_records(
            task_group,
            run_id,
            actor,
            owner_id=payload.get("run_owner_id"),
            project_id=payload.get("run_project_id"),
        )

        if run_dir.exists():
            delete_job_run(task_group, run_id)
            run_files_removed = True

        accepted_id = str(payload.get("accepted_command_id") or "")
        staging_removed = False
        if accepted_id:
            try:
                accepted_id = str(uuid.UUID(accepted_id))
            except ValueError as exc:
                raise ValueError("The accepted-input reference is invalid.") from exc
            accepted_root = (config.home() / "accepted-inputs").resolve()
            accepted_dir = accepted_root / accepted_id
            if accepted_dir.is_symlink():
                raise ValueError("The accepted-input path cannot be a symbolic link.")
            try:
                accepted_dir.resolve().relative_to(accepted_root)
            except ValueError as exc:
                raise ValueError("The accepted-input path is invalid.") from exc
            if accepted_dir.is_dir():
                shutil.rmtree(accepted_dir)
                staging_removed = True
            receipt = config.home() / "receipts" / f"{accepted_id}.json"
            if receipt.is_file() and not receipt.is_symlink():
                receipt.unlink()

        from mn_molecular_design_studio.job_catalog import reconcile
        reconcile(root)
        self._reconcile_purged_peer_jobs()
        return {
            "status": "purged",
            "run_files_removed": run_files_removed,
            "staging_inputs_removed": staging_removed,
            "database": database_result,
        }

    def submit(self, command):
        from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import submit_workflow_request, validate_workflow_request
        from mn_molecular_design_studio.workspace import mark_studio_job
        from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, write_json
        payload = command["payload"]
        needs = payload.get("requirements", {})
        snapshot = self.snapshot()
        from .store import Store
        if (payload.get("allow_code_revision_refresh") is True
                and command.get("node_id") == payload.get("origin_node")):
            # The coordinator permits this only for jobs explicitly pinned to
            # their submitting workstation. Validate against the current local
            # code while retaining strict image and checkpoint checks.
            needs = dict(needs, _allow_code_revision_refresh=True)
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
        metadata.update(cluster_job_id=command["job_id"], execution_node_id=self.node_id,
                        organizer_node_id=command.get("payload", {}).get("origin_node"),
                        organizer_output_directory=command.get("payload", {}).get("organizer_output_directory"))
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
        if not path.is_relative_to(run):
            raise ValueError("The requested artifact path is outside its source job.")
        if not path.is_file():
            raise FileNotFoundError("Missing source artifact: the requested file is absent from its source job.")
        if path.stat().st_size > bundles.MAX_BYTES:
            raise ValueError("The requested source artifact exceeds the 512 MB download limit.")
        response = transport.node_request("/blobs", binary=path.read_bytes(), method="POST", timeout=90)
        if response["blob_id"] != bundles.digest_file(path):
            raise ValueError("The downloaded artifact failed its integrity check.")
        return response | {"name": path.name}

    def collect_accepted_designs(self, command):
        """Transfer only explicitly accepted candidates and their structures."""
        from mn_molecular_design_studio.accepted_design_import import build_accepted_design_bundle
        from mn_molecular_design_studio.engine.workflow_engine.services.local_automation import resolve_job

        reference = command["payload"]["local_reference"]
        job = resolve_job(reference)
        run_dir = Path(job["run_dir"]).resolve()
        payload, manifest = build_accepted_design_bundle(
            run_dir, node_id=self.node_id, validation_job=command["job_id"],
        )
        response = transport.node_request("/blobs", binary=payload, method="POST", timeout=180)
        if response.get("blob_id") != hashlib.sha256(payload).hexdigest():
            raise ValueError("The accepted-design export failed its SHA-256 integrity check.")
        return response | {
            "name": f"accepted-designs-{job.get('run_id') or run_dir.name}.zip",
            "candidate_count": manifest["candidate_count"],
            "accepted_candidate_ids": manifest["accepted_candidate_ids"],
        }

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
            if once:
                self.heartbeat()
                self.poll()
                return

            def catalog_sync_loop():
                while not self._stop.is_set():
                    try:
                        self.heartbeat()
                    except Exception as exc:
                        log.warning("Catalog synchronization paused: %s", exc)
                    if self._stop.wait(self._heartbeat_interval()):
                        break

            # Liveness has its own small request path, so a full job/catalog
            # reconciliation cannot make a healthy workstation look offline.
            liveness_thread = threading.Thread(
                target=self._liveness_loop, name="studio-agent-liveness", daemon=True,
            )
            sync_thread = threading.Thread(
                target=catalog_sync_loop, name="studio-agent-catalog-sync", daemon=True,
            )
            liveness_thread.start()
            sync_thread.start()
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
                liveness_thread.join(timeout=5)
                sync_thread.join(timeout=5)
