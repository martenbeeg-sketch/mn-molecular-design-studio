from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import tempfile
import uuid
from pathlib import Path

from fastapi import APIRouter, Body, Depends, HTTPException, Query, Request, Response
from fastapi.responses import FileResponse

from . import bundles, config, service
from .store import Store, TERMINAL


def install(app, current_user, admin_user):
    router = APIRouter(prefix="/api/v1/cluster", tags=["workstations"])

    def coordinator():
        if config.load()["mode"] not in {"coordinator", "peer"}:
            raise HTTPException(503, "This installation is not the coordinator.")
        return Store()

    def call(function, *args, **kwargs):
        try:
            return function(*args, **kwargs)
        except LookupError as exc:
            raise HTTPException(404, str(exc)) from exc
        except PermissionError as exc:
            raise HTTPException(403, str(exc)) from exc
        except (ValueError, KeyError, TypeError) as exc:
            raise HTTPException(422, str(exc)) from exc

    def node(request: Request, store=Depends(coordinator)):
        auth = request.headers.get("Authorization", "")
        if not auth.startswith("Bearer "):
            raise HTTPException(401, "A workstation credential is required.")
        return call(store.authenticate_node, auth[7:])

    @router.get("/nodes")
    def nodes(user=Depends(current_user), store=Depends(coordinator)):
        return {"nodes": store.nodes(), "coordinator_node": config.identity()}

    @router.post("/placement-options")
    def placement_options(payload: dict = Body(...), user=Depends(current_user), store=Depends(coordinator)):
        if set(payload) != {"requirements", "origin_node", "affinity_node"} or not isinstance(payload.get("requirements"), dict):
            raise HTTPException(422, "Provide workstation requirements, origin_node, and affinity_node.")
        node_ids = {node["id"] for node in store.nodes()}
        origin_node = str(payload.get("origin_node") or "")
        affinity_node = payload.get("affinity_node")
        if (origin_node not in node_ids or affinity_node is not None
                and (not isinstance(affinity_node, str) or affinity_node not in node_ids)):
            raise HTTPException(422, "The origin or data-owning workstation is not registered.")
        candidates = call(
            store.placement_options, payload["requirements"],
            origin_node=origin_node, affinity_node=affinity_node,
        )
        return {"candidates": candidates}

    @router.post("/nodes", status_code=201)
    def enroll(payload: dict = Body(...), user=Depends(admin_user), store=Depends(coordinator)):
        return call(store.enroll, payload.get("name"), payload.get("node_id")) | {"coordinator_url": config.load()["coordinator_url"]}

    @router.patch("/nodes/{node_id}")
    def node_settings(node_id: str, payload: dict = Body(...), user=Depends(admin_user), store=Depends(coordinator)):
        if set(payload) - {"disabled", "name"} or "disabled" in payload and type(payload["disabled"]) is not bool:
            raise HTTPException(422, "Use name and a boolean disabled setting.")
        with store.connect() as db:
            if not db.execute("SELECT 1 FROM nodes WHERE id=?", (node_id,)).fetchone():
                raise HTTPException(404, "Workstation not found.")
            if "disabled" in payload:
                db.execute("UPDATE nodes SET disabled=? WHERE id=?", (int(payload["disabled"]), node_id))
            if payload.get("name"):
                db.execute("UPDATE nodes SET name=? WHERE id=?", (str(payload["name"])[:120], node_id))
        return {"status": "updated"}

    @router.get("/organization")
    def organization(user=Depends(current_user), store=Depends(coordinator)):
        from mn_molecular_design_studio.workspace import load_organization
        from mn_molecular_design_studio.auth import accessible_project_ids
        org = load_organization()
        allowed = set(accessible_project_ids(user))
        return dict(org, projects=[p for p in org.get("projects", []) if p["id"] in allowed])

    @router.post("/projects", status_code=201)
    def create_project(payload: dict = Body(...), user=Depends(current_user), store=Depends(coordinator)):
        from mn_molecular_design_studio.workspace import create_project
        return call(create_project, payload.get("name"), payload.get("description", ""), payload.get("color", "#ef5350"), owner_id=user["id"])

    @router.get("/jobs")
    def jobs(query: str = "", project_id: str | None = None, node_id: str | None = None,
             statuses: list[str] = Query(default=[]), tools: list[str] = Query(default=[]),
             limit: int = Query(100, ge=1, le=1000), offset: int = Query(0, ge=0),
             include_hidden: bool = False,
             user=Depends(current_user), store=Depends(coordinator)):
        return call(store.search, service.user_record(user), query=query, project_id=project_id, node_id=node_id,
                    statuses=statuses, tools=tools, limit=limit, offset=offset, include_hidden=include_hidden)

    @router.get("/jobs/{job_id}")
    def job(job_id: str, user=Depends(current_user), store=Depends(coordinator)):
        return call(store.job, job_id, service.user_record(user))

    @router.post("/jobs/{job_id}/actions", status_code=202)
    def action(job_id: str, payload: dict = Body(...), user=Depends(current_user), store=Depends(coordinator)):
        return call(store.action, job_id, payload.get("action"), payload.get("parameters") or {}, service.user_record(user))

    @router.delete("/jobs/{job_id}/purge")
    def purge_job(job_id: str, response: Response, user=Depends(current_user), store=Depends(coordinator)):
        """Permanently remove a finished run from its owning workstation."""
        actor = service.user_record(user)
        existing = call(store.job, job_id, actor)
        if existing.get("status") not in TERMINAL:
            raise HTTPException(409, "Only finished jobs can be permanently removed.")
        context = call(store.purge_context, job_id, actor)
        job = context["job"]
        local_reference = str(job.get("local_reference") or "")

        if local_reference and job.get("node_id") != config.identity():
            operation = call(
                store.action,
                job_id,
                "purge",
                {
                    "accepted_command_id": context.get("command_id") or "",
                    "actor": {
                        "id": actor.get("id"),
                        "role": actor.get("role"),
                        "_cluster_project_ids": actor.get("project_ids") or actor.get("_cluster_project_ids", []),
                    },
                    "run_owner_id": job.get("owner_id"),
                    "run_project_id": job.get("project_id"),
                },
                actor,
            )
            response.status_code = 202
            return {
                "status": "purge_queued",
                "cluster_job_id": job_id,
                "command_id": operation["command_id"],
                "node_id": operation["node_id"],
            }

        task_group = run_id = ""
        run_files_removed = False
        database_result = {"removed_sources": 0, "removed_members": 0, "removed_groups": 0, "preserved_saved_selections": True}
        staging_removed = False

        if local_reference:
            if job.get("node_id") != config.identity():
                raise HTTPException(409, "This job is stored on another workstation; purge it from that workstation.")
            parts = local_reference.split("/")
            if len(parts) != 2 or any(part in {"", ".", ".."} for part in parts):
                raise HTTPException(422, "The job has an invalid local run reference.")
            task_group, run_id = parts
            if (job.get("task_group") and job["task_group"] != task_group
                    or job.get("run_id") and job["run_id"] != run_id):
                raise HTTPException(409, "The job's local run reference does not match its indexed identity.")

            from mn_molecular_design_studio.molecular_database_results import remove_run_records
            database_result = call(
                remove_run_records, task_group, run_id, user,
                owner_id=job.get("owner_id"), project_id=job.get("project_id"),
            )
            from mn_molecular_design_studio.engine.workflow_engine.core.jobs import delete_job_run
            try:
                call(delete_job_run, task_group, run_id)
                run_files_removed = True
            except FileNotFoundError:
                # A prior interrupted cleanup may already have removed the
                # folder. Keep the purge idempotent and finish coordinator cleanup.
                pass

        command_id = str(context.get("command_id") or "")
        if job.get("node_id") == config.identity() and command_id:
            try:
                command_id = str(uuid.UUID(command_id))
            except ValueError:
                command_id = ""
            if command_id:
                accepted_root = (config.home() / "accepted-inputs").resolve()
                accepted_dir = (accepted_root / command_id).resolve()
                try:
                    accepted_dir.relative_to(accepted_root)
                except ValueError as exc:
                    raise HTTPException(422, "The accepted-input path is invalid.") from exc
                if accepted_dir.is_symlink():
                    raise HTTPException(422, "The accepted-input path cannot be a symbolic link.")
                if accepted_dir.is_dir():
                    shutil.rmtree(accepted_dir)
                    staging_removed = True
                receipt = config.home() / "receipts" / f"{command_id}.json"
                if receipt.is_file() and not receipt.is_symlink():
                    receipt.unlink()

        removed = call(store.purge_metadata, job_id, actor)

        # Results indexing is file based, so reconciling the local catalog drops
        # the now-deleted run immediately. Input blobs are removed only when no
        # remaining submission or command refers to their digest.
        if task_group and run_id:
            from mn_molecular_design_studio.engine.workflow_engine.core.jobs import runs_root
            from mn_molecular_design_studio.job_catalog import reconcile
            reconcile(runs_root())
        digest = removed.get("unreferenced_bundle_id") or ""
        input_bundle_removed = False
        if re.fullmatch(r"[0-9a-f]{64}", digest):
            blob = config.home() / "transfer-cache" / digest
            if blob.is_file() and not blob.is_symlink():
                blob.unlink()
                input_bundle_removed = True

        return {
            "status": "purged",
            "cluster_job_id": removed["cluster_job_id"],
            "run_files_removed": run_files_removed,
            "staging_inputs_removed": staging_removed,
            "database": database_result,
            "removed_coordinator_commands": removed["removed_commands"],
            "removed_submission": removed["removed_submission"],
            "input_bundle_removed": input_bundle_removed,
        }

    @router.get("/commands/{command_id}")
    def operation(command_id: str, user=Depends(current_user), store=Depends(coordinator)):
        return call(store.command, command_id, service.user_record(user))

    @router.get("/commands/{command_id}/file")
    def operation_file(command_id: str, user=Depends(current_user), store=Depends(coordinator)):
        command = call(store.command, command_id, service.user_record(user))
        if command["status"] != "completed" or not command["response"].get("blob_id"):
            raise HTTPException(409, "This file transfer is not complete.")
        path = call(bundles.blob_path, command["response"]["blob_id"])
        if not path.is_file():
            raise HTTPException(410, "This temporary download expired. Request the original file again from its workstation.")
        return FileResponse(path, filename=Path(command["response"].get("name", "result")).name)

    @router.post("/submissions", status_code=202)
    def submit(payload: dict = Body(...), user=Depends(current_user), store=Depends(coordinator)):
        return call(service.accept_submission, payload, user)

    async def receive_blob(request, *, owner_id=None):
        digest = hashlib.sha256()
        size = 0
        fd, temporary_name = tempfile.mkstemp(dir=config.home(), prefix="upload-")
        temporary = Path(temporary_name)
        try:
            with os.fdopen(fd, "wb") as handle:
                async for chunk in request.stream():
                    size += len(chunk)
                    if size > bundles.MAX_BYTES:
                        raise HTTPException(413, "Transfer exceeds the 512 MB staging limit.")
                    digest.update(chunk)
                    handle.write(chunk)
            value = digest.hexdigest()
            temporary.replace(bundles.blob_path(value))
            if owner_id:
                with Store().connect() as db:
                    db.execute("CREATE TABLE IF NOT EXISTS blob_grants(digest TEXT,owner_id TEXT,PRIMARY KEY(digest,owner_id))")
                    db.execute("INSERT OR IGNORE INTO blob_grants VALUES (?,?)", (value, owner_id))
            return {"blob_id": value, "size_bytes": size}
        finally:
            temporary.unlink(missing_ok=True)

    @router.post("/blobs")
    async def user_blob(request: Request, user=Depends(current_user), store=Depends(coordinator)):
        return await receive_blob(request, owner_id=user["id"])

    @router.post("/node/blobs")
    async def node_blob(request: Request, actor=Depends(node)):
        return await receive_blob(request)

    @router.get("/node/blobs/{digest}")
    def node_download(digest: str, actor=Depends(node), store=Depends(coordinator)):
        with store.connect() as db:
            rows = db.execute("SELECT payload FROM commands WHERE node_id=? AND status IN ('queued','claimed')", (actor["id"],)).fetchall()
        if not any(json.loads(r["payload"]).get("bundle_id") == digest for r in rows):
            raise HTTPException(404, "Input bundle not assigned to this workstation.")
        path = call(bundles.blob_path, digest)
        if not path.is_file():
            raise HTTPException(404, "Input bundle unavailable.")
        return FileResponse(path)

    def peer_replica():
        if not config.is_peer():
            raise HTTPException(409, "Enable independent workstation sharing first.")
        from .replication import Replica
        return Replica()

    @router.post("/node/peer/hello")
    def peer_hello(payload: dict = Body(...), actor=Depends(node), replica=Depends(peer_replica), store=Depends(coordinator)):
        if payload.get("node_id") != actor["id"] or actor["id"] == config.identity():
            raise HTTPException(403, "A peer can publish only its own resource snapshot.")
        call(store.heartbeat, actor["id"], payload.get("snapshot", {}), [])
        own = next((n for n in store.nodes() if n["id"] == config.identity()), None)
        snapshot = dict(own["snapshot"]) if own else {}
        if not own or own["availability"] != "online":
            snapshot.update(accept_jobs=False, scheduler_online=False)
        probes = payload.get("reference_probes", [])
        if not isinstance(probes, list) or len(probes) > 10000 or any(not isinstance(p,str) for p in probes):
            raise HTTPException(422, "Invalid reference inventory request.")
        config.atomic_json(config.home() / ("peer-probes-" + actor["id"] + ".json"), {"references":probes})
        return {"node_id":config.identity(), "snapshot":snapshot, "epoch":replica.epoch, "reference_probes":store.reference_probes()}

    @router.get("/node/peer/records")
    def peer_records(response: Response, since: int = Query(0, ge=0), actor=Depends(node), replica=Depends(peer_replica)):
        response.headers["Cache-Control"] = "no-store"
        return call(replica.export, since)

    @router.post("/node/peer/records")
    def receive_records(payload: dict = Body(...), actor=Depends(node), replica=Depends(peer_replica)):
        return {"accepted":call(replica.receive, payload.get("records"))}

    @router.post("/node/heartbeat")
    def heartbeat(payload: dict = Body(...), actor=Depends(node), store=Depends(coordinator)):
        call(store.heartbeat, actor["id"], payload.get("snapshot", {}), payload.get("jobs", []), payload.get("known_references"))
        from mn_molecular_design_studio.workspace import load_organization
        return {"status": "ok", "organization": load_organization(), "node_id": actor["id"], "reference_probes": store.reference_probes()}

    @router.post("/node/claim")
    def claim(actor=Depends(node), store=Depends(coordinator)):
        result = {"command": call(store.claim, actor["id"])}
        if config.is_peer():
            result["server_node_id"] = config.identity()
        return result

    @router.post("/node/commands/{command_id}/finish")
    def finish(command_id: str, payload: dict = Body(...), actor=Depends(node), store=Depends(coordinator)):
        call(store.finish, actor["id"], command_id, payload.get("lease_token"), payload.get("response"), payload.get("error"))
        return {"status": "ok"}

    app.include_router(router)
