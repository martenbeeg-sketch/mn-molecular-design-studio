from __future__ import annotations

import hashlib
import json
import os
import tempfile
from pathlib import Path

from fastapi import APIRouter, Body, Depends, HTTPException, Query, Request
from fastapi.responses import FileResponse

from . import bundles, config, service
from .store import Store


def install(app, current_user, admin_user):
    router = APIRouter(prefix="/api/v1/cluster", tags=["workstations"])

    def coordinator():
        if config.load()["mode"] != "coordinator":
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
             user=Depends(current_user), store=Depends(coordinator)):
        return call(store.search, service.user_record(user), query=query, project_id=project_id, node_id=node_id,
                    statuses=statuses, tools=tools, limit=limit, offset=offset)

    @router.get("/jobs/{job_id}")
    def job(job_id: str, user=Depends(current_user), store=Depends(coordinator)):
        return call(store.job, job_id, service.user_record(user))

    @router.post("/jobs/{job_id}/actions", status_code=202)
    def action(job_id: str, payload: dict = Body(...), user=Depends(current_user), store=Depends(coordinator)):
        return call(store.action, job_id, payload.get("action"), payload.get("parameters") or {}, service.user_record(user))

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

    @router.post("/node/heartbeat")
    def heartbeat(payload: dict = Body(...), actor=Depends(node), store=Depends(coordinator)):
        call(store.heartbeat, actor["id"], payload.get("snapshot", {}), payload.get("jobs", []), payload.get("known_references"))
        from mn_molecular_design_studio.workspace import load_organization
        return {"status": "ok", "organization": load_organization(), "node_id": actor["id"], "reference_probes": store.reference_probes()}

    @router.post("/node/claim")
    def claim(actor=Depends(node), store=Depends(coordinator)):
        return {"command": call(store.claim, actor["id"])}

    @router.post("/node/commands/{command_id}/finish")
    def finish(command_id: str, payload: dict = Body(...), actor=Depends(node), store=Depends(coordinator)):
        call(store.finish, actor["id"], command_id, payload.get("lease_token"), payload.get("response"), payload.get("error"))
        return {"status": "ok"}

    app.include_router(router)
