"""One discoverable job-control contract for local runs and workstation jobs."""
from __future__ import annotations

from pathlib import Path
from typing import Literal

from fastapi import APIRouter, Depends, HTTPException, Query, Request
from pydantic import BaseModel, ConfigDict, Field


class JobAction(BaseModel):
    model_config = ConfigDict(extra="forbid")
    action: Literal["inspect", "artifact", "cancel", "hide", "restore", "retry", "change_gpu", "md_source_options", "pause", "resume"]
    parameters: dict = Field(default_factory=dict)


ACTION_PARAMETERS = {
    "inspect": {}, "artifact": {"path": {"type": "string", "minLength": 1}},
    "cancel": {}, "hide": {}, "restore": {}, "retry": {},
    "change_gpu": {"gpu_id": {"type": "integer", "minimum": 0}},
    "md_source_options": {}, "pause": {}, "resume": {"gpu_device": {"type": "string"}},
}


def install(app, current_user, *, resolve_job, call):
    from mn_molecular_design_studio.cluster import config, service
    from mn_molecular_design_studio.cluster.store import Store, TERMINAL
    from mn_molecular_design_studio.engine.workflow_engine.services import local_automation as local
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs
    router = APIRouter(prefix="/api/v1/executions", tags=["job control"])

    def coordinator_request(request, payload=None):
        from mn_molecular_design_studio.cluster import transport
        path = request.url.path + ("?" + request.url.query if request.url.query else "")
        bearer = request.headers.get("authorization", "").removeprefix("Bearer ")
        return call(transport.request, path, payload, token=bearer, method=request.method)

    @app.get("/api/v1/workstations", tags=["job control"])
    def workstations(request: Request, user=Depends(current_user)):
        if config.is_worker():
            from mn_molecular_design_studio.cluster import transport
            token = request.headers.get("authorization", "").removeprefix("Bearer ")
            return call(transport.request, "/api/v1/cluster/nodes", token=token)
        if config.enabled():
            return {"nodes": Store().nodes(), "local_node": config.identity()}
        from mn_molecular_design_studio.cluster.capabilities import resources
        scheduler_online, snapshot = call(resources)
        return {"nodes": [{"id": config.identity(), "name": config.load()["name"], "availability": "online",
                           "snapshot": {"resources": snapshot, "scheduler_online": scheduler_online,
                                        "accept_jobs": scheduler_online and config.load()["accept_jobs"]}}],
                "local_node": config.identity()}

    @router.get("")
    def list_executions(request: Request, project_id: str | None = None, distributed: bool = True, query: str = "",
                        limit: int = Query(50, ge=1, le=1000), offset: int = Query(0, ge=0),
                        include_hidden: bool = False, user=Depends(current_user)):
        if distributed and config.is_worker():
            return coordinator_request(request)
        if distributed and config.enabled():
            page = call(Store().search, service.user_record(user), project_id=project_id, query=query,
                        limit=limit, offset=offset, include_hidden=include_hidden)
            page["items"] = [{**row, "reference": "cluster:" + row["cluster_job_id"]} for row in page["items"]]
        else:
            from mn_molecular_design_studio.auth import accessible_project_ids
            from mn_molecular_design_studio.job_catalog import search_jobs
            page = call(search_jobs, project_id=project_id, query=query, limit=limit, offset=offset,
                        include_hidden=include_hidden, access_owner_id=user["id"],
                        accessible_project_ids=accessible_project_ids(user), access_admin=user.get("role") == "admin")
            page["items"] = [{**row, "reference": "local:" + row["task_group"] + ":" + row["run_id"]} for row in page["items"]]
        return page

    def resolve(reference, user):
        if reference.startswith("cluster:"):
            if not config.enabled():
                raise HTTPException(503, "Workstation sharing is disabled.")
            store = Store()
            return call(store.job, reference[8:], service.user_record(user)), store
        return resolve_job(reference, user), None

    def actions_for(row, store):
        active = row.get("status") not in TERMINAL | {"completed", "succeeded", "paused", "stopped"}
        if store:
            accepted = bool(row.get("local_reference"))
            enabled = {"inspect": accepted, "artifact": accepted, "md_source_options": accepted,
                       "hide": not row.get("hidden"), "restore": bool(row.get("hidden")),
                       "cancel": active,
                       "retry": row.get("status") == "failed" and not accepted and row.get("error") ==
                           "The workstation's engine images, code, or checkpoints no longer match this submission.",
                       "change_gpu": row.get("status") == "waiting_for_workstation" and not accepted}
        else:
            resumable = (Path(row["run_dir"]) / "worker_request.json").is_file()
            enabled = {"inspect": True, "hide": not row.get("hidden"), "restore": bool(row.get("hidden")),
                       "cancel": active or row.get("status") == "paused", "pause": active and resumable,
                       "resume": row.get("status") in {"paused", "stopped"} and resumable}
        return [{"action": name, "enabled": bool(available), "parameters_schema": {
                    "type": "object", "properties": ACTION_PARAMETERS[name], "additionalProperties": False,
                    "required": list(ACTION_PARAMETERS[name]) if name in {"artifact", "change_gpu"} else []},
                 "description": "Retry workstation acceptance after a compatibility failure; does not rerun a scientific calculation."
                    if name == "retry" else "Resource and job-state eligibility is checked again when requested."}
                for name, available in enabled.items()]

    @router.get("/{reference}")
    def inspect(reference: str, request: Request, user=Depends(current_user)):
        if reference.startswith("cluster:") and config.is_worker():
            return coordinator_request(request)
        row, store = resolve(reference, user)
        details = {"job": row} if store else call(local.inspect_job, row["task_group"] + "/" + row["run_id"], include_hidden=True)
        return {"reference": reference, "distributed": store is not None, **details,
                "actions": actions_for(row, store)}

    @router.get("/{reference}/actions")
    def available_actions(reference: str, request: Request, user=Depends(current_user)):
        if reference.startswith("cluster:") and config.is_worker():
            return coordinator_request(request)
        row, store = resolve(reference, user)
        return {"reference": reference, "items": actions_for(row, store)}

    @router.get("/{reference}/results")
    def results(reference: str, request: Request, user=Depends(current_user)):
        if reference.startswith("cluster:") and config.is_worker():
            return coordinator_request(request)
        row, store = resolve(reference, user)
        if store:
            with store.connect() as db:
                latest = db.execute("SELECT id FROM commands WHERE job_id=? AND kind='inspect' ORDER BY created_at DESC LIMIT 1",
                                    (reference[8:],)).fetchone()
            if not latest:
                return {"reference": reference, "inspection_required": True,
                        "inspection_url": "/api/v1/executions/" + reference + "/results"}
            operation = call(store.command, latest["id"], service.user_record(user))
            return {"reference": reference, "pending": operation["status"] not in {"completed", "failed", "cancelled"},
                    "command": operation, "poll_url": "/api/v1/cluster/commands/" + latest["id"]}
        return {"reference": reference, "pending": False,
                **call(local.job_results, row["task_group"] + "/" + row["run_id"], include_hidden=True)}

    @router.post("/{reference}/results", status_code=202)
    def inspect_results(reference: str, request: Request, user=Depends(current_user)):
        if reference.startswith("cluster:") and config.is_worker():
            return coordinator_request(request, {})
        row, store = resolve(reference, user)
        if not store:
            return results(reference, request, user)
        operation = call(store.action, reference[8:], "inspect", {}, service.user_record(user))
        return {"reference": reference, "pending": True, **operation,
                "poll_url": "/api/v1/cluster/commands/" + operation["command_id"]}

    @router.post("/{reference}/actions", status_code=202)
    def action(reference: str, payload: JobAction, request: Request, user=Depends(current_user)):
        if reference.startswith("cluster:") and config.is_worker():
            return coordinator_request(request, payload.model_dump())
        row, store = resolve(reference, user)
        contract = next((a for a in actions_for(row, store) if a["action"] == payload.action), None)
        if not contract or not contract["enabled"]:
            raise HTTPException(409, "This action is unavailable in the job's current state.")
        params = payload.parameters
        schema = contract["parameters_schema"]
        if set(params) - set(schema["properties"]) or set(schema["required"]) - set(params):
            raise HTTPException(422, "Use the parameters from this job's action schema.")
        if payload.action == "change_gpu" and (type(params["gpu_id"]) is not int or params["gpu_id"] < 0):
            raise HTTPException(422, "gpu_id must be a non-negative integer.")
        if payload.action == "artifact" and (not isinstance(params["path"], str) or not params["path"].strip()):
            raise HTTPException(422, "Provide an artifact path from job inspection.")
        if store:
            return call(store.action, reference[8:], payload.action, params, service.user_record(user))
        folder = Path(row["run_dir"])
        canonical = row["task_group"] + "/" + row["run_id"]
        if payload.action == "inspect":
            return call(local.inspect_job, canonical, include_hidden=True)
        if payload.action == "cancel":
            return call(local.cancel_job, canonical)
        if payload.action in {"hide", "restore"}:
            from mn_molecular_design_studio.job_catalog import set_job_hidden
            call(set_job_hidden, row["task_group"], row["run_id"], payload.action == "hide")
        elif payload.action == "pause":
            call(jobs.pause_job, folder, reason="Paused through Studio API")
        elif payload.action == "resume":
            if "gpu_device" in params and not isinstance(params["gpu_device"], str):
                raise HTTPException(422, "gpu_device must be a string.")
            call(jobs.resume_job, folder, **params)
        return {"reference": reference, "action": payload.action, "status": "accepted"}

    app.include_router(router)
