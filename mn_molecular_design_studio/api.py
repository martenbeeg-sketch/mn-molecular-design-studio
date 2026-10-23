"""Local REST API for Molecular Design Studio workflows and jobs."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path
from typing import Any, Callable

from mn_molecular_design_studio.workspace import ensure_studio_workspace, mark_studio_job, validate_project_id

ensure_studio_workspace()

from fastapi import Body, Depends, FastAPI, HTTPException, Query, status
from fastapi.security import HTTPAuthorizationCredentials, HTTPBearer
from pydantic import BaseModel

from mn_molecular_design_studio.auth import (
    accessible_project_ids,
    authenticate,
    change_password,
    development_user,
    is_development_mode,
    require_job_access,
    require_project_access,
    resolve_session,
    user_input_root,
)
from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import (
    is_operational_path_key,
    is_portable_path,
    resolve_stored_path,
)
from mn_molecular_design_studio.engine.workflow_engine.runtime import reference_root
from mn_molecular_design_studio.job_catalog import (
    job_artifacts,
    job_list_facets,
    path_is_within_accessible_job,
    search_jobs,
)
from mn_molecular_design_studio.engine.workflow_engine.services.local_automation import (
    AutomationError,
    cancel_job,
    inspect_job,
    job_results,
    load_campaign_request,
    resolve_job,
    submit_design_campaign_request,
)
from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import (
    submit_workflow_request,
    validate_workflow_request,
    workflow_catalog,
    workflow_schema,
)


app = FastAPI(
    title="Molecular Design Studio API",
    version="0.1.0",
    description="Local API for discovering, validating, and submitting molecular design workflows.",
)
_bearer = HTTPBearer(auto_error=False)


@app.middleware("http")
async def refresh_workspace_settings(request, call_next):
    # The Streamlit settings page and API are separate processes. Reload the
    # configured results directory for every request so a saved change applies
    # without restarting the API server.
    ensure_studio_workspace()
    return await call_next(request)


class CancelRequest(BaseModel):
    reason: str = "Cancelled through Molecular Design Studio API"


class LoginRequest(BaseModel):
    username: str
    password: str


class PasswordChangeRequest(BaseModel):
    new_password: str
    current_password: str | None = None


def _session_user(credentials: HTTPAuthorizationCredentials | None = Depends(_bearer)) -> dict[str, Any]:
    from mn_molecular_design_studio.cluster import config as cluster_config
    if cluster_config.is_worker():
        user = resolve_session(credentials.credentials if credentials else None)
        if user is None:
            raise HTTPException(status_code=401, detail="Sign in with your coordinator account.")
        return user
    if is_development_mode():
        return development_user()
    if credentials is None or credentials.scheme.casefold() != "bearer":
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="Sign in is required.")
    user = resolve_session(credentials.credentials)
    if user is None:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="Session expired or invalid. Sign in again.")
    return user


def _current_user(user: dict[str, Any] = Depends(_session_user)) -> dict[str, Any]:
    if user.get("must_change_password"):
        raise HTTPException(status_code=403, detail="Change your temporary password before using the API.")
    return user


def _admin_user(user: dict[str, Any] = Depends(_current_user)) -> dict[str, Any]:
    if user.get("role") != "admin":
        raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Admin access is required.")
    return user


def _require_api_job_access(user: dict[str, Any], row: dict[str, Any]) -> None:
    from mn_molecular_design_studio.job_catalog import get_job_access_record

    record = get_job_access_record(str(row.get("task_group") or ""), str(row.get("run_id") or ""))
    try:
        require_job_access(user, record)
    except PermissionError as exc:
        # Return the same public response for an unknown and an inaccessible
        # job so references cannot be used to enumerate other users' runs.
        raise HTTPException(status_code=404, detail="Job not found.") from exc


def _resolve_api_job(reference: str, user: dict[str, Any], *, include_hidden: bool = False) -> dict[str, Any]:
    try:
        row = resolve_job(reference, include_hidden=include_hidden)
    except AutomationError as exc:
        if user.get("role") != "admin":
            if "is ambiguous" in str(exc):
                raise HTTPException(status_code=409, detail="Job reference is ambiguous. Use a full run ID from your accessible jobs.") from exc
            raise HTTPException(status_code=404, detail="Job not found.") from exc
        raise _http_error(exc) from exc
    _require_api_job_access(user, row)
    return row


def _http_error(exc: Exception) -> HTTPException:
    message = str(exc)
    from mn_molecular_design_studio.cluster.transport import ClusterUnavailable
    if isinstance(exc, ClusterUnavailable):
        return HTTPException(status_code=503, detail=message)
    if isinstance(exc, LookupError):
        return HTTPException(status_code=404, detail=message)
    if isinstance(exc, PermissionError):
        return HTTPException(status_code=403, detail=message)
    if isinstance(exc, AutomationError):
        if message.startswith("No job matches"):
            return HTTPException(status_code=404, detail=message)
        if "is ambiguous" in message:
            return HTTPException(status_code=409, detail="Job reference is ambiguous. Use a full run ID from your accessible jobs.")
        return HTTPException(status_code=422, detail=message)
    if isinstance(exc, (OSError, TypeError, ValueError)):
        return HTTPException(status_code=422, detail=message)
    return HTTPException(status_code=500, detail=message)


def _call(function: Callable[..., Any], *args: Any, **kwargs: Any) -> Any:
    try:
        return function(*args, **kwargs)
    except (AutomationError, OSError, RuntimeError, TypeError, ValueError, LookupError) as exc:
        raise _http_error(exc) from exc


def _with_request_file(
    payload: object,
    operation: Callable[[Path], Any],
    *,
    base_dir: str | None,
) -> Any:
    if not isinstance(payload, dict):
        raise HTTPException(status_code=422, detail="Request body must be a JSON object.")
    directory = Path(base_dir or Path.cwd()).expanduser().resolve()
    if not directory.is_dir():
        raise HTTPException(status_code=422, detail=f"base_dir is not a directory: {directory}")
    request_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            suffix=".json",
            prefix=".mn-molecular-design-studio-",
            dir=directory,
            delete=False,
        ) as request_file:
            json.dump(payload, request_file)
            request_path = Path(request_file.name)
        return _call(operation, request_path)
    except HTTPException:
        raise
    except OSError as exc:
        raise _http_error(exc) from exc
    finally:
        if request_path is not None:
            request_path.unlink(missing_ok=True)


def _request_directory(user: dict[str, Any], base_dir: str | None) -> Path:
    if user.get("role") == "admin":
        directory = Path(base_dir or Path.cwd()).expanduser().resolve()
    else:
        if base_dir:
            raise HTTPException(status_code=403, detail="Only an admin can choose a server-side request directory.")
        directory = user_input_root(str(user["id"])).parent / "requests"
        directory.mkdir(parents=True, exist_ok=True, mode=0o700)
        try:
            directory.chmod(0o700)
        except OSError:
            pass
    if not directory.is_dir():
        raise HTTPException(status_code=422, detail=f"Request directory is not a directory: {directory}")
    return directory.resolve()


def _path_is_within(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.expanduser().resolve())
        return True
    except (OSError, ValueError):
        return False


def _is_path_input_key(key: str) -> bool:
    lowered = key.casefold()
    if lowered == "source_job" or lowered in {"name", "type", "sequence", "chains", "chain"}:
        return False
    if lowered.endswith(("_name", "_type", "_sequence", "_chains", "_chain", "_id", "_code", "_goal", "_format", "_mode", "_count", "_strength")):
        return False
    return is_operational_path_key(lowered)


def _validate_user_path_values(value: Any, user: dict[str, Any], request_dir: Path, key: str = "") -> None:
    """Keep non-admin workflow inputs inside their uploads or readable data."""
    if user.get("role") == "admin":
        return
    if isinstance(value, dict):
        if key == "json_config":
            # RFdiffusion3 resolves each input by basename against input_files,
            # whose owned filesystem paths are checked separately. A spec's
            # input is a staged resource name, not a path under request_dir.
            for specification in value.values():
                if isinstance(specification, dict):
                    _validate_user_path_values({k: v for k, v in specification.items() if k != "input"}, user, request_dir)
            return
        for child_key, child_value in value.items():
            _validate_user_path_values(child_value, user, request_dir, str(child_key))
        return
    if isinstance(value, list):
        for child_value in value:
            _validate_user_path_values(child_value, user, request_dir, key)
        return
    if not isinstance(value, str) or not _is_path_input_key(key):
        return
    raw = value.strip()
    if not raw:
        return
    # These references are resolved by the cluster catalog/staging layer. They
    # are not local filesystem paths on the submitting workstation, so local
    # path authorization must not try to resolve a foreign run directory.
    if raw.startswith(("cluster://", "cluster-input:///")):
        return
    if key.casefold()=='model_weights':
        from mn_molecular_design_studio.generator_input import RF_MODELS
        if raw in RF_MODELS:return
    path_like = (
        is_portable_path(raw)
        or Path(raw).expanduser().is_absolute()
        or raw.startswith(".")
        or "/" in raw
        or "\\" in raw
        or Path(raw).suffix.casefold() in {".a3m", ".cif", ".csv", ".fa", ".fasta", ".json", ".jsonl", ".mmcif", ".msa", ".pdb", ".tsv", ".yaml", ".yml", ".zip"}
        or key.casefold().endswith(("_path", "_file", "_dir", "_directory", "_root", "_folder", "_pdb", "_cif", "_mmcif", "_fasta", "_fa", "_json", "_jsonl", "_csv", "_tsv", "_yaml", "_yml", "_a3m", "_msa", "_weights", "_repository", "_table"))
    )
    if not path_like:
        return
    if is_portable_path(raw):
        resolved = resolve_stored_path(raw, must_exist=True)
        if resolved is None:
            raise HTTPException(status_code=422, detail=f"{key} does not resolve to an existing managed file.")
        path = resolved.resolve()
    else:
        path = Path(raw).expanduser()
        if not path.is_absolute():
            path = request_dir / path
        path = path.resolve()
    allowed = (
        _path_is_within(path, user_input_root(str(user["id"])))
        or _path_is_within(path, reference_root())
        or path_is_within_accessible_job(
            path,
            access_owner_id=str(user["id"]),
            accessible_project_ids=accessible_project_ids(user),
        )
    )
    if not allowed:
        raise HTTPException(
            status_code=403,
            detail=f"{key} must refer to a file you uploaded, a shared reference file, or a job you can access.",
        )


def _guard_source_job(payload: dict[str, Any], user: dict[str, Any], *, project_id: str) -> None:
    parameters = payload.get("parameters")
    if isinstance(parameters, list):
        from mn_molecular_design_studio.batching import parameter_items
        for item in _call(parameter_items, parameters):
            _guard_source_job({**payload, "parameters": item}, user, project_id=project_id)
        return
    if payload.get("workflow") in {"sequence.mpnn", "sequence.abmpnn", "analysis.sabr", "sequence.esm_if1", "sequence.caliby", "sequence.boltzgen_inverse", "ligand.md.analysis"} and isinstance(parameters, dict) and parameters.get("source_run_dir"):
        from mn_molecular_design_studio.job_catalog import get_job_access_record
        run_dir = resolve_stored_path(parameters["source_run_dir"], must_exist=True)
        if not run_dir or not run_dir.is_dir():
            raise HTTPException(status_code=422, detail="The selected sequence-design source job is unavailable.")
        row = _resolve_api_job(f"{run_dir.parent.name}/{run_dir.name}", user)
        record = get_job_access_record(str(row.get("task_group") or run_dir.parent.name), run_dir.name)
        if not record or str(record.get("project_id") or "") != project_id:
            raise HTTPException(status_code=422, detail="The sequence-design source job must belong to the selected project.")
    sources = parameters.get("validation_sources") if isinstance(parameters, dict) else None
    if sources is not None:
        if not isinstance(sources, list) or not sources or len(sources) > 128:
            raise HTTPException(status_code=422, detail="Select between 1 and 128 design campaigns.")
        from mn_molecular_design_studio.job_catalog import get_job_access_record
        for source in sources:
            if not isinstance(source, dict):
                raise HTTPException(status_code=422, detail="Each design campaign must specify its run directory and design IDs.")
            run_dir = resolve_stored_path(source.get("source_run_dir"), must_exist=True)
            if not run_dir or not run_dir.is_dir():
                raise HTTPException(status_code=422, detail="A selected design campaign is unavailable.")
            row = _resolve_api_job(f"design/{run_dir.name}", user)
            if Path(row["run_dir"]).resolve() != run_dir.resolve():
                raise HTTPException(status_code=422, detail="Select a managed design campaign.")
            record = get_job_access_record("design", run_dir.name)
            if not record or str(record.get("project_id") or "") != project_id:
                raise HTTPException(status_code=422, detail="Every selected design campaign must belong to the selected project.")
    source_job = payload.get("source_job")
    if source_job is None:
        return
    if not isinstance(source_job, str) or not source_job.strip():
        raise HTTPException(status_code=422, detail="source_job must be a visible job ID or code.")
    row = _resolve_api_job(source_job, user)
    from mn_molecular_design_studio.job_catalog import get_job_access_record

    record = get_job_access_record(str(row.get("task_group") or ""), str(row.get("run_id") or ""))
    source_project_id = str(record.get("project_id") or "") if record else ""
    if source_project_id != project_id:
        raise HTTPException(
            status_code=422,
            detail="The selected results job must belong to the same project as this analysis job.",
        )


def _request_payload_and_project(payload: dict[str, Any], user: dict[str, Any]) -> tuple[dict[str, Any], str]:
    request_payload = dict(payload)
    try:
        project_id = validate_project_id(request_payload.pop("project_id", None))
    except ValueError as exc:
        if user.get("role") == "admin":
            raise HTTPException(status_code=422, detail=str(exc)) from exc
        raise HTTPException(status_code=404, detail="Project not found.") from exc
    if project_id is None:
        raise HTTPException(status_code=422, detail="Select a project before submitting a job.")
    try:
        require_project_access(user, project_id)
    except PermissionError as exc:
        raise HTTPException(status_code=404, detail="Project not found.") from exc
    return request_payload, project_id


@app.get("/api/v1/health")
def health() -> dict[str, str]:
    return {"status": "ok", "app": "mn-molecular-design-studio", "api_version": "v1"}


@app.post("/api/v1/auth/login")
def login(payload: LoginRequest) -> dict[str, Any]:
    from mn_molecular_design_studio.cluster import config as cluster_config
    if cluster_config.is_worker():
        from mn_molecular_design_studio.cluster import transport
        return _call(transport.request, "/api/v1/auth/login", {"username": payload.username, "password": payload.password}, method="POST")
    if is_development_mode():
        user = development_user()
        return {
            "access_token": "development-mode",
            "token_type": "bearer",
            "expires_at": None,
            "user": {key: user[key] for key in ("id", "username", "display_name", "role", "must_change_password")},
        }
    authenticated = authenticate(payload.username, payload.password)
    if authenticated is None:
        raise HTTPException(status_code=401, detail="Username or password is incorrect, or the account is unavailable.")
    user, token, expires_at = authenticated
    return {
        "access_token": token,
        "token_type": "bearer",
        "expires_at": expires_at,
        "user": {key: user[key] for key in ("id", "username", "display_name", "role", "must_change_password")},
    }


@app.get("/api/v1/auth/me")
def who_am_i(user: dict[str, Any] = Depends(_session_user)) -> dict[str, Any]:
    return {**{key: user[key] for key in ("id", "username", "display_name", "role", "must_change_password")},
            "project_ids": accessible_project_ids(user), "session_expires_at": user.get("session_expires_at")}


@app.post("/api/v1/auth/logout")
def logout(credentials: HTTPAuthorizationCredentials | None = Depends(_bearer), user: dict[str, Any] = Depends(_session_user)) -> dict[str, str]:
    from mn_molecular_design_studio.auth import revoke_session

    revoke_session(credentials.credentials if credentials is not None else None)
    return {"status": "signed_out"}


@app.post("/api/v1/auth/password")
def update_password(
    payload: PasswordChangeRequest,
    user: dict[str, Any] = Depends(_session_user),
) -> dict[str, str]:
    try:
        change_password(
            str(user["id"]),
            payload.new_password,
            current_password=payload.current_password,
            forced=bool(user.get("must_change_password")),
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    return {"status": "password_changed", "message": "All sessions were signed out. Sign in again with the new password."}


@app.get("/api/v1/workflows")
def get_workflows(category: str | None = None, user: dict[str, Any] = Depends(_current_user)) -> list[dict[str, Any]]:
    return _call(workflow_catalog, category=category)


@app.get("/api/v1/workflows/{workflow_id}/schema")
def get_workflow_schema(workflow_id: str, user: dict[str, Any] = Depends(_current_user)) -> dict[str, Any]:
    return _call(workflow_schema, workflow_id)


@app.post("/api/v1/workflows/validate")
def validate_workflow(
    payload: dict[str, Any] = Body(...),
    base_dir: str | None = Query(default=None, description="Directory used to resolve relative input paths."),
    user: dict[str, Any] = Depends(_current_user),
    credentials: HTTPAuthorizationCredentials | None = Depends(_bearer),
) -> dict[str, Any]:
    request_payload, project_id = _request_payload_and_project(payload, user)
    request_dir = _request_directory(user, base_dir)
    _validate_user_path_values(request_payload, user, request_dir)
    remote = _cluster_submission(request_payload, project_id, user, request_dir, credentials, validate_only=True)
    if remote is not None:
        return remote
    _guard_source_job(request_payload, user, project_id=project_id)
    return _with_request_file(request_payload, validate_workflow_request, base_dir=str(request_dir))


@app.post("/api/v1/workflows/placement-options")
def workflow_placement_options(
    payload: dict[str, Any] = Body(...),
    base_dir: str | None = Query(default=None, description="Directory used to resolve relative input paths."),
    user: dict[str, Any] = Depends(_current_user),
    credentials: HTTPAuthorizationCredentials | None = Depends(_bearer),
) -> dict[str, Any]:
    request_payload, project_id = _request_payload_and_project(payload, user)
    request_dir = _request_directory(user, base_dir)
    _validate_user_path_values(request_payload, user, request_dir)
    from mn_molecular_design_studio.cluster.service import placement_options_from_api
    token = credentials.credentials if credentials else "development-mode" if is_development_mode() else ""
    return _call(
        placement_options_from_api, dict(request_payload, project_id=project_id), user,
        request_dir=request_dir, token=token,
        authorize_path=lambda path: _validate_user_path_values(str(path), user, request_dir, "input_path"),
    )


@app.post("/api/v1/workflows", status_code=202)
def submit_workflow(
    payload: dict[str, Any] = Body(...),
    base_dir: str | None = Query(default=None, description="Directory used to resolve relative input paths."),
    user: dict[str, Any] = Depends(_current_user),
    credentials: HTTPAuthorizationCredentials | None = Depends(_bearer),
) -> dict[str, Any]:
    request_payload, project_id = _request_payload_and_project(payload, user)
    request_dir = _request_directory(user, base_dir)
    _validate_user_path_values(request_payload, user, request_dir)
    remote = _cluster_submission(request_payload, project_id, user, request_dir, credentials)
    if remote is not None:
        return remote
    _guard_source_job(request_payload, user, project_id=project_id)
    result = _with_request_file(request_payload, submit_workflow_request, base_dir=str(request_dir))
    return mark_studio_job(result, project_id=project_id, owner_id=str(user["id"]))


def _cluster_submission(payload, project_id, user, request_dir, credentials, *, kind="workflow", validate_only=False):
    from mn_molecular_design_studio.cluster.service import submit_from_api
    token = credentials.credentials if credentials else "development-mode" if is_development_mode() else ""
    response = _call(submit_from_api, dict(payload, project_id=project_id), user, request_dir=request_dir,
                     token=token, kind=kind, validate_only=validate_only,
                     authorize_path=lambda path: _validate_user_path_values(str(path), user, request_dir, "input_path"))
    payload.pop("execution", None)
    return response


@app.post("/api/v1/campaigns/validate")
def validate_campaign(
    payload: dict[str, Any] = Body(...),
    base_dir: str | None = Query(default=None, description="Directory used to resolve relative input paths."),
    user: dict[str, Any] = Depends(_current_user),
) -> dict[str, Any]:
    request_payload, _ = _request_payload_and_project(payload, user)
    request_dir = _request_directory(user, base_dir)
    _validate_user_path_values(request_payload, user, request_dir)
    return _with_request_file(request_payload, load_campaign_request, base_dir=str(request_dir))


@app.post("/api/v1/campaigns", status_code=202)
def submit_campaign(
    payload: dict[str, Any] = Body(...),
    base_dir: str | None = Query(default=None, description="Directory used to resolve relative input paths."),
    user: dict[str, Any] = Depends(_current_user),
) -> dict[str, Any]:
    request_payload, project_id = _request_payload_and_project(payload, user)
    request_dir = _request_directory(user, base_dir)
    _validate_user_path_values(request_payload, user, request_dir)
    result = _with_request_file(request_payload, submit_design_campaign_request, base_dir=str(request_dir))
    return mark_studio_job(result, project_id=project_id, owner_id=str(user["id"]))


def _pipeline_payload(payload: dict, user: dict, *, runnable=False):
    import re
    from mn_molecular_design_studio.pipeline_graph import clean_graph, validate_graph
    from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import (
        source_designs, source_generations, source_molecule_selection, source_project_designs,
    )

    request, project_id = _request_payload_and_project(payload, user)
    if set(request) - {"graph", "version", "reuse_completed", "gpu_device"}:
        raise HTTPException(status_code=422, detail="Unsupported pipeline request fields.")
    graph = _call(clean_graph, request.get("graph"))
    _validate_user_path_values(graph, user, _request_directory(user, None))
    for node in graph["nodes"]:
        if node["data"]["kind"] in {"design_source","generation_source","backbone_source","molecule_source","project_design_source"} and node["data"]["params"].get("source_run_dir"):
            params = node["data"]["params"]
            _guard_source_job({"parameters": {"validation_sources": [{"source_run_dir": params["source_run_dir"]}]}}, user, project_id=project_id)
            kind=node['data']['kind']
            loader=(source_molecule_selection if kind == 'molecule_source' else source_project_designs if kind == 'project_design_source' else source_generations if kind in {'generation_source','backbone_source'} else source_designs)
            extra={'stage':'generation.backbone' if kind=='backbone_source' else 'generation.backbone_sequence'} if kind in {'generation_source','backbone_source'} else {}
            _call(loader, params["source_run_dir"], project_id, params.get("selected_candidate_ids"),**extra)
    if runnable:
        plan = _call(validate_graph, graph, require_inputs=True)
        graph = plan["graph"]
        from mn_molecular_design_studio.rfdiffusion3_input import inventory
        for node in graph["nodes"]:
            if node["data"]["kind"] != "target":
                continue
            params = node["data"]["params"]
            target_path = resolve_stored_path(params["target_pdb"], must_exist=True)
            if target_path is None or not target_path.is_file():
                raise HTTPException(status_code=422, detail="The target structure is unavailable.")
            available = _call(inventory, target_path)["chains"]
            chains = [c for c in re.split(r"[,;\s]+", params["target_chains"]) if c]
            if not chains or not set(chains).issubset(available):
                raise HTTPException(status_code=422, detail="Select existing chains from the target structure.")
        from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
        from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import node_request
        nodes = {n["id"]: n for n in graph["nodes"]}
        target_contexts = {}
        for node_id in plan["order"]:
            current = nodes[node_id]
            current_kind = current["data"]["kind"]
            current_params = current["data"]["params"]
            parents = plan["incoming"][node_id]
            target = dict(target_contexts.get(parents[0], {})) if len(parents) == 1 else {}
            if current_kind == "target":
                target = {
                    "target_pdb": current_params.get("target_pdb", ""),
                    "target_chains": [c for c in re.split(r"[,;\s]+", current_params.get("target_chains", "")) if c],
                    "hotspots": current_params.get("hotspots", ""),
                    "coldspots": current_params.get("coldspots", ""),
                }
            elif current_kind == "binder_context":
                if current_params.get("clear_unselected"):
                    target["hotspots"] = current_params.get("hotspots", "")
                    target["coldspots"] = current_params.get("coldspots", "")
                else:
                    target["hotspots"] = current_params.get("hotspots") or target.get("hotspots", "")
                    target["coldspots"] = current_params.get("coldspots") or target.get("coldspots", "")
            target_contexts[node_id] = target
        for node in graph["nodes"]:
            kind = node["data"]["kind"]
            if NODE_CATALOG[kind]["category"] not in {"Design engine", "Design pipeline"} and kind != "surface_analysis":
                continue
            parents=plan["incoming"][node["id"]]
            target = dict(target_contexts.get(parents[0], {})) if parents else {}
            params = node["data"]["params"]
            if not target.get("target_pdb") and params.get("target_pdb"):
                target["target_pdb"] = params["target_pdb"]
                target["target_chains"] = [c for c in re.split(r"[,;\s]+", params.get("target_chains", "")) if c]
            target.setdefault("target_chains", [])
            target.setdefault("hotspots", "")
            child_request = _call(node_request, node["data"]["kind"], node["data"]["params"], {"target": target}, graph["name"])
            _with_request_file(child_request, validate_workflow_request, base_dir=str(_request_directory(user, None)))
    return graph, project_id, request


@app.get("/api/v1/pipelines/catalog")
def pipeline_node_catalog(user: dict[str, Any] = Depends(_current_user)):
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    return NODE_CATALOG


@app.post("/api/v1/pipelines/save")
def save_pipeline(payload: dict[str, Any] = Body(...), user: dict[str, Any] = Depends(_current_user)):
    from mn_molecular_design_studio.pipeline_graph import save_definition
    graph, project_id, request = _pipeline_payload(payload, user)
    if type(request.get("version", False)) is not bool:
        raise HTTPException(status_code=422, detail="version must be true or false.")
    return _call(save_definition, graph, str(user["id"]), project_id, version=request.get("version", False))


@app.get("/api/v1/pipelines")
def list_saved_pipelines(project_id: str, user: dict[str, Any] = Depends(_current_user)):
    from mn_molecular_design_studio.pipeline_graph import list_definitions
    _, project_id = _request_payload_and_project({"project_id": project_id}, user)
    return _call(list_definitions, str(user["id"]), project_id)


@app.post("/api/v1/pipelines/validate")
def validate_pipeline(payload: dict[str, Any] = Body(...), user: dict[str, Any] = Depends(_current_user),
                      credentials: HTTPAuthorizationCredentials | None = Depends(_bearer)):
    from mn_molecular_design_studio.pipeline_graph import validate_graph
    request, project_id = _request_payload_and_project(payload, user)
    request_dir = _request_directory(user, None)
    _validate_user_path_values(request.get("graph"), user, request_dir)
    remote = _cluster_submission(request, project_id, user, request_dir, credentials, kind="pipeline", validate_only=True)
    if remote is not None:
        return remote
    payload = dict(request, project_id=project_id)
    graph, _, _ = _pipeline_payload(payload, user, runnable=True)
    plan = _call(validate_graph, graph, require_inputs=True)
    return {"valid": True, "order": plan["order"], "node_count": len(graph["nodes"]), "connection_count": len(graph["edges"])}


@app.post("/api/v1/pipelines/run", status_code=202)
def run_visual_pipeline(payload: dict[str, Any] = Body(...), user: dict[str, Any] = Depends(_current_user),
                        credentials: HTTPAuthorizationCredentials | None = Depends(_bearer)):
    from mn_molecular_design_studio.pipeline_graph import save_definition
    from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import enqueue_graph_pipeline
    from mn_molecular_design_studio.engine.workflow_engine.services.local_automation import job_reference
    request, selected_project = _request_payload_and_project(payload, user)
    request_dir = _request_directory(user, None)
    _validate_user_path_values(request.get("graph"), user, request_dir)
    remote = _cluster_submission(request, selected_project, user, request_dir, credentials, kind="pipeline")
    if remote is not None:
        return remote
    payload = dict(request, project_id=selected_project)
    graph, project_id, request = _pipeline_payload(payload, user, runnable=True)
    reuse = request.get("reuse_completed", True)
    if type(reuse) is not bool:
        raise HTTPException(status_code=422, detail="reuse_completed must be true or false.")
    device = str(request.get("gpu_device") or "all")
    if device not in {"all", "0", "1"}:
        raise HTTPException(status_code=422, detail="Choose automatic GPU assignment, GPU 0, or GPU 1.")
    version = _call(save_definition, graph, str(user["id"]), project_id, version=True)
    run_dir = _call(enqueue_graph_pipeline, graph, pipeline_name=graph["name"], owner_id=str(user["id"]), project_id=project_id, reuse_completed=reuse, gpu_device=device)
    result = {"job": job_reference(run_dir), "revision": version["revision"]}
    return mark_studio_job(result, project_id=project_id, owner_id=str(user["id"]))


@app.post("/api/v1/pipelines/adopt-completed")
def adopt_completed_pipeline_groups(payload: dict[str, Any] = Body(...), user: dict[str, Any] = Depends(_current_user)):
    """Group already completed runs in Results without scheduling any work."""
    if set(payload) != {"project_id", "pipelines"} or not isinstance(payload.get("pipelines"), list):
        raise HTTPException(status_code=422, detail="Provide project_id and a pipelines list.")
    _, project_id = _request_payload_and_project({"project_id": payload.get("project_id")}, user)
    pipelines = payload["pipelines"]
    if not 1 <= len(pipelines) <= 8:
        raise HTTPException(status_code=422, detail="Adopt between 1 and 8 completed pipelines.")
    from mn_molecular_design_studio.job_catalog import get_job_access_record

    resolved_rows: dict[str, dict[str, Any]] = {}
    normalized: list[dict[str, Any]] = []
    all_refs: set[str] = set()
    for pipeline in pipelines:
        if not isinstance(pipeline, dict) or set(pipeline) != {"name", "steps"}:
            raise HTTPException(status_code=422, detail="Each pipeline needs exactly a name and three steps.")
        name = str(pipeline.get("name") or "").strip()
        steps = pipeline.get("steps")
        if not name or len(name) > 160 or not isinstance(steps, list) or len(steps) != 3:
            raise HTTPException(status_code=422, detail="Each pipeline needs a name and exactly three steps.")
        step_specs = []
        for index, spec in enumerate(steps):
            if not isinstance(spec, dict) or set(spec) != {"task_group", "run_id", "stage_label"}:
                raise HTTPException(status_code=422, detail="Each step needs task_group, run_id, and stage_label.")
            task_group = str(spec.get("task_group") or "").strip()
            run_id = str(spec.get("run_id") or "").strip()
            label = str(spec.get("stage_label") or "").strip()
            if task_group not in {"design", "refolding-validation"} or not run_id or "/" in run_id or ".." in run_id or not label or len(label) > 120:
                raise HTTPException(status_code=422, detail="A pipeline step has an invalid run reference or label.")
            ref = f"{task_group}/{run_id}"
            if ref in all_refs:
                raise HTTPException(status_code=422, detail="A completed job cannot be assigned to more than one adopted pipeline.")
            all_refs.add(ref)
            row = _resolve_api_job(ref, user)
            if str(row.get("status") or "").casefold() not in {"completed", "succeeded"} or row.get("success") is False:
                raise HTTPException(status_code=422, detail=f"{ref} is not a successful completed job.")
            record = get_job_access_record(task_group, run_id)
            if not record or str(record.get("owner_id") or "") != str(user["id"]) or str(record.get("project_id") or "") != project_id:
                raise HTTPException(status_code=422, detail=f"{ref} is not registered to the current user in the selected project.")
            resolved_rows[ref] = row
            step_specs.append({"task_group": task_group, "run_id": run_id, "stage_label": label})

        paths = [Path(resolved_rows[f"{step['task_group']}/{step['run_id']}"]["run_dir"]) for step in step_specs]
        inputs = [json.loads((path / "input.json").read_text(encoding="utf-8")) for path in paths]
        gen_meta = json.loads((paths[0] / "metadata.json").read_text(encoding="utf-8"))
        seq_meta = json.loads((paths[1] / "metadata.json").read_text(encoding="utf-8"))
        refold_meta = json.loads((paths[2] / "metadata.json").read_text(encoding="utf-8"))
        gen_tool = str(gen_meta.get("tool") or "")
        seq_params = inputs[1].get("params") if isinstance(inputs[1].get("params"), dict) else {}
        refold_params = inputs[2].get("params") if isinstance(inputs[2].get("params"), dict) else {}
        if gen_tool not in {"boltzgen_generator", "rfantibody"}:
            raise HTTPException(status_code=422, detail=f"{step_specs[0]['task_group']}/{step_specs[0]['run_id']} is not a BoltzGen or RFantibody generator run.")
        if str(seq_meta.get("tool") or "") != "ligandmpnn" or str(seq_params.get("model_type") or "") != "ab_mpnn" or str(seq_params.get("sabr_mode") or "") != "sabr":
            raise HTTPException(status_code=422, detail="The middle step must be an existing SAbR-annotated AbMPNN run.")
        if str(refold_meta.get("tool") or "") != "optimized_prediction_boltz2" or str(refold_params.get("engine") or "") != "boltz2":
            raise HTTPException(status_code=422, detail="The final step must be an existing Boltz-2 refold.")
        source_path = resolve_stored_path(seq_params.get("source_run_dir"), must_exist=True)
        if not source_path or source_path.resolve() != paths[0].resolve():
            raise HTTPException(status_code=422, detail="The AbMPNN run does not descend from the selected generator run.")
        selected_ids = seq_params.get("selected_candidate_ids")
        if not isinstance(selected_ids, list) or len(selected_ids) != 1:
            raise HTTPException(status_code=422, detail="The adopted AbMPNN run must select exactly one generated candidate.")
        gen_candidates = paths[0] / "artifacts" / "normalized_candidates" / "candidates.jsonl"
        try:
            gen_rows = [json.loads(line) for line in gen_candidates.read_text(encoding="utf-8").splitlines() if line.strip()]
        except (OSError, json.JSONDecodeError):
            gen_rows = []
        if not any(str(candidate.get("candidate_id")) in {str(value) for value in selected_ids} for candidate in gen_rows):
            raise HTTPException(status_code=422, detail="The selected AbMPNN candidate is absent from the generator outputs.")
        ab_candidates = paths[1] / "artifacts" / "normalized_candidates" / "candidates.jsonl"
        try:
            ab_rows = [json.loads(line) for line in ab_candidates.read_text(encoding="utf-8").splitlines() if line.strip()]
        except (OSError, json.JSONDecodeError):
            ab_rows = []
        binder_sequences = {str(candidate.get("binder_sequence") or "") for candidate in ab_rows if candidate.get("binder_sequence")}
        refold_inputs = inputs[2].get("inputs") if isinstance(inputs[2].get("inputs"), dict) else {}
        tasks = refold_inputs.get("tasks") if isinstance(refold_inputs.get("tasks"), list) else []
        protein_chains = tasks[0].get("protein_chains") if tasks and isinstance(tasks[0].get("protein_chains"), list) else []
        binder_sequences_in_refold = {str(chain.get("sequence") or "") for chain in protein_chains if chain.get("role") == "binder"}
        engine_options = refold_params.get("engine_options") if isinstance(refold_params.get("engine_options"), dict) else {}
        if not binder_sequences or not binder_sequences.intersection(binder_sequences_in_refold):
            raise HTTPException(status_code=422, detail="The Boltz-2 binder sequence does not match the AbMPNN output.")
        if engine_options.get("use_msa") is not True:
            raise HTTPException(status_code=422, detail="The adopted Boltz-2 run must have target MSA enabled.")
        if not tasks or tasks[0].get("has_target_structure") is not False or tasks[0].get("restraints"):
            raise HTTPException(status_code=422, detail="The adopted Boltz-2 run must have target structural templates disabled.")
        mapping_path = paths[2] / "artifacts" / "raw" / "optimized_prediction" / "input_mapping.json"
        try:
            mapping = json.loads(mapping_path.read_text(encoding="utf-8"))
            item_details = mapping["items"][0]["input_details"]
        except (OSError, json.JSONDecodeError, KeyError, IndexError, TypeError):
            raise HTTPException(status_code=422, detail="The Boltz-2 run is missing its input mapping needed to verify the template configuration.")
        if mapping["items"][0].get("raw_yaml") is not False or int(item_details.get("template_count") or 0) != 0:
            raise HTTPException(status_code=422, detail="The adopted Boltz-2 run must use sequence-only inputs with zero templates.")
        step_specs[0]["selected_candidate_ids"] = [str(value) for value in selected_ids]
        normalized.append({"name": name, "steps": step_specs})

    from mn_molecular_design_studio.engine.workflow_engine.workflows.pipeline_adoption import adopt_completed_pipelines
    return _call(adopt_completed_pipelines, normalized, owner_id=str(user["id"]), project_id=project_id, resolved_rows=resolved_rows)


@app.post("/api/v1/pipelines/adoptions/{run_id}/undo")
def undo_completed_pipeline_adoption(run_id: str, project_id: str, user: dict[str, Any] = Depends(_current_user)):
    _, selected_project = _request_payload_and_project({"project_id": project_id}, user)
    parent = _resolve_api_job(run_id, user)
    from mn_molecular_design_studio.job_catalog import get_job_access_record
    access = get_job_access_record(str(parent.get("task_group") or "campaign"), str(parent.get("run_id") or run_id))
    if not access or str(access.get("owner_id") or "") != str(user["id"]) or str(access.get("project_id") or "") != selected_project:
        raise HTTPException(status_code=404, detail="Pipeline not found.")
    from mn_molecular_design_studio.engine.workflow_engine.workflows.pipeline_adoption import undo_pipeline_adoption
    return _call(
        undo_pipeline_adoption,
        Path(parent["run_dir"]),
        owner_id=str(user["id"]),
        project_id=selected_project,
    )


@app.get("/api/v1/jobs")
def get_jobs(
    task_group: str | None = None,
    status: str | None = None,
    statuses: list[str] = Query(default=[]),
    tool: str | None = None,
    tools: list[str] = Query(default=[]),
    query: str = "",
    project_id: str | None = None,
    limit: int = Query(default=100, ge=1, le=1000),
    offset: int = Query(default=0, ge=0),
    cursor: str | None = None,
    include_hidden: bool = False,
    user: dict[str, Any] = Depends(_current_user),
) -> dict[str, Any]:
    selected_statuses = [*statuses, *([status] if status else [])]
    selected_tools = [*tools, *([tool] if tool else [])]
    return _call(
        search_jobs,
        task_group=task_group,
        statuses=selected_statuses,
        tools=selected_tools,
        query=query,
        project_id=project_id,
        access_owner_id=str(user["id"]),
        accessible_project_ids=accessible_project_ids(user),
        access_admin=user.get("role") == "admin",
        limit=min(limit, 1000),
        offset=offset,
        cursor=cursor,
        include_hidden=include_hidden,
    )


@app.get("/api/v1/jobs/facets")
def get_job_facets(user: dict[str, Any] = Depends(_current_user)) -> dict[str, list[str]]:
    return _call(
        job_list_facets,
        access_owner_id=str(user["id"]),
        accessible_project_ids=accessible_project_ids(user),
        access_admin=user.get("role") == "admin",
    )


@app.get("/api/v1/jobs/{reference}")
def get_job(reference: str, include_hidden: bool = False, user: dict[str, Any] = Depends(_current_user)) -> dict[str, Any]:
    _resolve_api_job(reference, user, include_hidden=include_hidden)
    return _call(inspect_job, reference, include_hidden=include_hidden)


@app.get("/api/v1/jobs/{reference}/results")
def get_job_results(
    reference: str,
    include_candidates: bool = False,
    include_hidden: bool = False,
    user: dict[str, Any] = Depends(_current_user),
) -> dict[str, Any]:
    _resolve_api_job(reference, user, include_hidden=include_hidden)
    return _call(
        job_results,
        reference,
        include_candidates=include_candidates,
        include_hidden=include_hidden,
    )


@app.get("/api/v1/jobs/{reference}/artifacts")
def get_job_artifacts(
    reference: str,
    include_hidden: bool = False,
    user: dict[str, Any] = Depends(_current_user),
) -> list[dict[str, Any]]:
    row = _resolve_api_job(reference, user, include_hidden=include_hidden)
    return _call(job_artifacts, str(row.get("task_group") or ""), str(row.get("run_id") or ""))


@app.post("/api/v1/jobs/{reference}/cancel")
def stop_job(
    reference: str,
    payload: CancelRequest | None = None,
    user: dict[str, Any] = Depends(_current_user),
) -> dict[str, Any]:
    reason = payload.reason if payload is not None else "Cancelled through Molecular Design Studio API"
    _resolve_api_job(reference, user)
    return _call(cancel_job, reference, reason=reason)


from mn_molecular_design_studio.cluster.routes import install as _install_cluster_routes
_install_cluster_routes(app, _current_user, _admin_user)

from mn_molecular_design_studio.molecular_database_api import install as _install_molecular_database_routes
_install_molecular_database_routes(app, _current_user)
from mn_molecular_design_studio.molecular_database_peer import install_routes as _install_molecular_peer_routes
_install_molecular_peer_routes(app)
