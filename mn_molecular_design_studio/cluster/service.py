from __future__ import annotations

import copy
import json
import tempfile
import uuid
from pathlib import Path
from urllib.parse import urlsplit

from . import bundles, capabilities, config, transport
from .store import Store


def user_record(user):
    from mn_molecular_design_studio.auth import accessible_project_ids
    return dict(user, project_ids=accessible_project_ids(user))


def lookup_job(reference, user, token=""):
    jid = str(reference).removeprefix("cluster:")
    if config.load()["mode"] == "coordinator":
        return Store().job(jid, user_record(user))
    return transport.request("/api/v1/cluster/jobs/" + jid, token=token)


def schema_check(payload, kind, request_dir):
    """Validate request types without pretending remote files exist locally."""
    from mn_molecular_design_studio.engine.workflow_engine.services import workflow_automation as workflows
    if kind == "pipeline":
        from mn_molecular_design_studio.pipeline_graph import clean_graph, validate_graph
        graph = clean_graph(payload.get("graph"))
        validate_graph(graph, require_inputs=False)
        return
    allowed = {"schema_version", "workflow", "parameters", "resources", "source_job", "priority"}
    if set(payload) - allowed or payload.get("schema_version") != 1 or type(payload.get("schema_version")) is not int:
        raise ValueError("Unsupported workflow request fields or schema version.")
    spec = workflows.WORKFLOWS.get(payload.get("workflow"))
    if spec is None:
        raise ValueError("Unknown workflow.")
    with workflows.deferred_path_validation():
        workflows._validate_generic_parameters(spec, payload.get("parameters", {}), Path(request_dir) / "cluster-request.json", source_job=payload.get("source_job"))
    workflows._validate_resources(payload.get("resources"))
    workflows._request_priority(payload)


def submit_from_api(payload, user, *, request_dir, token, authorize_path, kind="workflow", validate_only=False):
    """Return None for ordinary local execution; otherwise stage and federate."""
    settings = config.load()
    execution = payload.get("execution") or {}
    if not isinstance(execution, dict) or set(execution) - {"placement", "idempotency_key"}:
        raise ValueError("execution accepts placement and idempotency_key.")
    placement = execution.get("placement", settings["default_placement"])
    if settings["mode"] == "disabled":
        if placement not in {"local", None} and payload.get("execution"):
            raise ValueError("Configure workstation sharing in Settings before submitting remote jobs.")
        return None
    source = payload.get("source_job")
    has_cluster_source = bool(isinstance(source, str) and source.startswith("cluster:")) or bool(bundles.referenced_nodes(payload))
    if placement == "local" and not has_cluster_source:
        from .agent import Agent
        if not Agent().storage_online():
            raise transport.ClusterUnavailable("This workstation's Results storage is unavailable. Restore its disk or select another workstation.")
        return None
    request_payload = copy.deepcopy(payload)
    request_payload.pop("execution", None)
    project_id = request_payload.pop("project_id")
    if isinstance(source, str) and source.startswith("cluster:"):
        job = lookup_job(source, user, token)
        if job["project_id"] != project_id:
            raise ValueError("The source job must belong to the selected project.")
        if not job.get("local_reference"):
            raise ValueError("The source job has not been accepted by its workstation yet.")
        request_payload["source_job"] = "cluster://" + job["node_id"] + "/" + job["local_reference"]
    elif source is not None:
        from mn_molecular_design_studio.auth import require_job_access
        from mn_molecular_design_studio.job_catalog import get_job_access_record
        from mn_molecular_design_studio.engine.workflow_engine.services.local_automation import resolve_job
        row = resolve_job(source)
        record = get_job_access_record(row["task_group"], row["run_id"])
        require_job_access(user, record)
        if not record or record.get("project_id") != project_id:
            raise ValueError("The source job must belong to the selected project.")
        request_payload["source_job"] = config.artifact_uri(config.identity(), row["task_group"], row["run_id"])
    schema_check(request_payload, kind, request_dir)
    if validate_only:
        return {"valid": True, "validation_scope": "request_schema", "file_validation": "Performed on the owning workstation before a scientific job is created."}
    prepared, bundle_id, references = bundles.pack(request_payload, base_dir=request_dir, node_id=config.identity(), authorize_path=authorize_path)
    owner_nodes = bundles.referenced_nodes(prepared)
    if len(owner_nodes) > 1:
        raise ValueError("Selected inputs belong to multiple workstations. Keep each pipeline branch on its data-owning workstation and compare the resulting scores.")
    needs = capabilities.requirements(request_payload, references=references)
    envelope = {"kind": kind, "payload": prepared, "project_id": project_id, "bundle_id": bundle_id,
                "requirements": needs, "origin_node": config.identity(), "placement": placement,
                "idempotency_key": execution.get("idempotency_key") or str(uuid.uuid4())}
    if settings["mode"] == "coordinator":
        if bundle_id:
            with Store().connect() as db:
                db.execute("CREATE TABLE IF NOT EXISTS blob_grants(digest TEXT,owner_id TEXT,PRIMARY KEY(digest,owner_id))")
                db.execute("INSERT OR IGNORE INTO blob_grants VALUES (?,?)", (bundle_id, user["id"]))
        return accept_submission(envelope, user)
    if bundle_id:
        transport.request("/api/v1/cluster/blobs", token=token, method="POST", binary=bundles.blob_path(bundle_id).read_bytes(), timeout=90)
    return transport.request("/api/v1/cluster/submissions", envelope, token=token, method="POST", timeout=30)


def accept_submission(envelope, user):
    from mn_molecular_design_studio.workspace import validate_project_id
    from mn_molecular_design_studio.auth import require_project_access
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import is_operational_path_key
    store = Store()
    actor = user_record(user)
    request = {k: envelope.get(k) for k in ("kind", "payload", "project_id", "bundle_id")}
    if request["kind"] not in {"workflow", "pipeline"} or not isinstance(request["payload"], dict):
        raise ValueError("Unsupported cluster submission.")
    project_id = validate_project_id(request["project_id"])
    if not project_id:
        raise ValueError("Select a project.")
    require_project_access(user, project_id)
    origin = str(uuid.UUID(envelope["origin_node"]))
    if origin not in {n["id"] for n in store.nodes()}:
        raise ValueError("The originating workstation must first be enrolled.")
    schema_check(request["payload"], request["kind"], config.home())
    nodes = set()

    def validate(value, key="", *, named_input=False):
        if isinstance(value, dict):
            for k, v in value.items():
                validate(v, str(k), named_input=named_input or key == "json_config")
        elif isinstance(value, list):
            for v in value:
                validate(v, key, named_input=named_input)
        elif isinstance(value, str) and is_operational_path_key(key):
            if value.startswith("cluster://"):
                parsed = urlsplit(value)
                parts = Path(parsed.path.lstrip("/")).parts
                if len(parts) < 2 or ".." in parts:
                    raise ValueError("Invalid source artifact reference.")
                node_id = str(uuid.UUID(parsed.netloc))
                local_ref = "/".join(parts[:2])
                with store.connect() as db:
                    row = db.execute("SELECT id FROM jobs WHERE node_id=? AND local_ref=?", (node_id, local_ref)).fetchone()
                if not row:
                    raise ValueError("The source job is not in the shared catalog yet. Refresh its workstation agent.")
                job = store.job(row["id"], actor)
                if job["project_id"] != project_id:
                    raise ValueError("Every source job must belong to the selected project.")
                nodes.add(node_id)
            elif value.startswith(("reference:///", "cluster-input:///")):
                parsed = urlsplit(value)
                parts = Path(parsed.path.lstrip("/")).parts
                if not parts or ".." in parts or parsed.netloc or parsed.query or parsed.fragment:
                    raise ValueError("Invalid staged input or reference path.")
                if parsed.scheme == "cluster-input" and not request["bundle_id"]:
                    raise ValueError("The staged input bundle is missing.")
            elif named_input and key == "input" and Path(value).name == value:
                # RFdiffusion3 resolves these basenames against staged input_files.
                return
            elif key in {"model_weights", "advanced_settings_file"} and not any(c in value for c in "/\\"):
                return
            elif key == "source_job" or value.startswith(("/", ".", "runs:", "app:")) or "/" in value or "\\" in value or Path(value).suffix.lower() in {".pdb", ".cif", ".mmcif", ".a3m", ".fasta", ".fa", ".csv", ".tsv", ".json", ".jsonl", ".yaml", ".yml", ".pt", ".ckpt"}:
                raise ValueError("Remote submissions require staged inputs or machine-qualified source references.")
    validate(request["payload"])
    if len(nodes) > 1:
        raise ValueError("A pipeline must execute on the workstation that owns its source runs.")
    if request["bundle_id"]:
        path = bundles.blob_path(request["bundle_id"])
        if not path.is_file():
            raise ValueError("The staged input bundle is unavailable.")
        with store.connect() as db:
            db.execute("CREATE TABLE IF NOT EXISTS blob_grants(digest TEXT,owner_id TEXT,PRIMARY KEY(digest,owner_id))")
            granted = db.execute("SELECT 1 FROM blob_grants WHERE digest=? AND owner_id=?", (request["bundle_id"], user["id"])).fetchone()
        if not granted:
            raise PermissionError("The input bundle does not belong to this user.")
    requirements = envelope.get("requirements")
    if not isinstance(requirements, dict):
        raise ValueError("Workstation requirements are missing.")
    request["requirements"] = requirements
    return store.submit(request, requirements, actor, origin_node=origin, placement=envelope.get("placement", "auto"),
                        affinity_node=next(iter(nodes), None), idempotency_key=envelope.get("idempotency_key"))


def shared_jobs(user, *, token="", **filters):
    if config.load()["mode"] == "coordinator":
        return Store().search(user_record(user), **filters)
    from urllib.parse import urlencode
    return transport.request("/api/v1/cluster/jobs?" + urlencode({k: v for k, v in filters.items() if v is not None}, doseq=True), token=token)
