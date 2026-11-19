"""Fast, idempotent admission of one source-bound BinderScope cohort job."""
from __future__ import annotations

import hashlib
import json
from uuid import uuid4

from mn_molecular_design_studio.binderscope import validate_assessment_settings
from mn_molecular_design_studio.data_browser import ReadOnlyLibrary


def validate_submission(payload, user, project_id, *, database=None):
    if type(payload.get("schema_version")) is not int or payload["schema_version"] != 1 or payload.get("workflow") != "analysis.binderscope":
        raise ValueError("Choose the registered BinderScope workflow.")
    if set(payload) - {"schema_version", "workflow", "parameters", "resources", "project_id", "priority", "execution", "binderscope_collection"}:
        raise ValueError("Unsupported saved BinderScope request fields.")
    parameters = payload.get("parameters")
    if not isinstance(parameters, dict) or set(parameters) - {"cpu_cores", "use_recorded_settings", "bfactors_are_plddt", "component_format"}:
        raise ValueError("Unsupported saved BinderScope settings.")
    cpu = parameters.get("cpu_cores", 4)
    if type(cpu) is not int or not 1 <= cpu <= 64:
        raise ValueError("BinderScope requires between 1 and 64 CPU cores.")
    recorded = parameters.get("use_recorded_settings", False)
    if type(recorded) is not bool:
        raise ValueError("Recorded settings must be true or false.")
    if not recorded and (parameters.get("bfactors_are_plddt") is not True or parameters.get("component_format") not in {"linear protein", "linear peptide", "cyclic peptide"}):
        raise ValueError("Confirm the saved binder format and confidence encoding, or use recorded per-design settings.")
    priority = payload.get("priority", 0)
    if type(priority) is not int or not 0 <= priority <= 1000:
        raise ValueError("Priority must be between 0 and 1000.")
    resources = payload.get("resources", {})
    if not isinstance(resources, dict) or set(resources) - {"gpu_device"} or resources.get("gpu_device", "none") not in {"none", "cpu", "off"}:
        raise ValueError("Saved BinderScope analysis is a CPU-only job.")
    selection = payload.get("binderscope_collection")
    if not isinstance(selection, dict) or set(selection) - {"schema", "sources"} or selection.get("schema") != "studio.binderscope_selection.v1":
        raise ValueError("BinderScope requires a versioned saved-structure selection.")
    sources = selection.get("sources")
    if not isinstance(sources, list) or not 1 <= len(sources) <= 2000:
        raise ValueError("BinderScope requires between 1 and 2,000 saved structures.")
    identities = []
    for source in sources:
        if not isinstance(source, dict) or source.get("kind") != "database":
            raise ValueError("Saved BinderScope assessments require SQL design selectors.")
        identity = tuple(source.get(field) for field in ("molecule_id", "structure_id", "group_id"))
        if any(not isinstance(value, str) or not value or len(value) > 512 for value in identity):
            raise ValueError("Each selected design requires exact SQL identifiers.")
        binder, target = source.get("binder_chains"), source.get("target_chains")
        if not isinstance(binder, list) or not isinstance(target, list) or len(binder) != 1 or len(target) != 1 or not all(isinstance(value, str) and value for value in binder + target) or binder == target:
            raise ValueError("Each complex requires distinct binder and target chains.")
        if recorded:
            validate_assessment_settings(source.get("assessment"))
        identities.append(identity)
    if len(set(identities)) != len(identities):
        raise ValueError("Each saved structure may appear only once in a cohort.")
    database = database or ReadOnlyLibrary()
    groups = {group["id"]: group for group in database.groups(user, scope="all", collection="all")}
    for molecule, structure, group in identities:
        if group not in groups or str(groups[group].get("project_id") or "") != str(project_id or ""):
            raise PermissionError("Every selected design must be accessible in this analysis project.")
    with database.connection() as connection:
        count = connection.execute(
            "SELECT count(*) FROM json_each(?) selected JOIN molecule_structures s "
            "ON s.molecule_id=json_extract(selected.value,'$[0]') AND s.id=json_extract(selected.value,'$[1]') "
            "AND s.group_id=json_extract(selected.value,'$[2]') JOIN molecule_members m "
            "ON m.molecule_id=s.molecule_id AND m.group_id=s.group_id", (json.dumps(identities),)).fetchone()[0]
    if count != len(identities):
        raise PermissionError("A selected saved design or complex is no longer available.")
    return {"valid": True, "workflow": "analysis.binderscope", "designs": len(sources), "jobs": 1,
            "cpu_cores": cpu, "structure_validation": "Rechecked by the queued job before calculation."}


def queue_submission(payload, user, project_id):
    from mn_molecular_design_studio.auth import user_input_root, is_development_mode
    from mn_molecular_design_studio.collection_policy import collection_scope
    from mn_molecular_design_studio.engine.workflow_engine.core.central_scheduler import submission_priority
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import submission_tags
    from mn_molecular_design_studio.engine.workflow_engine.core.workflow_queue import enqueue_workflow_call
    from mn_molecular_design_studio.engine.workflow_engine.services.local_automation import job_reference
    from mn_molecular_design_studio.engine.workflow_engine.workflows.binderscope import run_saved_binderscope
    from mn_molecular_design_studio.workspace import mark_studio_job

    validate_submission(payload, user, project_id)
    content = json.dumps(payload["binderscope_collection"], sort_keys=True, separators=(",", ":")).encode()
    digest = hashlib.sha256(content).hexdigest()
    directory = user_input_root(user["id"]) / "binderscope_selections"
    directory.mkdir(parents=True, exist_ok=True, mode=0o700)
    path = directory / (digest + ".json")
    if not path.exists():
        temporary = directory / (digest + "." + uuid4().hex + ".tmp")
        with temporary.open("xb") as handle:
            handle.write(content)
        temporary.chmod(0o600)
        temporary.replace(path)
    parameters = payload["parameters"]
    cpu = parameters.get("cpu_cores", 4)
    with submission_tags(studio_owner_id=str(user["id"]), studio_project_id=project_id or "",
                         studio_development_mode=is_development_mode()), collection_scope("automatic"), submission_priority(payload.get("priority", 0)):
        run = enqueue_workflow_call(run_saved_binderscope, task_group="analysis", tool="binderscope", job_type="analysis",
            kwargs={"selection_file": path, "selection_sha256": digest, "cpu_cores": cpu,
                    "use_recorded_settings": parameters.get("use_recorded_settings", False),
                    "bfactors_are_plddt": parameters.get("bfactors_are_plddt", False),
                    "component_format": parameters.get("component_format", "unspecified")},
            cpu_cores=cpu, gpu_device="none")
    return mark_studio_job({"workflow": "analysis.binderscope", "job": job_reference(run)},
        project_id=project_id or None, owner_id=str(user["id"]), allow_private_source=not project_id)
