"""Preparation helpers and editable starter graphs shared with non-Streamlit clients."""
from __future__ import annotations

from pathlib import Path

from fastapi import APIRouter, Depends, HTTPException, Query
from fastapi.responses import FileResponse
from pydantic import BaseModel, ConfigDict, Field


class Contract(BaseModel):
    model_config = ConfigDict(extra="forbid")


class Preview(Contract):
    parameters: dict = Field(default_factory=dict)
    offset: int = Field(default=0, ge=0)
    limit: int = Field(default=100, ge=1, le=1000)


class PipelineDraft(Contract):
    template: str
    project_id: str
    parameters: dict = Field(default_factory=dict)
    source_job: str | None = None
    selected_candidate_ids: list[str] = Field(default_factory=list, max_length=10000)
    configuration: dict = Field(default_factory=dict)


class StructurePreview(Contract):
    file_id: str


def install(app, current_user, *, validate_paths, resolve_job):
    from mn_molecular_design_studio import api_files, workflow_setup
    from mn_molecular_design_studio.auth import require_project_access
    from mn_molecular_design_studio.workspace import validate_project_id
    router = APIRouter(prefix="/api/v1", tags=["workflow preparation"])

    def call(fn, *args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except (PermissionError, FileNotFoundError) as exc:
            raise HTTPException(404, "Input not found or inaccessible.") from exc
        except (ValueError, TypeError, KeyError, OSError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.get("/workflows/{workflow_id}/setup")
    def setup(workflow_id: str, user=Depends(current_user)):
        return call(workflow_setup.setup_contract, workflow_id)

    @router.get("/workflows/{workflow_id}/templates")
    def templates(workflow_id: str, library: str = "", offset: int = Query(0, ge=0),
                  limit: int = Query(100, ge=1, le=1000), user=Depends(current_user)):
        return call(workflow_setup.templates, workflow_id, library, offset=offset, limit=limit)

    @router.get("/workflows/{workflow_id}/templates/{template_id}/file")
    def template_file(workflow_id: str, template_id: str, library: str = "", user=Depends(current_user)):
        path = call(workflow_setup.template_file, workflow_id, template_id, library)
        return FileResponse(path, filename=path.name, media_type="chemical/x-pdb")

    @router.post("/workflows/{workflow_id}/templates/{template_id}/stage", status_code=201)
    def stage_template(workflow_id: str, template_id: str, library: str = "", user=Depends(current_user)):
        path = call(workflow_setup.template_file, workflow_id, template_id, library)
        staged = call(api_files.store_bytes, user, path.read_bytes(), name=path.name, media_type="chemical/x-pdb")
        parameters = {"template_library": "upload", "template_structure": staged["reference"],
                      "binder_chain": "A", "peptide_chain": "B"} if workflow_id == "design.kejia_idr" else {"scaffold_ids": [template_id]}
        return {"file": staged, "workflow": workflow_id, "parameters": parameters}

    @router.post("/workflows/{workflow_id}/preview")
    def preview(workflow_id: str, payload: Preview, user=Depends(current_user)):
        parameters = call(api_files.resolve_references, payload.parameters, user)
        validate_paths(parameters, user)
        allowed = {p["name"] for p in call(workflow_setup.setup_contract, workflow_id)["schema"]["parameters"]}
        if set(parameters) - allowed:
            raise HTTPException(422, "Unknown workflow parameters.")
        return call(workflow_setup.preview, workflow_id, parameters, limit=payload.limit, offset=payload.offset)

    @router.post("/workflows/{workflow_id}/examples/{example_id}", status_code=201)
    def example(workflow_id: str, example_id: str, user=Depends(current_user)):
        setup = call(workflow_setup.setup_contract, workflow_id)
        values = next((e for e in setup["examples"] if e["id"] == example_id), None)
        if values is None:
            raise HTTPException(404, "Example not found.")
        parameters = {k: v for k, v in values.items() if k not in {"id", "template_available"}}
        files = []
        if workflow_id == "design.kejia_idr":
            from mn_molecular_design_studio.kejia_idr_input import example_template
            path = example_template()
            content = call(path.read_bytes)
            record = call(api_files.store_bytes, user, content, name=path.name, media_type="chemical/x-pdb")
            files.append(record)
            parameters["template_structure"] = record["reference"]
        return {"workflow": workflow_id, "parameters": parameters, "files": files}

    @router.get("/pipelines/templates")
    def pipeline_templates(user=Depends(current_user)):
        from mn_molecular_design_studio.pipeline_graph import STARTER_TEMPLATES, LEGACY_STARTER_TEMPLATES, TEMPLATE_REFERENCES
        from mn_molecular_design_studio.scientific_tasks import catalogue
        return {"items": [{"id": key, "label": label, "reference": TEMPLATE_REFERENCES.get(key)}
                          for label, key in STARTER_TEMPLATES.items()], "tasks": catalogue(),
                "legacy_items": [{"id": key, "label": label} for label, key in LEGACY_STARTER_TEMPLATES.items()]}

    @router.post("/pipelines/drafts")
    def draft(payload: PipelineDraft, user=Depends(current_user)):
        from mn_molecular_design_studio.pipeline_graph import SUPPORTED_TEMPLATE_IDS, new_graph
        from mn_molecular_design_studio.scientific_tasks import TASKS, new_task_graph
        project = call(validate_project_id, payload.project_id)
        if not project:
            raise HTTPException(422, "Choose an existing project.")
        call(require_project_access, user, project)
        if payload.template not in SUPPORTED_TEMPLATE_IDS:
            raise HTTPException(422, "Choose a template from /pipelines/templates.")
        template = payload.template
        source = None
        if payload.source_job:
            if payload.parameters:
                raise HTTPException(422, "Edit the returned continuation graph to change parameters; source identity is fixed by source_job.")
            if template not in {"kejia_idr", "peptide_trap"}:
                raise HTTPException(422, "This starter does not provide a generated-results continuation helper.")
            source = resolve_job(payload.source_job, user)
            if source.get("project_id") != project:
                raise HTTPException(422, "Continue results in their original project.")
            template += "_generated"
        if payload.configuration and template not in TASKS:
            raise HTTPException(422, "Task configuration requires a scientific task ID; legacy graph settings remain explicit node parameters.")
        configuration = call(api_files.resolve_references, payload.configuration, user)
        validate_paths(configuration, user)
        graph = call(new_task_graph, template, configuration) if template in TASKS else call(new_graph, template)
        params = call(api_files.resolve_references, payload.parameters, user)
        validate_paths(params, user)
        if source:
            from mn_molecular_design_studio.generation_io import generated_candidates
            candidates = call(generated_candidates, Path(source["run_dir"]))
            valid = {c["candidate_id"] for c in candidates if c.get("source_tool") == payload.template}
            chosen = payload.selected_candidate_ids or sorted(valid)
            if not chosen or not set(chosen) <= valid:
                raise HTTPException(422, "Select generated candidates from a successful matching workflow.")
            graph["nodes"][0]["data"]["params"].update(source_run_dir=source["run_dir"], selected_candidate_ids=chosen)
        elif payload.selected_candidate_ids:
            raise HTTPException(422, "Candidate selection requires a source job.")
        if params:
            # Configure the first scientific node, as the standalone forms do.
            target = graph["nodes"][0]["data"]["params"]
            if set(params) - set(target):
                raise HTTPException(422, "Use parameters from the starter's first node; edit other nodes on the returned graph.")
            target.update(params)
        return {"project_id": project, "graph": graph, "saved": False}

    @router.post("/structures/preview")
    def structure_preview(payload: StructurePreview, user=Depends(current_user)):
        from mn_molecular_design_studio.rfdiffusion3_input import inventory
        from mn_molecular_design_studio.docking_structure_check import inspect_receptor_pdb
        path = call(api_files.resolve, user, payload.file_id)
        if path.suffix.lower() not in {".pdb", ".ent", ".cif", ".mmcif"}:
            raise HTTPException(422, "Choose a PDB or mmCIF structure.")
        if path.stat().st_size > 50_000_000:
            raise HTTPException(422, "Use a structure smaller than 50 MB for an interactive preview.")
        result = {"inventory": call(inventory, path), "file": call(api_files.describe, user, payload.file_id)}
        if path.suffix.lower() in {".pdb", ".ent"}:
            result["receptor_check"] = call(inspect_receptor_pdb, path.read_text())
        return result

    app.include_router(router)
