"""Shared data operations for Studio, student clients and LLM adapters."""
from __future__ import annotations

import base64
import hashlib
import shutil
import uuid
from pathlib import Path
from typing import Literal

from fastapi import APIRouter, Body, Depends, HTTPException, Query, Request
from fastapi.responses import FileResponse
from pydantic import BaseModel, ConfigDict, Field


class Contract(BaseModel):
    model_config = ConfigDict(extra="forbid")


class LibraryAccess(Contract):
    project_ids: list[str] = Field(default_factory=list, max_length=1000)


class LibraryRegistration(LibraryAccess):
    paths: list[str] = Field(min_length=1)
    name: str = Field(min_length=1, max_length=250)
    id_column: str = "id"
    smiles_column: str = "smiles"
    name_column: str = ""


class PreparationRegistration(Contract):
    manifest_path: str


class ScreeningQuery(Contract):
    offset: int = Field(default=0, ge=0, strict=True)
    limit: int = Field(default=100, ge=1, le=1000, strict=True)
    query: str = Field(default="", max_length=500)
    score_min: float | None = Field(default=None, allow_inf_nan=False)
    score_max: float | None = Field(default=None, allow_inf_nan=False)
    compound_ids: list[str] = Field(default_factory=list, max_length=1000)


class FinalImport(Contract):
    limit: int = Field(default=1000, ge=0, le=1000, strict=True)


class ExportRequest(Contract):
    molecule_ids: list[str] = Field(min_length=1, max_length=100_000)
    format: Literal["CSV", "FASTA", "SMILES", "SDF", "Structures ZIP"] = "CSV"


class BoundInputs(Contract):
    structure_id: str | None = None


class GroupMembers(Contract):
    molecule_ids: list[str] = Field(min_length=1, max_length=100_000)


class GroupProject(Contract):
    project_id: str


class ChainRole(Contract):
    role: str


class CropRequest(Contract):
    chains: list[str] = Field(min_length=1)
    structure_id: str | None = None
    name: str | None = None
    project_id: str = ""
    group_id: str | None = None


class ImportAnalysis(Contract):
    filename: str = "molecules.csv"
    content_base64: str = ""
    file_id: str | None = None
    smiles_column: str
    id_column: str = ""
    name_column: str = ""
    sheet_name: str = ""


class DockingPreparation(Contract):
    project_id: str
    workflow: str
    preparation_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    state: dict = Field(default_factory=dict)


def install(app, current_user, admin_user, *, validate_paths):
    from mn_molecular_design_studio import api_files
    from mn_molecular_design_studio.auth import accessible_project_ids, require_project_access
    from mn_molecular_design_studio.external_libraries import LibraryCatalog
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    from mn_molecular_design_studio.workspace import validate_project_id

    router = APIRouter(prefix="/api/v1", tags=["shared data"])

    def call(function, *args, **kwargs):
        try:
            return function(*args, **kwargs)
        except (PermissionError, FileNotFoundError, KeyError) as exc:
            raise HTTPException(404, "Record not found or inaccessible.") from exc
        except (TypeError, ValueError, OSError) as exc:
            raise HTTPException(422, str(exc)) from exc

    def projects(values, user):
        result = []
        for value in dict.fromkeys(values):
            project = call(validate_project_id, value)
            if not project:
                raise HTTPException(422, "Choose an existing project.")
            call(require_project_access, user, project)
            result.append(project)
        return result

    def readable_library(library_id, user):
        try:
            record = LibraryCatalog().get(library_id)
        except ValueError as exc:
            raise HTTPException(404, "Library not found or inaccessible.") from exc
        if (user.get("role") != "admin" and record.get("owner_id") != user["id"]
                and not set(record.get("project_ids") or []).intersection(accessible_project_ids(user))):
            raise HTTPException(404, "Library not found or inaccessible.")
        return record

    def editable_library(record):
        from mn_molecular_design_studio.cluster.config import identity
        if record.get("_replica_node") not in (None, identity()):
            raise HTTPException(409, "Update library metadata on its owning workstation; its sharing grants will then synchronize to peers.")

    def run(task_group, run_id, user):
        from mn_molecular_design_studio.runtime import runs_root
        from mn_molecular_design_studio.molecular_database_results import _authorized_run
        if any(not part or Path(part).name != part or part in {".", ".."} for part in (task_group, run_id)):
            raise HTTPException(404, "Studio run not found.")
        return call(_authorized_run, Path(runs_root()) / task_group / run_id, user)[0]

    @router.post("/files", status_code=201, openapi_extra={"requestBody": {"required": True, "content": {
        "application/octet-stream": {"schema": {"type": "string", "format": "binary"}},
        "text/plain": {"schema": {"type": "string"}}}}})
    async def upload(request: Request, filename: str = Query(default="input.dat", max_length=200), user=Depends(current_user)):
        """Stream a private input to disk. This does not import molecular records."""
        file_id = uuid.uuid4().hex
        folder = api_files.directory(user, file_id)
        folder.mkdir(parents=True, mode=0o700)
        digest, size = hashlib.sha256(), 0
        try:
            path = folder / api_files.storage_name(filename)
            with path.open("xb") as stream:
                path.chmod(0o600)
                async for chunk in request.stream():
                    size += len(chunk)
                    digest.update(chunk)
                    stream.write(chunk)
            return api_files.finish(user, file_id, name=filename, size=size, sha256=digest.hexdigest(),
                                    media_type=request.headers.get("content-type"))
        except BaseException:
            shutil.rmtree(folder)
            raise

    @router.get("/files/{file_id}/metadata")
    def file_metadata(file_id: str, user=Depends(current_user)):
        return call(api_files.describe, user, file_id)

    @router.get("/files/{file_id}")
    def download(file_id: str, user=Depends(current_user)):
        record = call(api_files.describe, user, file_id)
        return FileResponse(call(api_files.resolve, user, file_id), filename=record["name"], media_type=record["media_type"])

    @router.delete("/files/{file_id}")
    def delete_file(file_id: str, user=Depends(current_user)):
        call(api_files.describe, user, file_id)
        call(shutil.rmtree, api_files.directory(user, file_id))
        return {"status": "deleted", "id": file_id}

    @router.get("/libraries")
    def libraries(project_id: str | None = None, user=Depends(current_user)):
        if project_id:
            projects([project_id], user)
        allowed = set(accessible_project_ids(user))
        items = [r for r in LibraryCatalog().list() if
                 (user.get("role") == "admin" or r.get("owner_id") == user["id"] or allowed.intersection(r.get("project_ids") or []))
                 and (not project_id or project_id in (r.get("project_ids") or []))]
        return {"items": items, "total": len(items)}

    @router.post("/libraries", status_code=201)
    def register_library(payload: LibraryRegistration, user=Depends(admin_user)):
        project_ids = projects(payload.project_ids, user)
        catalog = LibraryCatalog()
        record = call(catalog.register, payload.paths, name=payload.name, id_column=payload.id_column,
                      smiles_column=payload.smiles_column, name_column=payload.name_column)
        editable_library(record)
        # Re-registering an existing data source preserves its sharing grants.
        record.setdefault("owner_id", user["id"])
        record["project_ids"] = sorted(set(record.get("project_ids") or []).union(project_ids))
        return catalog.database.save_external_library_record(record)

    @router.get("/libraries/{library_id}")
    def library(library_id: str, user=Depends(current_user)):
        return readable_library(library_id, user)

    @router.patch("/libraries/{library_id}/access")
    def library_access(library_id: str, payload: LibraryAccess, user=Depends(admin_user)):
        record = readable_library(library_id, user)
        editable_library(record)
        record["project_ids"] = projects(payload.project_ids, user)
        return LibraryCatalog().database.save_external_library_record(record)

    @router.get("/libraries/{library_id}/preparations")
    def preparations(library_id: str, kind: Literal["all", "fingerprints", "centroids"] = "all", user=Depends(current_user)):
        from mn_molecular_design_studio.external_libraries import FEATURE_SCHEMA, CENTROID_SCHEMA
        record = readable_library(library_id, user)
        expected = {"fingerprints": FEATURE_SCHEMA, "centroids": CENTROID_SCHEMA}.get(kind)
        items = [p for p in record.get("preparations", []) if expected is None or p.get("schema_version") == expected]
        return {"library_id": library_id, "items": items, "total": len(items)}

    @router.post("/libraries/{library_id}/preparations", status_code=201)
    def attach_preparation(library_id: str, payload: PreparationRegistration, user=Depends(admin_user)):
        record = readable_library(library_id, user)
        editable_library(record)
        return call(LibraryCatalog().attach_preparation, library_id, payload.manifest_path)

    @router.get("/submissions/status")
    def submission_status(kind: Literal["workflow", "pipeline", "campaign"], key: str = Query(min_length=1, max_length=200), user=Depends(current_user)):
        from mn_molecular_design_studio.api_idempotency import submission_status as inspect_submission
        return call(inspect_submission, user, kind, key)

    @router.get("/screening/{task_group}/{run_id}")
    def screening_summary(task_group: str, run_id: str, user=Depends(current_user)):
        from mn_molecular_design_studio.bulk_screening_results import external_results
        result = call(external_results, run(task_group, run_id, user), user)
        if result is None:
            raise HTTPException(404, "Screening results are not available yet.")
        return result

    @router.post("/screening/{task_group}/{run_id}/search")
    def screening_search(task_group: str, run_id: str, payload: ScreeningQuery, user=Depends(current_user)):
        from mn_molecular_design_studio.screening_queries import query_results
        return call(query_results, run(task_group, run_id, user), user, **payload.model_dump())

    @router.get("/screening/{task_group}/{run_id}/files/{kind}")
    def screening_file(task_group: str, run_id: str, kind: str, user=Depends(current_user)):
        from mn_molecular_design_studio.screening_queries import result_file
        path = call(result_file, run(task_group, run_id, user), user, kind)
        return FileResponse(path, filename=path.name)

    @router.post("/screening/{task_group}/{run_id}/final-selection")
    def screening_final(task_group: str, run_id: str, payload: FinalImport, user=Depends(current_user)):
        from mn_molecular_design_studio.bulk_screening_results import import_final_results
        return call(import_final_results, run(task_group, run_id, user), user, payload.limit)

    @router.post("/molecules/exports", status_code=201)
    def export(payload: ExportRequest, user=Depends(current_user)):
        from mn_molecular_design_studio.molecule_exports import export_molecules
        content, extension, media_type = call(export_molecules, MolecularDatabase(), user, payload.molecule_ids, payload.format)
        return call(api_files.store_bytes, user, content, name="studio-selection." + extension, media_type=media_type)

    @router.post("/molecules/{molecule_id}/bound-inputs", status_code=201)
    def bound_inputs(molecule_id: str, payload: BoundInputs, user=Depends(current_user)):
        from mn_molecular_design_studio.bound_docking import prepare_bound_inputs
        prepared = call(prepare_bound_inputs, MolecularDatabase(), user, molecule_id, structure_id=payload.structure_id)
        files = {}
        for key, media_type in (("complex_pdb", "chemical/x-pdb"), ("ligand_sdf", "chemical/x-mdl-sdfile")):
            path = Path(prepared[key])
            files[key] = call(api_files.store_bytes, user, path.read_bytes(), name=path.name, media_type=media_type)
        return {**prepared, "files": files, "workflow_inputs": {
            "complex_pdb": files["complex_pdb"]["reference"], "ligand_sdf": files["ligand_sdf"]["reference"],
            "ligand_key": prepared["ligand_key"]}}

    @router.post("/molecules/{molecule_id}/crop", status_code=201)
    def crop(molecule_id: str, payload: CropRequest, user=Depends(current_user)):
        from mn_molecular_design_studio.molecular_preparation import cropped_record
        database = MolecularDatabase()
        records = call(database.list_molecules, user, molecule_ids=[molecule_id])["items"]
        if not records:
            raise HTTPException(404, "Molecule not found or inaccessible.")
        record = records[0]
        structures = record.get("representations") or record["structures"]
        structure = next((s for s in structures if payload.structure_id is None or s["id"] == payload.structure_id), None)
        if payload.structure_id and structure is None:
            raise HTTPException(404, "Structure not found or inaccessible.")
        files = [{"path": str(call(database.structure_file, user, structure["id"]))}] if structure else []
        cropped = call(cropped_record, {**record, "files": files}, list(dict.fromkeys(payload.chains)))
        if payload.name is not None:
            cropped["name"] = payload.name
        return call(database.import_records, user, [cropped], name=payload.name or "Cropped molecule",
                    project_id=payload.project_id, group_id=payload.group_id,
                    origin={"type": "crop", "source_molecule_id": molecule_id, "structure_id": payload.structure_id})

    @router.post("/molecules/structures/{structure_id}/portable")
    def portable_structure(structure_id: str, user=Depends(current_user)):
        return call(MolecularDatabase().make_structure_portable, user, structure_id)

    @router.post("/molecules/groups/{group_id}/members/remove")
    def remove_members(group_id: str, payload: GroupMembers, user=Depends(current_user)):
        call(MolecularDatabase().remove_from_group, user, group_id, payload.molecule_ids)
        return {"removed_count": len(set(payload.molecule_ids))}

    @router.patch("/molecules/groups/{group_id}/project")
    def group_project(group_id: str, payload: GroupProject, user=Depends(current_user)):
        return call(MolecularDatabase().set_group_project, user, group_id, payload.project_id)

    @router.patch("/molecules/groups/{group_id}/molecules/{molecule_id}/chains/{chain}/role")
    def chain_role(group_id: str, molecule_id: str, chain: str, payload: ChainRole, user=Depends(current_user)):
        call(MolecularDatabase().set_chain_role, user, group_id, molecule_id, chain, payload.role)
        return {"status": "updated"}

    @router.post("/molecules/imports/analyze")
    def analyze(payload: ImportAnalysis, user=Depends(current_user)):
        from mn_molecular_design_studio import molecule_imports
        if bool(payload.file_id) == bool(payload.content_base64):
            raise HTTPException(422, "Supply either file_id or content_base64.")
        if payload.file_id:
            record = call(api_files.describe, user, payload.file_id)
            if record["size_bytes"] > molecule_imports.MAX_BYTES:
                raise HTTPException(422, "Use external library preparation for tables over 100 MB.")
            content = call(api_files.resolve, user, payload.file_id).read_bytes()
            filename = record["name"]
        else:
            if len(payload.content_base64) > molecule_imports.MAX_BYTES * 4 // 3 + 8:
                raise HTTPException(422, "Use external library preparation for tables over 100 MB.")
            content = call(base64.b64decode, payload.content_base64, validate=True)
            filename = payload.filename
        # No molecular records are written during analysis.
        result = call(molecule_imports.analyze_small_molecule_table, content, filename,
                      smiles_column=payload.smiles_column, name_column=payload.name_column,
                      id_column=payload.id_column, sheet_name=payload.sheet_name)
        result.pop("records", None)
        return result

    @router.get("/molecules/groups/{group_id}/imports")
    def import_batches(group_id: str, user=Depends(current_user)):
        return {"items": call(MolecularDatabase().import_batches, user, group_id)}

    @router.get("/molecules/imports/{batch_id}/rows")
    def import_rows(batch_id: str, limit: int = Query(100, ge=1, le=1000), offset: int = Query(0, ge=0), user=Depends(current_user)):
        rows = call(MolecularDatabase().import_batch_rows, user, batch_id)
        return {"items": rows[offset:offset + limit], "total": len(rows), "offset": offset, "limit": limit}

    @router.get("/molecules/imports/{batch_id}/source")
    def import_source(batch_id: str, user=Depends(current_user)):
        path = call(MolecularDatabase().import_source_file, user, batch_id)
        return FileResponse(path, filename=path.name)

    @router.get("/docking/preparations")
    def docking_preparations(project_id: str, workflow: str, user=Depends(current_user)):
        from mn_molecular_design_studio.docking_preparations import list_saved
        projects([project_id], user)
        return {"items": call(list_saved, user["id"], project_id, workflow)}

    @router.post("/docking/preparations")
    def save_docking_preparation(payload: DockingPreparation, user=Depends(current_user)):
        from mn_molecular_design_studio.docking_preparations import save
        projects([payload.project_id], user)
        state = call(api_files.resolve_references, payload.state, user)
        validate_paths(state, user)
        return call(save, user["id"], payload.project_id, payload.workflow, payload.preparation_id, state)

    app.include_router(router)
