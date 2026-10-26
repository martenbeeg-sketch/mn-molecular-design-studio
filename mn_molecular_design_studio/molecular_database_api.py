"""The same molecular library and import contracts for API clients and UI."""
from __future__ import annotations
import base64
from pathlib import Path
from typing import Any

from fastapi import APIRouter, Body, Depends, HTTPException, Query
from fastapi.responses import FileResponse
from pydantic import BaseModel, Field

from mn_molecular_design_studio import molecule_imports as imports
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio.molecular_database_ingest import sync_workspace


class GroupRequest(BaseModel):
    model_config = {"extra": "forbid"}
    name: str
    project_id: str = ""
    molecule_ids: list[str] = Field(default_factory=list)


class ImportRequest(BaseModel):
    model_config = {"extra": "forbid"}
    name: str = "Imported molecules"
    project_id: str = ""
    group_id: str | None = None
    method: str = "smiles"
    text: str = ""
    filename: str = "molecules.csv"
    content_base64: str = ""
    file_id: str | None = None
    sequence_type: str = "protein"
    chain_columns: dict[str, dict] = Field(default_factory=dict)
    name_column: str = ""
    id_column: str = ""
    sheet_name: str = ""
    review_small_molecules: bool = False
    include_ligands: bool = True
    records: list[dict] = Field(default_factory=list)
    identifier_source: str = "PDB"


class SelectionRequest(BaseModel):
    model_config = {"extra": "forbid"}
    molecule_ids: list[str]
    name: str = "Database selection"
    project_id: str = ""


class RunSelectionRequest(BaseModel):
    model_config = {"extra": "forbid"}
    molecule_ids: list[str]
    group_id: str | None = None
    name: str = "Saved results"
    project_id: str | None = None


class ComparisonRequest(SelectionRequest):
    evaluation_filters: dict = Field(default_factory=dict)
    settings: dict = Field(default_factory=dict)


class ComparisonPreviewRequest(BaseModel):
    model_config = {"extra": "forbid"}
    selection: dict = Field(default_factory=dict)
    settings: dict = Field(default_factory=dict)


def install(app, current_user):
    router = APIRouter(prefix="/api/v1/molecules", tags=["molecular database"])

    def call(function, *args, **kwargs):
        try:
            return function(*args, **kwargs)
        except (PermissionError, FileNotFoundError) as exc:
            raise HTTPException(404, "Molecular record not found or inaccessible.") from exc
        except (ValueError, TypeError, OSError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.get("/groups")
    def groups(project_id: str = "", scope: str = "organization", collection: str = "all", user=Depends(current_user)):
        call(sync_workspace, user)
        return {"items": call(MolecularDatabase().groups, user, project_id=project_id, scope=scope, collection=collection)}

    @router.get("/output-vocabulary")
    def output_vocabulary(user=Depends(current_user)):
        from mn_molecular_design_studio.output_vocabulary import vocabulary
        return vocabulary()

    @router.get("/metrics")
    def metrics(group_id: str = "", project_id: str = "", scope: str = "organization", collection: str = "all", user=Depends(current_user)):
        call(sync_workspace, user)
        return {"items": call(MolecularDatabase().metric_fields, user, group_id=group_id,
                              project_id=project_id, scope=scope, collection=collection)}

    @router.post("/groups", status_code=201)
    def create_group(payload: GroupRequest, user=Depends(current_user)):
        database = MolecularDatabase()
        if payload.molecule_ids:
            call(database.molecules_by_ids, user, payload.molecule_ids)
        identifier = call(database.create_group, user, payload.name, project_id=payload.project_id)
        if payload.molecule_ids:
            call(database.add_to_group, user, identifier, payload.molecule_ids)
        return call(database.group, identifier, user)

    @router.get("/groups/{group_id}")
    def group(group_id: str, user=Depends(current_user)):
        return call(MolecularDatabase().group, group_id, user)

    @router.post("/groups/{group_id}/members")
    def add_members(group_id: str, payload: SelectionRequest, user=Depends(current_user)):
        return {"added_count": call(MolecularDatabase().add_to_group, user, group_id, payload.molecule_ids)}

    @router.post("/imports", status_code=201)
    def import_molecules(payload: ImportRequest, user=Depends(current_user)):
        def operation():
            from mn_molecular_design_studio import api_files
            if len(payload.content_base64) > imports.MAX_BYTES * 4 // 3 + 8:
                raise ValueError("The import exceeds 100 MB.")
            if payload.file_id and payload.content_base64:
                raise ValueError("Supply either file_id or content_base64.")
            data = base64.b64decode(payload.content_base64, validate=True) if payload.content_base64 else b""
            filename = payload.filename
            if payload.file_id:
                descriptor = api_files.describe(user, payload.file_id)
                if descriptor["size_bytes"] > imports.MAX_BYTES:
                    raise ValueError("Use external library preparation for files over 100 MB.")
                data = api_files.resolve(user, payload.file_id).read_bytes()
                filename = descriptor["name"]
            audit = None
            if payload.method in {"fasta", "sequences"}:
                records = imports.parse_sequences(payload.text or data.decode("utf-8-sig"), kind=payload.sequence_type)
            elif payload.method == "smiles":
                records = imports.parse_smiles(payload.text or data.decode("utf-8-sig"))
            elif payload.method == "table":
                if payload.review_small_molecules:
                    mapped = payload.chain_columns
                    if len(mapped) != 1 or next(iter(mapped.values())).get("type") != "small_molecule":
                        raise ValueError("Reviewed tables require one mapped small-molecule SMILES column.")
                    audit = imports.analyze_small_molecule_table(data, filename, smiles_column=next(iter(mapped)),
                        name_column=payload.name_column, id_column=payload.id_column, sheet_name=payload.sheet_name)
                    records = audit.pop("records")
                else:
                    records = imports.parse_table(data, filename, chain_columns=payload.chain_columns,
                        name_column=payload.name_column, id_column=payload.id_column, sheet_name=payload.sheet_name)
            elif payload.method in {"sdf", "structure"}:
                records = imports.parse_structure(data, filename, include_ligands=payload.include_ligands)
            elif payload.method == "zip":
                records = imports.parse_archive(data)
            elif payload.method == "identifier":
                filename, content = imports.fetch_identifier(payload.text, payload.identifier_source)
                records = imports.parse_sequences(content.decode()) if payload.identifier_source == "UniProt" else imports.parse_structure(content, filename)
            elif payload.method == "records":
                records = [{"name": row.get("name"), "components": row["components"], "metadata": row.get("metadata") or {}} for row in payload.records]
            else:
                raise ValueError("Supported methods: table, fasta, sequences, smiles, sdf, structure, zip, identifier, records.")
            return MolecularDatabase().import_records(user, records, name=payload.name, project_id=payload.project_id,
                group_id=payload.group_id, import_audit=audit, source_content=data if audit else None)
        return call(operation)

    @router.post("/selections")
    def freeze_compounds(payload: SelectionRequest, user=Depends(current_user)):
        return call(MolecularDatabase().compound_snapshot, user, payload.molecule_ids, name=payload.name, project_id=payload.project_id)

    @router.get("/evaluations/options")
    def evaluation_options(group_id: str = "", project_id: str = "", collection: str = "all", user=Depends(current_user)):
        return call(MolecularDatabase().evaluation_options, user, group_id=group_id, project_id=project_id, collection=collection)

    @router.post("/evaluations/search")
    def evaluation_search(payload: dict = Body(...), user=Depends(current_user)):
        allowed = {"group_id", "project_id", "scope", "collection", "query", "molecule_type", "molecule_ids", "evaluation_filters"}
        if set(payload) - allowed:
            raise HTTPException(422, "Unsupported evaluation search fields.")
        rows = call(MolecularDatabase().evaluation_rows, user, **payload)
        return {"items": rows, "total": len(rows)}

    @router.post("/comparisons", status_code=201)
    def save_comparison(payload: ComparisonRequest, user=Depends(current_user)):
        return call(MolecularDatabase().save_comparison, user, payload.name, payload.molecule_ids,
                    project_id=payload.project_id, evaluation_filters=payload.evaluation_filters, settings=payload.settings)

    @router.post("/comparisons/preview")
    def compare_evaluations(payload: ComparisonPreviewRequest, user=Depends(current_user)):
        from mn_molecular_design_studio.comparisons import compare, compatibility_warnings, prepare_rows
        selection = payload.selection
        if set(selection) - {"group_id", "project_id", "molecule_ids", "evaluation_filters"}:
            raise HTTPException(422, "Unsupported comparison selection.")
        settings = payload.settings
        if set(settings) - {"fields", "directions", "weights", "thresholds", "method", "combine_same_subject", "evaluator_runs"}:
            raise HTTPException(422, "Unsupported comparison settings.")
        rows = call(MolecularDatabase().evaluation_rows, user, **selection)
        frame = call(compare, rows, **settings)
        prepared = call(prepare_rows, rows, combine_same_subject=settings.get("combine_same_subject", False), evaluator_runs=settings.get("evaluator_runs"))
        warnings = compatibility_warnings(prepared, settings.get("fields") or [])
        if warnings and settings.get("fields"):
            raise HTTPException(422, "; ".join(warnings))
        import json
        return {"items": json.loads(frame.to_json(orient="records")), "total": len(frame), "warnings": warnings}

    @router.post("/selections/proteins")
    def freeze_proteins(payload: SelectionRequest, user=Depends(current_user)):
        from mn_molecular_design_studio.database_selections import freeze_protein_selection
        if not payload.project_id:
            raise HTTPException(422, "A pipeline selection must belong to a project.")
        return {"source_run_dir": str(call(freeze_protein_selection, MolecularDatabase(), user, payload.molecule_ids, payload.project_id, name=payload.name))}

    @router.patch("/groups/{group_id}/schema")
    def update_schema(group_id: str, payload: dict = Body(...), user=Depends(current_user)):
        call(MolecularDatabase().set_schema, user, group_id, payload)
        return {"status": "updated"}

    @router.post("/refresh")
    def refresh(user=Depends(current_user)):
        return {"indexing_errors": call(sync_workspace, user, force=True)}

    @router.get("/jobs/{task_group}/{run_id}")
    def run_molecules(task_group: str, run_id: str, user=Depends(current_user)):
        from mn_molecular_design_studio.runtime import runs_root
        from mn_molecular_design_studio.molecular_database_results import list_run_molecules
        if Path(task_group).name != task_group or Path(run_id).name != run_id or task_group in {".", ".."} or run_id in {".", ".."}:
            raise HTTPException(404, "Studio run not found.")
        items = call(list_run_molecules, Path(runs_root()) / task_group / run_id, user)
        return {"items": items, "total": len(items)}

    @router.post("/jobs/{task_group}/{run_id}/selections")
    def save_run_molecules(task_group: str, run_id: str, payload: RunSelectionRequest, user=Depends(current_user)):
        from mn_molecular_design_studio.runtime import runs_root
        from mn_molecular_design_studio.molecular_database_results import save_run_selection
        if Path(task_group).name != task_group or Path(run_id).name != run_id or task_group in {".", ".."} or run_id in {".", ".."}:
            raise HTTPException(404, "Studio run not found.")
        return call(save_run_selection, Path(runs_root()) / task_group / run_id, user, payload.molecule_ids,
                    group_id=payload.group_id, name=payload.name, project_id=payload.project_id)

    @router.get("/structures/{structure_id}/file")
    def structure_file(structure_id: str, user=Depends(current_user)):
        path = call(MolecularDatabase().structure_file, user, structure_id)
        return FileResponse(path, filename=path.name)

    @router.post("/structures/{structure_id}/materialize")
    def materialize(structure_id: str, user=Depends(current_user)):
        return {"path": str(call(MolecularDatabase().structure_file, user, structure_id, materialize=True))}

    @router.get("")
    def molecules(group_id: str = "", project_id: str = "", collection: str = "all", query: str = "", molecule_type: str = "",
                  sort: str = "added", descending: bool = True, limit: int = Query(100, ge=1, le=1000), offset: int = Query(0, ge=0), user=Depends(current_user)):
        call(sync_workspace, user)
        return call(MolecularDatabase().list_molecules, user, group_id=group_id, project_id=project_id, collection=collection, query=query,
                    molecule_type=molecule_type, sort=sort, descending=descending, limit=limit, offset=offset)

    @router.post("/search")
    def search(payload: dict = Body(...), user=Depends(current_user)):
        allowed = {"group_id", "project_id", "collection", "query", "molecule_type", "sort", "descending", "limit", "offset", "metric_filters", "molecule_ids", "evaluation_filters"}
        if set(payload) - allowed:
            raise HTTPException(422, "Unsupported molecular search fields.")
        return call(MolecularDatabase().list_molecules, user, **payload)

    @router.get("/{molecule_id}")
    def molecule(molecule_id: str, user=Depends(current_user)):
        result = call(MolecularDatabase().list_molecules, user, molecule_ids=[molecule_id])
        if not result["items"]:
            raise HTTPException(404, "Molecular record not found or inaccessible.")
        return result["items"][0]

    @router.patch("/groups/{group_id}/molecules/{molecule_id}")
    def annotate(group_id: str, molecule_id: str, payload: dict = Body(...), user=Depends(current_user)):
        if set(payload) - {"name", "metadata"}:
            raise HTTPException(422, "Only the name and metadata can be changed. Molecular sequences are immutable.")
        call(MolecularDatabase().update_metadata, user, group_id, molecule_id, payload.get("metadata") or {}, name=payload.get("name"))
        return {"status": "updated"}

    app.include_router(router)
