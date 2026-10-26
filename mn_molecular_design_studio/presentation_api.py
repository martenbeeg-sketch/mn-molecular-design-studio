"""Structure viewer descriptions and saved figure/report exports for API clients."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Annotated, Literal
from urllib.parse import quote, urlencode

from fastapi import APIRouter, Depends, HTTPException, Query
from fastapi.responses import FileResponse
from pydantic import BaseModel, ConfigDict, Field, model_validator


class Contract(BaseModel):
    model_config = ConfigDict(extra="forbid")


class ViewerOptions(Contract):
    representation: Literal["cartoon", "ball-and-stick", "cartoon+ball-and-stick", "surface"] = "cartoon"
    color: Literal["chain-id", "uniform", "element", "plddt"] = "chain-id"
    background: str = Field(default="#ffffff", pattern=r"^#[0-9a-fA-F]{6}$")
    show_hydrogens: bool = False
    camera_position: list[Annotated[float, Field(allow_inf_nan=False)]] | None = Field(default=None, min_length=3, max_length=3)
    camera_target: list[Annotated[float, Field(allow_inf_nan=False)]] | None = Field(default=None, min_length=3, max_length=3)


class SavedView(Contract):
    file_id: str | None = None
    structure_id: str | None = None
    options: ViewerOptions = Field(default_factory=ViewerOptions)

    @model_validator(mode="after")
    def one_source(self):
        if bool(self.file_id) == bool(self.structure_id):
            raise ValueError("Choose either a private file or a database structure.")
        return self


class Report(Contract):
    format: Literal["CSV", "SVG", "HTML"] = "HTML"
    top_counts: list[Annotated[int, Field(ge=1, strict=True)]] = Field(default_factory=lambda: [100, 1000, 10000], min_length=1, max_length=20)
    metric: str = Field(default="mean_score", max_length=100)
    limit: int = Field(default=1000, ge=1, le=10000)


def install(app, current_user, *, resolve_job, call):
    from mn_molecular_design_studio import api_files
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path
    router = APIRouter(prefix="/api/v1", tags=["visualization and reports"])
    def call(function, *args, **kwargs):
        try:
            return function(*args, **kwargs)
        except (PermissionError, FileNotFoundError) as exc:
            raise HTTPException(404, "Result not found or inaccessible.") from exc
        except (ValueError, TypeError, OSError) as exc:
            raise HTTPException(422, str(exc)) from exc

    def selected_structure(database, user, molecule_id, structure_id):
        rows = call(database.list_molecules, user, molecule_ids=[molecule_id], limit=1, collection="all")["items"]
        if not rows:
            raise HTTPException(404, "Molecule not found or inaccessible.")
        record = rows[0]
        structures = record.get("representations") or record.get("structures", [])
        selected = next((s for s in structures if not structure_id or s["id"] == structure_id), None)
        if not selected:
            raise HTTPException(404, "Structure not found or inaccessible.")
        return record, selected

    @router.post("/visualizations", status_code=201)
    def save_view(payload: SavedView, user=Depends(current_user)):
        if payload.file_id:
            record = call(api_files.describe, user, payload.file_id)
            url = record["download_url"]
        else:
            call(MolecularDatabase().structure_file, user, payload.structure_id)
            url = "/api/v1/molecules/structures/" + quote(payload.structure_id, safe="") + "/file"
        manifest = {"schema_version": 1, "layers": [{"url": url, "requires_authentication": True}],
                    "options": payload.options.model_dump()}
        saved = call(api_files.store_bytes, user, json.dumps(manifest).encode(), name="structure-view.json", media_type="application/json")
        return {"view": manifest, "file": saved}

    @router.get("/molecules/{molecule_id}/viewer")
    def molecule_view(molecule_id: str, structure_id: str | None = None, user=Depends(current_user)):
        from mn_molecular_design_studio.bound_docking import is_bound_pose, receptor_file
        database = MolecularDatabase()
        record, selected = selected_structure(database, user, molecule_id, structure_id)
        source = call(database.structure_file, user, selected["id"])
        layers = [{"molecule_id": molecule_id, "structure_id": selected["id"], "label": record["name"],
                   "url": "/api/v1/molecules/structures/" + quote(selected["id"], safe="") + "/file",
                   "representation": "ball-and-stick" if is_bound_pose(selected) else "cartoon",
                   "format": selected.get("format") or source.suffix.lstrip("."),
                   "requires_authentication": True}]
        if is_bound_pose(selected):
            context = selected.get("context") or {}
            target_id = context.get("target_structure_id")
            if target_id:
                call(database.structure_file, user, target_id)
                url = "/api/v1/molecules/structures/" + quote(target_id, safe="") + "/file"
            else:
                receptor = call(receptor_file, database, user, selected)
                if receptor is None:
                    raise HTTPException(404, "The saved target structure is unavailable.")
                url = "/api/v1/molecules/" + quote(molecule_id, safe="") + "/viewer/receptor?" + urlencode({"structure_id": selected["id"]})
            layers.insert(0, {"molecule_id": context.get("target_molecule_id"), "structure_id": target_id,
                              "label": "Docking target", "url": url, "representation": "cartoon",
                              "requires_authentication": True})
        return {"schema_version": 1, "layers": layers, "options": ViewerOptions().model_dump(),
                "components": record.get("components", []), "coordinates_aligned": is_bound_pose(selected)}

    @router.get("/molecules/{molecule_id}/viewer/receptor")
    def viewer_receptor(molecule_id: str, structure_id: str, user=Depends(current_user)):
        from mn_molecular_design_studio.bound_docking import is_bound_pose, receptor_file
        database = MolecularDatabase()
        _, selected = selected_structure(database, user, molecule_id, structure_id)
        if not is_bound_pose(selected):
            raise HTTPException(404, "No saved target is linked to this structure.")
        path = call(receptor_file, database, user, selected)
        if path is None:
            raise HTTPException(404, "The saved target structure is unavailable.")
        return FileResponse(path, filename=path.name)

    @router.get("/jobs/{reference}/presentation")
    def job_presentation(reference: str, offset: int = Query(0, ge=0), limit: int = Query(100, ge=1, le=1000),
                         user=Depends(current_user)):
        row = resolve_job(reference, user)
        candidates = call(read_candidates, Path(row["run_dir"]))
        items = [{"candidate_id": c["candidate_id"], "stage": c.get("stage"), "metrics": c.get("metrics", {}),
                  "binder_chains": c.get("binder_chains", []), "target_chains": c.get("target_chains", [])}
                 for c in candidates[offset:offset + limit]]
        return {"items": items, "total": len(candidates), "offset": offset, "limit": limit,
                "metric_fields": sorted({key for c in candidates for key in c.get("metrics", {})}),
                "reports_url": "/api/v1/jobs/" + quote(reference, safe="") + "/reports"}

    @router.get("/jobs/{reference}/candidates/{candidate_id}/viewer")
    def candidate_view(reference: str, candidate_id: str, user=Depends(current_user)):
        job = resolve_job(reference, user)
        root = Path(job["run_dir"]).resolve()
        candidate = next((row for row in call(read_candidates, root) if row["candidate_id"] == candidate_id), None)
        if not candidate:
            raise HTTPException(404, "Candidate not found.")
        source = candidate.get("complex_pdb") or candidate.get("binder_pdb")
        path = resolve_stored_path(source, run_dir=root, must_exist=True) if source else None
        if path is None or not path.resolve().is_relative_to(root):
            raise HTTPException(404, "Candidate structure is unavailable in this job's artifacts.")
        url = "/api/v1/jobs/" + quote(reference, safe="") + "/artifact?" + urlencode({"path": str(path.resolve().relative_to(root))})
        return {"schema_version": 1, "candidate_id": candidate_id, "layers": [{"url": url, "requires_authentication": True,
                "binder_chains": candidate.get("binder_chains", []), "target_chains": candidate.get("target_chains", [])}],
                "format": path.suffix.lstrip("."), "options": ViewerOptions().model_dump()}

    @router.post("/jobs/{reference}/reports", status_code=201)
    def report(reference: str, payload: Report, user=Depends(current_user)):
        from mn_molecular_design_studio.bulk_screening_results import screening_results_payload
        from mn_molecular_design_studio.result_presentation import convergence_table, render_report
        job = resolve_job(reference, user)
        root = Path(job["run_dir"])
        screening = call(screening_results_payload, root, preview_limit=1)
        if screening is not None:
            rows = convergence_table(screening, payload.top_counts)
            if payload.metric not in {"mean_score", "boundary_score", "maximum_rank_score_change"}:
                raise HTTPException(422, "Choose mean_score, boundary_score or maximum_rank_score_change.")
            title = "Docking score convergence"
            x, y, series, units = "round", payload.metric, "top_n", screening.get("score_units") or ""
        else:
            candidates = call(read_candidates, root)
            rows = [{"candidate_index": index, "candidate_id": c["candidate_id"], **c.get("metrics", {})}
                    for index, c in enumerate(candidates[:payload.limit], 1)]
            if payload.format != "CSV" and not any(isinstance(r.get(payload.metric), (float, int)) for r in rows):
                raise HTTPException(422, "Choose a numeric metric from /presentation; CSV exports all available metrics.")
            title = "Candidate metrics"
            x, y, series, units = "candidate_index", payload.metric, "series", ""
        content, media, extension = call(render_report, rows, format=payload.format, title=title, x=x, y=y, series=series, units=units)
        file = call(api_files.store_bytes, user, content, name="results-report." + extension, media_type=media)
        return {"file": file, "row_count": len(rows), "scope": "saved convergence summaries" if screening else "bounded candidate selection",
                "truncated": screening is None and len(candidates) > payload.limit}

    app.include_router(router)
