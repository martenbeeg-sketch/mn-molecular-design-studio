"""REST projection of Projects and the saved Data browser."""
from __future__ import annotations

from typing import Literal

from fastapi import APIRouter, Depends, HTTPException, Query, Response
from pydantic import BaseModel, Field, StrictBool


class ProjectRequest(BaseModel):
    model_config = {"extra": "forbid"}
    name: str
    description: str = ""
    color: str = "#ef5350"


class ProjectArchive(BaseModel):
    model_config = {"extra": "forbid"}
    archived: StrictBool


class CleanupReview(BaseModel):
    model_config = {"extra": "forbid"}
    mode: Literal["clean", "delete"]
    selected_runs: list[str] = Field(default_factory=list, max_length=10000)
    include_collections: StrictBool = False


def visible_projects(user, *, include_archived=True):
    from mn_molecular_design_studio.auth import accessible_project_ids
    from mn_molecular_design_studio.workspace import load_organization

    allowed = set(accessible_project_ids(user))
    return sorted((project for project in load_organization().get("projects", [])
                   if isinstance(project, dict) and project.get("id") in allowed and not project.get("deleted_at")
                   and (include_archived or not project.get("archived_at"))),
                  key=lambda row: (row["name"].casefold(), row["id"]))


def browse(user, *, folder, ancestors, branch, query, historical, structures, offset, limit):
    from mn_molecular_design_studio.candidate_sets import json_safe
    from mn_molecular_design_studio.data_browser import ReadOnlyLibrary, expand_scope, listed_results, scope_path
    from mn_molecular_design_studio.data_browser_table import accessible_entries, result_table
    from mn_molecular_design_studio.database_browser_cache import BrowserCache
    from mn_molecular_design_studio.scientific_browser_tree import scientific_branches, scientific_scopes

    database = ReadOnlyLibrary()
    if not database.path.is_file():
        return {"folders": [], "branches": [], "items": [], "total": 0, "next_offset": None}
    projects = visible_projects(user)
    cache = BrowserCache(database.path)
    try:
        index = cache.get_index(database, user, scope="all", collection="all" if historical else "library",
                                hide_smoke_test=False, include_related=True, hide_pipeline_batches=True)
        index["project_names"] = {project["id"]: project["name"] for project in projects}
        scopes = scientific_scopes(index, projects)
        for parent in ancestors:
            if parent not in scopes:
                raise FileNotFoundError("Data browser folder not found.")
            expand_scope(scopes, index, parent)
        if not folder:
            return {"folders": [{"id": scope["id"], "label": scope["label"]} for scope in scopes.values() if not scope["parent"]],
                    "branches": [], "items": [], "total": 0, "next_offset": None}
        if folder not in scopes:
            raise FileNotFoundError("Data browser folder not found or inaccessible. Include its ancestor IDs when opening a nested folder.")
        scope = scopes[folder]
        expand_scope(scopes, index, folder)
        path = scope_path(scopes, folder)
        outcomes = cache.get_saved_results(database, index, scope["subject"], run_view=scope["run_view"],
            include_references=any(scopes[parent].get("route") == "collection" for parent in path),
            fallback_references=(scope.get("entity") or {}).get("kind") == "target" or scope.get("route") == "scientific-project")
        outcomes = [row for row in outcomes if row["outcome_kind"] not in {"Imported molecule", "Reference molecule"}]
        branches, selections = scientific_branches(index, scope, outcomes)
        if branch and branch not in selections:
            raise FileNotFoundError("Data browser result branch not found or inaccessible.")
        entries = accessible_entries(index, listed_results(outcomes, selection=selections.get(branch), query=query, structures=structures))
        page = entries[offset:offset + limit]
        table = result_table(database, index, page, cache=cache)
        items = [{"id": entry["id"], "label": entry["label"], "structure_id": entry["structure_id"],
                  **{key: entry["result"].get(key) for key in ("molecule_id", "group_id", "evaluation_id")},
                  "row": json_safe(row)} for entry, row in zip(page, table)]
        return {"path": [{"id": item, "label": scopes[item]["label"]} for item in path], "branches": branches,
                "folders": [{"id": row["id"], "label": row["label"]} for row in scopes.values() if row["parent"] == folder],
                "items": items, "total": len(entries), "next_offset": offset + limit if offset + limit < len(entries) else None}
    finally:
        cache.close()


def install(app, current_user, admin_user, call):
    router = APIRouter(prefix="/api/v1", tags=["workspace"])

    @router.get("/projects")
    def projects(include_archived: bool = True, user=Depends(current_user)):
        return {"items": call(visible_projects, user, include_archived=include_archived)}

    @router.get("/projects/{project_id}")
    def project(project_id: str, user=Depends(current_user)):
        selected = next((row for row in call(visible_projects, user) if row["id"] == project_id), None)
        if selected is None:
            raise HTTPException(404, "Project not found or inaccessible.")
        return selected

    @router.post("/projects", status_code=201)
    def create_project(payload: ProjectRequest, user=Depends(current_user)):
        from mn_molecular_design_studio.workspace import create_project_for_user
        return call(create_project_for_user, user, payload.name, payload.description, payload.color)

    @router.patch("/projects/{project_id}/archive")
    def project_archive(project_id: str, payload: ProjectArchive, user=Depends(current_user)):
        from mn_molecular_design_studio.project_lifecycle import archive_project
        project(project_id, user)
        return call(archive_project, user, project_id, payload.archived)

    @router.post("/admin/projects/{project_id}/cleanup/reviews", status_code=201)
    def cleanup_review(project_id: str, payload: CleanupReview, user=Depends(admin_user)):
        from mn_molecular_design_studio.project_lifecycle import ProjectCoordinator
        return call(ProjectCoordinator().prepare, user, project_id, **payload.model_dump())

    @router.get("/admin/project-cleanup/{operation_id}")
    def cleanup_status(operation_id: str, user=Depends(admin_user)):
        from mn_molecular_design_studio.project_lifecycle import ProjectNode
        try:
            return ProjectNode().read("operations", operation_id, user)
        except FileNotFoundError as exc:
            raise HTTPException(404, "Project cleanup not found or inaccessible.") from exc
        except PermissionError as exc:
            raise HTTPException(404, "Project cleanup not found or inaccessible.") from exc

    @router.get("/molecules/browser/results")
    def browse_results(folder: str = "", ancestor: list[str] | None = Query(None, max_length=12), branch: str = "",
                       query: str = "", historical: bool = False, structures: bool = False,
                       offset: int = Query(0, ge=0), limit: int = Query(25, ge=1, le=250),
                       export: Literal["csv"] | None = None, user=Depends(current_user)):
        try:
            result = browse(user, folder=folder, ancestors=ancestor or [], branch=branch, query=query,
                            historical=historical, structures=structures, offset=offset, limit=limit)
        except FileNotFoundError as exc:
            raise HTTPException(404, str(exc)) from exc
        if export == "csv":
            from mn_molecular_design_studio.data_browser_table import table_csv
            return Response(table_csv([item["row"] for item in result["items"]]), media_type="text/csv",
                            headers={"Content-Disposition": 'attachment; filename="data-browser-page.csv"'})
        return result

    app.include_router(router)
