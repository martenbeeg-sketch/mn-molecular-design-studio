"""Thin MCP adapter: all scientific and data operations go through Studio REST.

The optional MCP dependency is loaded only when this integration is started.
No Streamlit, database or workflow execution modules are imported here.
"""
from __future__ import annotations

import os
from typing import Any
from urllib.parse import quote

from mn_molecular_design_studio.api_client import StudioClient


def create_server(client=None):
    from mcp.server import MCPServer
    from mcp.server.mcpserver import Context
    from mcp_types import ToolAnnotations
    # Annotations are resolved by the SDK in module globals, including when this
    # optional import occurs inside the factory.
    globals()["Context"] = Context
    studio = client or StudioClient(os.environ.get("STUDIO_API_URL", "http://127.0.0.1:8766"),
        token=os.environ.get("STUDIO_API_TOKEN"),
        allow_development=os.environ.get("STUDIO_MCP_ALLOW_DEVELOPMENT", "").lower() in {"1", "true"})
    server = MCPServer("Molecular Design Studio", version="0.1.0", log_level="WARNING", instructions=(
        "Discover workflows and schemas, validate inputs, then submit scheduled jobs. "
        "Reuse the same idempotency key when retrying a submission. Long jobs continue "
        "independently; poll returned job/command IDs. Large libraries and complete "
        "screening results remain external. Use bounded result queries; final molecular "
        "Database imports are limited to an explicit 0–1,000 compounds. Target and "
        "ligand remain separate molecular records. Node endpoints are internal."))
    read = ToolAnnotations(readOnlyHint=True, destructiveHint=False, idempotentHint=True)
    write = ToolAnnotations(readOnlyHint=False, destructiveHint=False)

    async def call(ctx, method, path, *, body=None, query=None, headers=None, content=False, content_bytes=None):
        incoming = ctx.headers
        bearer = None
        if incoming is not None:
            authorization = next((v for k, v in incoming.items() if k.lower() == "authorization"), "")
            if not authorization.startswith("Bearer "):
                from mn_molecular_design_studio.api_client import StudioAPIError
                raise StudioAPIError("A Studio user bearer token is required.")
            bearer = authorization[7:]
        return await studio.request(method, path, body=body, query=query, headers=headers,
                                    bearer=bearer, http_client=incoming is not None, content=content, content_bytes=content_bytes)

    def segment(value):
        return quote(str(value), safe="")

    @server.tool(annotations=read)
    async def studio_identity(ctx: Context) -> dict[str, Any]:
        """Identify the authenticated Studio account and accessible projects."""
        return await call(ctx, "GET", "/api/v1/auth/me")

    @server.tool(annotations=read)
    async def studio_operations(ctx: Context) -> dict[str, Any]:
        """Discover available public REST operations; excludes internal node endpoints."""
        return await call(ctx, "GET", "/api/v1/capabilities")

    @server.tool(annotations=read)
    async def studio_operation_schema(method: str, path: str, ctx: Context) -> dict[str, Any]:
        """Get the full input/output contract for a public REST operation from studio_operations."""
        return await call(ctx, "GET", "/api/v1/operations/schema", query={"method": method, "path": path})

    @server.tool(annotations=read)
    async def list_workflows(ctx: Context, category: str | None = None) -> dict[str, Any]:
        """List scientific workflow IDs, descriptions and parameter descriptors."""
        return {"items": await call(ctx, "GET", "/api/v1/workflows", query={"category": category} if category else {})}

    @server.tool(annotations=read)
    async def workflow_schema(workflow_id: str, ctx: Context) -> dict[str, Any]:
        """Get a workflow's parameters, defaults and executable JSON Schema."""
        return await call(ctx, "GET", f"/api/v1/workflows/{segment(workflow_id)}/schema")

    @server.tool(annotations=read)
    async def workflow_setup(workflow_id: str, ctx: Context) -> dict[str, Any]:
        """Get form fields, pipeline node definitions and available preparation helpers/examples for any workflow."""
        return await call(ctx, "GET", f"/api/v1/workflows/{segment(workflow_id)}/setup")

    @server.tool(annotations=read)
    async def workflow_templates(workflow_id: str, ctx: Context, library: str = "", offset: int = 0, limit: int = 100) -> dict[str, Any]:
        """Browse installed Kejia IDR or Peptide Trap templates without submitting a calculation."""
        return await call(ctx, "GET", f"/api/v1/workflows/{segment(workflow_id)}/templates",
                          query={"library": library, "offset": offset, "limit": limit})

    @server.tool(annotations=write)
    async def stage_workflow_template(workflow_id: str, template_id: str, ctx: Context, library: str = "") -> dict[str, Any]:
        """Copy an installed design template to private input storage and return usable workflow parameters."""
        return await call(ctx, "POST", f"/api/v1/workflows/{segment(workflow_id)}/templates/{segment(template_id)}/stage",
                          query={"library": library})

    @server.tool(annotations=read)
    async def preview_workflow(workflow_id: str, parameters: dict, ctx: Context, offset: int = 0, limit: int = 100) -> dict[str, Any]:
        """Preview normalized Kejia inputs or matching peptide-trap sequence windows; no docking or design starts."""
        return await call(ctx, "POST", f"/api/v1/workflows/{segment(workflow_id)}/preview",
                          body={"parameters": parameters, "offset": offset, "limit": limit})

    @server.tool(annotations=write)
    async def stage_workflow_example(workflow_id: str, example_id: str, ctx: Context) -> dict[str, Any]:
        """Stage an available workflow example as private file handles and parameters; does not run it."""
        return await call(ctx, "POST", f"/api/v1/workflows/{segment(workflow_id)}/examples/{segment(example_id)}", body={})

    @server.tool(annotations=read)
    async def preview_structure(file_id: str, ctx: Context) -> dict[str, Any]:
        """Inspect an uploaded structure's chains, residues and ligands, plus the quick PDB receptor checks."""
        return await call(ctx, "POST", "/api/v1/structures/preview", body={"file_id": file_id})

    @server.tool(annotations=read)
    async def list_workstations(ctx: Context) -> dict[str, Any]:
        """Inspect all configured workstations, their CPU/GPU resources and availability."""
        return await call(ctx, "GET", "/api/v1/workstations")

    @server.tool(annotations=read)
    async def validate_workflow(request: dict, ctx: Context) -> dict[str, Any]:
        """Validate a full workflow request, including project_id, without submitting it."""
        return await call(ctx, "POST", "/api/v1/workflows/validate", body=request)

    @server.tool(annotations=read)
    async def workflow_placement(request: dict, ctx: Context) -> dict[str, Any]:
        """Preview eligible computers and resources without submitting a job."""
        return await call(ctx, "POST", "/api/v1/workflows/placement-options", body=request)

    @server.tool(annotations=write)
    async def submit_workflow(request: dict, idempotency_key: str, ctx: Context) -> dict[str, Any]:
        """Submit a validated workflow. Reuse this key for the same logical submission."""
        return await call(ctx, "POST", "/api/v1/workflows", body=request, headers={"Idempotency-Key": idempotency_key})

    @server.tool(annotations=read)
    async def submission_status(kind: str, idempotency_key: str, ctx: Context) -> dict[str, Any]:
        """Inspect a workflow, pipeline or campaign submission after a lost response, without resubmitting."""
        return await call(ctx, "GET", "/api/v1/submissions/status", query={"kind": kind, "key": idempotency_key})

    @server.tool(annotations=read)
    async def list_jobs(ctx: Context, project_id: str | None = None, distributed: bool = True,
                        limit: int = 50, offset: int = 0) -> dict[str, Any]:
        """List accessible jobs, optionally across all Studio workstations."""
        params = {"limit": limit, "offset": offset}
        if project_id:
            params["project_id"] = project_id
        params["distributed"] = distributed
        return await call(ctx, "GET", "/api/v1/executions", query=params)

    @server.tool(annotations=read)
    async def get_job(reference: str, ctx: Context) -> dict[str, Any]:
        """Inspect a job; use cluster:<job-id> for a distributed job, otherwise its local run ID."""
        return await call(ctx, "GET", "/api/v1/executions/" + segment(reference))

    @server.tool(annotations=read)
    async def job_results(reference: str, ctx: Context) -> dict[str, Any]:
        """Read local results or request remote inspection, returning a command ID to poll."""
        return await call(ctx, "POST" if reference.startswith("cluster:") else "GET",
                          "/api/v1/executions/" + segment(reference) + "/results")

    @server.tool(annotations=read)
    async def job_actions(reference: str, ctx: Context) -> dict[str, Any]:
        """Discover available job actions and their parameter schemas, including restricted compatibility retry."""
        return await call(ctx, "GET", "/api/v1/executions/" + segment(reference) + "/actions")

    @server.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=True))
    async def control_job(reference: str, action: str, ctx: Context, parameters: dict | None = None) -> dict[str, Any]:
        """Perform an action discovered through job_actions under the current user's job permissions."""
        return await call(ctx, "POST", "/api/v1/executions/" + segment(reference) + "/actions",
                          body={"action": action, "parameters": parameters or {}})

    @server.tool(annotations=read)
    async def command_status(command_id: str, ctx: Context) -> dict[str, Any]:
        """Poll an inspection, file transfer or other workstation command."""
        return await call(ctx, "GET", "/api/v1/cluster/commands/" + segment(command_id))

    @server.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=True, idempotentHint=True))
    async def cancel_job(reference: str, ctx: Context) -> dict[str, Any]:
        """Cancel an accessible active job; cluster IDs use the cluster: prefix."""
        return await call(ctx, "POST", "/api/v1/executions/" + segment(reference) + "/actions", body={"action": "cancel"})

    @server.tool(annotations=read)
    async def list_libraries(ctx: Context, project_id: str | None = None) -> dict[str, Any]:
        """List metadata for accessible external libraries; no compound tables are imported."""
        return await call(ctx, "GET", "/api/v1/libraries", query={"project_id": project_id} if project_id else {})

    @server.tool(annotations=read)
    async def library_preparations(library_id: str, ctx: Context, kind: str = "all") -> dict[str, Any]:
        """List fingerprint or centroid preparations linked to a source library."""
        return await call(ctx, "GET", "/api/v1/libraries/" + segment(library_id) + "/preparations", query={"kind": kind})

    @server.tool(annotations=read)
    async def screening_results(task_group: str, run_id: str, ctx: Context, offset: int = 0,
                                limit: int = 100, query: str = "", score_min: float | None = None,
                                score_max: float | None = None) -> dict[str, Any]:
        """Query a bounded page of external scores, retaining original cumulative docking rank."""
        return await call(ctx, "POST", f"/api/v1/screening/{segment(task_group)}/{segment(run_id)}/search",
                          body={"offset": offset, "limit": limit, "query": query,
                                "score_min": score_min, "score_max": score_max})

    @server.tool(annotations=read)
    async def screening_convergence(task_group: str, run_id: str, ctx: Context) -> dict[str, Any]:
        """Get round history, convergence summaries and the configured stopping criteria."""
        return await call(ctx, "GET", f"/api/v1/screening/{segment(task_group)}/{segment(run_id)}")

    @server.tool(annotations=write)
    async def import_final_shortlist(task_group: str, run_id: str, limit: int, ctx: Context) -> dict[str, Any]:
        """Explicitly promote 0–1,000 final ranked compounds after successful completion."""
        return await call(ctx, "POST", f"/api/v1/screening/{segment(task_group)}/{segment(run_id)}/final-selection", body={"limit": limit})

    @server.tool(annotations=write)
    async def export_molecules(molecule_ids: list[str], ctx: Context, format: str = "CSV") -> dict[str, Any]:
        """Create a private CSV/FASTA/SMILES/SDF/structures-ZIP export and return its file handle."""
        return await call(ctx, "POST", "/api/v1/molecules/exports", body={"molecule_ids": molecule_ids, "format": format})

    @server.tool(annotations=write)
    async def prepare_bound_pose(molecule_id: str, ctx: Context, structure_id: str | None = None) -> dict[str, Any]:
        """Prepare the linked receptor and mapped ligand as MD inputs without starting MD or merging DB molecules."""
        return await call(ctx, "POST", f"/api/v1/molecules/{segment(molecule_id)}/bound-inputs", body={"structure_id": structure_id})

    @server.tool(annotations=write)
    async def stage_text_input(filename: str, text: str, ctx: Context) -> dict[str, Any]:
        """Stage a small SMILES, FASTA, PDB, SDF or parameter file and return its private workflow input reference."""
        return await call(ctx, "POST", "/api/v1/files", query={"filename": filename},
                          headers={"Content-Type": "text/plain"}, content_bytes=text.encode())

    @server.tool(annotations=read)
    async def pipeline_templates(ctx: Context) -> dict[str, Any]:
        """List the same editable starter pipelines offered by the Studio canvas."""
        return await call(ctx, "GET", "/api/v1/pipelines/templates")

    @server.tool(annotations=read)
    async def create_pipeline_draft(template: str, project_id: str, ctx: Context, parameters: dict | None = None,
                                    source_job: str | None = None, selected_candidate_ids: list[str] | None = None) -> dict[str, Any]:
        """Create an editable starter graph or continue completed Kejia/Peptide Trap candidates; does not save or run it."""
        return await call(ctx, "POST", "/api/v1/pipelines/drafts", body={"template": template, "project_id": project_id,
                          "parameters": parameters or {}, "source_job": source_job,
                          "selected_candidate_ids": selected_candidate_ids or []})

    @server.tool(annotations=read)
    async def molecule_viewer(molecule_id: str, ctx: Context, structure_id: str | None = None) -> dict[str, Any]:
        """Get structure-view layers, including the separate saved target and bound ligand when available."""
        return await call(ctx, "GET", f"/api/v1/molecules/{segment(molecule_id)}/viewer",
                          query={"structure_id": structure_id} if structure_id else {})

    @server.tool(annotations=read)
    async def candidate_viewer(reference: str, candidate_id: str, ctx: Context) -> dict[str, Any]:
        """Get an authorized complex structure URL and binder/target chains for a local job candidate."""
        return await call(ctx, "GET", f"/api/v1/jobs/{segment(reference)}/candidates/{segment(candidate_id)}/viewer")

    @server.tool(annotations=read)
    async def result_presentation(reference: str, ctx: Context, offset: int = 0, limit: int = 100) -> dict[str, Any]:
        """Page through local result candidates and discover metric fields for plots or reports."""
        return await call(ctx, "GET", f"/api/v1/jobs/{segment(reference)}/presentation", query={"offset": offset, "limit": limit})

    @server.tool(annotations=write)
    async def save_visualization(request: dict, ctx: Context) -> dict[str, Any]:
        """Save a private structure-view description with representation, coloring, background and camera settings."""
        return await call(ctx, "POST", "/api/v1/visualizations", body=request)

    @server.tool(annotations=write)
    async def export_result_report(reference: str, ctx: Context, format: str = "HTML", metric: str = "mean_score",
                                    top_counts: list[int] | None = None, limit: int = 1000) -> dict[str, Any]:
        """Export convergence curves (top 100/1,000/10,000 together) or candidate metrics as HTML, SVG or CSV."""
        return await call(ctx, "POST", f"/api/v1/jobs/{segment(reference)}/reports", body={"format": format,
                          "metric": metric, "top_counts": top_counts or [100, 1000, 10000], "limit": limit})

    @server.tool(annotations=read)
    async def validate_pipeline(request: dict, ctx: Context) -> dict[str, Any]:
        """Validate an editable pipeline graph, its connections and inputs without running it."""
        return await call(ctx, "POST", "/api/v1/pipelines/validate", body=request)

    @server.tool(annotations=write)
    async def run_pipeline(request: dict, idempotency_key: str, ctx: Context) -> dict[str, Any]:
        """Submit a Studio pipeline; reuse this key when retrying the same submission."""
        return await call(ctx, "POST", "/api/v1/pipelines/run", body=request, headers={"Idempotency-Key": idempotency_key})

    @server.tool(annotations=ToolAnnotations(readOnlyHint=False))
    async def studio_api_request(method: str, path: str, ctx: Context,
                                  body: dict | None = None, query: dict | None = None,
                                  idempotency_key: str | None = None) -> dict[str, Any]:
        """Call another discovered public operation under the current user's permissions.

        Does not allow external URLs, arbitrary headers or internal node routes.
        For workflow/pipeline submission prefer the named tools with retry keys.
        """
        from mn_molecular_design_studio.api_client import public_path
        selected = public_path(path)
        if method.upper() == "POST" and selected in {"/api/v1/workflows", "/api/v1/pipelines/run", "/api/v1/campaigns", "/api/v1/cluster/submissions"}:
            execution = (body or {}).get("execution") or {}
            if not isinstance(execution, dict):
                raise ValueError("execution must be an object.")
            if not idempotency_key and not execution.get("idempotency_key") and not (body or {}).get("idempotency_key"):
                raise ValueError("A job submission needs an idempotency key. Prefer the named submission tools.")
            if selected == "/api/v1/cluster/submissions" and idempotency_key:
                if (body or {}).get("idempotency_key") not in {None, idempotency_key}:
                    raise ValueError("Submission idempotency keys must match.")
                body = {**(body or {}), "idempotency_key": idempotency_key}
        result = await call(ctx, method, selected, body=body, query=query,
                            headers={"Idempotency-Key": idempotency_key} if idempotency_key else None)
        return result if isinstance(result, dict) else {"items": result}

    @server.resource("studio://files/{file_id}")
    async def file_resource(file_id: str, ctx: Context) -> bytes:
        """Read a small private input/export file; large files use REST streaming."""
        return await call(ctx, "GET", f"/api/v1/files/{segment(file_id)}", content=True)

    @server.resource("studio://structures/{structure_id}")
    async def structure_resource(structure_id: str, ctx: Context) -> str:
        """Read an accessible molecular structure for visualization or inspection."""
        content = await call(ctx, "GET", f"/api/v1/molecules/structures/{segment(structure_id)}/file", content=True)
        return content.decode("utf-8")

    return server


def http_app(server=None, client=None):
    """Bearer-authenticated HTTP transport; no process-wide user token fallback."""
    from starlette.responses import JSONResponse
    from mcp.server.transport_security import TransportSecuritySettings
    studio = client or StudioClient(os.environ.get("STUDIO_API_URL", "http://127.0.0.1:8766"))
    server = server or create_server(studio)
    hosts = os.environ.get("STUDIO_MCP_ALLOWED_HOSTS", "127.0.0.1,127.0.0.1:*,localhost,localhost:*,[::1],[::1]:*")
    origins = os.environ.get("STUDIO_MCP_ALLOWED_ORIGINS", "http://127.0.0.1:*,http://localhost:*")
    application = server.streamable_http_app(stateless_http=True, transport_security=TransportSecuritySettings(
        allowed_hosts=[v.strip() for v in hosts.split(",") if v.strip()],
        allowed_origins=[v.strip() for v in origins.split(",") if v.strip()]))

    class AuthenticatedApp:
        async def __call__(self, scope, receive, send):
            if scope["type"] == "http":
                headers = dict(scope.get("headers") or [])
                authorization = headers.get(b"authorization", b"").decode()
                try:
                    if not authorization.startswith("Bearer "):
                        raise ValueError("A Studio bearer token is required.")
                    # This also rejects the local development admin identity.
                    await studio.request("GET", "/api/v1/auth/me", bearer=authorization[7:], http_client=True)
                except Exception:
                    # Never expose tokens, upstream connection details or traceback text.
                    response = JSONResponse({"detail": "Sign in to a production Studio API with a valid user token."},
                                            status_code=401, headers={"WWW-Authenticate": "Bearer"})
                    return await response(scope, receive, send)
            await application(scope, receive, send)

    return AuthenticatedApp()


def run(transport="stdio", host="127.0.0.1", port=8767):
    server = create_server()
    if transport == "stdio":
        server.run()
    elif transport == "streamable-http":
        import uvicorn
        uvicorn.run(http_app(server), host=host, port=port, access_log=False)
    else:
        raise ValueError("Choose stdio or streamable-http.")


if __name__ == "__main__":
    run()
