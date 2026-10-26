# REST and MCP integration

The Studio REST API is the shared backend for the current application, LLM tools
and a future student web interface. The MCP adapter calls REST; it does not import
Streamlit, open the molecular database or run scientific programs itself.

## Discovery and identity

- `/docs` and `/openapi.json`: REST documentation and complete OpenAPI, including
  internal machine APIs.
- `GET /api/v1/capabilities`: public operations available to the signed-in role.
  Machine endpoints are excluded.
- `GET /api/v1/operations/schema?method=POST&path=/api/v1/molecules/exports`:
  the operation contract and referenced schemas, filtered by role.
- `GET /api/v1/workflows` and `/workflows/{id}/schema`: scientific adapters and
  parameter descriptors, now also with `parameters_json_schema` and, where
  applicable, `options_json_schema`.
- `GET /api/v1/pipelines/catalog`: node/port/default-parameter definitions.

Other paths below are relative to `/api/v1`. Normal calls use
`Authorization: Bearer <Studio session token>`. Sign in through `/auth/login`;
the existing temporary-password change flow applies. Node credentials remain a
separate internal contract.

## Private file handles

Upload bytes to `POST /files?filename=input.pdb`. The upload streams to private
disk storage and returns a descriptor:

```json
{
  "id": "<file-id>",
  "name": "input.pdb",
  "reference": "studio-file:///<file-id>",
  "download_url": "/api/v1/files/<file-id>"
}
```

Use the reference as a path parameter in a workflow/pipeline request. Studio
resolves it for the authenticated owner and uses existing remote input staging.
The original extension is preserved. Handles remain private even when their
source molecules belong to a shared project.

`GET /files/{id}`, `GET /files/{id}/metadata` and `DELETE /files/{id}` download,
inspect and explicitly remove a private file. Retain inputs while jobs or saved
preparations still use them. Uploading never implicitly imports compounds.

For curated molecular inputs, `/molecules/imports` accepts `file_id`, ID/name/
chain mapping and `review_small_molecules=true`. Reviewed SMILES tables use the
same validation, duplicate report and retained source as the UI. Analyze without
importing through `/molecules/imports/analyze`; retrieve group `/imports`,
`/molecules/imports/{batch_id}/rows` and `/molecules/imports/{batch_id}/source`.

## External libraries and large screening results

`GET /libraries`, `/libraries/{id}` and `/libraries/{id}/preparations` provide
metadata and fingerprint/centroid lineage. Administrators register data-local
files with `POST /libraries`, grant project visibility with
`PATCH /libraries/{id}/access`, and can attach existing completed manifests with
`POST /libraries/{id}/preparations`. Normal library jobs attach their own outputs.
Libraries without an explicit owner/project grant are admin-only.

Edit replicated library metadata on its owning workstation; grants then
synchronize through existing replication. Compounds stay external. The scheduled
`library.prepare` and `library.centroids` workflows retain their existing
data-local and administrator path-registration rules.

For a campaign run:

- `GET /screening/{task_group}/{run_id}` returns a bounded preview, round history,
  convergence summaries, score units and stopping criteria.
- `POST /screening/{task_group}/{run_id}/search` accepts offset, limit (1–1,000),
  text, compound-ID and numeric score filters.
- `GET /screening/{task_group}/{run_id}/files/{kind}` serves complete `scores`,
  `failures`, `rounds`, `convergence`, `manifest` or `model` files.
- `POST /screening/{task_group}/{run_id}/final-selection` with `{"limit": 1000}`
  explicitly imports a final shortlist after successful completion. Only the
  campaign owner/admin may import, and the allowed range stays 0–1,000.

Search preserves original one-based cumulative docking `rank`. Unfiltered pages
skip whole Parquet row groups. Filtered pages scan bounded batches until filled;
`total` is null to avoid forcing a full count. `source_total`, `has_more` and
`next_offset` are returned. Late matches can still require long scans; this is
not a new full compound index or database import.

Local job artifacts also stream through
`GET /jobs/{reference}/artifact?path=artifacts/relative-file`. Both streaming
routes support HTTP ranges, check access/traversal, and have no 512 MiB file
ceiling. These serve files on the owning API host. The older generic remote
bundle/artifact protocol retains its 512 MiB limit. Collected campaign outputs
can be downloaded from the organizer.

## Exports, bound poses and MD

`POST /molecules/exports` accepts molecule IDs and `CSV`, `FASTA`, `SMILES`, `SDF`
or `Structures ZIP`, returning a private file descriptor.

`POST /molecules/{id}/bound-inputs` accepts an optional structure ID. It uses the
saved target/pose relationship and validated chemistry-to-coordinate map to
prepare a matching complex PDB and ligand SDF. Its `workflow_inputs` can be passed
to `ligand.md.openmm` or `ligand.md.gromacs`, plus desired simulation parameters.
Preparation does not start MD. Ambiguous atom maps are rejected; target and
compound remain separate database molecules.

Membership removal, group project reassignment, chain roles, chain cropping,
portable structure snapshots and saved docking preparations also have public
routes. Cropping uses the same function as Streamlit. Existing owner/write checks
are retained.

## Retries and revisions

Send `Idempotency-Key` on workflow, pipeline and design-campaign submissions.
Workflow/pipeline requests can instead use `execution.idempotency_key`; both
values must match when supplied together. Identical retries return the original
response, including direct local execution. Changed payloads return 409. Normal
cluster idempotency remains in place.

A lost reply may leave a pending/uncertain outcome; it is never blindly submitted
again. Inspect `GET /submissions/status?kind=workflow&key=<key>` and your jobs.
Use `kind=pipeline` or `campaign` as appropriate. Records are user-private and
contain a payload hash/response outside the molecular record tables.

Pipeline history is available through `/pipelines/{id}/versions` and
`/pipelines/{id}/versions/{revision}`, with the required `project_id` query.
Scope matches the editor's existing owner/project scope.

## Administration and browser clients

Admin routes cover account creation/list/detail, disable/restore, temporary
password reset, project-grant replacement and existing Studio settings. Reset
revokes sessions and forces a password change; the protections against disabling
oneself/the last admin remain. Workers direct account management to the coordinator.

For a separate-origin browser, set `STUDIO_API_CORS_ORIGINS` to explicit comma-
separated origins before API startup. Authorization, upload, retry and Range
headers are supported; download headers are exposed. A same-origin reverse proxy
also works. Enable production accounts before student network access; this change
does not silently convert or expose the current local development setup.

## MCP adapter

The optional `mcp` extra uses the official Python SDK, tested with 2.3.0. It
supports stdio and Streamable HTTP using the [SDK transport
interfaces](https://py.sdk.modelcontextprotocol.io/run/).

Configure an MCP host to launch the Studio environment's Python with arguments
`-m mn_molecular_design_studio.mcp_server`. Environment:

- `STUDIO_API_URL`: API origin, default `http://127.0.0.1:8766`.
- `STUDIO_API_TOKEN`: signed-in user's session token for stdio.
- `STUDIO_MCP_ALLOW_DEVELOPMENT=1`: explicit local stdio access to the development
  identity. HTTP never accepts that identity.

The application also provides `mn-molecular-design-studio mcp` with
`--transport stdio|streamable-http`, `--host` and `--port` (HTTP default 8767).
These launch the integration service; scientific jobs still use Studio's API,
scheduler and normal results pages.

HTTP checks each caller's Studio bearer token and forwards that identity to REST,
without falling back to a process-wide token. `STUDIO_MCP_ALLOWED_HOSTS` and
`STUDIO_MCP_ALLOWED_ORIGINS` configure host/origin checks; defaults are localhost.
This is Studio's bearer integration, not an OAuth server or a Codex registration
step. An external MCP host needs configured bearer headers or can use stdio.

Named tools cover discovery, operation/workflow schemas, validation, placement,
submission/status/cancellation, retry status, pipelines, library lineage,
screening/convergence, final shortlist import, exports, bound-pose preparation and
text-file staging. `studio_api_request` reaches other discovered public operations
under the same permissions. Job submissions require retry keys. Internal node
routes, external URLs and authentication-header overrides are excluded.

Small files are resources at `studio://files/{id}` and molecular structures at
`studio://structures/{id}`. Inline responses/resources are bounded to 1 MiB;
large files remain available through authenticated REST streaming. This bounds
LLM context, not dataset size.

## Rollout boundaries

This delivers the missing data API and working MCP adapter. The student web
frontend is a separate client still to be built on these contracts. Existing
Streamlit actions continue to use shared Python services and saved form state.
Public OAuth/SSO, per-user resource quotas and event streams are not introduced.

## Initial integration verification on 2026-10-26

- 330 affected API, access-control, pipeline, molecular database, import, docking
  handoff, replication and convergence tests passed.
- JSON Schemas for all 113 workflow adapters passed Draft 2020-12 validation.
- The running API exposes 117 HTTP operations, including its 14 internal machine
  operations; public discovery filters those internal routes out.
- A real stdio MCP subprocess exposed 25 tools and successfully queried the live
  library catalogue, docking schema, MD preparation contract and ranks 1–5 of
  the campaign's 1,590,790 external scores.
- HTTP MCP tests used production instructor/student accounts and verified that
  user identities stay separate. No scientific jobs were launched by these
  integration checks.

## Additional client coverage

Kejia IDR and Peptide Trap use the same discovery, validation, placement,
scheduled submission, results and cancellation contracts as other workflows:
`design.kejia_idr` and `design.peptide_trap`. Their schemas include the canvas
parameter ranges and allowed scaffold/template choices. The design stages use
CPU resources; their downstream prediction nodes use the existing GPU scheduler.

The following preparation and interaction contracts are also available, with
the signed-in user's existing permissions:

| Function | REST operation | Named MCP operation |
| --- | --- | --- |
| Scientific form fields and helper discovery | `GET /workflows/{id}/setup` | `workflow_setup` |
| Installed Kejia/trap templates (paginated) | `GET /workflows/{id}/templates` | `workflow_templates` |
| Template coordinates | `GET /workflows/{id}/templates/{template}/file?library=...` | `studio_api_request` |
| Private copy of a template, with usable parameters | `POST /workflows/{id}/templates/{template}/stage?library=...` | `stage_workflow_template` |
| Normalize Kejia inputs / preview trap matching windows | `POST /workflows/{id}/preview` | `preview_workflow` |
| Stage the Tamarind or C37 example without running it | `POST /workflows/{id}/examples/{example}` | `stage_workflow_example` |
| Structure chain/residue/ligand inventory and quick PDB checks | `POST /structures/preview` | `preview_structure` |
| Canvas starter catalogue | `GET /pipelines/templates` | `pipeline_templates` |
| Editable starter or Kejia/trap continuation graph | `POST /pipelines/drafts` | `create_pipeline_draft` |
| Standalone or shared workstation resources | `GET /workstations` | `list_workstations` |
| Common local/distributed job listing and inspection | `GET /executions`, `/executions/{ref}` | `list_jobs`, `get_job` |
| Available job controls and parameter schemas | `GET /executions/{ref}/actions` | `job_actions` |
| Apply a discovered job control | `POST /executions/{ref}/actions` | `control_job`, `cancel_job` |
| Read saved results / request remote inspection | `GET` or `POST /executions/{ref}/results` | `job_results` |
| Candidate metric fields and bounded tables | `GET /jobs/{ref}/presentation` | `result_presentation` |
| Candidate complex viewer layers and chain roles | `GET /jobs/{ref}/candidates/{id}/viewer` | `candidate_viewer` |
| Molecular viewer layers, including separate bound target and ligand | `GET /molecules/{id}/viewer` | `molecule_viewer` |
| Save representation, coloring, background and camera settings | `POST /visualizations` | `save_visualization` |
| HTML report, SVG figure or CSV table | `POST /jobs/{ref}/reports` | `export_result_report` |
| Collection metadata schema editing | `PUT /molecules/groups/{id}/schema` | `studio_api_request` |
| Local workstation connection/placement settings | `GET`, `PATCH /admin/workstation` | `studio_api_request` (admin) |
| Replication conflict inspection | `GET /admin/replication/conflicts` | `studio_api_request` (admin) |

Job listing returns an unambiguous `reference`: `local:<task-group>:<run-id>` or
`cluster:<job-id>`. A plain local run ID is also supported. Local execution
includes pause/resume when the run has a resumable worker request. Distributed
actions retain the workstation agent's existing support; `retry` specifically
means retrying a failed compatibility acceptance, not rerunning a calculation.
Action discovery is a preview; the backend rechecks current state/resources.
Permanent deletion remains the explicit existing cluster purge operation.

For remote results, `GET` reads the most recent inspection or explains how to
request one. `POST` queues a new inspection and returns its command/poll URL.
Legacy connected workers forward distributed execution controls to their
coordinator using the caller's identity. Standalone installations fall back to
local jobs and local resource discovery.

For example, prepare a Peptide Trap request through
`stage_workflow_example("design.peptide_trap", "c37")`, inspect its matching
windows with `preview_workflow`, then validate/submit the returned parameters
through the normal workflow API. Use `create_pipeline_draft("peptide_trap",
project_id, parameters=...)` to create the editable full pipeline instead.
For a continuation, pass `source_job` and optional `selected_candidate_ids`.
Draft creation neither saves nor starts a pipeline. Save/validate/run the returned
graph with the existing pipeline operations.

Reports use saved convergence summaries (top 100, 1,000 and 10,000 together by
default) or a bounded candidate metric selection. They retain recorded score
units and do not scan/import the complete external compound table. Report/view
files are private file handles. View descriptions supply authenticated structure
URLs and styling for a client renderer; they do not manipulate an open browser's
camera or reproduce every Streamlit chart/widget. Reports and candidate views
are requested from the workstation holding the run; distributed inspection and
artifact commands remain available for remote runs.

Workstation settings never return enrollment tokens or passwords. Connection
changes report that API/agent restart is required. Editing peer addresses does
not enroll a new trusted peer or generate credentials.

The MCP adapter now defines 41 named tools. Other public operations continue to
be discoverable through `studio_operations` / `studio_operation_schema` and
callable through `studio_api_request`. These extensions do not build the student
frontend, register an LLM client, add OAuth/SSO or add per-user quotas.

The subsequent full pytest run, including 25 new API/MCP cases, finished with
**1,334 passed, 15 skipped, zero failures/errors** in 107.88 seconds. See the
[full verification report](api-mcp-pytest-20261026.md) for scope, fixes and skips.
