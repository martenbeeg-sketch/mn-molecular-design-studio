# Molecular Design Studio functional audit

November 19, 2026 · baseline `79d301db644a61d5406223157ea1c37909b217ce`

## Assessment

The main weaknesses are at the boundaries between the canvas, REST, the job
runner, and saved scientific data. A stale editor can overwrite a newer draft,
an execution can inherit another execution's logical Run receipt when IDs
collide, and comparisons can rank results with incomplete execution provenance.
Peer naming replay also reproduces a previously reported transaction lock.
These findings warrant targeted fixes before treating the UI and REST as
interchangeable clients for the same work.

The baseline registers **123 REST workflows**, **79 workflows referenced by
the All tools allowlist**, and **96 canvas node kinds**, of which **49 distinct
workflow IDs** are directly referenced by a node. All allowlisted and
node-referenced workflow IDs resolved in the REST catalog. The numbers describe
different entities because the UI also groups modes and expands engine profiles
into cards, and the canvas includes input, merge and output nodes. Authenticated
administrator discovery returned **141 operations on 126 paths**. The inventory
found no duplicate OpenAPI operation IDs. Only two workflows advertised
specialized preparation helpers from the generic setup service.

## Reproduced findings

Priority P1 means a correctness or data-preservation issue with material user
impact. P2 means a functional gap, misleading status, or avoidable manual
workaround. The evidence comes from disposable fixtures and source review.

| ID | Priority | Evidence and impact | Recommended correction |
| --- | --- | --- | --- |
| STUDIO-01 | P1 | Two editors saved different edits from the same pipeline draft. Both saves returned 200 and the stale edit became the current draft. Attempting to submit the original revision returned 422. The canvas autosaves the whole graph, so two open tabs suffice. `pipeline_graph_ui.py:1091`, `api.py:934`, `pipeline_graph.py:2384`. | Require the editor's saved base revision on updates, compare and replace under one lock, and return 409 with a reload or Save as option. Keep immutable run snapshots separate from the editable draft. |
| STUDIO-02 | P1 | Two saved fixture evaluations had the same metric definition, target, chain scope and parameter digest. One had a recorded image digest, the other had none. REST comparison returned 200 with ranks 1 and 2 and zero warnings. The image check discards the empty value before comparing. This is an extension of the earlier provenance concern beyond the corrected check for differing target revisions. `comparisons.py:161`, `molecular_database_api.py:264`. | Use the same compatibility gate in preview, save and UI. Require verified equivalent execution provenance for pooled ranking or show unknown-provenance rows separately. Keep a descriptive table available. |
| STUDIO-03 | P2 | A saved-MD pipeline without a source returned 422 in local validation. With a simulated peer and automatic placement the same draft returned 200 with `valid: true` and `validation_scope: request_schema`. The response declares its scope, but labels a missing required selector as valid and gives clients no comparable readiness result. The worker checks required inputs later. `api.py:1060`, `cluster/service.py:31`, `cluster/service.py:132`, `cluster/agent.py:904`. | Check workstation-independent requirements in both modes. Return separate graph, source and workstation readiness states, with a clear deferred state when remote files must be checked on their owner. |
| STUDIO-04 | P2 | Submit the same incomplete MD draft with an idempotency key. The first response is a definitive 422 before queueing. A replay returns 409, and submission status reports `uncertain`. The ledger uses that state for every exception. A user can be sent looking for a job that never existed. `api.py:1080`, `api_idempotency.py:64`. | Distinguish pre-dispatch rejection from a failure that can occur after dispatch. Preserve the existing safeguard against retrying genuinely uncertain submissions, and give clients a definitive validation outcome. |
| STUDIO-05 | P2 | An authorized administrator could list a draft and fetch its saved version before archive. Both reads returned 422 after archive, while the saved files were still present. The GET routes use a helper that validates a project for new work. `api.py:1053`, `api.py:1350`, `workspace.py:173`. | Authorize reads on accessible archived projects. Keep save and execution subject to the writable-project check. Show archived drafts read-only in the editor. |
| STUDIO-06 | P2 | Appending records to an existing reusable dataset with a valid new project reference returned 201. The existing dataset retained only its prior reference. The import passes the reference in the incoming record origin, whereas the reusable link comes from the collection origin. `molecular_database_api.py:229`, `molecular_database.py:800`, `datasets.py:31`. | Reject project references on append with a pointer to the revision-checked project-reference route, or implement an explicit transactional link update and return the links actually saved. |
| STUDIO-07 | P2 | The dataset member UI returned an input record with zero observations. The authenticated generic REST list for the same group and input identity returned a scientific observation from another authorized group. That aggregate behavior is intentional for the generic molecule query. The gap is that REST offers no equivalent input-only member projection, so clients cannot reproduce what Datasets shows by applying the group filter. `datasets_ui.py:142`, `datasets.py:186`, `molecular_database.py:1273`, `molecular_database_api.py:322`. | Expose paginated dataset members through the same input-only service used by the UI. Keep whole-molecule aggregate search available under its current endpoint. |
| STUDIO-08 | P2 | **Create MD analysis pipeline** constructs a valid starter with three nodes in the canvas. `md_analysis` is absent from the template catalog, and the draft endpoint returns 422 for it. The tool workflow and graph node both exist. `pipeline_graph_ui.py:1025`, `pipeline_graph.py:817`, `pipeline_graph.py:2478`, `workflow_setup_api.py:104`. | Register the starter in the common template registry and check the UI and API draft paths against the same catalog. |
| STUDIO-10 | P1 | A disposable Run operation receipt identified one execution. Two jobs shared the same bare run ID, owner and project but had different task groups and workstations. Both inherited the same logical Run output. The fallback index uses the unqualified run ID, so older receipts can misattribute a result or affect the scope of downstream Run actions. No production collision was investigated. `logical_runs.py:226`. | Join receipts by cluster job ID or the verified node/task-group/run tuple. Keep unresolved legacy links visible when the producer is ambiguous. Add a cross-node and cross-task collision regression. |
| STUDIO-11 | P2 | An authenticated GET on logical Run results dispatched a remote inspection in the fixture. After a temporary inspection failure, subsequent GET requests polled the same failed command. The cache key only changes with the source timestamp. The individual execution API already distinguishes GET from an explicit inspection POST. `api.py:798`, `run_results.py:49`, `job_control_api.py:141`. | Make GET return cached state and explicit refresh requirements. Start or retry inspections through an authorized action with a bounded retry policy. Preserve the completeness check on exports. |
| STUDIO-12 | P2 | The first Results review write on a fresh database fails with `cannot start a transaction within a transaction`. The fresh connection reaches the caller with a transaction already active after schema setup. A subsequent connection has no active transaction. The failing review test reproduced on the current baseline and on `89e6d6b`. `molecular_database.py:367`, `:1934`, `browser_revision.py:13`. | Commit initialization and schema backfill in a distinct transaction before handing the connection to revision-checked writers. Test the first write on a new and an upgraded database, then confirm that stale reviews still fail safely. |

### Previously recorded peer replay defect

**STUDIO-09 · P1.** The previous release audit recorded SQLite lock failures
when applying peer project metadata alongside molecular and name records.
`Replica._apply` holds the molecular write context, then the project branch
calls `workspace.save_organization`, which opens another connection to project
the names in the same SQLite file. This can block applying otherwise valid
project/name updates. See `cluster/replication.py:374` and `:621`,
`workspace.py:126`, and `entity_names.py:72`. This is documented against
`026547e` in `docs/testing/2026-11-19-database-workspace.md`. Five naming
tests reproduced the lock in the focused peer run. The broad run also
reproduced it during shortlist review synchronization. Use the replay
connection for the SQL projection or move it to an explicit post-commit step
with retry semantics.

## Missing or uneven public API coverage

These are API parity findings from the authenticated operation inventory and
the actual UI call sites. A nearby API with a different data model is listed
where it already exists.

| Area | Existing coverage | Missing or uneven |
| --- | --- | --- |
| Projects | Cluster organization, create and rename. Streamlit uses project services for inventory and lifecycle. | Standalone project list/detail and archive/restore. Project cleanup review, progress and resume are UI operations. The peer lifecycle route uses node credentials and is excluded from public discovery. Expose read-only review and operation status first and carry through the UI's confirmation and recovery rules for any future mutation. |
| Pipelines | Catalog, compatibility, templates, draft creation, save, immutable versions, validate, run and completed-run adoption. | No single-draft read/update contract with revision checking. History reads fail for archived projects. MD starter drift. A successful local enqueue returns a saved revision while the remote worker's submission receipt returns the job reference without that revision. Establish one revision-bearing run contract. |
| Datasets | Discover, register, import, edit project references and modify group membership. | Input-scoped member list/detail with the same view as Datasets. The existing generic metadata patch and group schema operations have a different purpose from revision-checked canonical Properties edits. |
| Data browser | Molecule, evaluation, campaign and model search. Job reports and individual structure files. | No public equivalent of the browser's folder selection, exact result/model table and whole-folder CSV projection with linked analyses. The generic molecule export uses molecular IDs. Per-job reports use a job reference. Neither expresses the browser's exact result selection. |
| Results and execution | `/jobs`, `/executions`, `/runs` and cluster commands cover local jobs, workstation executions and logical runs. Per-artifact reads, job reports and logical Run CSV/JSONL exports exist. | The authenticated scientific bundle produced by `build_run_zip` is UI-only. Logical Run aggregation and individual execution inspection disagree on GET behavior. Document when to use each reference type and which representation applies to a selected result. |
| Tool tile preparation | REST lists workflows, parameter schemas and generic validation. Many specialized UI forms assemble requests or use direct service calls. | A machine-readable tile-to-workflow/profile mapping and preparation capability labels. The separate tile allowlist, grouping rules, graph registry and setup helper registry let discoverability drift. Specialized gaps need to be verified by journey rather than inferred from missing form fields. |

## Tool and workflow usability

- **Clarify how a tool can be used.** The inventory found a dedicated BinderScope
  tool and a Data browser action, with no BinderScope workflow mapping in the
  canvas registry. MD Analysis has a canvas node and a UI starter. Vina has a
  standalone workflow, with no directly mapped Vina canvas node. Label these
  entry points at the tile and in REST discovery. A missing canvas node alone
  does not establish a broken workflow or justify making every tool a node.
- **Make saved pipelines recoverable across sessions.** Canvas autosave, Save
  version, Restore version and Run are separate operations that currently share
  draft writes. Show the saved draft revision, a durable save error, and the
  immutable revision that produced each run. Audit evidence supports consolidating
  save and readiness behavior at these specific seams.
- **Preserve the difference between inputs and results.** The Datasets UI uses
  input-only membership, including sequence/SMILES entries without coordinates.
  Data browser separately derives saved outcomes from SQL and links available
  analysis evidence. The user-facing separation is backed by tests. Explicit
  dataset links are project references, whereas project use in a calculation also
  needs retained run provenance. The comparison screen currently labels
  completed evaluations `Pass` even with zero cutoffs and explains this in a
  caption. Use a comparison-specific label such as **Meets selected cutoffs**
  and display **No cutoffs applied** when the rules are empty
  (`comparison_ui.py:22`, `comparisons.py:106`).
- **Name the scope of exports.** The authenticated UI ZIP includes retained final
  candidate evidence through a verified scientific bundle. The original
  parent-only ZIP finding was partially addressed. The user must distinguish
  final scientific export from an archive of every raw intermediate and child
  run. Keep the verified serving boundary when adding API parity.
- **Avoid a broad rewrite.** `run_app.py` still contains the retired legacy
  Database rendering and compatibility handoffs, while current navigation calls
  the Data browser and Datasets modules. Retain verified compatibility routes
  and move only the active read, validation and presentation contracts needed
  by both UI and API into reusable services. Size alone is insufficient
  justification for removing a module.

## Verification and limits

- The detached checkout was inspected at the deployed commit. The main UI and
  API reported that same release directory. Audit requests used authenticated
  in-process HTTP against fixture accounts. Streamlit regressions exercised page
  state, controls and SQL handoffs. The signed-in production browser and
  each worker's installed engine inventory were outside this audit.
- The test launcher mounted the host filesystem read-only outside the audit
  directory, hid `/mnt` and `/run`, unshared networking and the process
  namespace, and gave each suite its own Studio, Results, reference and cluster
  roots. GPU inference, model downloads and opt-in engine containers were outside
  this run. Real scientific accuracy, GPU readiness, and mixed-version or live
  fleet behavior require their own acceptance evidence.
- Reproduction uses disposable HTTP and SQLite fixtures. Outputs and exact test
  totals are listed in the evidence ledger below. A probe that passes while
  asserting a defect confirms its observed behavior.
- The workflow overlap fixture spends time in synchronous result collection
  after its mocked execution stages finish. A timed stack sample showed
  PDBFixer/OpenMM called from candidate serialization
  (`graph_pipeline.py:2197`, `core/candidates.py:823`).
  This demonstrates a test-time cost and helps locate delays in finalization.
  Production timing requires representative workload measurements.

### Evidence ledger

| Check | Outcome | Local evidence |
| --- | --- | --- |
| Authenticated inventory and defect reproduction | 10 passed. The probes assert the observed behavior, including the defects. | `final.log`, `final.xml`, `test_functional_probes.py` |
| Datasets, Data browser, comparison, project lifecycle and Run results | 183 passed | `ui-data.log`, `ui-data.xml` |
| Pipeline control, scheduler, cluster, catalog, MD and result presentation | 284 passed | `operational.log`, `operational.xml` |
| Focused end-to-end fixture regression across the browser, canvas, native execution contracts, REST and the prior audit corrections | 523 passed, 4 warnings | `focused.log`, `focused.xml` |
| Peer naming and molecular replication | 18 passed, 5 failed. The five naming replay failures share the SQLite lock in STUDIO-09. | `peer-check.log`, `peer-check.xml` |
| Broad committed suite, interrupted at about 97% during slow end-of-suite workflow fixtures | 3,832 passed, 12 failed, 64 skipped, 103 warnings before interruption | `suite.log`, `suite.xml` |
| Broad-failure triage on the current and previous release | The three temporary-file failures passed in a sandbox with its own writable `/tmp`. The two Complex Analysis assertions reproduced on both releases, and the stale sidebar-order assertion reproduced on the current release. | `triage.log`, `prior-triage.log` |

The 12 broad failures break down into six peer replay cases (STUDIO-09), one
first-write review failure (STUDIO-12), three sandbox temporary-file failures,
and two test expectations. The sidebar test expects the retired **Database**
and **Complex Analysis** entries. The reversed-chain ingestion test
successfully attaches the evaluation to the exact structure, then indexes the
omitted `target_id` field directly. The ingestion code deliberately drops the
old target identity when the selected target chain scope changes. Update the
test to assert that missing target identity explicitly. These two expectations
are separate from the application findings.

The evidence directory is `/tmp/mds-functional-audit.VrC6Hu`. The launcher is
`isolated.sh` and the detached source is in its `source` worktree. To rerun
the authenticated probes on this installation, use:

```bash
audit=/tmp/mds-functional-audit.VrC6Hu
bash "$audit/isolated.sh" replay -m pytest -q -s "$audit/test_functional_probes.py"
```

## Repair order

1. Fix stale draft writes, logical Run receipt identity, the mixed-provenance
   ranking gate and the peer transaction boundary. Retest concurrent editor
   saves, cross-workstation ID collisions and mixed replay.
2. Fix first-write review initialization. Align local and remote validation, separate definite submission rejections
   from uncertain dispatch, and fix archived pipeline access and MD template
   discovery. Give failed remote inspections an explicit recovery path. Verify
   each through normal UI and HTTP entry points.
3. Add input-scoped dataset reads, project lifecycle visibility and exact
   saved-result/export parity. Provide one reference model per read surface
   and test ACL, history, and revision behavior before adding write operations.
4. Maintain a small acceptance matrix per tile and recipe covering discovery,
   preparation, validation, enqueue/recovery, result registration, inspection,
   and export. Record which workstation and runtime were exercised separately
   from the offline contract results.
