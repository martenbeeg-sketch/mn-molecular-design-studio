# Database simplification P0 review corrections F1–F5

**Date:** 5 November 2026

**Base reviewed:** `cb8a4ae3d734734183e320a9266cc3fcb0de70f3`
**Scope:** Correct the five findings from the independent P0 review. No Phase 2B, P1, deployment, synchronization or production-data work.

## Verification before changes

I read the independent `review.md` and the existing P0 implementation report, then ran the review's temporary probes against the app's Conda environment:

```text
/home/user/mambaforge/envs/mn-molecular-design-studio/bin/python -m pytest -q \
  /tmp/mds-p0-independent-review.BLno9T/test_review_probes.py \
  /tmp/mds-p0-independent-review.BLno9T/test_legacy_graph_probe.py \
  --basetemp=/tmp/mds-p0-independent-review.BLno9T/verified-state-ours --tb=short
```

**Before correction: 8 failed, 11 passed.** The failures independently reproduced all five findings:

| Finding | Before-correction evidence |
|---|---|
| F1 | An empty API body returned HTTP 200 and imported the requested default 1,000 (8 memberships and 11 sources in the fixture); exact-ID input was rejected as an extra field. |
| F2 | A saved graph limit of 1,000 caused the stable-stop case to write 10 sources/observations and the iteration-limit case to write 6, without a current Results-page selection. |
| F3 | Both audited original groups were hidden from Library, while both real-adapter `peer-*` replicas became visible there. |
| F4 | An inaccessible private observation for the same molecule/run/compound suppressed the local selected Evaluation; the fixture had zero local observations. |
| F5 | With the docking import choice omitted, children used their Results-only default while the queued batch parent retained `automatic` policy. |

## Corrections

### F1 — Public API selection boundary

**Root cause.** `FinalImport.limit` defaulted to 1,000, and the HTTP contract did not accept exact compound IDs. A valid request with no curation choice therefore authorized an import.

**Changed:** `integration_api.py`, `mcp_server.py`, `bulk_screening_results.py`; regressions in `test_api_integrations.py`, `test_mcp_server.py`, `test_bulk_screening_collection.py` and `test_bound_docking.py`.

**Correction.** The API and MCP tool now default to zero. Both accept an optional nonempty exact-ID list; an empty list, unknown ID, or simultaneous positive limit and IDs is rejected before scientific writes. A positive ranked top-N remains available when the caller explicitly supplies it. The publisher also rejects a positive limit unless its caller marks an explicit curation action. Exact selections reuse canonical IDs and stable group/source keys, so replay is idempotent.

**Evidence.** The permanent HTTP tests prove `{}` creates no source or Evaluation and leaves the Parquet Results intact; exact IDs import only the requested records and replay without adding rows. Empty, unknown and conflicting selections return HTTP 422. The MCP schema test confirms default limit `0` and optional exact IDs. The pre-fix API probes had demonstrated the default import and rejected exact-ID contract described above.

**Remaining limitation.** None observed within the public HTTP/MCP and trusted publisher paths covered here.

### F2 — Historical graph limits

**Root cause.** Ordinary campaign completion read `final_database_import_limit` from the saved graph and passed that historical value directly to the real publisher.

**Changed:** `engine/workflow_engine/workflows/active_learning_controller.py`, `pipeline_graph.py`; regressions in `test_active_learning_stopping.py`.

**Correction.** Completion always publishes with limit zero. It preserves the graph snapshot's saved value and records that value as `saved_import_limit_preserved_but_not_applied` in the Results manifest. The graph field's label and help text describe its legacy compatibility role and direct users to make a separate Results selection.

**Evidence.** A real `publish_final_selection` call against a temporary `MolecularDatabase` is exercised for stable-stop and iteration-limit completion with both saved values `0` and `1000`. All four cases keep the graph value, retain the complete Parquet row count and create no molecule, membership, source or observation rows. Before the fix, the same real-publisher probes wrote 10 and 6 observations for the two positive-limit cases.

**Remaining limitation.** None observed for ordinary completion and the two tested final-publication paths. Explicit Results-page imports remain available after completion.

### F3 — Replicated shortlist visibility

**Root cause.** The audited eligibility rule matched only the original collection primary key. The existing peer adapter creates a `peer-*` key and stores the original group ID and source workstation in replica metadata.

**Changed:** `molecular_database_visibility.py`; real-adapter regressions in `test_database_simplification_p0.py`.

**Correction.** The filter resolves a peer group to `origin._replica_group` only when `_replica_node` is present, then applies the same exact audited owner and run checks. It filters collection membership only; it does not blacklist canonical molecules or alter replication state.

**Evidence.** Parameterized tests cover both audited source IDs using the real `apply_record` adapter. Source and replica groups are both excluded from Library; an independent local membership restores the shared molecule to Library without making the shortlist group eligible. All indexed records, the run-scoped Results reader, direct molecule/group history and private-group ACL checks retain their expected behavior. Before correction, both replicas were visible in Library. Snapshots confirm the source/replica scientific rows and replication `records`, `state` (including a tombstone) and `conflicts` are unchanged by visibility reads.

**Remaining limitation.** None observed for the two audited owner/run pairs or replica metadata produced by the existing adapter.

### F4 — Evaluation deduplication identity

**Root cause.** The duplicate check matched molecule, run basename, compound ID and namespace globally. It did not require the same owning group, source workstation or source provenance, so unrelated private or peer observations could suppress a local Evaluation.

**Changed:** `bulk_screening_results.py`; regressions in `test_bulk_screening_collection.py`.

**Correction.** Duplicate suppression now joins the exact observation to its source and group and requires the campaign owner, original audited final-group identity (or its replica metadata), source workstation, run/compound context, and matching final source key/provenance. Foreign, private and other-workstation sources remain distinct. Historical automatic observations are not rewritten.

**Evidence.** Tests prove a matching historical local automatic Evaluation is retained unchanged and not cloned; replaying the same selected source is idempotent. A foreign private observation and a same-run observation from another workstation each allow the local selected Evaluation to be created. The local user cannot read the private foreign group or see its score, and the other workstation's score and authorship remain associated with its own replica. The pre-fix private-source probe produced no local observation.

**Remaining limitation.** None observed for the owner/source/run/compound provenance emitted by the current publisher and replication adapter.

### F5 — Batch docking parent policy

**Root cause.** Batch submission marked the parent external-only only if an item explicitly contained `import_to_database: false`. Omitted values were not resolved from the workflow function's default before the parent was queued.

**Changed:** `engine/workflow_engine/services/workflow_automation.py`; request/queue/aggregate regression in `test_bulk_screening_collection.py`.

**Correction.** Before parent creation, submission reads the workflow parameter's declared default, resolves each child's effective choice and validates that each choice is boolean. A Results-only child makes the parent policy manual. The existing input-batch and completion hooks then inherit that policy.

**Evidence.** A real workflow request creates a two-item batch with the choice omitted, executes the batch completion path with deterministic fixture children, and asserts manual policy for parent and both children. The ingest hook is never called and all scientific tables in the isolated database remain empty. The pre-fix probe observed `automatic` on the parent.

**Remaining limitation.** None observed for the omitted-choice multi-item docking path exercised here.

## Full regression and state boundary

The previous P0 regression selections plus the new tests were run together in the app environment:

```text
/home/user/mambaforge/envs/mn-molecular-design-studio/bin/python -m pytest -q \
  tests/test_database_simplification_p0.py tests/test_bulk_screening_collection.py \
  tests/test_molecular_database_results_ui.py tests/test_active_learning_graph.py \
  tests/test_active_learning_stopping.py tests/test_api_integrations.py \
  tests/test_ligand_tools.py tests/test_molecule_imports.py \
  tests/test_historical_results_import.py tests/test_accepted_design_import.py \
  tests/test_docking_results.py tests/test_molecular_database_peer.py \
  tests/test_molecular_database_policy.py tests/test_molecular_database.py \
  tests/test_molecular_database_integration.py tests/test_screening_convergence.py \
  tests/test_peer_workstations.py tests/test_candidate_peer_transfer.py \
  tests/test_mcp_server.py tests/test_bound_docking.py \
  tests/test_native_collection_execution.py tests/test_workflow_queue_lifecycle.py \
  tests/test_complex_analysis_dashboard.py tests/test_complex_analysis_collection.py \
  tests/test_protein_library.py
```

**Result: 569 passed in 145.15 seconds.** `git diff --check` also passed.

All new database and replication checks use temporary fixture databases, temporary application/Results paths and an isolated replication adapter. No production Database/Results files, live peer state, workdirs or service were changed. No peer synchronization, deployment, service restart, push or reconstruction was performed. Existing untracked audit documents were left untouched. The correction is ready for the requested independent review; that separate review is not represented as having been run here.
