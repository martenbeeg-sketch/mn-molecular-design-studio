# P0 closure review and dispatch-race correction

Date: 2026-11-05

**Verdict: APPROVED WITH LIMITATIONS for the corrected working tree. P0 can be closed at the code-review level. Deployment remains pending.** No blocking scientific-selection defect remains in the inspected and tested paths. The limitations are live-workstation/version verification and production computation/transport, which were deliberately not exercised. Commit `1d4c8fb` alone does not include the additional dispatch-race fix described below.

## Reviewed versions and working-tree changes

HEAD: `1d4c8fb7b74b0494700a60d200dc452a9c2b3417` — `Close P0 F4 and F5 findings`.

Earlier implementation commits:

- `82ef6fd591443071a47213faafa320c7dec013bc` — `Fix P0 database findings F1-F5`.
- `cb8a4ae3d734734183e320a9266cc3fcb0de70f3` — `Default screening curation to explicit selection`.

The closure review found one additional F5 dispatch race. After the user explicitly changed the review-only instruction to authorize a fix, local changes were made to:

- `mn_molecular_design_studio/engine/workflow_engine/workflows/input_batch.py`: accept validated, authenticated submission tags when the parent Studio marker has not yet been written.
- `tests/test_api_integrations.py`: add 14 real-worker timing/policy/retry cases and seven public-parameter authorization cases.
- This new closure report.

No commit, staging or push was performed. Pre-existing untracked audit documents were left untouched: the three directories under `docs/audits/2026-11-05-database-{phase2a-decisions,phase2a,simplification}/`, and the four `docs/database-phase1-{closure-verification,corrective-review,final-independent-verification,independent-review}.md` files.

## Additional F5 failure and bounded correction

The API persists `studio_owner_id` and `studio_project_id` in job metadata before dispatch, but writes `molecular_design_studio.json` after dispatch returns. Previously, `run_input_batch` read only the latter file and cached absent lineage. A worker starting during this interval therefore created children without authorized Studio markers. The API could mark the parent before either child completed, yet both explicitly requested Vina imports were silently skipped; retry skipped the completed children and still produced zero Evaluations.

The correction reads only the existing server-owned metadata tags when the parent marker is absent. It requires a string project tag and the existing nonempty-string owner validation. Generic ownership fields, item parameters, nested lineage objects and `input.json` cannot authorize a child. Existing valid marker behavior and all per-item import policies remain unchanged. The non-coalesced parent stays manual to prevent aggregate duplication. No replication adapter, protocol, schema, ACL, tombstone, historical observation or source identity code was changed by this follow-up.

The new regression starts `local_worker.run_worker_job` in a fresh thread before the API writes its marker, holds the first Docker substitute with synchronization events, and releases it after the API has marked the parent. This models separate-worker ContextVars and deterministically recreates the actual ordering. Only Docker output is substituted; authenticated HTTP submission, validation, queue persistence/restoration, batch/Vina adapters, staging and completion ingestion execute normally.

Restoring only HEAD's original lineage helper in an isolated pytest plugin makes the permanent true/true separate-child case fail with **0 observations versus 2 requested**. The corrected helper passes that same test. The control did not alter repository files.

## Verification results

Interpreter: `/home/user/mambaforge/envs/mn-molecular-design-studio/bin/python`.

| Verification | Result |
|---|---:|
| Same 25-module combined regression selection documented in the F1–F5 correction report, including new permanent cases | **607 passed in 169.85 seconds** |
| Original independent scenarios, legacy completion probes, previous final review probes and expanded closure probes | **97 passed in 51.71 seconds** |
| Final focused permanent API timing/authorization tests | **21 passed in 25.21 seconds** |
| Permanent regression with original HEAD lineage helper restored | **1 expected failure: 0 versus 2 Evaluations** |
| `git diff --check` | Pass |

Evidence: [combined log](/tmp/mds-p0-closure-review.5cZ0XR/fixed-combined.log), [combined XML](/tmp/mds-p0-closure-review.5cZ0XR/fixed-combined.xml), [independent log](/tmp/mds-p0-closure-review.5cZ0XR/fixed-independent.log), [independent probes](/tmp/mds-p0-closure-review.5cZ0XR/test_closure.py), [focused log](/tmp/mds-p0-closure-review.5cZ0XR/fixed-api.log), [baseline failure](/tmp/mds-p0-closure-review.5cZ0XR/baseline-permanent.log).

### F4

The original independent same-group/canonical-molecule collision and foreign-authoritative replica collision pass. The current duplicate check binds the expected original observation ID and exact historical source key, group, owner, namespace, workstation, run and compound. Replica reuse additionally requires the authoritative node, derived replica IDs and original group/observation/source markers. A matching unrelated source row cannot substitute for the selected observation.

Extended probes cover ten historical formats/provenance cases and seven replica-marker cases. Missing or inconsistent identity evidence does not suppress a valid selected local Evaluation. Exact historical same-source replay avoids cloning its observation; selected-source replay stays idempotent. Distinct local/foreign sources, foreign private observations and ordinary campaigns outside the audited shortlists are covered. Original scientific rows are retained, including measurement IDs and values; foreign-private access remains denied.

### F5

Both normal restored-queue execution and dispatch-before-marker execution were tested with compatible and differing Vina settings:

| Per-item choices | Expected and observed imported Evaluations |
|---|---:|
| omitted / omitted | 0 |
| false / false | 0 |
| true / true, coalesced or separate | 2 |
| true / false; false / true | 1, belonging to the true item |
| true / omitted; omitted / true | 1, belonging to the true item |

Assertions verify distinct exact scores, CSV/pose bytes, pose digests and structure references, observation counts, reuse of an intentionally imported canonical ID, real authenticated owner and nonempty project, settings, compound/run/source attribution, parent-completion stability and full serialized-worker retry. Coalesced retry reruns its scalar adapter deterministically; completed separate children resume without rerunning Docker. Neither creates duplicate Evaluations. Seven forged public parameter fields return HTTP 422 before dispatch.

### F1–F3 and preservation

- Omitted HTTP/MCP selections create no imported scientific rows; exact IDs and explicit positive top-N selections validate and replay correctly. The independent MCP cases use the actual SDK tool callable, StudioClient and isolated ASGI API.
- Stable-stop and iteration-limit completion retain saved positive legacy limits and Results while publishing zero imports. No additional positive automatic final-publication entry was found in the inspected callers.
- Both exact SpiroChem shortlists remain excluded on source and real-adapter replicas. Scientific records, direct lookup, independent curated memberships, private ACLs and tombstones remain available/effective. Development filtering remains restricted to the eight audited owner/run pairs.
- Intentional import, binder, historical Results/reference, structure, accepted-design and owner-qualified replication regressions pass. Preservation fixtures include an intentional compound import and genuine binder membership alongside excluded shortlists.
- Ordinary Results/Database navigation does not initiate scientific ingestion in the exercised paths. Existing derived Evaluation/metric cache rebuilding may write derived projections; repeated reads remain stable.

## Precisely bounded limitations

All scientific writes were confined to temporary fixture databases and run directories. The real replication adapter was used with isolated source/replica state; no live peer synchronization or workstation deployment parity was tested. Docker inference, live MCP transport and a new full inspection of HY-L126's production supplier dataset were not performed. The preservation conclusion combines unchanged relevant production code with regression/fixture evidence; it is not an inventory or checksum verification of all production scientific records.

No production data, service, deployed process or peer state was changed. No deployment, restart, synchronization, reconstruction or Phase 2B work occurred.

## Smallest safe deployment proposal — not executed

1. Commit the reviewed working tree, including the dispatch correction and tests, and record that exact SHA. Deploying only `1d4c8fb` would omit the reproduced race correction.
2. Arrange a quiet job window on **pc002975, PC004032 and PC004033**. Take consistent, restorable backups/snapshots of each workstation's scientific database and associated Results/source files, application configuration and replication state, including tombstones. Include SQLite WAL state through a supported consistent backup method. Record current code versions for rollback.
3. Install the same reviewed SHA and compatible application environment on all three workstations through the normal approved deployment process. Verify both checkout SHA and the code version actually loaded by API/workers. This change requires no reconstruction, scientific migration, reimport, campaign replay or forced peer synchronization.
4. Check read-only Library/Results/direct-history access, both excluded shortlist groups, independent intentional imports/HY-L126, existing binders, ACLs and tombstones. Use isolated canaries for omitted and explicit mixed batch imports, including retry. Compare scientific row counts before/after ordinary navigation and confirm no unexpected ingestion.
5. If checks fail, stop new submissions and restore the prior code version. Preserve any newly produced scientific outputs and database changes. Use backups only through a deliberate recovery plan that accounts for post-backup data; do not overwrite new scientific records blindly.

P0 implementation review can close for the corrected code; deployment acceptance requires the recorded SHA, backups and workstation checks above.
