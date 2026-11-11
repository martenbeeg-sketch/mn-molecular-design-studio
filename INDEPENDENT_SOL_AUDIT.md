# Independent SOL audit: Molecular Design Studio after Phase 2

**Audit date:** 2026-11-11, Europe/Berlin.  
**Committed baseline:** `e7cd3a95c5163f6e973187a24aab002f5e5c351b` (`e7cd3a9`).  
**Checkout:** detached, separate managed worktree at `/home/user/.codex/worktrees/independent-sol-audit/mn-molecular-design-studio`.  
**Deliverable scope:** technical assessment and audit evidence; application source, committed tests, authentication configuration and the active development branch are unchanged.

## Assessment

The two phases make substantial, independently reproducible improvements. New automatic observation attachment is conservative about scientific identity; the central Database structure reader verifies retained bytes and returns a private snapshot; strict finalization rejects undeclared dangling scientific outputs; updated review writers and editors enforce optimistic concurrency; the batch job lookup finds older exact-source work without scanning Results during coverage rendering; and Tag/Fusion accepts the supported multi-chain partner assembly.

The application is **not yet verified as consistently correct across its scientific workflows**. The most consequential remaining problems are version checks bypassed outside the batch planner, target-revision-blind comparison, a receptor-viewer integrity bypass, and a confirmed Phase-2 regression that makes a newly completed native interface-chemistry Evaluation appear incompatible with the identical request. Scientific parameter identities also differ between assessment planning and comparison. Concurrent API sessions still create duplicate scientific work.

The focused verification run passes **191 tests**. Additional adversarial probes reproduce the defects below; passing the implementation regression suites does not establish the untested cross-module contracts. The broader committed-suite result is recorded in the testing section and is not represented by the earlier reported 663-test subset.

## Method, independence and operational record

I first inspected the committed architecture, runtime boundaries, scientific descriptors, source adapters, job projection, review transaction and peer application code **before reading `AUDIT.md`**. I then compared those observations with the original findings and both implementation ledgers. Commit history was used to distinguish residual defects from regressions. The implementation ledgers were treated as claims to verify, not as independent evidence of correctness.

The worktree was explicitly created at the full recorded SHA; no new branch, commit, merge, deployment, service restart or authentication change was performed. Existing untracked development documents in the primary checkout were left alone. Tests used disposable SQLite databases, structures, private inputs, queue records and simulated peers. Native interface chemistry was calculated only on a tiny disposable two-chain fixture. GPU and installed-container tests were not enabled.

### Isolation incident and immediate correction

The initial environment-only test isolation was insufficient. API import/request middleware calls `ensure_studio_workspace()`, which reloads Settings and can replace `RUN_DIR` with the default shared Results directory. A custom API probe consequently created two **unexecuted audit queue directories** in production Results: `analysis/20261111-172829-e9cd2cd2` and `analysis/20261111-172856-874abe10`. Their metadata remained queued, their logs were empty, and neither contained a result. This breached the requested no-production-write constraint; it must not be described as a wholly isolated initial run.

I immediately interrupted the initial broad run and moved only these two audit-created directories to `/tmp/mds-independent-sol-20261111/escaped-results/`. A subsequent read-only check found neither ID in the production job catalog or molecular-ingestion table and found no remaining Results directory created since the initial run began. No pre-existing production scientific file was edited as part of that correction. Production services and authentication were not changed.

All subsequent test runs used an OS mount/network sandbox: the host filesystem is read-only except `/tmp` and this audit worktree; `/mnt` and `/run` are hidden; networking is disabled; and the test runtime has isolated Studio, Results, reference and cluster roots. A direct sandbox check confirmed production Results was invisible. Custom fixture Settings also explicitly retain the isolated Results/reference paths. The sandbox launcher is preserved with the evidence. The interrupted initial broad run is excluded from completed-suite totals.

## Independent verification of the eight implemented findings

| Original ID | Independent status | Evidence and remaining boundary |
|---|---|---|
| **AUD-001** | **Verified for new automatic attachments** | `molecular_database.py:34`, `:61`, `:787`: all supplied candidate/model/pose/digest selectors must match a unique representation in the source group/run. The committed identity suite covers insertion order, ambiguity, nested selectors, source versus calculation digests, zero values and peer-qualified source groups. Historical bad links are deliberately not repaired. Numeric-zero identity is still lost later in the review projection: SOL-007. |
| **AUD-002** | **Verified at the shared Database structure reader; incomplete across presentation routes** | `molecular_database.py:1540`, `:1564`: retained digest verification precedes a private, atomic, read-only snapshot. Changed sources, corrupt replica/cache data, downloads, materialization and analysis are covered. My additional source-change probe confirms the old snapshot remains exact and the next read fails. The legacy bound-receptor viewer bypasses its existing frame verifier; raw Results downloads also do not enforce recorded digests: SOL-006. |
| **AUD-003** | **Verified for the stated finalization contract** | `engine/workflow_engine/core/artifact_collection.py:319` and `core/jobs.py:606`, `:649`: scientific outputs and exact retained aliases take precedence over temporary-origin heuristics; both successful-completion entry points use the validator. All 31 committed artifact tests pass in the focused run. Wholly absent, undeclared, unreferenced regular files remain outside what the validator can infer. |
| **AUD-004** | **Partially resolved; one confirmed regression** | The batch coverage helper now compares protocols and complete scientific parameters, prefers successes deterministically and retains exact source/model/chain scope. Its version/permutation tests pass. Inspector/Results matching still bypasses the protocol contract (SOL-001); a real native chemistry result fails reuse against its own request (SOL-002); comparison uses a different parameter signature (SOL-004). |
| **AUD-005** | **Verified for updated writers/editors, with the documented mixed-version limit** | `molecular_database.py:1828`, `:1851`: content/generation are read together, `BEGIN IMMEDIATE` makes batch compare-and-set atomic, and insert/update/delete generations survive restart. The 18 committed concurrency cases and the simulated-peer review convergence/conflict tests pass. My deletion/restart probe also rejects both stale pre-creation and pre-deletion revisions. An old binary can still perform a blind SQL write; triggers invalidate new forms but cannot prevent that old writer. |
| **AUD-006** | **Lifecycle classification corrected; admission remains partial** | `complex_analysis_batch.py:47`: collecting, transferring, finalizing, paused, unknown and future nonterminal states block batch submission; completion blocks pending ingestion. The UI rechecks the index before callbacks. There is no atomic server-side reservation of scientific work: two independent request keys still produce two jobs (SOL-005). Actual SQLite failures also escape the intended unavailable-state handler (SOL-008). |
| **AUD-011** | **Indexed lookup correction verified; migration/error limitations remain** | `job_catalog.py:289`, `:351`: source filtering happens before 500-row paging, with an expression index; Database coverage uses one shared lookup and no broad Results traversal. Committed tests cover >1,000 unrelated and >1,000 same-source jobs, ACLs and read instrumentation. A separate migration probe verifies old projections block lookup, then ordinary reconciliation restores an exact collecting job. Inspector still has separate Results matchers with a first-1,000 limit; SQLite errors are not contained (SOL-001, SOL-008). |
| **AUD-012** | **Verified** | `complex_analysis_batch.py:103`, `:199`: Tag/Fusion uses the component-plus-partners profile, accepting one assessed chain and disjoint multiple partner chains. The five applicability cases and existing Tag/Fusion adapter/workflow tests pass; Rosetta/ipSAE retain their single-pair restriction and construct declarations remain explicit. |

“Verified” is scoped to the contract and exercised paths in that row. It is not a retrospective certification of production data or of an unchanged old deployment.

## Evidence-backed findings

Severity reflects scientific correctness or operational impact. SOL identifiers are local to this independent audit; original AUD identifiers remain stable.

| ID | Severity | Classification | Finding |
|---|---|---|---|
| SOL-001 | High | Residual AUD-004; additional cross-surface evidence | Inspector and Results reuse bypass the expected evaluator protocol. |
| SOL-002 | Medium | **Confirmed Phase-2 regression** | A real native chemistry result cannot be reused for an identical request. |
| SOL-003 | High | Newly discovered, predates these phases | Binding comparisons ignore target structure revisions. |
| SOL-004 | Medium | Newly discovered, predates these phases; AUD-014 evidence | Subject-specific paths and IDs split otherwise identical scientific protocols. |
| SOL-005 | Medium | Residual admission risk, acknowledged but not solved by Phase 2 | Independent/concurrent API submissions create duplicate scientific jobs. |
| SOL-006 | High for receptor interpretation; Medium for raw download integrity | Newly discovered cross-module extension of AUD-002/AUD-014 | Receptor visualization bypasses frame validation; Results downloads bypass saved digest validation. |
| SOL-007 | Medium | Newly discovered residual identity gap | Numeric model zero disappears between retained structures and review/assessment. |
| SOL-008 | Medium | Newly discovered residual error-handling defect | Real SQLite lookup failures crash the planner instead of producing the promised unavailable state. |

### SOL-001 — Protocol-aware reuse is not shared by Inspector and Results

**Locations:** `complex_analysis_ui.py:2397`, `:2454`, `:3901`, `:904`, `:1002`; `complex_analysis_batch.py:20`.

The indexed Database matcher enforces version compatibility only when `expected_protocol` is supplied. `_render_database_workflow()` receives a registered schema but calls that matcher without `schema["assessment_protocol"]`. A completed job can therefore be displayed/reused even when its recorded protocol is explicitly obsolete. This affects the ordinary individual Database Inspector route, not merely an unused helper.

**Reproduction:** The same exact-source indexed row with `method_version="obsolete"` is rejected when the matcher receives `expected_protocol={"method_version":"current"}`, but returned when called with the Inspector's actual default arguments. The final independent probe fails on `len(actual) == 1`. Source inspection confirms the omitted argument at the real caller. The Results Biotite/contextual matchers separately compare selected settings and drop/ignore runtime protocol revisions; they also retain their independent first-1,000-job searches.

There is a related **bounded runtime-identity limitation**: registered QC expectations contain the local runtime-source hash and image tag, but not an immutable image digest. A saved Evaluation with an arbitrary different recorded digest is still classified compatible when those expected fields match. A versioned tag is useful provenance, but by itself cannot prove identical installed container bytes across machines. This audit did not inspect or claim an actual mismatch between production image builds. Rosetta/ipSAE's unknown digest is correctly conservative in the batch helper; it does not make the Inspector's omitted protocol check safe.

**Impact:** Which route a researcher uses determines whether historical or unknown-version evidence is treated as compatible. The Phase-2 claim is valid for the batch helper, not for the whole application.

**Recommendation:** Pass the registered protocol to every reuse decision and share the same match/status policy across individual Database and Results routes. Represent unknown runtime identity explicitly; when equivalence depends on the container, compare a verified digest or a documented immutable compatibility contract. Preserve old Evaluations as historical evidence.

### SOL-002 — Fixed native chemistry output is incompatible with its own request

**Locations:** `engine/workflow_engine/workflows/analysis.py:1477`; `complex_analysis_batch.py:13`, `:279`; `job_catalog.py:342`; `complex_analysis_collection_ui.py:815`; `complex_analysis_ui.py:2438`.

The UI requests `analysis.binder_interface_chemistry` with `{}` because this version has no editable scientific cutoffs. The evaluator persists its effective fixed parameters plus `structure_digest` and `model_id`. The new exact parameter comparison compares that expanded map with `{}`. The index matcher does the same. Consequently even a newly calculated, correct, current-version result has no matching request according to the planner.

**Reproduction:** Build an authorized Database source, run the actual registered native chemistry function on a tiny two-chain PDB, mark the disposable run as a Studio run and ingest it. The result reports `success=true`; ingestion creates one Evaluation on the original representation. With the same model and chain groups:

```text
requested parameters: {}
persisted scientific parameters include:
  model_id: "0"
  structure_digest: a9f286ce...
  contact_cutoff_angstrom: 4.5
  sasa_probe_angstrom: 1.4
  sasa_sphere_points: 960
coverage: Existing Evaluation · different settings
indexed matching jobs: 0
```

The exact same retained sample, passed to the committed pre-Phase-2 `evaluation_coverage()` loaded from `2680bf0^`, returns **Already evaluated · compatible settings**. This is a confirmed behavioral regression, not a synthetic missing-version fixture. The earlier helper was too permissive generally, but this particular effective-default equivalence was lost when replacing subset comparison with raw-map equality.

**Impact:** Users are told to permit a new Evaluation despite unchanged settings. Enabling that option can repeatedly calculate the identical fixed method; the individual Inspector also misses the completed indexed job. This invalidates the claimed round-trip reuse contract for this native method.

**Recommendation:** Canonicalize **effective** scientific settings at both ends, including registered fixed defaults. Keep subject/model/digest identity in the subject key rather than the protocol map. Do not restore unrestricted subset matching. Add a real submit/calculate/ingest/review/index round-trip contract for every offered assessment method; the current version tests predominantly construct already-normalized synthetic settings.

### SOL-003 — Different target revisions are silently ranked together

**Locations:** `comparisons.py:139`; callers `comparison_ui.py:66`, `molecular_database_api.py:218`; producer `evaluations.py:65`, `:104`.

The descriptor preserves `target_revision`, but `compatibility_warnings()` checks only `target_id`, metric definitions, parameter digests and container digests. Two repaired or otherwise revised target structures can retain the same molecular identity and receive different binding measurements. The documented target-revision distinction is therefore lost at the pooling gate.

**Reproduction:** Two completed rows have the same target ID, method and parameter signature, different 64-character target revisions and a numeric binding field. The compatibility helper returns `[]`; comparison assigns ranks `[1.0, 2.0]`. The UI enables ranking when that warning list is empty, and the preview API uses the same gate. Source history shows this omission predates Phase 1.

**Impact:** A seemingly comparable binding-score ranking can span different receptor structures without warning. The subject-composition key's inclusion of target revision does not protect comparisons between different subjects.

**Recommendation:** For target-dependent fields, gate pooling on target identity **and** relevant structural revision/chain scope. Treat known-versus-unknown context explicitly. Keep sequence-only measurements exempt where scientifically justified. Exercise both the API and UI with repaired receptor revisions.

### SOL-004 — Scientific signatures include subject-specific paths and selectors

**Locations:** `evaluations.py:18`, `:92`; `complex_analysis_batch.py:8`, `:13`; `complex_analysis_design_review.py:160`; `comparisons.py:151`.

The assessment matcher removes `candidate_ids`, `selected_structure_path` and `selected_model_id` before comparing scientific settings. The Evaluation parameter digest does not. Those fields are persisted by the real source adapters and analysis workers, so otherwise identical calculations on two designs produce different protocol digests solely because they refer to different designs.

**Reproduction:** Two observations have identical method, pH/contact cutoff and expected protocol, with only the selected structure path, selected model and candidate ID changed. `scientific_parameters()` returns equal maps. `evaluation_descriptor()` returns two different `parameters_digest` values. `review_measurements()` then creates **two disjoint columns** for the same measurement, each populated for only one model. The comparison compatibility gate likewise reports different parameters for pooled ranking.

**Impact:** The workbench cannot reliably compare the common per-design assessment profile. This is a demonstrated correctness/maintenance reason for AUD-014: separate scientific-signature implementations already disagree. File size itself is not evidence of a defect.

**Recommendation:** Use one pure effective-settings projection in planning, descriptor/index generation, review-series grouping and comparison. Keep the exact selected subject separately bound. Version any derived projection migration; do not silently rewrite canonical historical observation/Evaluation identities.

### SOL-005 — Final UI rechecks do not make admission atomic

**Locations:** `complex_analysis_collection_ui.py:888`; `api.py:672`, `:729`; `api_idempotency.py:31`; `job_catalog.py:351`; cluster request-key uniqueness in `cluster/store.py`.

The batch UI checks a catalog snapshot, then independently calls submission. The API does not reserve a unique scientific input/protocol/settings identity. Request-level idempotency only protects the same owner's same request key, and keys from different sessions are independent. Each workstation's indexed coverage query is local as well.

**Reproduction:** Two identical serial authenticated TestClient submissions, with actual source authorization/staging and real queue creation but worker dispatch disabled, each return HTTP 202 and different job IDs. A separate controlled concurrent test uses two independent request keys and a barrier immediately before admission: both return 202 and create distinct queue jobs. For that race test the already-authorized source bundle is prepared once and reused to keep staging races out of the admission experiment; the real admission/idempotency/queue code still runs. A positive control using the **same** idempotency key returns the same queue job twice.

**Impact:** Simultaneous researchers or separate peer instances can create duplicate computation even after both UI checks pass. The tests prove the local concurrent case. The distributed risk follows from local-only lookup and per-request-key uniqueness; it was not exercised on live workstations.

**Recommendation:** Add atomic admission/reservation using a normalized scientific fingerprint in the existing submission/scheduler infrastructure, qualified by canonical owner/source identities and access scope. Preserve explicit reassessment/retry semantics. Keep request idempotency for transport retries, but do not present it as scientific deduplication. The ledger already acknowledges this limitation; that acknowledgement is accurate and remains operationally material.

### SOL-006 — Presentation integrity differs from analysis integrity

**Locations:** `bound_docking.py:70`, `:191`; `presentation_api.py:112`; `api.py:1366`; `molecular_database.py:1540`.

**Bound receptor viewer:** `_resolve_pose_receptor()` verifies the saved pose's receptor digest/coordinate frame. The viewer endpoint instead calls `receptor_file()`, whose legacy `target_path` branch checks existence only. Thus the source can be rejected for an MD/analysis preparation and still be displayed as the pose's receptor.

**Reproduction:** Create an accessible disposable ligand pose with a recorded receptor digest and legacy target path, then change that receptor file. `_resolve_pose_receptor()` correctly raises a coordinate-frame mismatch. `GET /api/v1/molecules/{id}/viewer/receptor?structure_id=...` returns HTTP **200** and the changed bytes for the same pose. This uses the real presentation router and Database authorization. The shared `structure_file()` fix is not on this legacy branch.

**Raw Results download:** An additional fixture has a required `artifacts.json` structure digest. After changing the structure, `validate_run_artifacts()` rejects it, while the authorized job-artifact endpoint still returns HTTP **200** and the changed bytes. That endpoint enforces path containment and job authorization but does not enforce retained scientific identity. General raw-file access and verified scientific retrieval currently provide different guarantees without a common integrity gate.

**Impact:** A docking pose may be interpreted against the wrong receptor coordinates. Recorded finalization integrity also does not imply subsequent Results download integrity. The main Database structure download remains correctly verified; these findings do not invalidate its positive tests.

**Recommendation:** Reuse the existing frame-validating receptor resolver for visualization and serve verified snapshots where a saved scientific identity is asserted. For retained Results structures, validate available recorded digests before inspection/analysis/download, while handling genuinely unversioned historical artifacts explicitly. Do not derive a replacement identity from changed bytes and silently attach old scores to it. Add router-level tests comparing what the viewer and analysis accept.

### SOL-007 — Numeric model zero still vanishes in review

**Locations:** `complex_analysis_design_review.py:80`; `complex_analysis_collection_ui.py:779`, `:815`; related fallbacks in `complex_analysis_sources.py:266`, `complex_analysis_shortlist_export.py:179`.

Phase 1 correctly canonicalizes zero on **new observations**, but the retained-structure review projection still uses truthiness fallbacks. A valid stored `model_id=0` or `model_index=0` becomes `""`. String `"0"` works, so the issue is easy to miss with normalized fixtures.

**Reproduction:** Retain a valid single-model PDB, store numeric zero in its structure context, then query `design_review_page()`. The selected model is `''`, not `'0'`. The batch planner then regards exact structure/model identity as incomplete. No coordinate corruption is required.

**Impact:** Legacy or peer-origin zero-valued models can lose assessment eligibility or clear model provenance at presentation boundaries even though new observations retain zero correctly.

**Recommendation:** Use explicit `None`/empty checks in all model-selection projections. Test one zero-valued retained model through review, Inspector, batch submission, export and peer remapping. Keep unknown model values unknown; do not backfill historical labels by guesswork.

### SOL-008 — SQLite errors escape the new coverage failure boundary

**Locations:** `job_catalog.py:57`, `:351`; `complex_analysis_collection_ui.py:749`, `:889`; `complex_analysis_ui.py:3900`.

The UI catches `OSError`, `RuntimeError`, `TypeError` and `ValueError`, but SQLite raises `sqlite3.Error` subclasses for locks, malformed databases and some migration failures. The catalog layer propagates those exceptions. The committed failure tests simulate `RuntimeError`, so they do not test the exception the storage layer actually emits.

**Reproduction:** A real Streamlit AppTest of `_render_batch_assessment()` injects `sqlite3.OperationalError("database is locked")` at the catalog boundary. The application has an uncaught exception and no usable unavailable/refresh state. The same omission exists at the fresh pre-submit lookup and individual Inspector caller.

**Impact:** This is fail-closed with respect to dispatch, but it can crash the assessment surface during contention or upgrade instead of preserving the intended recovery UI. It violates the claimed post-Phase-2 failure contract. Inspection of the preceding Results-search caller shows the same incomplete exception family, so this audit does **not** classify its origin as a newly introduced Phase-2 regression.

**Recommendation:** Translate storage failures into a typed catalog-unavailable error at the catalog boundary, preserving diagnostics; keep submission disabled and expose retry/refresh without treating failure as absence. Test actual SQLite exception classes and a disposable locked/malformed catalog rather than only synthetic generic exceptions.

## Remaining original findings: AUD-007–010

All four remain independently reproducible. Their recommendations remain valid in the post-Phase-2 architecture.

| ID | Reproduction at this baseline | Recommendation assessment |
|---|---|---|
| **AUD-007** | Import an accessible molecule with a retained structure, create a curated group, and use the public `add_to_group()`. The new group has one molecule membership but `design_review_page()` returns **0 models**, while the original source group returns **1**. The query still filters structure ownership by the selected group. | Resolve saved representation/Evaluation membership rather than copying coordinates. For molecule-only membership, offer accessible compatible models for explicit selection. Revalidate source ACLs and honor frozen comparison IDs. The new job index does not solve collection membership. |
| **AUD-008** | `remember_review_selection()` intersects the existing selection with visible-page IDs. A selected page-one model becomes `[]` on page two and its hydrated cache is removed. Source inspection confirms acceptance/search/target filtering occurs after the 250-model page query. | Push cohort/search constraints before pagination and retain stable selected IDs independently of the visible page; revalidate access/cohort when acting. Bound page controls. This remains necessary despite the improved bounded target-scope query. |
| **AUD-009** | Metadata-only review contains explicit B/A roles and saved component sequences `KLMNPQ`/`ACDEFG`; `build_shortlist_record()` exports empty binder/target sequence maps. The light record also lacks the full source lineage consumed by export. Opening a hydrated path changes what is available. | A minimal shared metadata projection must include saved sequence/provenance fields, using explicit chain correspondence rather than guessing roles. Coordinate hydration is unnecessary. Test actual Review → shortlist → JSON/CSV output before opening 3D. |
| **AUD-010** | Six harmless formula-shaped rationale strings (`=1+1`, `+1+1`, `-1+1`, `@SUM(1,1)`, tab-prefix and CR-prefix variants) survive as exact CSV cells. No spreadsheet formula or external command was executed. | Apply a documented spreadsheet-literal policy to untrusted string cells, including names/IDs as well as rationale. Preserve real numeric negatives and exact JSON text. CSV quoting alone is insufficient. |

These are previously known unresolved defects, not regressions attributed to the two completed phases.

## Operational assessment: AUD-013

Read-only service inspection confirmed both local primary services are active from deployment directory `/home/user/programs/deployments/mn-molecular-design-studio/p11-25e003c/source`, **not this audit's commit/worktree**. Both processes explicitly use `MN_MOLECULAR_DESIGN_STUDIO_AUTH_MODE=development`. An unauthenticated read of `GET http://127.0.0.1:8766/api/v1/auth/me` returned HTTP 200, `username=developer`, `role=admin`, `development_mode=true`. Socket inspection confirmed the API listener at **127.0.0.1:8766**. This does not establish public network exposure. The other workstations and their bindings/auth modes were not rechecked.

The operational risk remains: researchers who reach this instance share a synthetic admin identity, bypassing per-user/project checks and coalescing human authorship. The source defaults to production authentication; this is an explicit deployment mode. The CAS correction prevents stale writes among updated clients but does not turn a shared synthetic account into independent author identities.

**Recommendation:** Plan the account/ownership transition and rollback on staging before any separately authorized deployment/authentication change. Verify unauthenticated rejection, account/password-change flows, owner/member/outsider/admin permissions, peer identity propagation and historical developer-owned access. Preserve localhost isolation meanwhile. No authentication change was attempted by this audit.

## AUD-014: only demonstrated coupling warrants action

The metadata review module imports a private candidate projection from a Streamlit module (`complex_analysis_design_review.py:41`), while signatures and model defaults are reconstructed independently across planner, descriptors, review, Inspector and export. SOL-002, SOL-004, SOL-006 and AUD-009 demonstrate actual divergent behavior. These are sufficient reasons for a **small shared pure projection/compatibility contract**, not a broad UI rewrite or a new scientific data model. Large file sizes alone are not findings.

## Migration and deployment implications

1. **Index readiness is a real gate.** The additive catalog schema/projection marks old rows dirty. The independent migration probe confirms coverage refuses an incomplete index and succeeds after ordinary reconciliation. Run and verify that existing reconciler before making the upgraded assessment surface available; do not interpret the refusal as missing scientific work. Refreshing Settings can also change the Results root, so verify the catalog's storage-root identity.
2. **Old scientific records remain historical.** The automatic-attachment fix does not repair previously incorrect links. Legacy external structures without retained digest evidence are intentionally refused by verified access; legacy CAS filenames supply independent evidence where supported. Inventory actual legacy records in a separately authorized migration assessment before rollout rather than inventing hashes from current files.
3. **Review generations are local recovery metadata.** Back up the existing SQLite database as a consistent unit, including its review-generation table/WAL recovery state. Peer review wire identity is unchanged. The peer tests show stale forms reject replayed edits/clears after restart, while offline concurrent rationales remain in existing conflict history rather than being merged automatically.
4. **Mixed binaries remain a compatibility hazard.** An old writer can bypass form-level CAS; its write invalidates new editors but may already have overwritten the newest rationale. An old catalog writer can also leave projections requiring refresh. Update all review-writing processes together or make the old clients read-only during a controlled transition. No live mixed-version rollout was tested.
5. **Immutable snapshots change caller assumptions.** Some broad-suite fixtures still expect the original filesystem path or try to mutate the returned read-only snapshot. Those expectations must be reconciled with the new API contract; weakening integrity to make them pass would be the wrong remedy.
6. **The committed baseline is not the running release.** Findings about this SHA and operational observations about `p11-25e003c` must remain separate. This audit does not certify deployment compatibility on the three workstations, container equivalence, model checkpoints, production data shape or a completed production migration.

## Test execution and evidence

**Runtime:** `/home/user/mambaforge/envs/mn-molecular-design-studio/bin/python`, Python 3.11.17; pytest 9.1.1, DuckDB 1.5.6, RDKit 2026.09.1, Streamlit 1.65.0 as reported by the installed environment. No package upgrades or new dependencies were installed.

| Run | Result | Interpretation |
|---|---|---|
| Focused committed correction/peer suite, 10 modules | **191 passed, 26 warnings, 35.43 s** | Independently reruns identity/integrity, artifact finalization, assessment planning/index/lifecycle, CAS editors, peer review convergence and Tag/Fusion contracts. Warnings are existing Biopython alignment-attribute deprecations. |
| Entire committed offline suite in the corrected OS sandbox | **3,006 passed, 14 failed, 55 skipped, 103 warnings, 934.24 s** | 3,075 JUnit test records, zero collection/setup errors; full invocation, not the reported 663-test subset. Individual failures and scope limits are discussed below. |
| Final independent adversarial harness, 25 parameterized cases | **20 failed, 5 passed, 21.24 s** | Failures are assertions of the intended behavior that expose the defects in this report, including six CSV-prefix cases. Positive controls cover source integrity, delete/restart CAS, failed-QC evidence blocking reuse, stable-key idempotency and catalog upgrade/reconciliation. These are not counted as passing application regressions. |

Focused modules:

```text
tests/test_molecular_database_identity_integrity.py
tests/test_artifact_collection.py
tests/test_assessment_planning_audit.py
tests/test_assessment_submission_audit.py
tests/test_assessment_job_index_audit.py
tests/test_review_concurrency_audit.py
tests/test_tag_fusion_applicability_audit.py
tests/test_peer_workstations.py
tests/test_molecular_database_peer.py
tests/test_complex_analysis_tag_fusion.py
```

The implementation's earlier 663-pass result is a selected-module regression result. I did not treat it as a claim that every committed repository test passes. The focused independent run confirms the principal corrected contracts while the broader run challenges additional surfaces.

### Broader committed-suite failures

All 14 failures reproduced in two subsequent isolated targeted runs: 10 failed in 6.01 s and the remaining four failed in 1.54 s. They are not explained away as cross-module test ordering. Several test assumptions have visibly drifted from the changed contracts, so the 14 failures must also not be equated with 14 confirmed application defects.

| Tests / count | Observed failure and assessment |
|---|---|
| `test_accepted_design_import.py` / 1 | The test imports the removed `_render_accepted_design_promotion` helper. This is a stale test/API reference; the intended replacement user flow needs an explicit maintained contract. |
| `test_analysis_container_ownership.py` / 1 | A fixture labels a derived projection as v4 and forbids calling `normalize_metrics`; the current v6 projection rebuild calls it. This is a projection-upgrade expectation disagreement, not evidence that a scientific engine recalculated historical data. |
| `test_bound_docking.py` / 1 | Asserts that the returned receptor pathname equals the original CAS pathname. The new verified-reader contract intentionally returns a private snapshot. Assert verified bytes/identity rather than weakening that contract. This failure is distinct from the actual legacy viewer bypass in SOL-006. |
| `test_complex_analysis_dashboard.py` / 2; `test_results_board_identity.py` / 1 | Tests still expect run-first navigation/state, while the current UI initializes Database-first state. The code/test behavior disagreement predates these fixes. Refresh explicit navigation/deep-link acceptance tests around the intended current architecture. |
| `test_complex_analysis_placement_workflow.py` / 2 | One test attempts to mutate the newly read-only snapshot and fails before its stale-reference assertion. The other expects a launch button without satisfying the new indexed-lookup readiness precondition. Correct fixture setup and verify the actual refusal/recovery behavior. |
| `test_inverse_folding.py` / 1 | A candidate's run-relative coordinate path is opened as a current-working-directory path. Verify the stored-path resolver contract before deciding whether the test or a consumer is wrong; this audit does not claim that every downstream consumer is safe. |
| `test_ligand_tools.py` / 1; `test_pipeline_graph.py` / 1; `test_rfd4_proteina.py` / 1 | Mocked execution is expected to finish successfully despite incomplete declared artifacts. The rerun's shared-branch result explicitly reports missing `child-rfd4/designed.cif`; the RFD4 queued-fold result explicitly reports missing declared `candidates.jsonl`. The GROMACS stub also receives `success=false`. These fixtures must retain real declared outputs when testing success. The strict finalization gate must not be weakened merely to restore their old expectations. |
| `test_measurement_semantics.py` / 1 | The legacy default ranking call supplies a nonzero selection limit even with no rules/ranking measurement; `evaluate_selection()` treats that as an explicit policy and records `pass`. The test expects `not_evaluated`. This is a real unresolved selection-semantics disagreement requiring a product/scientific contract decision. Current final-acceptance projection does not automatically promote this `selection` field to final computational acceptance. |
| `test_scientific_tasks.py` / 1 | Starter-template count is 14 while the test expects 13. Update the intended catalogue contract rather than inferring a scientific failure from the count alone. |

The complete test IDs, assertions and stack traces are retained in the broad JUnit/log and the two targeted triage logs. No committed test was edited, filtered out of the full invocation or newly marked skipped. Targeted triage invocations name the failed tests explicitly.

The 55 skips comprise 18 opt-in GPU cases, 18 unavailable image-side batch-API-source cases, 10 opt-in QC-container cases, three opt-in RFD4-container cases, three missing `mdtraj` cases and three missing `chaperone` cases (including two module-level skips). The 103 warnings comprise 80 Biopython alignment-attribute deprecations, 15 DuckDB batch-reader deprecations and eight `pydssp` numerical warnings on fixtures. A redundant timing rerun of tail modules was stopped after the full suite completed; it is not counted as a completed run or substituted for the complete results.

### Reproduction and retained evidence

All final tests use the corrected launcher:

```bash
cd /home/user/.codex/worktrees/independent-sol-audit/mn-molecular-design-studio
bash /tmp/mds-independent-sol-20261111/run-isolated-tests.sh replay -q
bash /tmp/mds-independent-sol-20261111/run-isolated-tests.sh replay-probes \
  -q -s --tb=short /tmp/mds-independent-sol-20261111/probes/test_post_phase2_audit.py
```

The second command is expected to report the negative contract assertions until the defects are corrected. It uses disposable fixture data and disables queue dispatch in submission probes.

Preserved local evidence under `/tmp/mds-independent-sol-20261111/`:

- `run-isolated-tests.sh`: exact mount/network/environment isolation and Python command.
- `probes/test_post_phase2_audit.py`: final independent reproductions and positive controls.
- `focused.log`, `focused.xml`: the 191-test verification run.
- `full-suite-sandboxed.log`, `full-suite-sandboxed.xml`: completed broad run.
- `probes-final.log`, `probes-final.xml`: final 25-case adversarial run.
- `native-probe.log`: real chemistry round trip plus pre-Phase-2 comparison.
- `deep-probes.log`: concurrent admission and frame-invalid receptor HTTP responses.
- `failure-triage.log`, `failure-triage.xml`, `extra-triage.log`, `extra-triage.xml`: independent reruns of all 14 broad-suite failures.
- `full-suite.log`: interrupted initial run, retained for transparency and excluded from completed totals.
- `escaped-results/`: the two quarantined audit-created, unexecuted queue directories.

An accompanying `INDEPENDENT_SOL_AUDIT_EVIDENCE.tar.gz` preserves the final scripts, logs and JUnit records outside temporary-directory retention. It contains fixture/audit evidence, not a copy of production scientific data. Its checksum and final baseline/isolation checks are recorded below.

Archive SHA-256: `5e5378d59548d4eabf410c1812aa2b05de146c3e8696031c58f3abf30f657cdb`.

`final-state.json` confirms the detached audit HEAD and primary development HEAD still equal the full recorded baseline, no tracked application file differs, the primary checkout retains exactly its original seven untracked documentation entries, both accidental audit IDs are absent from Results and the two checked production tables, and both primary service PIDs remain unchanged (`3653828`, `3654759`). The archive includes the read-only verification script and its output.

## Coverage limitations

- This was a code/contract audit and an offline fixture test exercise, not a new biological validation campaign. No experimental affinity, stability, expression or tag-performance inference is warranted.
- No GPU inference, installed QC/RFD4 containers, live distributed scheduler submission, live remote peer rollout or browser-to-production mutation was intentionally performed. The initial isolated-queue-path mistake and its correction are disclosed above.
- Peer convergence, replay, clears, restarts and conflict preservation were tested through the repository's two-peer temporary fixtures, not by changing real peers. The concurrent admission proof uses one API instance with real queue records and dispatch disabled.
- No production database migration, legacy-link repair, digest backfill, historical re-finalization, authentication transition or active-branch modification was performed. Production data was not used as a mutable scientific fixture; the two accidental queue directory creations are separately disclosed.
- The indexed timing/instrumentation tests establish exact-source lookup behavior and zero broad Results reads in batch coverage; they do not benchmark a full production page, cold remote owners or an entire legacy catalog migration.
- Unknown historical method/image versions remain unverified. The QC image-digest probe demonstrates what the comparison ignores; it does not establish that installed production image bytes actually differ.
- Existing skipped dependencies and GPU/container cases are listed with the test results. No passing-test count should be interpreted as coverage of those unavailable runtime paths. Native CPU fixture calculations are exercised; these do not constitute a production scientific validation campaign.

## Prioritized recommendations

1. **Restore scientific boundary consistency before relying on pooled conclusions:** enforce protocol checks on all reuse surfaces (SOL-001), guard target revisions during comparison (SOL-003), and use the verified receptor/frame resolver for visualization (SOL-006).
2. **Correct native reuse and catalog recovery:** canonicalize effective fixed/default settings to fix the confirmed native chemistry regression (SOL-002), and contain actual SQLite failures with a recoverable unavailable state (SOL-008). Verify real round trips rather than only normalized synthetic maps.
3. **Unify the small scientific presentation contract:** separate subject selectors from scientific protocol signatures (SOL-004), preserve numeric-zero models (SOL-007), and make metadata-only exports complete (AUD-009). Reuse existing modules and retain stable historical IDs.
4. **Finish campaign review continuity and safe export:** resolve saved membership models, filter before pagination, preserve selections across pages, and escape spreadsheet string cells (AUD-007–010).
5. **Make duplicate-prevention guarantees accurate:** preserve conservative lifecycle checks, but add an atomic scientific admission key in the existing queue/coordinator if cross-session/peer deduplication is required (SOL-005). Continue supporting explicit reassessment.
6. **Stage deployment and authentication separately:** verify catalog readiness, legacy compatibility, snapshot consumers and a coordinated client upgrade. Reconcile the broad committed-suite failures. Only then perform a separately authorized release/account transition; the current live development identity remains an operational exception (AUD-013).

The implementation should retain the verified Phase-1 integrity and Phase-2 CAS protections while closing these specific gaps. No broad rewrite is justified by this audit.
