# Molecular Design Studio — independent application audit

Audit date: 2026-11-11. Baseline: `25e003cdb1a3198f15a4f425342ed3671288bbbb`, branch `codex/p0-production-visibility`.

This is a bounded, risk-prioritized technical audit, not a certification of the entire application. Findings below come from source inspection, newly executed isolated checks, and a read-only observation of the running application. Earlier development approvals were not used as evidence of correctness. **No application code, production scientific records, services, or replication state were changed.**

## History and maintenance of this document

- 2026-11-11: Initial root `AUDIT.md`; no previous root audit existed. Existing audit documents under `docs/` and unrelated untracked cleanup material were left untouched. This document does not supersede their conclusions or claim to have re-audited them.
- IDs `AUD-001` through `AUD-014` are stable. Preserve them when implementing fixes. Append resolution evidence, commit, and date rather than removing the original finding. A changed conclusion must explicitly explain the new evidence.
- Line numbers refer to the baseline above. Function names remain the primary locator after edits.
- 2026-11-11: Implementation phase 1 completed for AUD-001/AUD-002/AUD-003. The original audit text describes the original audit operation; later implementation changes and verification are recorded separately below. Production services, data, authentication and replication state remain unchanged.

## Overall assessment

The application has substantial useful infrastructure: project authorization, immutable Evaluation identities, explicit acceptance stages, source-qualified peer metadata, portable artifacts, and a functional production review page. The current tests cover many successful workflows and important negative cases.

However, **scientific identity and lifecycle guarantees are not enforced uniformly at subsystem boundaries**. The most important reproduced problems are incorrect model attachment during observation ingestion, unchecked local structure bytes, an artifact-finalization classification bypass, unsafe assessment-reuse decisions, and lost human rationale through a stale form. These warrant focused corrections before treating the review/assessment workflow as dependable for unattended or multi-user use. This audit found no evidence that production records have already been corrupted; the defect demonstrations used disposable fixtures.

The storage architecture does not need replacement to address these findings. Use the existing Molecular Database, Evaluation, scheduler, Results, and peer contracts. The recurring weakness is inconsistent application of those contracts across older and newer paths.

## Priority index

| ID | Severity | Status | Finding |
|---|---|---|---|
| AUD-001 | High | Confirmed | Observation auto-attachment ignores model identity |
| AUD-002 | High | Confirmed | Local structure retrieval does not enforce the recorded digest |
| AUD-003 | High | Confirmed | Temporary-path heuristic permits dangling scientific outputs |
| AUD-004 | High | Confirmed | Assessment reuse ignores evaluator version and depends on row order |
| AUD-005 | High | Confirmed | Stale shortlist forms overwrite newer rationale silently |
| AUD-006 | Medium | Confirmed | A collecting analysis is offered for duplicate submission |
| AUD-007 | Medium | Confirmed | Membership-only saved collections open with no review models |
| AUD-008 | Medium | Confirmed | Pagination clears selections and applies filtering after paging |
| AUD-009 | Medium | Confirmed | Metadata-only review export omits available sequence/provenance fields |
| AUD-010 | Medium | Confirmed | CSV export preserves spreadsheet formula input as executable text |
| AUD-011 | Medium | Confirmed | Assessment rendering scans Results files repeatedly and truncates job lookup |
| AUD-012 | Medium | Confirmed | Batch Tag/Fusion applicability contradicts the registered source adapter |
| AUD-013 | Medium | Confirmed | Running production interface uses a shared unauthenticated admin identity |
| AUD-014 | Low | Improvement | Scientific projection logic depends on large UI modules |

No Critical finding was established. “Confirmed” describes the specific evidence and scope stated below; it does not assert that every caller or historical record is affected.

## Findings

### AUD-001 — Observation attachment can select the wrong prediction model

**ID:** AUD-001  
**Severity:** High  
**Category:** Bug  
**Status:** Confirmed  
**Location:** `mn_molecular_design_studio/molecular_database.py`, `MolecularDatabase.add_observation`, lines 712–732; downstream identity consumers in `evaluations.py`, lines 70–108.  
**Problem:** When an observation lacks an explicit `representation_id`, the fallback chooses the most recent structure with the same run and optional candidate ID, ignoring an explicitly supplied `model_id`. It can persist contradictory model and representation identities.  
**Evidence:** An isolated Database contained one candidate with `model-1` and a later `model-2` in the same run. Calling `add_observation(..., context={"candidate_id":"c1", "model_id":"model-1"})` attached the observation to the `model-2` structure. The probe printed different expected and attached structure IDs. Existing review code may reject this contradictory record, producing missing evidence instead of displaying it; this does not make the stored association correct.  
**Root cause:** A legacy convenience fallback treats run/candidate equality as sufficient structure correspondence.  
**Expected behavior:** Evidence is attached only to its exact recorded model/representation. Ambiguous evidence remains explicitly unresolved or is rejected.  
**Recommended fix:** Match all supplied model, pose, candidate, and digest constraints. Require a unique compatible result before auto-attachment; never choose by recency among scientifically different models. Prefer explicit representation bindings from ingestion adapters. Do not retroactively rewrite historical observations as part of this fix.  
**Dependencies:** None. Coordinate with AUD-002 and AUD-004.  
**Verification:** Two models of one candidate in one run; reversed insertion order; explicit and absent model IDs; ambiguous candidates; pose IDs; zero-valued model indices; peer-qualified representations. Confirm stored context, Evaluation descriptor, review, and Inspector all agree, and repeated ingestion remains idempotent.

**Implementation — 2026-11-11:** Corrected automatic attachment to require exactly one source-group/run representation matching every supplied candidate, model, pose and source-digest selector. Ambiguous/unknown selectors remain unattached. Recorded calculation digests are preserved separately from source identity; zero model indices are canonicalized only on newly inserted observations. No historical observations or projection schema were rewritten. Added 14 isolated regressions, including insertion order, nested selectors, repeated ingestion and peer-qualified source groups. Eleven functional regression cases failed against the original implementation (one additional initial failure was a corrected replica-fixture setup error). Verification: 74 tests passed across identity, Evaluation/comparison, Database and peer suites. Fix commit recorded in the phase-1 ledger below. Status: implemented; production data reconciliation is deliberately outside scope.

**Follow-up — 2026-11-11:** An additional failing-before/passing-after regression covers numeric `model_id=0` as well as zero model indices: newly recorded IDs are string-normalized so the existing descriptor preserves zero. This does not rebuild or change historical Evaluation identities. The combined identity/integrity file now has 28 passing cases.

### AUD-002 — Local artifact reads can violate immutable structure identity

**ID:** AUD-002  
**Severity:** High  
**Category:** Bug  
**Status:** Confirmed  
**Location:** `molecular_database.py`, `structure_file`, lines 1450–1471; `molecular_database_api.py`, structure file/materialize routes, lines 266–273; `complex_analysis_sources.py`, Database source resolution, lines 262–290. All paths are under `mn_molecular_design_studio/`.  
**Problem:** Local retrieval checks authorization and file existence, but not whether current bytes match the representation's stored SHA-256. Materialization hashes the current bytes into a new input copy without rejecting a mismatch against the original representation. The analysis adapter also computes a current digest; its explicit comparison here applies to a selected Evaluation, which is optional.  
**Evidence:** Attach a local file, retain its recorded content digest, then change only that fixture file. `structure_file(user, structure_id)` returned the changed bytes without error: recorded digest `32429cc2…`, returned digest `3a2e0ca9…`. The direct REST download route delegates to this method. No production file was changed.  
**Root cause:** Integrity validation is distributed across callers and remote-transfer helpers instead of guaranteed by the common authorized read boundary.  
**Expected behavior:** An exact structure reference returns matching bytes or an explicit integrity error; it must not silently become another scientific representation.  
**Recommended fix:** At the shared retrieval/materialization boundary, verify the appropriate recorded content/portable digest before handing bytes to download, Inspector, or analysis. Handle legacy digest-less records explicitly without silently inventing certainty. Preserve the distinction between original, normalized, and calculation structures. Address concurrent replacement while reading, not just a preliminary `stat`.  
**Dependencies:** None; apply before expanding reuse or export reliance on exact identity.  
**Verification:** Changed local file, corrupted cached replica, missing file, authorized/unauthorized request, unchanged file, and legacy no-digest representation. Test direct download and Database analysis submission. Existing remote digest checks and Results ownership protection must remain intact.

**Implementation — 2026-11-11:** Authorized structure reads now validate the recorded portable/content SHA-256 and return an atomically written, read-only private snapshot of those exact bytes for both download and analysis. The original path cannot be swapped between validation and FileResponse opening it. Portable conversion uses the same verifier, new file imports retain a content digest, and existing corrupt content-addressed files are rejected rather than relabeled. Peer download integrity failures return an explicit 422 response. Thirteen added regression cases cover changed sources, cached peer copies, source replacement, real REST download/materialize, analysis without a selected Evaluation, missing/unauthorized access, cache repair from a verified source, and legacy handling. Eight regressions failed against pre-fix code; a separately corrected adapter-test import then reproduced the missing integrity gate as a downstream parser error. Verification: 191 tests passed across identity, Database/policy/integration, peer, Evaluation/comparison, collection and hydration suites. Fix commit recorded in the phase-1 ledger. Status: implemented with an explicit legacy limitation: a digest-less external legacy path cannot establish immutable identity and is refused for verified access; existing CAS filenames provide independent digest evidence. No historical metadata is backfilled or changed automatically.

### AUD-003 — Optional-path heuristics bypass scientific-output finalization

**ID:** AUD-003  
**Severity:** High  
**Category:** Bug  
**Status:** Confirmed  
**Location:** `mn_molecular_design_studio/engine/workflow_engine/core/artifact_collection.py`, `validate_run_artifacts`, lines 399–426; `core/jobs.py`, finalization calls around lines 618–619 and 668–670.  
**Problem:** A dangling coordinate file in an output directory is classified as optional if either the link name or target contains `/tmp/`, `/processed/`, or another optional marker. This check runs before the scientific output suffix/location check.  
**Evidence:** Fixture link `artifacts/raw/engine/output/model.cif -> /tmp/removed-container-output/model.cif`, with no candidate/manifest declaration, returned `status="complete"`, `required_artifact_count=0`, and `classification="optional_execution_artifact"` under the default strict validator. Explicitly declared required artifacts still fail correctly; the defect concerns undeclared outputs that the strict fallback is intended to catch.  
**Root cause:** Substrings describing execution locations are treated as authoritative evidence of dispensability. Required-output declaration coverage is incomplete, so the fallback matters.  
**Expected behavior:** Output semantics take precedence over a temporary container origin. Future successful collection cannot retain an indispensable inaccessible output merely because it came from a cache or temporary path.  
**Recommended fix:** Prioritize required manifests and scientific-output classification; require explicit optional-artifact evidence or a narrowly validated workflow classification for exceptions. Preserve historical links as links and retain the existing non-strict historical validation mode. Do not repair or delete historical Results.  
**Dependencies:** None.  
**Verification:** Output `.cif`, `.pdb`, and `.sdf` links into `/tmp`, `/processed`, `/workspace`, and `/boltz-context`; declared optional caches; a valid canonical copy plus obsolete alias; missing declared output; historical non-strict inventory. Test the actual finalization gate as well as the helper.

**Implementation — 2026-11-11:** Scientific coordinate files and entries in output/pose/prediction directories are classified before execution-cache heuristics. A target containing `/tmp` or `/processed` cannot declare its link optional. Optional manifests remain supported; a container alias may survive only when the exact resolved canonical file is already retained and validated. Unknown links still block strict finalization. Added 21 regression cases covering container origins, both normal finalization entry points, explicit optional aliases, exact canonical correspondence, missing canonical data, cached coordinates/confidence outputs, and preservation of historical links. Fifteen initial cases reproduced failures; two additional cache-payload cases also failed before the final classification correction. The previous manifest test now explicitly declares its optional link instead of relying on its target name. Verification: all 31 artifact cases passed; the combined targeted identity/artifact/collection/shortlist/comparison/catalog run passed 84 tests. Two pre-existing AppTest fixtures now restore their module replacements, preventing fake Database/comparison implementations leaking into subsequent regression tests. Status: implemented; final combined verification and fix commit appear in the phase-1 ledger below. Limitation: this is not a claim that every producer declares every possible scientific payload; wholly absent regular files with no manifest/candidate/downstream reference cannot be inferred. Workflows with opaque optional links must declare them explicitly. No historical Results are repaired, deleted, or re-finalized.

### AUD-004 — Compatible Evaluation reuse is incomplete and order-dependent

**ID:** AUD-004  
**Severity:** High  
**Category:** Consistency  
**Status:** Confirmed  
**Location:** `mn_molecular_design_studio/complex_analysis_batch.py`, `evaluation_coverage`, lines 182–230; `complex_analysis_collection_ui.py`, parameter construction and coverage, lines 777–813; `evaluations.py`, version/image/settings fields, lines 91–118.  
**Problem:** Reuse compares the requested parameter subset and identity fields, but has no expected evaluator version/image/protocol input. It also returns on the first exact matching record, even if that record failed and a compatible successful Evaluation exists later.  
**Evidence:** A matching descriptor with `engine_version="obsolete-engine-version"` was reported “Already evaluated · compatible settings”. Given the same two records, `[failed, completed]` returned “Previously failed”, while `[completed, failed]` returned “Already evaluated · compatible settings”. Thus enumeration order changes behavior. These were direct function probes.  
**Root cause:** Reuse is a handwritten subset comparison plus first-match status selection rather than an explicit scientific compatibility decision over all matching evidence.  
**Expected behavior:** Reuse is deterministic and based on the requested scientific protocol, exact input, and a documented version compatibility rule. A failed attempt cannot obscure an existing compatible success.  
**Recommended fix:** Supply expected method/protocol identity and normalized scientific settings from existing workflow contracts. Separate runtime settings from scientific parameters. Evaluate all matches; prefer an eligible success, then active work, then failed/mismatched evidence. Unknown version compatibility must be visible, not silently assumed. Keep all original Evaluations.  
**Dependencies:** AUD-001/AUD-002 establish reliable source identity. Group implementation with AUD-006 and AUD-011, without replacing the scheduler.  
**Verification:** Permute successful/failed/in-progress records; change method version, image/protocol, cutoffs, chain roles, model and digest; distinguish runtime-only changes; test missing version metadata and multiple successful conflicting records. No repeated calculation should occur solely because SQL returned rows in another order.

### AUD-005 — Stale human-review forms silently overwrite newer decisions

**ID:** AUD-005  
**Severity:** High  
**Category:** Bug  
**Status:** Confirmed  
**Location:** `complex_analysis_shortlist_export.py`, `render_shortlist_review_export`, lines 398–435; `molecular_database.py`, `save_complex_review`, lines 1707–1733; `molecular_database_peer.py`, capture base revisions, lines 529–540.  
**Problem:** The form reads the latest annotation but initializes widget state only once. A subsequent save writes all selected form values unconditionally, even if another session or peer has updated the annotation since this form was opened.  
**Evidence:** Isolated Streamlit AppTest: open “Original rationale”; save “Newer rationale saved by another session” directly through the existing Database method; rerun the form. Database held the newer rationale while the text area still showed the old one. Clicking Save restored “Original rationale” with no error. No real user annotation was edited.  
**Root cause:** There is no optimistic version check between form load and save. Peer conflict detection cannot protect a stale form saved after the new peer value has already been received: capture treats that save as a new sequential edit.  
**Expected behavior:** Researchers see intervening edits and explicitly reconcile them. An untouched stale field must not silently erase a newer rationale.  
**Recommended fix:** Bind the editor to the saved review version/base value, refresh untouched fields on revision changes, and perform compare-and-set at the write boundary. Save only intentional changes. Reuse existing conflict/history semantics; do not create a second shortlist store. Include user identity in editor state scope.  
**Dependencies:** None. Exercise existing peer synchronization after the fix.  
**Verification:** Two local sessions and two staged peers; edit versus edit, edit versus clear, refresh with unsaved typing, saving a multi-design selection with one stale member, and repeated unchanged saves. Verify newer rationale survives or a visible conflict is raised; exports must reflect persisted decisions.

### AUD-006 — Collecting jobs do not block duplicate analysis submission

**ID:** AUD-006  
**Severity:** Medium  
**Category:** Bug  
**Status:** Confirmed  
**Location:** `mn_molecular_design_studio/complex_analysis_collection_ui.py`, `_render_batch_assessment`, lines 816–846; `complex_analysis_ui.py`, `_matching_database_analysis_jobs`, lines 2508–2522.  
**Problem:** The job-state classifier recognizes several active states but omits `collecting`. An exact matching collecting job leaves the row classified as “Missing evaluation” and eligible to submit again.  
**Evidence:** An isolated AppTest mocked only the job lookup to return an exact `collecting` match. Coverage remained “Missing evaluation”, the submit button was enabled, and clicking it invoked the stub submission callback. `collecting` is an existing application state, including `cluster/candidate_batches.py` and distributed campaign collection. No scheduler job was actually submitted.  
**Root cause:** UI-local status enumerations differ from scheduler/collection lifecycle states; unrecognized states fall through to the previously computed missing state.  
**Expected behavior:** Collection/transfer/finalization in progress prevents equivalent resubmission. Unknown states require clarification or refresh rather than being treated as absence.  
**Recommended fix:** Reuse one existing lifecycle classification for active execution and collection. Make unknown/matching states conservative and visible. Keep explicit retry controls for terminal failures and retain per-design completion status.  
**Dependencies:** Group with AUD-004/AUD-011.  
**Verification:** Queued, running, collecting, transfer/finalization, failed, cancelled, unknown and completed-but-not-ingested states; inspect coverage and assert submission callback counts. Preserve queued jobs and do not force-finalize records.

### AUD-007 — Saved collections with membership references appear empty in review

**ID:** AUD-007  
**Severity:** Medium  
**Category:** Bug  
**Status:** Confirmed  
**Location:** `molecular_database.py`, `add_to_group`, lines 1312–1324 and `save_comparison`, lines 1292–1296; `complex_analysis_design_review.py`, `design_review_page`, lines 49–54. Paths under `mn_molecular_design_studio/`.  
**Problem:** Normal saved groups add molecule memberships without copying structures. The new review query reads only structures whose owning `group_id` is the selected group, so membership-only curated/comparison collections cannot supply their existing accessible models.  
**Evidence:** Create one source molecule with two retained models; create a saved group and call `add_to_group`. Source review returned two models, saved-group review returned zero, despite the new group containing the molecule. This uses normal public Database operations.  
**Root cause:** Selection membership and structure ownership are conflated in the review query.  
**Expected behavior:** A saved collection resolves its intended accessible retained representations without copying coordinates or silently selecting unrelated models of the same molecule.  
**Recommended fix:** Resolve explicit saved representation/Evaluation membership through existing relationships. Where a molecule-only collection does not identify a particular model, expose the compatible accessible models for explicit choice. Honor frozen comparison Evaluation IDs. Never broaden ACLs merely because a molecule appears in another group.  
**Dependencies:** AUD-001/AUD-002; coordinate query changes with AUD-008/AUD-009.  
**Verification:** Membership-only curated group, frozen comparison, source-owned and remote structures, multiple models/targets of one molecule, revoked source access, and repeated import. Keep exact structure IDs and avoid duplicate molecule/structure creation.

### AUD-008 — Paging discards selection and makes filtering page-local

**ID:** AUD-008  
**Severity:** Medium  
**Category:** UX  
**Status:** Confirmed  
**Location:** `complex_analysis_design_review_ui.py`, `render_design_review`, lines 65–84, 100–118 and 151–157; `complex_analysis_design_review.py`, `remember_review_selection`, lines 183–196. Paths under `mn_molecular_design_studio/`.  
**Problem:** The Database is paged before target, acceptance, and search filters are applied. The selection is then intersected with only the current page's filtered rows. Moving pages loses selected models and their hydrated state, preventing comparison across pages and making a search miss valid records on other pages.  
**Evidence:** Calling the existing selection helper with a selected page-1 structure and page-2 allowed IDs returned an empty selection and removed `complex-analysis-workflow-samples`. The caller derives those allowed IDs from the already limited page. Target filtering at lines 69–81 likewise occurs after the query.  
**Root cause:** Page visibility is treated as the authorization/selection universe. Filtering and pagination are applied in the wrong order for a campaign-wide workflow.  
**Expected behavior:** Search/cohort filters cover the authorized collection, counts describe that query, and paging does not clear valid selections.  
**Recommended fix:** Filter indexed metadata before paging, retaining bounded hydration. Keep stable selected model references independently of the visible page; revalidate their ACL and cohort eligibility when acting. Make any intentional selection reset explicit. Bound the page control to available pages.  
**Dependencies:** AUD-007; coordinate with AUD-011.  
**Verification:** At least 251 models; selected models on two pages; an accepted/search-matching model only on page 2; target filter; sorting; revoke access; change campaign. Verify analysis, comparison, shortlist and return from Inspector use the same intended selection.

### AUD-009 — The new lightweight review projection loses export context

**ID:** AUD-009  
**Severity:** Medium  
**Category:** Consistency  
**Status:** Confirmed  
**Location:** `complex_analysis_design_review.py`, record/sample construction, lines 120–142; `complex_analysis_shortlist_export.py`, `_sequence_by_chain`, lines 48–57, `_source_lineage`, lines 127–145, and `build_shortlist_record`, lines 181–186 and 278–284. Paths under `mn_molecular_design_studio/`.  
**Problem:** The metadata-only review passes `components` but no parsed `tracks`/numbering. The exporter reads chain sequences only from tracks/numbering, so already recorded target/binder chain sequences are omitted. The review record also omits `sources`, which the lineage exporter expects.  
**Evidence:** A review sample with saved target chain A=`ACDEFG` and binder B=`KLMNPQ` exported empty `target.sequences_by_chain` and `binder.sequences_by_chain`. The candidate-level binder sequence survived, showing this is a projection mismatch, not missing source data. Source lineage cannot be emitted from this path because its record has no `sources` field.  
**Root cause:** The exporter was designed for hydrated samples; the new lightweight review reuses it without satisfying its metadata contract.  
**Expected behavior:** JSON and CSV preserve available recorded sequences and source relationships without requiring coordinate loading. Missing data must be identified as missing only when it is actually unavailable.  
**Recommended fix:** Add an exact chain-aware fallback to recorded components, preserving author/source chain mapping. Read only the selected authorized source relationship metadata needed by the export. Define the minimal shared review/export projection explicitly; do not hydrate all coordinates.  
**Dependencies:** Coordinate with AUD-007/AUD-008 and AUD-010.  
**Verification:** Export from the metadata-only Review & assess path before opening 3D; compare JSON and CSV against saved components and source records. Cover renamed chains, remote provenance, multiple prediction models, unknown fields, and differing acceptance stages. Verify actual downloads in a staging browser after serializer tests.

### AUD-010 — Human text can become a formula in CSV exports

**ID:** AUD-010  
**Severity:** Medium  
**Category:** Security  
**Status:** Confirmed  
**Location:** `mn_molecular_design_studio/complex_analysis_shortlist_export.py`, `manifest_csv`, lines 322–360.  
**Problem:** User-entered rationale and other text fields are written to CSV without neutralizing spreadsheet formula prefixes. Quoting CSV fields does not make them literal spreadsheet text.  
**Evidence:** A fixture rationale `=1+1` emerged from `manifest_csv` as the same cell value `=1+1`, with no literal-text prefix. No spreadsheet formula or external command was executed in the audit.  
**Root cause:** Serialization handles CSV syntax but not the convenient spreadsheet-consumption use case.  
**Expected behavior:** Opening the laboratory-planning CSV cannot execute researcher/imported text as a formula. JSON should preserve the original exact text.  
**Recommended fix:** Apply a documented spreadsheet-safe escaping policy to string cells beginning with dangerous formula/control prefixes, including imported IDs/names and rationale. Preserve real numeric values and JSON fidelity. Do not treat ordinary CSV quoting as mitigation.  
**Dependencies:** Group with AUD-009.  
**Verification:** `=`, `+`, `-`, `@`, tab/CR prefixes, quotes, delimiters, newlines, Unicode and legitimate negative numeric measurements. Test the raw CSV and literal display in a spreadsheet; do not execute a malicious payload for validation.

### AUD-011 — Assessment coverage reintroduces Results scans into page reruns

**ID:** AUD-011  
**Severity:** Medium  
**Category:** Architecture  
**Status:** Confirmed  
**Location:** `complex_analysis_collection_ui.py`, `_render_batch_assessment`, loops at lines 749–823; `complex_analysis_ui.py`, `_matching_database_analysis_jobs`, lines 2417–2488; `complex_analysis_design_review_ui.py`, target relationship enumeration, lines 69–80. Paths under `mn_molecular_design_studio/`.  
**Problem:** Each design × missing-method pair searches up to 1,000 jobs and reads candidate job metadata/manifests from Results during rendering. The same search is repeated on subsequent widget reruns. The hard first-1,000 limit also excludes older matching jobs. Target-scoped review separately enumerates all relationship pages before filtering the visible model page.  
**Evidence:** Direct code trace: nested sample/method loop → `_matching_database_analysis_jobs` → `search_jobs(limit=1000)` → per-row `metadata.json` and analysis manifest reads. There is no paging of the job search or shared per-render lookup. This audit did not benchmark the full production latency or assert a particular number of seconds.  
**Root cause:** Existing Results-driven deduplication was reused in a Database-driven page; indexed execution identity/status is not exposed to the coverage planner.  
**Expected behavior:** Ordinary review and coverage rendering uses bounded Database/catalog queries; exact selected Results access occurs only when needed. Older matching work remains discoverable.  
**Recommended fix:** Query the existing indexed job/Database source identity and settings once per batch, or add a compact index/projection to that existing catalog. Avoid a new service or data model. Reuse per-render results, push target constraints into queries, and preserve a final exact-source check at submission. Failed lookup must not be interpreted as confirmed absence.  
**Dependencies:** Group with AUD-004/AUD-006; coordinate pagination with AUD-008.  
**Verification:** Instrument file/network accesses on render; zero broad Results traversal. Use >1,000 historical jobs with an older exact match, multiple selected models and methods, remote owner offline, and catalog errors. Record cold/warm render and rerun times before/after; assert no duplicate submissions.

### AUD-012 — Tag/Fusion batch rules reject supported partner assemblies

**ID:** AUD-012  
**Severity:** Medium  
**Category:** Consistency  
**Status:** Confirmed  
**Location:** `complex_analysis_batch.py`, method profile lines 42–44 and `assessment_applicability`, lines 169–177; `complex_analysis_sources.py`, lines 415–417. Paths under `mn_molecular_design_studio/`.  
**Problem:** Tag/Fusion is categorized as `single_ppi`, enforcing exactly one binder and one target chain. The existing Database source adapter explicitly permits one assessed component and one or more disjoint partner chains.  
**Evidence:** `assessment_applicability(sample, "analysis.tag_fusion_suitability", ["B"], ["A", "C"])` returned false with “exactly one chain in each group”. The adapter's condition is `len(binder_chains) != 1 or not target_chains`, so the same partner assembly passes that specific prerequisite.  
**Root cause:** Reuse of the Rosetta/ipSAE single-pair category over-constrains another evaluator.  
**Expected behavior:** UI eligibility agrees with the registered workflow's scientific input constraints.  
**Recommended fix:** Give Tag/Fusion its correct capability constraint using the existing profile mechanism. Preserve the single-pair restrictions of Rosetta and ipSAE; do not weaken termini/topology checks.  
**Dependencies:** None; test with AUD-004.  
**Verification:** One component with one and multiple partner chains, overlapping chains, absent partner, and multiple assessed components. Confirm UI, adapter and workflow validation agree without executing scientific calculations.

**Phase-2 implementation — 2026-11-11:** Tag/Fusion now uses the existing profile mechanism with a one-component/multiple-partner capability. Non-overlap and PPI checks remain, Rosetta/ipSAE retain their single-pair restriction, and construct/termini validation is unchanged. A multi-partner regression failed before the correction. Five applicability cases plus the existing registered Tag/Fusion adapter/workflow suite passed: 22 tests, 26 existing Biopython deprecation warnings. Commit and group-level verification are recorded in the phase-2 ledger. No scientific calculation or production operation was performed.

### AUD-013 — Production currently relies on localhost isolation, not user authentication

**ID:** AUD-013  
**Severity:** Medium  
**Category:** Security  
**Status:** Confirmed  
**Location:** `mn_molecular_design_studio/api.py`, `_session_user`, lines 97–111; `auth.py`, development identity, lines 34–56; running primary `studio-api.service` and `studio-ui.service` environment.  
**Problem:** The operational application uses `MN_MOLECULAR_DESIGN_STUDIO_AUTH_MODE=development`. Requests receive the synthetic administrator and human decisions share the developer identity. Per-user/project authorization is bypassed for clients able to reach this endpoint.  
**Evidence:** Read-only `/proc/<MainPID>/environ` inspection confirmed that setting for both primary services. Unauthenticated `GET http://127.0.0.1:8766/api/v1/auth/me` returned HTTP 200, `username=developer`, `role=admin`, `development_mode=true`. The browser visibly showed “Development mode · sign-in disabled”. Socket inspection showed API/UI bound to `127.0.0.1`, **not a public listener**. No unauthenticated mutation was attempted; remote workstation auth configuration was not rechecked.  
**Root cause:** An explicit development configuration remains the production operating mode. This is an operational trust-model limitation, not an accidental default in `auth.py` (whose default is production).  
**Expected behavior:** If the application is used by multiple researchers, authenticated identities and project permissions must govern access and authorship. A deliberately single-user trusted-local deployment must state that limitation explicitly.  
**Recommended fix:** Plan a controlled transition using existing account/identity facilities, with verified account access and ownership continuity first. Do not flip production auth blindly or rewrite existing developer-owned records. Preserve localhost binding and peer authentication. Record any intentional single-user exception.  
**Dependencies:** Operational decision/account readiness; AUD-005 should include real user identities.  
**Verification:** Staged production-auth configuration: unauthenticated 401, password-change enforcement, owner/member/outsider/admin access, peer user/project propagation, and review author identity. Check all three hosts before a separately authorized operational change. Existing development tests must remain explicitly isolated.

### AUD-014 — Scientific presentation projections are coupled to UI internals

**ID:** AUD-014  
**Severity:** Low  
**Category:** Maintainability  
**Status:** Improvement  
**Location:** `complex_analysis_design_review.py`, line 41 imports `_database_candidate_view` from `complex_analysis_collection_ui.py`; `complex_analysis_collection_ui.py` contains 2,830 lines, `complex_analysis_ui.py` 4,417 lines at this baseline.  
**Problem:** A metadata-query module imports scientific identity/acceptance reconstruction from a Streamlit module. Evidence and eligibility reconstruction is spread across review, hydrated collection, Inspector, assessment and export paths. Changes can fix one path while leaving another inconsistent.  
**Evidence:** The import and the boundary mismatches in AUD-004, AUD-009 and AUD-012 are directly observable. File size alone is not a defect and does not justify a rewrite.  
**Root cause:** Successive UI extensions reused private helpers rather than a small shared presentation contract.  
**Expected behavior:** UI components render a common, testable projection while source identity, applicability and evidence binding remain independently usable.  
**Recommended fix:** During the relevant fixes only, extract the minimal pure projection/compatibility helpers into existing non-UI modules. Avoid a broad UI rewrite, new application or replacement scientific model. Keep adapters for current callers until equivalent behavior is tested.  
**Dependencies:** Let AUD-004/AUD-009/AUD-012 determine the smallest extraction; do not delay high-severity fixes for cleanup.  
**Verification:** Contract tests compare metadata review, hydrated Inspector, assessment planning and export for the same exact model, including remote aliases and unknown fields. Importing projection helpers should not require constructing Streamlit UI state.

## Implementation ledger — phase 1 (2026-11-11)

The first implementation phase is complete. Original findings and their evidence remain above; these resolution entries apply to the corrected code and covered contracts, not to retrospective production-data reconciliation.

| Finding | Resolution | Implementation commits |
|---|---|---|
| AUD-001 | Resolved for automatic attachment of new observations: unique exact candidate/model/pose/source-digest match; ambiguity remains unattached; zero-valued model selectors retained | `2a2534014aece5ec46cf5849126f15504b001f4b`, follow-up `cc984c51e89400304cc68640842d094bc2e8aff2` |
| AUD-002 | Resolved at authorized local/remote file read, download, materialization and portable-copy boundaries, with explicit handling of unverifiable legacy paths | `ab91b40b1fbab330f66cb5dd6b4eea6bdf45a6d5` |
| AUD-003 | Resolved for strict finalization of referenced/declared artifacts and dangling coordinate/output paths; historical non-strict inventory preserved | `2f20f9b1a6d165914c5e089e5472fabdb3308c4e` |

**Final verification:** 482 passed, 4 warnings in 187.58 seconds, in one combined invocation of 22 test modules. This includes 49 newly added parameterized regression cases (28 identity/integrity and 21 artifact cases), the existing Database/policy/integration, Evaluation/comparison, Explorer/hydration/navigation, shortlist, independent-peer, scientific-lineage, authenticated API, remote Results and scheduler-worker suites. No tests were filtered out for this run. The four warnings are the existing `pydssp` invalid-divide warnings in tiny validation-acceptance fixtures.

Test roots (Studio home, app home, Results, references and cluster) were redirected to `/tmp/mds-audit-fixes-20261111/verified/`; GPU tests were disabled and scheduler mode set to local. Authentic production scientific records were not used as mutable fixtures. Final combined log: `/tmp/mds-audit-fixes-20261111/phase1-verified.log`. Targeted final run: 84 passed in 6.08 seconds; log `phase1-targeted-final.log` in the same temporary evidence directory. Reproduction and intermediate logs remain there; source regression tests are committed in the repository.

**Test-isolation corrections:** The initial combined run had 428 passing tests, 30 failures and 21 setup errors after a pre-existing AppTest left `MolecularDatabase = FakeDatabase` installed in the shared module. Registering cleanup removed those downstream failures. A second run had 477 passing tests and four failures caused by a separate existing AppTest leaving `render_comparison` replaced by a stub. That fixture now uses scoped monkeypatching. The targeted order-sensitive run and the final 482-test run passed after both cleanup corrections. These were test harness repairs, not changes to application behavior or test exclusions.

**Remaining risks and scope:**

- A legacy external structure without a recorded digest cannot be scientifically verified. Verified access now reports this explicitly; no historical hash or identity is invented. Legacy CAS paths remain verifiable from their recorded content-addressed filenames. Any retention/re-import decision for external legacy sources must be explicit and preserve the old record.
- Artifact finalization cannot infer a wholly absent regular file that no producer declares or references. Opaque optional execution links must be declared in the existing artifact manifest. This phase did not certify every engine's manifest coverage.
- Historical observations are not rewritten. Discovery or correction of any already incorrect historical association requires a separate, evidenced operation.
- AUD-004 through AUD-012 remain open for later implementation phases; AUD-014 remains a limited consolidation recommendation. AUD-013's operational authentication transition is explicitly deferred pending separate authorization/account readiness.
- No push, deployment, service restart, authentication change, production migration or scientific-data modification was performed. The read-only service check still identifies both primary production services with release directory `p11-25e003c`; these fixes are in the main application repository for a subsequent authorized release.

Next phase: assessment compatibility/lifecycle and indexed job lookup (AUD-004/AUD-006/AUD-011, with AUD-012), followed by human-review concurrency (AUD-005), following the guidance below. Phase-1 completion is not a declaration that the remaining High findings are resolved or that a production rollout is authorized.

## IMPLEMENTATION_GUIDANCE

### Recommended order and grouping

1. **Scientific identity and completion:** AUD-001, AUD-002, AUD-003. Separate small commits are appropriate; verify explicit and legacy inputs before merging. These guard persistent evidence and successful job finalization.
2. **Assessment planning and lifecycle:** AUD-004, AUD-006, AUD-011 together. Share the exact source/protocol match and active-state interpretation; measure query/file access and remove the first-1,000 lookup limitation. AUD-012 is a small companion capability fix.
3. **Human decision safety:** AUD-005 independently, including concurrent sessions and staged peers. Preserve existing Lamport conflict handling and annotation identities.
4. **Selection and export continuity:** AUD-007/AUD-008 followed by AUD-009/AUD-010. Test the complete campaign → review → assessment/3D → shortlist → export path with multiple pages and saved collections.
5. **Operational authentication:** AUD-013 through a separately controlled configuration change after account/ownership readiness. Do not use this audit as permission to interrupt production or change credentials.
6. **Limited consolidation:** AUD-014 only where needed to keep the above corrections consistent.

### Constraints that must remain intact

- Keep one Molecular Database and its existing Evaluation and synchronized human-review contracts. Avoid a parallel shortlist or candidate database.
- Keep Results on their source workstation. Synchronize compact campaign/job metadata and retained Database records; do not replicate rejected populations or whole Results folders for UI convenience.
- Project/owner ACLs must be revalidated before coordinates, exports, analyses and annotations are accessed. Molecular equality does not grant access to another source group.
- A molecule, candidate occurrence, sequence variant, predicted model, docking pose and structure digest are distinct identities. Do not “fix” missing evidence by joining names, run labels, or sequence alone.
- Native/stage/final computational acceptance and human shortlisting are separate. Never promote an earlier pass to final acceptance.
- Preserve method definitions, units, versions, source/calculation digests and zero-versus-missing distinctions. No universal quality score.
- No automatic calculations when opening a page. Preserve queued/collecting work; retries require deliberate selection and reliable idempotency.
- Preserve Results → Scientific Inspector, portable `runs:///`/`cluster://` access, and the existing target-based overlay/contact map contracts.
- Do not alter historical Results, existing scientific observations, the SpiroSpace campaigns, or unrelated cleanup files. Any historical correction requires its own evidenced plan.
- Keep symlink-preserving backups and narrowly documented historical exceptions. AUD-003 changes future finalization, not historical execution environments.
- No deployment, push, service restart or production migration is implied by this document.

### Regression gates

Each finding's Verification field is its completion gate. In addition:

- Re-run the existing isolated suites listed below after relevant changes, plus new failing-before/passing-after regressions for the reproduced defects.
- Use two independent staged peers for review conflicts, structure aliases and Evaluation identity. Test offline replay and source access revocation.
- Browser-check a multi-engine collection and a >250-model synthetic cohort, actual JSON/CSV downloads, and exact Inspector navigation. Use existing authentic structures for a small final smoke test; no new design or prediction campaign is required.
- For digest checks, distinguish source files from evaluator-normalized outputs and handle a legitimate new revision as a new representation, not a retroactive digest update.
- For UI performance fixes, assert bounded query counts and no broad Results reads, not merely that widgets render.
- Record the fixing commit, exact test results, remaining limitations and resolution date under the existing finding ID.

## Audit coverage and evidence

### What was actually checked

| Subsystem | Work performed | Limitations |
|---|---|---|
| Repository/change baseline | Git status/history; relevant recent modules and tests; no root `AUDIT.md` or applicable ancestor `AGENTS.md` found | Not every historical commit reviewed |
| Database identity/ingestion | Read selected `MolecularDatabase` methods, Evaluation projection, exact-source adapters; isolated multi-model and file-integrity probes | No production-wide record reconciliation |
| Assessment/UI | Read review, collection, batch and export paths; existing AppTests and new isolated stale-form/collecting-state probes | Not every registered evaluator executed |
| Results/artifacts | Read candidate materialization and finalization validator/call sites; existing artifact tests and dangling-link probe | Not every generator's manifest coverage audited |
| Authorization/API | Read auth/API routes and file handle validation; isolated API tests; unauthenticated read-only primary `/auth/me` check | No penetration test, remote configuration audit, or credential changes |
| Peer replication | Read capture/apply and review conflict code; existing independent-peer tests executed | No live peer mutations or crash-injection campaign |
| Scientific lineage | Selected source/count/query inspection and lineage suite, including compact publication/retained links and offline metadata cases | No live backfill or verification of all historical campaign totals |
| Structures/MD/catalog navigation | Comparison, Inspector, catalog and relationship test suites; selected query/helper inspection | No fresh real Mol* GPU/browser overlay run; MD engines not executed |
| Production UI | Observed existing `8502` page without changing the user's 7-model selection; saw Tau additional-35 native-passed cohort, linked Biotite saved result/digest, shared developer login; captured page and console state | This was read-only observation, not a new end-to-end submission/export certification |

The production page was responsive to inspection. Its visible saved Biotite result was `20261103-224940-7fd7d198`, Evaluation `evaluation-a64bfcc95315cec7034fb481a3ff584c`, with a displayed exact structure digest. Browser error/warning log query returned no entries at that observation. These facts do not establish correctness of every metric or workflow.

### Executed checks

Python: `/home/user/mambaforge/envs/mn-molecular-design-studio/bin/python`.

All test invocations redirected Studio home, application home, Results, references and cluster home into `/tmp/mds-audit-20261111/`. GPU tests were disabled; the integration/navigation runs additionally set scheduler mode `local`. API and peer execution used the suites' isolated fixtures/mocks. No production calculation was submitted.

1. **161 passed in 16.31 s**:

   ```text
   tests/test_complex_analysis_design_review.py
   tests/test_complex_analysis_collection.py
   tests/test_complex_analysis_system_navigation.py
   tests/test_complex_analysis_database_hydration.py
   tests/test_artifact_collection.py
   tests/test_molecular_database_peer.py
   tests/test_scientific_lineage.py
   tests/test_evaluations_and_comparisons.py
   tests/test_complex_analysis_shortlist_export.py
   ```

2. **144 passed in 90.12 s; 4 warnings**:

   ```text
   tests/test_peer_workstations.py
   tests/test_api_integrations.py
   tests/test_api_extended.py
   tests/test_central_scheduler_worker.py
   tests/test_remote_results_evaluations.py
   ```

   The warnings were `pydssp` invalid-divide warnings in small validation-acceptance fixtures. No failed assertion was observed. An initial invocation also named a nonexistent `tests/test_complex_analysis_sources.py`; pytest exited before running tests. That path was removed for the successful invocation above. No tests were excluded with `-k` or marked away by this audit.

3. **41 passed in 4.01 s**:

   ```text
   tests/test_complex_analysis_comparison.py
   tests/test_complex_analysis_inspector.py
   tests/test_complex_analysis_catalog.py
   tests/test_complex_analysis_relationships.py
   ```

**Total: 346 existing tests passed.** Passing tests do not refute the additional reproduced boundary cases. The probes deliberately observed current incorrect behavior; they were not added as application tests or counted as passing regressions.

### Reproduction evidence retained locally

Temporary evidence directory: `/tmp/mds-audit-20261111/` (not production and not guaranteed durable across cleanup):

- `focused-tests.log`, `integration-tests.log`, `navigation-tests.log`: the runs above.
- `probes.py`, `probes.log`: disposable Database fixtures for AUD-001/002/003/004/007/008/009/010.
- `ui-probes.py`, `ui-probes.log`: real Streamlit AppTest plus temporary Database demonstrating AUD-005.
- `batch-probes.py`, `batch-probes.log`: applicability and stub-submission AppTest for AUD-006/AUD-012.
- `production-ui.txt`, `production-ui.png`: read-only production UI observation.

The finding descriptions provide the reproduction conditions and observed outputs independently of these temporary files. To rerun probes, set all Studio/cluster roots to a fresh temporary directory before importing application modules; never run them against the production Database. Temporary form tests emitted expected bare-run/deprecated-Streamlit warnings, without exceptions.

## Systemic patterns and remaining limits

- **Identity fallback is too permissive in some low-level paths**, despite stronger checks in newer callers. Fix boundaries rather than adding more caller-specific exceptions.
- **Multiple partial definitions of compatibility and lifecycle coexist.** This explains version-blind reuse, collecting-state duplication, and Tag/Fusion disagreement.
- **The lightweight UI projection is not yet a complete contract.** Its pagination, saved-group resolution and export assumptions differ from the hydrated workspace.
- **Peer conflict handling does not substitute for form-level concurrency control.** The stale human editor needs an explicit saved revision.
- **The test suite is valuable but dominated by small happy-path cohorts at these boundaries.** Add model ambiguity, multi-page, changed-byte, stale-form and old-job cases; do not simply add more copies of the existing success tests.

Unverified: full production historical data correctness; all scientific numerical definitions; real container/tool availability; long-running scheduler recovery; storage/backup restoration; every workflow and UI page; all remote services; browser download behavior after proposed corrections; independent live Mol* overlay in this audit; adversarial security beyond the stated observations. Existing deployment/backup investigations were not repeated. No unsupported claim of data loss, public network exposure, or broken replication is made.

**Recommended next step:** Sol should implement the narrow fixes in the order above, starting with the five High findings, and append verification evidence to this document. A redesign or Database replacement is not justified by the evidence collected in this pass.
