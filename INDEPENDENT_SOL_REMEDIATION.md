# Post-audit corrections and release gate

Date: 2026-11-11. Audited baseline: `e7cd3a95c5163f6e973187a24aab002f5e5c351b`.

The independent assessment and its original evidence remain unchanged in
`INDEPENDENT_SOL_AUDIT.md` and `INDEPENDENT_SOL_AUDIT_EVIDENCE.tar.gz`.
This document describes subsequent implementation, not an independent
certification of that implementation. The user subsequently authorized fixes,
commit/push and distribution to the two `compute` workstations, choosing to fix
release blockers before deployment rather than activate the audited build as-is.

## Implemented corrections

| Finding | Correction and scope |
|---|---|
| SOL-001 | Database Inspector defaults to registered protocol checking; its rendering gate refuses an unverified completed protocol. Results Biotite/contextual matching checks recorded protocols rather than dropping versions. Results searches filter exact upstream identities in SQL before paging, without a first-1,000 limit. Saved image/checkpoint digests without a verified expected counterpart do not become automatic compatibility claims. |
| SOL-002 | One shared effective-settings projection supplies native interface-chemistry fixed defaults and separates subject/model/digest selectors from scientific cutoffs. The real calculate/ingest/review/index round trip is a permanent regression. |
| SOL-003 | Binding comparison gates include target structure revision and target-chain scope, including known/unknown differences; sequence-only fields retain their existing exemptions. |
| SOL-004 | Evaluation signatures use the shared scientific-settings projection rather than private file paths and candidate selectors. Protocol provenance remains distinct. Derived projections advance to `evaluation.v7`; canonical observation and Evaluation IDs are not rewritten. Review ignores stale cached descriptor versions and projects authorized observations afresh. |
| SOL-005 | Exact Database assessments have durable, owner-local admission in the existing request store, independent of transport request keys. File locks cover concurrent processes; project sessions share scientific admission while private owners remain separate. Pending/uncertain outcomes are not blindly repeated. Existing indexed work is checked under the admission lock. Explicit `execution.reassessment_of` identifies the exact previous job, is retry-idempotent, and cannot duplicate active work. |
| SOL-006 | Receptor visualization uses the same pose-coordinate-frame verifier as analysis and serves a verified snapshot. Recorded Results manifest identities are enforced for downloads, Inspector entry and analysis path resolution; changed recorded coordinates are rejected rather than assigned a replacement identity. Unversioned historical artifacts are not falsely described as verified. |
| SOL-007 | Numeric-zero model selection survives descriptor, metadata review, Database source authorization, Inspector and shortlist-export projections. Unknown values remain unknown. |
| SOL-008 | Assessment planning, fresh pre-submit lookup and Inspector handle actual SQLite exceptions without interpreting failure as missing work. Submission-state storage failures return HTTP 503; uncertain submission state is retained. |
| AUD-007 | Molecule-only curated membership can offer accessible source models without copying scientific structures. Groups that own their own source models remain restricted to those models rather than silently acquiring another campaign's representations. This does not certify every historical frozen comparison collection. |
| AUD-008 | Stable model selection survives nonempty page changes; off-page models are re-read by exact IDs and reauthorized against source membership and selected system before actions. Cohort/search/measurement filtering still follows the bounded model page; see remaining limitations. |
| AUD-009 | Metadata review includes accessible saved source provenance; shortlist sequence maps use recorded chain components without coordinate hydration. |
| AUD-010 | CSV text cells beginning with spreadsheet formula prefixes after leading space/control whitespace receive a literal apostrophe. Native numeric values and JSON exports retain their original semantics. |

### Concurrent input preparation

Preparation is serialized per account, source bundles are qualified by workflow,
and newly admitted work receives a separate sealed input copy. JSON/JSONL source
references are rewritten to that copy and coordinate SHA-256 is rechecked before
queue admission. Later preparation cannot replace the queued request's source
manifest or redirect its coordinate file. Duplicate requests reuse the durable
response rather than creating another input copy or queue job.

The added regressions cover real concurrent preparation and admission, different
method bundles, later source-cache changes, effective-default and numeric-format
equivalence, project/private admission scope, crash uncertainty, and explicit
failed/completed-job reassessment.

One additional admission regression was reproduced during final hardening:
invalid resource settings returned 422 but reserved the scientific key, so a
corrected request received 409 without ever having queued work. API schema,
resource and source preflight now run under the preparation lock **before**
scientific reservation. True post-admission uncertainty remains conservative.

### Distributed admission policy

A replicated Database structure remains viewable through its authorized source
adapter, but API submission on the non-owning workstation is explicitly refused.
Submit on the source owner. This is a deliberate fail-closed fence, not a claim
that a new distributed reservation or remote Database-dispatch protocol exists.
All participating writers must be upgraded together. An older binary can bypass
the new admission and review controls; production mixed-version writing remains
unsupported during cutover.

## Regression contract alignment

The broad audit run exposed 14 failures, not 14 independently established engine
defects. Maintained tests now use the current explicit passing-collection UI,
Database-first navigation, versioned derived reprojection, verified snapshot
bytes, recorded-path resolution and initialized job lookup. Successful GROMACS,
RFD4 and shared-pipeline fixtures now retain their declared scientific artifacts;
strict finalization is not weakened. The starter-template contract includes the
existing fourteenth template. Legacy ranking without a scientific criterion no
longer invents a pass solely from the default output limit; explicit metric/rule
selection remains separate.

## Validation

- The maintained adversarial release module passes **40 tests**, including input
  sealing and invalid-resource admission recovery.
- Earlier focused runs passed 160 and 72 relevant cases; those are supporting
  intermediate results, not substitutes for the final broad run.
- Broad suite: **3,077 passed, 37 skipped, 103 warnings in 943.82 seconds**.
  All 3,006 baseline passing cases still pass. Eleven original failures pass under
  the same case names; three obsolete navigation tests were replaced by explicit
  Database-first navigation contracts. Thirty-nine new release cases ran here.
  Eighteen image-adapter launcher tests that skipped in the detached audit worktree
  passed because their external source fixture exists beside the primary checkout;
  these are deterministic launch fixtures, not live image/GPU certification.
- Final API/planning/index/peer integration run: **248 passed in 95.06 seconds**,
  including all 40 release regressions. This run validates the last API preflight
  guard and its additional regression; the broad run started before that guard.
- Tests run with `scripts/test_release_isolated.sh`: read-only host root, hidden
  production Results mount, private application/cluster homes and disabled network.
- GPU/container computation remains opt-in and was not enabled. No dependency
  upgrades, authentication changes or production scientific writes are part of
  these corrections.

Reproduce with an existing Studio environment:

```bash
MDS_TEST_PYTHON=/path/to/studio/python bash scripts/test_release_isolated.sh -q
```

The launcher requires Linux Bubblewrap. Interrupted broad development runs are
excluded from final totals. Retain final JUnit/logs separately from the original
audit evidence.

Final local JUnit/logs, projection-rehearsal statistics, source checksums,
non-sensitive cleanup counts and per-host staging receipts are packaged in
`INDEPENDENT_SOL_REMEDIATION_EVIDENCE.tar.gz` (SHA-256
`12a55c559bd40ae266f02eca676d4409b45997f6ceabf32c759dd94adfa0bed5`).
The package excludes scientific Databases, raw structures, complete private
cleanup inventories and credentials. It includes the blocked remote launcher
logs rather than claiming remote regression passes.

### Actual-data projection migration rehearsal

A consistent SQLite backup was made from a read-only connection to the primary
production Molecular Database, then tested only in a private, network-disabled
scratch copy. **41,213** derived Evaluations advanced from v6 to v7 in **83**
batches, in 82.54 seconds. Every non-derived table retained its exact row-stream
digest, including scientific records and review state; canonical Evaluation IDs
were unchanged and SQLite `integrity_check` returned `ok`. No scientific engine
was invoked. This validates a primary Database projection upgrade, not a full
three-host recovery rehearsal or correctness of historical scientific links.
The private 1.7 GB Database copy is not included in Git or public evidence archives.

## Activation gate: not cleared

No production service is restarted, no scheduler is stopped, and no scientific
queue or data is edited by this release preparation.

Read-only checks found the running source on all three hosts is still
`p11-25e003c/source`. At the recorded preflight, MN Compute Scheduler reports no
queued/running/paused jobs for the Studio app on any host. This is **not** proof
that every legacy Studio queue is drained: the persistent Results catalog still
records 1 queued, 16 running and 19 collecting entries on pc002975, and 3 collecting
entries on PC004033. Those records were not cancelled, repaired or relabelled.
Two long-lived local-worker processes were also inspected: their configured
application/run roots are isolated scratch previews, not the production Results
root. They were left alone and are not represented as active production jobs.

The user noted that old jobs may be test artifacts. A read-only inspection of
their metadata/request/result files found no obvious fixture markers and no
`result.json` for all 36 primary and 3 PC004033 nonterminal records. That does not
prove either scientific importance or live execution; their disposition remains
unresolved. Absence from the central scheduler is not permission to delete them.

The user subsequently suggested considering workdir/Database cleanup. A private,
read-only primary-copy inventory found 1,050 molecular groups and four explicit
import batches. Its 15 lineage campaign records are not a complete inventory of
older design campaigns: retained campaign data also exist in molecular groups.
All inventory items remain `UNREVIEWED_PRESERVE`. Cleanup is a separate operation,
requiring user-confirmed protected IDs, a cross-host source/target/Evaluation/file
dependency closure, restore-tested recovery sets, and an exact reviewed quarantine
manifest. No workdir was emptied and no Database rows, replica cursors, epochs or
tombstones were reset. A clean release must not be conflated with a data reset.

Do not activate this candidate until an operator has established the disposition
of those exact jobs, closed admission/peer writers in one maintenance window,
verified current per-host recovery sets and isolated restores, and verified
uniform API/agent/UI/worker provenance. Follow the existing
`docs/database-simplification-p0-deployment.md` procedure. The previous release's
backup receipt is not proof of a current quiet-point backup for this candidate.
The additive job projection advances to `database_analysis:v2`; its ordinary
reconciler must finish before upgraded analysis admission is reopened.

Read-only environment inspection found Python 3.11.17, RDKit 2026.9.1, NumPy
2.4.6, SciPy 1.17.1, DuckDB 1.5.6, Pydantic 2.14.1 and Biopython 1.88 on all
three hosts. The primary has FastAPI 0.143.0 / Streamlit 1.65.0; both compute
hosts have FastAPI 0.143.3 / Streamlit 1.66.0. No packages were changed to mask
that difference; isolated remote regression results must be recorded separately.

Authentication remains unchanged. Local production services explicitly use
development authentication; this candidate does not turn that shared identity
into per-researcher authorization. Keep the existing loopback restrictions and
plan the separately reviewed account/ownership transition.

## Remaining limitations

- AUD-008's pre-pagination cohort/search/measurement filtering is not fully
  redesigned; visible-page counts must not be described as whole-cohort totals.
- Historical invalid scientific links and unknown identities are not repaired.
  Frozen collection semantics and legacy caller projections need further review.
- Unknown container/checkpoint equivalence is conservative, potentially requiring
  explicit reassessment rather than automatic reuse. No live container parity
  certification or GPU scientific validation is claimed.
- Atomic scientific admission covers exact Database assessments on their source
  owner, not every generic Results/campaign/pipeline request.
- Updated CAS writers remain required on every review-writing process; an old
  blind writer is not made safe by this release's local generation triggers.
- Sealed inputs and verified snapshots are required application state. Include
  them in recovery backups; do not delete them while queued or retained jobs
  reference them.

Commit identity and per-host source/checksum receipts belong outside the immutable
source artifact. Staging the same candidate on `compute@pc004032.marionegri.it`
and `compute@pc004033.marionegri.it` is distinct from activating their services.

## Source distribution and remote validation

Application correction commit: `8dd52dbf502198146d2e4fdc975da2be1279f36a`,
pushed to `origin/codex/p0-production-visibility`. The same tracked source archive
is staged on the primary and both `compute` hosts under
`~/programs/deployments/mn-molecular-design-studio/sol-8dd52dbf5021/source`.
Its SHA-256 is
`a1d85ef6de968d50b4f450179cf788953297dd886dc883418fddfaa7d101cb17`.
Every application/test/launcher file also passed its per-file SHA-256 check on
each host. An incremental Git bundle imported the pinned commit without changing
the compute workstations' active `main` checkouts.

Remote isolated regression launches on **both** compute hosts stopped before
pytest collection with `bwrap: loopback: Failed RTM_NEWADDR: Operation not permitted`.
An explicit user-namespace probe also failed; filesystem-only namespace setup
was denied. These are unavailable remote test runs, not passing regressions or
evidence of an application failure. No host security policy was weakened and no
unisolated fallback test was run. Runtime parity on those hosts remains an
activation blocker. Staging receipts retain the nonzero launcher exit status.

Service PIDs, active states and source directories were identical before and
after distribution. All active production APIs/agents still use the previous
`p11-25e003c/source`; the primary UI is unchanged and peer UIs remain inactive.
The source artifact is a **staged candidate, not a production-ready activation**.
This later documentation/evidence update does not change its application code.

### Cleanup dry run

Using the user's suggested smoke/testing-name and single-entry heuristics, the
primary inventory flags 284 name matches and 223 tool-output groups with exactly
one member and at most one structure. These sets overlap. Group membership is
not a verified whole-campaign result count: a one-result child may belong to a
valuable multi-job campaign, and shared molecules may occur in several groups.
No flag is deletion approval.

The primary union is 436 flagged groups, not 507 independent campaigns. All
other data remain preserved too; lack of a flag is not an importance assessment.

Separate read-only transactions inventoried both peer Databases. They each have
1,050 groups and 15 lineage campaign records at capture, but no local import-batch
rows; the primary has four. Host-specific import history and replica state must
not be treated as interchangeable. Inventories are private local review artifacts,
not new scientific imports; all records remain `UNREVIEWED_PRESERVE`. Counts across
replicas must not be added together to estimate unique scientific data.

## Subsequent student-readiness acceptance

The user authorized proceeding with recovery, isolated teaching and student
acceptance. `STUDENT_READINESS.md` records that later work and its current GO/NO-GO
boundary. New acceptance found and corrected default admin request staging in
immutable deployments and post-ingestion descriptive-role changes in scientific
admission identity. Admission v2 includes conservative v1 reservation reconciliation.
Remote container isolation now permits actual peer regression execution without
changing the hosts' security policies. A separate, authenticated, bounded teaching
pilot is running; research production remains unchanged and is not activated.
The original independent report remains the baseline assessment, not a retrospective
independent certification of these additional implementation changes.
