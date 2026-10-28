# Native execution-job ledger — Phase 6C.2

## Result and scope

Phase 6C.2 adds durable execution jobs to new campaign ledgers, independently of
scientific attempts. Deterministic fixtures can represent several independently
accounted attempts in one job, including partial results, missing receipts,
operational retries and late validation. **Production generation still reserves
and executes exactly one native attempt per job.** Multi-member reservations use
an accounting-only API and are refused by campaign dispatch and reconciliation.
No native batching capacity, persistent worker or new engine adapter is enabled.

Baseline: `f4586b1` and `41ff81b`. The design follows sections 9–10 and 14 of
[native-execution-batching-architecture.md](native-execution-batching-architecture.md)
and preserves the measurement limitations in sections 17–19. No protein generation,
validation inference, Docker command, peer deployment or GitHub push was performed
for this phase. The Phase 6C.1 timing run was not repeated.

## Implemented contract

New ledgers use campaign `schema_version=2`. Existing schema-v1 files are read and
recovered in place with their original interpretation; opening them does not
migrate, rewrite or synthesize job membership. Unknown campaign versions are refused.
`CampaignLedger.create(..., schema_version=1)` remains available for legacy fixtures.
The creation-version argument applies only to new files.

`execution_job_ledger.py` implements the accounting contract. All mutations use the
campaign's existing exclusive file lock, temporary file, file fsync, atomic replace
and directory fsync. Job and attempt authorization are one transaction. A failed
reservation writes neither a partial manifest nor an incremented global ordinal.
Exceptions and conflicts leave the previous ledger intact.

| Record | Durable fields and interpretation |
| --- | --- |
| Execution job | Nested `record_version=1`; `execution_job_id`; reservation `request_key` and digest; immutable ordered `manifest` and digest; immutable `configuration` and digest; authorization time; adapter mode; execution history/current execution; verified `output_receipts` |
| Configuration | Normalized campaign settings, declared native context and requested resource requirements. Production dispatch records the parent input payload; its bytes' SHA-256 and parent run reference identify the declared source input. Resolved scheduler requests and actual allocations are observations on an execution. |
| Manifest member / attempt | Global `attempt_id` and ordinal; zero-based job position; requested seed and explicit RNG mapping declaration; source input reference; job membership; authorization time; legacy worker state; separate `native_status` and tri-state `native_started`; result execution/evidence, output IDs, receipts and candidate IDs |
| Operational execution | Stable `execution_id`; retry key; stable accounting submission key; submission/launch/contact states; scheduler identifiers including child and peer IDs; worker identity; resource allocation; resolved resource request when observed; lease-duration receipts; terminal evidence when supplied |
| Attempt receipt | Job/execution/attempt identity, manifest digest, outcome, native state, native output identifiers, complete returned candidate payloads, verification evidence and canonical content digest |
| Candidate | Attempt/job/execution/receipt ownership, native identifiers, original sequence/structure/ancestry payload, separate native/validation/diversity acceptance and selection disposition |

A job identity is `campaign:j<first-global-ordinal>`; an attempt remains
`campaign:a<global-ordinal>`. Executions append `:e0001`, `:e0002`, and so on.
Identities are stable once reserved. Historical `batch_id` remains **dispatch
grouping**. It is neither used as a job ID nor reconstructed into native membership.
Several new single-attempt jobs can share the same dispatch group.

Digests are SHA-256 over sorted, compact JSON with nonfinite numbers refused.
Membership/configuration digests and the correspondence between manifest entries
and attempt records are checked on reads and writes. Worker transitions cannot
change authorization, seed, source input, RNG declaration or membership. Candidate
updates cannot replace verified ownership, native acceptance or scientific payloads.
These checks detect contract conflicts; they are not signatures against an actor
who can rewrite the entire file and recompute its digests.

The seed derivation is unchanged from Phase 6B. A manifest records independent
**requested** seeds with `native_mapping_verified=false`. This declaration does not
claim that one future batched invocation implements independent per-attempt seeds.
Actual native state and existing candidate provenance remain attached to verified
receipts. Replicate/RNG equivalence within a multi-attempt invocation is still a
prerequisite for enabling such an adapter.

## Authorization and submission

`reserve_execution_job(request_key, source_inputs, native_configuration=...,
resource_requirements=...)` atomically authorizes the entire ordered manifest.
It refuses a manifest exceeding the remaining global budget; it does not silently
truncate a requested size. Repeating a key with the same declaration returns the
same job. A changed declaration under that key raises `LedgerConflict`, including
after controller reconstruction. Cancelled or not-started reservations still consume
their authorization; budget is never silently recycled.

This explicit API always creates `mode=accounting_only`. It performs no native,
scheduler or network action. `prepare_execution_submission()` persists accounting
intent and a stable submission key, but likewise performs no submission.
`observe_execution()` binds observations to that execution/key. A different
scheduler/worker identity under an existing binding is a conflict.

The existing `authorize(n, batch_id=...)` API keeps its legacy budget-clamping
behavior. In v2 it creates **n separate single-member jobs** atomically, rather than
one n-member native job. Production still calls it with `n=1`. The controller,
dispatcher and child creators refuse accounting-only ledgers before launching work.

The new accounting submission key is not substituted into peer protocols. Existing
single-attempt remote submission keys, child tags, scheduler calls, authenticated
placement, image/checkpoint compatibility and transfer verification retain their
existing roles and are recorded as observations. The ledger is not an exactly-once
scheduler: after an intent/acknowledgement gap, a future multi-attempt adapter must
query the same external identity before considering any new action. Tests exercise
durable identity and conflicting acknowledgements without dispatching.

## States and recovery

Submission states are `reserved`, `submission_pending`, `queued`, `running`,
`completed`, `failed`, `interrupted` and `cancelled`. Launch knowledge is separate:
`not_exposed`, `may_have_started`, or `not_started_confirmed`. Even a submission
intent conservatively changes launch knowledge to `may_have_started` because the
external side effect could occur before the next durable acknowledgement.
Contact state is independently `unknown`, `connected` or `disconnected`.

Native attempt states are `authorized`, `running`, `completed`, `failed`,
`not_started`, `cancelled` and `unresolved`. The old worker-facing `status` field
remains available. A worker PID/running observation does not prove that a particular
native attempt started. `native_started=true` requires explicit native evidence or
a verified per-attempt result. Missing post-launch evidence remains unknown;
it is not converted to `false`.

| Boundary / event | Implemented accounting behavior |
| --- | --- |
| Crash before atomic authorization replacement | Old file remains valid; no partial members or budget debit. A subsequent reservation starts at the original ordinal. |
| Crash after reservation, before submission | Job, ordered members, seeds and configuration are recoverable by the reservation key. No new authorization is needed. |
| Duplicate submission observation | Same execution/key and scheduler identity is idempotent; conflicting identity fails atomically. No submission method dispatches work. |
| Worker/peer disconnect and reconnect | Explicit accounting observations retain the same job/execution/membership; contact status changes without reserving or resubmitting. Existing peer transport is not replaced. |
| Partial completion | Each verified attempt result is retained immediately, with all candidates. Other members remain independently pending/running. |
| Process exit, including exit zero | Terminal execution evidence settles the process, not missing native attempts. Missing per-attempt results become `unresolved`, regardless of process success. |
| Native rejection | A verified completed attempt may return zero candidates and authoritative rejected-trajectory evidence. Missing files or an incomplete process are not rejection evidence. |
| Interrupted execution | Complete verified members remain complete; remaining members become unresolved. Previously known native-start evidence is retained. |
| Late verified receipt | An unresolved member can resolve after process termination; unrelated completed members and their candidates are unchanged. |
| Cancel before submission | Reservation becomes cancelled with native-start false. Repeated cancellation is idempotent. Authorization remains consumed. |
| Cancel after submission / while queued | Record a cancellation request. Only explicit no-native-launch terminal evidence can mark members cancelled/not started; a queue label alone is insufficient. |
| Retry before native launch | Explicit settled no-launch evidence permits a new operational execution under the same job, manifest, attempt IDs and seeds. Retry keys are idempotent. No scientific authorization is added. |
| Ambiguous or post-launch retry | Refused. Partial results, known native starts and ambiguous launch state prohibit replay. A same-seed rerun is not treated as proof of physical-attempt preservation. |
| Budget exhausted, validation pending | Generation exhaustion is distinct from finalization. Candidates and validation remain reconcilable; late accepted validation can promote a budget outcome to target reached. |

Native-result and terminal-execution receipts have strict identity/content checks.
An identical result under an alternate receipt ID is an idempotent no-op; it adds
neither a new candidate nor a second result record. Changed content, candidate ID
collisions across attempts, foreign membership, wrong manifests, conflicting
terminal evidence and results from superseded executions require explicit
reconciliation. Confirmed no-launch evidence cannot be erased by later native results.
The API trusts a verified adapter's evidence object; it does not itself open or
validate arbitrary scientific artifact paths supplied by a caller.

Production mirrors existing child/peer observations into each single-member job.
The existing native candidate verification remains the receipt boundary. Native
output IDs come from `raw_metadata.native_candidate_id`, rather than relabeling
global normalized candidate IDs as native IDs. When a legacy API fixture supplies
no native identifier, the native-ID list is empty; candidate identity is still retained.

RFD4's existing `command.json` retry guard remains in force. Controller recovery now
also requires pre-native-launch evidence for other engines: an unavailable child or
recorded Docker/native invocation cannot justify automatic stochastic replay.
This conservative automatic-recovery check applies to old ledgers too, without
changing their stored interpretation. The legacy v1 explicit `retry()` API retains
its original record format. V2 prelaunch retries create a distinct execution and
persist submission state before resuming the existing scheduler job. Their guard
record includes the child and command-record hash when available. A recovered
unlinked child with native launch evidence is not automatically relaunched.

## Scientific and compute accounting

Old counters keep their original meaning. In particular, `attempts_started` is a
worker-start observation and `failed_interrupted_attempts` includes operational
retry history; neither is relabeled as an exact native count.

V2 adds total authorized jobs, submitted-job intents, execution retries, unsuccessful
execution events, total native authorizations, pending native attempts, known native
starts, separate native-state counts, selected/returned-unselected candidates and
pending validation candidates. `execution_job_summary()` reports the receipt/state
frontier per job independently of process submission state. A submitted-job count
reflects durable intent, not proof that an external process launched.

Multiple candidates from one receipt remain one completed attempt. Explicit
`native_pass_filters=true/false` can establish native acceptance/rejection for an
engine with authoritative acceptance; missing information is unavailable. RFD4
native acceptance stays unavailable even if a process succeeded or a fixture
includes a true filter flag. A trajectory-rejection count requires explicit native
trajectory and acceptance counts; absent acceptance data is not inferred as rejection.

Validation and diversity acceptance are separate decisions. Selection supports
`pending`, `selected` and `returned_unselected`, with reasons through candidate
updates. Every returned candidate payload is retained, including surplus outputs.
The contract does not impose a new selection policy on existing campaigns.

Native GPU-duration receipts belong to an **execution**, under its job. They are
not copied onto every member. Replaying the same lease duration is a no-op; a
changed duration conflicts. The explicit job API refuses one lease identity owned
by two executions. Production bridges its existing child-run duration observations
to the current single-member execution. Validation durations remain independent
attempt/candidate-work records, explicitly tagged by the validation call sites.
Totals sum native execution receipts once plus validation receipts. Existing live
compute estimation and scheduler resource placement remain on the single-attempt
production path. Durations are leased-worker estimates, not GPU utilization or
native lifecycle timings.

## Read-only execution planner

The existing pure planner adds `execution_job_manifest_preview`. It displays proposed
job IDs, ordered attempt IDs/ordinals, job positions, requested seeds, input membership,
resource requirements and manifest digests. It shares the pure membership/seed
builder with durable reservations. Its synthetic namespace is derived from planning
inputs; preview IDs and seeds are not the IDs/seeds of a submitted campaign.

For example, global budget 13 across two inputs with requested native size 5 produces
hypothetical memberships of 5/2 for input 0 and 5/1 for input 1, exactly 13 attempts.
The result remains unsupported, with no supported execution plan. Engine-specific
rounding/output multiplicity remains in the existing planner fields; a membership
preview does not authorize rounded extra native samples.

Previews are bounded to 16 job examples and 100 member examples. Exact aggregate
counts, omitted-job/member counts and truncation flags remain visible. A partial
job has no complete manifest digest. Native configuration/input identity and native
RNG mapping are explicitly unresolved/unverified. The preview's configuration digest
covers its supplied campaign/resources and empty native context; it is not a resolved
native invocation digest. All previews have `authorization=false`, and the parent
plan retains `submission_authorized=false`.

Choices 1/5/10/25/50/100 remain descriptive. Supported campaign choices remain `[1]`
for BindCraft2/RFD4 and empty for other engines. The existing form shows the additive
JSON through its read-only preview; saved pipeline parameters and execution controls
are unchanged. Planning creates no ledger, reserves no attempt, probes no device and
calls no create/enqueue/dispatch operation.

## Verification and retained evidence

The new deterministic tests use fixed ledger clocks, fixed campaign IDs/seeds and
synthetic records. They cover atomic whole-manifest authorization across controller
instances, failure before atomic replacement, reservation replay/conflict, submission
intent/duplicate acknowledgement, foreign or malformed membership, configuration
and order tampering, partial and late results, missing/conflicting receipts, failed
native attempts versus worker failures, exit-zero ambiguity, interrupted execution,
restart, disconnect/reconnect, queued cancellation, prelaunch retries, post-launch
refusal, single-job lease accounting, candidate lineage, surplus selection, pending
validation/finalization, legacy reading, production guards and non-executing bounded
planner arithmetic. Existing mocked controller/peer tests remain in the sweep.

Focused regression command:

```text
/home/user/mambaforge/envs/mn-molecular-design-studio/bin/python -m pytest -q \
  tests/test_execution_job_ledger.py tests/test_native_execution_contract.py \
  tests/test_execution_timing.py tests/test_generator_capabilities.py \
  tests/test_distributed_campaign.py tests/test_rfd4_proteina.py \
  tests/test_rfdiffusion3.py tests/test_native_generators.py \
  tests/test_bindcraft_boltzgen_contracts.py tests/test_distributed_campaign_ui.py \
  tests/test_candidate_peer_transfer.py tests/test_workflow_queue_lifecycle.py \
  tests/test_execution_planning.py tests/test_rfd4_native_diagnostics.py
```

Final sweep: **366 tests passed in 21.40 s**, including 45 new job-ledger cases.
Implementation/test commit: `2fb4cdc` (`feat: add durable native execution-job and
attempt accounting`). The report and retained evidence are committed separately.
Changed Python files also pass syntax compilation and `git diff --check`.
No full broad suite or GPU test suite was run.

[Retained checks](evidence/phase6c2-ledger-checks.json) record read-only inspection of
the original P4/P5 ledgers, normalized indexes and referenced scientific artifacts.
Recovery/finalization ran only on temporary ledger copies. Original schema-v1
configuration normalization remains equal; reads preserve original contents and
job membership remains unavailable. Finalization on copies preserves scientific
attempts, candidates and counters, including the existing late-validation promotion.
Original historical termination reasons are not rewritten.

| Retained record | Ledger SHA-256, before and after | Candidate validation identities |
| --- | --- | --- |
| P4 `20261028-101057-a07ef0a4` | `0f575ca3e563181423524207c34ebc5fbcf1c6c9b6ff8c5bf71409ae10bc1460` | One same-sequence validation, identity pass |
| P5 `20261028-124448-350f7294` | `205efddd4c18cea8f8ca9d9ee22a9b5ad193964f5e059854631bf2465e891b39` | Two same-sequence validations, identity pass |

All 16 checked ledger/index/structure/prediction file observations retain their
hashes. Existing saved-pipeline settings, deployed source trees and peer services
were not changed. Unrelated edits in `run_app.py` and untracked `tool_references.py`
are excluded from both commits.

## Requirements before enabling multiple native attempts

This phase verifies the accounting contract with deterministic evidence. It does
not verify a live multi-member native adapter, distributed partial artifact recovery,
or a batched scheduler/deployment path. Remaining gates include:

1. Verify engine-specific native input/replicate/trajectory-to-global-attempt and RNG
   mapping, including actual native starts, incremental completion and every output's
   sequence/structure ancestry. The contract's requested seeds alone are insufficient.
2. Implement and test a trusted incremental receipt producer with artifact hashes,
   completeness checks and recovery of partial files; distinguish explicit no-launch
   proof from missing progress. Existing single-attempt scientific verification is
   not a multi-member receipt adapter.
3. Verify the real scheduler's idempotent submission/query boundary, authenticated
   peer reconnect/artifact recovery and lease release for a multi-member job. The
   pure accounting methods do not supply exactly-once external execution, terminate
   workers, or perform cross-origin controller takeover.
4. Replace one-active-child assumptions in dispatch headroom, target backpressure,
   cancellation and live compute estimation with manifest-aware behavior before any
   multi-member job is allowed into those paths. Exact remainder sizes, output
   multiplicity and validation draining need adapter-level verification.
5. Separately verify scientific equivalence, cancellation/resume boundaries and
   engine/hardware limits for each proposed size. Keep ambiguous post-launch replay
   disabled. Flat-file ledger scalability has not been benchmarked for large manifests.

Phase 6C.1 measured per-invocation import/model work on one bounded local RFD4 run;
it did not measure multi-attempt reuse savings or throughput. Its uncertainties
remain unchanged. This phase makes no new performance claim and does not use CPU
container timing as evidence of model-reuse speedup.
