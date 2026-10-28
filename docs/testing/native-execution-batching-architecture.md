# Phase 6C — Native execution batching, Docker lifecycle and campaign scheduling

Date: 2026-10-28. Starting branch: `fix/scientific-pipeline-validation`, HEAD
`a25a570` (after `726a9a4`). Evidence labels in this report distinguish **measured**,
**verified in source**, and **proposed**. Source support is not throughput evidence.

Sections 1–16 preserve the Phase 6C report. Sections 17–20 record Phase 6C.1,
starting from `d8a1b69` (including `d88e8c1`), and supersede the RFD4 timing
unknowns where measurements are now available. Phase 6C.1 ran one bounded native
PD-L1 request; the earlier statement that Phase 6C ran no inference remains historical.

## 1. Result and scope

**Verified in source:** the distributed BindCraft2 and RFD4 adapters create one
scheduler child, container and native trajectory/request per authorized attempt.
Their ordinary native engines can reuse model objects for multiple attempts in a
process. RFdiffusion and RFdiffusion3 also reuse initialized sampler/trainer objects
within an invocation. BindCraft2 reuse is per native worker, not necessarily per
container. Native vectorized batching and sequential process reuse are different.

**Measured:** 24 CPU-only no-op container launches on `pc002975` took 0.174–0.332 s
across the four installed engine images. P5's retained RFD4 worker intervals remain
39 s and 51 s. The no-op runs exclude GPU runtime initialization, scientific Python
imports, checkpoint loading, compilation and inference. They are not an estimate
of those phases inside the P5 workers. Docker launch alone is inexpensive in this
local experiment; a meaningful avoidable model/startup cost in production has
**not been quantified**. No speedup or preferred native batch size is established.

This phase adds host lifecycle timing, descriptive capabilities, compatibility
validation and clearer UI explanations. It leaves campaign native execution at
one attempt per job. It does not add persistent inference workers, enable larger
batches, change scientific thresholds or targets, deploy source, or push commits.
No generation or validation inference was rerun. PD-L1 retained evidence supplied
the scientific compatibility checks; CPU-only containers supplied the new timing
baseline. A larger inference comparison would not resolve per-phase costs without
native phase instrumentation, and is deferred rather than inferred from totals.

## 2. Actual revisions and evidence provenance

Read-only SSH inspection of both peers, using the documented `compute` account,
confirmed their checkout HEADs at `4f9ba772fcb48c1961d74c99f2f5704f34c06197`.
They are not at local `a25a570`. Each peer has an unrelated modified `run_app.py`
and untracked `tool_references.py`, as does the local checkout; these were left
untouched. No service restart, checkout update or deployment occurred.

| Inspected local image | Docker image ID (prefix) | Native files inspected |
| --- | --- | --- |
| `mn-bindcraft2:v1.0.4-cu13-pae-v1` | `673ec7cf961e` | `/opt/bindcraft/bindcraft/{campaign,design_workers,campaign_output,trajectory,af2}.py` |
| `mn-rfd4-proteina:0.0.3-cu130` | `c34efb5afca7` | `/opt/rfd4/src/rfproteina/{generate,proteina}.py` and supporting modules |
| `mn-rfdiffusion:latest` | `8b40d74ca019` | `/opt/RFdiffusion/scripts/run_inference.py`, `rfdiffusion/inference/model_runners.py` |
| `mn-foundry:cu128` | `bd959dfe52dc` | `/app/foundry/models/rfd3/src/rfd3/{engine,run_inference}.py`, `/app/foundry/src/foundry/inference_engines/base.py` |

The active RFD3 tile selects `mn-foundry:cu128`, not the installed `:cu128-batch2`
tag. Both were checked; their inspected RFD3 engine files are identical. RFdiffusion's
mutable `latest` tag does not provide a trustworthy upstream revision; the image ID
and file checksums identify this inspection.

Both peers' RFD4 image IDs are `9ffb161fd399…`, different from the local image ID.
Direct SHA-256 checks inside those images matched the local `generate.py` and
`proteina.py` files exactly. This supports the lifecycle findings for these inspected
files, not complete binary/image equivalence. It does not replace the authenticated
scheduler's image-content fingerprint and checkpoint compatibility requirements.

Raw launch samples, full image IDs, native source checksums and peer revision
observations are committed in [container evidence](evidence/phase6c-container-startup.json).
Read-only retained P4/P5 config, index, same-sequence and identity checks are in
[retained campaign evidence](evidence/phase6c-retained-campaign-checks.json).
Local source snapshots and the no-op measurement script are retained under
`/mnt/data/RESULTS/mn-molecular-design-studio/verification/phase6c-native-batching/`.
Third-party source and weights are not copied into this repository.

## 3. Current application execution lifecycle

The CPU campaign controller reserves attempt IDs and deterministic seeds in a
locked, atomically replaced ledger. It stages one child per attempt before making
that child scheduler-visible. For peers, the existing authenticated submission
contract carries exact runtime/reference requirements and SHA-addressed inputs;
results pass through durable artifact-transfer receipts. No new transport exists.

A child obtains a physical GPU allocation, invokes Docker synchronously and waits
for exit. Commands use `docker run --rm`. RFD4 uses the shared `run_docker_job()`;
BindCraft2 and RFD3 have dedicated subprocess calls; RFdiffusion uses
`design._run_shell_steps()`. Output normalization follows native exit. Parent
reconciliation verifies artifacts, provenance and candidate identity before recording
attempt receipts. Downstream validation creates independent scheduler jobs.

Each native invocation loses its GPU-resident models at process/container exit.
Mounted checkpoints and disk caches can remain warm. A warm filesystem cache is
not a persistent GPU model. One workflow can invoke several Docker commands;
one container can contain several Python processes; neither is a scheduler job count.

## 4. Engine-specific native behavior

All statements in this section are **verified in inspected source**, unless marked
otherwise. Multi-attempt performance and restart behavior have not been measured.

| Engine | Native count controls and loop | Initialization scope | Currently verified campaign job size |
| --- | --- | --- | --- |
| BindCraft2 | `max_trajectories` bounds claimed trajectories; `number_of_final_designs` is an accepted-design stopping target, not a trajectory count. `run_campaign_arm()` loops sequentially. | `run_campaign()` constructs design/MPNN/validation models before the loop. Native worker dispatch can create multiple Python processes inside one container. Each worker has its own models. | 1 trajectory, 1 final-design target, 1 worker/GPU |
| RFD4-Proteina | `num_replicates` repeats **each input**; `batch_size` controls dataloader/GPU batching. The wrapper sends at most 100 requested replicates per native command, chunking larger requests into more Python subprocesses within the outer container. | `run_generation()` calls `load_ckpt_n_configure_inference()` once before dataset/dataloader construction and `trainer.validate(model, ...)`. One generator model is reused across loader batches. | 1 design request, native dataloader batch 1 |
| RFdiffusion | `inference.num_designs`; the Python entrypoint loops over design indices. This is sequential sampling, not a vectorized GPU batch. | Sampler selection occurs before the design loop. Sampler initialization loads checkpoints/model; `sample_init()` prepares a design without constructing a new sampler each iteration. | No campaign adapter enabled |
| RFdiffusion3 / Foundry | `diffusion_batch_size` controls vectorized samples; `n_batches` expands each specification into sequential batch examples. Potential native output count is specifications × batches × diffusion samples. | `RFD3InferenceEngine.run()` calls cached `BaseInferenceEngine.initialize()` before `_run_multi()`; one trainer/model is used across the loader. | No campaign adapter enabled |

### BindCraft2 details

`campaign.py:run_campaign_arm()` claims a trajectory through `CampaignProgress`,
folds the campaign JAX key by trajectory number and runs design with the existing
model objects. Model parameters and JIT caches can be reused within a worker.
Different lengths/settings may compile new graphs. Desperation settings may rebuild
the validation pool, so “all models load exactly once for every configuration” is
not established. Native workers are separate Python processes, not vectorized
trajectories sharing one GPU model allocation.

Trajectory names, drawn trajectory seed, metrics tables, per-trajectory structures,
filter decisions and `.campaign_state.json` provide incremental evidence. Native
state increments `trajectories` when a trajectory is **claimed**, before it finishes.
After interruption that count alone does not prove completed trajectories. The
current one-attempt adapter requires successful native-state/output reconciliation.
Native filter acceptance is authoritative under that reconciliation; process exit
alone is not scientific acceptance. One trajectory can undergo native MPNN redesign
and native validation, potentially producing several sequence-level rows in ordinary
workflows. Trajectory-only mode requires an explicit downstream sequence handoff.
No blanket external ProteinMPNN step should be inserted into ordinary full BindCraft2.

### RFD4 details

The Studio runtime's `range(0, count, 100)` creates distinct chunk output roots and
uses `(seed + chunk_start) % 2**32` for chunk seeds. The 100 boundary is a wrapper
chunk size, not a verified safe GPU batch or campaign batch limit. Multiple prepared
inputs multiply the native replicate count. Campaign support is limited to its
existing single selected target/guided binder path.

Dataset example IDs distinguish replicates (`_rep…`); generated samples may also
receive sample suffixes. The native seed is assigned to the invocation/job split,
and the same `inference_seed` is written into output metrics. Multiple replicates
therefore do **not** establish one independently persisted campaign seed per
attempt. Repacking a job, changing vectorized batch size or retrying a subset can
change the RNG stream; equivalence must be demonstrated before enabling batching.

Native structure/metric-pool work is scheduled as batches complete; final metrics
CSV aggregation occurs at epoch end. The Studio wrapper creates its portable
result index only after successful native subprocess completion and metrics-table
checks. An interrupted container can leave useful CIFs without a complete metrics
CSV/index. Those files are forensic partial results, not current verified campaign
receipts. Native task filters/evaluations are not an authoritative campaign
acceptance flag. RFD4 co-design provides a sequence and complex: the verified
adapter needs same-sequence downstream validation, not automatic ProteinMPNN redesign.

Ordinary RFD4 options can also start a separate folding-oracle container/process.
Its lifecycle differs from generator checkpoint reuse. P5 disabled refolding; the
one-generator-container observation does not generalize to oracle-enabled tasks.

### RFdiffusion details

The inspected script seeds by design index when deterministic mode is enabled;
otherwise designs consume a process RNG stream. A top-level Studio seed must not
be advertised as independently applied to every trajectory without verifying its
mapping to native indices. Output-prefix/index PDB and TRB files are written per
design; optional trajectory PDBs follow. The TRB's `time` starts inside the design
loop, after sampler initialization, and includes sampling plus output work. It is
not isolated GPU inference time. A killed design can leave an incomplete PDB/TRB
pair; previous complete pairs can be inspected but campaign recovery is unverified.
Backbones need a real ProteinMPNN sequence lineage step before binder validation.
There is no authoritative binder acceptance flag from backbone generation.
Custom scaffold preparation can add a separate container before inference.

### RFdiffusion3 details and count hazard

`BaseInferenceEngine.initialize()` loads the checkpoint and constructs the trainer
once while its `initialized_` flag is false. `_run_multi()` performs a model forward
per specification/batch and dumps each sample's CIF/JSON output. Metadata retains
the invocation seed and example/sample IDs; it does not supply distinct independent
campaign seeds for every vectorized sample. Partial batch files are not atomic
campaign receipts. `_model_forward()`'s native elapsed interval includes data movement
and trainer work, lacks an explicit CUDA-synchronized kernel boundary and must not
be relabeled pure inference time.

The current tile computes `batch_size=min(diffusion_batch_size,num_designs)` and
`n_batches=ceil(num_designs/batch_size)`. For example, 5 requested designs with a
native batch of 2 can produce 6 outputs per specification. `_normalize()` selects
only the requested number. **This existing ordinary-workflow behavior was not
changed in Phase 6C.** A campaign adapter must retain/account for all six native
attempts and outputs, or schedule an exact remainder invocation; it must not reuse
that truncation behavior as campaign accounting. Sequence-head/optional Foundry-MPNN
outputs still need verified backbone-to-sequence lineage. No campaign acceptance
can be inferred from generation or optional sequence generation alone.

### Cancellation and failures across engines

The existing scheduler can stop a run and its mounted Docker containers; remote
cancellation uses the existing authenticated peer mechanism. Operational cancellation
is available. None of these inspections verifies a graceful scientific checkpoint
at an arbitrary cancellation instant. Successful process exit, completed native
attempt, usable candidate, native rejection, and downstream acceptance remain
separate facts. Preserve logs/files on failure; incomplete output is not a negative
scientific result and not an automatically recoverable successful attempt.

## 5. Other installed generators and capability restrictions

The existing campaign registry also describes BoltzGen, PXDesign and Genie3, all
with campaign support deferred. Their new native-execution fields default to unknown,
not unsupported-by-the-engine, and cannot enable campaign dispatch. App-level
inspection shows BoltzGen has `num_designs` and `diffusion_batch_size`, PXDesign
uses `N_sample`, and Genie3 sets dataset `n_sample`; none proves attempt-level seeds,
acceptance, complete output cardinality, partial recovery or model reuse.

Protpardelle's ordinary adapter passes `--num-samples`, `--batch-size` and disables
external MPNN with `--num-mpnn-seqs 0`. Proteina-Complexa has an ordinary generator
manifest. These are not verified distributed campaign adapters/registry entries.
Their native lifecycle remains unverified for this contract. Image installation or
a count option is insufficient to expose native batching. Other specialized design
tools require the same engine-specific audit before registry inclusion.

`NativeExecutionCapabilities` now describes multi-attempt invocation support,
semantics, maximum verified campaign job size, maximum measured native invocation
size where available, multiple candidate records, seed provenance, partial recovery,
source model reuse, runtime-measured reuse and evidence. `None` means unknown;
0 verified campaign attempts means no enabled adapter, not a native engine limit.
BindCraft2 and RFD4 retain campaign maximum 1; RFdiffusion/RFD3 retain 0. Source
model reuse is true for the four inspected engines; runtime-measured reuse is false.
A future planner must also require explicit supported choices and a complete
attempt-level adapter; it must never derive authorization from source reuse alone.

## 6. Measured timing evidence and unknown phases

P4/P5 evidence is reused from the original reports and retained metadata:

| Observation | Interval | Meaning / limit |
| --- | --- | --- |
| P4 RFD4 child `20261028-101059-dc2accd7` | Worker 09:11:01–09:12:58 UTC, 117 s | Whole worker, including setup/inference/output; no model-phase split |
| P4 child creation to worker start | 09:10:59–09:11:01 UTC, about 2 s | Local scheduling/start interval, second-resolution timestamps |
| P5 RFD4 on `pc004032` | 11:45:03–11:45:42 UTC, 39 s | Whole worker on RTX 5090 |
| P5 RFD4 on `pc004033` | 11:45:08–11:45:59 UTC, 51 s | Whole worker on RTX 5090 |
| P5 parent | 11:44:48–11:48:53 UTC, 245 s | Includes distributed generation, transfer, serial validation and reconciliation |
| P5 Boltz validation | 66 s and 68 s | Separate complete validation workers |
| P5 follow-up | Both allocations overlap 12:27:38–12:29:01 UTC, about 83 s | Two independent Boltz jobs on two local 4090s; no native generator batching |

These are retained event timestamps, not new inference. P4/P5 allocation GPU-seconds
are elapsed leased seconds × leased GPUs; they are not utilization or kernel timings.
Cross-host clocks/polling can also affect end-to-end interval interpretation.

New **CPU-only measured** Docker client intervals, three samples per cell:

| Installed image | `/bin/true` median (range), s | `python3 -c pass` median (range), s |
| --- | --- | --- |
| BindCraft2 | 0.215 (0.199–0.216) | 0.265 (0.249–0.313) |
| RFD4 | 0.211 (0.201–0.228) | 0.264 (0.245–0.332) |
| RFdiffusion | 0.230 (0.174–0.286) | 0.207 (0.186–0.227) |
| Foundry | 0.209 (0.207–0.211) | 0.260 (0.227–0.260) |

Commands used `--pull=never --network=none --runtime=runc`, an overridden entrypoint,
no GPU request, no inputs/weights mounts and no engine imports. Images were already
local; daemon and filesystem caches were warm/unspecified. CPU test activity overlapped.
Three samples are descriptive, not a statistical performance study. The variation
also explains why subtracting these independent measurements cannot estimate Python
startup, model loading or any phase in P4/P5. No ratio-based speedup is claimed.

| Lifecycle component | Available evidence after this change |
| --- | --- |
| Scheduler queue delay | Existing submitted/allocated/worker timestamps; compute only same-clock, semantically matched intervals. Creation-to-remote-start also includes probes, staging and transfer, so is not pure queue delay. |
| Docker command launch through return | New paired `docker_command` host spans, exit codes or exception types |
| Native process initialization / imports | Unknown; no-op Python is only a CPU baseline |
| Checkpoint loading / GPU model initialization | Repeated scopes verified in source; durations unknown |
| Compilation | Shape-dependent source paths verified; duration unknown |
| Native inference | Unknown as an isolated synchronized interval |
| Native serialization | Source writes located; duration unknown |
| Candidate normalization | New paired host spans; no retrospective P4/P5 measurements |
| Artifact transfer | Existing durable requests, bytes and checksums; isolated transfer latency unknown |
| Container termination / removal | Included in synchronous Docker command interval; isolated stop/removal time unknown |

## 7. Implemented instrumentation and compatibility corrections

`execution_lifecycle.jsonl` is a best-effort append-only host diagnostic. A span has
schema version, unique span ID, phase, UTC timestamp, PID, start/finish events and
monotonic elapsed seconds. Docker spans include return code or exception type.
Retries and multi-command workflows append distinct span IDs. An unmatched start
means interrupted/incomplete telemetry, not a completed measurement. No sequences,
commands or exception text are copied into timing events. File-write errors do not
change scientific execution or mask an original exception. Telemetry is not a ledger
or acceptance authority, is not fsync-durable, and may have a truncated tail after a
hard kill. Consumers must ignore malformed/unpaired records.

Instrumentation covers the shared Docker runner, shared design command steps,
BindCraft2's dedicated Docker call, and RFD3's dedicated Docker call. Normalization
spans cover BindCraft2, RFD4, RFD3 and generic generator outputs. The BindCraft2
normalizer's span covers that normalizer function; RFD4/RFD3 also include their
candidate-index writes. These are host phase scopes, not identical native operations.
Generic normalization is timed only when execution succeeded. Existing commands,
mounts, GPU selection, seeds, output selection and acceptance rules are unchanged.
Shared command steps use `docker_command` for direct Docker calls and
`execution_command` for other commands, including shell wrappers.

Timing files currently remain in the worker run directory. They are not added to
peer scientific receipt schemas or automatically fetched by the parent. The origin
and peers must receive a separately reviewed deployment before peer jobs emit the
new spans. Original runs are never backfilled with invented timings.

The campaign config now accepts `dispatch_batch_size` as an input alias for legacy
`batch_size`. Conflicting values fail validation; successful normalization persists
the same old `batch_size` config. A supplied campaign `native_attempts_per_job` is
accepted only as integer 1, then removed from the normalized legacy config. Larger
values, floats, booleans and `auto` fail explicitly. This is a guard against silently
ignoring a requested behavior, not a larger-batch implementation. Existing absent
fields produce byte-compatible config values and permit original ledger recovery.

The UI calls the budget “Global maximum native attempts” and explains the effective
one-attempt/job execution, generation capacity and independent validation concurrency
before submission. Existing widget keys, workflow/API fields and defaults remain.
No universal tile batch option or incompatible schema migration was introduced.

## 8. Recommended configuration hierarchy — proposed

```text
Workflow
  → Generator tile / campaign configuration
  → Native execution jobs (scheduler resource leases)
  → Native attempts / trajectories (global authorization)
  → Candidate normalization (zero, one or several records per attempt)
  → Downstream validation jobs
  → Global campaign ledger and acceptance/diversity selection
```

Four independent concepts:

| Scope | Setting | Unit / authority |
| --- | --- | --- |
| Campaign | `accepted_target`, `max_attempts`, acceptance/filter/diversity policies | Global candidates and authorized native attempts across all workers |
| Generator tile | `native_attempts_per_job` | Requested attempts within one native invocation; engine-specific capability gate |
| Execution | `max_concurrent_generation_jobs` | Independent scheduler jobs; separate from attempts in a job |
| Execution | `max_concurrent_validations` | Independent downstream validation jobs; separate from generation |

The tile must display native semantics and the resolved invocation count, input
multiplicity, seeds, vectorized microbatch, expected maximum native outputs, resource
requirements and any remainder job. Native GPU microbatching remains a distinct
engine option (for example RFD4 dataloader batch or RFD3 diffusion samples).
A tile may create several jobs, and one job may require one GPU for many sequential
attempts. Tile count is not a job count and attempt count is not a GPU count.

Proposed selections: `1, 5, 10, 25, 50, 100`, with **Auto** as a later optional policy.
These are UI candidates, not supported values today. Enabled campaign choices remain
`[1]` for BindCraft2/RFD4; other engines expose no campaign choices. Disabled choices
must state the missing evidence. Never silently clamp 25 to 1 or transform it into
25 independent jobs while calling that native batching.

Auto should initially resolve to 1 and later choose only an already verified size,
persisted before authorization. It should consider remaining global budget, target
gap, unresolved validations, eligible hardware/memory, measured initialization cost
and stopping/cancellation granularity. Auto must not learn from unbounded exploratory
generation or assume linear speedups. Job membership must not change after launch.

## 9. Existing ledger limitations and proposed job/attempt contract

Today `batch_id` is dispatch grouping. There is no first-class native multi-attempt
execution job in the ledger: one attempt points to one child; retries retain attempt
ID/seed and record retry history. **Do not reinterpret historical `batch_id` as a
container identity.** Authorization can reserve several attempts atomically, and
candidate receipts are attempt-based, but this alone is not a batched adapter.

A separately versioned ledger extension should add:

| Record | Required content |
| --- | --- |
| Execution job | Stable job ID; immutable ordered attempt manifest and digest; invocation/config version; scheduler/cluster IDs; resource lease; worker/GPU identity; start/end; retry executions; output receipt frontier |
| Attempt | Global ordinal/ID; authorized seed or explicit native seed/key mapping; execution-job membership; native input/replicate/trajectory index; actual-start evidence; completion/failure/not-started/unknown state; output candidate IDs |
| Candidate | Stable attempt-scoped ID; native output ID/hash; sequence/structure ancestry; native acceptance or unavailable; validation stages; diversity outcome; final selection disposition |
| Receipt | Campaign/job/attempt/execution identity; manifest digest; output checksums; stable idempotency key; terminal/partial marker and completeness evidence |

Reserve every job's attempt manifest and budget atomically **before** submission.
A request for 5 attempts with only 3 authorized slots remaining may resolve to a
3-attempt job only if the engine supports that exact size and the UI/manifest records
it; otherwise create verified smaller jobs. Never authorize 5 against a budget of 3.
A crash between reservation and submission must be recoverable by manifest identity.

Counters must separately expose authorized attempts; submitted/launched execution
jobs and execution retries; actual started, completed, failed, interrupted, cancelled,
not-started and unresolved native attempts; native accepted/unavailable candidates;
downstream validated candidates; diversity-qualified candidates; pending validations;
selected candidates and returned-but-unselected candidates. Count a job's GPU lease
once, not once per contained attempt. Report failed execution events separately from
unique failed attempts: the current `failed_interrupted_attempts` includes retry
history and is not a unique-attempt counter. Current `attempts_started` is based on
scheduler child start, not a native per-trajectory start marker. New detailed counters
must not falsely relabel these historical observations.

Every produced output must have a recorded disposition. Native rejection with zero
candidate rows is still a completed native attempt. A candidate absent because a
process died is not a scientific rejection. Multiple candidates from one attempt
must not be counted as multiple attempts, and candidates must never silently move
between attempts because their native IDs collide.

## 10. Global stopping, partial failure, retry and cancellation — proposed

The existing controller already subtracts active attempts and pending validation
candidates from dispatch headroom, reconciles newly queued children before budget
finalization, and drains late validation before target finalization. **Preserve these
fixes.** It can retain surplus valid candidates. A full future job's authorized
attempt count must replace the old assumption that each active child represents one
attempt. Candidate multiplicity and pending validation load must also be explicit.

Near the target, reduce the next **unlaunched** job to an engine-verified exact size,
prefer one-attempt jobs, and pause generation while sufficient candidates await
validation. For already-running jobs, request “stop before next native attempt”
only where an engine boundary is verified; otherwise let the bounded job drain.
Native accepted-design stopping cannot substitute for the global validated/diversity
target. If native rejection or diversity reconciliation reduces the selected count,
resume only while the globally authorized budget permits it.

Potential overshoot is bounded in attempts by the remaining authorized in-flight
manifest, not by concurrency alone. Candidate overshoot additionally depends on a
verified per-attempt candidate-output bound; without that bound do not promise a
numerical candidate overshoot limit. Retain surplus outputs, acceptance evidence and
returned-but-unselected reasons. Do not discard them to force counters to the target.
The accepted target can be reached before all validation completes; finalization
must still reconcile relevant in-flight attempts, receipts and validation results.

| Event | Required semantics |
| --- | --- |
| Partial native success then process failure | Verify and receipt complete individual attempts; retain incomplete outputs with unresolved status; never turn the entire job into a fresh attempt |
| Retry before native launch | Reuse the same reservation, ID, seed and manifest; a new execution ID records operational retry |
| Retry after launch | Resume only with trustworthy attempt completion/RNG state. RFD4 currently forbids automatic retry after a Docker command is recorded; preserve that guard. A base seed alone is insufficient for subset replay. |
| Duplicate receipt | Same attempt/content is idempotent; conflicting content or ownership is an explicit reconciliation error, not a second acceptance or silent overwrite |
| Controller restart | Reconcile durable job submission, remote identity, manifests, artifacts and receipt frontier before reserving more work. Origin-controller recovery remains supported; cross-origin takeover is separate work. |
| Queued cancellation | Persist not-started/cancelled members and cancel scheduler work; do not count them as completed or silently recycle authorization |
| Active cancellation | Request native boundary stop if supported; otherwise scheduler Docker stop, mark unknown/incomplete members, retain files and drain verified receipts |
| Worker/peer disconnect | Keep unresolved reservations; reconnect and query the same authenticated job. Do not resubmit merely because polling timed out. |

Strict physical-attempt budgets need an explicit policy for post-launch retries:
re-executing an unfinished stochastic trajectory can consume another physical attempt
even when its logical seed is unchanged. Separate logical reservation and execution
cost; do not claim a hard physical-attempt bound from logical-ID reuse. Keep
post-launch multi-attempt retries disabled until resume proof and accounting exist.

The deterministic simulated-group tests use existing per-attempt ledger primitives
and diagnostic `execution_job_id` metadata. They exercise reservation limits,
separate attempt/candidate cardinality, surplus retention, partial failure, stable
seeds, restart, duplicate receipts and cancelled/not-started members. They do not
implement scheduler submission, native progress, partial output verification or
safe batched retries. A schema-v2 implementation must add transaction/crash tests
for those boundaries, malformed manifests and conflicting receipts.

## 11. Four-GPU scheduler implications

The verified pool remains `pc002975` with two RTX 4090s, and `pc004032`/`pc004033`
with one RTX 5090 each. P5 exercised remote RFD4 on both 5090s; its follow-up exercised
two concurrent Boltz validations on the local 4090s. It did not establish simultaneous
four-device utilization or multi-GPU native inference.

A sequential native job can reserve one GPU while processing several attempts.
Vectorized microbatching can increase per-device VRAM requirements. Neither implies
that the scheduler should allocate several GPUs. BindCraft2's internal worker fanout
must remain constrained to the scheduler lease and the verified one-worker setting
unless a separate resource audit permits more workers.

Keep generation and validation concurrency independent, but use the same physical
lease inventory so they cannot allocate the same exclusive GPU twice. When users
request four generation jobs plus two validations, explain that the four-device pool
may queue some work. Retain workstation ID, GPU UUID/device mapping, allocated CPU,
VRAM requirements, exact image/reference compatibility, native attempt membership,
remote scheduler IDs and artifact hashes for each job. Heterogeneous hardware may
need different verified sizes; there is no universal safe VRAM estimate from count.
No scheduler compatibility bypass or artifact-transfer contract change is needed here.

## 12. Backward-compatible API and P6 interface plan

Implemented now: input-only `dispatch_batch_size` alias and explicit one-attempt
campaign guard; additive capability metadata; descriptive UI captions. Normalized
old config remains unchanged, including the product
`batch_size * max_concurrent_batches` for legacy generation dispatch capacity.
The alias is scoped to the campaign object; it does not rename RFD4/Protpardelle
native `batch_size` elsewhere. Ordinary saved workflows keep their native options.

Proposed next API version: optional tile-level `native_execution` settings and an
execution-level explicit `max_concurrent_generation_jobs`. Preserve requests with
neither new field exactly. When new concurrency is supplied, define conflict rules
against legacy product fields and show the effective value; do not multiply the new
job count by the old dispatch size. Keep legacy serialized names readable and expose
a read-only effective execution plan before submission. Additive capability responses
must distinguish ordinary native support, verified campaign choices and measured
hardware limits. Unsupported requests fail preflight with engine-specific reasons.

The interface should group global target/budget/policy under Campaign, native
attempt count and engine semantics on Generator, and job/validation concurrency plus
eligible workers/GPUs under Execution. Display requested versus resolved job sizes,
expected number of jobs, possible native output multiplicity, resource leases,
partial recovery/cancellation limits and budget/target headroom. Do not put one
universal batch-size slider on every workflow tile. Auto stays unavailable until
its resolution policy and bounded behavior are reviewable.

## 13. Estimated benefits and unverified assumptions

A verified sequential engine could avoid repeated process/model construction by
packing attempts in one invocation. If an unchanged per-invocation cost S were
measured, N attempts in jobs of K could avoid approximately
`(N - ceil(N/K)) * S` of that cost before accounting for parallel scheduling,
shape recompilation, output work, early stopping and validation backpressure.
This is a conditional cost model, not a measured benefit. S is currently unknown.
The measured no-op Docker intervals cannot stand in for S.

Potential benefits include fewer Docker launches, fewer imports/checkpoint reads,
retained per-process models/JIT caches and reduced job/transfer overhead. Potential
costs include longer resource leases, less flexible heterogeneous placement,
stranded attempts near the target, validation queues, larger failure domains,
VRAM growth with vectorized batches and changed deterministic RNG behavior.
Native model-object reuse does not prove full checkpoint/compilation amortization,
scientific equivalence, better candidate yield or faster campaign completion.

## 14. Concrete staged follow-up plan

1. **Measure native phases without changing science.** Add opt-in native markers
   around process entry/import completion, model/checkpoint construction, GPU setup,
   compilation, attempt start/end and serialization in one verified engine. Separate
   CPU wall time from synchronized GPU timing. Review the instrumented entrypoint,
   then use one bounded PD-L1 attempt with the retained settings and scheduler lease;
   preserve native command, target/checkpoint/image hashes and logs. Only if needed,
   authorize a small explicitly bounded multi-attempt comparison to test reuse.
2. **Add a read-only execution planner.** Resolve tile settings, exact attempt counts,
   input multiplicity, hardware requirements and legacy concurrency. Reject unverified
   values. Add deterministic plan tests for budget remainder, multiple inputs and
   RFD3's rounded native outputs. Do not activate native batching yet.
3. **Version the job/attempt ledger.** Add immutable manifests, idempotent submission,
   actual native-start evidence, attempt receipts, partial/unknown states and unique
   counters. Provide explicit schema-v1 read compatibility; migrate only with a saved
   copy and no reinterpretation of historical batch IDs. Test crashes at every durable
   boundary, duplicate/conflicting receipts, restart and validation drain.
4. **Enable one engine-specific sequential path behind a capability gate.** BindCraft2
   needs globally mapped trajectory keys, one worker and incremental-state reconciliation;
   RFD4 needs verified replicate-to-attempt seed/provenance and partial index recovery.
   Select the first engine based on those proofs, not assumed speed. Start at the
   smallest separately verified size; compare full outputs and accounting to single
   jobs without changing acceptance criteria. RFD4 stays dataloader batch 1 initially.
5. **Introduce conservative stopping/backpressure.** Shrink unlaunched jobs near budget
   or target, account for pending validation, retain all surplus outputs and require
   verified boundary cancellation. Validate resource leases and same-origin restart.
6. **Extend only after evidence.** RFdiffusion needs ProteinMPNN lineage; RFD3 additionally
   needs exact native cardinality and sample provenance. Verify vectorized microbatching
   separately from sequential attempts. Promote proposed sizes 5/10/25/50/100 and Auto
   only individually with engine/hardware evidence. Persistent GPU workers are a
   separate architecture and remain outside this plan's initial implementation.

## 15. Focused verification

The project Conda interpreter was used:

```text
/home/user/mambaforge/envs/mn-molecular-design-studio/bin/python -m pytest \
  tests/test_native_execution_contract.py tests/test_execution_timing.py \
  tests/test_generator_capabilities.py tests/test_distributed_campaign.py \
  tests/test_rfd4_proteina.py tests/test_rfdiffusion3.py \
  tests/test_native_generators.py tests/test_bindcraft_boltzgen_contracts.py \
  tests/test_distributed_campaign_ui.py tests/test_candidate_peer_transfer.py \
  tests/test_workflow_queue_lifecycle.py -q
281 passed in 20.52s
```

An earlier run had 194 passes and one old-label assertion failure; the UI
assertion was updated to the clearer budget label. Final targeted checks passed:
51 native-generator/UI/timing tests in 14.37 s, then all four timing tests in 0.16 s
with explicit coverage of early failure and shell-wrapper phase labels.
`compileall` for changed modules and `git diff --check` passed.
No broad suite, GPU pytest suite, native multi-attempt performance test or deployment
acceptance campaign was rerun. New tests cover timing failures/nonzero exits,
telemetry write failure, alias conflict/type validation, rejection of unverified
choices, old-config recovery equivalence and simulated multi-attempt ledger grouping.
Existing controller tests retain late-validation finalization, peer transfer,
resource selection and RFD4 post-launch retry protection.

Retained P4/P5 checks read the original ledgers and normalized indexes; all three
candidate sequences match their validation evidence and retain
`identity_check.sequence_chain_identity=pass`. Config normalization equals saved config.
P5 ledger SHA-256 is still
`205efddd4c18cea8f8ca9d9ee22a9b5ad193964f5e059854631bf2465e891b39`;
P4 is `0f575ca3e563181423524207c34ebc5fbcf1c6c9b6ff8c5bf71409ae10bc1460`.
Original ledgers, predictions and structures were not edited.

## 16. Files and commits

Implementation files:

- `mn_molecular_design_studio/engine/workflow_engine/core/execution_timing.py`
- `mn_molecular_design_studio/engine/workflow_engine/core/docker_runner.py`
- `mn_molecular_design_studio/engine/workflow_engine/workflows/{bindcraft2,design,generators,rfd4_proteina,rfdiffusion3}.py`
- `mn_molecular_design_studio/engine/workflow_engine/workflows/distributed_campaign.py`
- `mn_molecular_design_studio/engine/workflow_engine/workflows/generator_capabilities.py`
- `mn_molecular_design_studio/distributed_campaign_config_ui.py`
- `tests/test_execution_timing.py`, `tests/test_native_execution_contract.py`, `tests/test_distributed_campaign_ui.py`

Documentation/evidence files:

- `docs/testing/native-execution-batching-architecture.md`
- `docs/testing/evidence/phase6c-container-startup.json`
- `docs/testing/evidence/phase6c-retained-campaign-checks.json`

`d88e8c1` — `feat: describe native execution contracts and record lifecycle timings`
contains the verified implementation and tests.

`d8a1b69` — `docs: record native batching architecture and bounded timing evidence`
contains the original Phase 6C report and its two evidence JSON files.
The unrelated `run_app.py` and `tool_references.py` are excluded from all commits.

## 17. Phase 6C.1 — Opt-in native instrumentation

The installed `rfd4-generate` console entrypoint is
`script_utils.entrypoints:generate`. It captures import warnings, imports
`rfproteina.generate.main_hydra`, queues those warnings and calls the Hydra entrypoint.
The isolated `rfd4_native_diagnostics.py` launcher follows that same warning policy
and entrypoint. It temporarily wraps known functions, forwarding the original
arguments and return objects, and restores the hooks on return or exception.
It does not edit the installed image, checkpoint, sampling code or filtering rules.

Enable explicitly with `options.native_phase_timing=true` in an **ordinary local**
RFD4 generation request. The existing workflow signature and default requests are
unchanged. The launcher is staged only for opt-in runs. A container-local timing
path activates it, preserving native Hydra argv; both the original and diagnostic
commands are retained in `native_command.json`. Campaign/scoring requests with this
option are rejected before job creation. The runtime also refuses to redirect an
`rfd4-evaluate` command through a generator timing launcher. No peer source was
updated, and diagnostics are not an enabled distributed campaign feature.

The launcher verifies the installed `generate.py` and `proteina.py` SHA-256 values
against the reviewed Phase 6C sources before invoking generation. An unreviewed
engine revision fails diagnostics explicitly rather than silently applying stale
hooks. A CPU-only `--help` invocation verified imports, source hashes and hook
installation against the actual image without loading a model or generating a design.
Optional CUDA-library warnings in that no-GPU check are retained in its logs; it
exited successfully. Those CPU import timings are not used as GPU timing evidence.

Native events are written to `artifacts/rfd4_proteina/native_lifecycle.jsonl`, with
unique span IDs, UTC timestamps, PIDs, monotonic elapsed wall intervals and exception
types. Writes are best-effort; a missing/truncated or unpaired record is incomplete
telemetry, never a successful native attempt receipt. Native exceptions and return
objects are preserved. No synchronization, Torch/NumPy RNG operations, worker-count
changes or scientific parameter overrides are added. Timing inevitably adds small
host overhead; this was not quantified by a paired inference comparison.

| Hook / boundary | Interpretation |
| --- | --- |
| `native_entrypoint_imports` | Native package/import work after Python enters the diagnostic script; interpreter startup is outside this interval |
| `load_ckpt_n_configure_inference` | Combined model/checkpoint loading and native inference configuration, including any work called by that function; not isolated disk-read time |
| `get_dataset_from_cfg` | Native dataset/input preparation |
| `Trainer.__init__` | Trainer object construction; not a complete GPU setup measurement |
| `torch.cuda._lazy_init` | Existing PyTorch lazy-initialization calls, observed only while `is_initialized()` is false; includes recursive calls and does not represent every possible CUDA/library initialization path |
| `setup_run_compile`, `reserve_shape_cache` | Compile configuration/cache setup if called; actual deferred compilation can remain inside prediction |
| `Trainer.validate` | Inclusive framework setup, inference, metrics, output and drain work |
| `Proteina.predict_step` | Native prediction plus strategy/decoding/atom-array work; not synchronized GPU-kernel time |
| `MetricsPool.drain_generated_structure` | Waiting for already submitted metrics/structure work, potentially in other processes; not total worker CPU time or isolated CIF serialization |
| `Proteina._save_and_log_csv` | Final metrics CSV serialization/logging call |
| `native_entrypoint_main` | Script imports and native main until return; interpreter atexit/teardown is outside this interval |
| `native_subprocess_launch_through_exit` | Adapter-side process launch/wait through actual native subprocess completion |

Spans are inclusive and may nest or overlap. In particular, recursive CUDA-init
spans must not be summed, and prediction/drain/serialization are already inside
`Trainer.validate`. Process completion is established by the parent wait/exit status,
not merely an entrypoint-return event. No timing event changes ledger accounting.

### One bounded PD-L1 measurement

Run `20261028-152136-81e2ed87` completed on `pc002975` GPU 0, RTX 4090 UUID
`GPU-195a1966-1148-340e-5fe2-bac2d7a9821a`. The existing shared scheduler leased
one GPU and eight CPU cores with a 16 GiB VRAM request; scheduler job
`a01f96704fd8422fb740e86d855af437` completed and released its reservation.
No worker service was restarted and no runtime source was deployed to peers.

The run reused P5 attempt 1's retained PD-L1 chain-A target and configuration:
one native replicate, dataloader batch 1, 60 steps, seed 1185903269, length range
60–120, unchanged hotspots, deterministic mode and native condition filtering,
refolding disabled, reward search off and compilation disabled. The complete
resolved native YAML was equal to the retained reference YAML. The diagnostic
option is not translated into a scientific Hydra setting. Target SHA-256 remained
`7da1c888011cdbd0937e3ff6b73429aa069e4c496cfa0cb8d8b53cc06ca73b57`.

The image remained `sha256:c34efb5afca78c0341cf996c9dadc8102a7fb3760e1792697a8779803e0e44f5`.
HF checkpoint revision was `a06e2509c7033b010fc43614e9fd840303588502`;
`model.safetensors` and `config.json` hashes matched P5's retained requirements:
`5b549dd61adb2f6c6b5cb95bda13f40f5e31556b4e6c755edc0ca6f011c68fb0` and
`154ea927ab8005b7eac2d047a792147fb4606ab4878e496aa8414b6b40dad368`.

| Measured host wall scope | Seconds |
| --- | ---: |
| Whole worker, second-resolution metadata | about 120 |
| Docker command through return/removal | 118.586 |
| Native subprocess launch through exit | 95.432 |
| Native entrypoint imports | 11.641 |
| Model/checkpoint loading and configuration | 5.610 |
| Dataset preparation | 16.460 |
| Trainer construction | 0.014 |
| Outermost observed CUDA lazy-init call | 0.018 |
| Inclusive native validation loop | 56.592 |
| Native prediction plus decoding | 17.250 |
| Metrics/structure future drain | 30.679 |
| Final metrics CSV serialization | 0.017 |
| Native entrypoint through return | 90.625 |
| Host candidate normalization/index writing | 1.701 |

These values are **not additive**. Three CUDA-init call pairs were observed, with
small recursive spans inside the outer call. `compile_lengths=null` disabled the
requested Torch compile path; no compile-setup spans were expected or observed.
Other lazy library/kernel work is not isolated. Detailed checkpoint disk I/O,
model GPU transfer, CIF serialization in metrics workers, interpreter startup,
atexit teardown, and pure GPU inference/compilation remain unmeasured separately.
The gap between whole-container and native-process intervals also includes outer
adapter imports, preparation and indexing; it is not isolated Docker startup.

This run demonstrates nontrivial per-invocation import/model work distinct from
prediction and whole-worker time. It does not establish the savings achievable by
packing attempts. Dataset and metrics costs have different reuse/parallelism scopes;
they must not be called model-loading time or assumed to amortize. Host CPU activity and filesystem-cache warmth were not controlled, and focused CPU
test work ran during the measurement session. The single local
4090 sample cannot be compared directly with P5's 5090 workers or used to predict
multi-attempt throughput. No uninstrumented paired GPU repeat or batch-size benchmark
was run. Runtime multi-attempt model reuse therefore remains unmeasured in capabilities.

One 116-residue co-designed candidate was normalized. Its target sequence matches
the retained 115-residue chain A; its designed chain B matches the normalized binder
sequence. The new binder sequence differs from the original P5 worker's sequence;
this was a separate execution on a different GPU. Current local and peer image IDs
also differ, and no controlled bitwise-equivalence comparison was performed. Transparent hook tests verify forwarding of arguments,
result identity and exceptions; this run verifies configuration and provenance, not
numerical equivalence across platforms. The fresh candidate is a timing artifact,
not a replacement P5 result or a new validated campaign acceptance. No Boltz run
or downstream acceptance was added.

An initial harness submission, `20261028-152050-33730633`, used the default run root
outside configured Results. The scheduler rejected it before Docker/native launch.
It is retained as an excluded setup failure with zero native attempts, and its local
metadata was marked failed. The corrected submission above is the only new native
generation attempt. Neither submission altered original campaign ledgers.

Raw paired timings, allocation, config comparison, command records, checkpoint/output
hashes and exclusions are in [Phase 6C.1 timing evidence](evidence/phase6c1-native-timing.json).
Full logs, the staged diagnostic script, scientific outputs and CPU import check are
retained under the run directory and
`/mnt/data/RESULTS/mn-molecular-design-studio/verification/phase6c1-native-timing/`.
The timing evidence identifies the exact staged launcher hash. A final defensive
command-builder guard only adds rejection of evaluation commands; its generation
argv equivalence is unit-tested without another inference run.

## 18. Phase 6C.1 — Read-only execution planner

`workflows.execution_planning.plan_campaign_execution()` is a deterministic data
projection. It does not probe GPUs/peers, reserve attempts, create jobs/ledgers,
transfer files, modify config, or dispatch work. Existing campaign validation is
reused for acceptance/resource compatibility; count arithmetic does not bypass it.
The result is never a submission request or authorization manifest.

The planner is available in the generator campaign form's **Read-only native
execution plan** preview and as a module CLI. Preview controls are separate from
saved tile/API configuration. Selecting 5/10/25/50/100 only changes the hypothetical
projection; it cannot change saved `distributed_campaign`, execution readiness or
actual native job size. The actual campaign settings remain displayed above it.
The existing APIs and workflow signatures require no schema migration or deployment.

The plan reports:

- Engine capabilities, provenance/recovery limitations and verified campaign choices.
- Global accepted target and maximum authorized native-attempt budget.
- Requested count, descriptive choices, supported choices and effective count.
- Compact job-size histograms and job-count upper bounds if the entire budget is
  used; actual acceptance-driven stopping and runtime/compute limits can reduce work.
- Native counts, Python invocation counts, nominal output multiplicity, rounding
  and proposed exact-remainder arithmetic.
- Legacy dispatch-product generation concurrency and independent validation concurrency.
- Resource defaults matching the existing adapters: BindCraft2 four CPU cores and
  RFD4 eight unless the parent tile/campaign overrides them; RFD4 campaign minimum
  VRAM 16 GiB. Validation defaults to four CPU cores unless its stage overrides them.
- Workstation/GPU selectors, caller-supplied inventory eligibility explanations,
  explicit exclusions, and unresolved freshness/identity/compatibility/resource data.
- Unsupported configuration explanations with `supported_execution_plan=null` and
  no effective native count. It never silently replaces an unsupported request by 1.

`submission_authorized` is always false. A supported configuration is not proof of
current resource availability or a guarantee of acceptance. Descriptive native
choices remain `[1,5,10,25,50,100]`; only `[1]` is supported for current BindCraft2/RFD4
campaigns. RFdiffusion/RFD3 and other deferred engines have no campaign choices.
Auto remains unavailable. Proposed `max_concurrent_generation_jobs` is diagnosed as
unsupported and cannot override the existing legacy dispatch product.

### Global budget and exact remainders

For a global budget of 12 and hypothetical size 5, the plan shows two jobs of 5
and one remainder job of 2, totaling 12 attempts. This arithmetic does not enable
those jobs. Unsupported requests have no runnable plan.

For several inputs, the hypothetical planner distributes the **global** budget
across input ordinals, with at most one extra attempt for earliest inputs, and
uses one input per hypothetical job. For example, budget 7 across three inputs
becomes input budgets 3/2/2, not seven replicates for each input. This convention is
explicitly hypothetical: current campaign adapters reject multiple-input planning
as an execution setting. Inputs with no budget receive no jobs. Job histograms are
compact even for large budgets; no attempt records are materialized.

The independent `describe_native_invocation()` helper describes ordinary native
multiplicity without implying campaign support:

- RFD4 with five replicates and three inputs requests 15 attempts. Dataloader batch
  2 produces eight loader batches with an exact final partial batch, not 16 attempts.
  Wrapper chunking above 100 replicates is counted separately from GPU batching.
- RFdiffusion accepts one input per native invocation; multiple inputs require
  separate invocations. Its design loop is sequential.
- BindCraft2 trajectory count is an upper bound; early native target stopping and
  sequence-level output multiplicity remain explicit, with unknown candidate bound.
- RFD3 with five requested outputs, three inputs and diffusion batch 2 requests 18
  native samples through the current rounded count. The proposed exact alternative
  is two full batches of 2 plus a remainder of 1 per input, 15 samples total, requiring
  a separate remainder invocation with the existing native interface. That alternative
  is labeled `implemented=false`. If rounded hypothetical physical attempts exceed
  the global budget, the planner reports the overrun and produces no supported plan.

Nominal outputs are not guaranteed usable structures, normalized candidates,
scientifically accepted candidates or validated candidates. Search/oracle count
expansion is outside the simple no-search ordinary-count helper.

### GPU inventory contract

Inventory can be supplied through `eligible_gpus` records with workstation `node_id`,
`gpu_id`, optional GPU UUID/model, local identity, online/compatibility status,
`observed_at`, `snapshot_current`, capacity and free/allocated fields. The pure planner
never invents missing evidence or derives exact peer compatibility from a GPU model.
Freshness, exact image/checkpoint/runtime compatibility and required resources must
be explicit for snapshot eligibility. A GPU can be eligible while allocated and
therefore require queuing; immediate-start availability requires positive free-device
and free-VRAM evidence. Missing inventory is `unknown_not_probed`, not an empty or
verified available pool. The current tile preview has no supplied inventory and
states that limitation; the existing authenticated placement/scheduler checks retain
authority over real submissions and physical leases.

A standalone JSON request contains `campaign` plus optional planner arguments such
as `engine_id`, `native_attempts_per_job`, `input_count`, `native_batch_size`,
`designs_per_input`, `generation_cpu_cores` and `eligible_gpus`. Run it with:

```text
python -m mn_molecular_design_studio.engine.workflow_engine.workflows.execution_planning request.json
```

The CLI reads the request and writes JSON to stdout. It performs no scientific or
scheduler operation. Deterministic requests and resulting examples are committed in
[planning evidence](evidence/phase6c1-planning-examples.json). The retained
`planner_request.json` in the verification directory was exercised through the CLI.

## 19. Phase 6C.1 — Compatibility and tests

Campaign execution and the Phase 6B ledger are unchanged. Multi-attempt execution,
new adapter support, persistent GPU workers, post-launch batched retry and new
acceptance thresholds remain disabled. Native diagnostics are opt-in and isolated;
ordinary default argv and sampling configuration retain their existing behavior.
The scheduler lease, late-validation finalization, artifact transport and exact peer
compatibility contracts are unchanged. No source was pushed or deployed and no
additional Codex agents were used.

The original P4/P5 ledger hashes and candidate evidence are checked again after
measurement/planning. No original scientific file was edited, promoted or replaced.
The unrelated modified `run_app.py` and untracked `tool_references.py` are excluded.

Focused tests cover global budgets and exact remainder arithmetic across multiple
inputs; native count multiplication versus GPU batching; RFD4 100-replicate wrapper
chunks; RFD3 rounded counts and proposed exact alternatives; unsupported/unknown
engines and choices; resource defaults and inventory uncertainty; legacy aliases and
config immutability; explicit proof that planning cannot call create/enqueue/dispatch;
preview controls leaving saved execution unchanged; transparent diagnostic return/
argument/exception forwarding; source-drift refusal; write failures; opt-in type checks;
and campaign/scoring diagnostics rejection before job creation. Existing controller,
peer transfer, lineage, queue and late-validation regressions are included.

The focused regression sweep passed **320 tests in 21.68 s**, covering the Phase 6C
regression set in section 15 plus `tests/test_execution_planning.py` and
`tests/test_rfd4_native_diagnostics.py`. Final planner/diagnostic checks, including an
additional invalid-null-resource case, passed **40 tests in 1.00 s**. These sets overlap;
they are not 360 distinct tests. Changed Python files passed syntax parsing,
`git diff --check` passed, both evidence JSON files parsed, the saved supported plan
matched the current pure planner exactly, and the CLI exercised that same request.
No full broad suite or additional native/validation inference was run.

The installed-image CPU import/source/hook check exited 0. The one GPU measurement
returned exit code 0, one candidate, matching sequence/structure roles and unchanged
resolved native config/checkpoint requirements. Both original P4/P5 ledger hashes,
legacy config normalization and retained same-sequence validation identity checks
remained unchanged after the measurement. The diagnostic scheduler allocation was
empty after completion, independently confirming release of GPU 0.

## 20. Phase 6C.1 — Files and commits

Implementation changes:

- New `workflows/execution_planning.py` and `workflows/rfd4_native_diagnostics.py`.
- Scoped opt-in staging/launching in `workflows/rfd4_proteina.py` and `workflows/rfd4_runtime.py`.
- Boolean diagnostic-option validation in `rfd4_proteina_input.py`.
- Read-only preview in `distributed_campaign_config_ui.py`, with parent CPU context
  from `pipeline_node_ui.py`.
- New `tests/test_execution_planning.py` and `tests/test_rfd4_native_diagnostics.py`.
- Updated architecture report and new timing/planning evidence JSON files.

Paths under `workflows/` above are relative to
`mn_molecular_design_studio/engine/workflow_engine/`. The standalone staged native
launcher and measurement outputs remain recoverable in the diagnostic run.
`f4586b1` — `feat: add native phase timing and read-only campaign planning`
contains the reviewed implementation and tests. The following documentation commit,
`docs: record native phase measurements and execution planning`, contains this updated
report and the timing/planning evidence. Its actual hash is reported in the completion
message. Only these scoped source, test, report and evidence files are committed.
