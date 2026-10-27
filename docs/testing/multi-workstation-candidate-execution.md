# Phase 6: candidate execution and scheduling — 27 October 2026

Implementation is on `fix/scientific-pipeline-validation`, starting from clean
revision `65904db`. Nothing was pushed. No peer installation, scheduler service,
account configuration, scientific threshold or native docking controller was changed.

Local execution is verified with real scientific outputs on two GPUs. Peer candidate
submission/transfer/recovery is implemented and regression tested, but actual
cross-workstation inference is **not verified**. Online peers have a different code
revision. Their eligibility checks were retained; no incompatible job was forced
onto them.

## Infrastructure reused

The graph controller still submits ordinary allowlisted workflow jobs. It uses the
existing run catalog, worker processes, central scheduler bridge or local
`WorkerScheduler`, shared CPU/GPU/RAM leases, cancellation, normalized candidate
collections and content-based graph cache. Peer work uses the existing authenticated
Store submission/action queues, command leases, execution receipts, agent, input blobs,
artifact commands and workstation heartbeats. There is no additional execution
service, scheduling database or distributed-computing dependency.

The initial checkout passed **1,552 CPU tests, 30 skipped**. The final broad sweep
passed **1,575 tests, 33 skipped, 17 warnings, 272.40 s**. An additional pair of
explicit partial/default fan-in regressions passed after collection of that sweep.
The final provenance/graph/cache subset passed **159 tests**, and the scheduler
retry fence passed **16 focused tests** (including its new regression). Previous
scientific contracts, task catalogue, two-stage MPNN/cache behavior, peer APIs and
active-learning docking remain covered by that suite.

## Implemented behavior

- Native collection batching remains the default for saved and new graphs.
  Explicit `candidate_execution.mode=sharded` freezes ordered input units and
  effective native requests in job-specific SHA-256 manifests before submission.
- Configurable batch size, in-flight bound, retries, CPU/RAM/VRAM requirements,
  local/automatic/explicit enrolled workstation destination, and GPU selection.
  The existing worker-count setting continues to bound workstation processes.
- Independent generation requests preserve constraints and record derived,
  non-repeating batch-start seeds; a diffusion trajectory stays on one worker.
- Candidate batches preserve sequences, chain roles, input IDs and parent links.
  New result IDs use a stable node/batch/ordinal namespace. Native IDs remain in
  provenance. Duplicate/ambiguous identities and missing input lineage are errors.
- Aggregation follows input batch order, verifies actual artifact contents and
  reports missing candidates. Progress counts represented completed inputs, rather
  than enqueued batches. Partial output never enters the completed-pool cache.
- Each result retains the input-manifest hash, effective request, worker/GPU,
  software/native checkpoint provenance, hardware inventory, output checksums and
  retry history. Completed equivalent results remain attributable across workers;
  no bitwise reproducibility across devices is promised.
- `auto` now remains an automatic **single-GPU** request through resource parsing.
  The previous conversion to `all` reserved both GPUs for a shard. Local workers
  now receive concrete allocated selectors before starting their native workflow.
  Admission limits are recorded before dispatch, rather than after a possible claim.
- Device reports include GPU name, UUID, free/total VRAM, utilization, compute
  capability, driver and driver CUDA API version. Installed PyTorch images can be
  inspected without GPU allocation/model loading to obtain compiled CUDA targets.
  Matching follows [NVIDIA's CUDA binary/PTX compatibility rules](https://docs.nvidia.com/cuda/cuda-programming-guide/01-introduction/cuda-platform.html).
- Peer eligibility checks retain exact code/container requirements and now compare
  explicitly required checkpoint file SHA-256 identities. Reference directories
  are never recursively copied. Worker reference probes hash named files.
- Portable input bundles recursively rebase selected candidate JSONL structures
  and file references; they include per-file hashes and retain the existing
  512 MB bound. A source directory represents the selected snapshot, not its
  original run/database. Upstream execution/model records remain opaque, origin-tagged
  provenance; they do not require installing prior-stage models or copying prior
  collections. Transfers do not require identical absolute paths.
- Remote outputs arrive through existing authorized artifact commands. Each blob
  is checksum verified and atomically published locally. Interrupted/expired
  transfers resume without rerunning inference. Foreign non-output paths remain
  machine-qualified provenance rather than being opened locally.
- Failed/lost batches retain attempts, seeds and native run artifacts. Retries
  operate at the configured batch unit; batch size one provides individual retries.
  Durable creation tags recover accepted local submissions after lost replies.
  A controller retry waits while the original scheduler retries/reconciles or its
  process is still alive, preventing competing attempts.
- A disconnected remote attempt retains its identity. After the configured
  recovery window it is quarantined as lost and cancellation is queued. It is
  **not** blindly resampled elsewhere while the original worker might still run.
  Other branches/batches can continue. Queued and running children can be cancelled.
- All-upstream dependency barriers remain the default. Explicit partial policy
  permits successful branches/pools to proceed and carries missing-branch/batch
  information into Results. No scientific selection threshold was changed.
- Optional incremental CPU analysis consumes completed upstream batches while the
  remaining GPU batches run. Only single-parent collection analyses are enabled;
  selection, Results, multi-parent and other complete-pool stages retain barriers.
- The existing pipeline editor/run view exposes settings, actual worker/GPU
  allocations, retries, errors, waiting reasons and completed candidate-stage
  progress. GPU selectors use discovered devices rather than hard-coded 0/1 labels.

## Adapter coverage and boundaries

Eligible collection adapters include ProteinMPNN/SolubleMPNN/shared Foundry MPNN,
AbMPNN, AF2, AF3, RF3, Boltz-2, Protenix and other existing optimized-prediction
profiles, ipSAE, Rosetta metrics and compatible collection QC workflows. The optional
interface-preserving soluble stage remains an ordinary lineage-preserving MPNN node.

Independent generation batches are enabled for RFdiffusion, RFD3, RFD4, Genie3 and
PXDesign. Native complete pipelines and generators outside this allowlist retain
normal execution. ESM-IF/Caliby and other unlisted sequence adapters are not newly
sharded. No model inference is partitioned across computers.

Only **shared ProteinMPNN, Boltz-2 and CPU sequence composition** were freshly
executed for this phase. Other listed adapters have scheduling/contract coverage,
not a new real-GPU validation claim. JAX images and custom CUDA extensions do not
have a universal compiled-target probe; their native runtime checks still apply.
Live streaming CPU/GPU overlap and disconnection/retry failure cases are verified
with deterministic workers; the real CPU check ran after GPU prediction.

## Real execution and benchmarks

Two local NVIDIA GeForce RTX 4090 devices, each 24,564 MiB, driver 595.91.07,
were used with separate real worker processes and shared resource leases. The
existing long-running central service was neither restarted nor reconfigured.
This is **one workstation with two GPUs**, not two workstations. No local RTX
5090 was present. The configured peers were online but code-ineligible; their
architecture/runtime compatibility was not experimentally established.

Installed images include `mn-ligandmpnn:latest`,
`mn-boltz2-opt:2.2.1-cu130-batch2` and `mn-protein-qc:1.0.0`.
The first two expose CUDA 12.8 and 13.0 builds respectively, including `sm_120`
and targets compatible with the observed CC 8.9 devices. Matching an image's
compiled targets does not certify all third-party kernels or every model mode.

Each GPU trial used the same four candidate inputs for its workload, fixed seed
20261027 and equivalent scientific settings. Single-GPU native batching used one
batch of four. Sharding used two batches of two. ProteinMPNN generated one sequence
per backbone at temperature 0.4; Boltz used one sample, one recycle, 20 diffusion
steps, single sequences, no target template and no automatic ipSAE. No checkpoints
or reference databases were downloaded.

| Workload (4 inputs) | One GPU, native batch | One GPU, two shards | Two GPUs, two shards |
|---|---:|---:|---:|
| ProteinMPNN, final-code retest | 23.75 s; 606.4/h | 23.02 s; 625.6/h | 13.34 s; 1,079.3/h |
| Boltz-2, optimized collection | 99.77 s; 144.3/h | 200.82 s; 71.7/h | 114.49 s; 125.8/h |
| Actual cross-workstation execution | Not run | Not run | Not verified |

ProteinMPNN improved **1.78×** versus its native single-GPU batch in the final run.
Boltz's small collection was **14.8% slower** across two GPUs and about twice as
slow when split on one GPU. Its single native batch remains preferable here.
Do not generalize these four-input smoke measurements to larger proteins/cohorts.
Desktop activity and a concurrent CPU test suite affected these non-exclusive
measurements; the raw device telemetry records that background load.

The first complete Boltz comparison took 129.42 / 228.42 / 161.48 seconds.
Inspection found approximately 4.5 GiB of installed checkpoint/CCD files included
in every candidate output checksum traversal. Output hashing now checks scientific
artifacts while retaining model identity separately in native/runtime provenance.
The optimized comparison took 99.77 / 200.82 / 114.49 seconds. Content checks on
scientific inputs/outputs remain ordinary SHA-256 checks; no stat-only digest cache
was retained. Timing improvements include background-load variation and should not
be treated as an isolated causal benchmark of hashing alone.

The final ProteinMPNN controller times were 2.93 / 0.94 / 0.85 seconds; scheduler
loop/admission times were 1.91 / 2.12 / 1.27 seconds. These are recorded elapsed
controller/scheduler work, not native model initialization timers. Local transfer
volume/time and real failure/retry cost were zero. Device utilization samples and
peak **whole-device** VRAM are in the JSON; they include display/other activity and
are not per-process VRAM measurements. Native model initialization and native CPU
preparation are included in wall time but **not separately instrumented**.

A further real CPU sequence-QC collection processed all four retained Boltz results
in **4.95 s**, using two CPU workers with four cores each and **zero GPU devices**.
All selected GPU/CPU jobs have completed native results and candidate collections;
execution evidence is not merely successful scheduler submission.

See [machine-readable benchmark evidence](candidate-execution-benchmarks-2026-10-27.json).
Raw manifests, requests, worker logs, native outputs, structures, metrics, checksums,
telemetry and complete attempt ledgers remain in the artifact paths in that file.
Pytest transcripts are retained under `.tmp/scientific-phase6-2026-10-27/logs/`.

## Rerunning

Use the existing `mn-molecular-design-studio` Conda environment. CPU tests never
start model inference unless GPU tests are explicitly enabled.

```bash
python -m pytest -q
python -m pytest tests/test_candidate_execution.py tests/test_candidate_peer_transfer.py -q
MDS_GPU_TESTS=1 \
MDS_CANDIDATE_BENCHMARK_SOURCE=/absolute/retained/backbone/run \
MDS_CANDIDATE_BENCHMARK_PREDICTION_SOURCE=/absolute/retained/sequence/run \
python -m pytest tests/test_candidate_execution_gpu.py -v -s
python scripts/summarize_candidate_benchmarks.py /absolute/benchmark.json --output summary.json
```

The optional CPU-only execution check accepts
`MDS_CANDIDATE_BENCHMARK_CPU_SOURCE=/absolute/retained/prediction/collection`.
No remote configuration changes or downloads are part of these commands.

A saved graph can opt in without changing native scientific parameters:

```json
{"candidate_execution":{"mode":"sharded","placement":"local","batch_size":4,
 "max_in_flight":2,"retries":1,"failure_policy":"all","min_vram_gb":8,
 "ram_gb":0,"cpu_cores":4,"stall_timeout_seconds":1800,"incremental_cpu":false}}
```

`placement` also accepts `auto` or an enrolled workstation UUID. Per-node
`data.execution` overrides the graph policy. The controller stays local for sharded
UI submissions. For single-candidate retries choose batch size one. Explicit GPU
selectors still apply; engine modes needing multiple GPUs retain their native
multi-device allocation. The 8 GiB default is an admission floor, not a universal
model memory estimate: set measured limits for large inputs. No reliable generic
sequence-length-to-VRAM estimator was invented.

## Changed files, commits and next steps

- Candidate planning, manifests, aggregation, recovery and CPU streaming:
  `candidate_execution.py`, `pipeline_graph.py`, `workflows/graph_pipeline.py`.
- Scheduler admission/receipts: `core/jobs.py`, `core/scheduler.py`,
  `core/workflow_queue.py`, `services/workflow_automation.py`,
  `workflows/optimized_prediction.py`.
- Peer portability/eligibility: `cluster/candidate_batches.py`, `bundles.py`,
  `capabilities.py`, `gpu_runtime.py`, `agent.py`, `store.py`.
- Monitoring: `pipeline_graph_ui.py`.
- Regressions: `test_candidate_execution.py`, `test_candidate_peer_transfer.py`,
  `test_cluster_capabilities.py`, `test_workflow_queue_lifecycle.py`.
- Real checks/benchmark export: `test_candidate_execution_gpu.py`,
  `scripts/summarize_candidate_benchmarks.py`, this guide and its JSON evidence.

Logical commits: `f6f6884` candidate sharding; `c2364be` resource admission;
`7d0e436` peer execution; `3cc396f` recovery/partial/incremental analysis;
`a0f26f9` checksum optimization/CUDA matching; `fb0c7a9` monitoring;
`22b5854` scheduler retry fencing; `02c9fe2` partial-pool propagation and explicit
fan-in regressions. `94e536d` contains real execution tests, benchmark artifacts
and this report. The subsequent portability fix bounds upstream provenance staging.

Before claiming multi-workstation readiness, deploy the same reviewed revision to
trusted peers through the normal authorized release process, refresh capability
reports, and execute a two-candidate cross-machine test including an interrupted
artifact transfer. That deployment was not authorized here. Further priorities are
real heterogeneous-GPU/JAX/custom-kernel validation, longer cohorts with repeated
uncontended measurements, and separate native initialization/preparation timers.
Keep native batching for small Boltz collections; use measured sharding for larger
independent cohorts and the verified ProteinMPNN case.
