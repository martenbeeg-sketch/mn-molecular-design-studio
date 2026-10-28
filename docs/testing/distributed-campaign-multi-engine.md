# Phase 6B P4 — Multi-Engine Campaigns and Execution Lifecycle

## Result

P4 adds a bounded RFD4-Proteina adapter to the existing P0–P3 campaign ledger and completes one real candidate handoff through Boltz-2 and campaign selection. BindCraft2 campaign behavior remains supported. No persistent inference workers, larger native batches, throughput increase, P5 peer support or P6 dashboard work were added.

The real run used the packaged PD-L1 target as a harmless benchmark target, source chain A, RFD4-Proteina on GPU 0, and one authorized native design request. It produced one 75-residue co-designed sequence/complex. RFD4 does not provide an authoritative campaign acceptance flag, so the ledger recorded native acceptance as **unavailable**, rather than inferring acceptance from a successful process or output file.

Boltz-2 then predicted the retained sequence in the target complex. Candidate identity and binder sequence matched the RFD4 source; the validation evidence records `identity_check.sequence_chain_identity=pass`. The measured outputs were iPTM 0.42641, pTM 0.72618 and mean pLDDT 86.701. The configured policy required those three measurements to be present and valid; it set no score thresholds. The candidate therefore passed the configured **measurement-presence** policy. This is a campaign selection result, not evidence of experimental binding or a claim that a particular interface score is strong.

| Campaign outcome | Count |
|---|---:|
| Authorized / started / completed native attempts | 1 / 1 / 1 |
| Native accepted | 0 |
| Native acceptance unavailable | 1 |
| Validated accepted candidates | 1 |
| Rejected candidates | 0 |
| Failed or interrupted attempts | 0 |
| Cumulative GPU allocation time | 185 seconds |

The validation worker used GPU 0, Boltz-2 2.2.1, one sample, one recycle, 20 sampling steps and no MSA. The campaign's configured attempt budget remained 1 throughout validation; rerunning the controller after recovery did not authorize another generator attempt. The final campaign result reports `target_met=true` and `campaign_success=true`; termination remains `budget_exhausted` because the single authorized attempt was consumed.

The normalized campaign row and same-sequence validation evidence are in the parent run's `artifacts/normalized_candidates/candidates.jsonl` and `campaign_result.json`. The durable attempt, validation and recovery events are in `artifacts/distributed_campaign/ledger.json`. The native RFD4 complex is retained in the child run's `artifacts/rfd4_proteina/output/inference_samples/protein_binder/protein_binder_generated.cif`; the Boltz-2 complex is retained in the validation run. The exact local parent run ID was `20261028-101057-a07ef0a4`, the RFD4 child was `20261028-101059-dc2accd7`, and the validation run was `20261028-101807-b37312d9`.

## Campaign adapters and scientific boundaries

| Engine | P4 campaign support | Attempt and acceptance semantics |
|---|---|---|
| BindCraft2 | Supported; existing P0–P3 adapter preserved | One independently seeded native trajectory per child. Native acceptance comes from authoritative `.campaign_state.json` reconciled against normalized candidates and native filters. Native, validated and diversity-qualified stopping remain available. |
| RFD4-Proteina | Supported for guided protein-binder co-design from one selected target structure | One native design request per child, which can yield zero or one normalized candidate. Seed, chain roles, target mapping, sequence, structure and source provenance are verified. Native campaign acceptance is unavailable; validated or diversity-qualified stopping requires explicit downstream criteria. No automatic ProteinMPNN redesign is inserted. |
| RFdiffusion | Inspected; deferred | Emits backbones. Campaign accounting must include a real ProteinMPNN sequence handoff and verify trajectory-to-sequence ancestry before a binder can be counted. |
| RFdiffusion3 / Foundry | Inspected; deferred | Exposes `diffusion_batch_size` and `n_batches` inside one native invocation, but campaign attempt-to-output accounting and ProteinMPNN/Foundry sequence lineage are not verified. |

BoltzGen, PXDesign and Genie3 remain capability-discovered but campaign-deferred because their native attempt, seed, output or acceptance semantics are not verified. No universal score threshold was added. The supported validation workflows remain optimized complex prediction and ipSAE, with stage-qualified measurements and same-sequence evidence.

RFD4 candidate handoff checks explicit target and binder chain roles, verifies binder sequence against its complex structure, verifies target chains against the staged target through source-to-label mappings and exact target sequence crops, and records seed, attempt, worker, GPU, native candidate ID and component mappings. Downstream prediction takes the existing RFD4 binder sequence and structure ancestry; it does not redesign the binder.

## Container and inference-process lifecycle

The shared `run_docker_job()` writes a command record and uses a synchronous `subprocess.run()` for a Docker command. The container is launched with `--rm`; model memory cannot persist between separate workflow processes. Model images, checkpoint directories and caches are mounted from local storage, so filesystem cache warmth may persist, but no app-owned GPU-resident model process or persistent inference worker exists.

| Engine | Current app lifecycle | What remains unknown |
|---|---|---|
| BindCraft2 | Each ordinary design job launches one Docker command. A distributed attempt is explicitly bounded to one trajectory and one final design, so every campaign attempt is a separate child/container/process. | Native multi-trajectory execution is available in ordinary runs, but model reuse across its internal trajectories was not established from app-side code or phase timing. |
| RFD4-Proteina | Each workflow job launches one Docker container and staged Python adapter. `rfd4-generate` receives `num_replicates`; up to 100 requests are sent in one native command. The adapter chunks requests above 100 into separate `subprocess.run()` calls. P4 campaign children request one output each, so the campaign starts a new container/process for each attempt. | Whether the native RFD4 loop retains an initialized model across multiple replicates in one command was not verified. The campaign adapter uses one replicate and does not test or alter that behavior. |
| RFdiffusion | The app builds a Docker command that runs the RFdiffusion Python entrypoint once with `inference.num_designs`; custom scaffold preparation can add a separate Docker command. | Reuse of the loaded model over multiple designs in one native process was not instrumented. |
| RFdiffusion3 | The app starts one Docker process with a `design` command and passes `diffusion_batch_size` and `n_batches` for that run. Checkpoints are mounted from `/weights`. | Internal model reuse across Foundry batches was not instrumented. No P4 campaign or GPU run was made. |

The campaign `batch_size` setting continues to bound controller dispatch concurrency; it is not a native per-container trajectory batch size. This work does **not** implement a default batch of 50, adaptive native batches, additional in-flight generation, or persistent workers. The existing RFD4 and RFdiffusion3 native batch options remain part of their ordinary workflow calls, but their multi-output model retention is unverified. No new GPU benchmark was run solely to compare batch sizes.

## Measured and inferred overhead

| Evidence | Measurement | What it establishes |
|---|---:|---|
| P4 RFD4 campaign child | 117 seconds from scheduler allocation to completion; 115.5 GPU allocation seconds | Total child allocation/process duration for one request. It does not separate Docker startup, model initialization, checkpoint reads, sampling or writeout. |
| P4 Boltz-2 validation | 69 seconds from worker creation to completion; about 69.5 GPU allocation seconds | Total worker wall/lease duration for one sample. It does not separate container startup, model initialization, checkpoint reads, prediction or normalization. |
| P4 combined campaign | 185 cumulative GPU allocation seconds | Sum of recorded generator and validation allocation time, not GPU utilization time. |
| Prior P0–P3 BindCraft2 campaign | 288 GPU allocation seconds for two one-trajectory attempts | Existing allocation evidence in [the P0–P3 report](distributed-acceptance-campaigns.md); no startup/model-load breakdown. Its referenced raw run directories are not currently retained. |
| Prior RFD4 + Boltz-2 integration run | 112 seconds for RFD4 and 77 seconds for Boltz-2, end-to-end job wall times | Existing evidence recorded in `scientific-pipeline-validation-2026-10-27.json`; it also has no initialization-phase split. |
| Candidate execution benchmark | `model_initialization_seconds: null` | Existing benchmark explicitly did not measure initialization, and is for sequence/docking tools rather than these binder generators. |

Repeated process and container launches are directly visible in the app code. Therefore, **the loss of in-memory model state between separate attempts is inferred from process boundaries**; the native packages are responsible for loading the models required by each new process. Docker startup duration, checkpoint-load duration, GPU model-initialization duration and compilation duration were not measured independently. The RFD4 test used `compile_model=false` (native `compile_lengths=null`), so it provides no JAX-compilation timing. We did not confirm whether JAX compilation occurs in the tested RFD4 inference path. Boltz-2 timing likewise contains no phase-level model-load or compile measurement.

The measured totals show the cost of complete jobs, not the savings a larger native batch would achieve. Establishing that savings requires phase instrumentation or comparable existing multi-output execution evidence; this P4 work deliberately did not start a large GPU batch experiment.

## Restrictions and compatibility

No target-specific check, checksum gate or special PD-L1 production restriction remains. PD-L1 was only the selected benchmark target. Ordinary RFD4 jobs and saved pipelines retain their default behavior because the optional `distributed_campaign` parameter defaults to absent. BindCraft2 campaigns retain their prior adapter and stopping semantics. The RFD4 workflow entrypoint accepts this new optional request field so it can invoke the existing controller; no existing field became required and no saved pipeline migration was made. The internal `enqueue_workflow_call()` helper also gained a default-on `dispatch` argument so a controller can attach durable parent/attempt metadata before the child becomes scheduler-visible. Existing queue callers keep dispatching immediately.

The new RFD4 campaign adapter is currently limited to guided protein-binder co-design against one selected target structure and chain set. It rejects campaign requests with secondary specificity structures, separate protein-sequence inputs or other RFD4 task modes because this adapter has not verified their attempt and candidate mappings. This limits only new RFD4 distributed-campaign requests; the ordinary RFD4 tasks remain available. RFdiffusion and RFdiffusion3 are not newly enabled for campaigns.

RFD4 retries are permitted only when the child run exists and no native Docker command has been recorded. A failed child after native launch, or a missing child whose launch state cannot be recovered, is not retried automatically because a stochastic output may exist and the engine does not expose an authoritative resume receipt. This prevents a retry from silently creating a second trajectory under the same attempt and seed. The effect is that these interrupted distributed RFD4 attempts are marked failed and consume their authorized attempt; ordinary RFD4 workflows and BindCraft2 retry behavior are unaffected.

The scheduler now treats a CPU-only coordinator allocation as `cpu` when rewriting nested GPU selectors. That preserves a child validation stage's own GPU request instead of replacing it with the parent's `none` allocation. The public top-level request is unchanged; regression coverage exercises the selector preservation.

## Verification

Focused regression set:

```text
conda run -n mn-molecular-design-studio python -m pytest \
  tests/test_distributed_campaign.py tests/test_generator_capabilities.py \
  tests/test_rfd4_proteina.py tests/test_workflow_queue_lifecycle.py -q
143 passed in 6.31s
```

Full CPU suite:

```text
conda run --no-capture-output -n mn-molecular-design-studio python -m pytest -q
1617 passed, 33 skipped, 17 warnings in 258.69s
```

Skips were the opt-in GPU/container checks and tests requiring `mdtraj`. Focused coverage includes capability discovery, attempt and candidate accounting, RFD4 result verification and retry safety, same-sequence validation requests, resource-selector preservation and workflow queue lifecycle.

No broad throughput benchmark, multi-GPU scale-up, remote deployment or P5/P6 work was performed.
