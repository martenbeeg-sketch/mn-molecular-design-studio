# Distributed Acceptance Campaigns — Phase 6B P0–P3

P4 multi-engine adapter, lifecycle and validation results: [Phase 6B P4 report](distributed-campaign-multi-engine.md).

## Delivered commits

- P0 ledger: `a6b10c3` — `feat: add durable distributed campaign ledger`
- P1 BindCraft2 attempt adapter: `0e2f95c` — `feat: add bounded BindCraft2 attempt adapter`
- P2/P3 integrated controller, validation and diversity: `b614a46` — `feat: execute distributed acceptance campaigns`
- GPU-discovered P2 staging and stopping fix: `a4d94b9` — `fix: reconcile queued campaign children before stopping`
- Native rejected-trajectory counter: `cbcb538` — `feat: report native rejected trajectories`

P2 dispatch and P3 acceptance share one controller and persisted candidate state, so those implementation changes landed together. P0 and P1 remain separate milestones.

## Controller and BindCraft2 behavior

An ordinary BindCraft2 design tile can enable `distributed_campaign` with a global `accepted_target`, `acceptance_level`, `max_attempts`, `batch_size`, `max_concurrent_batches`, GPU selector, resource requirements, retry limit and optional runtime/compute budgets. Validation stages and diversity rules are optional. Distributed mode is disabled by default, preserving ordinary BindCraft2 execution and saved graphs.

The parent is a CPU scheduler job. It creates isolated one-trajectory BindCraft2 child jobs, stages a separate target/settings/output directory for each, assigns a reproducible seed and explicit GPU, and leaves resource leasing to the existing scheduler. Children use `number_of_final_designs=1`, `max_trajectories=1`, and one worker per GPU. The parent validates native `.campaign_state.json` against normalized candidate rows and candidate structures before recording a receipt. A successful worker process and a scientifically accepted design remain separate states.

The durable ledger is stored at `artifacts/distributed_campaign/ledger.json` beside the parent run. A file lock protects attempt authorization and receipt updates; writes use atomic replacement and fsync. Attempt IDs and seeds survive retries, duplicate receipts are idempotent, and a restarted origin controller reconciles child jobs and results before dispatching more work. Native and validation retries use the configured bound. Cancellation stops queued work and requests stops for active child jobs.

## Acceptance and stopping semantics

- **Native:** count verified candidates whose normalized `native_pass_filters` flag is true, after reconciling that count with the native campaign state. `native_rejected_trajectories` separately counts completed trajectories minus native accepted designs. `rejected_candidates` counts rejected normalized candidate rows; a screen rejection with no normalized candidate row appears in the trajectory counter and retained native diagnostics.
- **Validated:** count a native-accepted candidate only after every configured validation stage passes its explicit required measurements and rules. Measurements must be finite and valid, and evidence must retain the same binder sequence, source ancestry, target/binder chains and sequence identity check. Selectors are stage-, tool- and output-candidate-qualified. An `analysis.ipsae` stage must follow an optimized complex prediction.
- **Diversity-qualified:** count native-accepted candidates after validation stages, if configured, and global diversity reconciliation. Supported filters include exact-sequence deduplication, sequence identity cutoff, target-aligned binder Cα RMSD cutoff, stage-qualified interface-cluster measurement and per-cluster cap. Candidates are processed in stable attempt order. Rejected rows and reasons remain in the parent collection when a normalized row exists.

The selected level controls global target completion. The controller stops authorizing attempts at the target or configured budget, waits for active atomic jobs, reconciles their outputs and retains any surplus valid candidates. If diversity rejection reduces the final count, generation continues while attempt budget remains. Budget, runtime, compute, cancellation and failure have explicit termination states. `max_compute_seconds` is cumulative GPU allocation lease time (elapsed leased worker seconds multiplied by leased GPU count); it is not measured GPU utilization, and atomic active jobs can create bounded overshoot.

The supported downstream workflow adapters are `refolding.optimized_prediction` and `analysis.ipsae`. Their GPU and CPU resource requests use the regular scheduler. Other validation workflows and arbitrary stage graphs are outside this P0–P3 adapter.

## Test results

The one full CPU sweep ran in the project Conda environment with:

```text
conda run -n mn-molecular-design-studio python -m pytest -q
1602 passed, 33 skipped, 17 warnings in 262.34s
```

Skips were opt-in GPU/container checks and checks requiring `mdtraj`. The full sweep preceded the GPU-discovered staging/stopping follow-up. After those follow-ups and the rejection-counter change, the focused regression set passed **73 tests**:

```text
python -m pytest tests/test_distributed_campaign.py tests/test_candidate_evidence.py \
  tests/test_candidate_selection.py tests/test_bindcraft_boltzgen_contracts.py \
  tests/test_cluster_progress.py -q
73 passed in 3.24s
```

Deterministic tests cover atomic attempt authorization, unique seeds, receipt idempotency, retry identity, asynchronous scheduler reconciliation, native versus validated counts, stage-qualified evidence, validation retry, target stopping, diversity rejection and generation continuation. The campaign controller’s validated-stop test uses a controlled prediction worker with normalized same-sequence evidence; it does not claim a new GPU prediction.

## Bounded local two-GPU campaign

One fresh native campaign ran on the local RTX 4090 pair, using the retained readiness PD-L1 target (`target.pdb`, chain A), the installed `mn-bindcraft2:v1.0.4-cu13-pae-v1` image and the short readiness overrides. The parent run was `20261028-085821-36356a14`; its children were `20261028-085821-84acdebe` and `20261028-085821-c5f60f7f`.

| Attempt | GPU | Seed | Native trajectories | Native accepted | Native rejection |
|---|---:|---:|---:|---:|---|
| `a00000001` | 0 | 1854686050 | 1 | 0 | screen: pLDDT 0.32; iPTM 0.14 |
| `a00000002` | 1 | 880427935 | 1 | 0 | screen: pLDDT 0.30; iPTM 0.11 |

Both scheduler child jobs completed with exit code 0, four CPU cores and their requested GPU IDs. Each authoritative native state recorded one trajectory and one `terminated.screen` rejection. The parent reconciled **2 authorized, 2 started, 2 completed; 0 accepted; 2 native-rejected trajectories; 0 failed/interrupted attempts**. No normalized candidates were emitted, so the consolidated candidate collection is empty and has no duplicate IDs. The campaign ended `budget_exhausted` at 0/2; the parent workflow completed computationally while `campaign_success` was false. GPU allocation accounting recorded 288 GPU-seconds. Docker containers exited and scheduler jobs completed, releasing their leases.

The first post-code-fix GPU launch exposed a child staging bug: a path inside BindCraft2’s working directory was interpreted relative to the parent run. Both children failed before native trajectories started. The controller could also finalize immediately after dispatching a still-queued child. Both cases are covered by the follow-up regression tests and fixed in `a4d94b9`. The successful bounded campaign above ran after those fixes. Two earlier scheduler-root preflight launches were cancelled while still queued; neither ran inference.

The retained readiness inference and this campaign are computational checks. They do not establish binder yield or predictive validity. The unchanged native filters were used throughout.

## Retained validation artifacts

The Results root contained real retained prediction/analysis outputs, including an AF2 initial-guess/refolding/Rosetta chain from 2026-10-26. A read-only scan of 333 candidate JSONL files found no retained complex-refolding candidate that met this campaign’s required same-sequence lineage and `identity_check.sequence_chain_identity=pass` contract. Those legacy measurements were not promoted into validated campaign evidence. No Boltz-2/ipSAE GPU validation run was launched: the bounded native campaign yielded no candidate to validate, and no compatible retained artifact was available. Validated stopping and diversity stopping were exercised by deterministic controlled-worker tests.

## Boundaries and next interfaces

- Execution in this implementation is on the local scheduler. `workstation_selection` is persisted and validated, but cross-workstation dispatch, remote compatibility verification, deployment and controller takeover are not implemented or tested. No remote machine was changed.
- GPU compute budgets use allocation leases, not device utilization measurements.
- P4 adds an RFD4-Proteina adapter with isolated one-request attempts, explicit non-authoritative native acceptance, verified sequence/structure rows and same-sequence downstream validation; see the report above.
- P5 can connect eligible peers to the existing scheduler/receipt transport after source, capability and artifact checks; same-origin controller recovery remains the supported recovery model.
- P6 can expand monitoring around the existing campaign progress and ledger counters. This work adds no separate dashboard.

No repository push, remote deployment or broad benchmark was performed.
