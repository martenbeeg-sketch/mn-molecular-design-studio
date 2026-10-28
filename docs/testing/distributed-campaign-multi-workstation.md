# Phase 6B P5 — Multi-workstation deployment and execution

**Date:** 2026-10-28
**Branch:** `fix/scientific-pipeline-validation`
**Purpose:** reproduce the existing PD-L1 RFD4-Proteina → Boltz-2 benchmark across enrolled workstations. This is computational workflow verification; it does not establish experimental binding, affinity, specificity, or activity.

## 1. Initial and final workstation revisions

| Workstation | Pre-deployment revision | Deployed runtime source revision |
|---|---|---|
| `pc002975` | `fix/scientific-pipeline-validation` at `d7317f0` | `fix/scientific-pipeline-validation` at `4f9ba772fcb48c1961d74c99f2f5704f34c06197` |
| `pc004032` | Existing `main` at `d56667628cc708b960c13824d8406f8f6a8bcd35` | `fix/scientific-pipeline-validation` at `4f9ba772fcb48c1961d74c99f2f5704f34c06197` |
| `pc004033` | Existing `main` at `d56667628cc708b960c13824d8406f8f6a8bcd35` | `fix/scientific-pipeline-validation` at `4f9ba772fcb48c1961d74c99f2f5704f34c06197` |

The two peer `main` branches were preserved. The P5 source fingerprint reported by the live agents on all three machines was `9ae11faca1ac4b2ada7e68f2c0c1ddc637dbd4d3a04d66bba0f63c53890f7cda`. The report-only commit described in section 15 was made after execution on the local checkout; it changes no Python or native source, and does not change that live source fingerprint.

## 2. Deployment actions

The reviewed source was transferred to the two peers as complete Git bundles and fast-forwarded onto `fix/scientific-pipeline-validation`. Bundle contents and commit ancestry were verified before deployment. The peer `main` refs, workstation identities and credentials, Results directories, databases, caches, reference files, and model installations were preserved. The local modified `run_app.py` and untracked `tool_references.py` were preserved and not included in the P5 commits. No GitHub push was made.

## 3. Services restarted

The peer `studio-api.service` and `studio-agent.service` services were restarted after their source updates, while no scientific inference was active. Both returned online and published fresh heartbeats. The local API and agent were active for submission and coordination. No service was restarted during the benchmark.

## 4. Live capabilities and compatibility

All three agents reported online at the end of the run, with the same source fingerprint. The peers reported one RTX 5090 each, compute capability 12.0, driver 595.99.02 and CUDA driver API 13020. `pc002975` reported two RTX 4090 GPUs, compute capability 8.9, driver 595.91.07 and CUDA driver API 13020.

Both peers reported the required RFD4 image content fingerprint `sha256:a0b8ed6e9c418b6ee7435ff633d77f1248abdab5b8d7b249b18b3a2dc5d3f82a`. Its CUDA 13.0 runtime profile lists `sm_120` among its compiled targets. Both peers published the requested RFD4 reference files with matching SHA-256 values before dispatch: config `154ea927ab8005b7eac2d047a792147fb4606ab4878e496aa8414b6b40dad368`, model `5b549dd61adb2f6c6b5cb95bda13f40f5e31556b4e6c755edc0ca6f011c68fb0`. Both immutable submissions retained those exact image and checkpoint requirements; neither used a compatibility bypass.

The peer snapshots also listed Boltz-2 images, but Boltz validation in this run was scheduled locally. Its `execution_inputs.json` recorded `mn-boltz2-opt:2.2.1-cu130-batch2` and image ID `sha256:3054688e91908156071238adee167ad4963cfbca05c8a4cf0e55a634ba5a2a55`, which matched local `docker image inspect`.

## 5. Reused PD-L1 benchmark configuration

The run reused the packaged PD-L1 structure and target chain A, 115 residues. The transferred target CIF was 67,676 bytes with SHA-256 `7da1c888011cdbd0937e3ff6b73429aa069e4c496cfa0cb8d8b53cc06ca73b57`; its chain-A sequence hash was `c434254957dd7a7d80ce4ca859f9a6d101a1205f4ff822177dd0b3c31f02e968`.

P4 reference run `20261028-101057-a07ef0a4` used one RFD4-Proteina design and then Boltz-2 on that same sequence. Its RFD4 child was `20261028-101059-dc2accd7`, and its Boltz validation was `20261028-101807-b37312d9`. It generated a 75-residue binder; Boltz reported iPTM 0.42641, pTM 0.72618 and mean pLDDT 86.701. The configured policy checked that measurements were present and valid; it set no score threshold. P5 retained that policy and the same PD-L1 chain identity, hotspots `A40,A99,A107`, RFD4 guided-binder task, 60 diffusion steps, batch size 1, one design per attempt, 60–120 residue length range, deterministic mode, native filters enabled, and RFD4 refolding disabled.

P5 used one logical controller campaign: accepted target 2, validated acceptance level, maximum 2 attempts, up to 2 concurrent generation attempts, and one validation at a time. Boltz used 1 sample, 1 recycle, 20 sampling steps, no target MSA or template, and required valid `iptm`, `ptm`, and `plddt_mean` measurements. No threshold was changed to obtain acceptance.

## 6. Actual RTX 5090 inference

Both authorized RFD4 attempts completed on different physical peer workstations. Each native job requested one GPU, used GPU 0, and returned one candidate.

| Attempt | Peer / remote run | Seed | Allocation and result |
|---|---|---:|---|
| 1 | `pc004033`, job `e8dceb04-b682-4fff-b968-442b7ae83898`, run `20261028-124506-6c44be8e` | 1185903269 | GPU 0, 8 CPU cores; one 116-residue candidate; worker interval 11:45:08–11:45:59 UTC |
| 2 | `pc004032`, job `038d8301-e0c7-4380-b729-52e32e105c3c`, run `20261028-124502-4c887566` | 1639588292 | GPU 0, 8 CPU cores; one 94-residue candidate; worker interval 11:45:03–11:45:42 UTC |

The RFD4 invocation records show the native `rfd4-generate` command, one replicate, one-item native batch, and the intended seed and 60-step configuration. The `pc004033` heartbeat sampled 100% GPU utilization during its active job. The `pc004032` heartbeat sampled the allocated GPU at 0% utilization near the end of its job, with about 3.5 GB less free VRAM than idle; its job nevertheless completed on the GPU-backed scheduler allocation and returned a native structure. Both peers reported GPU 0 free again after completion. This distinguishes sampled utilization from scheduler allocation and does not infer continuous utilization from a heartbeat sample.

## 7. Cross-workstation campaign result

The local CPU controller run was `20261028-124448-350f7294`, campaign ID `911d1cc8-4f5a-4a66-8ba1-c46200098805`, cluster submission `5af10b31-163a-45fd-88b0-2f8ceb291fd6`. Its two RFD4 children ran on `pc004033` and `pc004032`, with distinct seeds, worker IDs, cluster job IDs, scheduler allocations, and native output structures. The local controller then ran one Boltz validation for each candidate on `pc002975` GPU 0:

- Attempt 1 validation: `20261028-124743-07cd45e1`, iPTM 0.46578, pTM 0.69548, mean pLDDT 81.473.
- Attempt 2 validation: `20261028-124630-a7d1645d`, iPTM 0.65414, pTM 0.82405, mean pLDDT 93.981.

Both validation records retained `identity_check.sequence_chain_identity=pass` and the corresponding source candidate ID. These values are model measurements under the configured policy, not claims of binding quality.

## 8. Attempt and acceptance accounting

The parent ledger finalized `completed` with termination reason `target_reached` at 2/2 validated candidates. It recorded 2 attempts authorized, 2 started, 2 completed, 0 failed or interrupted, 0 rejected candidate rows, 2 validated accepted candidates, and 2 RFD4 native acceptances unavailable. RFD4 does not emit an authoritative campaign acceptance flag, so those native results were recorded as unavailable rather than inferred from successful process exit. Candidate IDs were unique, and the two binder sequence hashes differed. No attempt was counted twice.

## 9. Structure and sequence identity

For each downloaded RFD4 complex, the parsed chain-A sequence exactly matched the 115-residue input target and the RFD4 candidate metadata; chain B exactly matched its candidate binder sequence. The retrieved RFD4 complex files were:

| Candidate | Binder length | Binder sequence SHA-256 | Complex CIF SHA-256 |
|---|---:|---|---|
| `dc_00000001_d68026ec8a3e62d5` | 116 | `ad6a651791941cbe4fbeebf027fb6c46b470034ebc4c3a2c66d80598285cdf68` | `f773917a1d42cb47ad7cff1d956c4f033c85463a002e67f6c0e7d5ac3a5e3744` |
| `dc_00000002_d68026ec8a3e62d5` | 94 | `da0ba51511da10c21f61dfc4e86e686907a2e733d2e62ae19a39f1ea32e6acee` | `45d0102568de7d1ca3d3de576a3a843ef0b4fcb5b8b98b7da7938aae7b104524` |

Both Boltz output complexes were parsed and their target and binder chains matched the respective RFD4 source sequences exactly. The campaign evidence independently recorded the sequence-chain identity check as `pass`. No binder sequences are reproduced in this report.

## 10. Artifact transfer and recovery

Each remote attempt received the same SHA-addressed input bundle, with 67,676 input bytes recorded. Each peer then completed 11 unique artifact requests; there were no duplicate paths. The received artifact sets summed to 909,738 bytes for `pc004033` and 864,623 bytes for `pc004032`. SHA-256 blob receipts matched the retrieved generated complex CIF, normalized structure CIF, and raw normalized candidate JSONL for both workers.

The campaign recovered both native structures and both candidate rows without resubmitting inference. The repaired transfer ledger persisted command IDs across polling. Earlier duplicate `request.json` requests from a pre-fix local run were drained before this run; they did not alter that run's saved scientific outputs. Interrupted transfers, duplicate completion receipts, disconnection/reconnection, cancellation and late validation are also covered by the focused deterministic regression tests. No active scientific job was stopped to simulate a fault.

## 11. Scheduler and GPU allocation

The controller ran locally as a CPU job with no GPU allocation. The two RFD4 jobs reserved one GPU each on separate peer schedulers, GPU 0 with 8 CPU cores per job. The validation jobs each reserved local GPU 0, one at a time. The controller ledger records 224 cumulative GPU allocation-seconds across the two RFD4 jobs and two Boltz jobs; this is elapsed lease time multiplied by allocated GPU count, not measured kernel time or utilization. Final peer heartbeats showed both 5090s free, and the local scheduler released the validation allocation after each run.

## 12. Execution lifecycle and timings

The complete parent run lasted about 245 seconds, from 11:44:48 to 11:48:53 UTC. Measured RFD4 worker intervals were 51 seconds on `pc004033` and 39 seconds on `pc004032`; their scheduler allocations began at 11:45:08 and 11:45:03 UTC. Boltz worker intervals were 68 seconds and 66 seconds, run serially on the local GPU. Each attempt invoked a separate RFD4 worker/container process; Boltz also ran as separate workflow jobs. Startup, model-load, and sampling phases were not instrumented independently, so the worker intervals include them. No persistent inference worker or native batch larger than one was used.

## 13. CPU and GPU test outcomes

Focused CPU regressions passed during implementation: 147 tests in 39.32 seconds before the final reference-probe follow-up; then 44 campaign and peer-transfer tests in 4.89 seconds, plus 3 targeted cancellation, transfer, and placement tests. `git diff --check` passed. The full broad CPU suite was not rerun after these focused changes. The bounded P5 run itself supplied the real GPU evidence: two successful remote RFD4 inferences on RTX 5090 peers and two local Boltz validation runs on the RTX 4090 workstation. No opt-in pytest GPU suite was substituted for these executions.

The P5 concurrency follow-up also passed 203 focused CPU tests in 16.76 seconds, including a new controller regression that limits dispatch to two validation jobs and keeps the third candidate pending. The real two-GPU Boltz execution is recorded in section 16.

## 14. Remaining limitations and excluded runs

The campaign acceptance policy establishes valid measurements and exact lineage only; it makes no experimental affinity, specificity, or activity claim. Native RFD4 acceptance is unavailable. The `pc004032` heartbeat did not sample high utilization while its allocation was active, although the GPU-backed scheduler job completed and the GPU was released. Container startup and checkpoint loading were not separated from worker wall time. Fault-recovery paths were covered by deterministic tests rather than by interrupting these successful jobs.

The following P5 setup attempts were not counted toward the verified campaign: `de90d8a1` stopped before inference because the campaign actor was missing; `4413088c` and `76014463` failed request-schema checks before inference; `bc8d08b7` ran two RFD4 attempts on the local RTX 4090s and was excluded as a wrong-placement run, with both outputs retained; and `24704453` dispatched a controller to the wrong workstation and was canceled after checks confirmed zero native attempts had started. Its children were canceled before inference. None is evidence of cross-workstation inference. No extra designs or larger batches were run to compensate.

## 15. Commits and changed files

P5 source commits on the deployed branch:

- `724aa75` — dispatch RFD4 campaigns through the peer scheduler.
- `2f60413` — preserve the campaign actor in peer workers.
- `ce00445` — omit local queue metadata from remote RFD4 requests.
- `742ad05` — strip the internal GPU flag from peer requests.
- `b5e9300` — persist peer placement and artifact-transfer command IDs.
- `4f9ba77` — wait for exact peer model-reference probes before dispatch.

P5 source changes are limited to `mn_molecular_design_studio/api.py`, `mn_molecular_design_studio/cluster/candidate_batches.py`, `mn_molecular_design_studio/cluster/context.py`, `mn_molecular_design_studio/engine/workflow_engine/workflows/distributed_campaign.py`, `tests/test_candidate_peer_transfer.py`, `tests/test_cluster.py`, and `tests/test_distributed_campaign.py`. This report is the only documentation file added for P5. The pre-existing local `run_app.py` modification and untracked `tool_references.py` were not staged or committed.

## 16. P5 follow-up — concurrent Boltz validation and four-GPU scheduling

### Deployment and live pool

At the follow-up check, `pc002975`, `4032`, and `4033` were online, accepting jobs, and reported online schedulers. Their live runtime fingerprint remained `9ae11faca1ac4b2ada7e68f2c0c1ddc637dbd4d3a04d66bba0f63c53890f7cda`. The local checkout was at `a72a66a`; the peers remained on runtime source commit `4f9ba772fcb48c1961d74c99f2f5704f34c06197`. The checkout difference is documentation-only relative to the deployed runtime source. This follow-up changed no runtime code, so no remote deployment or service restart was needed or performed.

| Workstation | GPU | UUID | Compute capability | VRAM | Follow-up allocation |
|---|---:|---|---:|---:|---|
| `pc004032` | 0, RTX 5090 | `GPU-b7f8339b-157f-e1f2-f136-bd0ad68bc780` | 12.0 | 31.84 GiB | Idle during local validation; used by P5 RFD4 attempt 2 |
| `pc004033` | 0, RTX 5090 | `GPU-370f52fb-5676-7316-104d-5a9bf016fae1` | 12.0 | 31.84 GiB | Idle during local validation; used by P5 RFD4 attempt 1 |
| `pc002975` | 0, RTX 4090 | `GPU-195a1966-1148-340e-5fe2-bac2d7a9821a` | 8.9 | 23.99 GiB | Candidate `dc_00000002_d68026ec8a3e62d5` |
| `pc002975` | 1, RTX 4090 | `GPU-c0569fbb-42d6-594e-a304-0a84a0c60a21` | 8.9 | 23.99 GiB | Candidate `dc_00000001_d68026ec8a3e62d5` |

All machines reported CUDA driver API 13020. At final inspection, each scheduler had zero active reservations; the local GPUs were idle. The four-GPU pool was exercised across the distributed generation and validation phases: the two RTX 5090 peers ran the existing RFD4 campaign, and the two local RTX 4090s ran the concurrent validation. The four devices were not all active at the same time.

### Concurrent validation results

The original P5 campaign remains unchanged: its ledger SHA-256 is `205efddd4c18cea8f8ca9d9ee22a9b5ad193964f5e059854631bf2465e891b39`, matching before and after the follow-up. A separate validation-only ledger, campaign `af8ce7ac-b56f-4610-a26c-ec0765b2f7ec`, reused the same two retained candidates and PD-L1 chain-A target. It dispatched with `max_concurrent_validations=2`, `gpu_device=auto`, one GPU per job, 8 CPU cores, and 16 GiB minimum VRAM. No additional RFD4 generation ran and no global default changed.

Boltz-2 used the P5 settings: one sample, one recycle, 20 sampling steps, optimization off, no MSA or target template, and required `iptm`, `ptm`, and `plddt_mean` measurements. Both jobs used `mn-boltz2-opt:2.2.1-cu130-batch2`, image ID `sha256:3054688e91908156071238adee167ad4963cfbca05c8a4cf0e55a634ba5a2a55`. The checkpoint hashes were `boltz2_conf.ckpt` `090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1`, `boltz2_aff.ckpt` `dcc5cd3722b1c9eaa34267e4ae32f55cbbf1963f4c19319381ccfa30fdd2ca9e`, and `ccd.pkl` `2d3b2f03a3c5665944adba51e33263511e51b21c9cd05d902f9c4b7c1e58d2f4`. The follow-up inputs retained the target file with SHA-256 `7da1c888011cdbd0937e3ff6b73429aa069e4c496cfa0cb8d8b53cc06ca73b57` and the two checksum-verified P5 complexes from section 9.

| Source candidate | Follow-up Boltz run | GPU | Seed | iPTM | pTM | Mean pLDDT | Identity / lineage |
|---|---|---:|---:|---:|---:|---:|---|
| `dc_00000001_d68026ec8a3e62d5`, 116 residues | `20261028-132737-b243cd10` | 1 | 1203626143 | 0.64732 | 0.77642 | 80.499 | Sequence-chain and source-candidate checks passed |
| `dc_00000002_d68026ec8a3e62d5`, 94 residues | `20261028-132734-d4512561` | 0 | 1801939714 | 0.43823 | 0.72855 | 91.590 | Sequence-chain and source-candidate checks passed |

Both scheduler allocations overlapped from 12:27:38 to 12:29:01 UTC, about 83 seconds. A 145-sample `nvidia-smi` log captured four samples with both local GPUs above zero utilization; the first such sample showed 4% on GPU 0 and 38% on GPU 1, with 3,152 MiB and 4,166 MiB in use. The output metrics differ from the original P5 predictions because these are separate seeded inferences; the repeat results do not replace or modify the original campaign evidence. As in P5, the acceptance policy checks measurement presence and validity, not score thresholds or experimental binding.

The follow-up ledger ended `completed/target_reached` with two attempts completed, two candidates validated accepted, and zero rejected or interrupted attempts. Replaying the controller after completion was idempotent. A first isolated harness pass did complete two Boltz predictions on distinct GPUs, but its parent run lacked the standard normalized-candidate index, so the strict lineage reconciler rejected those outputs. That pass is marked failed and preserved separately; it is excluded from acceptance counts. The corrected pass above included the index and passed both lineage checks. No original P5 results were overwritten, and no active job was interrupted.

The full follow-up evidence is in `/mnt/data/RESULTS/mn-molecular-design-studio/verification/phase6b-p5-concurrent-validation/20261028T122730Z-f047554f`, including source snapshots, checksum manifest, request files, execution summary, scheduler metadata, and GPU samples. The initial harness pass and its failure disposition are preserved under `/mnt/data/RESULTS/mn-molecular-design-studio/verification/phase6b-p5-concurrent-validation/20261028T122312Z-06e45a21`.

## 17. P6 interface and retained-results verification

The Results interface now projects the existing campaign ledger and normalized candidate index without changing either record. It shows global campaign acceptance separately from valid prediction measurements, and labels missing index or lineage evidence as unresolved metadata rather than a scientific candidate rejection. The first isolated validation harness is therefore visible as an index/lineage failure with no current-run measurements; the corrected follow-up remains a separate validation-only campaign linked to the original P5 campaign. The panel does not replace original measurements with follow-up values.

The campaign panel also reports controller dispatch batch size and dispatch capacity separately from the generator's native inference `batch_size`. It shows generation dispatch capacity separately from downstream validation concurrency, and includes attempt and validation worker/run/GPU provenance. For the retained P5 records it reads controller batch size 1, two concurrent generation batches (capacity 2), native RFD4 inference batch size 1, and the follow-up's two concurrent validation slots. The original P5 campaign used one validation slot.

Workstation settings now show detected GPUs, currently available GPUs, and currently allocated GPUs from the latest scheduler snapshot, with per-device model and scheduler state. Offline nodes are labeled as last reported. The inventory fixture was checked against the verified pool above: `pc002975` has two RTX 4090s, `pc004032` one RTX 5090, and `pc004033` one RTX 5090. The interface does not claim all four GPUs are active at once.

The retained PD-L1 Results views rendered successfully for P4 campaign `686e0a0a-54fb-4468-b211-d1d2159c78fa` (one accepted candidate and one valid validation measurement), original P5 campaign `911d1cc8-4f5a-4a66-8ba1-c46200098805` (two accepted candidates and two valid measurements), validation-only follow-up `af8ce7ac-b56f-4610-a26c-ec0765b2f7ec` (two accepted candidates and two valid repeat measurements on local GPUs 0 and 1), and the preserved index-missing harness `b8428038-8533-476a-acd3-30c47662ce89` (two unresolved metadata outcomes, zero scientific rejections, and zero current-run measurements). The native inference settings and candidate rows were read from retained files; no generation or validation inference was rerun for UI verification.

Verification: `tests/test_distributed_campaign_ui.py` plus `tests/test_cluster_ui.py` passed, 18 tests total. `py_compile` passed for the changed Python modules. Streamlit AppTest rendered all four retained Results contexts and the four-device workstation inventory. The Results projection is read-only. No source deployment, worker change, generation capability, saved pipeline, scientific threshold, accepted/rejected record, or worker provenance was changed.
