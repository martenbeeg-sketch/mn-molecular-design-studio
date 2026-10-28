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

## 16. P6 readiness

P5 verified that the existing campaign ledger can place two distinct RFD4 attempts on compatible enrolled peers, retain one shared attempt budget and seeds, recover and checksum native artifacts, validate those exact sequences downstream, and finalize after late validation. The backend is ready for P6 to expose these existing worker, validation, accounting, and transfer states. P6 UI work was not part of this verification and remains untested.
