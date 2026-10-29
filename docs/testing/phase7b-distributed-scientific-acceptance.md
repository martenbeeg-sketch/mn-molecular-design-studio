# Phase 7B — bounded live execution completed; partial acceptance

> **29 October follow-up:** The earlier partial-acceptance summary below describes the previous run. A later complete native run produced a genuine passing candidate and completed independent comparison. The active three-workstation revision is now `e1ebb2c`, fingerprint `dd1c7388c8426c1045f381064081ecf628af5ad2cc2abbdd1225b16635ec57a7`. See section 32 for the new evidence and remaining limitations.

Updated 29 October 2026, Europe/Berlin. **Final status: bounded live execution and Protenix RTX 5090 compatibility recovery completed; overall acceptance remains partial.**

All three workstations now run corrected revision `32609863cfb99faebb7c3dd5b1c4ba97acb0f9f0`, with matching active fingerprint `24de4cb1fdff303ec0fbd092f78af1d1222e5ed19b6c82a0db322021d3ceddb7`. The 28 October scientific run used `2edf8de`; its results remain preserved. BindCraft2 completed one trajectory with a genuine screen rejection and zero eligible candidates. RFD4 produced one fresh candidate; it and the two retained P5 candidates each now have independent Boltz-2 and Protenix V2 results. Seven predictor input attempts produced six successful predictions; one Protenix failure on RTX 5090 was recovered on RTX 4090 without repeating completed inference. Results/ZIP and exact sequence/target identity passed verification.

All four physical GPUs executed real work across the run. **Four simultaneous native inference processes were not demonstrated.** Positive complete-native comparison and downstream internal BindCraft stages remain unverified because the sole trajectory was rejected. The original Protenix V2 stock image failed on the 5090. The 29 October image fix now passes real stock inference on 4090 and 5090, but the resulting poses differ substantially and numerical equivalence is not certified. The remaining scientific and concurrency limits prevent an unqualified full acceptance claim.

The user authorized correcting problems, rebuilding the faulty image and deploying working fixes to the other machines. The 28 October corrections and 29 October image selection/preflight correction were tested, committed locally and deployed through the reviewed idle-only service override/rollback pattern. Original checkouts, local edits, old images and 27 protected scientific files remain unchanged; no shared scheduler/tunnel service, model weight, scientific threshold or generation budget was changed. No GitHub push or merge occurred. Sections 1–28 preserve the 28 October history; **sections 29–31 record the authorized 29 October image correction and current status**. Portable summary: [phase7b-live-acceptance.json](evidence/phase7b-live-acceptance.json).

Evidence labels used below are **Real GPU execution**, **Retained scientific artifact**, **Deterministic execution fixture**, **Source inspection**, and **Unverified**. CPU-only container imports/probes are Source inspection, not Real GPU execution.

Raw evidence: `/mnt/data/RESULTS/mn-molecular-design-studio/verification/phase7b-distributed-acceptance`. Portable machine summary: [phase7b-preflight.json](evidence/phase7b-preflight.json).

## 1. Exact source and GitHub baseline

**Source inspection.** All three original checkouts were already at `abf0f2ef5b3d1117df72c6a74a56e05bb3f21f2f`, the Phase 7A.2 commit “feat: add three scientific workflow modes with distributed branch execution.” Both peers had moved beyond historical `4f9ba77`; this was observed over SSH, not assumed. `git ls-remote origin refs/heads/fix/scientific-pipeline-validation` returned the exact same full commit. No push was needed or performed.

The original local branch remains `fix/scientific-pipeline-validation`. All three original checkouts retain their modified `run_app.py` and untracked `tool_references.py`; their file hashes match across machines. No original checkout was reset, stashed, rebased, merged or updated. The four requested architecture/testing documents were read.

Fixes were reviewed and committed separately as `a628c33becb3258df39e51ff89cf10643649e59d` on isolated branch `codex/phase7b-acceptance`. No code or documentation from this task was pushed to GitHub. There was no merge into main.

## 2. Observed runtime and staged deployment

**Source inspection.** Existing API/agent services were active on all three machines, using the existing Studio Python 3.11.17 environments. Published heartbeat source fingerprints matched the original source trees, including the preserved paper-reference edits:

`5569f8a69023be09c6079d2608fb383d60ba7a1dbcca1bd29e96c92c7fa53d5b`.

This is source-on-disk and heartbeat evidence. Local API/agent processes predate the latest code; exact in-memory module revisions are **Unverified**. A matching heartbeat alone is not proof that every already-imported Python module was refreshed.

| Host | Original checkout HEAD | Clean staged HEAD | Activation |
| --- | --- | --- | --- |
| pc002975 | `abf0f2e` | `a628c33` | Not performed |
| pc004032 | `abf0f2e` | `a628c33` | Not performed |
| pc004033 | `abf0f2e` | `a628c33` | Not performed |

Local pinned source is `/home/user/programs/deployments/mn-molecular-design-studio/phase7b/source`; each peer has `/home/compute/programs/deployments/mn-molecular-design-studio/phase7b/source`. All three checkouts are clean, detached at the full staged commit and import the expected package path. Their common clean runtime fingerprint is `96814c66c8f1afdf5a7beae74a30f15f6452fddf32470e930bd73e4abc8f92fe`. The full-history source bundle was verified before the two remote clean clones. This is staged deployment, not an active-runtime claim.

The existing environments are `/home/user/mambaforge/envs/mn-molecular-design-studio` and `/home/compute/programs/miniforge3/envs/mn-molecular-design-studio` on the peers. No global installation or environment replacement was performed.

## 3. Activation boundary and rollback

**Source inspection; activation Unverified.** Reproducing the retained native configuration required a graph handoff fix. The new source must be activated consistently before fresh acceptance. Seven host-specific service drop-ins and an exact activation/rollback plan are prepared under `deployment/` in the evidence root. None was installed.

The approval action is to switch API/agent services on all three machines plus the controller UI to the pinned clean source, and restart only idle Studio application scheduler bridges. The existing shared MN Compute Scheduler and SSH tunnel services remain in use. Merely refreshing the API would leave old bridge processes spawning workers from the old import root.

The request's section 3 explicitly requires stopping before service modifications. Approval was requested; elapsed time is not approval. Activation must recheck Studio queues and leases first and must not interrupt useful scientific work. It causes a short UI/API reconnect and dispatch pause. It includes no model download, image rebuild/replacement or GitHub push.

Rollback removes only the newly introduced `phase7b.conf` drop-ins, restores the previous service environment, and restarts idle application bridges from the original source. Existing units/drop-ins must be backed up before any change. The original three checkouts, settings, credentials, images, weights and Results remain available; rollback requires no destructive Git action. Full procedure: `/mnt/data/RESULTS/mn-molecular-design-studio/verification/phase7b-distributed-acceptance/deployment/activation-plan.md`.

## 4. Four physical GPUs and live availability

**Source inspection.** UUIDs and telemetry came from `nvidia-smi` and the existing scheduler inventory.

| Host | Index | GPU | UUID | VRAM |
| --- | --- | --- | --- | --- |
| pc002975 | 0 | NVIDIA GeForce RTX 4090 | `GPU-195a1966-1148-340e-5fe2-bac2d7a9821a` | 23.99 GiB |
| pc002975 | 1 | NVIDIA GeForce RTX 4090 | `GPU-c0569fbb-42d6-594e-a304-0a84a0c60a21` | 23.99 GiB |
| pc004032 | 0 | NVIDIA GeForce RTX 5090 | `GPU-b7f8339b-157f-e1f2-f136-bd0ad68bc780` | 31.84 GiB |
| pc004033 | 0 | NVIDIA GeForce RTX 5090 | `GPU-370f52fb-5676-7316-104d-5a9bf016fae1` | 31.84 GiB |

All three schedulers were online and their snapshots showed zero GPU reservations at inspection. The local scheduler's global counts included 388 queued jobs, but a separate app-filtered query returned **zero queued/running Molecular Design Studio jobs**. Those global jobs were left untouched. GPU 1 had some external/display memory use and nonzero sampled utilization; a free Studio lease does not prove the physical GPU is otherwise idle. Exact memory/utilization snapshots are in `preflight/*-models.json`.

Four registered devices are not four active scientific allocations. New Phase 7B leases, UUID-linked process receipts and simultaneous four-device intervals are **Unverified**.

## 5. Authenticated workers and artifact transport

**Source inspection.** Authenticated read-only peer API requests succeeded to both configured peers using existing credentials; no credential was printed or changed. All three agents were online, accepting jobs, reporting accessible storage, and using the existing scheduler and transfer system. Snapshot timestamps and resource observations are preserved in raw JSON.

Results/reference roots and preserved service/configuration hashes are recorded per machine. Source bundles were transferred by SSH/SCP into new deployment directories. That transfer is not a scientific artifact-transfer acceptance test. Real fresh candidate upload/download, checksum receipts, process lifecycle transfer and live remote cancellation remain **Unverified**; retained scientific lineage and deterministic transfer/recovery coverage are recorded separately below.

## 6. Images and checkpoints

**Source inspection.** Existing images were inspected without pull/rebuild. The following are exact local IDs and portable content fingerprints; peer image IDs differ because their Docker image stores represent identity differently. The portable content fingerprints match on all three workstations, including all five rows.

| Image | pc002975 image ID | Portable content fingerprint |
| --- | --- | --- |
| `mn-bindcraft2:v1.0.4-cu13-pae-v1` | `sha256:673ec7cf961eb892f6977ed425920a2fae22eb0b0aae4510b9e7cbecc7fedbcb` | `sha256:14782c560a9c20271119332d22a497d4f88328cd6842261f809ae17d283b9113` |
| `mn-rfd4-proteina:0.0.3-cu130` | `sha256:c34efb5afca78c0341cf996c9dadc8102a7fb3760e1792697a8779803e0e44f5` | `sha256:a0b8ed6e9c418b6ee7435ff633d77f1248abdab5b8d7b249b18b3a2dc5d3f82a` |
| `mn-boltz2-opt:2.2.1-cu130-batch2` | `sha256:3054688e91908156071238adee167ad4963cfbca05c8a4cf0e55a634ba5a2a55` | `sha256:10b16971ffc0567920260f71b45dd5fe0041b9fbce2a9f39655bd92bad852f12` |
| `mn-protenix-opt:2.0.0-cu128-batch2` | `sha256:f9b86d5f136581795d0fd4ba968e6ae37a5fa256a63eebec5d896fe69182f5ca` | `sha256:5f3881201c380dc823cf4abb924c586b1dbdd52eb7f30cefc332b6d8fb2cb33b` |
| `mn-ligandmpnn:latest` | `sha256:7f85a354ebe11a8edf2722aa3cf2247189ff74f62985d4ac9e72e8a9a40e181d` | `sha256:1d8a6644e5c449dd6e164fe047f1cf7027e30a77afcfe06dc65bc1beb3503c96` |

Selected checkpoint/config hashes match on all three machines:

| Relative reference artifact | SHA-256 |
| --- | --- |
| `rfd4-proteina/huggingface/hub/models--nvidia--RFD4-Proteina/snapshots/a06e2509c7033b010fc43614e9fd840303588502/config.json` | `154ea927ab8005b7eac2d047a792147fb4606ab4878e496aa8414b6b40dad368` |
| `rfd4-proteina/huggingface/hub/models--nvidia--RFD4-Proteina/snapshots/a06e2509c7033b010fc43614e9fd840303588502/model.safetensors` | `5b549dd61adb2f6c6b5cb95bda13f40f5e31556b4e6c755edc0ca6f011c68fb0` |
| `boltz_models/boltz2_conf.ckpt` | `090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1` |
| `boltz_models/boltz2_aff.ckpt` | `dcc5cd3722b1c9eaa34267e4ae32f55cbbf1963f4c19319381ccfa30fdd2ca9e` |
| `boltz_models/ccd.pkl` | `2d3b2f03a3c5665944adba51e33263511e51b21c9cd05d902f9c4b7c1e58d2f4` |
| `protenix/checkpoint/protenix-v2.pt` | `8f931f9774a396b67033d0e58628e1834f4a1448165e04254b40a780b0c0d599` |
| `protenix/common/components.cif.rdkit_mol.pkl` | `d1cfb71f5993a3ebea7c47877022d7f597bbfbaf86e28a4770e957da6c50cd35` |
| `alphafold_models/params_model_1_multimer_v3.npz` | `611da8fc7478928f68de12e8b226260ef1f4ce62bcc29b008572e52f4f212959` |
| `alphafold_models/params_model_1_ptm.npz` | `5e564f79af5bcd54ccef6e2a6bb0ff01015d01650ebc41d4575e35f0de9ecc84` |
| `mpnn/v48_020_epoch300_hyper.pt` | `2635afa506cbad850af3fffa7629d80b72ee166d386e2d63a05435b5bd9b459f` |

The staged compatibility probe reports RFD4, Boltz-2 and Protenix V2 compatible with the four observed GPUs and driver APIs. Protenix reports PyTorch `2.7.1+cu128` and targets `sm_75 sm_80 sm_86 sm_90 sm_100 sm_120 compute_120`; this is a CPU import/compiled-target inspection, not inference. The original generic probe failed before importing torch because it did not reproduce `mn-predict`'s existing removal of legacy `PROTENIX_OPT_MODE`. The fix follows that launcher exactly, preserving real driver/architecture checks.

The torch probe does not establish BindCraft2's JAX runtime compatibility. Prior real native 4090 evidence is retained; fresh native placement still requires its ordinary preflight and process evidence. No assumption is made that every installed model runs on every GPU.

## 7. Integration fixes and review

**Source inspection + Deterministic execution fixture.** The complete BindCraft2 graph tile previously rejected the retained benchmark's `engine_overrides` and native worker controls, despite standalone support. It also coupled the accepted-design target and trajectory cap. The fix exposes a separate native trajectory limit, native workers/GPU and an advanced native-settings JSON object in the existing form. It maps these into the same complete `design.bindcraft2` request. A default zero trajectory override preserves old saved graphs' count-derived cap; auto worker/default settings remain omitted from the native request.

Advanced JSON is validated before submission and cannot change trajectory limits, accepted-design budget, workers/GPU or full/generator mode. Nonfinite/invalid JSON fails. The native operation remains one whole invocation. The selected complete tile continues to use the adapter's `trajectory_only=False` default. No scientific thresholds, algorithms, models or campaign packing limits were changed. A second fix makes the Protenix CPU compatibility probe use its installed launcher's startup convention.

Code and tests were reviewed using the staged diff. Validation includes saved/reopened native configuration, request preservation, backward compatibility, scientific cache-key changes, invalid-budget settings, the actual Streamlit tile form, and driver/architecture constraints. No runtime compatibility bypass was added.

## 8. PD-L1 reference and retained lineage

**Retained scientific artifact.** The target is `/mnt/data/RESULTS/mn-molecular-design-studio/workdir/runs/design/20261028-124448-350f7294/artifacts/rfd4_proteina/input/PD-L1.pdb.cif`.

Target SHA-256: `7da1c888011cdbd0937e3ff6b73429aa069e4c496cfa0cb8d8b53cc06ca73b57`. Chain A is 115 residues; sequence SHA-256 `c434254957dd7a7d80ce4ca859f9a6d101a1205f4ff822177dd0b3c31f02e968`. RFD4 hotspots are A40, A99, A107. Author residue mapping is stored in `retained/lineage.json`.

P4 ledger `20261028-101057-a07ef0a4`: `0f575ca3e563181423524207c34ebc5fbcf1c6c9b6ff8c5bf71409ae10bc1460`.
P5 ledger `20261028-124448-350f7294`: `205efddd4c18cea8f8ca9d9ee22a9b5ad193964f5e059854631bf2465e891b39`.

P4 predates the external Results-root migration. Its original run is under `/home/user/.local/share/mn-molecular-design-studio/workdir/runs/design/20261028-101057-a07ef0a4`; portable references were resolved against each recorded source root, following the existing retained-import convention. No original paths/files were rewritten. Original ledgers, candidate records, source inputs/results and available source/prediction structures were hashed before/after replay. All tracked hashes match. Parsed target/binder chains exactly matched the retained sequences for all three P4/P5 candidates (75/116/94 residues).

## 9. Mode 1 — complete native workflow

**Prepared configuration, Source inspection.** Target chain A → complete `bindcraft2` → native outcomes → Native Results; no external ProteinMPNN, Boltz or Protenix nodes. The saved bounded configuration uses binder length 60, native seed 20261027, one accepted-design stopping target, **one trajectory maximum**, one native worker/GPU and the original benchmark's one-step native stage overrides, model_1_multimer_v3 design/model_1_ptm validation, one sequence candidate, native filtering and failed-trajectory preservation. The retained one-trajectory cap is stricter than the request's maximum two. No generation was repeated to seek an accepted binder.

**Retained scientific artifact.** Original native run `20261028-085821-36356a14` records two historical rejected trajectories and zero normalized candidates. Its imported outcome ran through the normal graph and exported successfully, without external predictions. This establishes rejection/empty routing, not a fresh native execution. Actual gradient design/internal sequence design/AF2 evaluation/filter process receipts for a new Phase 7B invocation are **Unverified**.

## 10. Mode 2 — native workflow plus comparative evaluation

**Retained scientific artifact.** The genuine complete native rejection source routed automatically through CandidateSet, both enabled empty predictor branches, analyzer, ranking and Results. Both predictors recorded zero processing jobs. No RFD4 co-design source was relabeled as a complete internally evaluated native workflow.

**Unverified.** A positive complete-workflow candidate was not available in the inspected benchmark, and no fresh Mode 1 result exists yet. After activation, reuse the genuine Mode 1 output if eligible and replace the prepared Mode 2 retained reference accordingly. Preserve its original sequences/native filters and use independently configured Boltz-2 and Protenix V2 branches; no mandatory redesign. Positive Mode 2 contracts are deterministic fixtures only. The positive scientific acceptance boundary remains explicitly incomplete.

## 11. Mode 3 — modular pipeline

**Prepared configuration, Source inspection.** One fresh RFD4 guided binder request uses the retained P5 chain/hotspots, 60 steps, one replicate/design, native microbatch one, 60–120 aa length range, seed 20261028, deterministic/native-condition filtering enabled and native refolding disabled. The graph label `guided` and P5 label `guided_binder` enter the same existing native preparation branch. A second explicit retained RFD4 branch supplies the two original P5 candidates, without another generator launch.

Initial ProteinMPNN is conditional and sequence-complete RFD4 outputs bypass it. Optional redesign remains disabled for this prepared live graph. Its disabled path is covered by retained execution; interface-fixed child identity and parent-evidence isolation are covered by deterministic tests. No new redesign variant was generated.

Both predictors are enabled: same exact target/binder sequences, no MSA or target template, seed 42, one sample. Boltz retains P5's one recycle/20 sampling steps; Protenix retains its ordinary 10 recycle/200 step profile. Ranking explicitly uses `prediction/predictor-0/iptm` within recorded protocol cohorts; no threshold or zero fill is introduced, and no numerical equivalence to a Protenix measurement is assumed.

**Retained scientific artifact.** The two P5 variants bypassed sequence inference, joined their original Boltz evidence (iPTM 0.46578 and 0.65414), produced consolidated Results/selection and ZIP export. Protenix is explicitly disabled in this retained replay. RFD4 native acceptance remains unavailable. Retained campaign validated acceptance records measurement presence and lineage, not a native filter pass or experimental binding. Fresh RFD4 and both fresh predictor branches are **Unverified** pending activation.

## 12. Saved graph configurations and normal execution

**Source inspection.** Prepared graph definitions were saved/versioned through Studio's normal graph API, with validation and CandidateSet preflight. They submit no jobs. These are snapshots for opening in the existing editor after consistent activation:

| Mode | Graph ID | Saved revision |
| --- | --- | --- |
| native_complete | `9c6b0015ed11428f8677991122532249` | `b78b23a7b869fb32e4815bc107bc6097037d2683fc5e11f9ba27a3a9d24d0342` |
| native_comparative | `264de98ded9d4edaa608043289febeee` | `019c43621405cdd411a23357b7bb91cdd77c30a2dddd9d22aa643253bc65015b` |
| modular_design | `2c1cb60b631d44ccbe8e45109b780d07` | `f474b3a1129367e3ab3050145e273f0392cc9f953bfacabb15363a215bd388a9` |

Owner is the existing development user; project is `37130be6726a40d1a37652b0d7106c37` (Phase 6A deployment verification). The existing browser displayed the three-mode setup and a four-node complete-native graph. No live Run button was used. Native input editing, plan review, per-stage/GPU monitoring and final fresh Results interaction are not claimed complete. Streamlit setup/Results contracts and the new native form are covered by deterministic execution fixtures.

## 13. CandidateSet fan-out/fan-in and batch accounting

**Retained scientific artifact + Deterministic execution fixture.** Candidate normalization, exact variant IDs, original sequences and native provenance survive retained fan-in and analysis. The modular replay has two original sequence-ready candidates, zero redesign variants, two retained Boltz predictions, two analyzed rows and two selected rows. Native/filter absence remains explicit. Empty native sources retain their two rejection attempts even with zero normalized rows.

Planned fresh bound is at most three eligible candidates per predictor, at most **six candidate–predictor inputs total**, one seed/sample each. Retained Boltz candidates are grouped in one two-input chunk; retained Protenix candidates use two independent one-input chunks. This supplies useful work for placement while also exercising one actual multi-input interface when run. Actual scheduler jobs, containers, native invocations, model initializations, samples, returned predictions and accepted candidates must be read from the new receipts; forecast chunk counts do not establish them.

No campaign `native_attempts_per_job` larger than one was enabled. Deterministic collection/recovery tests preserve the existing campaign boundary and native multi-input behavior. No new live native batching claim is made.

## 14. Distributed parallel scheduling

**Deterministic execution fixture.** Existing tests exercise independent generator submission, downstream prediction while another generator remains active, four independent in-flight chunks, stale/unavailable machines, shared lease accounting and strict receipt reconciliation. Existing authenticated worker/image/checkpoint eligibility remains authoritative.

**Unverified.** No new physical lease or overlapping scientific allocation was created. The planned native and modular graphs expose useful independent tasks, but neither four-device occupancy nor placement of any specific new model is promised. Capture UUIDs, lease IDs, allocation/start/finish timestamps and actual process receipts during the authorized live continuation. Prior P5 evidence used all four devices across phases but did not show four simultaneous devices; it is not relabeled as Phase 7B concurrency.

## 15. Provenance and scientific identity

**Retained scientific artifact.** The read-only lineage check parses original complexes and verifies exact chain-A target and binder sequences. The normal graph preserves source/native candidate IDs, backbone IDs, sequence variant IDs, original sequence/reference structures, target hash, predictor ancestry and independent measurement namespaces. Missing historical trajectory/model metadata remains unavailable.

The analyzer's same-variant/same-sequence joins, distinct redesign child IDs, missing predictor entries, native rejection versus operational failure and parent-evidence isolation are exercised in deterministic tests. Complete fresh model/settings/seed/GPU/native receipts and remote artifact hashes remain **Unverified**. No confidence score is described as experimentally measured binding.

## 16. Results and exports

**Retained scientific artifact.** Three new lightweight retained graph runs completed through the normal engine and `build_run_zip`; exports contain CandidateSet workflow summary and native-source provenance. ZIP identities are:

| Mode | Retained graph run | Candidates | ZIP SHA-256 |
| --- | --- | --- | --- |
| native_complete | `20261028-185830-26e7d024` | 0 | `9177afa19956270cccf7f962556108d5e12643b2c276f8177b2c92cc19689018` |
| native_comparative | `20261028-185830-0eb629d1` | 0 | `1e3049a92bd96a013a97ea1a407d97a0b2e635372a786d09a183248d6dcd8ad5` |
| modular_design | `20261028-185831-33dafdbb` | 2 | `a0732e1c9401f0b82944984e9b82917fe45bd123a6bc40cafef0391cf518260b` |

ZIP paths are recorded in `retained/summary.json`. The modular ZIP includes consolidated measurements/selection from genuine retained Boltz outputs. Empty native mode ZIPs preserve source/rejection evidence. Fresh per-GPU progress, real partial failures, positive complete-native comparison Results and fresh scientific ZIP export remain **Unverified**. Analyzer-only upstream cache reuse is fixture/source evidence, not a new live inference claim.

## 17. Recovery and failure coverage

**Deterministic execution fixture.** Focused suites cover unavailable/disconnected workers; failed/retried submissions; corrupted/interrupted output transfer; partial native output and missing artifacts; duplicate/conflicting receipts; failed prediction siblings with preserved successful evidence; restart/configuration checks; completed-member preservation; queued cancellation and no-launch fences. Representative contracts are in `test_candidate_execution.py`, `test_candidate_peer_transfer.py`, `test_execution_job_ledger.py`, `test_workflow_modes.py`, `test_native_collection_execution.py`, `test_cluster.py` and `test_peer_workstations.py`.

No useful scientific inference was killed. Actual fresh remote cancellation/transfer interruption is **Unverified**. Replaying completed receipt records is not treated as a scientifically equivalent restart of a stochastic generator.

## 18. Runtime and utilization measurements

**Source inspection.** Captured GPU memory/utilization samples, scheduler resource snapshots, model/file sizes and exact source/image/checkpoint hashes are available. CPU-only container imports yielded compiled-target profiles. Test suite wall times below are software-test measurements, not model timings.

**Unverified.** New scientific queue time, leases, native process durations, model initialization, candidate/transfer/prediction/analyzer/end-to-end durations and concurrent GPU utilization were not measured because no new inference ran. Retained graphs' lightweight completion times are not scientific performance data. No startup/batching speedup is claimed from fewer scheduler jobs. Phase 6C.1's older measured RFD4 import/load phases retain their original scope.

## 19. Validation record and scientific outcome

**Deterministic execution fixture** (overlapping suites, not a sum of unique tests):

- Baseline: **163 passed in 233.99 s** — modes, CandidateSet, candidate execution/peer transfer, distributed campaigns and job ledger.
- Native-settings fix: **219 passed in 208.52 s** — new configuration tests, workflow modes, graph, candidate execution/peer transfer, native collections and ledger; before the separate Protenix probe follow-up.
- Compatibility follow-up: **85 passed in 10.81 s** — configuration/probe, peer-workstation and cluster tests.
- Final native configuration/UI: **15 passed in 0.53 s** — includes actual Streamlit form editing and apply.
- Source diff whitespace validation passed. Staged source imported at the exact pinned revision on all three machines. CPU-only probes used existing installed images.

Initial harness setup mistakes (a nonexistent test filename, pre-workspace path resolution, and unsupported raw graph settings) failed before inference and were corrected. They are not counted as scientific failures or successful live runs. The native graph rejection exposed the configuration gap fixed here. Earlier retained replay artifacts are preserved separately; the final summary points to the replay with explicit original-structure checks.

**Software outcome:** staged source, focused regressions, retained graph routing/lineage and retained export pass. Active runtime switch and fresh distributed execution are incomplete.

**Scientific outcome:** no new candidate or prediction. Genuine retained BindCraft2 rejection remains rejection; retained P4/P5 RFD4/Boltz evidence remains available. Positive-looking fixtures are never accepted as scientific evidence.

## 20. Limitations, remaining priorities and continuation

Phase A is complete for source/staging checks but active-runtime attestation awaits activation. Phase B retained checks pass. Phases C/D/E are not executed. Phase F passes retained export and deterministic recovery contracts; fresh UI/scientific recovery remains incomplete. Phase 7B's definition of done is **not met**.

The next required action is the separately approved, idle-only service/bridge switch in section 3. After activation, verify common source fingerprints and exact model/reference probes before using the normal Studio Run action. Run one bounded complete BindCraft2 invocation, then reuse its genuine outcome in Mode 2 where eligible. Run the prepared RFD4 modular graph with automatic downstream prediction and the six-input ceiling. Measure real physical leases/process receipts and four-device overlap where compatible work exists. Finish fresh Results/provenance/ZIP verification and safe bounded cancellation. Do not expand generation budgets to force acceptance or device occupancy.

If Mode 1 again yields no eligible complete-workflow candidate, preserve its native rejection and leave Mode 2's positive path explicitly incomplete unless an independently provenance-verified complete workflow source is available. RFD4 co-design does not satisfy that boundary. Models/images may not be replaced or services further modified beyond the reviewed activation scope without the request's stated approval boundary.


## 21. Approved activation attempt and verified rollback

### Current outcome

**Software execution failure before scheduler acceptance; scientific outcomes unverified.** The approved activation of `a628c33becb3258df39e51ff89cf10643649e59d` succeeded operationally on all three machines. A real submission of the prepared complete-native graph then failed with HTTP 422 before any scheduler job or model inference started. The native settings string was treated as a filename. Because correcting that defect changes the explicitly approved runtime revision, acceptance stopped and the reviewed rollback was executed. The staged source was not modified in place.

| Acceptance item | Actual evidence in this continuation | Status |
| --- | --- | --- |
| Activate API/agents/controller UI and idle application bridges | Seven services switched; exact process import paths and fresh common fingerprint verified | Performed, then rolled back |
| Shared scheduler and SSH tunnels | Existing infrastructure used; no configuration/restart applied to either | Preserved |
| Mode 1 fresh complete BindCraft2 | Normal Studio API/graph submission attempted; HTTP 422 before scheduler acceptance | Operational failure; no scientific invocation |
| Mode 2 positive complete-native comparison | No fresh Mode 1 candidates exist | Unverified; no RFD4 substitution |
| Mode 3 fresh modular workflow | Prepared graph passed normal API validation; not submitted after the rollback boundary was reached | Unverified |
| Four-GPU scientific work / actual batching | 348 host observations; zero running Studio jobs and zero leases | Unverified, not inferred from registration |
| Original P4/P5 and native benchmark artifacts | All 27 previously tracked source hashes rechecked after activation/rollback | Unchanged |
| Original runtime health | API health, authenticated peers, original fingerprints, service unit hashes and idle queues verified | Restored |

### Pre-activation safeguards and actual service actions

**Source inspection.** Immediately before activation, each scheduler had zero reservations and no queued/running Studio jobs. No scientific Docker container was active. The four UUIDs remained the ones in section 4. Local GPU 1 had approximately 1.5 GiB of external/display memory use and sampled 0% utilization at this check; no Studio allocation was inferred from that memory.

The seven prepared service drop-ins were installed only after private backups of original units and existing drop-in directories. API/agent services were switched on pc002975, pc004032 and pc004033; the controller UI was switched on pc002975. Idle application scheduler bridges were refreshed. The exact clean source root and pinned HEAD were checked before installation. The original Git checkouts, paper-reference modifications, installed images, checkpoint stores, Results roots, shared compute scheduler and SSH tunnels were preserved.

The three APIs returned `status=ok`; fresh authenticated heartbeats reported `96814c66c8f1afdf5a7beae74a30f15f6452fddf32470e930bd73e4abc8f92fe` on all three machines. `/proc` verification established that each affected service and application bridge used the pinned source through `PYTHONPATH`, rather than relying only on source-on-disk hashes. The restarted Studio UI rendered its tool catalog, scientific mode selector and saved graphs. Normal API file/graph validation passed for the saved Mode 1 and Mode 3 graphs with the CPU controller explicitly placed locally; per-stage placement remained the saved policy.

The prepared Mode 2 retained-import reference returned HTTP 404 at the public API ownership/catalog guard. It was therefore not run or silently replaced. The intended continuation after a successful fresh Mode 1 run remains to reference that genuine managed native output; the lightweight retained source used by the earlier direct graph replay is not evidence of a fully verified public-API retained-source journey.

### Actual times and GPU observations

**Source inspection; no Real GPU execution.** These are systemd lifecycle timestamps in Europe/Berlin, not model runtimes:

| Host | First activation stop | Staged services started | Rollback began | Original services restarted |
| --- | --- | --- | --- | --- |
| pc002975 | 19:14:47 | 19:15:04 | 19:22:14 | 19:22:21 |
| pc004032 | 19:15:02 | 19:15:04 | 19:22:21 | 19:22:27 |
| pc004033 | 19:15:03 | 19:15:05 | 19:22:27 | 19:22:28 |

Controller service reconnect time was approximately 17 seconds during activation and seven seconds during rollback. The peer service stop-to-start intervals were approximately 2.45/1.82 seconds for activation and 5.46/0.50 seconds for rollback. These intervals use the individual host's journal timestamps; `Started` is a service lifecycle event, with API health verified separately. The staged services remained active for approximately seven minutes while validation and diagnosis occurred. None of these timings measures scientific throughput or model initialization.

`live/gpu-scheduler-samples.jsonl` contains 348 workstation samples across the three machines from 19:17:55 to 19:22:30 Berlin time (18:17:55–18:22:30 UTC). There are zero running Studio job observations and zero reservation observations. No new GPU UUID-linked scientific lease, native process receipt, native batch, prediction, initialization event, or overlapping scientific interval exists for this continuation. GPU confidence/acceptance/ranking results and scientific elapsed times remain unavailable, rather than assigned zeros or borrowed from retained work.

### Exact submission failure and recovery disposition

At 19:19:24 Berlin time, the saved Mode 1 graph was sent to the regular `/api/v1/pipelines/run` endpoint with unchanged scientific settings, `reuse_completed=false`, the CPU graph controller placed locally, and idempotency key `phase7b-live-native-a628c33-20261028`. No standalone native tool was invoked.

The native `engine_overrides_json` value is 417 characters. `portable_paths` classified the new `*_json` field as an operational file path because it was missing from the existing inline-content allowlist. Reading the queued worker request during dispatch consequently raised `OSError: [Errno 36] File name too long`. Shortening the JSON is not a valid solution: the field would still be a file reference and its content would not reliably survive serialization.

The rejected submission created local controller directory `campaign/20261028-191924-818b1d58` and wrote its input/worker request, but obtained **no scheduler job ID or GPU lease** and launched no worker/native process. This partially created controller was explicitly marked failed with `scientific_execution_started=false`; the original request files and an operational disposition record are preserved. The request-level idempotency entry remains uncertain and was not cleared or retried with a new key. This is an operational setup failure, not a native scientific rejection.

Modes 2/3 were not launched after the explicitly stated stop/rollback boundary. No accepted binder, predictor measurement, complete-native success, comparative rank or fresh ZIP was fabricated. Prior retained exports in section 16 remain separate and unchanged.

### Rollback verification and bounded local correction

**Source inspection.** Only the seven newly installed `phase7b.conf` drop-ins were removed, after copying their installed bytes to the private rollback backup. The original service definitions were left intact. Idle application bridges were restarted from the original checkouts. The original fingerprint `5569f8a69023be09c6079d2608fb383d60ba7a1dbcca1bd29e96c92c7fa53d5b` again appeared in all three fresh online peer snapshots. The controller API was healthy, authenticated requests to both peers succeeded, Studio queues/reservations were empty, and all seven original service-file hashes matched their pre-activation values. All 27 protected scientific source hashes also matched after rollback.

The minimal code correction adds `engine_overrides_json` to the same inline-content declaration already used for `options_json` and related engine input text. No scientific setting, algorithm, threshold, generator budget, scheduler policy or checkpoint changed. This correction is committed **locally only** as `191c5a1` on `codex/phase7b-acceptance`; it is not deployed, pushed or merged.

**Deterministic execution fixture.** The added regression first reproduced two failures on the staged code (long JSON worker serialization and the normal non-admin API path guard), with 11 other checks passing. After the one-line correction, 49 focused portable/native/API tests passed in 19.70 seconds. A further exact normal `/pipelines/run` journey using the existing PD-L1 routing fixture queued one controller with the full long native-settings JSON preserved; repeating the idempotency key reused that queue entry. This journey plus inline-input regressions passed 14 tests in 1.61 seconds. Its worker launch was deliberately stubbed and is not GPU evidence. These suites overlap. Source diff whitespace validation passed.

### Final acceptance status and remaining boundary

**Phase 7B is incomplete and the original functioning runtime is restored.** The approved source could be activated consistently, but its normal native graph submission contract failed. The safe stop/rollback policy was exercised successfully without interrupting scientific work. This is a discovered integration defect, not a reason to change the benchmark, thresholds, workflow semantics or generation counts.

A continuation requires review and authorization of the corrected runtime revision instead of `a628c33`, followed by another idle-only activation and the same bounded scientific acceptance. It must also replace the unmanaged Mode 2 retained-import reference with the genuine new managed native result, when available. No new architectural implementation phase is required for the local correction, but the approved runtime pin cannot be silently changed. Fresh complete BindCraft2/RFD4 execution, independent Boltz/Protenix evidence, real batching, four-GPU work, fresh Results/ZIP and scientific elapsed times remain unverified.

Full continuation evidence: `activation/*-before.json`, `activation/*-active.json`, `activation/*-rollback.json`, `activation/*-service-events.json`, `activation/rollback-connectivity.json`, `live/native_complete-submission.json`, `live/native_complete-job.json`, `live/*validation*.json`, `live/gpu-scheduler-samples.jsonl`, `live/inline-json-reproduction.log`, `live/inline-json-fix-regression.log`, `live/native-api-fix-regression.log` and `live/original-source-hashes-after-activation-rollback.json` under the evidence root. Summary: [phase7b-activation-rollback.json](evidence/phase7b-activation-rollback.json).

## 22. Authorized correction, redeployment and live continuation

**Source inspection + Deterministic execution fixture.** Following the user's explicit instruction to fix problems and deploy working corrections to the other machines, revision `191c5a183656962b8311f2b587d6f486f7acb620` passed 29 focused tests and was activated through the same reviewed service drop-ins with only a new isolated source/backup directory (`phase7b-v2`). All three active fingerprints matched `136917928e508d86f8b710055a936e31272626fb2829730cb82914e45e807ed4`; exact service and idle application-bridge process environments were verified. All APIs were healthy.

Two normal `/api/v1/pipelines/run` submissions then exposed a second integration defect before any native inference: `The candidate pipeline account is unavailable.` The local API tagged ordinary workflow controllers durably, but graph controllers lacked the development-account tag, and candidate batch lifecycle calls did not forward the existing `allow_development` permission. This did not indicate a scientific rejection. The native graph `20261028-193447-f20dc418` failed without a native child. The modular graph `20261028-193506-1465487b` preserved its two retained P5 candidate rows and exported partial Results despite generator/predictor dispatch failures. The normal Results API returned these two rows. Original services were restored after both controllers reached terminal states; no scientific work was interrupted. Raw evidence: `activation-v2/`, `live-v2/`.

Revision `9668069435b4e9b1485163b2e9c8c9c74b074e73` applies the existing ordinary-workflow submission-tag pattern to the graph API and passes the durable controller development flag through peer submit/poll/cancel. It does not change authentication mode or infer admin access from a missing account. The existing actor still requires either a real enabled user or explicit permission for the single synthetic development identity, and checks project access. Production/development API identity, idempotent submission, peer lifecycle, cancellation, campaign and API integration tests passed: **91 passed in 23.56 s**. No model inference is represented by these fixtures.

The revised deployment uses isolated `phase7b-v3/source` directories, pinned cleanly at `9668069`, on all three machines. Immediately before activation, app-filtered Studio queued/running jobs and scheduler reservations were empty and all four GPU utilization samples were zero. Only the seven Studio services and idle application bridges restarted. Original checkouts, units, unrelated scheduler/tunnel services, models and Results were preserved. After readiness settled, all APIs were healthy and exact process environments and fresh peer snapshots agreed on fingerprint `c9a5fd5e7e2300a3dff5943cba5e41447004ae79e43da3296f7ccc4bccc728bc`. An immediate first readiness request reached pc004033 before its API was listening; the subsequent check succeeded without another restart or code change. Evidence: `activation-v3/*-verify_active.json`.

**Live execution in progress.** Normal API submissions, each with a new durable idempotency key after verifying prior attempts performed no inference:

| Mode | Parent graph run | Submitted UTC | Saved graph revision |
| --- | --- | --- | --- |
| Complete native | `20261028-193839-caa59d85` | 18:38:39 | `b78b23a7` |
| Modular RFD4 + retained P5 | `20261028-193858-4f33be74` | 18:38:58 | `f474b3a1` |

Both use the normal CPU graph executor and distributed scheduler with the prepared scientific settings and retry limits. No native execution scripts were used. UUID-linked compute processes, app jobs, leases and GPU utilization are sampled every two seconds across the three machines into `live-v3/gpu-scheduler-samples.jsonl`. Final scientific outcomes, elapsed times, native counts, comparative Mode 2 and four-device acceptance remain pending until the actual receipts are collected.

## 23. Final deployment and bounded recovery

**Source inspection + Deterministic execution fixture.** Real prediction dispatch at `9668069` revealed a third integration defect: `resources.cpu_cores conflicts with parameters.cpu_cores`. The graph reserved its configured eight cores while the optimized-prediction adapter embedded a four-core default. Five candidate chunks failed validation before creating a native predictor process. The generator work was allowed to finish; its genuine outputs survived in partial Results. No active science was stopped.

Correction `2edf8de82c9d43962bcb53f8a48f79080f9eebc3` binds the node's CPU choice into the effective parameters before freezing the immutable request. Ordinary resource conflict validation remains intact. Real Boltz/Protenix request construction through the normal validator, CandidateSet lifecycle, native collection and graph suites passed **153 tests in 97.32 s**. The focused CPU regression passed 16 tests in 1.88 s; these suites overlap. No model, numerical setting or acceptance rule changed.

After all graph controllers and native jobs became terminal, original services were restored, then the reviewed override pattern activated new clean `phase7b-v4/source` directories on all three machines. Each activation rechecked app-filtered queues, scheduler leases and GPU activity. Service files and rollback backups remain available under each deployment directory. Final service/bridge process environments point to the pinned source and every fresh peer fingerprint matches. The APIs are healthy. The local controller UI is healthy and was used to inspect and download final Results. Existing runtime environments/images/checkpoints were reused; no installation, pull or image rebuild occurred.

The modular continuation retained the newly completed RFD4 run plus the original two P5 candidates. The only graph changes were the genuine retained source and explicit placement of useful existing prediction chunks: fresh Boltz on local GPU 1, retained Boltz on pc004032, fresh Protenix on pc004033, retained Protenix on local GPU 0. CPU/VRAM, scientific settings, one seed/sample and retry limits were preserved. It validated and ran through `/api/v1/pipelines/run` and the ordinary distributed executor. No standalone scientific runner was used.

The pc004033 Protenix process loaded the pinned checkpoint, then failed during its first model operation with `CUDA error: no kernel image is available for execution on the device`; its native CLI returned zero but the collection adapter correctly detected zero structures and marked the task failed. No prediction was fabricated from that exit status. Its failed manifest/logs and zero output count are retained. The CPU compiled-target probe was therefore insufficient to certify this installed stock runtime on RTX 5090. The exact underlying kernel/library cause remains unresolved; replacing images or numerical semantics was outside this bounded recovery.

The same stock Protenix model succeeded on both retained candidates on a 4090. Recovery then changed only the failed fresh-candidate Protenix placement to local GPU 1 and enabled normal completed-step reuse. All relevant cache manifests were verified before submission. The recovery ran **one** new predictor child, reused **13** steps, and completed with zero failed steps. No generator or successful predictor was resampled. This used seven candidate–predictor input attempts in total, below the original eight-input ceiling, across six distinct input pairs. The generation counts remained one complete BindCraft invocation/one trajectory and one RFD4 request/one replicate; optional redesign generated zero variants.

Raw evidence: `activation-v3/cpu-resource-*.log`, `activation-v4/`, `live-v4/*-submission.json`, `live-v4/*-job.json`, `live-v4/recovery-cache-verification.json`, `final-evidence/execution-receipts.json` and `final-evidence/active-runtime-verification.json`.

## 24. Three live workflow outcomes and elapsed times

**Real GPU execution.** All times below are UTC on 28 October; Europe/Berlin is UTC+1 on this date. Parent elapsed times include graph scheduling, staging, artifact collection and analysis; they are not native model throughput.

| Operation | Managed run | Start → completion UTC | Elapsed | Outcome |
| --- | --- | --- | --- | --- |
| Mode 1 graph | `20261028-193839-caa59d85` | 18:38:39 → 18:40:46 | 127 s | Completed, 0 eligible candidates |
| Native BindCraft2 child, pc004033 | `20261028-193907-c9b8a0a1` | 18:39:09 → 18:39:59 | 50 s worker interval | 1 trajectory, rejected at screen |
| Initial Mode 3 graph | `20261028-193858-4f33be74` | 18:38:58 → 18:42:43 | 225 s | Fresh generation completed; predictor CPU conflict preserved partial Results |
| Native RFD4 child, pc002975 GPU 0 | `20261028-193936-433a03d6` | 18:39:37 → 18:41:32 | 115 s worker interval | 1 sequence-complete 106-aa co-design |
| Mode 2 graph, pc004033 | `20261028-194227-fc3fdda9` | 18:42:27 → 18:42:28 | 1 s | Genuine empty native source, 0 predictor jobs |
| Mode 3 comparative continuation | `20261028-194549-e217df15` | 18:45:49 → 18:52:21 | 392 s | 3 Boltz + 2 Protenix outputs; 1 explicit Protenix failure |
| Mode 3 recovery | `20261028-220218-57400312` | 21:02:18 → 21:05:36 | 198 s | 3 candidates, 6 predictions, 13 cached steps, 0 failures |

**Mode 1.** The normal complete `design.bindcraft2` request recorded `trajectory_only=false`, the approved one-trajectory maximum, one worker, length 60, seed 20261027 and the preserved one-step stage overrides. The native campaign state reports one screen termination, zero MPNN-eligible trajectories and zero accepted candidates. The native log reports screen iPTM **0.12**, pLDDT **0.31**, rejection reason **pLDDT**. Native filters were unchanged. The failed trajectory CIF, losses CSV, campaign state, settings, checkpoint metadata and native stdout are preserved on pc004033 and in `live-v3/native-rejection-evidence.tar.gz` / `final-evidence/pc004033/20261028-193907-c9b8a0a1/`. Internal ProteinMPNN, AF2 validation and later filtering were not reached; running the complete workflow entrypoint is not proof that these conditional stages executed. No extra trajectory was generated to force passage.

**Mode 2.** The graph uses the genuine completed managed native child on its execution workstation. Only the retained reference and its identical target path were changed from the prepared graph; target SHA-256 still equals `7da1c888011cdbd0937e3ff6b73429aa069e4c496cfa0cb8d8b53cc06ca73b57`. This avoids the earlier unmanaged-import API guard failure. Both enabled predictor branches recorded `batches=[]`, `processing_jobs=0`, and skipped empty input; analyzer, ranking and Results completed with zero rows while retaining the native rejection provenance. No RFD4 candidate was relabeled as a complete native output. Positive-candidate Mode 2 scientific evaluation remains **Unverified**, as required for zero eligible native candidates.

**Mode 3.** One RFD4 guided co-design used the established PD-L1 chain A, hotspots A40/A99/A107, seed 20261028, 60 steps and the retained bounded settings. Its output was sequence-complete and bypassed initial ProteinMPNN. Optional redesign was disabled. The two existing P5 candidates were reused without new generation. Both independently configured predictor branches, unified analysis, cohort ranking and export executed automatically. Recovery required a new normal graph submission after the empirical runtime failure; this is not a claim that the initial end-to-end run was flawless.

## 25. Physical GPU, scheduler lease and native process evidence

**Real GPU execution.** Scheduler job IDs below identify the application leases; full lease owner IDs, run paths, GPU UUIDs, first/last samples and receipts are in `final-evidence/gpu-execution.json`. Native process intervals are conservative first/last `nvidia-smi --query-compute-apps` observations, sampled approximately every two seconds. They establish device-linked process activity, not continuous forward-pass utilization.

| Workstation / GPU UUID | Useful workload | Scheduler job / lease | Native PID | CUDA process observations UTC |
| --- | --- | --- | --- | --- |
| pc004033 / `GPU-370f52fb-5676-7316-104d-5a9bf016fae1` | BindCraft2 native rejection | `2495488a005145938343cf57aec0c70d` | 142401 | 18:39:11.615–18:39:57.302 |
| pc002975 GPU 0 / `GPU-195a1966-1148-340e-5fe2-bac2d7a9821a` | RFD4 fresh co-design | `4005bb5ab916469fae7f391abe8441fe` | 995823 | 18:40:30.966–18:41:26.344 |
| pc004033 / same UUID above | Protenix failed first model operation | `f8c7e7fc2e074391a86886f630122492` | 161074 | 18:47:00.970–18:47:11.313 |
| pc002975 GPU 1 / `GPU-c0569fbb-42d6-594e-a304-0a84a0c60a21` | Boltz, fresh 106-aa input | `fcf13b7d0b13454998be753fcb409dcf` | 1015532 | 18:47:23.299–18:47:46.145 |
| pc004032 / `GPU-b7f8339b-157f-e1f2-f136-bd0ad68bc780` | Boltz, two retained inputs | `ac10d8766bcb42b5925d568ddc839832` | 106184 | 18:47:46.162–18:47:51.087 |
| pc002975 GPU 0 / same UUID above | Protenix, retained 94-aa input | `7b5db3f80b514d05bc1f550f4ffd7b67` | 1016592 | 18:49:06.051–18:49:35.140 |
| pc002975 GPU 0 / same UUID above | Protenix, retained 116-aa input | `092f39c7c96c49d6aa3bacc0008c2d41` | 1025788 | 18:51:40.920–18:52:09.927 |
| pc002975 GPU 1 / same UUID above | Protenix recovery, fresh 106-aa input | `0fb693cfed964fac96476c164b960674` | 1261764 | 21:04:56.536–21:05:25.574 |

All four UUIDs have successful useful scientific execution: native BindCraft rejection on pc004033, Boltz on pc004032 and local GPU 1, and RFD4/Protenix on local GPU 0. The failed Protenix job is additional evidence, not the only scientific activity on pc004033.

**Concurrency limit.** Maximum observed concurrent GPU leases was **three**, including preparation/loading: 18:46:39.000–18:47:13.735 and 18:47:28.178–18:47:51.161 UTC. The two native-design leases overlapped 18:39:37.969–18:40:02.117 UTC. The conservative CUDA-process observation intervals do not overlap sufficiently to establish simultaneous multi-device model forwards; four-process simultaneous inference is **Unverified**. In particular, the two Boltz observation intervals nearly abut and must not be promoted into proven overlap. Reference fingerprint publication delayed pc004032's eligibility until after pc004033's Protenix failure, and the local Protenix jobs spent substantial time initializing before CUDA work. No redundant design/prediction was added to manufacture occupancy. An unrelated RustDesk process was excluded from scientific counts.

The monitor stopped after final idle checks. No queued/running Studio jobs or GPU reservations remained on any workstation. Shared compute scheduler and SSH tunnel services were neither modified nor restarted.

## 26. Actual native batching, measurements and scientific integrity

**Real GPU execution.** Counts below come from returned `native_execution.json`, native batch manifests and logs, not execution-plan forecasts. There were **eight native scientific processes/allocation events**: one BindCraft2, one RFD4, two Boltz and four Protenix (including one failed attempt). Five initial predictor chunks plus one recovery chunk processed seven input attempts and returned six successful predictions. The five CPU-conflict submissions and earlier account/inline failures launched zero native inference and are excluded from these counts.

| Native operation | Candidate inputs / chunks | Native commands / Docker processes | Successful outputs | Native manifest elapsed |
| --- | --- | --- | --- | --- |
| BindCraft2 | 1 whole native invocation, 1 trajectory | 1 / 1 | 0 eligible; 1 genuine rejection | 50 s worker interval; native-only duration not separately measured |
| RFD4 | 1 design, microbatch 1 | 1 / 1 | 1 co-design | 115 s worker interval; generation log retained |
| Boltz on local GPU 1 | 1 input / 1 chunk | 1 / 1 | 1 | 81.455 s |
| Boltz on pc004032 | 2 inputs / 1 chunk | 1 / 1 | 2 | 24.046 s |
| Protenix failed on pc004033 | 1 input / 1 chunk | 1 / 1 | 0 | 44.337 s |
| Protenix local GPU 0 | 2 inputs / 2 chunks | 2 / 2 | 2 | 178.200 s + 151.584 s |
| Protenix recovery GPU 1 | 1 input / 1 chunk | 1 / 1 | 1 | 158.586 s |

Boltz's two-input native manifest records `native_yaml_directory_single_persistent_worker`, one native command, seed `[42]`, two distinct input IDs and two returned structures. Its stdout reports processing/predicting two inputs. This verifies a real multi-input native process, not a tensor forward batch size of two. The two Protenix retained inputs used two one-input native directory invocations and separate model/checkpoint loads; live multi-input Protenix model reuse remains **Unverified**. Protenix logs show one checkpoint load per process; no total model-initialization count is inferred for BindCraft's conditional stages or hidden model internals.

**Scientific identity verification.** Read-only parsing checked all original and successful predicted complexes against exact binder sequences and the PD-L1 115-aa chain-A sequence hash `c434254957dd7a7d80ce4ca859f9a6d101a1205f4ff822177dd0b3c31f02e968`. Every prediction joined the correct variant, with a distinct prediction identity, seed 42, image identity, recorded chain roles and source ancestry. Native measurements/filters remain separate from independent predictor measurements. No MSA/template, acceptance threshold or native filter was changed.

| Source / binder length | Binder sequence SHA-256 | Boltz iPTM | Protenix iPTM | Recorded cohort rank |
| --- | --- | --- | --- | --- |
| Fresh RFD4 / 106 | `ffff82c4ee4f91b3dd2d4ad345c02a4700b901882fe4cf61aacb2c1530aca310` | 0.59279 | 0.28279 | 1 in its own recorded cohort |
| Retained P5 / 116 | `ad6a651791941cbe4fbeebf027fb6c46b470034ebc4c3a2c66d80598285cdf68` | 0.48661 | 0.31441 | 1 in P5 cohort |
| Retained P5 / 94 | `da0ba51511da10c21f61dfc4e86e686907a2e733d2e62ae19a39f1ea32e6acee` | 0.41334 | 0.58751 | 2 in P5 cohort |

These are fresh computational confidence measurements, not experimental binding results. They differ from the original retained Boltz measurements and do not overwrite them. Ranking correctly follows the configured Boltz iPTM field descending within recorded comparable cohorts; it neither averages engines nor infers one global rank across the fresh/P5 cohorts. All three rows are selected by the configured unthresholded selection policy; this is not a native acceptance or experimental-binding claim. Fresh RFD4 native filter status remains `not_evaluated`; the retained P5 native filter status remains `unresolved`.

The original P4/P5/native source ledgers and 27 tracked scientific files retain all pre-run hashes. The three original Git checkouts remain at `abf0f2e`, with their exact modified/untracked paper-reference files preserved byte-for-byte.

**Resource receipt limit.** The BindCraft graph requested an 8-GiB minimum, but its native worker scheduler receipt reports `min_vram_gb=0.0`. CPU count, one GPU and the scientific trajectory cap were enforced, and the selected idle 5090 had ample memory, but end-to-end propagation of that admission floor is not verified. The actual receipt is retained rather than rewritten to the graph's requested value.

**Provenance limits.** BindCraft's original worker `native_execution.json` confirms one returned Docker/native invocation, but its collected graph `source_provenance.native_execution` is null. The root receipt was retrieved read-only separately and preserved with the native source snapshot; the report's actual count uses that receipt. The original native campaign records AF2 checkpoint hashes and its image-internal MPNN checkpoint hash (`b166276f80a45cb3de24bceb3b47b2f47a60ff654660f09837cd775846411ad1`), without claiming that MPNN/validation ran. Boltz candidate records contain full checkpoint hashes. Protenix candidate `checkpoint_provenance` and `effective_settings` fields are empty: the actual model/settings are instead traced through the read-only reference mount, native checkpoint-load logs, unchanged preflight/final file hashes and exact CLI/batch manifest. This export-field gap remains unresolved and is not described as complete per-row provenance.

## 27. Results, exports and partial-output preservation

**Real GPU execution + Retained scientific artifact.** Normal Results APIs returned the native, comparative-empty, partial modular and recovered modular outcomes. Studio's final Results page displayed three candidates, six completed predictions, the saved graph and downloadable data. The native Results page displayed zero candidates and preserved rejection/native-source counts. Both native and recovered modular ZIPs were downloaded through the UI. The recovered export contains 95 files, including 13 structure files, all three exact selected FASTA sequences, six predictor records, native-source provenance, execution graph and ranking decisions. ZIP CRC and all **25** export-manifest SHA-256 checks passed.

| Export | SHA-256 | Verification |
| --- | --- | --- |
| Native graph UI download | `19a280c91b5e7a67e5d13a3503caa2a36ddc9ca4c14d543b39962da7bdb3165e` | Valid ZIP, genuine screen rejection provenance |
| Mode 2 empty graph ZIP | `472c5ee9799b66a90d02fb8afbdf1d6e4b40f2c80820c8ff45b83ea40a13a97b` | Standard `build_run_zip` on preserved managed-run snapshot; empty native provenance |
| Partial modular graph ZIP | `243c60265b052e3eaf24d87e0359ee8d33c5a2448696da1c397e47c683ea44ba` | Standard ZIP helper; 3 candidates and 5 successful predictions survive failure |
| Recovered modular UI download | `1b2388e95c047613e4d7db97d9b91df94501bcee96ba6d7e3304739788d41ebe` | 3 candidates, 6 predictions, exact FASTA, all checksums valid |

Local UI downloads remain in `/home/user/Downloads/`; evidence copies are `live-v3/native-complete-ui-download.zip` and `final-evidence/modular-recovered-ui-download.zip`. The Mode 2/partial ZIPs are under `final-evidence/`. Raw native rejection artifacts are additionally preserved in their separate source-evidence archive; an empty-candidate graph ZIP is not claimed to contain every native compile/trajectory file.

The partial modular run remains immutable, including its failed 5090 prediction. Before recovery, independent identity/ranking checks passed for its three candidates and five successful predictions. After recovery, the same variants and five original successful prediction records survived; the sixth result was added through ordinary cache reuse and analysis. No completed native job was repeated. Deterministic recovery tests remain separately labeled; no live cancellation, node-disconnect or interrupted-transfer scenario was intentionally injected during active science.

## 28. Final acceptance disposition

**Bounded live execution is finished; Phase 7B acceptance is partial.** The deployed fixes work across all three machines and the bounded scientific runs, recovery, Results and exports have real execution evidence. There are no unfinished acceptance jobs. Corrected runtime `2edf8de` remains active and healthy; the original runtime/units and every deployment backup remain available for rollback. Code corrections and this report are local commits only; nothing was pushed or merged.

| Acceptance item | Final finding |
| --- | --- |
| Consistent active source/services and authenticated workers | Verified on all three machines |
| Four physical GPUs dispatch real useful science | Verified across the run with UUID/lease/process evidence |
| Four simultaneous native inference processes | Unverified; three overlapping leases include preparation, and conservative CUDA-process intervals do not demonstrate concurrent forwards |
| Bounded complete-native entrypoint and true rejection | Verified; downstream internal MPNN/AF2 stages not reached |
| Positive Mode 2 evaluation | Unverified because zero eligible complete-native candidates; empty routing verified |
| Fresh RFD4 + automatic downstream evaluation | Verified with documented integration/runtime failure and normal graph recovery |
| Exact identities, separate predictor measurements and cohort ranking | Verified for all three variants and six successful predictions |
| Native multi-input batching | Boltz two-input invocation verified; Protenix multi-input reuse unverified |
| Completed outputs survive partial failure and recovery | Verified live; 13 cached steps, one recovery prediction, no duplicate generation |
| Full per-row provenance export | Partial: Protenix checkpoint/settings fields and collected BindCraft native-receipt gap documented above |
| Results / downloadable ZIP | Verified; exact sequences and checksums passed |
| Original work, scheduler/tunnels, images and scientific limits preserved | Verified within recorded evidence; no GitHub push/merge |

Remaining priorities are narrowly identified: diagnose the installed Protenix stock RTX 5090 kernel mismatch and make eligibility reflect actual runtime support; complete missing per-row/native-receipt export coverage; demonstrate concurrent useful native forwards and multi-input Protenix reuse only when independent eligible work justifies it; evaluate a genuine positive complete-native output when one exists. These are not authorization to expand generation, replace images or start another architectural phase.

Machine-readable final summary: [phase7b-live-acceptance.json](evidence/phase7b-live-acceptance.json). Full original worker receipts/logs, API outcomes, GPU lease/process samples, identity audit and verified exports are under `final-evidence/` in the evidence root. The final audit helper only reads existing outputs and starts no scientific job.


## 29. RTX 5090 root cause and canonical image correction — 29 October

**Source inspection + Real GPU execution.** The user explicitly authorized fixing the Docker image and deploying working fixes to the peer machines. The canonical local repository is `mn-tool-containers` (remote `martenbeeg-sketch/mn-tool-containers`), matching the user's container-repository clarification. Dockerfiles remain there. Original dirty checkouts were preserved; changes live on isolated `codex/protenix-blackwell`.

The old image's PyTorch advertises SM120, but stock Protenix loads its own `fast_layer_norm_cuda_v2.so`. `cuobjdump` found native SM70/80/86/89/90/100 kernels and **no SM120 or PTX**. The pinned extension compiler hardcodes a target list ending at SM100. On pc004033, an independent 48-case PyTorch LayerNorm diagnostic passed; the combined stock fused-LayerNorm/reference diagnostic failed with CUDA's missing-kernel-image error. That error surfaced at the next PyTorch operation, despite synchronization, so its traceback alone did not identify the originating launch. Binary inventory and the successful rebuilt-extension test establish the missing stock-extension target as the compatibility defect. No GPU hardware failure was demonstrated.

`protein-design/protenix-opt/Dockerfile.sm120` appends SM120 to that compiler list and rebuilds only the stock extension. The original CUDA/PyTorch stack, Protenix source, model weights, optimizer kit, mn-predict adapter, entrypoint, environment and every base image layer are retained. Build took 181.6 seconds. Source recipe commit: `93e8aca5899ebcdd22ceefb6fe4bd28fa5c8b6a5`; the build/catalog helper is `./build.sh protenix-v2-opt-sm120`. The ordinary batch build also builds this follow-on image (integration/evidence commit `5d01e3c8ba0e0c15a59b711e614f14e76a39ad3f`). A container-source Git bundle is preserved locally and on both peers. No old image tag was replaced.

New image: `mn-protenix-opt:2.0.0-cu128-sm89-sm120-batch3`. Its portable content fingerprint matches on all three hosts: `sha256:d9d00bf5909b607c65b7d9aaeb57c4c7a50e5a40bd0d110da2ea5ccc266fecb4`. Local Docker ID is `7468360d…`; peer IDs are `bb9be7ff…` because their stores expose identity differently. The rebuilt extension SHA-256 is `5aad24cec735c5f7d7f4c1ab7eb1616ffb1ef5d30e3cc990737615864d620885`. One image includes both SM89 and SM120; separate Ada/Blackwell images are unnecessary for this compatibility fix. Separate optimized variants remain a future evidence-based choice.

All 48 diagnostic cases passed independently on pc002975 GPU0, pc004032 GPU0 and pc004033 GPU0: float32/float16/bfloat16, four widths and four scale/bias combinations. Largest absolute difference from the PyTorch reference was 0.0078125, in bfloat16 and within the stated diagnostic tolerance. Diagnostics used scheduler leases and are **not molecular inference**. An initial diagnostic command omitted interactive stdin and executed no probe; its empty successful shell exit is explicitly excluded from validation evidence.

Studio revision `3260986` selects the new image and inspects required Protenix extension targets as well as PyTorch targets. The old image now fails eligibility on SM120 while retaining SM89 eligibility; missing extension inventory fails closed. Architecture-specific PTX is not misclassified as portable generic PTX. Regression validation: 133 cluster/graph/native-configuration tests passed, plus 22 focused tests after image selection. Shell syntax and catalog checks passed.

## 30. Controlled activation and real PD-L1 execution with the fixed image

**Real GPU execution.** Before activation, all three hosts had no Studio queued/running work, no scheduler reservations and no running containers. GPU0 on each host was idle. Local GPU1 had external/display activity and was left alone. Seven Studio service overrides were replaced using exact backups of the active v4 overrides, with rollback restoring v4. Only idle Studio application bridges were refreshed. Shared schedulers and SSH tunnels were untouched. APIs, agents, controller UI and bridges now run clean v5 source at `3260986`, fingerprint `24de4cb1fdff303ec0fbd092f78af1d1222e5ed19b6c82a0db322021d3ceddb7`; all peer heartbeats agree and APIs report healthy.

Normal Studio graph `20261029-090448-18c4a6a4` reused the fresh retained RFD4 candidate; **zero new designs** were generated. Two independent Protenix branches ran through the normal graph executor and distributed scheduler, locally on a 4090 and remotely on pc004033's previously failing 5090. Both retained seed 42, 10 recycles, 200 diffusion steps, one sample, bfloat16, no MSA and no template. Target A is 115 residues and binder B is 106 residues, with exact hashes from section 8 and the fresh candidate record preserved. Native input fields equal the old 4090 input except for the descriptive name.

| Workstation / UUID suffix | Scheduler lease job / native host PID | GPU process observations (UTC) | Native forward | iPTM / pTM |
| --- | --- | --- | --- | --- |
| pc002975 / `bac2d7a9821a` | `b3684c821b0e43b78e1f25db372890f2` / `2231745` | 08:07:22.906–08:07:51.893 | 13.71 s | 0.25146 / 0.58199 |
| pc004033 / `5a9bf016fae1` | `1023e79007d94269b9497771f8f09b05` / `383789` | 08:05:31.445–08:05:50.777 | 7.39 s | 0.25138 / 0.58415 |

Each native receipt records one input chunk, one Docker process, one native command and one prediction output. Each log contains exactly one checkpoint load; model initialization counts are verified from logs rather than the receipt's forecast wording. Total: **two native invocations, two initializations, two successful predictions**. Parent elapsed time was 190 seconds, including unequal CPU initialization, staging and collection. The model times are a single-case observation, not an optimized-image benchmark. The native forward intervals did not overlap; no new concurrent-forward claim is made. pc004032 passed its real kernel test but did not receive redundant protein inference. The earlier four-physical-GPU scientific evidence remains unchanged.

The test graph initially omitted a standard shared-filter decision, so ranking correctly remained unresolved while both predictions survived. A second normal graph, `20261029-091307-c02a3682`, used the unchanged prepared filter policy and **cached both predictors**, then recomputed analysis/ranking/export in five seconds. The filtered source view feeds analysis; predictor inputs and scientific cache keys stay unchanged. The single variant ranks first by the explicitly declared 5090 comparative iPTM. This is relative selection, with `native_filter=not_evaluated`; it is not native acceptance. No prediction, design, threshold or budget was repeated/relaxed to obtain that selection.

## 31. Scientific integrity, remaining limits and current final status

**Retained scientific artifact + Real GPU execution.** Both collected predictions preserve the exact variant, source candidate, 106-aa binder SHA-256 `ffff82c4ee4f91b3dd2d4ad345c02a4700b901882fe4cf61aacb2c1530aca310`, and 115-aa PD-L1 chain mapping. Original native measurements remain separate from the two predictor namespaces. Results API and authorized CSV download returned HTTP 200. Studio's own Results ZIP builder produced a 671,304-byte archive with 49 entries, both prediction records, their structures and 17 verified export checksums. ZIP SHA-256: `a2164076f1c1e25ec54f4edd353219adc10a62ccde0cc33027aaca4825964a83`.

**Numerical equivalence remains Unverified.** Close iPTM values do not establish matching poses. Sequence-matched all-residue Cα fitting gives 6.368 Å complex RMSD between the new 4090 and 5090 outputs; target-fit binder pose RMSD is 14.158 Å. The old versus new 4090 outputs also differ (5.577 Å complex; 9.991 Å target-fit binder), with old iPTM 0.28279. This protocol did not request deterministic algorithms; exact inputs/settings and identity were checked, but the source of the pose variation has not been isolated. No structural-equivalence or scientific-quality certification, improved binding claim, or positive native candidate is inferred. No equivalence threshold was set or relaxed after seeing the results.

Checkpoint SHA-256 `8f931f9774a396b67033d0e58628e1834f4a1448165e04254b40a780b0c0d599` still matches on all three machines; weights remain external and read-only during inference. Each prediction retains its actual image ID, seed, identity check and native command. The existing per-row `checkpoint_provenance` field remains empty; this report's independent checkpoint hash does not silently fill that export gap. All 27 protected historical scientific file hashes and all old image IDs remain unchanged. No GitHub push or merge was performed.

**Current acceptance:** the stock Protenix RTX 5090 execution defect is corrected and deployed on all three machines, with real 4090/5090 molecular inference and end-to-end Results verification. Overall Phase 7B remains **partial**: positive complete-native Mode 2, downstream BindCraft stages after the genuine rejection, four simultaneous useful native forwards, per-row checkpoint/native-receipt export coverage, numerical equivalence and Blackwell `exact`/`fast`/`big` optimization-mode validation remain open. No new architectural phase or additional scientific generation was started.

Portable image/activation/execution summary: [protenix-blackwell-compatibility.json](evidence/protenix-blackwell-compatibility.json). Raw receipts, original failures, kernel checks, native logs, UUID samples, structures, audits, CSV and ZIP are under `/mnt/data/RESULTS/mn-molecular-design-studio/verification/protenix-blackwell-image-20261029`. v5 service backups and rollback scripts are under the existing Phase 7B evidence root's `activation-v5/` and each host's `phase7b-v5/rollback-units` deployment directory.

## 33. v14 deployment and objective-continuation verification — 29 October 2026

**Deployment.** Release `phase7b-v14` uses the actual v12 deployment line (`e1ebb2c`) with the already-active v13 controller fixes, then adds the CandidateSet/Results projection fix in local commit `24f3f43`. It was installed with the v13 idle-check/rollback procedure, not the stale v9 scripts. Before activation each host had no active Studio generation jobs or scheduler reservations. The local RTX 4090 display workload on GPU1 was left undisturbed. API, agents, UI and idle application bridges are healthy on pc002975, PC004032 and PC004033; all report source fingerprint `f791da2d62145d14762d16bbbf269de774c0f928daff2e1159a797422660067b`. Original checkouts, service rollback copies, scheduler and SSH tunnels were preserved. Deployment manifest and pre/post checks are in `activation-v14/`.

**Real normal-API/UI campaign.** The bounded PD-L1 run `20261029-103552-de38a6af` completed through the Studio graph executor and distributed campaign controller in 9m13s as shown in the Results UI (9m16s from API create/complete timestamps). It used the existing two-attempt maximum and the explicit Boltz-2 iPTM criterion `>= 0.7`; no budget or threshold was expanded. Attempt 1 (seed `2075037775`, native candidate `dc_00000001_e433786c54328396`, CandidateSet `cs_91004477c665cfc63d123ce58e9bbe722db2eb7bc22615bd1b61230fc8881b69`) was genuinely rejected at iPTM `0.67484`. The controller then generated attempt 2 automatically (seed `2031933765`, native candidate `dc_00000002_e433786c54328396`, CandidateSet `cs_286d76d3e052339cd2e6133e55bd7bc22cb2e3d06c0e2ad4c00a278bb5ae9cd2`); its Boltz-2 iPTM was `0.71954`, meeting the objective. The campaign terminated `target_reached` at 1/1 on attempt 2 of 2. Both candidates, actual per-attempt reasons, and prediction provenance are visible in the v14 Results table and API payload. Native acceptance remains `unavailable` for RFD4 and is displayed separately from validated acceptance and ranking. The v14 API/UI re-read the already-completed real campaign; no duplicate designs were generated solely for verification.

The Results page's downloadable job ZIP was verified readable (33 entries, no corrupt member; includes `campaign_result.json` and `candidates.jsonl`). API/UI snapshots and the revalidated payload are `activation-v14/api-job-real-campaign-v14.json` and `api-results-real-campaign-v14.json`; the downloaded ZIP is `PD-L1_RFD4_validated_acceptance_continuation_v13_20261029-103552-de38a6af.zip`.

**Shared multi-generator objective: not implemented.** The existing controller remains one generator per ledger. `validate_graph` in `pipeline_graph.py` attaches a campaign only to a BindCraft2 or RFD4 tile and overwrites `generator_id` from that tile. `graph_pipeline.node_request` submits that single workflow, while `_campaign_ledger_path`, `_create_attempt_child`, `_reconcile_attempts` and `_run_distributed_campaign` in `engine/workflow_engine/workflows/distributed_campaign.py` bind the ledger and attempt adapters to the one parent input and its generator-specific staged files. Independent generator tiles therefore produce separate parent jobs/ledgers; CandidateSet merge combines their outputs after execution but does not give them one shared target, attempt budget, acceptance counter or recovery state. The exact implementation boundary is a graph-level campaign coordinator that stages per-generator inputs under one durable campaign parent, assigns each authorized attempt to its adapter, and records every branch receipt in the existing ledger/CandidateSet path. The UI/config entry point is `distributed_campaign_config_ui.py` / `pipeline_node_ui.py`; graph validation, execution, controller and recovery changes belong in `pipeline_graph.py`, `engine/workflow_engine/workflows/graph_pipeline.py`, and `engine/workflow_engine/workflows/distributed_campaign.py`. A schema-only flag would falsely advertise functionality, so none was added.

**Status.** The rejection-triggered continuation loop is operational and verified via the normal v14 API/UI with genuine GPU generation and Analysis evidence. The shared multi-generator objective remains the next code change; no second scheduler or broad architectural rewrite is justified. Changes remain local; nothing was pushed or merged.


## 32. Production campaign follow-up — 29 October 2026

**Real GPU execution.** Complete native graph `campaign/20261029-092228-2f0240a5` used the PD-L1 benchmark, native settings and filters, a cap of four trajectories and one desired native result. The first trajectory failed the native annealing iPTM filter. The second reached gradient design, screening, refinement, annealing, hardening, mutation, ProteinMPNN and AF2 validation. Three of ten redesigns passed the native filters; the native workflow selected its best candidate. Native accepted iPTM was 0.80 and i_pDAE was 0.23. These are native measurements, not independent predictions.

The genuine accepted candidate was reused in graph `campaign/20261029-092753-33a4b24b` without repeating native generation. Its 60-residue binder sequence SHA-256 was `6b16f77d8743eb6fc864b019a05f6e4ae53662220b898be78eb62ee58c72dd8d`. The target was chain A, 115 residues, sequence SHA-256 `c434254957dd7a7d80ce4ca859f9a6d101a1205f4ff822177dd0b3c31f02e968`. Boltz-2 reported iPTM 0.6834 and pTM 0.86425; Protenix V2 reported iPTM 0.58092 and pTM 0.6102. The shared result ranked the candidate first using Boltz iPTM. The audit checked identity and source separation, and the Results ZIP checked 16 file checksums. Per-row Protenix checkpoint provenance remains empty; it must not be inferred from the prediction score.

The new campaign controller supports generated-count or accepted-count objectives, optional unique accepted sequences, bounded attempts and runtime/compute settings, and per-validation-stage AND/OR conditions with explicit missing-data behavior. Existing ledgers retain their accepted-objective default. BindCraft2 peer dispatch retains the complete native workflow. A race in graph ownership attribution was fixed before the retry. Scientific thresholds were not relaxed.

**Real GPU execution.** Generated-objective graph `campaign/20261029-093720-e6f21fe9` ran two separate RFD4 native generator jobs and reached two generated candidates. PC004033 and PC004032 ran one native process each with an observed overlap from 08:38:41 to 08:39:15 UTC. The initial downstream redesign failed because historical checkpoint paths were incorrectly treated as required worker files; both generated results survived. Revision `2a67b92` preserves historical measurements/provenance and queues work when distinct peers are busy. The two saved designs were selected for downstream recovery in graph `campaign/20261029-094343-58e50522`, without repeating generation. Earlier retry `campaign/20261029-093446-df50af61` failed before GPU execution due to the ownership race and is retained as failure evidence.

The sampled scheduler and GPU evidence is `/mnt/data/RESULTS/mn-molecular-design-studio/verification/production-campaigns-20261029/gpu-scheduler-samples.jsonl`. Four physical GPUs ran genuine native processes: PC004033 `GPU-370f52fb-5676-7316-104d-5a9bf016fae1`, PC004032 `GPU-b7f8339b-157f-e1f2-f136-bd0ad68bc780`, and pc002975 `GPU-195a1966-1148-340e-5fe2-bac2d7a9821a` and `GPU-c0569fbb-42d6-594e-a304-0a84a0c60a21`. Sampling records lease run IDs, timestamps, process IDs and GPU UUIDs. Four simultaneous processes are not claimed.

**Current limits.** The controller presently schedules one native attempt per execution job, even when a generator can run multiple designs natively. A deterministic 1,000-output test verifies controller bookkeeping, not native multi-design batching or a 1,000-protein experiment. Combined objectives spanning generators, arbitrary nested cross-stage boolean criteria, and streaming graph fan-out while generation continues are not implemented; the controller does overlap its own validation with generation. Generator-specific objective support currently covers BindCraft2 and RFD4, not every generator tile. These gaps prevent declaring the requested production definition of done satisfied.

Evidence, submissions, tests and the comparative ZIP are under `/mnt/data/RESULTS/mn-molecular-design-studio/verification/production-campaigns-20261029`. All code commits are local; nothing was pushed or merged. The original checkouts and scientific data were retained, and service activation used idle-only checks with rollback copies.

The following intervals are sample bounds (UTC), rather than exact process start and stop times. `run_id` identifies the scheduler reservation; process IDs are from the workstation's `nvidia-smi` sample.

| Workstation / GPU UUID | Scheduler reservation run ID | Native process | Observed interval |
| --- | --- | --- | --- |
| PC004033 / `GPU-370f52fb-5676-7316-104d-5a9bf016fae1` | `c9b35b1c4b444a529264fa580af5fb0b` | RFD4 Python PID 425262 | 08:38:32–08:39:15 |
| PC004032 / `GPU-b7f8339b-157f-e1f2-f136-bd0ad68bc780` | `adb1441622454dfb96bc29dcf804d227` | RFD4 Python PID 353848 | 08:38:41–08:39:17 |
| pc002975 / `GPU-195a1966-1148-340e-5fe2-bac2d7a9821a` | `dc495dce0a0545aba59605c3993b03a6` | LigandMPNN Python PID 2329110 | 08:44:15–08:44:22 |
| pc002975 / `GPU-c0569fbb-42d6-594e-a304-0a84a0c60a21` | `182d345ca1ef445d8545c8ce82c587b9` | LigandMPNN Python PID 2329747 | 08:44:22–08:44:29 |

The two RFD4 jobs used two independent native processes, each requesting one design with native batch size one. The two redesign jobs each processed one candidate. These are execution counts, not an extrapolation from collection size or a claim that generation used native multi-design batching.

**Complete native peer campaign.** Graph `campaign/20261029-094447-d4f389fb` used the new distributed BindCraft2 adapter, accepted target one, maximum two native trajectories, and unchanged native filters. Attempt one ran on PC004033 as cluster job `8a82d216-bd5a-4768-b3dd-830120a04c49`. Its real native log reached all design stages, ProteinMPNN and AF2 validation; three of ten redesigns passed and native selection chose candidate 8. The controller collected its authoritative native campaign state and recorded candidate `dc_00000001_9df225e6a61e5793` as native accepted. The campaign stopped at `target_reached` after one attempt. Native acceptance is distinct from comparative acceptance; no independent validation was requested for this particular candidate.


**Recovered modular result.** The original recovery graph completed both predictor branches (four Boltz-2 and four Protenix predictions) but its analyzer failed because the saved source target lacked a computed identity. The outputs remained intact. Revision `e1ebb2c` repairs the target handoff. Recovery graph `campaign/20261029-095604-8df0a975` reused all generator, redesign and predictor steps from the verified Studio cache and completed analysis, ranking and Results without rerunning inference. An audit checked exact binder and target sequences for all eight predictions and confirmed two original designs and two separately identifiable redesigns. Ranked by Boltz-2 iPTM, the four scores were 0.85530, 0.83795, 0.80955 and 0.64706; the corresponding Protenix iPTM values were 0.41445, 0.40633, 0.49136 and 0.41581. Ranking, rather than scientific binding, is what these scores establish. The ZIP was readable and its 29 export checksums matched. Full sequence hashes and predictor identities are recorded in `modular-audit.json`; the export is `modular-results.zip` in the evidence directory.

**Real validated objective.** Graph `campaign/20261029-094642-a66be082` requested one unique RFD4 candidate passing Boltz-2 iPTM >= 0.7, with at most two attempts. Its first native generation had native acceptance unavailable; independent validation returned iPTM 0.89697, so the validated objective reached its target after one attempt. This verifies genuine downstream acceptance and stopping, but not live automatic continuation after rejection. The validation job `refolding-validation/20261029-095105-6aa0a12e` incorrectly carried the internal labels `benchmark_engine=bindcraft2` and `parent_role=distributed_bindcraft2_validation`; its durable candidate lineage and prediction identify RFD4 and Boltz-2. Revision `631bf55` fixes those labels for subsequent jobs without rewriting the original receipt.

**Elapsed wall time.** From Studio job creation to completion: original complete native 296 s; native comparative graph 126 s; two-job RFD4 generation graph 173 s; modular redesign/predictions graph 484 s (ended in analyzer failure); corrected analysis/cache recovery 14 s; distributed BindCraft2 peer campaign 329 s; validated RFD4 objective campaign 379 s. These intervals include scheduler queues and artifact collection and are not native GPU runtime measurements.

**Verification and final status.** Focused regressions included 59 campaign, transfer and saved-source tests and 24 CandidateSet tests, plus earlier 241 integration tests and 91 workflow/transfer tests. The corrected revision `e1ebb2c` is active and healthy on all three workstations, with the same source fingerprint `dd1c7388c8426c1045f381064081ecf628af5ad2cc2abbdd1225b16635ec57a7`. The earlier Protenix numerical-equivalence and per-row checkpoint-provenance limitations remain. Overall acceptance remains **partial**: all three demonstrated bounded scientific modes now have completed Results, but native multi-attempt batching, combined multi-generator objectives, broader generator adapters, cross-stage boolean groups, full incremental graph fan-out and live rejection-triggered continuation are not established. The deployment and documentation commits are local; no Git push or merge was performed.
