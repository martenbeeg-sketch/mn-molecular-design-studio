# AF2 validation and Boltz-2 configuration — 27 October 2026

Implemented on `fix/scientific-pipeline-validation`, starting from clean commit `fd11ffb`. The existing F1/F2 component and sequence contracts and F3/F7 measurement/selection infrastructure were retained. Nothing was pushed.

The implementation commits are `4caa1f4` (AF2) and `6154f13` (Boltz-2). The integration tests and this report are in the subsequent test/documentation commit; obtain its hash with `git log -1 --format=%h -- tests/test_prediction_configuration_gpu.py`.

Seven initial CPU regressions failed before repair: a requested AF2 target-only mode was forced back to whole-complex initialization; legacy AF2 results lacked explicit conditioning metadata; Boltz aliases did not change sample/recycle requests; conflicting aliases were accepted; duplicate sample/recycle controls remained visible; and PAE handling did not match configuration. The transcript is `.tmp/phase3/logs/initial-reproductions.txt`.

**AF2 modes and interpretation**

| Setting | `conditioned_self_consistency` | `target_only` |
|---|---|---|
| Staged prediction PDB | Designed target and binder coordinates | Selected target coordinates only |
| Native ColabDesign protocol | `binder`, redesign of the staged binder | `binder`, new binder of the supplied sequence length |
| `use_initial_guess` | `True` | `False` |
| `prep_inputs` binder argument | Explicit designed binder chain | `binder_chain=None`, `binder_len=len(sequence)` |
| Binder sequence | Checked against the staged designed binder | Supplied as sequence without binder coordinates |
| Target template | Target coordinates, sequence and side chains retained | Target coordinates, sequence and side chains retained |
| Binder/interface templates | Optional explicit settings; default off | Rejected |
| Designed reference | Retained for post-hoc structural comparison | Retained separately for post-hoc comparison; never passed to `prep_inputs` |
| Meaning | Pose self-consistency conditioned on the design | Less-conditioned complex prediction; still conditioned on target conformation/template and any supplied MSA |

The existing installed ColabDesign implementation supports both native preparations; no substitute predictor or approximate removal of binder coordinates was introduced. Whole-complex initialization remains useful evidence of conditioned self-consistency and is not described as independent binding-orientation verification. Neither mode provides an experimental binding probability.

Canvas nodes expose mode selection and descriptions. Candidate-based saved AF2 graphs without a mode retain whole-complex initialization. Saved direct target-plus-sequence graphs retain target-only behavior. New canvas screening nodes default to conditioned evaluation; the selected-design validation page retains its existing target-only default and now offers both modes. A direct sequence input cannot request conditioned evaluation without a designed complex. The existing monomer/fixbb path remains separate and identifies its target-coordinate initialization.

A node runs one selected mode. Users can screen with conditioned evaluation, filter/select candidates, then run a target-only node and compare Results evidence. The original generated reference survives that staged handoff. Target-only result IDs have a separate suffix so the two observations cannot overwrite each other. No automatic second AF2 evaluation or universal acceptance threshold was added.

Both modes retain chain maps and compare predicted structures against the retained design using the existing comparison implementation: target-fitted binder C-alpha RMSD, target alignment RMSD, matched atom counts and coverage. Confidence metrics, model-specific structures and PAE remain available. Failed/missing measurements use the existing measurement-outcome contract. Output chain/sequence mismatches fail validation; missing output cannot be substituted with the input structure. The binder-alone auxiliary prediction remains part of the existing runner.

Native execution manifests record ColabDesign/JAX versions, AF2 parameter names and checkpoint SHA-256 hashes, effective seeds/model seeds, model and recycle counts, input PDB hash, exact sequences and chains, prepared binder coordinate atom counts, template masks and MSA row counts. Results also retain installed image identity, reference, source candidate, prediction protocol and conditioning description. Evidence rows and common-metric provenance expose conditioning and effective settings. Graph cache keys include mode/settings and configuration/native-runner implementation digests.

**Boltz-2 settings and migration**

| Configuration | Effective behavior |
|---|---|
| `samples` | `samples_per_item` in the batch manifest, then native `--diffusion_samples` |
| `num_recycles` | Batch option `recycling_steps`, then native `--recycling_steps` |
| Legacy `diffusion_samples` / `recycling_steps` in request-builder clients | Mapped to canonical names; unequal simultaneous values raise a conflict error |
| Old saved canvas with shipped inactive defaults (`diffusion_samples=3`, `recycling_steps=10`) alongside canonical settings | Removes the inactive controls, preserves the previously effective canonical values, and records visible migration notices |
| Old saved canvas with only inactive aliases | Preserves its historical effective defaults (one sample, three recycles) with notices identifying the discarded inactive values; use canonical controls to change execution explicitly |
| Old saved canvas with edited contradictory aliases | Rejected with an explanation of the historical behavior and the canonical setting to resolve |
| `write_full_pae` | Fixed `True` for the installed engine; no active toggle; new `False` requests are rejected with the limitation explained |
| Old saved canvas `write_full_pae=False` | Migrates to its historical actual output, full PAE, with a notice |

Migration is idempotent and saved with configuration version 2. Editing a new node uses version 2 immediately. Canonical sample/recycle changes alter both the request and cache identity; explicit zero values cannot silently restore defaults. Adopted run graphs describe effective settings without recreating duplicate controls.

Real GPU testing exposed a deeper PAE limitation after the initial adapter fix: omitting `--write_full_pae` still produced PAE in the installed Boltz-2 model path. The CLI option is not a reliable output switch for this image. The active toggle was therefore removed, rather than claiming that `False` worked. Full PAE is explicitly requested and retained for built-in or downstream ipSAE. The retained diagnostic test with `False` is an initial failed test, not a supported final configuration.

Results include effective native batch settings, actual native command lines/return status, engine version/image digest, exact input artifacts, source/chain assignments, sample-specific structures/confidence/PAE, and checkpoint/input content identities. The sampling and recycle values in result protocols are the effective values used in execution.

**Verification**

CPU full sweep: **1,415 passed, 26 skipped, 0 failed**, 119.12 seconds. After final control/identity and zero-value fixes, the affected regression set passed **115 tests**; checkpoint-provenance/queue/component regressions passed **141 tests**. GPU tests are skipped in the ordinary CPU suite. Optional missing dependencies/preparation prerequisites account for the other skips. Existing deprecation/numerical-fixture warnings remain.

Six real GPU pytest cases passed, using installed weights and two RTX 4090 devices (24 GB, driver 595.91.07). Inference seed was 20261027. No download, numerical-equality requirement, large campaign or benchmark was used.

| GPU test | Actual result |
|---|---|
| AF2 conditioned mode | ColabDesign 1.1.3, JAX 0.6.2, `model_1_multimer_v3`, one recycle/model. Native input retained T/V and **166 binder coordinate atoms**; binder template mask off. Predicted T/V sequences and target-aligned pose comparison verified. |
| AF2 target-only mode | Same parameter family/count/recycles. Native PDB contained T only; `binder_chain=None`, redesign false, **zero binder coordinate atoms**, binder template mask off; predicted T/V identities and target-aligned pose comparison verified. |
| Boltz-2 one sample / one recycle | Native `--diffusion_samples 1 --recycling_steps 1 --write_full_pae`; one predicted complex with sample-specific PAE. |
| Boltz-2 two samples / two recycles | Native `--diffusion_samples 2 --recycling_steps 2 --write_full_pae`; two predicted complexes with separate confidence/PAE artifacts. |
| RFD4 candidate schema → Boltz-2 → ipSAE → Rosetta | Passed real prediction and both analysis containers; source T/V → prediction A/B; existing experimental 7EOW fixture encoded in the RFD4 schema, not newly generated by RFD4. |
| Prediction → ipSAE → sequence QC → selection | Passed real prediction/analysis/QC; missing-structure control remained invalid and could not pass selection. |

The AF2/Boltz configuration smoke tests use 24-residue fragments from each protein chain of the existing 7EOW reference (48 residues total), retaining its source pose and explicit T/V roles. These deliberately small fragments test execution/conditioning, not binding-quality prediction. The compatibility tests use the existing full reference fixture. Boltz smoke tests use 20 diffusion steps and single-sequence mode without target templates. AF2 uses single-sequence inputs and also runs its existing binder-alone auxiliary fold. Low-confidence predictions and large pose deviations in this fragment fixture are retained, not filtered away.

Machine-readable [execution evidence](af2-boltz-configuration-2026-10-27.json) identifies retained jobs, GPU/image identities, native feature manifests/commands, chains/sequences, structures, effective settings and observed metrics. Full logs, raw/normalized structures, source references, PAE, selection/analysis outputs and execution manifests remain under `.tmp/gpu-tests/`; pytest/debug transcripts remain under `.tmp/phase3/logs/`.

**Changed files**

AF2 commit `4caa1f4`:

- `prediction_configuration.py`, `prediction_provenance.py`: shared conditioning semantics and content/container identities.
- `pipeline_graph.py`, `pipeline_node_ui.py`, `run_app.py`: mode selection/descriptions, compatibility defaults and effective controls.
- `engine/workflow_engine/workflows/graph_pipeline.py`, `design_validation.py`, `refolding.py`: requests, distinct staging, exact sequence/chain checks, retained references, comparison validity and manifests.
- `engine/workflow_engine/tools/af2_initial_guess_binder_eval.py`: native feature/checkpoint/seed provenance.
- `output_vocabulary.py`: conditioning/effective settings in result evidence and metric provenance.
- `tests/test_af2_configuration.py`, `tests/test_rfdiffusion_beta_pairing.py`: regressions and corrected labels.

Boltz-2 commit `6154f13`:

- `prediction_configuration.py`, `pipeline_graph.py`, `pipeline_node_ui.py`: canonical settings, saved-graph migration/notices, duplicate-control removal and fixed PAE explanation.
- `engine/workflow_engine/workflows/graph_pipeline.py`, `optimized_prediction.py`, `pipeline_adoption.py`: request/native parameter mapping, conflicts/PAE validation, effective settings, native execution and checkpoint provenance.
- `tests/test_boltz_configuration.py`, `tests/test_pipeline_adoption.py`: alias/conflict/cache/PAE/zero-value regressions and adoption contracts.

The subsequent test/report commit contains `tests/test_prediction_configuration_gpu.py`, this report, its JSON evidence, and the link from `scientific-pipelines.md`.

**Remaining boundaries**

AF2 complex evaluation supports one binder chain with explicit target roles. A target fragmentation that needs an additional residue correspondence map is reported as unavailable pose agreement rather than silently using an incorrect chain map. Template conformation, target sequence/side chains, parameter-family priors and optional MSA are residual conditioning in target-only mode. The computational tests do not establish prospective predictive validity or experimental affinity. ipTM, ipSAE, pLDDT and structural agreement are model-derived evidence, not binding probabilities.

Full PAE cannot be disabled in the installed Boltz-2 path. Changing that behavior would require a separately verified engine change. Stochastic seed zero still means random at execution; reusing a completed graph node reuses the recorded result rather than requesting fresh sampling. Use a recorded nonzero seed and rerun when reproducing an inference. No new full generator campaign was run; existing F1/F2/F3/F7 CPU regressions and the affected downstream GPU integrations were exercised.

**Repeat commands**

From the repository, using the existing Conda environment:

```bash
conda activate mn-molecular-design-studio
python -m pytest -q
python -m pytest tests/test_af2_configuration.py tests/test_boltz_configuration.py tests/test_pipeline_graph.py tests/test_pipeline_adoption.py tests/test_rfdiffusion_beta_pairing.py tests/test_design_validation.py -q
MDS_GPU_TESTS=1 MDS_TEST_GPU=0 python -m pytest tests/test_prediction_configuration_gpu.py -k af2 -v -s
MDS_GPU_TESTS=1 MDS_TEST_GPU=1 python -m pytest tests/test_prediction_configuration_gpu.py -k boltz -v -s
MDS_GPU_TESTS=1 MDS_TEST_GPU=1 python -m pytest tests/test_scientific_gpu.py::test_rfd4_candidate_boltz2_interface_analysis tests/test_scientific_gpu.py::test_real_prediction_to_qc_and_selection -v -s
```

Run the GPU commands sequentially on each selected device. `MDS_GPU_ARTIFACTS` can choose the retained artifact directory. Missing local images/weights are skips; launched inference or contract failures are failures. No model download is authorized by these tests.
