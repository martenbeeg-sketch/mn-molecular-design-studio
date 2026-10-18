# Proteina-Complexa implementation handoff — 2026-10-17

Latest update: the full native task smokes have completed and a separate generator tile/canvas node has been implemented. See [Generator implementation and verification](proteina-complexa-generator.md) for the current state, completed checks, and remaining GPU verification. The readiness-only notes below describe the earlier implementation handoff.

## Status

Ready for GPU smoke testing after the user switches to GPT-6 Luna. **No Proteina-Complexa GPU inference or design jobs have been started or queued by this implementation task.** The scheduler and paused V7 campaign were left untouched.

The existing Proteina-Complexa tile now opens a dedicated green Studio page modeled on [Tamarind’s input page](https://app.tamarind.bio/tools/proteina-complexa). It includes all three native tasks, example inputs, project selection, target upload/fetch, task-specific selections, sampling settings, independent job submission, scheduler controls, readiness checks, and dedicated result rendering.

## Installed engine and task defaults

Image: `mn-proteina-complexa:latest`, image ID begins `9adf3e7f75c4`. The actual installed source is **Proteina-Complexa 1.0.0**, commit `32b71ae1d9a8767414eae3921cf35969430db82f`. The Dockerfile ARG in the container repository names a different revision; the tile reports the version inspected inside the installed image. No image rebuild was needed. The UI does not offer Tamarind’s Claude optimization because it is absent from this image.

| Task / API value | Native config | Example | Length | Designs/job | Search |
| --- | --- | --- | --- | --- | --- |
| Protein Binder / `protein_binder` | `search_binder_local_pipeline` | PD-L1, chain A, hotspots A37,A39,A49,A98 | 70–150 | 4 | Best of N, 2 candidates |
| Ligand Binder / `ligand_binder` | `search_ligand_binder_local_pipeline` | 5SDV ligand IAI, residue A803 | 80–140 | 2 | Best of N, 2 candidates |
| Motif + Ligand / `motif_ligand` | `search_ame_local_pipeline` | M0024_1nzy_v3, ligand L:1 at A272, six motif residues | 120–220 | 4 | Single Pass |

Sampling defaults are 400 steps, guidance 1.0, noise 0.1, GPU batch size 1, deduplication on, and seed 0 (random, resolved value recorded). Number of jobs counts independent scheduler jobs; designs per job is the native requested design count. Native sequence redesign/evaluation settings are preserved per task. Motif + Ligand has no configured search reward and therefore only offers Single Pass.

The original input is preserved. A job-owned prepared copy normalizes ligand naming and chain IDs for the native ligand/motif pipelines, retains bonds, and remaps the selected motif residues while keeping the requested atom selections. Current address validation requires letter chain IDs and selections without insertion codes. PDB and mmCIF are accepted; fetched structures use the shared Studio conversion where needed.

## Shared resources

All six generator/autoencoder checkpoints are present under `/mnt/db/reference_files/proteina-complexa/ckpts`:

- `complexa.ckpt`, `complexa_ae.ckpt`
- `complexa_ligand.ckpt`, `complexa_ligand_ae.ckpt`
- `complexa_ame.ckpt`, `complexa_ame_ae.ckpt`

Other mounts:

- RF3: `/mnt/db/reference_files/foundry/rf3_foundry_01_24_latest_remapped.ckpt` mounted as `/rf3.ckpt`.
- AF2: `/mnt/db/reference_files/alphafold_models` (15 parameter NPZ files).
- Complete shared LigandMPNN model folder: `/mnt/db/reference_files/protpardelle-1c/model_params/LigandMPNN`.
- Cached ESM2/ESMFold models and persistent JAX/CUDA cache: `/mnt/db/reference_files/proteina-complexa/hf-cache`.
- Analysis tools: `/mnt/db/reference_files/proteina-complexa/tools`. Foldseek and MMseqs were copied from the already installed `mn-protpardelle-1c:cu128` and `mmseqs-gpu:latest` images. `provenance.json` records image IDs and source paths. Both execute their version commands successfully inside the Complexa image on CPU. The existing image’s bundled shape-complementarity executable is used directly.

The container helper writes only explicit runtime paths to a job-local `.env` for the native CLI, composes the installed task configuration, and launches `complexa design`. The existing `generator_only` API flag now launches `complexa generate`. No inherited Studio credentials are written into the job dotenv.

## Verification completed

- Full repository pytest: **305 passed**, one existing Starlette/httpx deprecation warning, 5.00 s.
- 29 Proteina-Complexa tests cover task inputs, malformed values, API preflight/schema, scheduler command construction, resource checks, input preservation, native task configs, environment isolation, native success-table handling, selected passing MPNN result structures (including multi-redesign list columns and joint motif criteria), ligand/protein chain roles, raw candidates never passing filters, and all three task forms.
- All three installed native configurations compose and resolve inside `mn-proteina-complexa:latest`.
- All three native conditional-feature readers instantiate and set up successfully on CPU using the bundled examples.
- `complexa validate design` passes for all three prepared configs (12, 13, and 9 checks respectively). These checks do not load models or run predictions.
- Live browser: tile, three task forms, example loading, ligand/motif selections, optional controls, and readiness verified. Preview: [input page](images/proteina-complexa-input.jpg). Submit remains disabled without a selected project.
- Both local API servers were reloaded with the updated schema. UI: `http://127.0.0.1:8502/`; API: ports 8766 and 8767. Do not restart the shared scheduler or resume V7.

CPU check artifacts and logs are in `/tmp/mds-complexa-preflight/{protein_binder,ligand_binder,motif_ligand}` (`prepare.log` and `validation.log`). `/tmp/complexa_config_preflight.py` captures the runner command with job-creation/status functions replaced by local scratch-folder stubs, removes Docker GPU arguments, and runs `--check-features`. The helper also supports `--prepare-only` and `--validate-only`. No jobs are created in Studio by this script.

## GPU smoke test to perform after the model switch

Use the ordinary Studio API and shared scheduler, workflow ID `design.proteina_complexa`, in a dedicated selected smoke project. Use the bundled example for each task. Prefer one independent job with one design, batch size 1, 400 sampling steps, seed 17. Keep native length ranges and the task’s native search strategy (Best of N with 2 candidates for protein/ligand; Single Pass for motif). Start with a single task before submitting the remaining tasks so runtime problems are caught early. GPU assignment can remain automatic; the scheduler controls GPU use.

Parameters can be built with `complexa_ui.example_parameters(task)` and overridden with `num_designs=1`, `seed=17`, `batch_size=1`, `cpu_cores=4`, a campaign name, and `gpu_device="all"`. Copy the packaged example into the submitting account’s private input folder before API submission; the UI’s Load example inputs action already does this. The endpoint is `POST /api/v1/workflows`; use the existing project/auth contract. `POST /api/v1/workflows/validate` validates without submitting.

Verify generation, native filtering, refolding/evaluation, final RAW/success CSVs, normalized candidates, result viewers/downloads, and final scheduler status. Native task-specific `all_successes_*.csv` tables are authoritative. A smoke can complete with zero filter-passing designs; do not fabricate passing designs. Generated structures remain in a separate pool and are never eligible for successful-design validation batches. A successful MPNN sequence must use its matching predicted structure and metrics. The normalizer reads the native success-criteria JSON to pick a passing redesign from `_all` list columns; it does not assume the engine’s default best candidate necessarily passed. If a passing redesign cannot be verified or its evaluated structure is missing, it is not eligible for successful-design validation batches. GPU runtime, model-loading compatibility, and design quality remain unverified until these smoke tests finish.

## Implementation locations

- `mn_molecular_design_studio/complexa_input.py`: validation, task defaults, resources, native config overrides, staging.
- `mn_molecular_design_studio/complexa_container.py`: job-local preparation and existing-image CLI launcher.
- `mn_molecular_design_studio/complexa_ui.py`: task input page and results.
- `mn_molecular_design_studio/engine/workflow_engine/workflows/design.py`: scheduler-compatible runner and result normalization.
- `mn_molecular_design_studio/engine/workflow_engine/services/workflow_automation.py`: schema and API input preflight.
- `mn_molecular_design_studio/run_app.py`: tile, version, setup/results routing.
- `mn_molecular_design_studio/pipeline_graph.py`: compatible protein-binder node sampling controls. The standalone tile exposes all three tasks; the existing canvas target port remains protein-focused.
- `tests/test_complexa.py`: CPU regressions.
- `mn_molecular_design_studio/assets/proteina_complexa`: packaged example coordinates and provenance.

Preserve this checkout; much of the existing workspace is untracked. Do not reset it.
