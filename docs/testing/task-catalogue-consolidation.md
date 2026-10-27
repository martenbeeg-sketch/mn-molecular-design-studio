# Task-oriented scientific catalogue — Phase 5

Implemented on 27 October 2026 (Europe/Berlin), continuing
`fix/scientific-pipeline-validation` from `47700a7`. Nothing was pushed.
The previous F1–F9 implementations were inspected in their committed source and
retained execution evidence, then 110 focused CPU regressions passed before edits.
The requested audit sections D, E, G2 and G3 informed the task routing.

## Catalogue and scientific methods

The old primary menu had 21 scientific recipes plus a blank canvas, including
all-generator and all-compatible-predictor presets. The primary menu now has eight
scientific objectives plus a blank canvas. Those objectives compose conditional,
editable graphs rather than forcing eight fixed workflows.

| Scientific objective | Methods available in task setup |
|---|---|
| Targeted Protein Binder Design — standard | RFdiffusion3, RFdiffusion, RFD4, BindCraft2 native pipeline or generator, BoltzGen native pipeline or generator, PXDesign generator or native pipeline, Genie3 generator or native pipeline, Proteina-Complexa co-design generator or native pipeline, ESMFold2 Binder Design, RFdiffusion beta pairing |
| Antibody CDR Redesign — standard | Existing candidates → SAbR → AbMPNN with selected CDRs/framework/cysteine protection; sequence QC and primary prediction |
| Template-Based Peptide / IDR Recognition — standard | Peptide Trap and Kejia IDR as distinct strategies; Kejia retains its template/threading controls and initial sequence stage |
| Molecular Structure Prediction — standard | Capability-filtered registered predictors for monomer, protein complex, ligand, DNA, RNA and cyclic partner modes; Boltz-2 is the single default |
| Small-Molecule Hit Prioritization — standard | Existing editable adaptive docking campaign, including its docking-engine controls; direct docking remains in Tools |
| New Protein Fold Design — standard | Unconditional Protpardelle generation → initial ProteinMPNN → monomer self-consistency; no target-binding claim |
| Ligand-Conditioned Protein-Pocket Design — experimental | RFD4 small-molecule binder co-design; exact ligand/atom constraints are supplied by the user |
| Catalytic Motif Scaffolding — experimental | RFD4 enzyme motif scaffolding; user-defined motif/active-site constraints, with no catalytic-activity claim |

The protected VHH format exposes BindCraft2 VHH generation, BoltzGen nanobody
generation and ESMFold2 VHH generation individually or as an explicit comparison.
SAbR annotation follows generation and precedes requested soluble redesign.
Unrestricted methods cannot be selected through the VHH task configuration.
Native complete workflows keep their own internal generation, inverse folding,
refinement and filtering. No external MPNN, predictor or Rosetta stage is appended
by default. Requested sequence variants receive subsequent prediction; additional
validation of an unchanged native output is explicit.

Prediction choices are derived from the existing molecular capability profiles:

- Monomer: RF3, Boltz-2, official AF3, AF3 PyTorch/OpenFold3 P2, AF3 JAX/OpenFold3 P2,
  Protenix V1/V2, ESMFold, ESMFold2 Full/Fast, Chai-1, AtlasFold-M,
  OpenFold3 OpenBind-0, OpenDDE and ColabFold.
- Protein complex: the above, plus AF2 target-only evaluation and BoltzGen Fold
  when a target PDB is supplied. ESMFold complex use remains explicitly exploratory.
- Ligand/DNA/RNA: RF3, Boltz-2, official AF3, both OpenFold3 P2 implementations,
  Protenix V1/V2, Chai-1 and OpenFold3 OpenBind-0. Chai-1 is omitted for CCD ligand
  notation; it accepts SMILES.
- Cyclic protein/peptide partner: Boltz-2, using the native cyclic-chain capability.

These are compatible registered adapters, not a claim that every checkpoint is
installed on every host. Existing engine preflight still controls installed
images/checkpoints, terms, supported chemistry and task-specific requirements.
The API catalogue also exposes per-mode and per-binder-format method discovery.
No engine was removed from Tools or the canvas node library.

## Configuration and execution

Task setup follows objective → molecular inputs → method → optional stages →
**Create editable pipeline**. Selecting a different objective does not overwrite
an existing graph. Source and engine nodes still open the existing detailed forms;
advanced task JSON is also available. Standard/experimental status, method
limitations, samples/seeds and encoded workload counts are visible. Runtime is
not estimated without reliable evidence.

New recipes use one attempt, one initial sequence and one prediction sample/model
by default, respecting engine-specific minimum sample counts. Three recycles are
used where the adapter exposes them. No global model settings or saved legacy
budgets were changed. Beta-pairing conditions, Kejia templates, Peptide Trap
registers, native acceptance budgets and sequence counts retain their distinct
meaning; equal attempt counts are not an equal GPU-time benchmark.

Comparisons require an explicit flag and at least two distinct selected methods.
Each method has separate input/context, generation, sequence, prediction,
optional selection and Results nodes. Settings stay on the individual nodes;
execution specifications retain effective requests, resolved seeds, runtime
identities and cache evidence. Candidate IDs are namespaced by Results branch.
Final candidates record the actual ancestor method nodes/settings, in addition
to the existing source-candidate lineage and graph snapshot.

A failed generator, prediction or method-specific request preparation blocks its
own descendants. Other branches continue. A comparison with a healthy complete
Results branch is completed with `comparison_status=partial` and explicit failed
node details; all-failed comparisons remain failed. Failed measurement rows can
reach review Results but do not qualify as a successful branch. Legacy/custom
graphs without explicit task comparison configuration retain their prior overall
failure rule. Cache lookup, seed resolution, scheduler submission and output
fingerprints use the existing mechanisms.

Optional selection is the existing evidence-qualified filter, with thresholds,
required finite metrics, manual IDs and diversity controls. Ranking remains
partitioned by recorded comparable protocols. Extra validation requires an
explicit shortlist stage. The default Boltz route reuses its native ipSAE output;
requesting a separate ipSAE node disables the redundant Boltz calculation.
Rosetta interface analysis is optional.

Boltz nodes now expose a seed and target-MSA control. Historical defaults retain
target MSA; an explicit single-sequence choice is submitted and recorded rather
than overridden. Samples/recycles/full-PAE semantics remain the Phase 3 controls.

## Sequence strategy and scientific constraints

Backbone-only methods receive the Phase 4 initial ProteinMPNN stage. Native
sequence–structure co-design outputs bypass generic redesign and preserve their
original sequences. The initial preparation adapter still preserves native
candidates if an engine emits a mixed collection.

SolubleMPNN is disabled by default. Enabling it adds the existing
`soluble_optimization` stage with mandatory interface fixing and parent retention.
The configured heavy-atom contact cutoff, fixed motifs, target sequences,
framework masks, cysteines and verified residue mappings remain enforced by the
Phase 4 implementation. Parents and children coexist through prediction,
selection, Results and exports. SolubleMPNN does not guarantee improved
experimental solubility.

New experimental RFD4 recipes omit the legacy tutorial's 5QA atom and A36/A41/A87
motif selectors. The user supplies constraints for the actual ligand/motif in
setup or the existing RFD4 form. Annotated CIF and native condition-spec routes
remain available. Native chemistry/condition validation and the existing strict
component handoff are unchanged.

## Saved workflows, IDs and APIs

The graph schema remains `molecular-design-studio.graph.v1`. Navigation does not
rewrite saved nodes, execution settings, source candidate references, immutable
versions or previous results. Only the already-established Phase 1–4 draft
migrations remain. Task configuration is creation/provenance metadata; the saved
editable nodes and connections determine execution.

All 34 former factory identifiers remain supported, including the deliberately
unconnected blank canvas:

- Former menu IDs: `unidock_active_learning`, `foundry`, `kejia_idr`, `peptide_trap`,
  `abmpnn_cdr_redesign`, `pxdesign`, `rfdiffusion`, `genie3`, `protpardelle`,
  `esmfold2_binder`, `all_generators_boltz2`, `vhh_design`, `rfd4_protein_binder`,
  `rfd4_small_molecule`, `rfd4_enzyme`, `cofold_protein`, `cofold_protein_protein`,
  `cofold_protein_ligand`, `cofold_protein_dna`, `cofold_protein_rna`,
  `cofold_protein_macrocyclic`, `blank`.
- Native/specialist factory IDs: `bindcraft2`, `boltzgen`, `complexa`,
  `rfdiffusion_beta_pairing`.
- Continuations: `complexa_generated`, `rfdiffusion_generated`, `genie3_generated`,
  `pxdesign_generated`, `protpardelle_generated`, `kejia_idr_generated`,
  `peptide_trap_generated`, `rfdiffusion_beta_pairing_generated`.

The eight new IDs are `targeted_protein_binder`, `antibody_cdr_redesign`,
`template_peptide_idr`, `molecular_structure_prediction`,
`small_molecule_hit_prioritization`, `new_protein_fold`,
`ligand_conditioned_pocket` and `catalytic_motif_scaffolding`.

`GET /api/v1/pipelines/templates` preserves its `items` shape for the primary
catalogue and adds `tasks` and `legacy_items`. `POST /api/v1/pipelines/drafts`
accepts old IDs and adds a `configuration` object for task IDs. The MCP
`create_pipeline_draft` tool forwards the same object. Project/file ownership
checks remain in place. Existing Kejia/Peptide Trap `source_job` helpers and all
UI continuation constructors remain; existing collections can be resumed without
rerunning generation. The new task configuration does not broaden the legacy
`source_job` helper to every generator.

Example unsaved API task configuration:

```json
{
  "template": "targeted_protein_binder",
  "project_id": "existing-project-id",
  "configuration": {
    "methods": ["rfd3", "rfd4"],
    "comparison": true,
    "num_designs": 4,
    "inputs": {"target_pdb": "authorized-target-path", "target_chains": "A"},
    "soluble_mpnn": true,
    "soluble_settings": {"sequences": 1, "interface_cutoff": 4.0},
    "method_settings": {"rfd3": {"seed": 42}, "rfd4": {"seed": 43}}
  }
}
```

## Validation

See [machine-readable evidence](task-catalogue-consolidation-2026-10-27.json).
CPU tests include every task and advertised method, all molecular modes,
engine-specific sampling minima, capability/input rejection, protected VHH
routing, conditional sequence design, legacy graph round trips, API ownership,
continuations, separate comparison identities, preparation/engine/measurement
failure handling, and Streamlit creation of every task and blank canvas.

The final full CPU sweep passed **1,549 tests**, with **30 expected skips**, no
failures, in **256.41 seconds**. The final task/UI suite passed **108 tests** in
54.73 seconds, including three additional comparison/UI cases added after the
full sweep's collection. A focused graph/API/configuration sweep passed 168 tests.
Skip causes were opt-in GPU/container tests and unavailable optional MD dependencies.
The commands and retained results are recorded in the evidence JSON.
The previous scientific-contract, antibody, AF2, Boltz, cache, measurement,
candidate-selection, Results and continuation regressions remain covered.

Two real GPU task integrations passed in 215.78 seconds. They reused the retained
real RFdiffusion and RFD4 generation outputs from Phase 4, then exercised the new
recipe nodes and production request adapters. They did not rerun generation:

| Route | Fresh inference and verified outcome |
|---|---|
| Retained RFdiffusion backbone → ProteinMPNN → SolubleMPNN → Boltz-2 → selection → Results | One initial sequence, original plus one soluble child; both exact sequences reached prediction and Results; mapped interface/cysteines stayed fixed |
| Retained RFD4 co-design → SolubleMPNN → Boltz-2 → selection → Results | No initial redesign node; exact native sequence retained alongside one child; both sequences reached prediction and Results with source lineage and protected interface/cysteines |

Installed RTX 4090 GPU 0 was selected, with driver 595.91.07. Both tests used seed
20261027, one soluble variant, temperature 0.4, and Boltz one sample/one recycle/
20 diffusion steps with target MSA and target templates explicitly disabled.
Installed images/checkpoints were used; no model downloads or all-engine campaign
were run. Commands, immutable runtime identities, parameters, candidate collections
and structures remain under the artifact paths in the evidence JSON.

The first GPU attempt failed before inference because the new test harness passed
an unsupported `cpu_cores` argument to the MPNN function. That harness error was
removed and both integrations passed on rerun. Two PDBIO warnings arose while
staging retained RFD4 nonfinite B-factors; they do not provide confidence evidence.
Initial CPU comparison tests also used the wrong execution-specification field
name; that assertion was corrected. A later capability probe exposed BoltzGen
Fold's two-sample minimum, which is now retained by recipe and UI defaults.

## Files and commits

- `3751557`: catalogue/routing and independent comparison execution —
  `scientific_tasks.py`, `pipeline_graph.py`,
  `engine/workflow_engine/workflows/graph_pipeline.py`.
- `24712cc`: task-first UI and optional stage controls — `scientific_tasks_ui.py`,
  `pipeline_graph_ui.py`, `pipeline_node_ui.py`.
- `f550799`: API/MCP compatibility — `workflow_setup_api.py`, `mcp_server.py`.
- Tests/documentation: the commit containing this report —
  `tests/test_scientific_tasks.py`, `tests/test_scientific_tasks_gpu.py`,
  `docs/visual-pipelines.md`, this report, the evidence JSON and the link in
  `docs/testing/scientific-pipelines.md`.

## Boundaries

No prospective binding, catalytic, solubility or experimental success validation
is claimed. Existing chain/chemistry, Rosetta one-target/one-binder and mapped
interface restrictions remain. The native complete methods were preserved and
CPU-tested as routes, not all rerun on GPUs. DNA/RNA/cyclic modes received capability
and graph coverage rather than new GPU benchmarks. Comparison budgets are explicit
attempt/sample counts rather than calibrated equal GPU time. No reliable universal
runtime estimator or cross-engine ranking calibration was introduced.

Initial sequence preparation uses ProteinMPNN; deliberate alternate inverse-folding
engines remain in Tools/the editor. Experimental chemistry and motif constraints
still require specialist input. Advanced setup settings use the existing node
forms or JSON. No docking, MD, screening implementation, distributed execution
framework, frontend framework or model algorithm was changed.

## Rerun commands

Use the existing Conda environment from the repository:

```bash
conda activate mn-molecular-design-studio
CUDA_VISIBLE_DEVICES='' python -m pytest -q
CUDA_VISIBLE_DEVICES='' python -m pytest tests/test_scientific_tasks.py -q
```

For downstream-only real GPU replay, set `MDS_TASK_GENERATION_RFDIFFUSION` and
`MDS_TASK_GENERATION_RFD4` to the retained successful generation directories in the
evidence JSON, then run:

```bash
MDS_GPU_TESTS=1 MDS_TEST_GPU=0 python -m pytest tests/test_scientific_tasks_gpu.py -v -s
```

Missing replay inputs or installed prerequisites skip; launched inference errors
fail. The GPU tests are automatically skipped during ordinary CPU sweeps.
