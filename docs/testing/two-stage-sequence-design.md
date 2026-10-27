# Two-stage sequence design and content-correct reuse — Phase 4

Validated on 27 October 2026 (Europe/Berlin), on `fix/scientific-pipeline-validation`.
Nothing was pushed. This implements the clarified F5 design and repairs the demonstrated F9 cache gaps in the [scientific audit](../predefined-pipeline-scientific-audit.md).

## Scientific behavior

The existing unified `sequence.mpnn` adapter already supported ProteinMPNN and SolubleMPNN, source-chain normalization, interface contact calculation, explicit residue selectors, and some fixed-motif protections. The missing parts were stage-aware routing, parent preservation, strict interface/mapping enforcement, general output verification, and complete cache input provenance. No new sequence-design engine was introduced.

The unified node now has three explicit modes:

| Mode | Behavior |
|---|---|
| `initial` | Backbone-only candidates receive ProteinMPNN sequences. Candidates already carrying sequence and structure retain their original sequence and ID, after correspondence checks. Mixed collections are partitioned, and only backbones go to inference. |
| `soluble_optimization` | Disabled by default. When enabled, the existing shared backend runs SolubleMPNN and retains every parent alongside its child variants. Interface fixing is mandatory. |
| `redesign` | Retains the existing deliberate, user-configurable model/redesign workflow. |

The all-generator comparison has an initial preparation node followed by the optional optimization node. Single-generator MPNN starters also use initial preparation, including co-design continuations. Untouched legacy v3 comparison defaults migrate to the two-stage graph; explicitly configured redesigns remain the user's choice. Existing RFD4 direct-prediction and VHH/SAbR pipelines remain available. To optimize either route in a custom graph, insert another unified MPNN node, select `soluble_optimization`, and enable it. Keep VHH optimization after its SAbR annotation.

**SolubleMPNN is trained on soluble proteins. Its output does not establish increased experimental solubility.**

The interface definition is a binder residue with any available non-hydrogen atom within **4 Å, inclusive**, of an available target non-hydrogen atom. The distance is configurable through `interface_cutoff` (greater than zero, at most 20 Å). Contacts are calculated on the candidate's source complex, before backbone completion or chain normalization, and mapped to the engine namespace. Ambiguous/missing chain assignments, mismatched sequences, unmappable protected positions, and unsupported native sequence conditions fail instead of falling back to unrestricted redesign.

All target chains remain fixed. Native fixed motifs, user selectors, native cysteines, sampled framework masks and SAbR framework/unselected-CDR constraints are retained. Constraints inherited through prediction are mapped using explicit binder chain order, residue counts and fixed identities. Native sequence conditions require a verified residue mapping unless they select only the already-fixed target.

For unified outputs, FASTA must agree with the designed structure; fixed binder positions and target sequences are checked after inference. LigandMPNN FASTA follows **sorted protein-chain IDs**, independently of PDB chain order; the adapter now follows that native convention and rejects inconsistent chain counts. Foundry output checks cover its explicit author/label chain mapping as well. RFD4 component metadata is rebuilt for each child, with non-protein context checked, so the prediction contract describes the redesigned sequence rather than the parent's components.

Each child retains `parents`, the complete `source_candidate`, effective model/seed/contact settings, fixed-residue evidence and chain-role metadata. Original and variant sequences coexist in candidate collections and Results. Existing stage-specific evidence selectors can filter on a variant's sequence-design evidence without borrowing scores from an incompatible parental sequence. Prediction, measurement and QC provenance remain separate.

## Cache contract

Cache contract 2 fingerprints actual file bytes, including cofold sequence/compound files, RFD4 input files, relative references in JSON/YAML conditions, source collections and their resolved structure artifacts. Same-size edits are detected without depending on timestamps. Output manifests verify contents as well as successful job status and artifact existence; entries without the new manifest are not trusted.

Execution specifications record the effective submitted workflow/parameters, local adapter/config identities, installed image IDs, and selected checkpoint contents/revisions. The relevant MPNN, RFdiffusion, RFD3/RF3, Boltz, RFD4, Complexa, ESM binder and AF2 identities are included where applicable. MPNN jobs additionally retain checkpoint and image identities in `artifacts/mpnn_settings.json`.

Unspecified stochastic seeds are resolved before lookup and submission. Reuse stores that resolved seed for the scientific specification and means reuse of a completed draw. `reuse_completed=False` requests fresh execution and draws new unspecified seeds; explicit seeds remain explicit. Per-node execution state records the seed and specification. Changing a sequence stage, fixed selector, cutoff, effective setting, checkpoint or upstream contents invalidates its descendants while leaving matching generation results reusable.

On resume, a completed or running child whose specification or output contents no longer match is superseded, never relabeled under the new key. Previous attempts remain recorded. Changed targets receive content-addressed staged filenames so older runs' target copies are preserved.

## Verification

Eight new regressions failed before implementation: stage routing/preservation, ambiguous interface fixing, configurable cutoff, and the three demonstrated input-content cache cases. The retained reproduction is `.tmp/phase4/logs/reproduction.log`.

| Check | Result |
|---|---|
| Full CPU suite | **1,443 passed, 28 skipped**, 195.68 s, no failures |
| Final focused sequence contracts, including native-cysteine mutation rejection | **16 passed**, 1.99 s |
| Resume and downstream variant-filter regressions | **2 passed**, 28.41 s |
| Final real GPU integration suite | **2 passed**, 111.85 s |

CPU skips cover opt-in GPU/container integrations and unavailable optional MD dependencies. Existing RFD4, VHH, prediction-configuration, measurement-provenance and selection regressions passed in the full sweep. One old expectation was updated deliberately: intermediate-tool fixed residues are now protected rather than discarded. A backbone-only test now requests sequence generation without interface fixing; requesting an interface lock for an unknown/mismatched sequence is explicitly tested as an error.

### Actual GPU executions

An installed **RTX 4090, 24 GB**, GPU index **0**, was used, with driver **595.91.07**. Two devices were inventoried; no 5090 was present. Models and images were already installed; no model downloads or all-generator campaign were run.

| Route | Actual work and checked outcome |
|---|---|
| RFdiffusion → ProteinMPNN → SolubleMPNN | One real RFdiffusion binder generation: 60-residue binder, PD-L1 A1–40 fragment, 15 diffusion steps. Real ProteinMPNN generated one sequence. Real SolubleMPNN generated two children, retaining the ProteinMPNN parent. All **8** interface residues and an additional fixed motif stayed unchanged; **27 and 25** permitted positions changed. |
| RFD4 co-design → original sequence → SolubleMPNN | One real RFD4 generation: 60-residue binder plus PD-L1 A1–40, 60 steps, batch 1, no native oracle/filter. Initial preparation retained the exact co-designed sequence and performed no MPNN inference. Real SolubleMPNN generated two children; all **16** interface residues and an additional fixed motif stayed unchanged; **40 and 38** permitted positions changed. |
| RFD4 parent + two variants → Boltz-2 → Results | All three were predicted with real Boltz-2, preserving their exact sequences, target/binder roles and parent links. The three candidates reached the Results collection and evidence export. |

Final MPNN sampling used seed **20261027**, ProteinMPNN temperature **0.1**, SolubleMPNN temperature **0.4**, one ProteinMPNN sequence and two SolubleMPNN sequences per parent, installed `v48_020` weights and `mn-ligandmpnn:latest` resolved to its immutable image ID. Boltz used the same seed, **1 sample, 1 recycle, 20 diffusion steps**, single-sequence mode and no target template. These small budgets test integration and data handling, not useful design quality.

Completed generation outputs from this phase were reused during downstream retests. The first RFdiffusion trial used custom-contig mode, which did not provide the intended binder-task roles and triggered a motif-mapping failure; the successful representative uses the explicit binder task. Other early failures exposed the missing preservation-job input argument and native FASTA ordering. They were fixed and inference rerun. A later RFD4 test-only assertion read a field from the wrong request level; its production redesign had succeeded. Logs and artifacts from these attempts were retained. Two PDBIO warnings report missing RFD4 B-factors during PDB staging; no confidence or binding claim is made from them.

Immutable image IDs, checkpoint digests, GPU inventory, final sequences, contact/fixed masks, effective settings, production run directories and transcript paths are in [the execution evidence](two-stage-sequence-design-2026-10-27.json). Native structures, FASTAs, commands, logs and Results remain under `.tmp/gpu-tests/`; CPU/GPU transcripts remain under `.tmp/phase4/logs/`.

## Files and commits

- Sequence contract/UI: `mpnn_input.py`, `mpnn_ui.py`, `pipeline_graph.py`.
- Unified/native sequence workflows: `engine/workflow_engine/workflows/mpnn.py`, `sequence_design.py`.
- Cache/orchestration: `cache_provenance.py`, `engine/workflow_engine/workflows/graph_pipeline.py`, and `pipeline_graph.py` fingerprints.
- Tests: `test_two_stage_mpnn.py`, `test_two_stage_mpnn_gpu.py`, `test_cache_provenance.py`, and updated `test_mpnn.py`, `test_pipeline_graph.py`, `test_rfdiffusion3.py`.
- Validation: this guide, the evidence JSON, and the link in `scientific-pipelines.md`.
- Sequence implementation commit: `688f9e7`.
- Cache implementation commit: `be5298c`.
- Final validation/documentation commit: the commit containing this report.

## Remaining boundaries

The contact mask uses atoms present in the candidate; it cannot infer missing sidechain contacts or establish the biological interface. Cysteine identity preservation does not validate disulfide geometry. Unverifiable method-specific sequence constraints and changed/unsupported non-protein context are rejected. ProteinMPNN and SolubleMPNN GPU checks used the shared installed backend; Foundry constraints/output handling received CPU regression coverage in this phase.

The configurable nine-generator comparison is retained, with heterogeneous native budgets/conditioning and a shared fan-in dependency. It is not an equal-budget method benchmark. Only the two representative generators were exercised here. Fresh execution currently applies to the graph; editing only a downstream sequence stage reuses matching upstream generation. An external MSA/database revision that a service does not report cannot be inferred from a container identity. No prospective binding benchmark, wet-lab validation or solubility measurement was performed.

## Rerun

```bash
conda activate mn-molecular-design-studio
python -m pytest -q
python -m pytest tests/test_two_stage_mpnn.py tests/test_cache_provenance.py -q
MDS_GPU_TESTS=1 MDS_TEST_GPU=0 python -m pytest tests/test_two_stage_mpnn_gpu.py -v -s
```

The GPU command generates one new representative sample per engine by default. For downstream-only replay, `MDS_PHASE4_RFDIFFUSION_GENERATION` and `MDS_PHASE4_RFD4_GENERATION` can point to the successful generation directories in the evidence JSON. A missing installed model/image is a skip; a launched inference or contract failure is a failure.
