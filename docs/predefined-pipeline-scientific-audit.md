# Scientific and architectural audit of predefined pipelines

**Project:** MN Molecular Design Studio  
**Audit date:** 27 October 2026 (Europe/Berlin)  
**Source revision:** `d56667628cc708b960c13824d8406f8f6a8bcd35`  
**Mode:** Read-only review of application code and pipeline definitions. This report is the only repository change. No scientific GPU jobs, model downloads, docking screens, or MD simulations were launched.

## A. Executive summary

The catalogue contains useful scientific building blocks, but it is not yet a uniformly validated collection of experiment-ready workflows. Its strongest feature is explicit separation of generation, sequence design, prediction, and analysis, with candidate lineage and scheduler integration. Its main weakness is that **a connected graph and a completed job often establish computational completion without establishing a defensible candidate-selection decision**.

The public starter registry contains **22 entries: 21 scientific starters and one blank canvas**. The graph factory also implements **12 additional named graphs**: four non-menu factory graphs and eight continuation graphs. These are inventoried separately below. Individual workflow adapters, native pipelines inside tools, configurable campaigns, and structural scaffold libraries are not counted as additional canvas templates.

The most important demonstrated findings are:

1. **The VHH starter fails graph validation as shipped.** Its BindCraft2 branch intentionally has an empty binder-length field for scaffold mode, but the graph validator rejects that value. After that is repaired, the downstream general ProteinMPNN step still needs antibody-specific protection: it has no CDR selection, fixes neither the framework nor the interface by default, and excludes cysteine at mutable sites. The validation failure was reproduced; sequence damage was not experimentally tested. [F1; E01, E04, E11]
2. **RFD4-to-Boltz-2 conversion loses biological chain roles.** RFD4 candidates carry `studio_molecule_id` and molecular components. The component conversion path serializes all proteins into `protein_sequences`, omitting target/binder partitioning, the designed binder sequence, and per-chain MSA controls. A CPU-only reproduction produced empty normalized target and binder role lists. This compromises the protein-binder starter's promised interface analysis and affects the interpretation/export of pocket and enzyme designs. Ligands are not simply discarded: they are carried as component identities, with separate chemical-fidelity limitations. [F2; E06, E09, E10]
3. **Most explicit starters do not select a final shortlist.** No public scientific starter includes the general `filter` node. ipSAE and Rosetta analysis append measurements; Results collects rows and may calculate pose RMSD, but does not establish an overall acceptance rule. Native end-to-end tools can have their own filtering, and Peptide Trap rejects unsuccessful native designs; those are important exceptions. [F3; E01, E05, E07, E13, E16]
4. **AF2 evaluation is more conditioned than several descriptions suggest.** Canvas candidate adapters force whole-complex initial-guess mode, despite disabling binder and interface templates. That can be useful for pose self-consistency, but is not equivalent to an unconstrained, independent search for the binding pose. [F4; E08]
5. **The “all generators” default spends computation without a controlled comparison.** Nine generators feed a shared redesign barrier, followed by prediction and two analyses. It mixes backbone generation and sequence/structure co-design, erases the original sequence path from final Results, and does not equalize budgets or supported site annotations. A failed generator can block their common downstream analysis. Retire this all-on default; retain engines and an explicitly configured comparison option. [F5; E01, E03, E11]
6. **Some displayed Boltz-2 settings are ineffective.** Changing legacy `diffusion_samples`, `recycling_steps`, or `write_full_pae` does not change the submitted request in the tested path. The request uses `samples` and `num_recycles`, and the runner forces full PAE. This is a reproducibility and presentation defect, rather than evidence that the underlying predictor is invalid. [F6; E04, E06, E10]
7. **The antibody starter ends in a measurement view that no longer exposes prediction metrics directly.** Sequence QC deliberately replaces top-level metrics and retains previous measurements in nested lineage. A downstream filter on top-level `iptm`, for example, cannot see the preceding Boltz-2 score. Preserve measurement provenance, but make the final selection explicitly join prediction and QC results. [F7; E12]
8. **The active-learning screen is a real adaptive campaign, but learns docking-score labels.** Its probabilities do not establish experimental binding. A changing, warm-started random validation split and score-rank stopping rule cannot demonstrate whole-library recall or independent generalization. The implemented multi-computer docking and external-library streaming are valuable and should be retained. [F8; E17–E19]

Recommended direction: six task-oriented catalogue entries, with conditional method variants, plus two clearly experimental ligand-conditioned design entries. Preserve the custom editor and native tool workflows. Do not add obligatory extra predictors, MD, relaxation, or consensus scoring as a substitute for validation. Repair chain/sequence contracts and selection semantics before increasing computational sophistication.

### Scope and strength of evidence

- **Demonstrated** means established by the inspected code, an existing test, or a small CPU-only reproduction. It does not mean an experimental biological result was reproduced.
- **Scientific concern** means an inference about methodological risk, with the relevant assumptions stated. No new wet-lab success rate is inferred from a method paper or a working adapter.
- **Verified in code** rates the visible graph/adapter implementation, not installed model weights, container internals, or predictive performance. **Partially verified** is used where the assessment depends materially on native execution inside an external tool.
- Scientific sources were checked online where accessible. Current upstream documentation can differ from the revision packaged by Studio; it supplies method context, not proof of local implementation equivalence.
- Existing tests exercised temporary files and simulated execution. Across the two selected test runs: **310 passed; 10 optional QC container tests skipped**. No GPU benchmark or prospective design screen was performed. Reproductions additionally checked every named graph's topology and the specific handoff/settings defects above. See H.

## B. Complete pipeline inventory

### B1. What is actually predefined

`STARTER_TEMPLATES` is the public canvas/API registry. Both the UI and `/pipelines/templates` use it. `new_graph()` constructs the graphs, with several additional continuation or legacy factory IDs. The API permits `_generated` continuation through its helper only for Kejia and Peptide Trap; other continuations are exposed by tool result pages. [E01, E02]

**Notation:** T = target structure; C = binder-design context; MPNN = the shared `sequence.mpnn` adapter, whose default model is **ProteinMPNN**, even where the saved node kind is named `ligandmpnn`; R = Results. Arrows indicate dependencies, not automatic filtering. Numerical ranges below are source lines at the audited revision.

Ratings in the table are **scientific validity / practical research utility / computational efficiency / implementation confidence / improvement priority**. Utility reflects the shipped workflow, including material execution or interpretation defects. Efficiency is a static assessment, not a measured runtime comparison.

### B2. Public menu and API starters

| ID | Purpose and actual stages | Factory location in `pipeline_graph.py`, `new_graph()` | Recommendation | Ratings |
|---|---|---|---|---|
| `unidock_active_learning` | Prioritize a large compound library: centroid sample → UniDock Pro → cumulative scores → stopping decision → changing score labels → dense classifier → library prediction → exploit/explore selection → feedback; external final results | 675–709; defaults 175–259 | **IMPROVE** | Moderate / High / Acceptable / Verified in code / High |
| `foundry` | Protein binder backbones: T → C → RFD3 → shared MPNN (2 sequences/backbone) → RF3 → R | 988; 1009–1025 | **CONSOLIDATE** | Moderate / High / Acceptable / Verified in code / Medium |
| `kejia_idr` | Bind an unstructured sequence fragment: template threading/flexible peptide docking → MPNN (5 sequences) → AF2 (5 models, 3 recycles) → R | 1000–1001; 1063–1070 | **IMPROVE** | Moderate / High / Poor / Verified in code / High |
| `peptide_trap` | Retarget beta-strand-cleft peptide traps: matching/threading + Rosetta design → AF2 (1 model, 3 recycles) → Rosetta interface analysis → R | 998–999; 1059–1062 | **KEEP** | Moderate / High / Acceptable / Verified in code / Medium |
| `abmpnn_cdr_redesign` | Existing antibody structures → SAbR CDR annotation → AbMPNN (selected CDRs) → Boltz-2 → sequence QC → R | 710–737 | **IMPROVE** | Moderate / High / Acceptable / Verified in code / High |
| `pxdesign` | Protein binder design: T → C → PXDesign backbone generation → MPNN (1 sequence) → parallel AF2 and Protenix → R | 1004–1008; 1026–1058 | **CONSOLIDATE** | Moderate / High / Poor / Verified in code / Medium |
| `rfdiffusion` | Protein binder design: T → C → RFdiffusion generation → MPNN (2 sequences) → AF2 → R | 1004–1007; 1026–1032 | **CONSOLIDATE** | High / High / Acceptable / Verified in code / Medium |
| `genie3` | Protein binder design: T → C → Genie3 generation → MPNN (1 sequence) → AF2 (5 models, 20 recycles) → R | 1004–1007; 1026–1032 | **CONSOLIDATE** | Moderate / High / Poor / Verified in code / Medium |
| `protpardelle` | Default is unconditional fold generation, not binding: Protpardelle-1C → MPNN (8 sequences) → ESMFold binder/monomer alone → R | 1004–1005; 1028; 1033–1036 | **IMPROVE** | Moderate / Moderate / Acceptable / Verified in code / Medium |
| `esmfold2_binder` | Gradient-guided binder sequence generation → MPNN (2 sequences by shared default) → RF3 → R; target supplies the generator | 995–996; 1016–1018 | **CONSOLIDATE** | Moderate / Moderate / Acceptable / Verified in code / High |
| `all_generators_boltz2` | T → C → nine generators in parallel → shared MPNN (1 sequence each) → Boltz-2 → ipSAE → Rosetta interface analysis → R | 873–935 | **RETIRE** | Moderate / Moderate / Poor / Verified in code / High |
| `vhh_design` | T → C → BindCraft2 VHH, BoltzGen nanobody, ESMFold2 VHH → shared general MPNN → Boltz-2 → ipSAE → Rosetta → R | 936–987 | **IMPROVE** | Low / Low / Unknown / Verified in code / High |
| `rfd4_protein_binder` | T → C → RFD4 sequence/structure co-design (8 attempts) → Boltz-2 → ipSAE → Rosetta → R | 782–871 | **IMPROVE** | Low / Moderate / Unknown / Verified in code / High |
| `rfd4_small_molecule` | T → RFD4 ligand-conditioned protein-pocket design (8 attempts, 100–155 residues) → Boltz-2 → R | 782–871; options 811–820 | **IMPROVE** | Uncertain / Moderate / Unknown / Verified in code / High |
| `rfd4_enzyme` | T → RFD4 enzyme motif scaffolding (8 attempts, 150–200 residues) → Boltz-2 → R | 782–871; options 821–832 | **IMPROVE** | Uncertain / Moderate / Unknown / Verified in code / High |
| `cofold_protein` | One protein sequence → all 15 compatible registered prediction branches → separate Results per engine | 738–781; tasks 640–647 | **CONSOLIDATE** | High / High / Poor / Verified in code / Medium |
| `cofold_protein_protein` | Target + partner sequences → all 17 compatible prediction branches → separate Results; AF2 and BoltzGen Fold also require a target structure | 738–781 | **CONSOLIDATE** | Moderate / High / Poor / Verified in code / High |
| `cofold_protein_ligand` | Protein + ligand identity/collection → 9 compatible prediction branches → separate Results | 738–781 | **CONSOLIDATE** | Moderate / High / Poor / Verified in code / High |
| `cofold_protein_dna` | Protein + DNA sequences → 9 compatible prediction branches → separate Results | 738–781 | **CONSOLIDATE** | Moderate / Moderate / Poor / Verified in code / Medium |
| `cofold_protein_rna` | Protein + RNA sequences → 9 compatible prediction branches → separate Results | 738–781 | **CONSOLIDATE** | Moderate / Moderate / Poor / Verified in code / Medium |
| `cofold_protein_macrocyclic` | Protein + one cyclic protein partner → Boltz-2 only → R | 738–781; cyclic marking 745–746 | **CONSOLIDATE** | Moderate / Moderate / Acceptable / Verified in code / Medium |
| `blank` | Unconnected Target and Results placeholders; user builds a graph | 988; 1025 | **KEEP** | Uncertain / High / Unknown / Verified in code / Low |

The blank canvas is an editing facility, not a scientific workflow; its scientific validity cannot be assigned independently of the graph the user creates. Its initial validation failure is expected because it is deliberately unconnected. The VHH validation failure is different: it is a supplied connected scientific starter.

### B3. Additional factory and continuation graphs

These are actual explicit constructors, but they are **not additional public menu entries**. Non-menu factory availability alone is not evidence that users can discover them from a current tool page.

| ID | Exposure and actual stages | Factory lines; relevant entry point | Recommendation | Ratings |
|---|---|---|---|---|
| `bindcraft2` | Non-menu constructor: T → C → complete native BindCraft2 pipeline → R | 989–990; native adapter E14 | **KEEP** | Moderate / High / Acceptable / Partially verified / Medium |
| `boltzgen` | Non-menu constructor: T → C → complete native BoltzGen pipeline → R | 991–992; native adapter E13 | **KEEP** | Moderate / High / Acceptable / Partially verified / Medium |
| `complexa` | Non-menu constructor: Proteina-Complexa generation/co-design → MPNN (2 sequences) → RF3 → R; target is configured on generator | 993–994 | **CONSOLIDATE** | Moderate / Moderate / Acceptable / Partially verified / Medium |
| `rfdiffusion_beta_pairing` | Non-menu constructor: beta-strand-conditioned RFdiffusion → MPNN (10 sequences) → AF2 → R | 1002–1003; 1071–1083; node fields 401–424 | **IMPROVE** | Moderate / Moderate / Poor / Verified in code / Medium |
| `pxdesign_generated` | Tool-results continuation: saved backbones → MPNN → parallel AF2/Protenix → R | 1004–1008; 1026–1058; `generator_ui.py:365–378` | **CONSOLIDATE** | Moderate / High / Poor / Verified in code / Medium |
| `rfdiffusion_generated` | Tool-results continuation: saved backbones → MPNN (2) → AF2 → R | 1004–1007; 1026–1032; same UI entry | **CONSOLIDATE** | High / High / Acceptable / Verified in code / Medium |
| `genie3_generated` | Tool-results continuation: saved backbones → MPNN (1) → AF2 (5 models, 20 recycles) → R | 1004–1007; 1026–1032; same UI entry | **CONSOLIDATE** | Moderate / High / Poor / Verified in code / Medium |
| `protpardelle_generated` | Tool-results continuation: saved backbones → MPNN (8) → ESMFold alone → R; can originate from different generation tasks | 1004–1005; 1028; same UI entry | **IMPROVE** | Moderate / Moderate / Acceptable / Verified in code / Medium |
| `complexa_generated` | Tool-results continuation: saved co-designed structures/sequences → MPNN (2) → RF3 → R | 993–994; `complexa_ui.py:207–222` | **CONSOLIDATE** | Moderate / Moderate / Acceptable / Partially verified / Medium |
| `kejia_idr_generated` | Tool/API continuation: saved threaded backbones → MPNN (5) → AF2 (5 models) → R | 1000–1001; 1063–1070; `kejia_idr_ui.py:49` | **IMPROVE** | Moderate / High / Poor / Verified in code / High |
| `peptide_trap_generated` | Tool/API continuation: saved accepted designs → AF2 → Rosetta → R | 998–999; 1059–1062; `peptide_trap_ui.py:15` | **KEEP** | Moderate / High / Acceptable / Verified in code / Medium |
| `rfdiffusion_beta_pairing_generated` | Tool-results continuation: saved beta-pairing backbones → MPNN (10) → AF2 → R | 1002–1003; 1071–1083; `rfdiffusion_beta_pairing_ui.py:168` | **IMPROVE** | Moderate / Moderate / Poor / Verified in code / Medium |

### B4. Boundaries: tools, subworkflows, custom graphs, and campaigns

| Object | Actual role | Audit treatment |
|---|---|---|
| `NODE_CATALOG`, `WORKFLOWS`, prediction `ENGINE_PROFILES`, inverse-folding/quality-tool registries | Reusable nodes or adapters, often with native collections | Inspected as implementation dependencies; not counted as predefined graphs. A `prediction_*` node does not create another pipeline template. |
| BindCraft2/BoltzGen/Proteina-Complexa complete native tools; hidden PXDesign/Genie3/Protpardelle pipeline node kinds | A single canvas node can wrap a substantial scientific workflow | Assess the scientific contract without counting its internal stages as separate public templates. A generator-only node has a materially narrower contract. |
| `campaign.validation_sequence`, `campaigns.py:169–212` | Runs a caller-supplied list of lineage steps | Configurable internal orchestration, not a separate fixed preset catalogue. |
| `design_campaigns.py:3781` and associated refinement/evaluation functions | Multi-engine campaign with configurable sequence refinement, screening, and evaluation | Separate campaign machinery. Its optional filters and staged refinement do not automatically exist in canvas starters. |
| Active-learning docking | A predefined graph executed by a dedicated feedback controller | Included as one adaptive scientific starter; do not interpret its loop as a generic DAG cycle. |
| Saved user graphs/versions | User-created or edited pipelines | No claim of scientific validation merely because graph validation succeeds. No private saved projects were enumerated or edited. |
| Kejia scaffold libraries, Peptide Trap C37/C104.1, structure templates, generator examples | Molecular input templates/examples | Not additional computational pipeline definitions. |
| MD, docking, interaction-analysis adapters in the vendored ligand runtime | Standalone scientific capabilities | No public canvas MD starter was found. These capabilities do not establish an implemented docking → MD predefined pipeline. |

## C. Detailed scientific assessment

### C1. Cross-cutting scientific and implementation findings

#### F1. VHH execution and antibody sequence constraints — demonstrated, critical

`new_graph('vhh_design')` sets the BindCraft2 VHH branch's `binder_length` to `""`. `validate_graph()` applies its positive-length/range regular expression to every string-valued binder length, including this scaffold-mode value. A direct call fails with `Binder length must be a positive length or range, such as 74 or 60-100.` [E01:943–945; E04:1467–1471]

The shared downstream MPNN defaults are `model_type='protein_mpnn'`, blank designed/fixed residue selectors, `fix_interface_residues=False`, and `omit_aas='C'`. The general sequence-design function also defaults `auto_detect_cdrs=False`. The constructor does not insert SAbR/AbMPNN or a framework mask. BindCraft2's generator annotations preserve its pre-MPNN gate, but do not export the scaffold edit mask as a general fixed-residue contract; ordinary BoltzGen generation normalization likewise does not supply a CDR/framework mask. [E11, E13, E14]

**Consequence:** correcting only the empty length would leave unrestricted antibody-chain redesign, including the possibility of replacing native disulfide cysteines at mutable positions. That is inconsistent with presenting the whole graph as protected VHH design. This is a demonstrated missing constraint, not a claim that every resulting sequence necessarily loses its fold.

**Required improvement:** either preserve each generator's VHH sequence through prediction, or explicitly annotate and redesign selected CDRs while fixing framework/disulfide residues. Validate antibody identity and sequence/structure correspondence after redesign. Do not use a generic no-cysteine policy as an antibody quality rule. Keep VHH mode experimental until both the validation defect and sequence-constraint contract are verified.

#### F2. RFD4 component conversion preserves molecules but loses roles — demonstrated, critical

RFD4 normalization parses each generated complex into components and sets `raw_metadata.studio_molecule_id`. `_optimized_prediction_tasks()` then takes the molecular-library branch. It writes `protein_sequences`, DNA/RNA strings, and ligand CCD/SMILES identities, but does not write `target_protein_sequences`, `binder_sequence`, or `molecule_chain_options`. This branch also bypasses target-template construction in the ordinary candidate path. [E06:343–363; E09:108–122]

`_prepare_prediction()` assigns target/binder roles only when its target/partner input fields are populated. For generic `protein_sequences`, the role lists remain empty. Prediction-result candidates copy those empty lists and omit the designed sequence when no binder chain is identified. A reproduction starting with a two-chain RFD4-style candidate confirmed:

```text
input candidate: target_chains=['A'], binder_chains=['B'], binder_sequence='LMNPQRSTV'
submitted task: protein_sequences='ACDEFGHIK:LMNPQRSTV'
normalized task: target_protein_ids=[], binder_protein_ids=[], binder_sequence=''
molecule_chain_options={}
```

No prediction, MSA fetch, or GPU operation was needed for that reproduction. [E10:2060–2078, 2625–2649]

For a protein binder, this also removes the intended target-MSA/binder-single-sequence partition. Downstream ipSAE may fall back to A/B chain assumptions; Rosetta needs exactly one identified binder and one target, and multi-chain biological targets are already outside that scorer's supported contract. Merely attaching the original reference does not restore prediction-chain roles. The full runtime consequences depend on chain ordering and scorer behavior, but the missing metadata is proven. [E07, E15]

For pocket/enzyme workflows, ligand CCD identities are carried into prediction. However, `parse_structure()` emits non-polymer components only for residues with more than one atom and reconstructs them by residue code; single-atom ions are omitted, and custom chemistry, covalent links, atom-level conditions, and catalytic constraints are not represented in this handoff. Whether a particular input is harmed depends on its chemistry. [E09, E25]

**Required improvement:** maintain an explicit mapping from original component/chain identity to prediction chain identity, preserve designed protein sequence/roles, propagate appropriate per-chain MSA policy, and validate chemical identity and relevant bonds. Preserve constraints as reference data for post-prediction assessment; do not force the prediction to satisfy them and call that independent validation.

#### F3. Analysis is not a shortlist decision — demonstrated design gap

The `filter` node supports explicit pass evidence, a finite numeric threshold, direction-aware sorting, and a limit. None of the 21 public scientific starters instantiates it. `results` concatenates candidate collections, namespaces identities, and calls `_add_target_aligned_pose_rmsd()` where a reference exists. It does not rank or threshold the collection. [E01, E05]

Standalone ipSAE and Rosetta analysis retain candidates and report score/error counts. `_finish_metric_analysis()` marks success according to the presence of normalized rows, even if measurements failed for some or all rows. Thus a completed analysis job is not proof that all candidates possess the promised metrics. [E07:298–343]

This does **not** mean all upstream engines lack filters. BindCraft2 generation exports only trajectories passing its pre-MPNN gate, complete native tools have broader selection, and Peptide Trap checks `accepted is True`. Conversely, a pre-MPNN trajectory gate is not final binding evidence. [E13, E14, E16]

**Required improvement:** show measured, failed, missing, rejected, and selected states separately. Define a task-specific shortlist decision using valid interface/geometry measurements and diversity. Keep controls and rejected candidates accessible. Avoid an arbitrary universal weighted sum or universal ipTM cutoff.

#### F4. AF2 initial guess changes the meaning of validation — demonstrated configuration; scientific inference

Canvas adapters set `use_initial_guess=True` for AF2 candidates. The refolding runner therefore does not pass `--target-template-only`; the evaluation script instantiates the binder protocol with an initial guess and reads the designed binder chain from the input complex. `use_binder_template=False` and `use_interface_template=False` disable separate template features, not that initial guess. [E08]

This is a reasonable **conditioned self-consistency test** and may match a particular author's protocol. It is not fully independent evidence that the sequence discovers the proposed binding orientation. Low target-aligned pose RMSD can partly reflect conditioning. AF2 and related predictors also share structural training information; agreement should not be multiplied into a probability of experimental success.

**Improvement:** label prediction mode precisely. Offer target-only-template plus partner sequence as a distinct validation mode where appropriate, and reserve a second method for a shortlist or a concrete disagreement. Do not require a second predictor for every sample simply to increase the number of scores.

#### F5. All-on generation and sequence replacement — demonstrated default; benefit unproven

The miniprotein graph runs RFD3, RFdiffusion, Genie3, PXDesign, Protpardelle, BindCraft2 Generator, BoltzGen Generator, Proteina-Complexa Generator, and ESMFold2 Binder Generator. Some generate backbones; others already optimize or co-design sequence and structure. All feed one general ProteinMPNN step with one sequence per input, with no original-sequence bypass into prediction/Results. [E01:873–935]

That is defensible as a deliberate sequence-redesign experiment, but not a neutral comparison of generator performance. The design engines do not receive equivalent constraints: Protpardelle and ESMFold2 do not receive hotspot conditioning through these controls; coldspots are supported only on the BindCraft2/BoltzGen branches. `node_request()` deliberately removes unsupported annotations. This is honest adapter behavior, but prevents interpreting the pooled graph as a controlled epitope-matched benchmark. [E03]

The merged MPNN node waits for all nine parents. A failed parent blocks this common descendant, although unrelated branches can otherwise complete. Namespaced IDs retain branch provenance, but do not solve this failure barrier or make raw scores comparable. [E05:1006–1029]

**Recommendation:** retire the all-nine default. A user may still construct a comparison with selected methods, separate per-branch validation, documented budgets and conditioning, preserved original sequences, and a joint diversity-aware selection at the end.

#### F6. Ineffective legacy prediction controls — demonstrated

The Boltz-2 node exposes both `samples`/`num_recycles` and older `diffusion_samples`/`recycling_steps` fields. `_optimized_prediction_request()` forwards the former. The native batch options force `write_full_pae=True`. In the CPU request probe, changing only `diffusion_samples` to 5, `recycling_steps` to 50, and `write_full_pae` to false still submitted **one sample, three recycles, and the same engine options**. [E04, E06:435–471; E10:1705–1720]

Expose only effective settings or map legacy aliases unambiguously. Preserve full PAE when requested by ipSAE, and explain why that output is needed. Reject incompatible saved settings rather than silently implying they were applied. This does not establish an incorrect prediction; it establishes an incorrect user-to-request contract.

#### F7. Final antibody QC obscures prediction measurements — demonstrated integration gap

`run_quality()` intentionally writes only newly measured QC values into the new candidate's top-level `metrics`, changes its stage/tool to analysis, and nests the source record in lineage. The existing test explicitly asserts that inherited filter flags are absent. This protects measurement provenance and should not be indiscriminately reversed. [E12]

However, `abmpnn_cdr_redesign` sends only this sequence-QC output to Results. `metric_value()` reads current metrics/common normalization, not the nested source's confidence values. A following numeric filter on `iptm` would therefore reject these rows as missing that metric. The final collection is not a flat combined prediction-and-QC selection table.

**Improvement:** retain provenance-qualified prediction and QC measurements in a joined final report, or branch cheap sequence QC before prediction and select explicitly. Sequence composition/charge/GRAVY and the sequence “instability index” are descriptive heuristics, not measured solubility, melting temperature, or antibody developability. [E12; S1, S2]

#### F8. Active learning approximates a docking-score oracle — verified protocol; generalization unverified

The implemented oracle is UniDock Pro with lower-is-better Vina scoring. The controller rejects alternative score directions/fields it cannot honor; the UI's configurable labels do not mean arbitrary engines are implemented. Defaults are 10,000 initial centroid compounds, later batches of 50% of that sample, 80% exploitation and 20% random exploration, a 1024-bit fingerprint network with three 3000-unit hidden layers, and a label percentile tightening from 10% to 0.01% over ten rounds. [E17, E18]

The classifier learns whether a compound is in a **docked score tail**, not whether it binds experimentally. Class weighting and changing label definitions also prevent treating its sigmoid outputs as calibrated experimental probabilities. OpenVS provides precedent for active-learning screening, but its published RosettaVS results do not validate this different docking engine/preparation/scoring protocol. [S6]

The model step makes a stratified random 80/20 split each round, seeds it by round, and may reload previous weights. A current validation compound can have been in an earlier training set. This is a demonstrated non-independent validation design; the amount of optimistic bias is not measured here. Analogue/scaffold similarity can further inflate random-split performance. There is no independent scaffold/time-split test or whole-library hit-recall estimate in this path. [E18:1344–1438]

At the default batch progression, reaching round ten would yield roughly 55,000 docked observations before losses/deduplication. A 0.01% positive tail then contains only approximately six compounds, making a 20% positive holdout extremely small. The code correctly fails when either class has fewer than two members, but that is a technical minimum, not a sufficient validation standard.

Stopping requires stable top-N identities/order and score changes within tolerance over the configured comparisons, or library exhaustion/an iteration limit. This is better than ignoring identity changes, but it demonstrates local ranking stability, not exhaustion of useful chemical space. An unchanged ranking can reflect a poor surrogate, restricted chemotypes, or oracle bias. Default unlimited iterations also leave the total cost unbounded by an explicit user budget. [E19]

**Improvement:** name labels “top docking-score class”; retain a never-trained audit set, report enrichment/precision at the actual purchase budget and improvement per docked batch, monitor diversity and preparation failures, and let a budget bound the campaign. Keep random exploration. Compare with a matched random-screen baseline before increasing network size or adding uncertainty ensembles.

#### F9. Reuse is useful but does not hash every scientific dependency — demonstrated cache-key gap

Graph reuse keys hash node parameters, upstream fingerprints, selected source JSONL, the target file in some paths, and several adapters. This is a substantial strength. However, the cofold source path only adds a structure content digest for `target_pdb`; the contents of `sequence_file` and `compound_set_file` do not enter that source key. Generic RFD4 input/condition files and every transitive helper/model/image dependency are not comprehensively fingerprinted here either. Cache validation checks status and file existence, not a complete content/version manifest. [E20]

Changing a collection file in place while preserving its path/settings can therefore leave a downstream key unchanged. Immutable database snapshots mitigate this for managed snapshots, but user-supplied sequence files remain an exposed case. Random seed zero is resolved by some engines only after the graph cache key is formed: “reuse completed” is repeat-output reuse, not fresh stochastic sampling. [E10:1927–1930]

**Improvement:** fingerprint the actual referenced input contents and the resolved scientific execution specification, including effective seeds, checkpoint/image identity, and relevant adapter versions. Preserve a “fresh sampling” operation. No new workflow engine is required to solve this concrete cache/provenance need.

#### F10. Scientific validation gaps should be addressed selectively

| Item | What exists | What is missing from the standard explicit graphs; scientific consequence |
|---|---|---|
| Sequence/structure consistency | Chain parsing, candidate roles, selected sequences, PDB/mmCIF conversion, motif checks in specific adapters | No single uniform preflight covers assemblies, cofactors, target accessibility, chain identity, and task-specific input assumptions. Perform this once at input; do not repeat expensive calculations. |
| Region selection | Target chains, hotspots/coldspots, branch-specific context; unsupported features are flagged/removed | A requested epitope is not guaranteed merely by connecting context. Measure actual contacts and accessible surface in the biological assembly. |
| Chemical context | Explicit ligand/nucleic-acid task fields, RFD4 components, some native ligand-context preservation | Component identity conversion is narrower than an atom/bond graph; not every predictor preserves metals, covalent ligands, modifications or water-mediated geometry. Assess only context relevant to the hypothesis. |
| Structural quality | Pose comparison and available interface/quality tools | No universal final clash/chemistry gate is wired into all starters. Use cheap clashes/geometry before interpreting a high-confidence pose. |
| Diversity | Native tools may select diversity; clustering nodes exist | No sequence/interface clustering is used by public explicit design starters. A shortlist of near duplicates wastes experiments. Cluster/select once before final export. |
| Confidence/uncertainty | pLDDT, PAE, ipTM, ipSAE and native engine outputs | These overlap in information content. They are not independent binding probabilities; pLDDT can be high for an incorrectly positioned partner. |
| Controls | Custom graphs and separate target inputs are possible | No preset requires unrelated/near-homologue targets, scrambled peptide controls, or known docking controls. Add a small control set when it discriminates the intended biological hypothesis. |
| Selection/export | Normalized candidates, parent IDs, structures and tool-specific tables | Most starters lack a provenance-qualified, diversity-aware synthesis/purchase list with explicit rejection reasons and planned assays. Existing artifact export is a useful base; no new file format is inherently necessary. |

### C2. Individual pipeline reviews

#### 1. Active-learning docking (`unidock_active_learning`) — IMPROVE

**Objective and fit:** efficient prioritization of a very large purchasable compound library against a prepared receptor/pocket. It is not protein binder design, ligand generation, or experimental affinity prediction.

**Strengths:** an independent centroid seed source, compatible fingerprint checks, stable source-library/compound namespaces, cumulative scoring, random exploration, explicit preparation failures, bounded pose retention, and working multi-computer sharding. The model does not require copying the entire library into Studio's molecule database. [E17–E19]

**Required changes:** F8; validate the receptor/protonation/cofactors and box with a known ligand where available. The default hybrid mode needs a reference-ligand SDF; it is not a valid default for every receptor. Tautomer enumeration is skipped by default, so targets/chemotypes sensitive to tautomer state require an explicit preparation choice. One retained pose per compound limits alternative-pose review. Keep preparation chemistry and docking method constant across cumulative labels unless scores are separated by protocol.

**Output/decision:** cumulative ranked compound scores, stable IDs/SMILES, failure records, round/model history, selected/retained poses, optional bounded database import. Rank only comparable scores from the same protocol. Purchase a chemically diverse, plausible subset after pose/strain/contact inspection, solubility/availability checks, and counter-screen planning; docking kcal/mol is not a measured free energy.

**Efficiency:** batch independent ligands; overlap CPU preparation and GPU docking within measured resource limits; train/infer between rounds. Whole-library prediction and a second streaming exploration pass can become I/O dominated. Early chemistry filters should remove only clear inadmissible inputs, avoiding unvalidated medicinal-chemistry exclusion rules. Multi-workstation benefit is high for large rounds, low for tiny ones. Fewer training epochs or a smaller model should be benchmarked against enrichment before adding computation.

#### 2. Foundry binder design (`foundry`) — CONSOLIDATE

**Objective and fit:** a target-conditioned backbone, a sequence compatible with it, then a predicted complex. The separation is scientifically sensible; the saved `foundry_mpnn` kind now uses the same default ProteinMPNN contract as the general node. It is not inherently ligand-aware just because some backend code is named LigandMPNN. [E01, E11]

**Weaknesses:** no wired final interface/quality decision, no diversity selection, and RF3 defaults to five diffusion samples per candidate with ten recycles. Additional samples can reveal ambiguity, but their default cost is not justified by a Studio benchmark. Target geometry/templates reduce the independence of validation without making it useless.

**Output/selection:** sequences, predicted complexes and available RF3 confidence. Select using interface confidence, design-versus-predicted pose agreement, steric/accessibility checks and sequence diversity; never score on global pLDDT alone. Use structural templates only for the intended target state and include relevant assembly context.

**Efficiency/catalogue:** retain a Foundry method profile within targeted protein binder design. Parallelize backbones and their sequence/prediction jobs; defer extra samples or a different predictor to a shortlist. It does not need a separate task tile purely because it uses RFD3/RF3.

#### 3. Kejia IDR (`kejia_idr`) — IMPROVE

**Objective and fit:** recognize a particular intrinsically disordered peptide fragment by compatible template binding modes, threading, flexible docking, sequence design and conditioned AF2 evaluation. Treating the peptide as flexible distinguishes it from a rigid folded-target binder task.

**Strengths:** the adapter verifies the reported target fragment against the output peptide sequence and preserves fragment ranges/template provenance. Defaults retain five sequences per accepted backbone and five AF2 models, rather than silently presenting the upstream threading stage as complete design. [E16]

**Weaknesses:** a fragment bound to a template is a proposed bound state, not evidence that the same state is accessible or specific in the full protein. The local starter omits the upstream iterative diffusion/redesign and clustering branches. The authors explicitly caution that strand-pairing templates can have high AF2 pass rates without corresponding experimental success; that is particularly relevant to overinterpreting AF2 here. [S7]

**Selection/efficiency:** filter threading strain/repulsion and impossible geometries first; cluster registers/binding modes; run initial AF2 evaluation on representative sequences, then expand only promising families. Five sequences × five models means up to 25 model evaluations per threaded backbone before considering template/register multiplicity. Final decisions need peptide contacts, pose agreement, specificity against related fragments, and assay-context accessibility. CPU threading/docking and independent GPU evaluation can overlap across cohorts; no cross-host model inference is needed.

#### 4. Peptide Trap (`peptide_trap`) — KEEP within its restricted scope

**Objective and fit:** retarget published beta-strand-cleft scaffolds to compatible peptide sequences. This is a distinct scaffold-based recognition problem, not general de novo binder generation. Available scaffold identities and sequence/burial matching constrain applicability.

**Strengths:** the adapter checks exact target-fragment and designed-sequence correspondence, retains only native accepted designs, and records scaffold/protocol provenance. The pipeline does not add an unnecessary general MPNN step after Rosetta has already designed the interface. [E16]

**Limitations:** the local workflow is an adaptation, not reproduction of all scaffold-generation and selection stages in the paper. AF2 is conditioned as in F4. Rosetta design followed by Rosetta scoring gives correlated, model-dependent energetic evidence; the second scoring stage can diagnose packing or unsatisfied contacts but is not an independent binding assay. No claim about anti-aggregation activity follows from the predicted complex alone. [S8]

**Selection/efficiency:** preserve native acceptance, verify peptide register and buried polar geometry, inspect AF2 disagreement, and choose diverse scaffold/register representatives with related-peptide controls. Run Rosetta analysis on the useful shortlist if it is expensive and does not alter a decision for obvious failures. CPU register/design jobs and independent GPU AF2 jobs can be farmed separately. Retain the task; improve selection/export without expanding it into a generic multi-engine graph.

#### 5. Antibody CDR redesign (`abmpnn_cdr_redesign`) — IMPROVE

**Objective and fit:** generate sequence variants of selected CDRs on an existing antibody scaffold/complex. It must not be interpreted as validated affinity maturation. Backbone-conditional sequence likelihood does not measure affinity or preserve epitope recognition by itself.

**Strengths:** SAbR annotation, explicit loop selection, specialized AbMPNN weights, and framework/unselected-CDR fixing form a coherent contract. The implementation reuses upstream CDR annotations where applicable and includes fixed-motif identity safeguards for supported design sources. Input eligibility allows structure-bearing completed designs, rather than claiming all were final filter-passing binders. [E11, E24]

**Required improvements:** fix final evidence integration in F7. Run cheap sequence liability/composition checks before expensive prediction; preserve parental controls, framework/disulfide identity, and target/chain mapping. Evaluate whether redesigned loops retain the intended interface and whether the model predicts an alternative orientation. A single prediction is a screen, not a quantified probability of retained binding.

**Output/selection:** parent-linked CDR variants, exact mutable/fixed positions, predicted complexes, confidence/pose metrics and QC. Select diverse CDR sequences with plausible loop geometry and no obvious expression liabilities, retaining the original antibody as a control. Independent candidates are parallel; SAbR is CPU work. Multi-workstation use is beneficial for large variant libraries, usually unnecessary for a few variants.

#### 6. PXDesign (`pxdesign`) — CONSOLIDATE, preserve optional dual evaluation

**Objective and fit:** target-conditioned backbone generation → ProteinMPNN → AF2 and Protenix. The two predictor branches are genuinely parallel; this is not sequential refolding of one predictor's coordinates by the other. Upstream PXDesign distinguishes preview and extended evaluation modes. [E01; S9]

**Weaknesses:** the canvas default runs both predictors for every sequence and simply collects their rows, rather than producing a matched, per-sequence disagreement/selection table. The branches are complementary evidence, not statistically independent votes. No final threshold/diversity choice is specified, and F4 applies to AF2.

**Improvement:** preserve original identity across branches, compare compatible pose/interface metrics, and make the second predictor a shortlist or explicit benchmark option. Retain PXDesign as a method profile under the same biological task as the other protein binder backbones. Generation and MPNN are sequential per candidate; independent AF2/Protenix jobs can overlap if memory/resources permit.

#### 7. RFdiffusion (`rfdiffusion`) — CONSOLIDATE without retiring the method

**Objective and fit:** target-conditioned backbone generation → ProteinMPNN → AF2 self-consistency screening is an established design strategy. Published successes support the strategy, not this application's targets, thresholds or success rate. ProteinMPNN without obligatory FastRelax is a legitimate option; absence of relaxation alone is not a defect. [S1, S3]

**Weaknesses:** final experimental selection is unwired; AF2's actual initial-guess mode needs explicit labeling. A score threshold borrowed from another campaign can fail for different target sizes or interface classes. Target cropping may expose an artificial epitope, so check the full assembly and glycan/accessibility context when relevant.

**Output/efficiency:** two sequences per backbone and AF2 complexes/confidence. Filter egregious backbone clashes/site misses before sequence expansion; select by valid interface/pose criteria and diversity. Batch independent backbones and predictions. Consolidate the tile into the protein binder task while retaining RFdiffusion's conditioning controls and exact workflow mode.

#### 8. Genie3 (`genie3`) — CONSOLIDATE, reduce unproven default validation cost

**Objective and fit:** an alternative all-atom/backbone generation method followed by sequence design and complex evaluation. The author-inspired AF2 settings are explicit; the application uses its Initial Guess adapter rather than assuming equivalence with the authors' ColabFold template protocol. [E01, E08; S10]

**Weaknesses:** five AF2 models at twenty recycles for each sequence are a large mandatory evaluation budget, with no comparative evidence that they improve the final experiment allocation. Model/recycle counts are not interchangeable with independent seeds, and repeated highly related predictions are correlated evidence.

**Improvement:** use a measured staged evaluation budget, retain disagreement rather than only a best-model summary, and route a shortlist to additional effort only if it changes selection. Preserve Genie3 as a separate method profile, not a separate scientific task. Independent samples can run on different GPUs; extending one recycle loop across workstations is not justified.

#### 9. Protpardelle (`protpardelle`) — IMPROVE and rename around fold generation

**Objective and fit:** the default is unconditional protein generation followed by eight ProteinMPNN sequences per backbone and ESMFold monomer self-consistency. That is a coherent designability test, not evidence of target binding. This use is distinct from the target-conditioned protein-binder family and should remain distinguishable. The upstream project also uses sequence-design/self-consistency evaluation. [E01; S11]

**Weaknesses:** monomer confidence alone does not measure agreement with the generated fold. Make designed-versus-refolded backbone agreement and sequence validity explicit in the decision. When a user changes generation to binder mode, ESMFold-alone evaluation remains a monomer screen and cannot validate the interface. The `_generated` continuation does not distinguish these source objectives.

**Efficiency/output:** deduplicate similar sequences/backbones before predicting all eight variants, then allocate more sequence samples only to promising folds. Export new sequences, generation/reference structures, refolded structures, self-consistency metrics and diversity representatives. Optional complex prediction belongs only to target-conditioned tasks. Multi-workstation gains are candidate-level and become worthwhile at larger sample counts.

#### 10. ESMFold2 binder (`esmfold2_binder`) — CONSOLIDATE with a sequence-preserving branch

**Objective and fit:** gradient-guided target-binding sequence optimization, followed by a different predictor, can be useful. Current Biohub guidance describes a complete design protocol, but does not establish that adding the application's general MPNN step improves it. [S12]

**Weaknesses:** the generator has already optimized sequences using folding/language-model objectives. MPNN replaces these sequences, so the final RF3 results test redesigned variants rather than the original design protocol. Two variants per design are produced by the shared MPNN default. Epitope controls from Binder Design Context are not supplied as hotspot conditioning here; the single-target-chain restriction is enforced. [E03, E06, E11]

**Improvement:** preserve original sequences into the independent prediction path; make redesign an explicitly compared, constrained option. Evaluate original versus redesigned candidates on the same criteria. Retain the method under targeted binder design, with engine-specific assumptions and single-chain limitations visible. Gradient optimization itself is sequential; designs/seeds and downstream predictions are independent GPU jobs.

#### 11. All generators (`all_generators_boltz2`) — RETIRE as a standard preset

The graph is useful for adapter stress testing or a deliberate broad method experiment. It is poorly suited to default experimental candidate design because of F5: unequal conditioning, mixed scientific contracts, unnecessary sequence replacement, all-parent failure barriers, and no controlled comparison/selection. A single design attempt per method also provides little evidence about method-level performance.

The analysis tail does not turn the pool into a calibrated consensus. Boltz-2 already has default built-in ipSAE generation; the explicit ipSAE node can repeat a related computation. Check exact chain pairs/cutoffs before reusing either value: equality of settings, not a shared metric name, determines redundancy. [E10:151; 1776–1811]

Retirement here means removing the all-on recommended starting point, while preserving custom graphs, individual generators and an explicit comparison recipe. Multi-GPU execution shortens wall time but does not remedy scientific incomparability or reduce total compute.

#### 12. VHH design (`vhh_design`) — IMPROVE; experimental until repaired

VHH generation is a meaningful distinct input/design regime: framework identity, CDR geometry, disulfides and antigen accessibility matter. It should be a protected format mode within binder design, not an unrestricted miniprotein MPNN workflow. F1 is a blocking defect; assigning high current utility or measured efficiency would be misleading.

Use one supported generator by default. Retain generator framework/CDR masks, preserve original sequences, and invoke SAbR/AbMPNN only when loop redesign is requested. Evaluate antigen contacts, clashes, loop/sequence fidelity, and sequence diversity. Additional predictors or generators belong to explicitly budgeted exploration. Independent designs/antibody variants benefit from separate GPU jobs, not cross-workstation splitting of a single gradient trajectory.

#### 13. RFD4 protein binder (`rfd4_protein_binder`) — IMPROVE, then consolidate as a method variant

The intended sequence/structure co-design → different complex predictor → interface analysis is sensible. Omitting obligatory MPNN respects the co-design output. However, F2 currently breaks the role/sequence contract before analysis, so this graph should not be offered as an established, complete binder-validation path.

Repair and test exact original-to-predicted chain mapping, target-only MSA policy, binder sequence export and correct interface pairs. Then compare predicted pose/interface and physical plausibility; report RFD4's native scores separately from prediction confidence. The eight generated candidates can run independently. Rosetta should operate on a valid shortlist and only on supported chain groupings. RFD4 remains a method choice under the protein-binder task after these repairs; no extra prediction layer is intrinsically required.

#### 14. RFD4 small-molecule pocket (`rfd4_small_molecule`) — IMPROVE; experimental

**Objective:** design a protein around a specified ligand, not dock a compound library into an existing receptor. It must remain separate from screening and from protein–protein binder design.

**Concerns:** restraints are prefilled for the tutorial's 5ENI/5QA atom names and hydrogen-bond roles. Those defaults are not generally meaningful for another ligand; validate them against the selected molecule before generation. The downstream handoff reduces chemistry to component definitions and loses designed roles/sequence as in F2. General Boltz-2 structure confidence does not test whether ligand atom restraints or pocket specificity were retained, and affinity prediction is not explicitly requested by this graph's ordinary task payload.

**Required outputs/decision:** exact designed sequence and ligand identity, full complex/reference atom mapping, pocket geometry, ligand strain/clashes, satisfied/buried polar interactions, predicted pose deviations and diversity. Predicting a plausible complex is only hypothesis support. Optional competition/analogue controls should reflect the desired ligand selectivity. Do not add MD to all candidates; use a chemically validated shortlist if dynamics answers a specific question.

**Efficiency:** independent designs/predictions parallelize. Atom-selection checks and obvious clashes should precede expensive follow-up. Current end-to-end efficiency is unknown because no model runs were performed and the handoff needs repair. Upstream RFD4 task documentation supports rich conditioning, not the validity of every local tutorial-derived combination. [S13]

#### 15. RFD4 enzyme motif scaffolding (`rfd4_enzyme`) — IMPROVE; experimental

**Objective:** embed selected catalytic residues/atom geometry around a ligand in a scaffold. This is not interchangeable with generic ligand binding: catalysis depends on reaction chemistry, protonation, stereoelectronics and the relevant substrate/transition-state representation.

**Concerns:** prefilled motif/active-site atoms and chlorine restraints are tutorial-specific. Boltz-2 prediction does not measure catalytic turnover, an activation barrier, or transition-state stabilization. The graph has no explicit post-prediction catalytic-atom correspondence/tolerance assessment; F2 further narrows chemical context and removes role/sequence metadata. Scaffolding the wrong motif or preserving the motif in the generator but losing it in prediction could otherwise pass an irrelevant global-confidence screen.

**Improvement:** require a named catalytic hypothesis, validated motif atom selection and chemistry; check motif backbone/sidechain and ligand-relative geometry after both generation and prediction, including required cofactors. Export motif mapping, sequences, structures, deviations and limitations. Reserve mechanistic energy calculations for questions where their assumptions and cost are justified. Candidate-level parallelism is useful; no generic mandatory QM or MD stage is warranted. Keep experimental until its scientific output contract is explicit. [S13]

#### 16–21. Cofolding starters — CONSOLIDATE the interface, preserve typed scientific modes

All six constructors dynamically enumerate compatible node profiles and create a separate Results branch for each. They do not perform a cross-engine consensus/rank aggregation. Capability matching is a useful strength; the default of running every compatible implementation is the main redundancy. “Compatible” is a declared input contract, not demonstrated accuracy for a specific molecular class. [E01:650–655, 738–781; E10]

| Starter | Individual assessment and appropriate decision |
|---|---|
| `cofold_protein` | Monomer folding is a clear task. Fifteen branches include related model families and both standard/fast implementations; this duplicates work before any failure or uncertainty warrants comparison. Use one fit-for-purpose engine, domain/local confidence and geometry assessment, with optional alternate predictions for unresolved regions. No binding claim follows from a monomer model. |
| `cofold_protein_protein` | Predict a specified stoichiometric complex. Seventeen branches include exploratory ESMFold-v1 multichain use and target-template-required engines. Without the optional target PDB, the supplied full graph is not runnable under `require_inputs=True` because AF2 and BoltzGen Fold require it. This contradicts treating sequence-only entry as sufficient for every default branch. Check biological assembly, partner/target MSA policy and interface geometry. Exclude exploratory monomer-trained engine use from the standard complex default. |
| `cofold_protein_ligand` | Pose prediction for a protein and known chemical identity; a compound collection expands into independent protein–compound tasks. Nine engines are excessive as a default screen. Chai-1 accepts SMILES rather than CCD-only inputs, so choosing CCD while keeping every branch can fail task-format validation. Select compatible branches from the actual notation. Assess stereochemistry, atom identity, strain, contacts and physical validity before ranking. A high confidence score is not an affinity estimate. Boltz-2 affinity outputs require explicit task/property support and have a separate interpretation. [S4, S5] |
| `cofold_protein_dna` | Typed nucleic-acid complex mode is scientifically distinct from protein-only prediction. Require correct DNA strand count, sequence, intended duplex state and relevant stoichiometry. Inspect DNA geometry, protein–DNA clashes, interface contacts and uncertainty. Nine alternative engines are an optional comparison, not validation by majority vote. Preserve relevant ion/cofactor limitations. |
| `cofold_protein_rna` | RNA secondary/tertiary structure ambiguity and ion/state dependence need specific interpretation. Check sequence/strand identity, RNA geometry and local/interface uncertainty. A few seeds are not a thermodynamic ensemble. Nine engines need not be automatic; use further methods only for a defined ambiguity. |
| `cofold_protein_macrocyclic` | Only Boltz-2 is advertised as supporting native cyclic partner chains in this shared profile, and the task explicitly marks the partner cyclic. This is appropriately narrower than ordinary peptide input. Verify closure geometry and the intended head-to-tail chemistry; do not imply support for arbitrary staples, disulfide topologies, noncanonical residues or covalent receptor attachment from a single boolean. Retain as a typed mode, not a standalone “multi-engine comparison” claim. |

The 15 monomer branches are RF3, Boltz-2, AF3 JAX, AF3 PyTorch, OpenFold3 JAX, Protenix v2, ESMFold, ESMFold2, ESMFold2 Fast, Chai-1, AtlasFold, OpenFold3 OpenBind-0, OpenDDE, Protenix v1, and ColabFold. Protein–protein adds AF2 Initial Guess and BoltzGen Fold. The nine ligand/DNA/RNA branches are RF3, Boltz-2, AF3, AF3 PyTorch, OpenFold3 JAX, Protenix v2, Chai-1, OpenFold3 OpenBind-0, and Protenix v1. These are the factory's actual selections at the audited revision, not a recommendation to run them all.

AF3 official-weight acceptance defaults false and correctly needs the user's terms acceptance before execution. Its inclusion in an all-engine starter must be conditional on that existing requirement and model availability; this is a preflight/usability issue, not a reason to bypass the requirement. [E10:1923–1926]

All cofold modes benefit from independent inputs, seeds or selected engines running on separate GPUs when volume warrants it. A single small input rarely justifies coordinating several workstations. Outputs should preserve model-specific metrics and per-input failures, without comparing different engines' raw `ranking_score` values as if calibrated on one scale.

#### 22. Blank canvas — KEEP

This is necessary for flexible custom workflows. It should carry no implied scientific endorsement. Topological/type validation protects execution contracts; scientific assumptions remain the user's responsibility, supported by task-aware input checks and interpretable outputs. An unconnected blank graph is expected to fail runnable validation. [E04]

#### 23–26. Non-menu factories

- **`bindcraft2` — KEEP:** wraps a native design/refinement/validation/filtering procedure. It is sufficiently different from the generator-only MPNN graph to preserve as a native workflow option. The canvas maps the same requested count to final designs and maximum trajectories, so it can exhaust attempts before obtaining the requested number of final accepted designs; label these budgets separately. Avoid automatically appending MPNN/prediction again. Native scoring remains model-dependent, and local wet-lab performance was not tested. [E06:661–664; E14; S14]
- **`boltzgen` — KEEP:** preserves its native generation/inverse-folding/refolding/filtering and diversity-selection contract; this differs materially from generation → generic MPNN. Keep a native method route rather than replacing it merely to simplify the software graph. Inspect final native candidate pools separately from intermediates and retain exact full-chain sequences. The final selection budget is meaningful; adding another full evaluation stack is optional and should change a decision. [E13; S15]
- **`complexa` — CONSOLIDATE:** the graph uses generation-only Proteina-Complexa plus general MPNN/RF3, not the complete native Proteina-Complexa search/filter/evaluation procedure. Preserve original co-designed sequences and make redesign optional. Ligand/motif tasks have different chemistry; do not merge them into protein–protein binder defaults. Keep a method route, not a second generic task tile.
- **`rfdiffusion_beta_pairing` — IMPROVE:** explicit target/binder secondary-structure adjacency and strand-register conditioning are a distinct scientific assumption, not equivalent to an ordinary hotspot. Retain as a specialist binder-design variant with explicit target strand ranges. Ten MPNN sequences per backbone and initial-guess AF2 can amplify apparently good strand-pairing solutions. Check strand register, backbone hydrogen bonds, sheet-edge aggregation risk, full-target accessibility and specificity; increase sequence expansion only after geometry/diversity filtering. [E01:401–424, 1071–1083]

#### 27–34. Continuation graphs

Continuations are valuable handoffs and should not be removed simply to reduce a menu count. Their scientific method follows their parent profile, with these individual qualifications:

| ID | Assessment beyond the parent review |
|---|---|
| `pxdesign_generated` | Same dual-predictor cost/selection issue as PXDesign. Reuse generation artifacts; do not rerun generation. Preserve matched identity between both prediction branches. **CONSOLIDATE** into that task's “start from generated backbones” input mode. |
| `rfdiffusion_generated` | Useful backbone-to-sequence/prediction continuation. Confirm target/site metadata and fixed residues before redesign; retain the same conditioning caveat. **CONSOLIDATE** as an input mode. |
| `genie3_generated` | Same five-model/twenty-recycle default cost. A continuation is not evidence that upstream backbones passed final filters. **CONSOLIDATE** as an input mode with staged evaluation. |
| `protpardelle_generated` | Source may be unconditional, motif-conditioned or target-conditioned. A fixed ESMFold-alone tail is insufficient for binding/multichain objectives. **IMPROVE** by carrying source task identity and choosing the correct self-consistency versus complex-assessment path. |
| `complexa_generated` | It already starts with a sequence/structure co-design output. Preserve that sequence; do not require redesign before it can be assessed. **CONSOLIDATE** with the same task-aware candidate-input route. |
| `kejia_idr_generated` | Retains fragment/template information and avoids repeating CPU threading. Still needs mode diversity, proper final acceptance and the costly 5×5 evaluation review. **IMPROVE** with the parent method; keep the continuation. |
| `peptide_trap_generated` | Starts from accepted generated trap sequences, appropriately skipping unnecessary inverse folding. **KEEP** as a useful continuation, with the same limited experimental claims as the parent. |
| `rfdiffusion_beta_pairing_generated` | Preserve the strand/register-conditioning context and evaluate sheet geometry/specificity. Ten sequences per backbone should be an intentional budget. **IMPROVE**, preserving the specialist input mode. |

## D. Redundancy and consolidation

| Overlap | Recommended consolidation | What must remain distinct |
|---|---|---|
| Foundry, RFdiffusion, Genie3, PXDesign; applicable Complexa/ESMFold2/RFD4 binder variants | One targeted protein-binder task with method profiles and backbone versus co-design contracts | Backbone generators need sequence design; co-design outputs may be validated directly. Native full tools retain their own filtering/search. Beta-pairing requires explicit sheet conditioning. |
| Nine-generator miniprotein preset and three-generator VHH preset | Remove all-on defaults; choose one method and optionally compare selected alternatives | VHH framework/CDR/disulfide constraints must be enforced. “Use every available engine” is not a scientific objective. |
| Six cofold starters | One molecular-structure-prediction entry with typed task selection and capability-aware inputs | Monomer, protein complex, ligand, DNA, RNA, and cyclic peptide keep different validation rules, supported chemistry and output interpretations. This is shared presentation, not a universal interchangeable model. |
| `_generated` variants and their parents | Shared task setup with “new design” versus “continue generated candidates” inputs | Handoffs preserve source stage, references and constraints. Do not replace a co-designed sequence with an unrequested inverse-folding step. |
| Shared `foundry_mpnn` and `ligandmpnn` kinds | Present the already shared ProteinMPNN/model/backend configuration consistently | Saved IDs can remain for compatibility. ProteinMPNN versus ligand-aware or antibody-specific weights are scientific choices, not just backend aliases. |
| Built-in Boltz ipSAE and explicit ipSAE node | Reuse validated metrics when pair selection, cutoffs and implementation match | Different chain groupings/cutoffs are different measurements; avoid silently substituting or double-counting them. |
| General MPNN and native tools' internal inverse folding/refinement | Avoid appending redundant generic steps to a complete native workflow | Native methods may use specialized masks, objectives or selection not reproduced by generic MPNN. |

**Do not consolidate away these scientific distinctions:** ligand-pocket design versus library docking; enzyme motif scaffolding versus ligand binding; peptide-template recognition versus rigid folded-target design; antibody loop redesign versus unrestricted sequence generation; unconditional fold design versus interface design. Kejia and Peptide Trap can share a task entry and common export infrastructure, but their assumptions and algorithms should remain explicit, separate strategy branches.

No new orchestration layer, universal scoring service, or distributed inference architecture is justified by the evidence. The concrete needs are narrower: repair contracts, select task-aware defaults, expose final evidence coherently, and use the existing scheduler/collection mechanisms.

## E. Proposed optimized catalogue

These are proposed task entries and conditional graphs, not source changes made by this audit. Six standard entries cover the existing useful work; two experimental entries preserve specialized ligand-conditioned research. Existing named recipes can remain method presets inside those entries.

### P1. Targeted protein binder design — standard task; VHH mode gated until repaired

- **Objective/users:** design proteins that recognize a defined folded protein surface; structural biology, assay reagents and binder discovery.
- **Required inputs:** biologically relevant target structure/assembly, chain/region selection, intended epitope or explicit exploratory binding-site mode, binder-format and experimental budget.
- **Recommended stages:** validate target/context → one selected generator → cheap geometry/site checks → sequence design only for backbones or explicit redesign → one appropriate complex prediction → pose/interface/quality assessment → diversity-aware shortlist → export.
- **Optional:** selected alternative generator, a second predictor for shortlist disagreements, native complete BindCraft2/BoltzGen route, protected VHH mode, beta-pairing conditioning, limited specificity counter-targets. Preserve original co-designed sequences alongside requested variants.
- **Outputs/ranking:** exact sequences, original/predicted complexes, chain/motif mapping, interface confidence, valid pose RMSD, clashes/site-contact summaries, rejected/failed states, diversity clusters and selected experimental list. Rank within a comparable protocol, combining explicit quality gates and diversity rather than multiplying confidence scores.
- **Limitations:** computational interfaces do not establish affinity, selectivity, expression or function. Target templates/initial guesses reduce independence. Large biological context may be missing.
- **Parallelism/workstations:** independent backbones, sequences and predictions; CPU filtering can overlap GPU jobs. Multiple workstations help large campaigns after candidate sharding is made explicit; individual inference need not span machines.

### P2. Antibody CDR sequence redesign — standard after evidence/report repair

- **Objective/users:** explore selected loops of an existing antibody or VHH while preserving its framework; antibody engineers with a parental sequence/structure.
- **Required inputs:** parent complex or antibody structure, target/binder chain roles, exact mutable CDRs and fixed residues, parent control, desired number of variants.
- **Recommended stages:** verify antibody/numbering/structure identity → SAbR → AbMPNN with framework/disulfide/motif protection → cheap sequence checks and deduplication → complex prediction → interface/loop assessment → joined QC/prediction shortlist.
- **Optional:** additional loop sampling or predictor for ambiguous loop/pose cases, related-antigen counterprediction where scientifically motivated.
- **Outputs/ranking:** variants with parent/mutation masks, predicted structures, pose/interface metrics, sequence liabilities, diversity representatives and synthesis list. Prefer plausible retained interfaces and diverse valid variants; no sequence-model score is an affinity-improvement measurement.
- **Limitations:** fixed-backbone sequence design does not guarantee loop conformation or retained specificity. A designed antibody may express poorly despite confidence.
- **Parallelism/workstations:** CPU annotation/QC, independent GPU variant predictions. Multi-workstation benefit scales with variant count; avoid overhead for small libraries.

### P3. Template-based peptide/IDR recognition — standard specialist task with explicit strategies

- **Objective/users:** design recognition of a defined peptide fragment; peptide recognition and IDR research.
- **Required inputs:** target sequence and fragment boundaries, relevant biological/assay context; select a compatible Peptide Trap scaffold strategy or Kejia flexible-template strategy. These are not interchangeable defaults.
- **Recommended stages:** validate fragment/template compatibility → strategy-specific matching/threading and design/flexible docking → early geometry/strain/register filtering → MPNN only when required → conditioned self-consistency evaluation → mode/sequence diversity selection.
- **Optional:** expanded AF2 model budget for an unresolved subset; related-sequence controls; upstream iterative refinement only as an explicitly documented extension.
- **Outputs/ranking:** full requested and modeled fragment sequences/ranges, scaffold/register provenance, designed sequences/structures, interface/pose/strain metrics and representative shortlist. Rank binding-mode plausibility and specificity evidence, retaining diverse registers.
- **Limitations:** high-confidence beta pairing and a bound peptide conformation do not establish accessibility or affinity in the native disordered protein. Peptide-trap predictions do not establish anti-aggregation activity.
- **Parallelism/workstations:** CPU template/register tasks, independent GPU predictions; overlap across cohorts. Workstation farming helps large template/register searches. Keep Kejia's current workflow experimental within this task until its selection and evaluation budget are addressed.

### P4. Molecular structure prediction — one standard entry with typed modes

- **Objective/users:** produce structural hypotheses for specified molecular systems; structural biologists and modelers, rather than a universal binding-screen claim.
- **Required inputs:** sequence/stoichiometry plus typed molecule identities; optional relevant templates/MSAs; explicit cyclic/modified/bond chemistry where supported.
- **Recommended stages:** task-specific validation/capability/availability checks → one suitable engine → molecular geometry/confidence review → identity-preserving export. Keep separate monomer, PPI, protein–ligand, DNA, RNA and cyclic-partner modes.
- **Optional:** further seeds or selected alternative models where uncertainty changes the downstream decision; ligand affinity head only within its supported small-molecule scope and with its own documented units/training limitations.
- **Outputs/ranking:** per-input/sample structures, model-specific confidence, relevant pair/chain/atom mapping, validity checks, failure states and comparison views. Compare equivalent geometric conclusions; do not sort pooled raw engine ranking scores on one scale.
- **Limitations:** static models are not solution ensembles or catalytic/affinity measurements. Complex confidence can be high for an incorrect pose; unusual chemistry and ions may be unsupported.
- **Parallelism/workstations:** independent inputs/seeds/models. Shard collections at the input level where supported. A local two-GPU mode is an engine-specific memory/performance option, not proof of cross-workstation inference support.

### P5. Small-molecule hit prioritization — standard task; optional adaptive execution

- **Objective/users:** prioritize purchasable compounds for binding/functional assays against a defined receptor site; screening scientists.
- **Required inputs:** chemically prepared receptor and pocket/box, reference ligand when required, compound identities/stereochemistry, preparation protocol and purchasing/compute budget.
- **Recommended stages:** validate receptor/known-ligand behavior and inputs → one preparation/docking protocol → pose/chemistry/contact checks on ranked compounds → diversity/availability-aware purchase shortlist. For ultralarge libraries, enable the existing active-learning controller with exploration, a held-out audit set and bounded stopping.
- **Optional:** a justified second pose/scoring method on a shortlist; matched random baseline; selected-state ensemble when receptor-state ambiguity is relevant. Use the existing direct docking tool route for modest libraries; a separate new preset is not required merely to wrap a single tool.
- **Outputs/ranking:** stable compound IDs, parent/state mapping, protocol-specific docking scores, poses, preparation failures, round/enrichment history, chemical diversity and a purchase/assay table. Experimental affinity is not inferred from score units.
- **Limitations:** docking/surrogate errors, receptor rigidity/state, protonation/tautomer coverage, chemical applicability domain, and assay artifacts remain. Counter-screens and binding confirmation are part of experimental planning.
- **Parallelism/workstations:** strongest existing multi-workstation case: independent ligands/shards with local preparation/GPU docking, compact score return and selective pose transfer. Training/selection are sequential between rounds; library inference may be streamed/sharded if measured as a bottleneck.

### P6. New protein folds and self-consistency — standard exploratory design task

- **Objective/users:** generate new folds and sequences likely to realize them; protein-design method development and scaffold exploration. This is the existing unconditional Protpardelle use, renamed by task.
- **Required inputs:** length/topology or explicit unconditional mode, optional fixed motif with validated mapping, sample/sequence budget.
- **Recommended stages:** generation → geometry/fixed-motif checks → sequence design as required → sequence deduplication → monomer/refold self-consistency → structure/sequence diversity selection.
- **Optional:** expanded sequence sampling for promising folds; a second predictor for uncertain folds; target-complex evaluation only when the scientific objective actually includes a binding partner.
- **Outputs/ranking:** exact sequence, generated/refolded structures, backbone agreement, local confidence/geometry, motif retention and diverse representatives. High monomer confidence alone is insufficient.
- **Limitations:** self-consistency is not thermodynamic stability, solubility, function or novelty relative to all known proteins. Generalization of an upstream model is not benchmarked by adapter tests.
- **Parallelism/workstations:** independent structures/sequences, lightweight CPU checks between GPU stages. Workstations help large generation cohorts; no mandatory distributed single inference.

### P7. Ligand-conditioned protein-pocket design — experimental

- **Objective/users:** design a protein pocket for a specified small molecule; specialists in ligand-binding protein design.
- **Required inputs:** exact ligand chemistry/state and atom mapping, relevant cofactors, desired interaction constraints and design region/length. Tutorial-specific atom names must not silently transfer.
- **Recommended stages:** chemical/constraint validation → RFD4 co-design → pocket/motif geometry check → independent complex prediction with preserved component/role identity → atom-level ligand/pocket comparison → diverse shortlist.
- **Optional:** constrained redesign with an appropriate ligand-aware model; selected analogue controls; focused relaxation/dynamics after a chemically valid shortlist.
- **Outputs/ranking:** designed sequences and complexes, ligand correspondence, strain/clashes/polar interactions, predicted pose agreement, diversity and experiment-ready constructs. Rank on pocket plausibility and task constraints, not global pLDDT or an unrequested affinity score.
- **Limitations:** binding/affinity/selectivity remain unproven; unusual chemistry/covalent context may not survive the current handoff. Repair F2 before presenting as an end-to-end validated path.
- **Parallelism/workstations:** independent designs and predictions; CPU chemistry analysis between stages. Multi-workstation gain is moderate/high at scale, contingent on explicit data/atom mapping.

### P8. Catalytic motif scaffolding — experimental, separate from pocket design

- **Objective/users:** position a specified catalytic motif in a compatible scaffold; enzyme-design specialists with a defined reaction hypothesis.
- **Required inputs:** reaction/substrate or appropriate mechanistic representation, catalytic residue/atom definitions, protonation/cofactors and acceptable geometry tolerances.
- **Recommended stages:** verify the catalytic input specification → motif scaffolding/co-design → fixed-residue and ligand-relative atom-geometry checks → independently predicted structure → repeat catalytic geometry assessment → shortlist/export with mechanistic limitations.
- **Optional:** justified sequence redesign, reaction-specific calculations or focused simulations on a few candidates; catalytic-dead motif controls. These are not universal mandatory stages.
- **Outputs/ranking:** constructs, residue/atom correspondences, exact motif retention/deviations, pocket geometry, uncertainty and a catalytic-assay plan. No claim of activity from fold confidence.
- **Limitations:** good scaffold geometry does not establish transition-state stabilization, turnover or selectivity. No automated generic quantum-chemistry stage is justified without reaction-specific evidence.
- **Parallelism/workstations:** independent scaffold/prediction jobs; optional mechanistic calculations parallelize over selected candidates. Do not distribute one inference by default.

## F. Computational and multi-workstation assessment

### F1. Per-family execution decisions

| Existing pipeline(s) | Expensive stage(s); avoidable work | Correct unit of parallelism and required sequencing | Multi-workstation assessment |
|---|---|---|---|
| Foundry | RFD3 generation; RF3's five samples per sequence; no early shortlist | Independent backbones/sequences; generation → sequence design → prediction per candidate | Useful for larger campaigns; generic graph currently submits collection jobs, not automatic cross-host candidate shards |
| RFdiffusion and continuation | Generation plus AF2; expansion of obvious site/clash failures | Backbones and AF2 inputs independent; local sequence design can batch | Useful at scale, especially GPU predictions; transfer normalized candidates/references once |
| Genie3 and continuation | Five AF2 models × 20 recycles per sequence | Backbones/sequence jobs; a model's recycle loop is sequential | High potential throughput benefit, but first reduce unjustified evaluation cost |
| PXDesign and continuation | AF2 + Protenix for every sequence | Parallel predictor branches from the same sequence collection | Multiple GPUs/workstations can compare shortlisted candidates; preserve matched identity and settings |
| Protpardelle and continuation | Eight sequence refolds per backbone | Independent sequence/structure pairs after generation | Useful for large cohorts; avoid predicting duplicates or calling monomer confidence interface validation |
| ESMFold2 binder | Sequential gradient loop, followed by redesign and RF3 | Independent design trajectories; predictions after each completed sequence | Separate GPUs useful; splitting one gradient loop across machines has no demonstrated benefit here |
| Kejia and continuation | Template/register flexible docking plus 5 sequences × 5 AF2 models | CPU registers/templates independent; GPU evaluation after design | Strong CPU/GPU overlap opportunity across cohorts; mode filtering should reduce expansion first |
| Peptide Trap and continuation | Rosetta design/relaxation, AF2, repeated Rosetta scoring | Independent scaffold/register jobs; scoring after valid prediction | Useful for larger register sets; compare cost of remote transfer with CPU design/scoring time |
| Antibody CDR redesign | Boltz prediction for every variant; late sequence QC | CPU annotation/QC first, then independent variants | Useful only beyond small libraries; preserve parent/framework constraints across shards |
| All generators | Nine generation stacks + pooled redesign + prediction + analysis | Generators independent but shared MPNN is a barrier | More hardware masks wall time, not waste; retire default before optimizing its distribution |
| VHH | Current validation failure; later three generators plus unrestricted redesign | First repair; then independent protected-format designs | Unknown current efficiency; no reason to distribute an invalid graph |
| RFD4 protein binder | Co-design + Boltz + two analysis stages | Independent designs, then correct role-preserving prediction | Potentially useful at scale after F2; CPU analysis of one cohort can overlap next GPU cohort |
| RFD4 pocket/enzyme | Co-design and complex prediction; optional expensive mechanistic work | Independent designs, validated atom mapping before follow-up | Useful for design cohorts; transfer component/constraint manifests with structures |
| Six cofold starters | 15/17/9/9/9/1 branches; duplicate implementation families | Inputs, seeds or selected engines independent | One primary model usually sufficient; sharding many inputs is more useful than running every engine once |
| Native BindCraft2/BoltzGen factories | Native iterative search/refinement/selection | Independent trajectories/design batches, subject to native controller state | Respect native pipeline batching/checkpoints; do not split one stateful search arbitrarily |
| Complexa factory/continuation | Co-design, redundant mandatory redesign, RF3 | Independent designs; preserve original sequence branch | Use native candidate-level batching; compare original vs redesigned only when requested |
| Beta-pairing factory/continuation | Ten sequence evaluations per backbone | Registers/backbones and sequence predictions independent | Good large-search potential after early strand/geometry filtering |
| Active-learning docking | Repeated docking, full-library inference, preparation and pose I/O | Ligand shards within a round; score accumulation → training → selection between rounds | **Already explicitly implemented**, with weighted assignments, worker batches and score/pose transfer; strongest supported case |

### F2. What the current architecture actually supports

The generic runner uses a CPU coordinator and scheduler child jobs. It can enqueue independent ready branches without waiting for each branch's computation to finish. Candidate collections are commonly passed to one native batch job. This is distinct from splitting a collection into independently scheduled cross-host shards. No automatic candidate-to-workstation sharding is established by `graph_pipeline._submit()`; the active-learning controller supplies that behavior explicitly for docking. [E05, E21]

Fan-in currently creates an all-parent barrier. Move scientifically independent per-generator filtering/prediction before final comparison where appropriate, or define explicit partial-pool behavior. This recommendation addresses a demonstrated failure coupling; it does not require replacing the scheduler.

Some prediction profiles advertise local two-GPU “big” modes. They concern a particular runner's inference/memory behavior, not distributed inference across ordinary networked workstations. Candidate-level jobs offer the clearest fault isolation and scientific provenance. There is no demonstrated need to introduce multi-host tensor parallelism for these presets. [E10, E21]

### F3. Data, storage and overlap

- **Protein prediction:** full PAE is quadratic in residue/token count. Keep matrices required for interface assessment and provenance, but avoid making unnecessary copies across every analysis node. Five models/samples multiply structures and sidecars; nested candidate lineage can also grow JSONL size. Preserve parent references and validated measurements without losing recoverability.
- **Large screening:** packed 1024-bit fingerprints require at least 128 bytes per compound before IDs/SMILES/container overhead. One billion fingerprints alone is approximately 128 GB decimal. Unpacking all docked features to float32 requires 4096 bytes per training compound before tensors/labels/model state. Streaming inference is essential; the current training preparation still materializes the cumulative training feature matrix. [E18]
- **Transfers:** the active-learning controller returns compact scores and selectively collects top poses with verification, then can remove bulk worker artifacts. This is scientifically and operationally preferable to transferring every pose when only a shortlist is reviewed. Preserve enough preparation/protocol/source identity to regenerate a discarded pose and audit failure bias. [E19]
- **CPU/GPU overlap:** prioritize CPU chemical preparation, annotation and cheap geometry/QC while separate GPU candidate jobs run. Do not occupy a GPU allocation with a long CPU-only analysis tail if the scheduler can run it separately; `_submit()` already distinguishes CPU analysis and CPU peptide/template tools. [E21]
- **Multi-workstation costs:** environment/checkpoint consistency, path portability, model availability, database snapshot identity, network bandwidth, job granularity and retries can outweigh gains for small jobs. Pin effective settings and input hashes per shard. A local batch is often more efficient than many tiny network jobs.

### F4. MD and energetic refinement

No MD stage occurs in the actual public starter graphs. Bound-ligand MD and related GROMACS/OpenMM capabilities exist as standalone adapters. Their availability is not evidence that a scientifically validated docking → MD pipeline is supplied. [E23]

Do not automatically attach short MD to every predicted binder or docked compound. A stable-looking trajectory started in a proposed complex is not evidence of binding affinity, catalytic activity, or pose uniqueness. Conventional MD must first have validated chemistry/protonation, force-field/cofactor treatment, equilibration and a sampling question. Independent replicas and uncertainty analysis are appropriate when dynamics is informative; trajectory frames are correlated observations. [S16]

Likewise, Rosetta interface energy is a model-dependent structural score, and MM/GBSA is not automatically a converged experimental free energy. Optional simulations or refinement belong on a small, chemically plausible shortlist when they can discriminate alternatives. Parallelize replicas/candidates; do not assume an ordinary trajectory must span workstations. The audit does not endorse or benchmark the separate MD implementation.

## G. Prioritized improvement roadmap

### G1. Critical scientific correctness and execution contracts

| Priority | Action | Evidence and completion criterion |
|---|---|---|
| P0 | Repair RFD4 component/role/sequence handoff before promoting the three RFD4 starters | F2. CPU contract tests must preserve exact target/binder identities, designed sequence, per-chain MSA policy and ligand identity through normalized prediction output; test multi-chain targets and nontrivial chain order. Validate supported cofactors/bonds explicitly. |
| P0 | Repair VHH runnable validation and enforce framework/disulfide/CDR protection | F1. The unedited scaffold-mode starter must validate after real inputs are supplied; redesigned framework residues and fixed motif identities must remain unchanged in test outputs. No unrestricted cysteine exclusion on a scaffold requiring disulfides. |
| P0 | Distinguish successful computation, valid measurement and selected candidate | F3. Results must show per-item missing/error states; analyses with zero scored candidates must not appear as scientifically complete. Selection must not treat unknown metrics as pass evidence. |
| P1 | State AF2 conditioning accurately and expose a distinct target-only validation choice | F4. Saved/effective settings and result labels must agree; an initial-guess run cannot be labeled an unconditioned independent pose prediction. |
| P1 | Make displayed Boltz-2 controls match effective requests | F6. Alias/default tests must demonstrate that changing a shown control changes the execution specification or is rejected/explained. |
| P1 | Produce a coherent final antibody evidence table without falsely relabeling lineage measurements | F7. A final selection must access valid preceding prediction metrics and QC with their original provenance, and preserve the parental comparator. |

### G2. High-value scientific usability improvements

1. **Define task-specific acceptance and experimental outputs.** Add explicit shortlist criteria for interface/fold/pocket/motif objectives, finite metric checks, diversity and clear failure reasons. Export exact synthesis sequences or purchasable compound IDs, structures, mutable/fixed positions, relevant chemistry, scores with units/directions, and planned controls. Do not use one universal weighted confidence score.
2. **Replace engine-combination tiles with the task catalogue in E.** Preserve native end-to-end workflows, specialist assumptions and continuation entry modes. Retire only the all-on generator default; keep custom comparisons possible.
3. **Add focused input preflight.** Validate the biological assembly, selected chains/regions and relevant ligands/cofactors once. Resolve template-required/SMILES-only/terms/model-availability requirements before launching any cofold branch. Tutorial atom selections must match actual chemistry before a costly job.
4. **Validate screening as a selection process.** Use a never-trained audit set or carefully documented holdout protocol, a matched random baseline, precision/enrichment at practical acquisition budgets, chemical diversity and preparation-failure accounting. Label surrogate probabilities by their docking-score target. F8.
5. **Complete reuse provenance.** Hash mutable collection/condition files and resolved scientific dependencies; expose fresh sampling separately from reuse. F9. Keep existing immutable run snapshots and parent links.

### G3. Computational efficiency improvements

1. Move cheap invalid-input, geometry, sequence and duplicate checks before prediction; apply costly interface/energetic analysis to a useful shortlist. Do not filter out novelty using unjustified heuristics.
2. Replace all-engine prediction defaults with a primary method and an optional escalation path. Measure whether extra samples/models/recycles improve experimental selection before preserving the Genie3/Kejia/beta-pairing expansion budgets.
3. Make MPNN conditional on scientific need; preserve original co-designed sequences as a direct evaluation route. Compare redesign as an experiment rather than silently replacing the parent design.
4. Reuse matching ipSAE measurements instead of recomputing them; preserve pair/cutoff/version provenance. Avoid repeated Rosetta analysis when the same structural question was already answered.
5. For large protein cohorts, add explicit candidate-level sharding using existing scheduler primitives only after profiling collection bottlenecks. Keep active-learning docking's current batching/retention strengths. Resolve fan-in failure behavior before promising broad multi-workstation speedups.
6. Monitor prediction/PAE storage, training RAM and external library scans. Tune batch sizes from measured throughput and memory, not advertised engine capability alone.

### G4. Optional future enhancements

- Prospectively benchmark the recommended task defaults on a small, representative target set with actual experimental outcomes, fixed budgets, novelty controls and recorded failure rates. Adapter tests cannot substitute for this evidence.
- Add selected near-homologue/related-peptide or chemical-analogue counterchecks where they test specificity. Broad proteome-wide counterprediction is not required for every early-stage design.
- Offer uncertainty-aware acquisition only if it improves screening enrichment per unit cost over the existing exploitation/random baseline. Multiple surrogate ensembles are not automatically beneficial.
- Add focused dynamics or reaction-specific calculations only for shortlisted questions that require them; do not create a universal MD certification step.
- Consider a reproducible comparison workspace with matched budgets and task/conditioning declarations. This is a custom research recipe, not justification for restoring every engine as a default.

## H. Evidence, verification and references

### H1. Repository evidence index

Line numbers refer to the source revision at the top of this report. Paths below are relative to the repository root; links open the audited local files at a representative line. References elsewhere cite these evidence IDs and, where useful, additional ranges.

| ID | Source and functions/lines | Relevance |
|---|---|---|
| E01 | [pipeline_graph.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/pipeline_graph.py:675), `new_graph`: 675–1090; `STARTER_TEMPLATES`: 1586–1608; `TEMPLATE_REFERENCES`: 617–638 | Authoritative inventory, stage graphs, counts and template defaults |
| E02 | [workflow_setup_api.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/workflow_setup_api.py:95): 95–127; [pipeline_graph_ui.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/pipeline_graph_ui.py:489): 489–490; [generator_ui.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/generator_ui.py:365): 365–378 | Menu/API exposure and continuation distinction |
| E03 | [pipeline_graph.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/pipeline_graph.py:425): 425–443; [graph_pipeline.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/engine/workflow_engine/workflows/graph_pipeline.py:611), `node_request`: 611–675 | Hotspot/coldspot support, inherited target context, single-chain ESMFold2 restriction |
| E04 | [pipeline_graph.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/pipeline_graph.py:1320), `validate_graph`: 1320–1509; node defaults 129–167, 327–340; `_prediction_task_support`: 551–569 | Type/cycle/input checks; VHH length rejection; prediction settings and capability constraints |
| E05 | [graph_pipeline.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/engine/workflow_engine/workflows/graph_pipeline.py:960), `run_graph_pipeline`: 960–1285; `filter_candidates`: 475–492; `_add_target_aligned_pose_rmsd`: 839–878 | Scheduler, failure barriers, fan-in, Results collection, RMSD and filtering |
| E06 | [graph_pipeline.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/engine/workflow_engine/workflows/graph_pipeline.py:190), `_optimized_prediction_tasks`: 190–427; `_optimized_prediction_request`: 435–471; `node_request`: 495–759 | Candidate → prediction data contract; effective settings; RFD4 component conversion; forced AF2 initial guess |
| E07 | [analysis.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/engine/workflow_engine/workflows/analysis.py:252), `_analysis_output_candidates`: 252–280; `_finish_metric_analysis`: 298–343; `run_ipsae_analysis_contract`: 1057–1103; `run_rosetta_metrics_contract`: 1106–1246 | Measurement versus acceptance, missing PAE, per-item errors and Rosetta chain restrictions |
| E08 | [refolding.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/engine/workflow_engine/workflows/refolding.py:4357), `run_af2_initial_guess_complex_refolding`: 4357–4400, 4526–4552; [af2_initial_guess_binder_eval.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/engine/workflow_engine/tools/af2_initial_guess_binder_eval.py:470): 470–496, 537–558 | Exact meaning of initial guess versus target-only template |
| E09 | [rfd4_proteina.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/engine/workflow_engine/workflows/rfd4_proteina.py:69), `normalize_outputs`: 69–123 | RFD4 chain roles, sequence, components and upstream provenance |
| E10 | [optimized_prediction.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/engine/workflow_engine/workflows/optimized_prediction.py:76), `ENGINE_PROFILES`: 76 onward; `_batch_options`: 1647–1740; `_prepare_prediction`: 1885–2218; `run_optimized_prediction`: 2254 onward, especially 2625–2678 | Effective model settings, random seeds, terms, input roles, batch execution, confidence and affinity output |
| E11 | [mpnn_input.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/mpnn_input.py:28), `defaults`: 28–36, `prepare`: 102–185; [mpnn.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/engine/workflow_engine/workflows/mpnn.py:158), `run_mpnn_sequence_design`: 158–291; [sequence_design.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/engine/workflow_engine/workflows/sequence_design.py:352): 352–363, 636–700 | Default model, mutable/fixed chains/residues, CDR protection, original sequence replacement and motif checks |
| E12 | [quality.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/engine/workflow_engine/workflows/quality.py:145), `run_quality`: 145–174; [quality_runner.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/engine/workflow_engine/workflows/quality_runner.py:28), `sequence_metrics`: 28–40; [output_vocabulary.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/output_vocabulary.py:310), `metric_value`: 310–314; [test_quality_tools.py](/home/user/programs/git-projects/mn-molecular-design-studio/tests/test_quality_tools.py:59): 59–93 | Intentional QC measurement isolation, lineage preservation and final filtering consequence |
| E13 | [design.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/engine/workflow_engine/workflows/design.py:618), `_normalize_boltzgen_generation_candidates`: 618–697; `_normalize_boltzgen_vanilla_candidates`: 936–1010; `run_boltzgen`: 3600 onward | Generator/native distinction; final full-chain identity, native pools and ligand/task controls |
| E14 | [bindcraft2.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/engine/workflow_engine/workflows/bindcraft2.py:842), `_normalize_bindcraft2_trajectory_candidates`: 842–979; `build_bindcraft2_settings`: 199–403; source identity 44–45 | Native workflow versus pre-MPNN gate; scaffold modes and masks |
| E15 | [analysis.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/engine/workflow_engine/workflows/analysis.py:693), `_chain_role_groups`: 693–727; Rosetta validation 1173–1182; [design_validation.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/design_validation.py:79), `candidate_passes_filters`: 79–106 and `candidate_description`: 109–131 | Explicit pass evidence and chain role semantics; fallback risks |
| E16 | [kejia_idr.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/engine/workflow_engine/workflows/kejia_idr.py:44), `normalize_output`: 44–69; [peptide_trap.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/engine/workflow_engine/workflows/peptide_trap.py:29), `normalize_output`: 29–62 | Fragment/sequence validation, provenance and native acceptance |
| E17 | [pipeline_graph.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/pipeline_graph.py:175): 175–259; [active_learning_controller.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/engine/workflow_engine/workflows/active_learning_controller.py:183), `_docking_params`: 183–229 | Actual oracle, preparation/box/reference-ligand choices, sampling/training defaults |
| E18 | [active_learning_controller.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/engine/workflow_engine/workflows/active_learning_controller.py:1287), `_OpenVSDenseModel`: 1287 onward; `run_active_learning_model_step`: 1313–1523 | Labels, splitting, warm starts, feature memory, training and streamed exploitation/exploration |
| E19 | [active_learning_controller.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/engine/workflow_engine/workflows/active_learning_controller.py:497), `_ranked_rows`: 497–501; `_top_vector`: 504–516; `_collect_best_poses`: 554–625; campaign loop 723–1090 | Score accumulation, identity-aware stopping, output/retention and transfer |
| E20 | [graph_pipeline.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/engine/workflow_engine/workflows/graph_pipeline.py:881), `_cache_valid`: 881–912; key construction 1030–1089, reuse 1118–1122 | Existing reuse strengths and unhashed scientific dependencies |
| E21 | [graph_pipeline.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/engine/workflow_engine/workflows/graph_pipeline.py:762), `_submit`: 762–788; `enqueue_graph_pipeline`: 943–958; [active_learning_controller.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/engine/workflow_engine/workflows/active_learning_controller.py:1093), `_dock_round`: 1093–1238, `_dispatch_docking_shard`: 265–379 | Local child submission versus explicit remote docking sharding |
| E22 | [campaigns.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/engine/workflow_engine/workflows/campaigns.py:169): 169–212; [design_campaigns.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/engine/workflow_engine/workflows/design_campaigns.py:3781): 3781 onward, refinement 1722/2610 and evaluation 2984 | Configurable lineage/campaign boundary; not a second public fixed-template inventory |
| E23 | [bound_ligand_md.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/engine/workflow_engine/vendor/ligand_runtime/mn_ligand/workflows/bound_ligand_md.py:1826), `run_ligandx_md`: 1826 onward; MD manifests under `engine/workflow_engine/manifests/` | Standalone MD capability, absent from public starter graphs |
| E24 | [test_sabr_abmpnn.py](/home/user/programs/git-projects/mn-molecular-design-studio/tests/test_sabr_abmpnn.py:58): 58–211; [test_pipeline_graph.py](/home/user/programs/git-projects/mn-molecular-design-studio/tests/test_pipeline_graph.py:22): 22 onward; [test_native_generators.py](/home/user/programs/git-projects/mn-molecular-design-studio/tests/test_native_generators.py:159): 159–195 | Graph contracts, CDR/motif handling and hidden native workflow coverage |
| E25 | [molecule_imports.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/molecule_imports.py:529), `parse_structure`: 529–579 | Component reconstruction, modified-residue treatment and omission of single-atom non-polymers |

### H2. CPU-only verification performed

The repository was initially clean. Tests were run with `PYTHONDONTWRITEBYTECODE=1`, `CUDA_VISIBLE_DEVICES=''`, Python `-B`, and pytest `-p no:cacheprovider`, using `/home/user/mambaforge/envs/mn-molecular-design-studio/bin/python`.

| Run | Test files | Outcome |
|---|---|---|
| 1 | `test_pipeline_graph.py`, `test_active_learning_graph.py`, `test_active_learning_stopping.py`, `test_active_learning_failures.py`, `test_bindcraft_boltzgen_contracts.py`, `test_sabr_abmpnn.py`, `test_native_generators.py`, `test_kejia_idr.py`, `test_peptide_trap.py`, `test_rfdiffusion_beta_pairing.py` | 157 passed in 17.78 s |
| 2 | `test_rfd4_proteina.py`, `test_mpnn.py`, `test_quality_tools.py`, `test_bulk_screening_collection.py`, `test_cofolding_results.py`, `test_active_learning_execution.py` | 153 passed, 10 skipped in 5.31 s; skipped tests require `MDS_QC_DOCKER_TESTS=1` |

Additional in-memory/temporary-directory probes:

1. Constructed and called `validate_graph()` on every public starter and all 12 additional explicit IDs. Public results: 20 passed topology validation, VHH failed the binder-length rule, blank failed its deliberately unconnected initial state. All 12 additional factory/continuation graphs passed topology validation. This did not supply real scientific inputs or prove runnable model availability.
2. Enumerated cofold branch kinds: 15 monomer, 17 PPI, nine each ligand/DNA/RNA, one cyclic-partner branch.
3. Built Boltz-2 requests before/after changing legacy sample/recycle/full-PAE controls. Submitted sample/recycle counts and option dictionaries were unchanged, as described in F6.
4. Built a minimal RFD4-style two-protein component candidate in a temporary directory, passed it through `node_request('boltz2', ...)`, then `_prepare_prediction(..., resolve_random_seed=False)`. Confirmed the empty role/sequence fields in F2. `_prepare_prediction()` normalizes/checks inputs without starting a job or fetching MSAs. No actual GPU allocation or computation occurred.
5. Supplied sequences to the unmodified PPI cofold graph and called `validate_graph(..., require_inputs=True)`: preflight rejected the missing target PDB required by its AF2 branch. Supplied an ATP CCD input to the unmodified ligand cofold graph: preflight rejected Chai-1's incompatible CCD notation. Both reproduce the need to choose branches from the actual inputs before execution.

These checks establish adapter/graph behavior. They do not validate docking accuracy, binder success rates, catalytic activity, native container internals, current remote computer availability, or scientific threshold calibration. No tests were added and no source was repaired during this audit.

### H3. Scientific interpretation of metrics

| Measurement | Defensible interpretation | Interpretation to avoid |
|---|---|---|
| ProteinMPNN/AbMPNN likelihood or confidence | Sequence compatibility under that model and backbone/conditioning | Binding affinity, improved affinity, solubility or experimental success probability |
| pLDDT | Predicted local structural accuracy, with engine-specific scaling | Thermodynamic stability or correctness of partner orientation |
| PAE, ipTM, ipSAE | Predicted relative/interface structural confidence under specified chain pairing/cutoffs | Independent binding assays or a product of independent experimental probabilities |
| Target-aligned binder RMSD | Agreement with a reference pose under a stated mapping and conditioning mode | Accuracy against experiment when the reference is itself a generated hypothesis |
| Rosetta interface score/energy | Model/protocol-dependent packing/interaction assessment | Directly measured binding free energy or calibrated cross-target affinity |
| Docking score | A ranking proxy for a fixed receptor/preparation/scoring protocol | Measured kcal/mol binding free energy or selectivity across unrelated protocols |
| Active-learning sigmoid | Surrogate prediction of membership in the chosen docking-score class | Calibrated probability of experimental binding |
| Boltz-2 affinity outputs | Separate trained small-molecule affinity/classification predictions in their supported domain | Protein–protein affinity, catalytic activity, or equivalence to ipTM/pLDDT |
| Short MD stability/contact persistence | Behavior over the sampled trajectory with a specified force field and starting state | Binding proof, equilibrium occupancy, or converged affinity |

### H4. Literature and primary method sources

The sources below support method assumptions and limitations. They do **not** validate Studio's complete combinations or transfer published success rates to its targets. Author documentation is identified separately from published evidence. No unsupported citation metadata has been invented.

| ID | Source checked | What it supports and limits |
|---|---|---|
| S1 | Dauparas et al. (2022), **Robust deep learning-based protein sequence design using ProteinMPNN**, *Science* 378, 49–56, [DOI:10.1126/science.add2187](https://doi.org/10.1126/science.add2187); [author-hosted paper](https://www.bakerlab.org/wp-content/uploads/2022/09/Dauparas_etal_Science2022_Sequence_design_via_ProteinMPNN.pdf) | Backbone-conditional sequence design with experimental examples; does not make sequence likelihood an affinity/stability measurement. |
| S2 | [LigandMPNN author repository](https://github.com/dauparas/LigandMPNN) | Atom-context-aware sequence design and fixed/mutable residue controls. Using its runtime with ProteinMPNN weights does not itself make a graph ligand-conditioned. |
| S3 | Watson et al. (2023), **De novo design of protein structure and function with RFdiffusion**, *Nature* 620, 1089–1100, [DOI:10.1038/s41586-023-06415-8](https://www.nature.com/articles/s41586-023-06415-8); [author sequence-design guidance](https://github.com/RosettaCommons/RFdiffusion#sequence-design-for-binders) | Experimental support for backbone generation → sequence design → prediction/filtering; author guidance permits ProteinMPNN without mandatory FastRelax. Does not validate every target/threshold or the local initial-guess configuration. |
| S4 | Abramson et al. (2024), **Accurate structure prediction of biomolecular interactions with AlphaFold 3**, *Nature* 630, 493–500, [DOI:10.1038/s41586-024-07487-w](https://www.nature.com/articles/s41586-024-07487-w); [EMBL-EBI interpretation guidance](https://www.ebi.ac.uk/training/online/courses/alphafold/alphafold-3-and-alphafold-server/how-to-assess-the-quality-of-alphafold-3-predictions/) | Biomolecular structure prediction and documented geometry/state/dynamics limitations. Confidence is structural; sampled static predictions are not automatically solution ensembles. |
| S5 | [DunbrackLab/IPSAE author implementation](https://github.com/DunbrackLab/IPSAE) | ipSAE derives interface evidence from structure and pairwise PAE, with explicit cutoffs. It reuses prediction information. The cited February 2025 [preprint endpoint](https://www.biorxiv.org/content/10.1101/2025.02.10.637595v1) returned access denial during this audit; detailed quantitative performance claims from that paper remain externally unverified here. |
| S6 | Zhou et al. (2024), **An artificial intelligence accelerated virtual screening platform for drug discovery**, *Nature Communications* 15, 7761, [DOI:10.1038/s41467-024-52061-7](https://www.nature.com/articles/s41467-024-52061-7); [author OpenVS code](https://github.com/gfzhou/OpenVS) | Published precedent for docking-guided iterative screening with exploration and experimental follow-up. Its RosettaVS oracle differs from Studio's UniDock Pro oracle. The catalogue's [Novartis/OpenVS reference](https://github.com/Novartis/OpenVS) returned 404; update that attribution. |
| S7 | [Kejia peptide/IDR binder author repository](https://github.com/drhicks/Kejia_peptide_binders) | Author-described flexible-fragment pipeline and the explicit caution about high AF2 passing rates for strand-pairing templates. Used as primary method guidance; no local experimental equivalence established. |
| S8 | Sahtoe et al. (2024), **Design of amyloidogenic peptide traps**, *Nature Chemical Biology* 20, 981–990, [DOI:10.1038/s41589-024-01578-5](https://doi.org/10.1038/s41589-024-01578-5); [publication record](https://pubmed.ncbi.nlm.nih.gov/38503834/) | Identity and scoped experimental rationale of the peptide-trap task. Direct publisher full-text fetch failed; the complete paper's selection protocol was not independently rederived here. Local code itself labels its protocol an adaptation. |
| S9 | [PXDesign author repository](https://github.com/bytedance/PXDesign) | Distinct generation-only, preview and extended AF2/Protenix evaluation modes. Author recommendations are not evidence that mandatory double prediction improves every Studio campaign. |
| S10 | [Genie3 author repository](https://github.com/aqlaboratory/genie3) | Task scope and author method context. The Studio AF2 adapter is not assumed to reproduce every upstream evaluation detail. |
| S11 | [Protpardelle-1C author repository](https://github.com/ProteinDesignLab/protpardelle-1c) | Unconditional/conditioned generation, sequence design and ESMFold self-consistency context. Monomer evaluation does not establish interface binding. |
| S12 | [Biohub ESM get-started and binder-design protocol](https://biohub.ai/esm/protein/get-started) | ESMFold2 structure and sequence-search/binder-design context. The catalogue's GitHub tutorial directory could not be fetched; this current primary documentation was accessible. No published VHH-specific success claim is inferred from the displayed scFv/minibinder examples. |
| S13 | [RFD4-Proteina author repository](https://github.com/RosettaCommons/RFD4-Proteina), including its linked Python/enzyme tutorials | Rich conditioning/co-design task scope. Repository identifies **Biomolecular design through latent reasoning with RFD4-Proteina** as a 2026 preprint. Tutorial capability does not establish catalytic activity or chemical-fidelity of Studio's subsequent handoff. |
| S14 | [BindCraft2 author repository](https://github.com/PacesaLab/BindCraft2) | Native sequence optimization, ProteinMPNN redesign and AF2 validation are a compound protocol. Studio records a source commit; current upstream docs are method context, not verification of every packaged native stage. |
| S15 | [BoltzGen author repository](https://github.com/HannesStark/boltzgen) | Native final-budget/diversity selection and design/refolding pipeline, which should not be replaced blindly by generic sequence redesign. |
| S16 | Braun et al. (2019), **Best Practices for Foundations in Molecular Simulations**, *Living Journal of Computational Molecular Science*, [author/journal guide](https://livecomsjournal.org/index.php/livecoms/article/download/v1i1e5957/939/2729) | Need for appropriate preparation, sampling, uncertainty and reproducibility. The audit's recommendation against obligatory short MD is a methodological judgment, not a benchmark of the local MD adapters. |
| S17 | [Boltz author repository and affinity guidance](https://github.com/jwohlwend/boltz) | Separates small-molecule `affinity_probability_binary` from quantitative `affinity_pred_value`, with different supervision/use cases. Neither is a protein–protein binding confidence head. Use exact model/version documentation for units; do not substitute either output for a measured Kd. |

### H5. Remaining verification needs

The following remain open rather than being asserted as demonstrated scientific outcomes: prospective wet-lab utility of each default; equivalence of packaged native engines to current author protocols; accuracy/calibration of engine-specific thresholds; magnitude of active-learning validation optimism; optimal sample/recycle budgets; preservation of every unusual ligand/cofactor/bond topology; and measured cross-workstation throughput for protein candidate collections.

These uncertainties justify focused tests or small benchmark studies after contract repair. They do not justify mandatory additional prediction, MD, or scoring for every candidate.
