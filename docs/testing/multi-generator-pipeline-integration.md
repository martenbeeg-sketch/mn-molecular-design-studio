# Phase 7A — Multi-generator scientific graph integration

Date: 28 October 2026. Branch: `fix/scientific-pipeline-validation`.
Starting HEAD: `fcc4dcc` (including the newer `2fb4cdc` native execution-job ledger).
This report describes executable implementation and CPU/retained-artifact evidence.
It makes no new inference-throughput or prospective binding-quality claim.

## Implemented graph and entry point

**Multi-Generator Binder Design & Evaluation** is available in the existing
Pipelines setup panel and `STARTER_TEMPLATES`, under ID `multi_generator_binder`.
Select generator and predictor tiles there, create the editable graph, set the
Target Structure, configure native tools using their existing forms, and enable
branches under **Branch and processing settings**. The first predictor starts
enabled; additional branches are explicit opt-ins. The normal save, validate,
run, monitor and Results actions apply. No second canvas, workflow engine,
campaign controller or scheduler is introduced.

```text
Target Structure ─┬─ RFD4-Proteina (ordinary complete co-design workflow) ─┐
                 └─ BindCraft2 (ordinary complete native pipeline) ─────┤
                                                                      v
                                                            Candidate collection
                                                                      |
                                                 Initial ProteinMPNN when needed
                                                                      |
                                                 Optional SolubleMPNN variants
                                                      (disabled by default)
                                                                      |
                                                         Shared scientific filters
                                                                      |
                                                ┌─────────────────────┴────────────┐
                                             Boltz-2                         Protenix V2
                                            (enabled)                    (enable explicitly)
                                                └─────────────────────┬────────────┘
                                                                      v
Candidate collection + filter review records ─────────────> Unified candidate analyzer
                                                                      |
                                                           Ranking and diversity
                                                                      |
                                                          Existing Results + ZIP
```

The original collection also enters the analyzer, so original backbones and
native-rejected/unresolved records remain visible after routing. The filter review
input supplies original and redesigned variants plus every shared-filter decision.
The graph is tool/target independent. Production code contains no PD-L1 filename,
sequence, hotspot, binder-length or acceptance-threshold special case.

The conceptual composition was inspired by Figure 2A of
[Claude Science et al., *De novo design of protein binders against scorpion toxins
with Claude*, 20 October 2026](https://www-cdn.anthropic.com/31ac578806cffd8b3d92d8affe6edc4c04d502e1.pdf).
The figure was inspected. Its campaign, targets, biological screening, thresholds
and agent are not implemented here. This is an Anthropic technical report; “Science”
in its author name is not a claim of publication in the journal *Science*.

## Reused generators and actual output capabilities

Inventory comes from the registered graph contracts and existing normalizers in
`workflows/generators.py`, `design.py`, `bindcraft2.py`, `rfd4_proteina.py` and the
dedicated Complexa/ESMFold2 adapters. Registration/source verification is distinct
from new GPU verification. Campaign-adapter eligibility is not a prerequisite.

| Existing tile / workflow | Actual ordinary output contract | Pipeline handling |
|---|---|---|
| BindCraft2, `design.bindcraft2` | Native gradient design, MPNN, validation and filters; sequence-complete evaluated rows, possible rejected outcomes, or no normalized rows after trajectory rejection | Whole native workflow reused. Original sequences/scores retained. No automatic external redesign |
| RFD4-Proteina, `design.rfd4_proteina` | Guided binder co-design emits sequences and complexes with component/chain mappings, conditions and native metrics; authoritative campaign acceptance unavailable | Original sequence/reference/constraints retained. Explicit task options, including refolding, are respected |
| RFdiffusion, `design.rfdiffusion_generate` | Backbone proposals; residue identities in coordinates are not a designed-sequence contract | Initial ProteinMPNN required |
| RFdiffusion3 / Foundry, `design.rfdiffusion3` | Generation-only backbone records; existing standalone/native modes have separate controls | Initial ProteinMPNN; no new campaign adapter |
| Genie3, `design.genie3_generate` | Generation-only backbones | Initial ProteinMPNN |
| Genie3 native saved pipeline, `design.genie3` | Native inverse folding and evaluation; raw/intermediate and evaluated candidate pools where available | Native pipeline preserved; intermediates retained for review |
| PXDesign, `design.pxdesign_generate` | Backbone proposals, even where structures contain residue letters | Initial ProteinMPNN |
| PXDesign native saved pipeline, `design.pxdesign` | Native MPNN and AF2-IG; extended preset also Protenix; raw/evaluated pools | Native pipeline preserved |
| BoltzGen, `design.boltzgen` | Generation, inverse folding, refolding, evaluation/filtering; multiple raw/evaluated pools | Complete workflow preserved |
| BoltzGen Generator | Native design-stage sequence/structure records; later sequence/evaluation stages are separate | Preserve supplied sequences; independent validation remains separate |
| Proteina-Complexa | Native generation/filtering/evaluation, including verified native redesign choices and unresolved acceptance when evidence is absent | Retain native chosen sequences and all raw evidence |
| Proteina-Complexa Generator | Backbone and sequence, explicit target/binder roles | Preserve supplied sequences |
| Protpardelle generator | Generation-only backbone proposals | Initial ProteinMPNN |
| Protpardelle native saved pipeline | Several MPNN variants per scaffold and native ESMFold measurements; raw samples retained | Preserve explicit native design/sequence indices and evaluations |
| ESMFold2 Binder Generator | Gradient-guided sequence design with native predicted complex | Preserve original sequence/native design-loop evidence |
| RFdiffusion β-pairing / Kejia IDR | Structure-bearing backbone proposals | Their existing constraints/preparation remain required; initial MPNN supplies sequences |
| Peptide Trap | Sequence-complete native Rosetta-threaded/refined designs | Retain native sequence/evidence; existing target-sequence input still required |

The compatibility builder lists registered tiles, not a promise that every task
mode accepts a shared protein target. Existing tool preflight checks required
scaffolds, task settings and chemistry. A generator's separate target override
must match the common target. Unsupported shared hotspot/coldspot constraints
fail preflight; an explicit Binder Design Context can change a branch's constraints.
Legacy standalone workflows retain their existing behavior.

## CandidateSet contract and scientific identities

`studio.candidate-set.v1` is additive metadata in existing candidate-v1
`raw_metadata.candidate_set`, accompanied by a hashed JSONL collection manifest.
No existing candidate schema or saved graph schema is migrated.

Identities distinguish target, generator, original source run, native attempt
(when supplied), backbone/design, sequence variant, prediction, analysis and
selection. Original native IDs, parent IDs and artifact references remain present.
Backbone identity uses explicit native backbone/design/parent identifiers where
available; otherwise the source candidate is its own design identity. Missing
native trajectory identity remains unavailable, rather than inferred from row count.

Target identity records the original file SHA-256, selected chains, their exact
sequences and author residue identifiers, and interface constraints. The existing
native adapters perform their own preparation and chain mapping. Target staging
copies bytes; it does not renumber, rebuild or change constraints.

Each variant retains its original sequence/structure, native filtering, original
source settings and worker allocation, artifact hashes, and parent-variant ancestry.
Native and shared filters are distinct. Prediction evidence retains the complete
normalized worker row, raw seed/protocol/model/checkpoint information when supplied,
input membership, worker/GPU attempts, output hashes and independent prediction ID.
Unrecorded model/attempt information is not manufactured.

Collection merge uses a temporary SQLite index. Policies are **keep_all** or
**mark_duplicates**. Marking creates an exact-sequence duplicate group/count while
retaining each generator/run history. Conflicting records for one stable identity
are errors. Selection can choose one exact-sequence representative per comparable
protocol without deleting other histories. Shared filtering, partitioning, sequence
variants and predictor fan-out/fan-in all operate on this same transport.

## Collection batching and resources

The existing `candidate_execution` controller is reused. Each predictor defaults to
chunks of five, at most two in flight, with explicit partial-result handling.
Generation nodes explicitly retain ordinary native mode. Existing native generation
controls and one-attempt campaign dispatch are unchanged.

These settings mean different things:

| Level | Meaning / implementation |
|---|---|
| Native generator batching | Existing ordinary engine design/trajectory/replicate controls; unchanged |
| Candidate collection chunking | Bounded sets of normalized inputs for downstream collection adapters |
| Scheduler jobs | Existing adapter workflow requests, resource leases, attempts and workstation dispatch |
| GPU microbatching | Existing engine controls; never equated to candidate chunk size |
| Acceptance campaign attempts | Current supported native job size remains one; Phase 6C.2 accounting does not authorize larger native dispatch |

Per-step controls expose chunk size, maximum concurrent jobs, destination,
GPU selector and CPU/VRAM requirements. Native sampling/sample multiplicity stays
in the tool form. The original pipeline-wide controls remain available. A per-step
policy overrides the pipeline policy. Unsupported explicit sharding and invalid
selectors fail validation; actual workstation/image/checkpoint eligibility and
physical GPU reservations remain the existing scheduler's responsibility.

The known hardware pool remains pc002975 (two RTX 4090s), pc004032 and pc004033
(one RTX 5090 each). This phase neither deploys source to these peers nor claims
fresh cross-workstation inference. Generation and predictor concurrency remain
independent. No requirement to occupy all GPUs is added.

Sharding freezes exact candidate-to-job membership and effective requests in
checksummed manifests. A 12-input collection with chunk size five is **5 + 5 + 2**.
Optimized predictors additionally record expected sample/seed multiplicity, returned
and valid output counts, and per-input state. A returned row is not automatically
all expected samples. Failed chunks do not erase successful chunks. Existing durable
submission tags, scheduler retry fences, peer receipts, recovery and cancellation
are reused. Stable CandidateSet shard output IDs derive from scientific identity
and native sample identity, not completion order. Conflicting prediction receipts
are refused and duplicate identical receipts reconcile idempotently.

CandidateSet JSONL reads, merges, filtering, sequence routing, batch aggregation,
analysis joins and ranking are streamed or disk indexed. Requests materialize at
most one processing chunk. Existing legacy native adapters and the final normal
Studio result.json/candidate writer can still materialize their input/final rows;
this is not an unbounded-scale replacement for all old native interfaces.

## Conditional sequence handling and filters

The existing shared ProteinMPNN workflow supplies initial sequences for
`generation.backbone` records; coordinate residue letters alone do not make those
sequence-ready. Co-designed/sequence-complete records bypass initial inference.
An entirely sequence-ready collection passes through without an MPNN scheduler job.

Optional SolubleMPNN remains disabled until explicitly enabled. Existing interface,
fixed-residue and constrained-position handling is reused. Originals remain available;
child variants get independent IDs and ancestry. Parent prediction/filter evidence
is retained as history and is not assigned as a child's measurements or acceptance.
Rejected, unresolved, intermediate or artifact-missing inputs are retained for review
and excluded from inappropriate expensive sequence/prediction requests.

Shared filters preserve **passed**, **failed**, **not_evaluated**,
**missing_evidence**, **invalid_metadata**, and **unresolved** states. Rules contain
exact measurement field, operator and threshold. Defaults add no new scientific
threshold. Native-rejected candidates are excluded from prediction by default,
but remain in the consolidated collection. Corrupt JSONL is an operational error;
a missing artifact is missing evidence, not scientific rejection. Nonfinite values
are retained as explicit invalid textual evidence in standards-compliant JSON.

## Predictors, analyzer, ranking and diversity

Boltz-2 reuses its verified optimized-prediction handoff and collection task adapter.
The additional default tile is **Protenix V2**, through the existing
`refolding.optimized_prediction` adapter, not the older Protenix wrapper.
Installed images were checked read-only:

| Adapter | Installed image | Verification in this phase |
|---|---|---|
| Boltz-2 | `mn-boltz2-opt:2.2.1-cu130-batch2`, ID `3054688e9190…` | Real retained P4/P5 evidence + deterministic chunk workers |
| Protenix V2 | `mn-protenix-opt:2.0.0-cu128-batch2`, ID `f9b86d5f1365…` | Existing collection adapter/request signature and installed image; deterministic worker fan-out/fan-in |
| Protenix V1 | `mn-protenix-v1-opt:1.1.0-cu130-batch2`, ID `d17846d1dfb6…` | Existing profile and installed image inspected; no fresh inference |
| ESMFold | `mn-esmfold:2.0.1-cu128-sm89`, ID `de83f59303fc…` | Existing multimer-FASTA profile/image inspected; no fresh inference |
| ESMFold2 / Fast | `mn-esmfold2-opt:0.1.0-cu130-batch1` / `…-fast-batch1`, IDs `f34fcd468c48…` / `e98b1fcdd8ee…` | Existing collection profiles/images inspected; no fresh inference |

ESMFold-family sequence-only predictions are not evidence of an MSA-based interface
confidence score. The binder-alone ESMFold tile remains distinct. Installed images
and request compatibility do not establish new GPU validation of all profiles.

The analyzer requires the source variants and joins predictions by variant ID,
then verifies exact sequence identity. It retains every independent seed/sample,
raw measurement validity, failure, disabled branch and missing predictor entry.
Existing analysis nodes can contribute separate analysis-result evidence. Optional
target-aligned pose RMSD uses the existing calculation only where an actual input
reference is recorded. DockQ/scDockQ, prediction agreement or clustering scores are
not fabricated when absent; supplied valid measurements remain available.

Measurement names are qualified:

```text
prediction/<branch>/<prediction-id>/<metric>
prediction/<branch>/<metric>                  # only one unambiguous valid output
analysis/<analysis-node>/<analysis-result-id>/<metric>
```

Several seeds do not create an implicit average/best-seed alias. The rank node
separates measurement collection/validity, thresholds, ordered multi-measurement
ranking, exact-sequence diversity and final selection. Ranking rules explicitly
declare higher/lower direction. Lexicographic ranking does not equate iPTM, ipSAE,
DockQ or native scores. Recorded protocol cohorts reuse the existing scientific
selector's comparability policy; missing protocols are conservative singleton
cohorts. Raw values and all unselected decisions remain available. Existing cluster
selection/filter nodes can still use externally supplied cluster measurements;
no new clustering/scientific algorithm is introduced.

## Results, export, cache and resume

The normal Results measurement/selection table is extended through its existing
`candidate_state_rows` projection. It shows generator/source run, attempt/backbone/
variant identity, original sequence/structure, native/shared filtering, independent
predictor outcomes and measurements, rank/selection, duplicate evidence and worker
allocation. Full candidate/analysis/model evidence stays in JSONL and Files.

The existing graph run monitor shows distinct stage counts and chunk membership,
completed inputs, retries, worker, GPU and errors. Native attempts/rejections are
reported only when source evidence provides them. Normalized candidates, sequence
states/variants, filter outcomes, processing jobs, predictions, analyzed candidates
and selected candidates are separate counters, not one “designs completed” number.

The regular run ZIP includes `artifacts/candidate_set/` with the collection manifest,
all candidate records, CSV identity index, predictor evidence, analysis/selection
records, graph/execution snapshot, selected FASTA and copied selected native/predicted
structures, source-to-export artifact mapping and checksums. Rejected/unresolved
records are retained. Existing selection/run export infrastructure is reused.

Graph content caching remains authoritative. Effective scientific inputs, settings,
implementation/runtime identities, variant sequences and artifact hashes determine
reuse. Analyzer/ranking edits do not invalidate upstream generation or predictions.
Completed independent branches survive partial failures and controller restart.
Inline RFD4 `options_json` is explicitly treated as content by portable-path handling,
fixing a restart bug that previously resolved `{}` as a filesystem location.

## Scientific evidence and regression results

[Retained execution evidence](evidence/phase7a-retained-graph.json) identifies the
actual saved graph/run, exact original source artifacts, sequence hashes,
Boltz measurements, export and stage counts. The replay uses normal
`run_graph_pipeline`, target fan-out, generator tile retained-source adapters,
collection merge, conditional sequence routing, shared filtering, predictor branch,
analyzer, ranking and Results. No hand-edited Results file substitutes for execution.

Unmarked older evidence is imported into new project-scoped snapshots with absolute
artifact references and original hashes; originals and their ownership files are
not modified. The reusable retained-source option checks exact managed-run access,
project membership, successful execution and matching generator/predictor identity.
Prediction replay requires original candidate ancestry plus sequence equality and
exact target-chain sequence correspondence; matching a sequence alone is insufficient.

Sources are the established 115-residue chain-A PD-L1 reference with hotspots
A40/A99/A107, P4 campaign `20261028-101057-a07ef0a4`, P5 campaign
`20261028-124448-350f7294`, and rejected BindCraft2 campaign
`20261028-085821-36356a14`. The BindCraft2 benchmark used a 60-residue binder setting
and produced zero accepted normalized candidates. Its two rejected trajectories
remain source evidence. No accepted BindCraft2 result is invented.

| Retained RFD4 binder | Retained Boltz iPTM | Evidence |
|---|---:|---|
| P4, 75 residues | 0.42641 | Existing P4 prediction, exact original candidate ancestry |
| P5, 116 residues | 0.46578 | Existing P5 prediction, source candidate `dc_00000001_d68026ec8a3e62d5` |
| P5, 94 residues | 0.65414 | Existing P5 prediction, source candidate `dc_00000002_d68026ec8a3e62d5` |

This retained run performs **zero new generation/prediction inference**, downloads
no models and submits no scheduler children. Additional predictor execution and
batch behavior use explicitly labeled deterministic workers through the real graph
and ordinary request/manifest/receipt controller interfaces.

Deterministic collection sizes **0, 1, 2, 5, 10, 25, 50, 100** exercise two
predictor branches, stable identities and exact five-member chunks. Additional
cases cover unequal generator sizes, empty/native-rejected outputs, duplicate
sequences from different sources, backbone-only routing, three child sequences per
backbone, 5+5+2 remainders, partial chunk failure, retries, restart/idempotence,
missing samples, different seeds, conflicting receipts, cross-sequence refusal,
missing/nonfinite evidence, scientific rejection, ranking/diversity, export,
disabled branches, analyzer-only cache edits, target/project preflight and old graphs.

UI verification uses Streamlit AppTest on the existing pipeline builder and Results
table projection, checking branch controls, stage counts and consolidated rows.
API verification uses the existing authenticated catalogue/save/validate endpoints,
checks collection execution-plan output and source guards, and asserts no job is
submitted during validation. Existing access-control regressions remain required.
Browser interaction and fresh native GPU inference are separate from these CPU tests.

Final verification ([machine-readable summary](evidence/phase7a-verification-summary.json)):

* Full CPU sweep: **1,768 passed, 33 skipped, zero failed**, 380.51 s. Skips are
  existing opt-in GPU/container or optional dependency checks. No GPU opt-in was set.
* Focused graph/candidate/peer/MPNN/task compatibility sweep: **211 passed**,
  190.73 s. This overlaps the full sweep.
* Final accounting/cache/Results/API/analysis follow-up: **20 passed**,
  15 deselected, 35.90 s; existing-selector preservation regression: **1 passed**,
  22 deselected, 3.55 s. These verify final small changes after the broad sweep;
  they are overlapping checks, not additional full-suite totals.
* Real retained graph: `campaign/20261028-164454-50d74667`, saved graph
  `934d533d05f043ffb51bc378d73f6c57` (the authoritative graph ID is in the JSON
  evidence). Three native source variants join three real retained Boltz predictions.
  Collection counts distinguish five historical native attempts, two rejected native
  attempts, three normalized/sequence-ready variants, three completed predictions,
  three analyzed variants and three configured export selections.
* Existing pipeline UI with the real retained run: zero AppTest exceptions;
  the builder, branch controls, plan, stage counts and three-row Results projection
  render. Actual managed-run/project API preflight also succeeds without mocked
  source access guards or job submission.
* `git diff --check`, changed-source syntax parsing, export checksums, exact P4/P5
  ledger hashes and original local-edit checks pass. No original scientific artifact
  was changed, and no new native inference, model download, service restart or
  remote deployment occurred.

Re-run the deterministic graph regressions with the existing Conda environment:

```bash
python -m pytest tests/test_candidate_sets.py tests/test_candidate_execution.py -q
```

`scripts/verify_multi_generator_pipeline.py --help` documents the retained replay
arguments. Its source inputs and project identity are recorded in the committed
retained-evidence JSON; the script imports new snapshots, executes the ordinary
pipeline, verifies original hashes and creates the normal scientific export.

## Limitations and next adapters

* Native output fidelity remains each existing normalizer's responsibility. A
  rejected native trajectory with no normalized candidate is an attempt/source
  record, not an invented candidate. Unknown attempt/backbone relationships remain
  unknown; source candidate IDs are used conservatively where no native mapping exists.
* Initial sequence generation and explicit redesign use existing MPNN scientific
  constraints. Tool-specific unsupported chemistry/chain roles fail the existing
  preparation checks. This phase does not add new scientific algorithms.
* Retained evidence reuse is an explicit provenance operation; its original native
  settings govern those results. Editing tile controls does not relabel retained
  output as a fresh computation. Clear the retained-run field to request fresh work.
* New collection operations use disk indexes/bounded iteration, while existing native
  adapters and final Results serialization retain their prior in-memory boundaries.
* Selection provides ordered measurements and exact-sequence diversity. Structural
  clustering, new agreement scores and richer ensemble aggregation need supplied
  measurements/adapters; they are not claimed as newly computed here.
* Only Boltz-2 retained scientific predictions are joined in the real replay.
  Protenix V2 fan-out is tested with deterministic workers, not new Protenix inference.
* Peers must pass their ordinary exact runtime compatibility checks before using
  newly changed source. No remote deployment or service restart was performed.

Adding a generator with ordinary normalized outputs requires its existing tile and
input-preparation contract, plus an output adapter/capability declaration if its
records lack sequence state, chain roles or native identity/evidence. It does not
require a Phase 6B campaign adapter. Adding a predictor requires an existing
collection-capable request adapter with exact source-candidate/sequence lineage,
model/seed/artifact provenance and real input cardinality; then connect it to the
analyzer. Unsupported explicit sharding is refused rather than silently emulated.

Phase 6C.2 still needs verified multi-attempt native dispatch, per-attempt seeds and
receipts, native partial recovery/retry semantics and performance measurements before
larger generator campaign jobs can be authorized. Its durable accounting ledger is
retained. Persistent GPU model workers and a new scheduler are outside this phase.

## Changed files and commits

Implementation: new `candidate_sets.py` and `candidate_pipeline.py`; scoped changes
to `pipeline_graph.py`, `pipeline_graph_ui.py`, `scientific_tasks_ui.py`,
`candidate_execution.py`, `result_presentation.py`, `api.py`, core
`portable_paths.py` and `workflows/graph_pipeline.py`.
Verification: `tests/test_candidate_sets.py`, updated catalogue expectation in
`tests/test_scientific_tasks.py`, `scripts/verify_multi_generator_pipeline.py`,
this report and `evidence/phase7a-retained-graph.json`.

Pre-existing `run_app.py`, untracked `tool_references.py` and paper-reference work
are excluded. Retained P4/P5 scientific artifacts, thresholds, native defaults,
standalone tools, prior templates and distributed campaigns are preserved.
No additional Codex agents, GitHub push or remote runtime deployment were used.

Implementation commit: **`9549be6` — `feat: compose batch-aware multi-generator scientific pipelines`**.
The following documentation commit contains this report and both evidence JSON files;
its hash is reported in the completion message and resolves with
`git log -1 --format=%h -- docs/testing/multi-generator-pipeline-integration.md`.

## Phase 7A execution-layer verification (follow-up)

The earlier graph partitioning tests did not establish native process grouping.
This follow-up exercises the existing scientific adapters at their process launch
boundaries and the installed Boltz/Protenix image-side batch adapters, replacing
execution with deterministic fixtures. No GPU inference, model download, remote
runtime deployment or push was performed. All retained source hashes still match.

### Generality and capability boundaries

`native_execution.tool_inventory()` enumerates **every registered generator,
sequence-design and prediction tile**, including hidden compatibility tiles,
MPNN checkpoint variants, Caliby models/ensemble modes, and RFD4's existing
inverse-folding, structure-prediction, sidechain-packing and diversification
operations. The inventory is exposed in the canvas execution preview and pipeline
validation API. Registry membership does not certify native model reuse. Inspected
contracts, fixture-tested paths and unverified lifecycles are explicitly different.
The generic graph, CandidateSet lineage, sequence eligibility, chunk controller,
fan-out/fan-in and receipt aggregation do not depend on the default two generators.

The pipeline factory accepts `sequence_designer` as well as generator/predictor
lists. Existing ESM-IF1, Caliby, BoltzGen inverse-folding and AbMPNN tiles can replace
MPNN; their own scientific eligibility/constraint checks still apply. Initial
sequence preparation filters backbone-only records for every sequence-design tile
and preserves native co-designed sequences. Explicit redesign creates variants.
Collection chunking now recognizes the existing inverse-folding adapters too;
ordinary workflows continue to pass the complete input list to their native runner.

| Tool/operation | Actual current execution contract | Model initialization / limits |
| --- | --- | --- |
| Ordinary RFD4 generation and existing-structure operations | One outer container; up to 100 replicates **per input** in each native command; larger replicate requests launch more native processes. Dataloader batch size is independent. | Generator checkpoint reused within each invocation; search/oracle/affinity work can add initialization. Native source and wrapper fixtures verified. |
| RFD4 inverse folding | Existing `existing_structure` + `operation=inverse_folding` supports several input files and replicates through the same native runtime. | Included in the inventory. A dedicated CandidateSet **sequence-stage** handoff is not verified. It is not interchangeable with the RFD4 structure-prediction tile. |
| Ordinary BindCraft2 | One campaign entrypoint/container handles the trajectory budget; accepted-design target remains distinct. | Reuse is per native worker. Auto worker count and validation-pool rebuilds prevent an exact model initialization forecast. Fixtures preserve ordinary budgets 1–100. |
| RFdiffusion / RFD3 | Existing native trajectory/vectorized sampling contracts retained, as documented in Phase 6C. | Source-verified lifecycle; no new multi-attempt campaign support. RFD3 rounding remains explicit in its existing planner. |
| Shared LigandMPNN backend, including ProteinMPNN/SolubleMPNN/AbMPNN/custom checkpoint selections | Compatible structures use native `pdb_path_multi`, `fixed_residues_multi`, and per-residue bias/omission maps in **one** model process. Different global settings split native groups. | One model load per compatible group. The native RNG stream spans structures; input order, seed, masks, parent IDs, group ordinal, all FASTA/backbone/packed/statistics artifacts and native output IDs remain recorded. No numerical-equivalence claim versus separate seeded processes. |
| Foundry MPNN backend | Existing one-structure command, one Docker/native process per structure within a scheduler chunk. | A chunk of five is **five native jobs**, five containers and repeated initialization. No verified multi-structure Foundry path is claimed. Backend/constraint-dependent plans show bounds until staging resolves them. |
| ESM-IF1 | One native runner/model process receives all inputs and sequences per input. | Model loads before the input loop; deterministic library fixtures verify one load for 1, 5, 10, 25, 50 and 100 inputs. |
| Caliby variants | One native runner receives all inputs; main Caliby model is constructed before the input loop. | Single/provided ensembles reuse the main model. Generated ensembles can initialize additional Protpardelle/MPNN models per input; that extra lifecycle is not counted as one verified model load. |
| BoltzGen inverse folding | One native Python process receives all inputs but instantiates/runs a native task for each input. | Initialization repeats inside the process. Combining native task configurations is a separate unverified ordinary-workflow optimization. |
| Boltz-2 optimized tile | One container per chunk; image-side adapter sends the complete YAML directory to one model entrypoint **per seed**. Samples per input are another axis. | Initialization repeats per seed/chunk. Protein complex fixture path verified against installed source. Multi-GPU/affinity legs can have additional model lifecycle costs. |
| Protenix V2 optimized tile | One container and native directory invocation per chunk, all JSON items included. | Installed source constructs the runner before its JSON loop. Native launch grouping verified with the installed image-side adapter. |
| Other selectable generators/cofolders and legacy prediction paths | Existing adapter contract preserved; inventory shows unverified internal lifecycle rather than inferring one model process from a batch API/container. | Deeper per-engine native lifecycle verification remains required before a model-reuse claim or optimization. No silent per-candidate workflow substitution is introduced. |

The RFD4/BindCraft2 lifecycle findings reuse retained Phase 6C source evidence.
Additional installed sources were copied read-only from stopped, **never started**
inspection containers. Their hashes are in the evidence JSON; third-party source
and weights are not committed. Boltz stock source loads a model before
`trainer.predict` over its data module; the optimized worker consumes the full
input directory. Protenix source constructs `get_default_runner` before iterating
JSON inputs. These are source findings, not throughput measurements.

### Accounting, forecast and recovery

`native_execution.json` records process launch windows, returned process status,
commands and exact input IDs separately from scheduler chunks. Model-call counts
for prediction use returned image-side receipts only where the native lifecycle
is verified. Unverified wrappers retain commands without inventing model counts.
RFD4's inner runtime records each native command and completion separately.
Unknown launch windows remain unknown. Allocation receipts are deduplicated by
allocation/worker identity; five containers under one scheduler allocation count
as one allocation. Scheduler reallocation is counted separately. A fixture verifies
three containers across two allocations, independently of native invocation count.

Per-node execution state/export/monitor reports candidates, chunks, scheduler
attempts, native invocations, Docker processes, GPU allocation events/slots and
prediction outputs separately. Missing receipts stay unknown. Read-only forecasts
include partial chunks and expected model initialization scope; counts depend on
eligible inputs and resolved backend/constraint groups. Preview settings neither
submit jobs nor mutate tile scientific controls. Example: 25 candidates, chunk 5,
Boltz three seeds and two samples gives **5 scheduler jobs, 5 containers, 15 native
invocations and 150 requested prediction outputs**, with concurrency still 2.

Failed local/remote jobs retain already normalized, artifact-backed member outputs
and do not rerun those members. Missing members and incomplete samples remain
explicit and the stage remains partial. LigandMPNN's grouped files are recovered
by exact filename boundaries; completed FASTA/backbone pairs can be normalized
without launching a model. ESM-IF1/Caliby/BoltzGen inverse runners atomically
checkpoint completed sequence records between samples/inputs; failed/interrupted
jobs can normalize that durable prefix. Fixed target/motif checks remain in the
existing normalizers. Remote partial outputs use existing checksum-verified artifact
transfers, including native receipts. An interrupted native job with unnormalized
raw artifacts is quarantined for reconciliation instead of blindly resampled.

### Execution fixture results

A deterministic suite drives real Studio request builders and scientific adapters,
intercepts actual Docker/native `subprocess.run` boundaries, and invokes copied
installed Boltz/Protenix batch adapters. It checks native multi-input files,
per-candidate constraints, output ancestry and returned sample counts, not just
Python list partitioning. MPNN native constraint maps and prefix-colliding IDs have
separate regressions. Generators preserve ordinary attempt budgets. The RFD4
runtime fixture also crosses its 100-replicate wrapper boundary at 101.

| Inputs | Scheduler chunks (size 5) | Shared MPNN / Protenix / one-seed Boltz native calls | Foundry MPNN native calls | Three-seed Boltz native calls |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 1 | 1 | 1 | 3 |
| 5 | 1 | 1 | 5 | 3 |
| 10 | 2 | 2 | 10 | 6 |
| 25 | 5 | 5 | 25 | 15 |
| 50 | 10 | 10 | 50 | 30 |
| 100 | 20 | 20 | 100 | 60 |

The same sequence/prediction launch matrix was repeated using copies of the three
retained PD-L1 candidate structures. Larger collections cycle those **three real
structures** under distinct deterministic fixture IDs; they are not 100 newly
generated scientific candidates. Fixtures echo sequence/structure outputs to test
mapping and process grouping; their scores/results are not new scientific evidence.
The earlier retained replay remains the real scientific prediction evidence. GPU
allocation metadata in these fixtures is simulated, not a real GPU reservation.

Machine-readable launch observations, full tile inventory, source hashes,
read-only forecasts and retained-source hash checks are committed in
[evidence/phase7a-native-execution.json](evidence/phase7a-native-execution.json).
Reproduce from the repository using the existing Python environment:

```bash
python scripts/verify_native_collection_execution.py \
  --image-source /path/to/copied-image-api-sources \
  --evidence-dir /tmp/phase7a-native-observations \
  --output /tmp/phase7a-native-verification.json
# Add --retained for the local retained PD-L1 coordinate fixtures.
python -m pytest tests/test_inverse_folding.py tests/test_candidate_peer_transfer.py -q
```

The copied-source directory contains `boltz/` and `protenix/`, each copied from
`/opt/mn_batch_api` in its pinned image. Without copied image sources, the default
fixture discovery uses the sibling container-source repository when available;
missing source skips those image-side tests explicitly. CI without these images
can still run Studio adapter, recovery, capability and model-library fixtures.

### Remaining work and campaign boundary

Phase 6B distributed acceptance campaigns still use one verified native attempt
per child. Their validation and Phase 6C.2 durable ledger are unchanged. Enabling
multi-attempt campaign dispatch still needs persisted per-attempt seeds, exact
attempt receipts, partial recovery/retry semantics and measured equivalence. No
ordinary workflow batching setting authorizes that change.

Separate ordinary-workflow work remains for the dedicated RFD4 inverse-folding
CandidateSet stage, verified Foundry multi-structure execution, BoltzGen inverse
model reuse, generated-Caliby-ensemble lifecycle accounting, and deeper verification
of other selectable cofolders. Those are listed as limitations rather than being
silently certified by the generalized orchestration layer.

Follow-up verification: full regression suite **1,861 passed, 33 skipped,
17 warnings (436.05 s)**. Final focused adapter/MPNN/inverse-library/controller/peer
suite **205 passed (58.18 s)**; subsequent forecast/inventory checks **2 passed**,
managed-path interruption recovery **4 passed**, and cache/sequence-routing checks
**2 passed**. Retained PD-L1 native-launch matrix **30 passed (383.64 s)**.
The skips are opt-in GPU/container tests and unavailable optional dependencies;
they do not establish live inference equivalence or throughput. The native
matrix observations include **61 deterministic cases** and **30 retained-coordinate
cases**, plus the library/recovery/regression checks described above.

Follow-up implementation commit: **`642f986` — `feat: verify native collection
batching and preserve completed members`**. This includes the native contract and
process receipts, shared MPNN grouping, generic inverse-folding collection handling,
partial recovery, read-only inventory/forecast, native launch fixtures and verification
script. The following documentation commit contains the updated report and execution
evidence. Pre-existing `run_app.py` and `tool_references.py` checksums/diff remain
unchanged and are excluded from both commits.
