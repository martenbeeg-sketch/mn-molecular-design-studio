# Phase 7A.2 — Three scientific workflow modes

28 October 2026. Starting HEAD `5c4da02`, branch
`fix/scientific-pipeline-validation`. This change extends Phase 7A's existing
Studio graph and candidate collection controller. It adds no scheduler, editor,
acceptance campaign or scientific scoring framework.

## Mode definitions and graphs

Open **Pipelines → Binder design scientific workflow modes**. Select the mode,
one or more generator operations and, where applicable, independent co-folders.
**Create scientific workflow pipeline** creates a draft; it submits no work.
The same modes are registered starter-template IDs:

| ID | Definition | Graph |
| --- | --- | --- |
| `native_complete` | Complete Native Workflow | Target → complete native workflows → native outcome collection → Native Results |
| `native_comparative` | Complete Native Workflow + Comparative Evaluation | Target → complete native workflows → CandidateSet → shared eligibility filter → independent predictors → unified candidate analyzer → comparative ranking → Results |
| `modular_design` | Modular Design Pipeline | Target → generator-only operations → CandidateSet → conditional initial sequence design → optional redesign → shared filter → independent predictors → unified candidate analyzer → comparative ranking → Results |

For modes 2/3 with several generators, **each generator owns its collection,
sequence/filter stages, predictor branches and analyzer**. Its downstream work can
start when that generator's normalized outputs are ready, while other generators
remain running. A final collection consolidates analyzed variants before the one
comparative ranking and Results node. There is no generator-wide collection
barrier. Canonical predictor evidence namespaces keep the configured measurement
names consistent across generator branches, while unique execution-node IDs,
prediction IDs and source provenance preserve independent histories.

Native workflow internals remain inside the existing native node. Mode 1 appends
no sequence design, external predictor, analyzer or common ranking. Multiple native
engines consolidate for inspection without comparing unrelated scores. Mode 2
adds no sequence design: native sequences, evaluations and filters remain intact.
The initial sequence node in mode 3 submits only compatible backbone-only inputs;
sequence-ready originals bypass it. **Redesign existing sequences · fix interface amino acids** is an explicit optional
setup checkbox. It exposes variants per design and a heavy-atom interface distance
(default 4 Å); mapped target-contacting binder residues and existing protected
positions are fixed through the established sequence-design implementation.
Original amino-acid identities are checked against the structure before redesign;
a mismatch fails rather than allowing unrestricted changes.

Optional SolubleMPNN starts disabled. When
explicitly enabled, it redesigns prepared sequences, preserving originals and
creating independently identified child variants. A retained original backbone is
not designed a second time by the optional soluble stage. General explicit
redesign remains available through existing sequence-design forms.

Boltz-2 and Protenix V2 are the initially reviewed independent co-folder choices.
The first selected branch starts enabled; additional branches require explicit
opt-in under **Branch and processing settings**. Each branch retains its own
sampling, seeds, MSA/template choices and resource settings. For a comparison,
configure the same prediction protocol across the candidates being compared.
The setup panel provides each selected co-folder's sampling/seed controls once
and applies them to every generator branch; new presets use explicit seed 42 by
default. Other prediction settings remain editable through the existing node forms.

## Reviewed ordinary tool capabilities

These are registered ordinary workflow operations, not assertions of newly tested
GPU inference or distributed acceptance-campaign eligibility.

| Complete native operation | Internal scientific stages / limitation |
| --- | --- |
| BindCraft2 `bindcraft2` | Native gradient design, MPNN, validation and filters; recorded ranked-table order where available |
| BoltzGen `boltzgen` | Generation, inverse folding, refolding, filtering and final selection |
| PXDesign `pxdesign_pipeline` | Native sequence design and AF2-IG; extended preset includes Protenix |
| Genie3 `genie3_pipeline` | Native inverse folding and structure evaluation |
| Protpardelle `protpardelle_pipeline` | Native MPNN and binder ESMFold evaluation; no implied independent interface confidence |
| Proteina-Complexa `complexa` | Native generation, filtering and evaluation; sequence/evaluation selection can remain unresolved |

The existing hidden saved-workflow forms for PXDesign, Genie3 and Protpardelle
are explicitly selectable here. Their generator-only tiles are separate operations;
no global catalog visibility or historical graph semantics are changed.

Modular operations include RFD4 guided co-design, RFdiffusion, RFdiffusion3,
Genie3/PXDesign/Protpardelle generator tiles, BindCraft2 Generator, BoltzGen
Generator, Proteina-Complexa Generator, ESMFold2 Binder Generator and RFdiffusion
β-pairing. Only registered operations accepting a shared structure target are
listed. Specialist conditioning, scaffolds and chain mappings still pass through
existing tool forms and execution preflight. Available operation IDs and limitations
are exposed in catalog `scientific_workflow` metadata.

RFD4 supplies sequence/complex co-design; native refolding depends on task options.
It is not offered as a complete evaluated native workflow, even when an output has
a sequence. Explicit native refolding options remain respected. RFdiffusion,
RFdiffusion3, Genie3, PXDesign and Protpardelle generation require initial sequence
assignment. Coordinate residue letters alone do not establish sequence readiness.

## Scientific evidence and ranking

Native measurements remain in original candidate metrics and an explicitly named
original-source measurement snapshot. They remain separate from qualified
`prediction/<branch>/<prediction-id>/<metric>` and `analysis/...` measurements.
One unambiguous prediction can expose `prediction/<branch>/<metric>`; several
samples do not silently produce a mean or best-seed score.

Native Results show the recorded native order only when a finite native rank
exists. It is scoped to the original engine/run, called **recorded native order**,
and does not claim new validation. Otherwise ranking is explicitly **unavailable**.
A generic score or row position is never turned into a new common native rank.

Comparative ranking starts with no implicit scientific scoring rule. Configure
ordered measurements in **Ranking and diversity**, for example
`prediction/predictor-0/iptm`, direction `higher`. Existing protocol cohorting,
finite-value checks, missing values, native/shared filter outcomes and exact
sequence diversity remain authoritative. Missing or incompatible evidence is not
assigned zero or promoted to a rank. Child variants do not acquire parent
prediction scores as new measurements. Original-source measurements are labeled
as original history, including when displayed for child variants.

## Distributed execution and batching

The new presets use automatic placement through the existing trusted-workstation
submission database and shared compute scheduler. A generator is submitted as
**one whole ordinary invocation**, using an immutable manifest with
`unit_kind=native_invocation`. Its native multi-design/trajectory settings, seed
and GPU microbatch controls are not partitioned or rewritten. This unit means one
job, not one scientific attempt or one accepted design. Existing process receipts,
partial recovery, transfer checksums, idempotency keys and cancellation apply.
RFD4's preset requests eight CPU cores and 16 GiB VRAM; resource controls remain
editable. Enabled acceptance campaigns retain their established controller and
one-attempt dispatch limits; this phase enables no multi-attempt campaign packing.

Independent generators are queued without waiting for one another. Prediction
chunks independently use automatic workstation placement, five eligible candidates
per chunk and at most two in flight **per predictor node** by default. Two enabled
predictor branches can therefore expose four simultaneous independent jobs when
there is enough input. Several generators can expose additional ready work.
The existing scheduler decides eligibility and physical leases using exact
image/reference compatibility, GPU architecture, memory and configured limits;
no controller pins jobs to four device numbers or reserves idle GPUs speculatively.
Offline or incompatible devices are not promised useful work. Automatic placement
falls back to the existing local scheduler when peer mode is disabled.

The existing pipeline monitor shows scheduler/heartbeat GPU snapshots, active or
reserved devices, available devices, stale/offline states, queued/pending stage
and chunk counts, worker/GPU allocations and per-stage counts. Queued work has no
GPU allocation until leased. Snapshot availability is separate from per-model
eligibility. Refresh uses the existing run-status action.

Native design counts, native GPU microbatching, candidate chunk size, maximum
in-flight jobs, predictor seeds/samples, returned outputs and scientific acceptance
counts remain separate. There is no universal batch-size control. Successful
branches/chunks survive a failed sibling through existing explicit partial pools.

## Results, export and compatibility

The builder reports the starting preset separately from actual enabled operations
and downstream stages. Results include the effective mode, generator operations,
candidate counts, per-stage counts and an execution plan with actual job/output
counts alongside forecasts. Editing a preset can yield a custom graph; the preset
does not restrict editing. Empty native outcomes also render the mode and native
rejection counts in the normal Results panel.

Normal save/version/reopen, validation, run, monitoring and ZIP export are reused.
Exports add `workflow_summary.json` and `native_sources.json`; the latter retains
source provenance even with zero candidate rows. Native exports include available
structures and `native_sequences.fasta` without calling those designs selected by
a new comparative filter. Native artifact absence remains missing evidence. The
native source runs retain their original raw artifacts; source references and
provenance are preserved, rather than recursively copying complete native workdirs
and model caches into the consolidated export. Comparative exports keep the
existing selected sequence/structure and prediction evidence files.

The graph schema is unchanged. `workflow_mode` and predictor `evidence_namespace`
are additive authored metadata preserved through canvas edits and serialization.
Saved graphs without modes, the original `multi_generator_binder` constructor,
legacy templates, branch controls, cache reuse and standalone tool contracts remain
supported. Namespace changes participate in scientific cache keys. No deployment
or service restart is required for source verification, and none was performed.

## PD-L1 verification

All new target-dependent integration tests use the established PD-L1 chain A with
115 residues. The committed CA-only routing fixture derives from original P4
coordinates and its 75-residue binder; provenance and the original structure hash
are in `tests/fixtures/pdl1_p4_ca.json`. This reduced deterministic fixture is not
scientific evaluation evidence. Synthetic native outcomes, MPNN children and
predictor workers are explicitly fixtures, including hypothetical positive native
outcomes absent from the real BindCraft2 benchmark.

`tests/test_workflow_modes.py` covers all presets, saved round trips, actual normal
graph execution, empty/single/multiple candidates, native rejection provenance,
zero extra native-mode predictors, conditional sequence assignment, explicit
redesign ancestry, missing/failed predictor evidence, failed generators, Streamlit
setup/Results rendering, normal catalog/save/list/validate API routes, independent
generator dispatch, generator/prediction overlap, four independent in-flight
chunks, stale four-GPU inventory and whole-invocation parameter preservation.
Existing Phase 7A batching, peer-transfer, graph, sequence and native-adapter
regressions are included in the recorded validation below.

[Retained replay evidence](evidence/phase7a2-retained-modes.json) records three
completed runs with **zero new designs, predictions or scheduler submissions**:

- Mode 1: original native BindCraft2 benchmark, zero normalized candidates and
  two rejected trajectories retained as provenance; no additional predictors.
- Mode 2: the same real complete native outcome passes through CandidateSet,
  empty predictor stages and analysis automatically. Positive complete-workflow
  comparison branches use deterministic fixtures; no accepted BindCraft2 result
  is fabricated.
- Mode 3: original P4/P5 RFD4 outputs and retained same-sequence Boltz-2 evidence,
  with lengths 75/116/94 and iPTM 0.42641/0.46578/0.65414. Protenix remains disabled
  in this retained replay and is tested using deterministic workers separately.

Original P4/P5 candidate collections, ledgers and benchmark target files match
all previously recorded SHA-256 values. Replay artifacts live under
`/mnt/data/RESULTS/mn-molecular-design-studio/verification/phase7a2-workflow-modes/`.

## Known limits and remaining GPU acceptance

No new biological design, model inference, four-GPU throughput measurement or
cross-workstation deployment was performed. Four-GPU queue exposure and overlap
are deterministic controller tests; they are not observed four-device utilization.
Existing installed runtime compatibility still governs real submissions. Native
ranking is unavailable where the adapter does not record it, and recorded native
order alone does not establish common scientific comparability.

A native invocation publishes its downstream collection after normalization.
Unpublished raw intermediate PDB/CIF files are not treated as valid independent
prediction inputs. Overlap currently uses completed normalized generator branches
while other generators run. Live per-design streaming inside a still-running
native invocation requires an engine-specific durable publication contract and is
not claimed here. Configured native invocation sizes are never reduced merely to
produce earlier downstream work.

Remaining acceptance work: run authorized PD-L1 complete workflows and co-folding
on eligible installed images across the four physical GPUs; observe leases,
concurrent generation/prediction and peer queue placement; exercise real remote
cancellation/partial recovery and Protenix V2 inference; verify positive native
complete-workflow comparison outputs and their exact scientific protocols. These
are separately authorized GPU acceptance runs, not requirements to fabricate new
scientific evidence during preset/routing verification.

## Validation record

All commands used the existing
`/home/user/mambaforge/envs/mn-molecular-design-studio/bin/python` environment.

- **346 passed in 203.25 s**: graph, Phase 7A collection, scientific task,
  candidate execution/peer transfer, inverse folding, execution planning,
  pipeline hierarchy, native generator and BoltzGen normalization regressions.
  Two existing NumPy warnings came from reduced-coordinate inverse-folding fixtures.
- **76 passed in 126.35 s**: new workflow-mode integration tests plus the existing
  MPNN interface-preservation tests. This includes the optional interface-fixed
  redesign setup and request controls, failed generators and four in-flight jobs.
- **46 passed in 125.42 s**: final candidate-controller, peer-transfer and
  CandidateSet regression repeat after native identity/provenance adjustments.
- Frontend `npm run build` completed. Vite reported the upstream ReactFlow
  `use client` directive warning; the packaged editor bundle was rebuilt.
- Changed Python source parsed, and `git diff --check` passed.
- Retained replay completed all three modes; previous scientific source hashes
  matched before and after. New test workers performed no scientific inference.

The focused repeats overlap the larger sets and are not extra distinct tests.
Final targeted native-identity/concurrency checks passed **10 tests in 64.48 s**;
final setup/serialization/UI checks passed **7 tests in 38.55 s**.
Unrelated pre-existing `run_app.py` paper-link edits and untracked
`tool_references.py` remain outside the commit; only this phase's Results hook is
included from `run_app.py`.
