# Complex Analysis workspace

## Operational flow

Complex Analysis is the contextual workspace for selecting a molecular
structure, choosing applicable existing analyses, configuring the selected
representation and component mapping, invoking Studio's registered workflow,
and viewing the resulting evaluation. All Tools remains the standalone catalog
and setup surface.

The workspace loads structures from Database or completed Results, uses the
existing PDB/mmCIF component parser and Mol* viewer, and retains per-user
persistent human review notes separately from scientific observations.
Database and Results reads use their current ACLs and stable
record/evaluation or job/candidate/artifact identifiers. Database has an
**Open in Complex Analysis** route for a selected record, evaluation, and
structure representation. Analyses can launch from either source. Saved
Database structures are adapted into private, deterministic input bundles for
registered Studio workflows; this creates no synthetic prediction job.
Evaluations are ingested back onto the exact Database structure using the
existing evaluation contract.

### First execution slice: Biotite Interface

For completed prediction Results with normalized candidates, and for an
authorized Database structure, the workspace:

1. Selects a specific Results job and structure artifact/model.
2. Detects the linked candidate and its recorded model metadata.
3. Lets the user select two disjoint protein chain groups, including multiple
   chains per group. The labels do not imply a binder or experimental binding.
4. Checks existing Biotite Interface analysis jobs for the same source job,
   candidate, structure digest, chain mapping, and contact distance.
5. Reuses an exact completed match, shows an active job, or offers an explicit
   **Run analysis** / **Retry analysis** action.
6. Submits `analysis.biotite_interface` through the existing authenticated
   workflow API and scheduler. Results inputs retain the source candidate and
   job lineage. Database inputs use the selected authorized record, structure,
   model, project, and chain mapping directly; the adapter stages a private
   input bundle and does not manufacture a Results prediction.
7. Shows saved measurements, progress, failure details, contact maps, and the
   persistent Evaluation identifier in this workspace. Results evaluations
   use their existing ingestion path; Database-source evaluations attach to
   the existing molecular record and exact structure representation.

Opening a structure never starts a calculation. Results launches require a
completed accessible job, an associated normalized candidate, a local
accessible PDB/mmCIF artifact, and that source job's project assignment.
Database launches require an authorized record, representation, and matching
project/private scope. The API rechecks those stable IDs and ACLs at validation
and submission. Multi-model coordinate files are treated as separate saved
Database representations when available; the workspace does not split an
unrepresented internal model from one coordinate file into selectable records.

Biotite measurements are structural/geometric observations. Hydrogen bonds
are available only when explicit hydrogens are present; missing hydrogen-bond
inputs remain unavailable. Neither contacts nor model confidence establish
experimental binding.

### Additional contextual protein–protein evaluations

Rosetta Interface Analysis, ipSAE, and PEP-Patch are now available from the
same Results context through their registered workflows. Each card checks for
an accessible matching evaluation before offering a launch and shows the saved
measurements, job state, Evaluation ID, and provenance when available.

- Rosetta Interface Analysis and ipSAE currently accept one selected protein
  chain in each of the two interface groups. Their exact selected model,
  structure digest, groups, source job, and method settings are saved in the
  analysis artifacts. Rosetta SASA and interface terms retain Rosetta's
  definitions. ipSAE uses its existing saved-structure and saved-PAE paths;
  PAE-derived values remain unavailable when a matching PAE input is missing.
- PEP-Patch accepts a selected protein surface chain group and explicit pH.
  It uses the existing QC batch workflow and displays its saved tables and
  residue-patch artifacts. Surface measurements are not affinity or specificity
  estimates.
- These integrations use the existing authenticated workflow submission,
  scheduler, evaluation ingestion, and results renderers. They add no new
  execution service or scientific calculation. All four workflows accept a
  Database source through the shared adapter, subject to each evaluator's own
  input constraints. Saved outputs flow back to the exact selected Database
  representation and are shown in the Scientific Inspector and Evaluation
  Board with method-specific measurements and provenance.

## Studio and BinderScope capability status

| Scientific or inspection capability | Studio status |
| --- | --- |
| Results/Database structures, access control, Mol*, saved Evaluations and provenance | **Supported.** Contextual actions use registered workflows, the existing scheduler and Evaluation ingestion; selected structure/model IDs and digests remain attached to results. |
| Protein–protein interface evaluation | **Supported through existing tools.** Biotite Interface, Rosetta Interface Analysis, ipSAE and PEP-Patch retain their own definitions, units, parameters and missing-input states. Interface Clustering reuses its registered workflow. |
| Protein–ligand evaluation | **Supported through existing tools.** PLIP and PandaMap report interactions; PoseBusters performs its own pose/chemistry validation. A reconstructed docking complex is attached to the existing Database molecular identity with receptor/pose digests, atom mapping, ligand isomeric SMILES and source provenance. |
| Scientific Inspector | **Implemented.** Mol* and saved evaluations are linked through author-residue selectors, including insertion codes. Interface, Confidence, Surface/Chemistry, Findings and optional Reference Context panels appear when their inputs exist. Evidence rows can select their mapped residues in Mol*. No universal quality score is calculated. |
| Collection Map and Evaluation Board | **Implemented.** Authorized Database structures can be paged and selected beyond the 2–5 Compare panel limit. Compatible protein–protein maps show saved Biotite contact footprints and per-structure frequencies; protein–ligand rows preserve pose and ligand identity without applying protein–protein metrics. Saved Interface Clustering and Foldseek assignments remain distinct. Missing evaluations are not zero-contact results. |
| Contact measurements | **Supported with distinct meanings.** Saved residue-pair rows are pair counts; each contacted residue is counted once per structure; contact frequency is the number of evaluated structures contacting that residue divided by the number with compatible saved maps. |
| Optional reference association | **Implemented as component-specific context.** It preserves source/reference structure IDs, digests, model IDs, chains and author numbering. A unique exact subsequence correspondence covering the full source construct is valid even when it covers only part of the reference; ambiguous or non-identical correspondence is not enabled for reference-dependent evaluation. |
| Paired reference visualization | **Supported when an authorized association exists.** Before placement is evaluated, the selected complex and reference display in separate original coordinate frames. A matching completed placement Evaluation enables the transformed source complex beside the original reference in one Mol* viewer. |
| Exact-correspondence structural placement | **Supported as an optional registered Evaluation.** A unique contiguous exact source subsequence is fitted to the explicitly associated reference chain using mapped Cα atoms and the shared Bio.PDB.Superimposer Kabsch fitter. Results retain both structure identities, digests, models, chains, residue map, transformation, aligned-pair count, source/reference coverage, RMSD, method/version and parameters. |
| Non-identical correspondence and placement | **No persisted Evaluation.** Mol* has a separate transient TM-align superposition for selected polymer chains; it does not save a Database-bound correspondence, placement Evaluation or provenance. No Studio US-align workflow is registered. |
| Full-length spatial context, glycans and advanced steric context | **Unavailable as validated Complex Analysis evaluations.** No full-length clash/occlusion, glycan accessibility, or context score is claimed. These require separately validated calculations and provenance. |
| Protein–DNA/RNA interactions, generalized specificity and fusion/tag accessibility | **Unavailable as general contextual analyses.** Existing prediction-time constraints or protein-only tools are not presented as post-hoc evaluations for these cases. |

BinderScope's Inspector, map, board and query concepts informed these views.
Studio retains its Database identities, ACLs, Mol* viewer, registered
evaluators, and per-user human review annotations; campaign-specific labels,
browser-local review state, and inferred rankings are not imported.

### Optional reference-frame placement

`analysis.reference_placement` uses the saved, user-scoped component reference
association and the authorized Database source adapter. Both structures must be
the exact digest-current, single-model representations recorded by that
association. The source component must map uniquely and exactly to a contiguous
subsequence of the selected reference chain. The evaluator checks each mapped
author residue number, insertion code and amino-acid identity, then fits the
corresponding Cα atoms with `Bio.PDB.Superimposer`; the same fitter is used by
the existing target-aligned benchmark path. At least three non-collinear pairs
with finite coordinates are required. Alternate Cα atoms use a unique highest
occupancy; tied occupancy is rejected.

The saved transformation maps selected-complex coordinates into the associated
reference frame, using a row-major column-vector matrix for Mol*. The exact
source/reference IDs, groups, models, chains, digests, mapping, method/version,
parameters, aligned-pair count, coverage and RMSD are stored in the normal
analysis job and Database Evaluation. RMSD has no acceptance cutoff. Mol* applies
the matrix to the source complex only in the viewer; saved coordinates and
existing evaluations are unchanged. If placement is absent, stale, inaccessible
or invalid, the two structures remain visible in separate original frames.

A good fit does not establish experimental binding or biological accessibility.
The evaluator does not calculate clashes, glycans, occlusion, or accessibility,
and a partial reference chain does not become full-length context by inference.
TM-align/US-align and validated non-identical structural correspondence remain
separate future work.

Mol*'s existing TM-align control is useful for interactive structural
comparison of selected polymer chains, including non-identical structures. Its
viewer transform is not the exact-reference workflow above: it does not
establish an explicit Database reference association or persist the aligned
residue correspondence, source/reference digests, method parameters and
transformation as a Studio Evaluation. It therefore remains a separate visual
comparison capability until an Evaluation-backed workflow with validated
inputs and provenance is implemented. No US-align integration was found.

The documented validation fixture uses a five-residue source construct and its
unique exact subsequence in a seven-residue reference chain. A known rigid
rotation and translation are recovered with five aligned pairs, source coverage
1.0, reference coverage 5/7 and near-zero Cα RMSD. The Mol* adapter test checks
that the original reference is loaded first and the source complex receives
the saved source-to-reference matrix. The local Database contained no saved
reference associations during this verification, so no genuine biological
association was available for a live paired-view demonstration.

BinderScope's family board, inspector, and binding-mode map informed the
compatible-collection board and side-by-side inspector below. Its family
labels, campaign filters, and browser-local review state are not portable
generic structure metadata.

### Compare saved complexes and Collection Map

**Compare saved complexes** keeps its focused two-to-five structure view. It
accepts authorized Results and Database structures, including compatible
protein–protein and protein–ligand representations, and preserves their exact
model, pose and component identities. The shared protein–protein map requires
matching selected target sequence and author numbering. Protein–ligand
comparison keeps ligand identity and saved pose measurements separate from
protein–protein contact metrics.

The **Collection Map & Evaluation Board** is the larger Database-centered
workspace. Its ACL-scoped structure discovery is paged and returns lightweight
record/representation summaries; evaluation observations and provenance are
loaded only for selected records. The map supports more than five selected
structures (up to 1,000 per interactive view), with row and residue
pagination. Compatible protein–protein rows use matching component sequence
and residue correspondence. Protein–ligand collections are handled separately
by pose and ligand identity.

For protein–protein collections, contact-frequency values are derived from
saved Biotite contact artifacts. Pair rows remain residue-pair observations;
unique contacted residues are deduplicated within each structure; frequency
counts how many evaluated structures contact each reference residue. A missing
artifact is distinct from an evaluated map with zero contacts. Saved Interface
Clustering assignments and representative structures are shown by run and
cluster. Foldseek assignments are labeled as structural similarity, not
binding-mode similarity. The Evaluation Board filters saved method-specific
measurements and provenance without merging scientific findings, prediction
confidence or human review into a score.

Before showing a shared target-residue map, the panel requires identical
selected target chain-group sequences and author residue numbering, including
insertion codes. Different saved target labels are surfaced as a note when the
sequence and numbering agree. Missing numbering, mismatched targets, missing
groups, and selecting the same representation twice block shared comparison.
The panel does not perform sequence alignment or infer residue correspondence.

Structures use separate Mol* viewers by default. The shared-frame overlay is
available when the source files have identical coordinates or both carry the
same explicit coordinate-frame ID. For different frames, an operator can
provide a row-major 4×4 moving-to-reference rigid transform and an explicit
target-residue correspondence in the **Overlay with supplied transform**
panel. The record is bound to the selected stable structure ID, model ID, and
structure digest;
residue keys must exist in the selected target chains and must be one-to-one.
The UI reports the transform source and method, mapped-pair count, omitted
selected target residues, and any user-supplied mapping note. The viewer applies
the transform to the moving structure in memory; source coordinates and saved
evaluations remain unchanged. Transform form validation rejects non-finite,
non-rigid, or reflective matrices. It does not verify fit accuracy, calculate
an alignment RMSD, infer missing residue mappings, or estimate mapping
uncertainty. Multi-model coordinate files are not transformed unless the source
selection is a single-model representation.

An **External JSON record** input mode now accepts the versioned
`mnds.external-alignment/v1` envelope. The downloadable template is bound to
the currently selected structure identities, model IDs, and coordinate-file
digests. Imported records must declare the moving-to-reference direction,
transform method, source tool/version, and explicit target-residue mapping;
Studio validates the envelope against current selections before applying the
matrix. The displayed source/version/record ID is user supplied and is not
independently authenticated. The record is not stored as an evaluation, and
Studio still does not derive or score an alignment.

The existing refolding target-aligned RMSD path now requires the same ordered
chain IDs and exact author residue numbers, insertion codes, and residue names
on both sides. It reports the measurement as unavailable when those lists do
not match; it no longer truncates unequal arrays or treats a positional prefix
as a residue correspondence. It does not remap numbering schemes.

When every selected structure has a matching completed Biotite evaluation
for the exact model, digest, chain groups, and 5 Å contact setting, the panel
shows a residue-presence map transcribed from the saved contact artifacts.
If target chain identifiers do not match, the existing author-numbered
contact lists are shown separately. The saved-map view itself transcribes
existing contact artifacts. The registered Interface Clustering launch below
is enabled only when target atom coordinates, author chain IDs, sequence, and
numbering match; its outputs remain the existing workflow's measurements.
Foldseek continues to be reached through its established collection tool.

#### Contextual Interface Clustering

For a compatible selection, Compare saved complexes exposes the existing
`analysis.interface_clustering` workflow and its Jaccard threshold and contact
distance parameters. The API rechecks access to every selected Results or
Database representation, then stages a private candidate JSONL bundle with
stable source identities, model IDs, structure digests, chain groups, and
associated evaluation IDs. The registered workflow, scheduler, Results
ingestion, and Evaluation display handle execution and outputs; the contextual
adapter does not calculate contacts or cluster assignments. This entry point
keeps all sources in one project access boundary. The workflow retains its
2,000-item ceiling, while the comparison panel currently selects two to five.

The saved-evaluation table includes Results candidate metrics and
structure-associated Database evaluations, with method and available
provenance filters. It keeps one row per existing record and preserves native
measurement definitions; it does not normalize methods, rank structures,
merge scores, or treat human review decisions as scientific measurements.
Missing saved evaluations are reported explicitly.

#### Live comparison verification (2026-11-02)

The comparison panel was opened against two existing completed Results
structures from separate runs and engines:

- Boltz-2: `refolding-validation/20261029-094704-dcc0cd11`, model
  `i189ea24a074b59c6_model_0.pdb`.
- Protenix v2: `refolding-validation/20261029-094613-7f418205`, model
  `idc0c4d94b21517ef_sample_0.cif`.

The UI showed the exact selected artifacts, target chain A and partner chain B,
a matching 115-residue target sequence and author residue numbering, and both
Mol* viewers side by side. Existing candidate evaluation rows from both
engines remained separate and displayed their own run, candidate, stored
measurements, and artifact provenance. The Boltz-2 sample had an exact saved
Biotite contact artifact, which the UI displayed as author-numbered target
residues. No matching Biotite map was saved for the selected Protenix v2
representation, and the UI stated that limitation. The structures had no
shared coordinate-frame metadata or identical coordinate digest, so no
overlay was offered. This was a read-only comparison; no analysis or
prediction workflow was submitted.

### Remaining comparative extensions

The current panel does not calculate structural superpositions, verify the
scientific fit of a supplied transform, align non-identical target constructs,
or infer residue mappings across numbering schemes. Compare saved complexes
can launch the registered Foldseek collection workflow for compatible
single-chain protein inputs; the Collection Map displays saved Foldseek
assignments but does not rerun that calculation. The supplied-transform view
is an inspection aid and does not produce scientific evaluation measurements.

#### Structural integrity and supplied-transform inspection (2026-11-02)

The existing refolding pose-RMSD helper no longer pairs unequal Cα arrays by
truncating them. It emits an unavailable status when chain lists, residue
counts, author residue numbers, insertion codes, or residue names do not match
in the same order. It does not infer a residue map. Tests cover matching,
missing, renumbered, and residue-name-mismatched inputs, plus binder-list
mismatches.

The Mol* comparison UI accepts a user-supplied rigid transform and target
correspondence, validates their binding to the selected structure/model IDs and
coordinate-file digests,
checks that each mapped residue exists in the selected target group, and
rejects duplicate mappings and non-rigid matrices. A partial map is reported
with counts of unpaired selected target residues. Transform provenance and a
user-entered mapping note are shown with the overlay. The viewer applies the
matrix to the whole moving structure for display only. Studio does not derive
or fit the transform, verify its scientific fit, calculate alignment RMSD, or
estimate mapping uncertainty. Transformation-source and mapping-note text are
user supplied and are not independently authenticated.

Focused tests cover the validation path and a Streamlit UI test confirms that
Mol* receives the validated transform. The running UI was inspected with the
saved Boltz-2/Protenix pair; both remain in their original frames until an
external transform and explicit map are supplied. No real aligned overlay was
claimed in this verification.
BinderScope's campaign-specific family clustering, embedding, rejection
filters, specificity map, validated full-length spatial context, glycan/steric
context, sequence liabilities, and fusion/tag accessibility remain planned
extensions. The optional Database reference association remains independent
from binder-design provenance and does not claim reference identity unless a
user or source workflow explicitly records it.

### Delivery status and remaining work

The Results/Database contextual workflows, optional component reference
associations, Scientific Inspector, protein–protein and protein–ligand
comparison paths, Collection Map, Evaluation Board, and saved clustering
display are implemented. They reuse Studio's registered tools, scheduler,
Database ACLs, Evaluation ingestion and Mol* viewer. Protein–ligand docking
representations remain attached to their original Database molecular record
and retain source pose, receptor, ligand chemistry, coordinate and mapping
provenance.

Remaining scientific work is limited to separately validated calculations:

1. Derive and validate structural correspondence/fit across non-identical
   constructs or numbering schemes, including any TM-align/US-align extension;
   verify supplied transforms against structures rather than treating their
   declared provenance as authenticated.
2. Add validated spatial context measurements for justified full-length
   references before exposing them as evaluations; exact-component placement
   alone does not provide these measurements.
3. Add advanced interface context such as all-atom full-length clashes,
   glycan occlusion, and construct/tag accessibility only through suitable
   registered evaluators with explicit inputs and provenance.
4. Add general protein–DNA/RNA interaction and counter-target specificity
   analyses only when their scientific methods and input requirements are
   defined and validated.

BinderScope campaign-specific family semantics, browser-local review state,
and composite rankings are outside this workspace. Predicted confidence,
scientific findings and human review remain separate throughout.

The All Tools catalog remains available for standalone workflow access.
Complex Analysis is an additional contextual entry point and does not create a
second scheduler, job store, Evaluation store, or scientific implementation.

## Comparison rules and limits

One saved observation remains one model/pose/calculation. The workspace does
not pool distinct poses or models. A reusable evaluation match is scoped to its
source collection, selected structure digests, chain groups, method, and
scientific parameters. Interface Clustering additionally requires identical
target heavy-atom coordinates, author chain IDs, residue numbering, and at most
2,000 items. Missing structures,
references, confidence values, or unsupported checks are displayed as
unavailable rather than guessed. Geometric and predicted evidence is not an
experimental binding result.

## BinderScope source and licensing

The inspected BinderScope source is MIT licensed. No BinderScope code was
copied or vendored for this implementation, so no BinderScope source notice was
modified or removed. Future source reuse must retain its copyright and MIT
permission notice and cite the upstream source. BinderScope's MIT license does
not grant separate rights to BindCraft 2; its campaign launcher and scheduler
integration are excluded from Studio.

## Verification

### Database collection-discovery profile (2026-11-03)

The local Database contained 21,134 molecules, 30,265 structure
representations and 38,093 observations. No evaluation projections were stale.
The detailed `list_molecules(..., limit=1000)` call took 13.38 seconds before
the discovery-path change; evaluation projection rebuilding accounted for
about 0.16 seconds. The profiled work included full-list selection and exact
count queries, plus loading and normalizing 1,543 observations and 1,530
structure rows for the first 1,000 molecules.

Complex Analysis structure pickers now use the ACL-filtered,
`created_at`-indexed `list_structure_records` page query. It returns only
record/component and minimal representation summaries, supports `offset` /
`next_offset`, and does not rebuild or hydrate evaluations. On the same local
Database, measured 1,000-record pages took 0.89, 0.82 and 2.64 seconds at
offsets 0, 1,000 and 2,000 (median 0.89 seconds; cache/page variation is
visible). Loading the full saved record after selecting one entry took 0.24
seconds and retained its two saved observations. The general detailed
`list_molecules` API remains available where callers need complete evaluation
history; this optimization narrows the picker query rather than changing
evaluation semantics.

The focused consolidation regression selection passed 117 tests, including
Complex Analysis, docking-to-complex, Database listing and the relevant API
ACL/submission cases. Mol* residue-selection tests passed 3/3 and TypeScript
type-checking passed. In the 2026-11-03 maintenance follow-up, the ipSAE API
assertion was traced to a stale mocked CSV: the Studio parser's 71-column
schema exactly matches the registered ipSAE runner, while the mock omitted
seven target-aligned pose columns and shifted later values. The fixture now
matches the runner field order; no production calculation changed. The
relevant Complex Analysis, Database, workflow, reference-integrity, ipSAE and
extended API selection passed 216 tests; one API test was skipped because the
optional `mcp` Python module is unavailable in this environment.

### Verified by real execution

On 2026-11-02, the Results route launched all four contextual protein–protein
evaluators on the same benign Boltz-2 model from
`refolding-validation/20261029-094704-dcc0cd11`, candidate
`20261029-094704-dcc0cd11_boltz2_01_001`, structure
`i189ea24a074b59c6_model_0.pdb`, with selected chain groups B/A:

- Biotite Interface job `analysis/20261102-173608-96121bb5` completed and
  ingested as Evaluation `evaluation-b0de78f9e962514fb9ba5f478ab012aa`.
  It recorded 373 heavy-atom pairs, 56 residue-pair contacts, 20 interface
  residues on B, 25 on A, and one geometric salt bridge. Hydrogen bonds were
  unavailable because the selected structure has no explicit hydrogens.
- Rosetta Interface Analysis job `analysis/20261102-175724-c0343bd3`
  completed as Evaluation `evaluation-8ee7a48a9a46945cea9d28c8b39b133a`.
  Its saved outputs included Rosetta interface energy -372.2 and surface
  hydrophobicity 0.13.
- PEP-Patch job `analysis/20261102-175854-97e2f675` completed as Evaluation
  `evaluation-bf77663ffe93f3c721b53e46359f7d3f` with its saved patch tables.
- ipSAE job `analysis/20261102-175952-2002ec7f` completed as Evaluation
  `evaluation-f530273afc25c35b45d38f3d40cb1e9e`; saved ipSAE, ipTM, and
  interface pLDDT scores were 0.218506, 0.384524, and 0.376324 respectively.

For the Results Biotite job, the workspace rendered the saved residue contact
map and its source chain/residue identifiers. Reopening the same structure,
model, chain groups, and settings reused the completed Evaluation without
rerunning it. The refreshed Results page also reopened the saved Rosetta,
ipSAE, and PEP-Patch evaluations for that same exact artifact and chain groups;
their measurement tables and Evaluation IDs were visible, with no new jobs
launched. The saved Biotite map showed 56 residue pairs. Mol* highlighting
reported 20 residues on the selected binder side, 25 on the target side, and
45 across both sides.

A direct Database launch was then verified against an existing generated
three-component complex already stored in Database: molecule
`mol-edbc300706a703d1dcb9789eeb2d7201`, structure
`structure-be9ac776347f2d009eb748cf6e2f925d`, exact structure SHA256
`4f9f1fb9d9eff2e4ab543bb54399475a1617ac30f11da718da5bda0289aae2f7`, model
`20261102-175754-f9172f27_localfold_01_001`, and chain groups A versus B,C.
Biotite job `analysis/20261102-181951-5366c951` completed and ingested as
Evaluation `evaluation-bd31e5e18a415be47ba26433cbce7e0e`. Its saved outputs
reported 544 heavy-atom pairs, 65 residue-pair contacts, 27 residues on A, 23
on B/C, and zero geometric salt bridges. Hydrogen bonds were unavailable
because explicit hydrogens were absent. The evaluation points to the exact
existing Database representation and digest; it contains source kind,
molecule/structure/group/evaluation IDs, model, chain groups, method, settings,
and source job metadata without fabricating a Results prediction job. Reopening
the same selection showed the saved completed result instead of a Run/Retry
control.

The local Complex Analysis page was inspected in the browser. Interface
selection feeds the existing main Mol* viewer, preserving author chain and
residue numbering. The Scientific Inspector also supports insertion-code
selectors; focused frontend and Streamlit tests verify those exact selectors.
An earlier duplicate on-demand viewer was removed after it collapsed to zero
height inside an expander.

### Verified by focused tests

The project-environment focused run passed 29 tests: all tests in
`tests/test_complex_analysis.py`, the Database analysis API validation and
submission/ACL test, and the saved pipeline scope ACL test. Coverage includes
imported Database source selection and submission, owner and project boundary
denials, stable structure/model/group IDs, no synthetic prediction job,
Database evaluation ingestion onto the exact representation, exact-result
reuse and access filtering, missing metadata, and selected residue numbering.
Container-backed execution was not part of these mocked/API tests. The three
additional Results evaluations and the Database Biotite run above are real
scheduler executions, separate from the mocked submission tests. The focused
pytest invocation reported 29 passed and no skipped/deselected tests. The
repository's collection hook excludes tests marked `gpu` or `gpu_e2e` unless
`MDS_GPU_TESTS=1`; those marks were not selected here. Ten container tests were
not included because the command targeted the workspace/API path. There was no
GPU generation, deployment, push, or commit.

A live Database launch was verified for a generated structure saved in
Database. The same imported-Database adapter and ACL behavior are covered by
focused tests, but a second real scheduler execution from an imported-only
record was not performed. Database contextual execution is limited to Biotite,
Rosetta Interface Analysis, ipSAE, and PEP-Patch; Rosetta and ipSAE retain their
one-chain-per-group constraints, PEP-Patch needs a protein surface chain, and
PAE-derived ipSAE outputs remain unavailable without a compatible saved PAE.
