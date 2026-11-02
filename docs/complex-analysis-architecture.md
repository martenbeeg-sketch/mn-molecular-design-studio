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
  input constraints. Database UI result display currently has a dedicated
  Biotite contact-map renderer; the other evaluators reuse the shared saved
  metrics/progress renderer.

## Studio and BinderScope capability matrix

| Capability | Studio analysis already available | BinderScope capability | Decision and planned integration |
| --- | --- | --- | --- |
| Structure parsing, components, chain IDs, provenance | `molecule_imports.parse_structure`, Database representations, Results artifacts, Mol* viewer, and candidate/evaluation lineage. | `structure.py` parses structures and assigns component/chain roles. | **Already reusable.** Keep Studio parser, access control, viewer, and provenance contracts; do not copy the BinderScope parser. |
| Heavy-atom contacts and residue maps | Biotite Interface computes contacts and residue-pair maps for declared protein chain groups. Results artifacts preserve source residue identifiers. | `interface.py` generates interface contacts/maps and 3D review layers. | **Already reusable for protein–protein.** Complex Analysis launches Biotite, renders its existing contact-map artifact, and can highlight the saved interface, binder side, or target side in Mol* using source author chain/residue IDs. Insertion-code selectors are not yet supported by the viewer selector syntax. |
| Buried solvent-accessible area (SASA) | Rosetta Interface Analysis reports `dSASA`, interface SASA fraction, interface energy, and related measures through `analysis.rosetta`. | `interface.py` reports interface buried area. | **Already implemented and now contextually invocable.** Reuse the existing Rosetta workflow and metrics; retain its Rosetta method and units. |
| Hydrogen bonds | Biotite Interface counts geometric hydrogen bonds when explicit H atoms are present. Rosetta Interface Analysis reports its own interface H-bond metrics. | `interface.py` reports interface hydrogen bonds and geometry. | **Already partially implemented.** Reuse the two existing evaluators with their method-specific definitions and explicit missing-input states. Do not merge their counts. |
| Unsatisfied polar atoms / interface H-bonds | Rosetta Interface Analysis reports interface unsatisfied H-bond metrics. Biotite does not calculate unsatisfied polar atoms. | `interface.py` reports unsatisfied polar atoms. | **Already implemented and now contextually invocable through Rosetta.** Identify as Rosetta output, not a universal polar-satisfaction score. |
| Salt bridges and contact geometry | Biotite Interface reports charged-side-chain pairs within its declared 4 Å geometric cutoff. It does not infer energetic strength. | `interface.py` includes polar geometry and salt-bridge context. | **Already reusable.** Present as geometric contacts with chain/residue numbering and cutoff provenance. |
| Hydrophobic surface patches | PEP-Patch is an existing workflow for selected protein chains with hydrophobic and pH-specific electrostatic surface patches, per-residue artifacts, areas, and fractions. | `interface.py` and `liabilities.py` report surface patches/chemical properties. | **Already implemented and contextually invocable.** Reuse saved patch artifacts and declared pH; patch area is not affinity or specificity. |
| Full-length target context and residue mapping | Studio retains target/reference IDs and source residue numbering. ipSAE can use saved structures, PAE, and design-reference evidence when available; this is not a general full-length target-frame mapper. | `target_frame.py` aligns cropped targets to optional full-length references and maps contacted residues. | **Partially present, requiring extension.** First reuse explicit references and existing ipSAE evidence. Plan a separate reusable context-mapping workflow for full-length placement, target residue mapping, and missing-reference reporting. Do not assume residue numbers match across constructs. |
| Glycans and construct accessibility | Some prediction/design inputs can carry glycan-related structure context, but Studio has no general glycan-occlusion evaluation in Complex Analysis. | `target_frame.py` checks glycan context and clashes. | **Missing and worth integrating.** Add only as an explicit reference-context workflow with supported glycan atoms/links and a transparent unavailable state. Do not infer glycan occupancy or biological accessibility from an absent structure. |
| Steric clashes | Existing measures include ipSAE Cα clash counts, selected validation Cα geometry, Rosetta/other workflow outputs, and PoseBusters ligand plausibility checks. They are distinct methods and do not provide a universal full-length all-atom clash test. | `target_frame.py` checks clashes in full-length context. | **Partially implemented.** Connect existing measurements with source/method labels. A full-length all-atom context check is a later extension; do not substitute a Cα threshold for it. |
| Binding-mode clustering | Interface Clustering compares contacted target residues by Jaccard similarity for one identical target structure and author numbering, with a 2,000-item workflow limit. Sequence Clustering and Foldseek are separate existing collection workflows. | `modes.py` and `families.py` compare contacted-residue footprints and families. | **Contextually invocable.** Compare saved complexes submits compatible saved Results/Database structures to the registered `analysis.interface_clustering` workflow. Studio records each source structure, model, chain group, digest, and associated evaluation IDs; the existing workflow computes contacts and Jaccard clusters. |
| Structural comparison across engines/runs | Studio stores distinct model/pose evaluations and supports selected metric comparisons. Foldseek is available for chain-level structure clustering. Existing refolding RMSD metrics are retained with strict residue-list validation. | `modes.py` compares complex poses and binding modes. | **Partially implemented.** Complex Analysis selects authorized saved structures from Results or Database, presents side-by-side Mol* viewers, and can display an overlay from a user-supplied transform plus explicit target-residue map. Studio validates the selected identities/models, residue references, one-to-one mapping, and rigid transform form; it does not fit or score that alignment. |
| Interface specificity | Some prediction/design workflows accept explicit counter-target or specificity panels. Studio has no general saved-pose homolog/paralog specificity evaluation. PEP-Patch is not a counter-target assay. | `specificity.py` evaluates ortholog/paralog panels. | **Partially present at the workflow level; missing as a general analysis.** Reuse existing explicit-panel workflows where compatible, then design a reusable post-prediction panel adapter. Do not claim proteome-wide specificity from one target model. |
| Sequence/structural liabilities | Sequence Composition, Protein-Sol, ESM-1v/ESM-IF1, PEP-Patch, PoseBusters, and some engine-specific developability filters exist. They measure different properties and are not interchangeable. | `liabilities.py` includes sequence and structural liability screens. | **Partially implemented.** Surface existing Studio workflows with their evaluator and provenance. Plan a generic plugin interface for additional well-defined liability checks; no composite pass/fail verdict or invented measurement. |
| Fusion/tag accessibility | Some BindCraft2 modes constrain termini during generation. Studio has no generic post-hoc tag/fusion accessibility analysis for arbitrary structures. | `termini.py` evaluates N/C terminus and fusion/tag accessibility. | **BinderScope-specific behavior requiring generalization.** Do not reuse a generation-time constraint as a post-hoc measurement. Plan a separately named, structure-context workflow that accepts explicit tags/fusions and reports geometric assumptions. |
| Protein–ligand interactions | PLIP and PandaMap provide saved-pose interaction analysis; PoseBusters checks pose/chemistry plausibility. The existing workflows retain pose and source compound provenance. | BinderScope is primarily an antibody/protein-complex analysis tool. | **Already implemented in Studio.** Connect these workflows for the selected ligand pose, reusing their renderers and compound identity mapping. |
| Protein–DNA/RNA interactions | Structure parsing and visualization support DNA/RNA components. Current Biotite Interface requires protein chains; no registered contextual protein–nucleic-acid interaction evaluator is available in this slice. | BinderScope's protein-chain assumptions do not make these analyses applicable. | **Not currently supported by the contextual analysis tools.** Keep structures selectable and visible, state unsupported checks explicitly, and add only scientifically suitable existing or separately reviewed workflow modules later. |
| Source viewer and human review | Mol* and authenticated Database/Results selection are reused; per-user persistent notes remain separate from observations. | BinderScope includes an offline review bundle, interactive family board, per-design inspector, binding-mode/species map, and rule-based query/filter view. | **Partially adapted.** The comparison panel provides a compatible structure board, side-by-side inspector, saved contact-map view, and evaluator filter. Values retain their original evaluation provenance. Human review remains separate from scientific measurements and is not part of filtering or scoring. Campaign family semantics and browser-local review storage are not imported. |

BinderScope's family board, inspector, and binding-mode map informed the
compatible-collection board and side-by-side inspector below. Its family
labels, campaign filters, and browser-local review state are not portable
generic structure metadata.

### Phase 3 comparative workspace (implemented foundation)

The **Compare saved complexes** panel provides two to five independently
selected structure slots. Each slot uses the existing authenticated Results
job/artifact or Database molecule/representation selectors, preserves the
exact model or artifact identity, and lets the user map target and
protein-partner chain groups. Compatibility is checked against Structure A.
The first comparison view is limited to protein–protein structures;
protein–ligand and protein–DNA/RNA structures remain supported by the common
single-complex model and viewer but are not yet compared in this panel.

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
structure-associated Database evaluations, with an evaluator filter. It keeps
one row per existing record and includes measurement and provenance payloads;
it does not normalize methods, rank structures, merge scores, or treat human
review decisions as scientific measurements. Missing saved evaluations are
reported explicitly.

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
or infer residue mappings across numbering schemes. It does not invoke
Foldseek from its collection selection. The supplied-transform view is an
inspection aid and does not produce scientific evaluation measurements.

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
filters, specificity map, full-length target frame, glycan/steric context,
sequence liabilities, and fusion/tag accessibility remain planned extensions.

### Planned implementation phases

1. **Complete Biotite Interface vertical slice** (implemented and verified in
   live execution): exact Results and Database structure/model, chain groups,
   evaluation reuse, workflow submission, measurements/contact map, and
   provenance.
2. **Complete protein–protein contextual connections**: Rosetta Interface
   Analysis, ipSAE, and PEP-Patch are connected to Results and Database
   submission; all three were live-run from Results. Database execution uses
   the shared adapter. Compare saved complexes also launches Interface
   Clustering for collections with identical target coordinates, author chain
   IDs, and numbering. Biotite saved residues can also highlight in Mol*.
3. **Database as a first-class source** (implemented): direct authorized source
   selection, stable structure/model/group IDs, ACL/project validation at API
   boundaries, private staging, scheduler submission, and evaluation
   ingestion onto the same structure representation. Imported-source behavior
   has focused API/adapter tests; the real Database launch used a generated
   structure already saved in Database.
4. Connect existing PLIP, PandaMap, and PoseBusters workflows for selected
   protein–ligand poses.
5. Add reusable full-length context/reference mapping, glycan context, and
   compatible steric-clash checks only where a validated existing workflow or
   separately reviewed module supports them.
6. Add compatible cross-run/model binding-mode comparison, explicit-panel
   specificity, advanced structural inspection, and protein liability
   integrations while keeping predicted, measured, geometric, and empirical
   quantities separate.
7. Add scientifically appropriate protein–DNA/RNA analyses when a registered
   supported workflow exists. Keep unavailable capabilities explicit.
8. Adapt BinderScope's family board, inspector, map, and query interface ideas
   to compatible Studio evaluation collections. Do not carry over campaign
   family semantics, browser-only review storage, or inferred pass/fail
   rankings.

The All Tools catalog remains available for standalone workflows throughout.
Complex Analysis should become their preferred contextual entry point; it must
not create a second scheduler, job store, evaluation store, or copied scientific
algorithm.

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

The local Complex Analysis page was inspected in the browser. The first
on-demand duplicate Mol* viewer collapsed to zero height inside an expander;
that attempt was removed. Interface selection now feeds the existing main
Mol* viewer, preserving the contact map's author chain and residue numbering.
Current Mol* selectors accept numeric author residue IDs; insertion-code
selectors remain visible in contact output but are not highlighted by this
viewer syntax.

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
