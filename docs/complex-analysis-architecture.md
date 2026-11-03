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
| Paired-reference visualization | **Supported when an authorized association and current placement Evaluation exist.** The complete saved reference coordinate file and transformed source share one Mol* view. Author-numbered sequence tracks distinguish mapped/present, resolved-but-omitted, and sequence-declared-but-unresolved reference residues. Other saved reference chains/components remain visible. Studio does not infer that the saved file is biologically full-length. |
| Exact-correspondence structural placement | **Supported as an optional registered Evaluation.** A unique contiguous exact source subsequence is fitted to the explicitly associated reference chain using mapped Cα atoms and the shared Bio.PDB.Superimposer Kabsch fitter. Results retain both structure identities, digests, models, chains, residue map, transformation, aligned-pair count, source/reference coverage, RMSD, method/version and parameters. |
| Non-identical correspondence and placement | **No persisted Evaluation.** Mol* has a separate transient TM-align superposition for selected polymer chains; it does not save a Database-bound correspondence, placement Evaluation or provenance. No Studio US-align workflow is registered. |
| Optional reference-frame steric/environment context | **Supported as a registered, optional Evaluation.** Requires a current explicit component association and matching saved Reference Placement Evaluation. Reports static heavy-atom distances, residue evidence and configured van-der-Waals envelope overlap against resolved omitted reference residues, explicitly selected reference chains and recognized modeled sugar residues. Geometry against the associated source component and other complex components is reported separately. A selected assembly ID is accepted only where the file confirms identity-operator membership; non-identity operators are not expanded. No biological accessibility or binding-quality score is claimed. |
| Tag/fusion terminal geometry and clearance advisory | **Supported as a registered post-hoc Evaluation.** Flexible tags use grid-path reach estimates; folded fusions use spherical exit-cone envelopes. Measured envelope clearance is retained separately from suitability; a folded fusion at an interface-contacting terminal residue gets an explicit advisory even if its sampled sphere clears. Incomplete termini, cyclic/multichain inputs, missing DSSP and unverified B-factors remain explicitly limited. No tag dynamics, biological accessibility or experimental SPR/BLI result is claimed. |
| Full-reference tag/fusion context | **Supported only as a separate optional Evaluation.** It requires an explicit validated component association and a saved exact Reference Placement Evaluation. It checks resolved reference heavy-atom geometry only; it does not calculate biological assemblies, glycan accessibility or experimental binding. |
| Protein–DNA/RNA interactions and generalized specificity | **Unavailable as general contextual analyses.** Existing prediction-time constraints or protein-only tools are not presented as post-hoc evaluations for these cases. |

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
and a partial reference file does not become full-length context by inference.
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

The small unit fixture uses a five-residue source construct and its unique
exact subsequence in a seven-residue reference chain. The real-coordinate
validation fixture uses the complete tracked 5O45 mmCIF file (a partial PD-L1
structural construct, not the full-length biological protein): the explicitly
associated chain-A source is the unique author-numbered 45–105 construct, and
the original chain-B partner is retained. A known rigid transform was applied
to the source atoms before Database registration. The registered placement
workflow recovered 61 pairs, source coverage 1.0, reference coverage 61/129,
and sub-millangstrom Cα RMSD. Independent coordinate checks verified that both
the source target and partner return to the reference coordinates under the
saved Mol* matrix. The records and association live only in an isolated test
Database; this does not create a biological association in production data.

The paired Mol* adapter loads the complete saved reference first, untransformed,
then the original selected complex with the saved source-to-reference matrix.
The reference is grey and the source remains chain-colored. Reference sequence
buttons select exact author-numbered residues in the first Mol* structure.
Mapped source residues, resolved reference-only residues and unresolved
declared sequence positions are distinct. CIF-declared branched entities and
other parsed components are shown as file annotations; they are not treated as
glycan occupancy, accessibility or context measurements. Without a current,
authorized placement, both structures stay in their original coordinate
frames, and ordinary Inspector panels remain available.

### Optional reference structural-context inspection and measurement boundary

This is an optional visualization layer for any explicitly associated Database
reference, not a binder-design feature. It presents the existing source complex
in the selected reference coordinate frame, the entire saved reference file
(including additional chains and resolved non-protein components), and the
selected reference chain's declared sequence/residue map. This shows the full
saved reference representation, not necessarily a full-length biological
protein. Mol* does not synthesize missing reference coordinates. A sequence
tile without coordinates is labelled unresolved; regions absent from the
source construct are not presented as measured occlusion or accessibility.
When reference identity, access, digest, mapping, or placement checks fail, the
ordinary source-complex Inspector remains usable and the two structures are
not superposed.

BinderScope's `target_frame.py` builds a campaign-aware construct/reference
frame, optionally weights placement around the binder interface, and then
calculates partner-subunit, rest-of-chain, glycan, cut-end and membrane-context
measurements. Its Inspector renders those outputs as distinct scene layers,
target-frame strips and evidence-linked findings. Studio reuses the presentation
idea of displaying the explicitly associated reference around the placed
complex and linking reference residues to Mol*; it does not port BinderScope's
BindCraft run numbering, interface-weighted fit, verdict thresholds or context
calculations. Studio's registered placement uses the existing Kabsch fitter on
the validated mapped Cα pairs. Studio now has a registered optional static
reference-context evaluator for resolved omitted reference regions, explicitly
selected saved reference chains, and explicitly modeled carbohydrate
coordinates. It uses the existing reference association, placement Evaluation,
Database access checks and normal Evaluation ingestion. It does not synthesize
a full assembly or infer biology.

The existing registered interface and ligand Evaluations retain their own
scopes; reference context is measured by a separate method-specific Evaluation:

| Existing Studio output | Reuse for reference context |
| --- | --- |
| Biotite Interface contacts, hydrogen bonds and residue artifacts | Provides evidence for the saved source-complex interface under its existing chain selections. It does not evaluate omitted reference domains, other reference subunits or glycans. Hydrogen bonds require explicit hydrogen atoms. |
| Rosetta Interface Analysis | Provides its registered interface energy and SASA terms for the selected complex chains. Those definitions do not automatically include the remaining reference assembly. |
| ipSAE / saved confidence fields | Retains its model/PAE-derived confidence outputs where inputs exist. It does not measure full-length steric fit or biological accessibility. |
| PLIP and PandaMap | Provide saved ligand–receptor interaction records for the selected pose/receptor. They do not report ligand contacts with reference-only domains or glycans. |
| PoseBusters | Validates ligand pose geometry/chemistry under its own inputs. It is not a receptor-context evaluator. |
| PEP-Patch, Interface Clustering and Foldseek | Retain their surface, contact-mode and fold-similarity meanings respectively. None supplies a full-length context score; Foldseek similarity is not a binding-mode measurement. |

The registered `analysis.reference_environment` workflow reports static
geometry; it does not calculate solvent exposure, conformational accessibility,
membrane orientation, glycan occupancy, or a binding-quality score. For each
context residue it stores closest atom-pair evidence (source and context atom,
element, author residue number and insertion code), distance, atom-pair counts
within 8 Å, and summed configured van-der-Waals envelope overlap above 0.4 Å.
It reports geometry separately for the associated source component and the
other selected-complex chains, so chain continuity within the mapped target is
not confused with partner proximity.
Bondi radii are used for H, C, N, O, F, P, S, Cl, Br, I and Se; the remaining
listed common elements use explicit configured values. Unsupported heavy atoms
are counted by element and excluded from pair calculations, with a limitation
shown in the Evaluation. Bondi values follow A. Bondi, “van der Waals Volumes
and Radii,” *J. Phys. Chem.* 68 (1964), 441–451,
[doi:10.1021/j100785a001](https://doi.org/10.1021/j100785a001); added values
are Studio-configured and are recorded individually in the saved parameters.
The Evaluation also stores total residue/atom-pair counts, the placement fit
coverage and RMSD, source/reference IDs and digests, reference coordinate
origin where explicit file metadata supports it, selected chains, assembly ID,
method/version, and limitations. Positive overlap is reported as potential
geometric overlap, not an experimental clash.

Resolved omitted protein residues in the associated reference chain are
measured separately from explicitly selected other chains and recognized
modeled sugar residues. Sugar components are labeled by coordinate residue
identity; covalent attachment/linkage is not validated. When an assembly ID is supplied, only membership whose saved
operator is the identity transform is confirmed; non-identity assembly
operators and operator expressions are rejected because Studio does not
expand them. A selected chain without an assembly ID is labelled a selected
reference component, not a biological partner. PDB REMARK 350 and mmCIF
assembly metadata are read where available. Sequence-only N-X-S/T motifs are
listed in a separate field and never treated as modeled sugar coordinates.
Unsupported elements without a configured radius, unresolved atoms, alternative
conformers not selected by Bio.PDB, and missing coordinates limit the static
atom inventory. Experimental and predicted reference labels are shown only
when stated by the coordinate file; otherwise origin is unknown.

The following capabilities remain unavailable and need separate methods,
provenance and validation before they can be presented as measurements:

1. **Non-identity biological-assembly expansion.** This needs exact assembly
   operator composition, chain mapping, atom transformation, operator
   provenance and fixtures covering identity, non-identity, compound and
   ambiguous expressions. The current evaluator fails closed for unsupported
   operators.
2. **Full-context solvent exposure or occlusion.** Inputs: a chemically
   complete placed complex, the selected reference assembly, atom/residue
   correspondence, protonation and missing-atom policy, and explicit
   partner/context chains. Existing Rosetta SASA can be evaluated on declared
   selections, but its existing interface output is not this full-context
   comparison. A new evaluator must state solvent model, radii, probe size and
   whether exposure differences are calculated per residue. Artifacts should
   contain per-residue values in named units, exact selections, software/method
   version and input digests. Validate against analytic/simple structures and
   an independent implementation, with tests for missing loops, chain changes
   and atom preparation.
3. **Membrane-relative placement.** Inputs: an explicitly selected membrane
   plane/orientation source (for example a saved membrane-frame structure or
   curated plane annotation) in the same validated reference frame, plus
   residue/atom selections. Without that orientation evidence, Studio cannot
   infer membrane accessibility. Any future output should report signed
   distances/orientation with frame provenance and validate against known
   membrane structures and deliberate frame changes.
4. **Construct-end context.** Existing contact artifacts and validated
   residue correspondence can identify the source construct's mapped termini
   and source-partner residue evidence. A future residue-level view can surface
   contacts near those termini while retaining the Biotite method/parameters;
   it cannot estimate the effect of absent sequence. Tests must include
   insertions, author numbering, missing residues and alternate mappings.

These outputs should become separate residue-linked Findings with evaluator,
definition, units, parameters, exact residue references, artifacts and source
provenance. They must not be merged into a universal quality score or conflated
with predicted confidence or human review.

BinderScope's family board, inspector, and binding-mode map informed the
compatible-collection board and side-by-side inspector below. Its family
labels, campaign filters, and browser-local review state are not portable
generic structure metadata.

### Tag and Fusion Suitability

Studio now has a post-hoc registered Database Evaluation for existing protein
complexes. BindCraft2's terminal design constraints remain design-time inputs
and are not treated as post-hoc evidence. The intrinsic evaluator requires one
selected protein chain, one or more disjoint protein interface-partner chains,
an exact single-model Database representation, and a user declaration about
whether the observed chain ends are the complete intended construct termini.
If the structure declares unresolved terminal polymer residues, that evidence
overrides the user's confirmation for that end. Cyclic peptides and
multichain constructs return explicit not-applicable statuses; a single-chain
component may still be selected independently where scientifically justified.

The registered `analysis.tag_fusion_suitability` workflow records, for both
termini, Cα exit direction from up to six inward residues, distance to the
selected interface residues, interfacial contact status, Bio.PDB
ShrakeRupley residue SASA (240 sphere points), DSSP Sander relative
accessibility and secondary structure when DSSP is available, and terminal
window B-factors. B-factors are interpreted as pLDDT only when the caller
explicitly confirms that source encoding. An explicitly selected, saved
confidence Evaluation may also be retained as separate Evaluation-wide
context; those scalar outputs are not localized to either terminus and do not
grade it.

Flexible His6 and Twin-Strep tags use a shortest clearance path on a 3 Å grid,
with stated exposed-motif and per-residue reach assumptions. GFP/SNAP-sized and
Fc/MBP-sized fusions use spherical envelopes sampled in a fixed 40° exit cone,
with tested linker increments. These are documented geometry approximations
for clearance and linker comparison. They do not simulate the chemical tag,
polymer conformations, tag self-avoidance, tag–receptor binding, expression,
folding, or experimental SPR/BLI behavior. The end comparison and
immobilization-orientation recommendation are advisory summaries of those
measurements, never a universal score or experimental prediction.

The inspected BinderScope `termini.py` and `docs/termini.md` explicitly block
folded fusions when the terminal residue itself contacts the selected target,
while allowing flexible tags to be assessed by path clearance. Studio originally
reported an open sampled fusion cone without carrying that attachment-site
limitation into the summary. Studio now preserves the measured cone/linker
status as geometric evidence and separately marks a folded-fusion advisory as
not recommended at a contacting terminus; it excludes that end from the
advisory end preference. This retains Studio's own clearance measurement while
making its suitability summary consistent with the documented physical
limitation. Flexible-tag clearance and terminal-interface evidence remain
separate.

The separate `analysis.tag_fusion_context` workflow is offered only when the
same assessed component has an explicit, current Database reference
association and a completed, ingested Reference Placement Evaluation matching
the source/reference IDs, models, chains, digests, method, and transformation
direction. The association may be attached to a different protein chain in
the same saved complex; the explicitly selected chain anchors placement of the
entire coordinate model. It applies the saved source-to-reference transform
to a transient in-memory copy and includes resolved reference heavy atoms
outside the mapped source component in the same clearance geometry. Original coordinates and
prior Evaluations remain unchanged. It does not infer a reference association,
refit the structures, assert biological assembly, or calculate glycan,
full-length accessibility, or experimental binding measurements. Without a
reference, the intrinsic result is complete within its selected-complex scope
and lists the external checks not performed.

The Tag & Fusion Inspector panel appears only when a matching saved Evaluation
is available. It displays N/C evidence, native clearance findings, advisory
comparison, method parameters and provenance. Each terminal button selects its
exact author-numbered residue, including insertion codes, in the existing Mol*
viewer. Results are saved as normal Evaluation records against the existing
Database molecular identity. BinderScope method attribution and the MIT notice
are in `docs/third-party-notices.md`.

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

The inspected BinderScope source is MIT licensed. Studio's Tag and Fusion
evaluator independently implements and attributes selected geometric and
presentation concepts from `binder_review/termini.py` and `docs/termini.md`;
no BinderScope Python source was copied or vendored. The upstream copyright and
MIT permission notice are retained in `docs/third-party-notices.md`. The
license does not grant separate rights to BindCraft 2; its campaign launcher
and scheduler integration are excluded from Studio.

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

The full Complex Analysis, Compare and Collection pickers now request 100
structure summaries initially instead of rediscovering the first 1,000 on
every Streamlit rerun. A local call-path profile on 2026-11-03 measured the
initial 100-record page at 0.146 s, ACL-revalidated rerenders at 0.128 s and
0.109 s, loading the next page at 0.271 s, and the subsequent 200-ID
ACL-revalidated rerender at 0.173 s. Cached state contains stable IDs only;
each rerun resolves labels and minimal representation summaries through the
current ACL-filtered Database query. Page navigation remains offset-based, and
revoked IDs trigger a current authorized-prefix rebuild so changes in visible
rows do not skip records. Coordinate files and Evaluation histories are loaded
only after selecting a structure.

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

### Reference environment Evaluation verification (2026-11-03)

The registered `analysis.reference_environment` workflow was run in an isolated
Database against the real experimental 1IEP X-ray structure in
`assets/vina/1iep.pdb`. The saved source representation was a documented
controlled truncation of chain A to author residues 280–450 while retaining
chain B; it was associated only inside the fixture Database with chain A of
the original file. This is not a production biological-reference association.
The exact saved placement used `Biopython PairwiseAligner affine semiglobal
sequence mapping`, aligned 171 Cα pairs, source coverage 1.0, reference
coverage 0.6241, and RMSD `4.61e-15 Å`.

The Evaluation measured 103 resolved omitted-reference residues and 13,325
heavy-atom pairs within 8 Å. It counted 17 configured-radius envelope overlaps
greater than 0.4 Å across six context residues; these are potential static
overlaps, not experimental clashes. The role-separated totals showed 16
overlaps against the associated target component itself and one against the
other source-complex component. The closest other-component pair was source
chain B residue GLN 447 OE1 to reference chain A residue TYR 456 OH at 2.6074 Å
(0.4326 Å configured-radius overlap). The all-source minimum of 1.3278 Å was
at the truncated target-chain boundary and is therefore shown separately from
the other-component minimum. No unsupported heavy atoms or modeled sugars
occurred in this structure. A controlled unit fixture separately verifies a
modeled NAG residue and sequence-only N-X-S/T motif remain distinct; the
evaluator does not validate covalent sugar attachment or glycan occupancy.

Both placement and context Evaluations were ingested into the same fixture
Database. Reopening the Database record restored the saved measurements and
the Inspector generated source-residue Mol* selectors for the nearest other
component evidence. The workflow also verified authorized association checks,
stale-digest rejection, duplicate residue-map rejection, rigid-transform
invariance, and unavailable context behavior. The full relevant regression
selection passed 166 tests after final edits. One Starlette `TestClient`
deprecation warning came from the existing Database policy test; no tests were
skipped.

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
