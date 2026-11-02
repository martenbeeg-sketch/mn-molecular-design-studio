# Complex Analysis workspace

## Operational flow

Complex Analysis is the contextual workspace for selecting a molecular
structure, choosing applicable existing analyses, configuring the selected
representation and component mapping, invoking Studio's registered workflow,
and viewing the resulting evaluation. All Tools remains the standalone catalog
and setup surface.

The workspace currently loads structures from Database or completed Results,
uses the existing PDB/mmCIF component parser and Mol* viewer, and retains
per-user persistent human review notes separately from scientific
observations. Database and Results reads use their current ACLs and stable
record/evaluation or job/candidate/artifact identifiers.

### First execution slice: Biotite Interface

For completed prediction Results with normalized candidates, the workspace:

1. Selects a specific Results job and structure artifact/model.
2. Detects the linked candidate and its recorded model metadata.
3. Lets the user select two disjoint protein chain groups, including multiple
   chains per group. The labels do not imply a binder or experimental binding.
4. Checks existing Biotite Interface analysis jobs for the same source job,
   candidate, structure digest, chain mapping, and contact distance.
5. Reuses an exact completed match, shows an active job, or offers an explicit
   **Run analysis** / **Retry analysis** action.
6. Submits `analysis.biotite_interface` through the existing authenticated
   workflow API and scheduler. The registered workflow validates the parameters,
   reads the selected candidate from the source Results collection, stages the
   selected structure into the analysis job, and records the exact model and
   chain mapping in its existing manifest and normalized candidate provenance.
7. Shows saved measurements, progress, failure details, contact maps, and the
   persistent Evaluation identifier in this workspace. The existing Results
   ingestion and evaluation projection remain authoritative.

Opening a structure never starts a calculation. This first contextual execution
path currently requires a completed Results job, an associated normalized
candidate, a local accessible PDB/mmCIF artifact, and a project assignment.
Database-only imported structures remain viewable, but direct analysis launch
from such a Database record is a later adapter phase. A multi-model coordinate
file is treated as one selected artifact; the workspace does not yet split its
internal models into separate selectable records.

Biotite measurements are structural/geometric observations. Hydrogen bonds
are available only when explicit hydrogens are present; missing hydrogen-bond
inputs remain unavailable. Neither contacts nor model confidence establish
experimental binding.

## Studio and BinderScope capability matrix

| Capability | Studio analysis already available | BinderScope capability | Decision and planned integration |
| --- | --- | --- | --- |
| Structure parsing, components, chain IDs, provenance | `molecule_imports.parse_structure`, Database representations, Results artifacts, Mol* viewer, and candidate/evaluation lineage. | `structure.py` parses structures and assigns component/chain roles. | **Already reusable.** Keep Studio parser, access control, viewer, and provenance contracts; do not copy the BinderScope parser. |
| Heavy-atom contacts and residue maps | Biotite Interface computes contacts and residue-pair maps for declared protein chain groups. Results artifacts preserve source residue identifiers. | `interface.py` generates interface contacts/maps and 3D review layers. | **Already reusable for protein–protein.** Complex Analysis now selects groups, launches Biotite, and renders the existing contact-map artifact. Residue-driven Mol* highlighting remains a follow-up. |
| Buried solvent-accessible area (SASA) | Rosetta Interface Analysis reports `dSASA`, interface SASA fraction, interface energy, and related measures through `analysis.rosetta`. | `interface.py` reports interface buried area. | **Already implemented in Studio, not yet contextualized.** Connect the existing Rosetta workflow and renderer; do not reproduce its SASA calculation. Keep its Rosetta method and units visible. |
| Hydrogen bonds | Biotite Interface counts geometric hydrogen bonds when explicit H atoms are present. Rosetta Interface Analysis reports its own interface H-bond metrics. | `interface.py` reports interface hydrogen bonds and geometry. | **Already partially implemented.** Reuse the two existing evaluators with their method-specific definitions and explicit missing-input states. Do not merge their counts. |
| Unsatisfied polar atoms / interface H-bonds | Rosetta Interface Analysis reports interface unsatisfied H-bond metrics. Biotite does not calculate unsatisfied polar atoms. | `interface.py` reports unsatisfied polar atoms. | **Already implemented by Rosetta for its protocol.** Expose through the existing Rosetta workflow; identify as Rosetta output, not a universal polar-satisfaction score. |
| Salt bridges and contact geometry | Biotite Interface reports charged-side-chain pairs within its declared 4 Å geometric cutoff. It does not infer energetic strength. | `interface.py` includes polar geometry and salt-bridge context. | **Already reusable.** Present as geometric contacts with chain/residue numbering and cutoff provenance. |
| Hydrophobic surface patches | PEP-Patch is an existing workflow for selected protein chains with hydrophobic and pH-specific electrostatic surface patches, per-residue artifacts, areas, and fractions. | `interface.py` and `liabilities.py` report surface patches/chemical properties. | **Already partially implemented.** Connect PEP-Patch as a contextual workflow. Keep its patch definition and declared pH; do not relabel patch area as affinity or specificity. |
| Full-length target context and residue mapping | Studio retains target/reference IDs and source residue numbering. ipSAE can use saved structures, PAE, and design-reference evidence when available; this is not a general full-length target-frame mapper. | `target_frame.py` aligns cropped targets to optional full-length references and maps contacted residues. | **Partially present, requiring extension.** First reuse explicit references and existing ipSAE evidence. Plan a separate reusable context-mapping workflow for full-length placement, target residue mapping, and missing-reference reporting. Do not assume residue numbers match across constructs. |
| Glycans and construct accessibility | Some prediction/design inputs can carry glycan-related structure context, but Studio has no general glycan-occlusion evaluation in Complex Analysis. | `target_frame.py` checks glycan context and clashes. | **Missing and worth integrating.** Add only as an explicit reference-context workflow with supported glycan atoms/links and a transparent unavailable state. Do not infer glycan occupancy or biological accessibility from an absent structure. |
| Steric clashes | Existing measures include ipSAE Cα clash counts, selected validation Cα geometry, Rosetta/other workflow outputs, and PoseBusters ligand plausibility checks. They are distinct methods and do not provide a universal full-length all-atom clash test. | `target_frame.py` checks clashes in full-length context. | **Partially implemented.** Connect existing measurements with source/method labels. A full-length all-atom context check is a later extension; do not substitute a Cα threshold for it. |
| Binding-mode clustering | Interface Clustering compares contacted target residues by Jaccard similarity for compatible target revisions/numbering. Sequence Clustering and Foldseek are separate existing collection workflows. | `modes.py` and `families.py` compare binding modes and families. | **Partially implemented.** Reuse Interface Clustering for compatible same-target collections and Foldseek/sequence clusters for their existing question. Plan to contextualize collection selection and outputs; never cluster incompatible targets or residue frames together. |
| Structural comparison across engines/runs | Studio stores distinct model/pose evaluations and supports selected metric comparisons. Foldseek is available for chain-level structure clustering. There is no generic multichain complex superposition/comparison mode in this workspace. | `modes.py` compares complex poses and binding modes. | **Missing and worth integrating.** Add a comparison workflow only for explicitly compatible structures, target revisions, chain mappings, and metric definitions; keep samples/poses independent. |
| Interface specificity | Some prediction/design workflows accept explicit counter-target or specificity panels. Studio has no general saved-pose homolog/paralog specificity evaluation. PEP-Patch is not a counter-target assay. | `specificity.py` evaluates ortholog/paralog panels. | **Partially present at the workflow level; missing as a general analysis.** Reuse existing explicit-panel workflows where compatible, then design a reusable post-prediction panel adapter. Do not claim proteome-wide specificity from one target model. |
| Sequence/structural liabilities | Sequence Composition, Protein-Sol, ESM-1v/ESM-IF1, PEP-Patch, PoseBusters, and some engine-specific developability filters exist. They measure different properties and are not interchangeable. | `liabilities.py` includes sequence and structural liability screens. | **Partially implemented.** Surface existing Studio workflows with their evaluator and provenance. Plan a generic plugin interface for additional well-defined liability checks; no composite pass/fail verdict or invented measurement. |
| Fusion/tag accessibility | Some BindCraft2 modes constrain termini during generation. Studio has no generic post-hoc tag/fusion accessibility analysis for arbitrary structures. | `termini.py` evaluates N/C terminus and fusion/tag accessibility. | **BinderScope-specific behavior requiring generalization.** Do not reuse a generation-time constraint as a post-hoc measurement. Plan a separately named, structure-context workflow that accepts explicit tags/fusions and reports geometric assumptions. |
| Protein–ligand interactions | PLIP and PandaMap provide saved-pose interaction analysis; PoseBusters checks pose/chemistry plausibility. The existing workflows retain pose and source compound provenance. | BinderScope is primarily an antibody/protein-complex analysis tool. | **Already implemented in Studio.** Connect these workflows for the selected ligand pose, reusing their renderers and compound identity mapping. |
| Protein–DNA/RNA interactions | Structure parsing and visualization support DNA/RNA components. Current Biotite Interface requires protein chains; no registered contextual protein–nucleic-acid interaction evaluator is available in this slice. | BinderScope's protein-chain assumptions do not make these analyses applicable. | **Not currently supported by the contextual analysis tools.** Keep structures selectable and visible, state unsupported checks explicitly, and add only scientifically suitable existing or separately reviewed workflow modules later. |
| Source viewer and human review | Mol* and authenticated Database/Results selection are reused; per-user persistent notes remain separate from observations. | BinderScope has interface review functionality. | **Already reusable/adapted.** Human review is secondary to workflow execution and cannot change scientific measurements. |

### Planned implementation phases

1. **Complete Biotite Interface vertical slice** (current): exact Results
   candidate/model, chain groups, reusable evaluation detection, workflow API
   submission, job state, contact map/measurements, and source/evaluation
   provenance.
2. Connect existing protein–protein analyses: Rosetta Interface Analysis,
   ipSAE, PEP-Patch, and compatible Interface Clustering collections.
3. Connect existing PLIP, PandaMap, and PoseBusters workflows for selected
   protein–ligand poses.
4. Add reusable full-length context/reference mapping, glycan context, and
   compatible steric-clash checks only where a validated existing workflow or
   separately reviewed module supports them.
5. Add compatible cross-run/model comparison, explicit-panel specificity, and
   protein liability integrations while keeping predicted, measured,
   geometric, and empirical quantities separate.
6. Add scientifically appropriate protein–DNA/RNA analyses when a registered
   supported workflow exists. Keep unavailable capabilities explicit.

The All Tools catalog remains available for standalone workflows throughout.
Complex Analysis should become their preferred contextual entry point; it must
not create a second scheduler, job store, evaluation store, or copied scientific
algorithm.

## Comparison rules and limits

One saved observation remains one model/pose/calculation. The workspace does
not pool distinct poses or models. A reusable evaluation match is scoped to its
source Results job, source candidate, selected structure digest, chain groups,
method, and scientific parameters. Interface Clustering remains limited to a
compatible target revision and residue numbering. Missing structures,
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

Focused tests cover component detection for protein–protein, protein–ligand,
protein–DNA and protein–RNA structures; model and chain-group selection;
workflow schema validation; mocked workflow staging/output provenance;
existing-evaluation matching; unavailable results; and job/Database access
boundaries. Workflow tests mock the container execution. No live analysis job,
GPU generation, deployment, or push is performed as part of the test suite.
