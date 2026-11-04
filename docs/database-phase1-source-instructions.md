# MN Molecular Design Studio — Phase 1: Database Simplification, Scientific Data Integrity and Historical Import

## Mission

Implement Phase 1 of the Molecular Database restructuring in the existing MN Molecular Design Studio repository.

The objective is to establish a **simple, scientifically reliable and maintainable Database ingestion and Evaluation system** before considering a clean rebuild of the current Database index.

I want the Database to function as a curated molecular library containing important reusable molecular identities, selected structures and scientifically meaningful results—not as a copy of every computation performed by Studio.

This phase must also establish a reliable way to import historical results, especially native BindCraft1 and BindCraft2 campaigns.

**Use the existing architecture wherever possible. Do not create another molecular Database, another Evaluation store or another scheduler.**

Work in the current repository and preserve all existing work, including uncommitted changes.

## 1. Start by reviewing the forensic audit

Read the completed audit documents:

- `database-forensic-audit.md`
- `database-inventory.csv`
- `data-classification.csv`
- `ingestion-root-causes.md`
- `reconstruction-feasibility.csv`
- `architecture-options.md`

Locate the original audit files in the current Studio workspace. Do not substitute an older audit or assume a path without checking.

Also review:
- `docs/protein-library.md`
- `docs/generator-output-vocabulary.md`
- `docs/molecular-database-behavior-recheck-20261020.md`
- `docs/audits/2026-11-03-complex-analysis/README.md`
- Relevant existing scientific and ingestion regression tests.

The audit recommends **Option B only after implementing the ingestion-policy corrections associated with Option A**.

We are implementing those prerequisites now.

Do not rebuild production data in this phase.

### Important audit findings

The November 4 audit reported approximately:

- 21,473 molecular identities.
- 21,537 structure records.
- 41,198 raw observations and corresponding Evaluation projections.
- 1,050 groups, including 377 empty groups.
- Approximately 15,000 automatically retained UniDock Pro screening input compounds.
- 3,848 initially identified reusable molecular identities under the proposed retention policy.
- 517 generated memberships lacking explicit acceptance status.
- 2,230 structures lacking straightforward representation-type classification.
- A 128-observation/source-link discrepancy between primary and peer-visible Database snapshots.
- A large replication-tombstone history that must not be reset casually.

These findings are a snapshot, not instructions to delete records.

Do not assume every unselected screening compound is scientifically worthless or that the initial 3,848 identities constitute a complete migration set.

## 2. Preserve the existing scientific architecture

Review the current implementations, especially:

- `molecular_database.py`
- `molecular_database_ingest.py`
- `molecular_database_peer.py`
- `molecular_database_results.py`
- `molecule_imports.py`
- `collection_policy.py`
- `protein_library.py`
- `evaluations.py`
- `output_vocabulary.py`
- `accepted_design_import.py`
- `complex_analysis_sources.py`
- `complex_analysis_dashboard.py`
- Relevant Database, Results, and analysis APIs.
- Existing workflow manifests and output adapters.

Identify which mechanisms already work correctly before modifying anything.

The existing logical concepts should remain:

**Molecule identity**
- A canonical protein, peptide, nucleic acid, compound or defined multicomponent system.

**Structure representation**
- Experimental/supplied PDB, predicted structure, generated structure, complex, docking pose or other coordinate representation.

**Scientific Evaluation**
- A scientifically defined measurement associated with exact inputs, parameters, model and method.

**Collection**
- A reusable grouping or explicitly curated selection.

**Results provenance**
- The original run, candidate, workflow, workstation and scientific artifacts.

These concepts may use multiple tables but must remain part of one consistent molecular library.

Do not introduce independent databases for proteins, compounds, structures, references or Evaluations.

## 3. Establish explicit Database inclusion rules

The principal architectural correction is that completing a computational job must not automatically mean that every input, output and intermediate becomes a reusable Database entry.

Create an explicit, understandable retention and ingestion policy for every workflow category.

### 3.1 Source proteins and structures

Automatically retain genuinely reusable source material when it is intentionally supplied to Studio:

- Imported protein sequences.
- User-supplied or imported PDB/mmCIF structures.
- Target protein sequences and source structures.
- Source structures used for sequence extraction.
- Explicitly imported compounds and ligand collections.
- Other recognized reusable biological or chemical source material.

Preserve accession, source file, provenance, composition, chain relationships and the exact way the source was used.

Do not turn a sequence-only target into a fictitious experimental or native structural reference.

Do not treat an idealized sequence-carrier CIF as an experimental structure.

Retain the recently implemented Protein Library functionality, including source collection during job submission and its ability to preserve source material even if a subsequent prediction fails.

However, source capture must not automatically promote all derived validation candidates or generated outputs.

The Protein Library should be a reference/source view within the common molecular library, not a parallel identity system.

### 3.2 Designed proteins and binders

For BindCraft2, RFD4, RFdiffusion-derived designs, BoltzGen and other generators:

- Keep complete generation history in My Results.
- Retain original native acceptance/filter evidence.
- Promote explicitly accepted designs, with their scientific policy and required structures, to reusable Database collections.
- Support explicitly selected additional designs even when their original status was rejected, unknown or unresolved, but retain that status.
- Do not treat unknown status as accepted.
- Do not promote all raw trajectories or unsuccessful refolds.
- Do not assign a new molecular identity merely because the same sequence appears in another run.

A validated target–binder complex may be stored as a representation associated with the relevant designed binder identity while preserving its relationship to the exact modeled target.

Do not silently equate a cropped target sequence with a complete target construct.

### 3.3 Protein structure prediction and cofolding

Protein structure prediction should generally:

- Reuse the existing appropriate molecular identity.
- Attach new predicted structure representations to that identity.
- Preserve model, engine, sequence, parameters and source run.
- Preserve multicomponent relationships where applicable.
- Avoid creating unnecessary new molecular collections for every prediction attempt.

Unselected models and intermediate attempts may remain exclusively in Results.

If there is no valid Database identity, use the established source-registration or explicit promotion process rather than creating arbitrary duplicate records.

### 3.4 Small molecules and docking

- Reuse canonical chemical identities.
- Preserve stereochemistry and scientifically significant protonation/tautomer distinctions.
- Record docking poses as representations or result evidence, not automatically as new compounds.
- Retain exact receptor identity, receptor preparation, ligand preparation, docking engine, scoring protocol and pose identity.
- Attach new measurements to the corresponding existing compound and pose where scientifically valid.

For large screening campaigns, keep the complete screening library and results accessible in Results/source storage without automatically placing every screened compound into the reusable molecular library.

Support explicit selection and promotion of hits.

An explicitly imported compound library is different from automatically indexed screening inputs; preserve that distinction.

### 3.5 Analysis-only workflows

Biotite, Rosetta, ipSAE, PEP-Patch, sequence liabilities, tag/fusion, reference context, Foldseek, clustering, PLIP, PandaMap, PoseBusters and other registered scientific analysis workflows should not normally create new molecular identities simply because an analysis was executed.

Instead, attach their Evaluations to the appropriate existing molecular and structure identities.

Do not create a new reusable collection for every analysis job.

If the input identity or structure cannot be resolved unambiguously, preserve the result in Results and report that Database association is unresolved.

Never silently attach an Evaluation to an incorrect molecule or pose.

## 4. Correct ingestion triggers and redundant collection creation

Investigate all ingestion entry points.

The audit identified overlapping paths involving:
- Job completion hooks.
- Workspace synchronization/scans.
- Database API GET requests.
- Explicit run imports.
- Existing source-Protein Library capture.
- Peer replication.

Implement one coherent contract for each operation.

Requirements:

1. Database reads must not unexpectedly trigger ingestion or mutation.
2. Explicit refresh/reindex operations must have controlled, understandable semantics.
3. Completion hooks must respect the declared workflow collection policy.
4. Repeated ingestion of unchanged runs must be idempotent.
5. Source capture and output promotion must remain distinguishable.
6. Empty groups should not be created by unsuccessful ingestion.
7. Analysis jobs must reuse existing identities and appropriate groups.
8. Import/promote operations must not silently rewrite unrelated curated collections.
9. Historical acceptance decisions and source scientific data must remain unchanged.
10. Existing ACL and cross-workstation ownership boundaries must remain enforced.

Do not interpret `automatic` as the universal default for unclassified jobs.

Use explicit workflow purpose and safe defaults.

Carefully preserve compatibility with historical runs whose policy was never recorded.

A missing legacy policy must not automatically become evidence that the corresponding output should be retained.

## 5. Verify automatic attachment of Evaluations

This is an essential requirement.

**When I select a molecule or saved structure in Database and analyze it using another Studio tool, I expect the resulting scientific measurements to become available under the original molecular record automatically, provided the analysis completed successfully and exact source identity is established.**

I should not have to manually import the measurement afterward.

Inspect every relevant execution path:

- Complex Analysis.
- Standalone analysis tools.
- Analyses launched from Results.
- Analyses launched from Database.
- Pipeline/canvas analysis nodes.
- Independently submitted cofolding or validation analyses.
- Protein–ligand docking and pose-analysis workflows.
- Historical normalized outputs.

Determine which already perform correct automatic association and which do not.

### Required Evaluation identity

Preserve and verify, as applicable:

- Canonical molecular identity.
- Structure representation ID.
- Exact structure digest.
- Model number or pose identity.
- Binder and target chains.
- Receptor/ligand identity.
- Source run and candidate.
- Scientific evaluator and version.
- Effective calculation parameters.
- Calculation stage.
- Source workstation.
- Evaluation status.
- Original native measurements and provenance.

One molecule may have many structures.

One structure may have many Evaluations.

One Evaluation must refer to its exact scientific inputs.

An analysis of structure A must not appear as though it were performed on structure B merely because both represent the same protein sequence.

An analysis of one ligand pose must not be presented as an evaluation of all poses of that compound.

### Existing versus new measurements

When a user requests analysis:

- Discover compatible completed Evaluations.
- Reuse them when their exact input identity and relevant parameters match.
- Keep scientifically distinct measurements separately.
- Identify stale or incompatible Evaluations.
- Do not infer that missing means failed.
- Preserve failures and unavailable inputs distinctly.
- Do not automatically rerun computationally expensive analyses because a page was opened.

Add regression tests for these cases.

## 6. Standardize scientific metric names and values

Studio already has `output_vocabulary.py` and `evaluations.py`.

Reuse them.

Do not create a second normalization dictionary.

Audit how complete the existing normalization actually is, including its application to raw Results, Database observations, Complex Analysis, tables, filters, comparisons and exports.

I want scientifically identical measurements to be searchable and comparable under consistent names, while differences in scientific definition must remain explicit.

For every supported normalized measurement, retain:

- Original field name.
- Original value.
- Source engine and version.
- Canonical field name.
- Scientific definition.
- Unit and valid conversion.
- Chain/component scope.
- Stage: generator, validation, independent analysis, etc.
- Method and estimator.
- Effective parameters.
- Conversion provenance.
- Any alias conflict or invalid measurement status.

Examples requiring particular care:

- pLDDT versus binder pLDDT, target pLDDT and complex pLDDT.
- Native generator iPTM versus independent validation iPTM.
- BindCraft2 and RFD4 interface PAE definitions.
- PAE reported as Å versus normalized by 31 Å.
- ipSAE versus BindCraft2-definition i_pDAE.
- Different RMSD definitions and alignment references.
- Interface contact counts at different cutoffs.
- Buried surface area versus differently defined interface SASA.
- Rosetta energy versus native predictor confidence.
- PLIP/PandaMap interaction counts and protocol differences.
- PoseBusters evaluations with distinct input conditions.

Do not normalize measurements based only on similar names or numerical ranges.

If two values cannot be scientifically compared, preserve their distinct normalized identities.

Avoid numerical alias conflicts and ambiguous metric provenance.

The existing vocabulary already has mechanisms for alias conversion, context checks, and conflict reporting. Extend those rather than duplicating them.

## 7. Normalize workflow parameters too

This is related to, but different from, metric normalization.

Investigate whether effective input parameters are consistently represented across versions of the same engine and across scientifically equivalent workflows.

For example:
- Contact distance thresholds.
- Interface definition settings.
- PAE aggregation and cutoff settings.
- Structural alignment references.
- Chain and residue selections.
- Number of prediction models.
- Recycles, seeds and sampling settings.
- Acceptance criteria and thresholds.
- Filters and model versions.
- Docking scoring and preparation settings.

Preserve the original native parameter names and values.

Add common parameter definitions only when equivalence is established.

Do not rewrite historical settings to current defaults.

Missing historical parameters must remain missing or unresolved unless a documented, version-specific default can be established.

The normalized measurement identity must not conceal differences in scientifically important calculation settings.

Expose both readable common labels and exact native settings in the user interface.

## 8. Implement a controlled historical results importer

This should become a general capability of Studio, beginning with **BindCraft1 and BindCraft2**.

Study the relevant open-source BinderScope readers:

https://github.com/schnickn/binderscope

In particular, inspect its BindCraft1/BindCraft2 reader implementation, documented compatibility limits and native run-folder structure.

Use it as a reference for file discovery, scientific reconciliation and provenance—not as a reason to replace Studio's existing Results or Database models.

Inspect Studio's existing `accepted_design_import.py` and normalized-candidate infrastructure.

Do not assume that importing a Studio-accepted candidate bundle is equivalent to importing an original external BindCraft campaign.

### Historical import workflow

Design and implement a safe, initially non-mutating workflow:

**Select folder/archive → Detect format and version → Inspect native results → Reconcile evidence → Preview Database changes → Explicit import authorization.**

During Phase 1, implement the reader, validation and dry-run preview with isolated fixtures.

Production import remains disabled until a later explicit authorization.

### BindCraft1

Support the documented native statistics and structure files where compatible.

Detect unsupported or incomplete layouts and report them clearly.

### BindCraft2

Account for its native stages:

- `1_Trajectories`
- `2_Refolded`
- `3_Ranked`

Read the actual native statistics tables, available campaign metadata, final structures and source target information.

Reconstruct trajectory, refolded candidate and ranked design relationships using exact native identifiers, trajectory hashes and sequences.

Do not use ambiguous filename substring matching.

The final ranked designs are not necessarily all designs that passed the native refolding filters.

Preserve separately:
- Started trajectories.
- Completed/refolded candidates.
- Native filter passes.
- Ranked designs.
- Unknown or incomplete outcomes.

Do not count a ranked representation and its corresponding refolded candidate as two different molecular identities.

Preserve original native acceptance/filter rules and distinguish them from subsequently calculated Studio evaluations.

### Historical source targets

Capture and validate:

- Target protein sequence.
- Original source PDB/mmCIF when available.
- Construct and crop identity.
- Chain mapping.
- Target numbering and hotspots.
- Whether coordinates were directly supplied or used only to derive sequence.
- Source accessibility and exact file digest.
- Any uncertainty about correspondence.

Do not claim that a predicted or sequence-derived target has an experimental native structure.

### Historical metrics

Preserve native BindCraft results unchanged.

Map compatible metrics through the existing normalized vocabulary.

Retain original filters and parameters.

Do not recompute or replace historical acceptance decisions as a side effect of import.

### Import preview

Before importing, display:

- Detected engine and version.
- Target identity and source-structure status.
- Total trajectories/refolds/passes/ranked designs.
- Existing matching molecular identities.
- Proposed new identities.
- Structures to attach.
- Evaluations to attach or reuse.
- Proposed collections.
- Duplicates and conflicts.
- Missing or ambiguous evidence.
- Unrecognized metrics or unsupported formats.
- Explicit retention selection.

Import previews must have no production Database side effects.

Use deterministic source content identities and checksums so a repeated import, renamed archive or relocated folder does not silently duplicate results.

The full historical computational run should remain associated with Results/source provenance. Only selected reusable molecular content should be promoted to Database.

## 9. Test against existing real scientific examples

Use read-only data inspection and isolated test fixtures.

Representative cases:

**R3 BindCraft2**
- The runs with 65 and 35 ranked designs.
- Verify native trajectory/refold/pass/ranked counts are not mixed.
- Check exact accepted binder sequence and target identity.
- Check that existing Database Evaluations are recognized instead of duplicated.

**R3 and Tau RFD4**
- Verify correct target sequence despite potentially misleading historical group names.
- Preserve accepted, rejected and unresolved state.
- Preserve original target/binder chain and modeled target provenance.

**PD-L1**
- Verify the structured-target record and associated source PDB/provenance.
- Preserve the distinction between its final acceptance policy and native PPI diagnostic.

**Protein Library**
- Verify the PDB-derived 8JYR and 7QF4 examples.
- Preserve source PDBs separately from derived sequence identities.
- Avoid duplicate protein identities when the same input is reused.

**Small molecules**
- Verify an existing compound appearing in multiple docking runs.
- Confirm that independent scientific evaluations attach appropriately without generating duplicate compounds.
- Verify that large screening input libraries are not automatically treated as selected hits.

**Scientific analyses**
- Check Biotite, Rosetta, ipSAE, PEP-Patch and protein–ligand analysis output association.
- Confirm exact model/pose/digest/chain and parameter matching.
- Confirm that available, missing, failed, incompatible and unavailable values are handled distinctly.

Do not launch new scientific generation or expensive GPU analysis to perform these checks.

## 10. Address multi-workstation considerations

The system operates on:

- pc002975
- PC004032
- PC004033

The current audit reports small but unresolved differences in materialized scientific records and a substantial replication-tombstone history.

Phase 1 must not reset replication, discard tombstones or reconstruct live peer state.

Review how the proposed ingestion and Evaluation changes interact with:
- Owner-authoritative records.
- Replicated metadata.
- Remote source structures.
- Access permissions.
- Exact file digests.
- Durable removals.
- Source workstation availability.
- Idempotent replay.

Provide a precise, read-only reconciliation plan for the 128-observation/source-link discrepancy.

Do not mistake temporary peer unavailability for missing scientific evidence.

A clean Database rebuild will be considered only after these boundaries are understood and verified.

## 11. Testing and implementation discipline

Work in the current Git working tree.

Before editing:
1. Inspect current branch and commits.
2. Check all uncommitted changes.
3. Determine whether earlier Luna work already implements part of this task.
4. Identify overlapping implementations and avoid duplicating them.
5. Write a concise implementation plan based on actual code.

Then implement Phase 1 incrementally.

Prefer small changes to existing modules.

Create focused regression tests for:
- Source capture and reuse.
- Automatic versus manual collection policies.
- Unknown/rejected/accepted designs.
- Prediction and analysis-only ingestion.
- Repeated ingestion.
- Prevention of empty groups.
- Automatic exact Evaluation attachment.
- Cross-run metric normalization.
- Parameter/alias conflicts.
- Historical BindCraft1/BindCraft2 parsing.
- Multiple candidates per trajectory.
- Repeated historical import.
- Missing or incompatible source structures.
- Cross-workstation ownership boundaries.

Do not weaken existing scientific input validation or remove historical provenance.

Do not create synthetic acceptance decisions merely to make tests pass.

Do not use production scientific records as mutable test fixtures.

## 12. Strict production safety

This phase is **code and isolated-test work**, not an authorization to migrate or clean the production Database.

Do NOT:
- Delete or purge production data.
- Drop or recreate production tables.
- Reclassify historical production records.
- Recompute acceptance decisions.
- Bulk import historical runs into production.
- Replace the live Database.
- Reset replication or remove tombstones.
- Restart production services.
- Enable bulk peer sharing.
- Synchronize or deploy to the compute workstations.
- Execute new expensive scientific workflows.
- Push code to GitHub.
- Perform an automatic production reindex or workspace refresh as a test.

Use explicitly isolated temporary databases and fixtures for all write tests.

A dry-run report is permitted; it must not persist scientific changes in production.

Preserve the existing original Results and scientific artifact files.

## 13. Deliverables

Finish with:

### A. Implementation summary
- Changes to ingestion behavior.
- Existing systems reused.
- New historical import capabilities.
- Automatic Evaluation attachment corrections.
- Normalization extensions.
- Preserved compatibility boundaries.

### B. Workflow inclusion-policy matrix
For each supported workflow:
- Inputs retained?
- Generated identities retained?
- Structures retained?
- Evaluations attached?
- Explicit user selection required?
- Original Results retained?
- Acceptance evidence required?

### C. Evaluation attachment and normalization audit
- Which tools already attached Evaluations correctly.
- Which paths were corrected.
- Which metric aliases are scientifically equivalent.
- Which similar-looking metrics must remain distinct.
- Which input parameter mappings are verified.
- Which remain unsupported or ambiguous.

### D. Historical importer report
- BindCraft1 and BindCraft2 layouts supported.
- Version/modality restrictions.
- Counts and identity reconciliation on real historical examples.
- Dry-run proposed Database changes.
- Unresolved scientific/source evidence.

### E. Test report
- Tests executed and results.
- Tests not executed.
- Relevant regression coverage.
- Any remaining scientific uncertainty.

### F. Phase 2 readiness report
Explain:
- What is now safe to reconstruct.
- What requires additional evidence or manual curation.
- How to preserve Database-only information.
- Whether new ingestion rules would reproduce unwanted records.
- Remaining replica inconsistencies.
- What an offline shadow Database would contain.
- What exact validation must precede a production cutover.

Commit verified Phase 1 code changes locally in logical commits, preserving unrelated pre-existing modifications.

**Do not push, deploy, synchronize peers or modify production data.**

## Final expected outcome

At the end of Phase 1, Studio should have a reliable, tested ingestion and Evaluation architecture in the codebase.

Future runs should follow explicit retention rules.

Existing molecular identities should be reused.

New scientific analyses should automatically attach to the correct saved molecules/structures when exact identity is established.

Measurements should be consistently named and normalized only where scientifically valid, while retaining raw values and full provenance.

Historical BindCraft1 and BindCraft2 campaigns should be inspectable through a clean, validated dry-run import path.

The system should be prepared for a selective reconstruction of the existing Database, but **the actual reconstruction and production migration belong to Phase 2 and must not be performed now.**


---

**Source-integrity note (editorial; not part of the original request):** The supplied original attachment ended after the first Section 9 small-molecule bullet with an incomplete `-`. The user then supplied the missing continuation recovered from the complete original prompt. This document removes that dangling marker and inserts the supplied continuation verbatim. The continuation was checked byte-for-byte against the user-provided recovery text. The original attachment SHA-256 is `52157c0c386358a5119acd6cd900f385081b3ea7d1dc1aca1497a9e800505360`; the recovered continuation SHA-256 is `a2408349dd31ff9451f8a994025767d02a849f61f2cccbec067c0ca3d2accbad`. The two subsequent clarifications follow unchanged.

---

# Subsequent clarification: independent local ingestion and peer replication

Additional clarification for the Phase 1 implementation currently underway:

An essential architectural requirement is that **pc002975, PC004032 and PC004033 each maintain their own local Molecular Database and independently generate, import and analyze scientific data**.

The existing peer-replication architecture must synchronize these locally authored records into one logically consistent molecular library.

Please incorporate this requirement into your current implementation and tests, without restarting the work or expanding into Phase 2.

In particular:

- Preserve independent local ingestion, even when other workstations are offline.
- Ensure that canonical molecular identities converge across machines without duplication.
- Preserve distinct structures and Evaluations from different workstations.
- Allow analyses performed on one workstation to attach to the exact existing molecular identity and structure originating on another.
- Preserve source ownership, provenance, access permissions and remote structure retrieval.
- Make legacy BindCraft imports idempotent across machines, even when archives or folders are copied or relocated.
- Verify offline/reconnect behavior, replay, duplicate prevention and deletion tombstones using isolated multi-peer tests.
- Account for the unresolved peer discrepancies identified in the forensic audit.

Do not centralize ingestion on pc002975.

Do not modify production Database contents, restart services, synchronize production peers or perform the Database rebuild.

Please include distributed-ingestion and replication correctness explicitly in your Phase 1 verification report. If additional implementation is necessary, integrate it into your existing plan rather than starting a separate project.

---

# Subsequent clarification: structural file retention and associations

Additional Phase 1 clarification — Structural file retention and associations

Please explicitly incorporate the following requirement into the ongoing Database implementation.

For every curated molecular record, retain or associate its actual structural files **whenever available and scientifically valid**.

Required representations:

- Proteins and targets: original PDB/mmCIF structures, supplied templates, and separately identified predicted structures.
- Designed binders: individual binder structures and complete target–binder complex structures.
- Small molecules: original chemical identity and coordinate-bearing SDF when available.
- Docking: original ligand SDF, prepared ligand state, docked ligand pose (preferably exportable as SDF), exact receptor PDB/mmCIF, and original PDBQT where relevant.
- Protein–ligand complexes: complete bound-complex PDB/mmCIF, linked to the exact receptor and ligand pose.
- Historical imports: preserve original source files, exact checksums, chain/atom mappings, model/pose identifiers and provenance.

The Database must distinguish original supplied structures, experimental structures, prepared structures, predicted structures and reconstructed/exported complexes.

A ligand, its different poses, the protein receptor and the bound complex must maintain correct relationships. Different docking poses must not become duplicate compound identities.

Reuse Studio's existing bound-docking and structure-export infrastructure. Do not create a new structure storage system.

When a ligand pose is stored as PDBQT, support validated reconstruction/export to SDF using the original chemical identity and atom mapping. Reject ambiguous conversions rather than silently creating chemically incorrect structures.

Retain selected scientifically valuable structures in the curated Database. Keep unselected high-volume docking poses and intermediate structures in Results.

For the distributed architecture, structure files must remain retrievable from their owning workstation through authorized content-addressed storage. Metadata synchronization alone is insufficient if a saved structure cannot subsequently be opened or exported. Preserve checksums and distinguish locally available from remotely retrievable or unavailable files.

Do not invent structures when no valid coordinates exist. Sequence-only proteins are legitimate molecular records.

Verify these requirements for protein design, cofolding, docking, historical BindCraft imports and existing small-molecule results using isolated tests.

In the Phase 1 report, distinguish what is already implemented from what was corrected and what remains unsupported. Do not modify production data or initiate Phase 2.
