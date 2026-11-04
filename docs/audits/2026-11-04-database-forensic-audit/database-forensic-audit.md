# MN Molecular Design Studio — Forensic Database Audit and Rebuild Feasibility

Audit date: 2026-11-04. Production data was inspected read-only. No records, schemas, services, workflows, or peer state were changed.

## Executive summary

The audit finds **a usable core data model with an overly broad and partly implicit ingestion boundary**. Molecular identity, structures, observations, and collections are distinct concepts in the current SQLite schema; the worst volume is from automatic retention of job inputs and outputs, not from a need to redesign scientific identity from scratch.

1. **Is the Database fundamentally too complicated?** The current logical model is not fundamentally too complicated. It already has separate molecule identities, structure representations, raw observations, collections, source links, and Results provenance. Its implementation has accumulated complexity around heuristic run classification, legacy fallback types, derived evaluation tables, and three-peer replication history. Those are real maintenance costs, but do not yet justify a second schema.
2. **Are identities and Evaluations generally correct?** They are generally well separated and the direct local consistency checks found no orphan memberships, structures, observations, sources, or evaluation projections. Canonical protein identity is composition/sequence based; small molecules use RDKit isomeric canonical SMILES. One cross-run compound ID (`mol-45cb9598df1fec323a416015bbecefc7`) appears in 15 runs and three docking methods, as intended. Exact Complex Analysis associations were verified for the audited 537 accepted R3/Tau structures. This is not a blanket scientific validation of all 21,537 structure rows; 2,230 structures lack a simple representation-type key.
3. **What appears scientifically valuable under the proposed rule?** The first-pass reusable-identity set is **3,848 distinct identities (17.9% of 21,473)**: 3,080 curated, 43 reference, and 725 explicitly accepted generated identities. These groups are disjoint in the current membership query. This estimate counts identities, not the scientific value of every Evaluation. Selected observations and exact structures should also be retained for those identities and other explicitly chosen studies.
4. **How much clutter is automatic?** At least 15,000 input compounds from three 5,000-item UniDock Pro screens are recorded as persistent workflow inputs, alongside three receptor entries. That is about 70% of all current molecule identities, and roughly 67% of membership rows. Of 1,050 groups, 920 have job provenance; 377 groups are empty (286 generated and 91 intermediate under the application’s fallback rules). All 15,000 compounds may be useful screening inputs, but they were not an explicitly curated hit selection. The current rule keeps them in the Database by default.
5. **What is Database-only?** No saved manual Complex Analysis annotations or complex-reference-context rows were present. However, 3,080 curated identities, exact collection memberships, imported-row audit metadata, source-only proteins, and policy evidence may encode user decisions not recoverable from Results alone. The exact count of facts that exist only in the Database is **unknown** because every source artifact was not exhaustively matched against all local and remote Results stores.
6. **Can the Database be safely reconstructed from Results and source structures?** A selected index is feasible, but a full replacement is not yet proven safe. The independent November 3 audit reconstructed and checked exact target/structure relationships for 537 accepted R3/Tau binders and checked rejected-membership cleanup on all three nodes. It also identified historical runs without recoverable ledgers. Imported records, user collections, sequence-source links, raw observations, and replica ownership need an explicit migration manifest.
7. **What would recreate clutter after a rebuild?** Default `automatic` inclusion, completion-hook ingestion plus workspace scans, input auto-registration, evaluation/intermediate group creation, weak empty-group prevention, substring run classification, and missing acceptance status. Rebuild without changing these contracts would repopulate the same class of records.
8. **Repair, rebuild, or redesign?** Evidence supports **Option B after first correcting the Option A ingestion rules**: keep the existing logical model, create a clean index only after inclusion policy and group-creation behavior are fixed, and migrate curated/reference/accepted/policy-specific information. Option C is not justified by the evidence found so far.
9. **Safest sequence if authorized later?** Preserve all current databases and replica catalogs; freeze explicit keep rules; checksum and inventory all Results/source artifacts; construct a candidate index offline; validate identity, policy, representation, metrics, and provenance crosswalks; rehearse rollback; then cut over from one owner manifest and one replication epoch across the three nodes. Do not delete old stores until convergence and review are complete.
10. **What remains uncertain?** The exact 128-observation/source-link peer discrepancy, physical database and structure-store sizes on PC004032/PC004033, complete source artifact coverage, exact classification of 2,230 structures, per-residue source-PDB missing-residue mappings, and the full set of Results that can still be read from each owner workstation.

## Scope, methods, and limits

The working directory was `main` at commit `d5df4d6` (“Add interactive analysis dashboard and reusable protein sources”). Five files had pre-existing uncommitted changes: the ingest file, Complex Analysis dashboard/view, and two test files. They were preserved. In particular, the local uncommitted ingest diff adds `validation_acceptance` run classification and defers accepted designs without coordinates; it was not assumed to be running on production peers. Peer runtime code fingerprints differed from the local checkout.

The primary database was read through SQLite URI `mode=ro`, at `/home/user/.local/share/mn-molecular-design-studio/database/molecular_library.sqlite3`. I used `PRAGMA user_version`, `sqlite_master`, `PRAGMA table_info`, direct aggregate queries, and JSON metadata extraction. I did **not** open it through `MolecularDatabase.connection()`: that routine executes `CREATE TABLE/INDEX/TRIGGER IF NOT EXISTS` and can modify schema as a side effect. The local database’s latest group update was `2026-11-04T12:16:55Z`; the reported counts are a dated snapshot, not a live guarantee.

PC004032 and PC004033 were inspected using authenticated, read-only paginated peer-record exports, 500 records per page (1,322 pages each). The remote API exposes molecular metadata, not SQLite PRAGMA output, filesystem statistics, or a full remote content-store listing. Therefore remote physical DB size, schema version, and local structure-store size are **unavailable**, not zero. Direct SSH inspection was unavailable; no peer synchronization or cleanup was triggered by this audit. The full snapshots show the two peer views agree on current live molecular entry/group counts, while differences against the primary remain as listed below.

The group-category counts in this report use the same rule as `MolecularDatabase._collection_kind()`: use persisted `origin.collection_kind`; otherwise classify old job groups as intermediate and old non-job groups as curated. Empty-group counts use direct group/member joins. Source-owned replica counts are reported separately from materialized views and are not summed as independent molecule libraries.

The only source-Results audit relied on here is the independent [2026-11-03 Complex Analysis audit](../2026-11-03-complex-analysis/README.md) and its CSV/JSON evidence. It inventoried 2,410 deduplicated source runs; 591 runs had no recoverable target sequence in that audit. It evaluated 192 runs for 780 pending candidate/run pairs. This does not prove all unrelated historical jobs or external imports are recoverable.

## Actual contents and three-workstation status

### Primary local database — pc002975

| Measure | Snapshot |
|---|---:|
| Database file | 1,715,109,888 bytes (about 1.60 GiB) |
| `molecular_files` store | 1,947,575,687 bytes; 35,312 files |
| Other reference/import source stores | 1,646,710 bytes; 6 files |
| Tables | 15; `PRAGMA user_version=0`; no explicit migration/version table found |
| Molecule identities | 21,473: 17,825 small molecule; 2,814 protein; 832 complex; 1 DNA; 1 RNA |
| Structure rows | 21,537: 12,936 PDBQT; 7,727 CIF; 854 PDB; 20 gzip |
| Observations / evaluation projections / metric index | 41,198 / 41,198 / 746,644 |
| Provenance links / ingestion signatures | 33,604 / 1,707 |
| Groups / memberships / empty groups | 1,050 / 22,306 / 377 |

Current group categories are 132 curated, 577 generated, 16 evaluation, 20 inputs, 287 intermediate (including legacy fallback), and 18 reference. Membership rows by category are curated 3,117 / 3,080 distinct identities; generated 1,779 / 1,199; evaluation 1,287 / 1,179; inputs 15,794 / 15,782; intermediate 278 / 121; and reference 51 / 43. Distinct-ID counts overlap between categories except for the explicit first-pass retention union above.

The 15,000-item input-screen footprint is concentrated in three groups of 5,000, 5,001, and 5,002 members. The source input data shows about 15,000 small molecules and three receptor entries. These are not 15,000 selected hits. Direct input provenance links identify UniDock Pro as the largest input workflow (15,001 identities in source context, including receptor/input associations). A current observation namespace count includes 15,001 UniDock Pro observations, 3,693 UDP, 3,626 Vina, 3,623 GNINA, 2,397 AlphaFold 3, 2,394 Boltz-2, 2,377 Nesso-1, 1,110 `validation_acceptance`, 876 RFD4 Proteina, 837 ipSAE, and 544 Biotite interface observations. These namespaces are preserved as separate methods; aliases and different output types mean they should not be combined as a universal score.

All 18,274 distinct non-remote structure path references existed during inspection: 18,098 pointed inside `molecular_files`, and 176 pointed elsewhere (primarily Results-owned locations). There were 139 `cluster-library://` reference rows. The store contains 17,214 files (417,294,177 bytes) not directly named by a current structure path. This is only a **candidate orphan-file** count: files may be shared, non-structure artifacts, or stale blobs, so none should be deleted on this evidence. Of 21,537 structure rows, 2,230 had no `representation_kind`, `structure_type`, `source_type`, or `tool` key in their context. Local relational checks found zero missing molecules/groups for memberships, structures, observations, and sources; zero source links without a matching same-group membership; and zero missing observation/evaluation counterparts.

### Peer materialization and source ownership

| View or source owner | Groups | Memberships | Distinct source IDs | Structures | Observations | Source links | Physical file/schema data |
|---|---:|---:|---:|---:|---:|---:|---|
| pc002975 direct local DB | 1,050 | 22,306 | 21,473 | 21,537 | 41,198 | 33,604 | Observed directly; schema version 0 |
| PC004032 peer-visible materialization | 1,050 | 22,303 | not exposed | 21,534 | 41,070 | 33,476 | Not exposed by read-only metadata API |
| PC004033 peer-visible materialization | 1,050 | 22,303 | not exposed | 21,534 | 41,070 | 33,476 | Not exposed by read-only metadata API |
| pc002975 authoritative source records | 684 | 10,196 | 9,379 | 6,802 | 11,342 | 13,183 | By source author; not an additive unique cluster count |
| PC004032 authoritative source records | 209 | 5,347 | 5,334 | 430 | 5,357 | 5,435 | By source author; not an additive unique cluster count |
| PC004033 authoritative source records | 157 | 6,760 | 6,717 | 14,302 | 24,371 | 14,858 | By source author; not an additive unique cluster count |

Source-author type counts were: pc002975 7,819 small molecules, 2,309 proteins, 66 complexes, one DNA and one RNA; PC004032 5,001 small molecules and 346 proteins; PC004033 5,838 small molecules, 122 proteins, and 800 complexes. Author-distinct IDs may overlap across workstations and should not be added to estimate a global unique total.

The confirmed current gap between the primary direct DB and peer-export aggregation is 3 memberships, 3 structure rows, 128 observations, and 128 source links. The peer exports agree with each other. The snapshot does not isolate the missing relation rows, so the root cause is **unresolved**, not assigned to any one peer. The peer job catalog also differs: PC004032 exported 6,186 live job rows, PC004033 5,564, and the primary local snapshot 6,212. The local replication catalog contains one conflict.

The local `replication.sqlite3` file is 1,179,967,488 bytes and has 653,406 molecular-entry keys, of which 631,100 (96.6%) are tombstones. Per author, entry-key tombstones were pc002975 120,422; PC004032 210,387; PC004033 300,291. Retaining those tombstones is necessary until a coordinated epoch/reconciliation plan exists. A rebuild of the SQLite library alone would not reset or correct these replicated ownership records.

Direct local provenance attribution found 12,381 source links for pc002975, 5,486 for PC004032, 14,903 for PC004033, 810 without workstation attribution, and 24 carrying a node UUID instead of a hostname. Groups span 25 project IDs plus two unscoped groups. The supporting [inventory CSV](database-inventory.csv) includes per-project counts under SHA256 pseudonyms, and method-namespace counts by source workstation; raw project IDs are not included. Unexposed remote values are marked instead of treated as zero.

Top observation namespaces by group-source workstation were pc002975: `validation_acceptance` 840, `ipsae` 810, `rfd4_proteina` 737, `visual_pipeline` 559, `candidate_cohort` 400, and `biotite_interface` 324; PC004032: `validation_acceptance` 189, `biotite_interface` 150, `bindcraft2` 109, `rfd4_proteina` 76, and `ipsae` 9; PC004033: `AlphaFold 3` 2,397, `udp` 756, `vina` 459, `gnina` 324, `Boltz-2` 258, and `Nesso-1` 248. A large additional block has unknown workstation attribution: `unidock_pro` 15,000, `gnina` 3,298, `vina` 3,165, and `udp` 2,937 observations. Namespace spelling is preserved, so differently named variants are not merged in this audit.

## Actual lifecycle and identity contracts

The application has the requested logical concepts, but they are not all equally authoritative:

1. **Input:** a job contains sequence/structure/compound data in Results. Protein Library input collection records authorized source protein inputs before or during run completion. The general ingester later scans completed jobs and may also add inputs.
2. **Identity:** protein/RNA/DNA identities use component type and uppercase sequence. Small molecules use RDKit canonical isomeric SMILES. A molecule's fingerprint is a sorted component multiset; chain labels and target/binder roles are stored as context/member metadata but do not define identity. Repeated identical chains/components remain part of the composition.
3. **Representations:** PDB/CIF/PDBQT and predictions/complexes/poses are separate `molecule_structures` attached to an identity and group with source context. Content-addressed local files are deduplicated by blob path/digest; a peer may instead have a `cluster-library://` reference and retrieve/cache it from the owner.
4. **Generation and prediction:** `_run_mode()` identifies generation/evaluation from labels and substring matches. The run's result candidates are normalized; generated protein binders can be represented by binder identity while a target-binder complex remains its structure representation.
5. **Acceptance:** raw acceptance policy and decision metadata are stored on member/source records. Rejected and unresolved decisions remain in Results; current generated memberships contain no explicit rejected rows. “No explicit status” remains a separate state.
6. **Database ingestion:** a completion hook invokes the Protein Library hook and forced ingestion; GET `/groups` and `/metrics` trigger workspace scans. `molecular_ingestions` signatures help avoid unchanged replay. This combination of hooks is redundant but not always duplicative because signatures and stable source keys provide idempotency.
7. **Evaluation:** raw `molecule_observations` store engine namespace, run, metric payload, and context. `molecule_evaluations` and `molecule_metric_values` are query projections; `evaluation.v5:studio.output.v1` is set on all current projection rows.
8. **Collection and replication:** groups/memberships express reusable sets; Results remains the complete execution archive. Peer replication uses owner-attributed group/entry metadata and durable tombstones; structures are retrieved from owners or cached when marked portable.

Canonical SMILES preserves explicit stereochemistry and does not neutralize, desalt, or collapse tautomers. Protonation/tautomer forms that canonicalize differently remain separate identities, but there is no general parent/variant relation that unifies them. That is a representation choice rather than a demonstrated defect in the sampled screening data. Different poses should remain structures or evaluations; the identity model supports this, and the sampled cross-method compound confirms reuse.

## Confirmed ingestion and scientific-data issues

The detailed path-by-path mechanisms and smallest corrections are in [ingestion-root-causes.md](ingestion-root-causes.md). The material findings are:

- `collection_policy` supports `automatic`, `manual`, and `none` with parent inheritance, but unmarked jobs resolve to `automatic`. A rebuild without a policy change recreates the current input/output accumulation.
- Input molecules not already associated with an accessible source group are placed in persistent `inputs` groups. This preserves reusable inputs but has no screen-size or “save hits only” rule.
- The current Database has 377 empty groups, despite UI filtering them from ordinary group listings. The group records continue to occupy DB and peer metadata.
- 1,231 generated memberships / 725 distinct identities have explicit accepted metadata; 31 are unresolved; 517 memberships / 468 identities have no explicit acceptance status. Those unknowns must not be classified as accepted.
- `complex_review_annotations` and `complex_reference_contexts` are both empty in the current DB. There is no human-review note corpus in those two tables to lose during a migration, but this does not establish that no informal notes exist elsewhere.
- 4 imported small-molecule batches contain 807 rows, of which 801 are valid and 6 invalid; the valid rows link to 792 distinct imported identities. Keep the import mappings, duplicate/stereochemical-variant decisions, warnings, and source files.
- 17,214 content-store files are not directly referenced by structure paths. They are a forensic review set, not safe cleanup targets.

The independent November 3 Complex Analysis audit provides a useful cohort-level validation: 537 accepted R3/Tau structures have exact saved molecule/structure/target associations and Biotite interface geometry. Target chain A and binder chain B were used for 537/537. ipSAE was measured for 125 BindCraft2 structures with compatible saved PAE. It was unavailable—not failed—for 412 RFD4 structures without compatible pairwise PAE. These measurements do not modify acceptance policy.

The source audit evaluated 780 pending RFD4 candidate/run pairs from 192 runs: 35 accepted, 727 rejected, 18 unresolved. The 35 newly accepted are 29 R3 and 6 Tau. Together with existing accepted structures, the audited cohort totals 406 R3 and 131 Tau: R3 BindCraft2 100, R3 RFD4 306, Tau BindCraft2 25, Tau RFD4 106. An additional target-identity check found 636 historically R3-labeled RFD4 identities split across 286 actual R3-target and 350 actual Tau-target records. The exact target IDs/sequence hashes are more reliable than campaign names. The separate PD-L1 saved-policy cohort records 162 accepted identities out of 300 candidates; the remaining 138 were rejected or deferred under that policy, not all affirmatively rejected. A native PPI diagnostic reported 118 pass and 44 fail among the 162 accepted, but all 162 remained accepted by the separate saved policy. These tests must not be merged or treated as equivalent acceptance policies. The audit verified removal of 727 newly rejected and 248 prior rejected memberships from all three source views while preserving Results; durable tombstones were used.

## Representative identity and structure checks

| Example | Direct evidence | Audit interpretation |
|---|---|---|
| R3 sequence-defined target | `mol-d2b7e1fd31ce750c4b75fe0b20119197`, one 31-aa protein component | Reused target identity; role belongs in run context. |
| Tau sequence-defined target | `mol-5b8294fd0a0b697c318184f64286a427`, one 70-aa protein component | Distinct from R3 even where historical campaign labels contain both names. |
| PD-L1 target | `mol-16a1556642b28bf0c9bee6739bf37cb4`, one 115-aa protein component | Structured target identity appears in RFD4 source runs; its 162 accepted designs retain a separate saved policy. |
| R3 BindCraft2 | `mol-0305630ebda977f02dcd95b74de61beb`, run `20261019-181522-21f9a96c`, source PC004032; exact analysis target reference `target-842ce3a9a69d510e4e0cbead7903afef`; chains A/B; ipSAE 0.209636 | Representative exact structure/analysis association from the 537-structure audit. |
| R3 RFD4 | `mol-57cb7d558e0f1d88ac59eb683ef2336c`, structure `structure-d6cc59901b09530590581429726a6be9`, source run `20261029-210422-9af4aeaa`, group `accepted-review-cd4bac69933895474d40272de91029c9`, target `mol-d2b7e1fd31ce750c4b75fe0b20119197` | Exact identity, structure digest, and target ID are preserved; one of the newly accepted RFD4 examples. |
| Tau RFD4 | `mol-1904dbb2787d38d96e68a50c62e24fac`, candidate `rfd4_000045_idr_peptide_binder_rep49`, source run `20261031-112821-9b884998` | Audit decision is tied to a specific candidate/run and Tau sequence target, not just a campaign label. |
| PD-L1 saved policy | group `job-32ac2af5edecceead7dd5d01b2fa0292`, run `20261102-081647-0ff0921a`, 162 explicit accepted member records | Kept as its own policy decision; native PPI filter evidence is separate. |
| 7QF4 sequence source | reference group `protein-library-f666120b02c99aca3f78f7aa151ec707`; sequence identity `mol-2bd2da0a8c340938d54f67c8e48dba64` has one 128-aa component and no coordinates in that group; source complex `mol-8b638e2ff35bd09b8fd25342cf26b985` is linked as `sequence_source_reference` and has 44 structure rows | The source PDB/system and matched reusable sequence identity are distinct. This prevents treating a sequence-only target as if its original full coordinates were the same representation. |
| 8JYR sequence source | reference group `protein-library-398a94501208de2b0b0f8a399241e74a`; source complex `mol-13f31a4d276a7e62bf017d993a1ebc7c` has 21 structures; protein components are also represented by IDs `mol-edbc300706a703d1dcb9789eeb2d7201` and `mol-ebbebc2dacd6423b9d5e73e61acd75ba` | The shared sequence-source relationship is retained alongside source structures and individual protein component identities. |
| Multi-run docking compound | `mol-45cb9598df1fec323a416015bbecefc7`; 15 runs; `gnina`, `udp`, `vina` | Repeated scoring reused a canonical compound identity. Pose/prepared structures remain separate representations. |

The current Protein Library documentation states the correct high-level contract: sequence/composition is identity; experimental, supplied, and predicted coordinates are distinct; sequence-only IDR targets do not receive invented templates; and PDB sequence matching requires full-chain composition/copy-count matches. The audit confirms those relationships for the examples above. Per-residue missing-residue mappings are not uniformly evidenced as a first-class field in the records reviewed; the original source coordinate files and sequence digest should be preserved so this can be checked during migration.

## Scientific retention and reconstruction

Keep raw, method-specific observations for retained designs, source proteins, compounds, and selected screening hits. An Evaluation should identify the exact molecule, target, structure digest/model, evaluator/version, scientific parameters, and Results run. The current primary has `rosetta_metrics` 14, `pep_patch` 4, `sequence_liability` 6, `tag_fusion_suitability` 9, and `rosetta_ligand` 7 observations. There is no namespace containing `experiment` among the 41,198 observations; this does not rule out experimental values stored under generic imported metadata. Native scores, Biotite, Rosetta, ipSAE, PEP-Patch, liability, tag/fusion, ligand pose validity, and experimental measurements remain distinct; no method-independent score should replace them.

Keep `molecule_observations` as the source record when retained. Rebuild `molecule_evaluations` and the metric-value index from those observations. Keep in Results the complete rejected trajectory, model attempts, logs, checkpoints, temporary conversions, unselected poses, large simulation trajectories, and failed outputs. Keep or migrate the exact selection and provenance for reusable collections. Keep imported source payloads and import mappings even when the compound identity is deduplicated into an existing molecule.

Reconstruction quantities and gaps are tabulated in [reconstruction-feasibility.csv](reconstruction-feasibility.csv). In brief, the audit supports high coverage for the 537 accepted R3/Tau complex-analysis cohort and 162 PD-L1 acceptance entries. The 4 imported batches can be reconciled from 807 import-row records, but all source-file checksums and paths still need a migration manifest. Current raw observations may overlap Results, but an exhaustive exact comparison is outstanding. The 517 generated membership rows without an explicit decision and 2,230 unclassified structure rows are not safe to automatically reclassify. User-selected groups and source-only records require migration even when their identities can be rediscovered elsewhere.

## Simpler long-term organization

Keep one canonical identity index. **My Proteins** and **My Compounds** should filter that index by molecule type. **Structures** should filter representations by experimental/supplied/predicted/complex/pose kind. **Reference Structures** should be a source/reference view with accession and chain correspondence. **Collections** should be explicitly curated memberships. **My Molecules** can provide the general identity view. These should mostly be views over shared storage, not independent databases.

Protein Library should become the reference/source collection view inside the general molecular library. It should share sequence identities with generated designs and evaluations, but preserve its special source accession, source file, chain, and structure provenance. My Results should continue to own the full run history and artifacts.

## Recommendation and safe sequence

Use Option B only after correcting intake rules. Option A changes are mandatory before rebuilding because current defaults and hooks would recreate screen inputs and group clutter. Keep the existing logical schema for now; only consider explicit schema versioning and a formal migration once an offline candidate index reveals a data concept that the current model cannot represent without ambiguity. The trade-offs are detailed in [architecture-options.md](architecture-options.md).

If a rebuild is later authorized, preserve every current DB, file store, peer catalog, and Results source; do not drop tombstones or overwrite the source index. Build an offline candidate from a signed keep-rule manifest. Reconcile identities, exact structure hashes, target/binder chains, source PDB relationships, import batches, policy-specific acceptances, raw observations, and source-owner links. Prove restoration of old snapshots. Then replace only a single isolated index, verify all checks, and roll out one consistent source-author manifest/epoch to the three workstations. Keep the old databases read-only until the user accepts the reconciliation and the cluster passes a normal sync cycle.

## Supporting artifacts

- [Per-workstation and deduplicated inventory](database-inventory.csv)
- [Data classification and retention](data-classification.csv)
- [Ingestion root causes and code paths](ingestion-root-causes.md)
- [Reconstruction feasibility matrix](reconstruction-feasibility.csv)
- [Architecture options and staged sequence](architecture-options.md)

No tests, scientific workflows, schema changes, deletion, peer sync, migration, or service restarts were performed for this report.
