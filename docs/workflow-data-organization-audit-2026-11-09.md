# Workflow and scientific data organization audit — 2026-11-09

## Finding and recommendation

The concern is justified. The current collection-first Complex Analysis page does not present the researcher's experiment reliably: one campaign can have many jobs on different workstations, several generated sequences and several predicted structures per sequence. Database groups are not an adequate substitute for that hierarchy.

The application already uses relational databases (SQLite). The demonstrated problem is incomplete normalization and inconsistent presentation of relationships, together with overlapping automatic ingestion and explicit collection paths. Replacing SQLite, or creating separate scientific databases for design, folding and docking, would not by itself fix those problems.

**Recommended direction:** one coherent relational scientific catalog with explicit campaign/job, design/variant, prediction/pose, acceptance and artifact relationships. Keep the existing operational stores and scientific files; keep the current molecular, structure and Evaluation identities. Restore a campaign-oriented entry point with separate design, folding and docking views. Make user-created collections optional saved cohorts, rather than a required transfer step.

A new organization should first be built and validated in an isolated copy. This audit does not authorize or perform migration, historical reconciliation, production changes or deployment.

## Scope and evidence boundaries

- Checkout: `2c30eb10435d96317c4c319e7ce3207fbdf3882b`.
- Production UI service points to the tracked snapshot `p7-2c30eb1` on pc002975.
- Molecular Database inspected with raw SQLite connections using `mode=ro` and `PRAGMA query_only=ON`: `/home/user/.local/share/mn-molecular-design-studio/database/molecular_library.sqlite3`.
- Results root from the existing application settings: `/mnt/data/RESULTS/mn-molecular-design-studio/workdir/runs`.
- Read-only inspection of the local job catalog, existing coordinator job metadata, selected PD-L1 metadata/candidate/ledger files, and relevant current source code.
- No scientific calculations, remote directory inventories, Database writes, ingestion, peer synchronization, service changes or historical cleanup. The paused Database cleanup files were not edited.
- Counts below are a primary-workstation snapshot, including its existing peer replicas. This is not an independent completeness check of every scientific file on all three workstations. Missing metadata does not prove missing coordinates; an existing path does not establish the recoverability of an entire campaign.

## What exists today

| Responsibility | Existing storage/implementation | Current boundary |
| --- | --- | --- |
| Original inputs and outputs | Results job directories, metadata/input/result JSON, normalized candidate JSONL, native tables and manifests | Durable source evidence, including data not collected into the library |
| Job search and execution hierarchy | `database/job_catalog.sqlite3`, coordinator metadata, `logical_runs.py`, distributed campaign ledgers | Can associate pipelines, controllers, batches and workers; some legacy associations use qualified name grouping |
| Molecular identity | `molecules` | Canonical component/sequence/chemistry identity; not one row per generation attempt |
| Membership and access | `molecule_groups`, `molecule_members` | A group can be a generated job, curated cohort, reference, input, intermediate or evaluation group |
| Structures | `molecule_structures`, persistent structure files and authorized peer references | A structure row includes a source group/context; row count is not necessarily independent-model count |
| Measurements and Evaluations | `molecule_observations`, `molecule_evaluations`, metric-value projection | Existing persistent evidence and indexed descriptors, with method and scientific subject context |
| Scientific source links | `molecule_sources.provenance`, member metadata, structure context | Much candidate, target, campaign and model information is embedded JSON rather than queryable relationships |
| Human shortlist/review | `complex_review_annotations`, existing peer synchronization | Separate from computational acceptance; preserve this mechanism |
| Peer state and security | Existing replication store, coordinator and access database | Operational state should not be replaced by a scientific-data reorganization |

Having separate files for operational responsibilities is reasonable. Requiring researchers to understand those files and move between their representations is the usability problem.

### Current snapshot

| Stored object | Count | Meaning/limitation |
| --- | ---: | --- |
| Canonical molecules | 21,473 | Includes 17,825 small molecules, 2,814 proteins, 832 complexes, one DNA and one RNA; not a binder-design count |
| Groups | 1,050 | Mixed purposes; 577 generated, 128 curated, 43 intermediate, 20 input, 18 reference, 16 evaluation, 248 without a recorded collection kind |
| Memberships | 22,306 | A molecular identity can be in multiple groups |
| Structure records | 21,537 | Includes repeated source/context references and several representation types |
| Sources | 33,604 | Includes design, folding, docking, analysis and campaign provenance |
| Observations / Evaluation projection rows | 41,198 / 41,198 | These are existing evidence records and their projection, not two independent result populations |
| Explicit `results_design_collection` groups | 0 | The new collection action has not produced persisted groups of this type in this inspected production database |

Of the 1,050 group origins, 1,001 contain a source run ID, but only 14 contain a top-level campaign/parent-campaign ID and 122 a campaign name. None contains the new top-level `generated_candidate_count` field. These are coverage counts for those exact fields, not assertions that all equivalent information is absent elsewhere.

Source-workstation labels occur for pc002975, PC004032 and PC004033. Some labels are absent or differ in casing; use node identities for joins, with display names as labels.

Historical explicit top-level `model_id`/`model_index`/`model` fields are absent from the inspected structure contexts. Some model evidence may remain in nested candidates, source artifacts and filenames. Current ingestion code does preserve parsed model selectors for new attachments. Historical model identity must be recovered from validated evidence or remain unknown, not invented during migration.

## Demonstrated problems

### 1. A campaign is not a molecular collection

The PD-L1 BindCraft2 example has this evidence-backed execution chain:

1. Parent visual pipeline: pc002975, `campaign/20261029-094447-d4f389fb`.
2. Campaign controller: pc002975, `design/20261029-094457-fb2fd80a`; its metadata explicitly names the parent and distributed campaign `d95dfc71-fe3e-4ca7-91a7-d2f1257b0ceb`.
3. The controller ledger's attempt records cluster job `8a82d216-bd5a-4768-b3dd-830120a04c49`.
4. Existing coordinator metadata maps that exact cluster job to PC004033 and `design/20261029-094502-6ad69fc5`.

The controller also records an execution/proxy source `runs:///design/20261029-094458-5f7061f2`. This is another execution identity to retain, not a reason to substitute or merge jobs by filename.

The Database has two generated groups: the controller group and the PC004033 worker group. Both point to canonical molecule `mol-6789314c240186b94a6ba805a84cc7a6`, but their group origins do not carry this complete parent hierarchy. The current collection picker therefore cannot produce the correct campaign view solely from its group rows.

The worker has an empty campaign ID in its coordinator payload. The ledger provides the explicit missing edge. Grouping only on nonempty child `campaign_id` values would still be insufficient.

### 2. Structure references are not independent predictions

The two BindCraft2 groups contain three structure records: one worker representation and two controller representations. All three record digest:

`a6ebc954057e53aa4907f277c13bec13e22fcf324a5272762a03b8669f3fea9b`

The controller's two files are available locally. Their contexts differ in the chain-role information they retain. These three references must not be presented as three independent prediction outputs. Conversely, equal coordinates alone do not justify deleting source records or merging scientific executions.

Separate the coordinate artifact, coordinate-model selector, scientific model/pose occurrence and its source associations. Preserve every original association.

### 3. The useful run overview was removed from the entry path

The earlier dashboard's `_result_sources()` and `_table_rows()` used `logical_runs` to show logical runs, participating jobs, workstations, engine names, structures and recorded acceptance counts. That matches the researcher's mental model much better than listing each Database group.

The current Database-only selector obtains library groups and `collection_statistics()`. It excludes groups without members or structures and lists the remainder separately. It does not join them to the logical-run hierarchy. Its counts therefore describe saved molecular membership, not campaign production.

Restoring the old grouping concept is appropriate. Restoring its synchronous whole-catalog paging and local candidate-file reads on every landing-page load is not necessary. Persist/index compact run summaries and query those with pagination; retrieve exact artifacts only on demand.

### 4. “Collect” overlaps existing automatic ingestion

`collection_policy.py` and completion ingestion already support automatic collection for known workflow purposes. Generated designs require a recorded applicable acceptance pass in the ingestion path; policy-less historical jobs remain Results-only. Backbones, predictions, analyses, docking outputs and curated imports have different handling.

The newer Results action creates an explicit curated snapshot of passing candidate/model rows. It can preserve a deliberate selection, but it creates another group even when accepted outputs already have generated-group memberships. There are currently no groups of the new snapshot type in this production database, yet the landing page tells users to start from one.

A separate collection action should not be mandatory just to open an already retained passing campaign. Successful ingestion and promotion need visible receipts/counts. A saved collection should mean a deliberate cross-run cohort or selection, not a repair step users must discover.

Automatic ingestion may fail while a scientific job remains completed: `record_completed_run()` logs indexing failures without changing scientific completion. That is an appropriate distinction, but the UI must show “execution completed; Database indexing incomplete,” rather than making the output appear absent.

### 5. Run totals and acceptance scopes cannot be inferred from retained rows

In the PD-L1 controller, saved metrics report one generated normalized candidate and one native accepted candidate, with acceptance level `native` and zero validated accepted candidates. The native evidence separately records eight scored rows and five rejected rows. These are producer-specific populations/stages; they must not be added or treated as a single total of unique designs.

The current acceptance inventory counts normalized candidate/model rows in one run. That is useful for selecting those rows, but is not necessarily the total number of native trajectories, sequence variants, attempted predictions or scored compounds.

Required summary fields need a counting unit, stage, source, completeness and authoritative aggregation scope. Final pass, native pass, unresolved and no decision stay distinct. Never sum parent aggregates and their worker outputs twice.

### 6. Database output retention differs by scientific operation

Current prediction ingestion selects `best_candidates_by_complex()` for automatic structure retention, labels those structures `best_prediction`, and prunes stale indexed prediction structures for the affected candidates. Original outputs remain a separate Results responsibility. The Database snapshot has 934 rows with `best_prediction` policy.

Current docking retention includes `best_docking_pose` (12,949 stored structure rows in this snapshot), receptor representations and a small number of bound-complex records. Scores may be associated with canonical compounds while a pose structure has its own receptor context.

Therefore “all models/poses in the run” and “models/poses retained in the Database” are different sets. The UI must display the difference explicitly and permit deliberate additional-model retention where supported.

### 7. Unavailable structures can disappear from the collection list

The R3 Database group for `20261030-014030-a12be5a3` has one member, zero structures and `source_run_available=false`. The current picker omits such groups. This is a known Database availability gap, not proof that the entire historical Results directory has been lost.

Other groups without structures include intentional sequence-only/reference/input records, so absence of coordinates is not automatically a data-loss finding. Availability needs an explicit status rather than silently dropping the experiment from the overview.

## Required scientific distinctions

| Research operation | Identity and relationships needed | What the researcher should see |
| --- | --- | --- |
| Binder/design campaign | Campaign → jobs/attempts → generated design → sequence variants | Generated/native-scored counts by stage; passing, rejected, unresolved; retained models and analysis coverage |
| Design-associated prediction/cofolding | Exact sequence variant + target construct/components → predictor execution → separately identified prediction outputs | Original generator and predictor separately; generated complex versus independent cofolded model; exact chain roles, sample/model and digest |
| Standalone folding/cofolding | Explicit input molecular assembly → prediction run → models | Valid prediction outputs and model confidence; no invented binder-design or final-acceptance status |
| Docking/screening | Receptor construct/model + prepared ligand variant → campaign/round/job → pose → score/evidence | Unique compounds, scored cases and poses as separate counts; chemistry/atom mapping; original score and pose availability |
| Post-design assessment | Exact design/model/pose + method/version/settings → existing Evaluation | Compatible coverage, results, failures and residue evidence; reused Evaluations remain distinguishable from new execution |
| Human selection | Existing design/model/structure + researcher decision and rationale | Experimental shortlist/exclusion, separate from computational acceptance and experimental outcome |

A binder's cofolded structure must remain directly reachable from its design. Running another predictor does not generate a new binder identity. Sequence redesign creates an explicit child variant; an identical sequence can reuse canonical molecular identity while preserving the new design occurrence and source.

For standalone folding and docking, a valid saved output is not equivalent to an accepted binder. Library admission and scientific acceptance must be recorded with their appropriate operation-specific policy.

## Proposed organization to validate in staging

This is a normalized extension/mapping of existing identities, not a competing identity system or execution engine.

```text
Project
  Experiment / logical campaign (design | folding | docking | assessment)
    Existing jobs, attempts, rounds and worker receipts
    Exact targets / input assemblies / receptor and preparation context
    Design occurrences -> canonical molecule -> sequence variants
      Prediction outputs / coordinate model selections
    Compound occurrences -> preparation variants -> docking poses
      Existing artifact references + digests + source locations
      Recorded acceptance decisions, scoped to their design/model/pose
      Existing Evaluations and their exact scientific subjects
      Existing human review / shortlist annotations

Optional saved cohort: explicit references to any compatible chosen subjects
```

Important implementation boundaries:

- Persist explicit experiment-to-job and source-to-output edges from existing metadata, ledgers and receipts. Preserve node-qualified job identity; do not use a workstation name or a physical SQLite row ID as scientific identity.
- Separate canonical molecular identity from a design occurrence in a campaign. A deduplicated sequence can occur in several campaigns without losing either source.
- Retain native IDs and existing normalized candidate IDs as source-qualified aliases. Preserve the existing CandidateSet backbone/variant/prediction lineage where present.
- Distinguish predictor output/sample identity from a PDB/mmCIF coordinate-model selector. Store both when known; unknown remains unknown.
- Preserve existing structure and Evaluation IDs. Add relationships to those objects, rather than reissuing scientific identities.
- Store acceptance decisions with policy, stage, thresholds, evidence and exact subject. Preserve multiple historical decisions and their original scopes.
- Keep all available raw history in Results; index its metadata and counts even when rejected coordinates are not admitted to the reusable molecular library. Do not bulk-promote failed or unfiltered designs.
- Keep coordinate and residue-level files in existing content-addressed/Results storage with authorized `runs:///` and `cluster://` references. A metadata reorganization does not require copying every file or creating another transport.
- A digest establishes content identity, not ownership or permission. Access and project checks apply to every source association and retrieval.
- Extend existing peer synchronization only for required new relationship metadata, preserving its identity/conflict/tombstone policy. Do not replace the scheduler, replication cursors, access store or shortlist system.

### User workflow

1. Open a project and choose **Designs**, **Folding**, or **Docking**.
2. See one row per logical experiment/campaign, aggregating its explicit jobs across workstations. Show target, methods, date, job count, scientific counts by stage, retained structures, analysis coverage and availability warnings.
3. Open the run to see its stages, sequence variants and exact models/poses. A passing design cohort opens directly when already retained. If its outputs have not been indexed/retained, show the precise action and missing count in that run.
4. Use **Analyze → Compare/Visualize → Shortlist** on the selected exact subjects. Cross-run selections can be saved as an optional cohort.

This preserves the useful old campaign table while avoiding slow live filesystem/remote discovery during Complex Analysis rendering.

## Database-engine options

| Option | Assessment |
| --- | --- |
| Separate scientific databases for designs, folding and docking | Not recommended for this problem. It separates strongly related objects and adds cross-database provenance, access and synchronization work. Use typed tables/views within a coherent scientific catalog instead. |
| Normalize the existing SQLite scientific catalog | Recommended first implementation target. It reuses current persistence and services while making missing relations explicit. Prove it in staging before choosing a production migration. |
| Central client/server relational database, such as PostgreSQL | A separate operational decision if a central always-online writer, transaction model or demonstrated concurrent-write workload is required. This audit establishes no engine-level need to change it. Changing engines without repairing relationships would preserve the current confusion. |

SQLite's own guidance distinguishes local/application-server use from many direct network writers and high write concurrency. Do not share one writable SQLite file between workstations as a substitute for the existing authorized API/replication path. [SQLite deployment guidance](https://www.sqlite.org/whentouse.html).

## Safe legacy mapping and migration approach

1. Define the entity/counting rules above and expected campaign screens before changing production. Build a versioned relational mapping in an isolated copy using existing IDs.
2. Read existing Database records, job catalogs, CandidateSets and selected manifests. Map explicit parent/child edges and aliases; do not derive scientific parentage from a name, identical sequence or filename alone.
3. Mark each mapped relationship as exact, validated legacy correspondence, ambiguous, unavailable or unsupported. Retain unknowns and original payloads. A mapping report is required for unresolved records.
4. Validate representative cases: the distributed PD-L1 BindCraft2 campaign; an RFD4 campaign with independent cofolded models; standalone folding; a docking campaign with selected and score-only poses. Preserve all source stages and distinguish aggregate/native counts.
5. Reconcile molecular IDs, structure IDs/digests, Evaluation IDs/settings, acceptance subjects, human review identities and project ACLs. Verify retries/imports are idempotent and no parent/worker double counting occurs. Do not deduplicate historical structure rows by hash alone.
6. Compare the new read-only views with existing Results and Database views, including unavailable-artifact behavior and page latency. Keep original storage and production reads intact until the mapping is accepted.
7. Only then plan an authorized additive production migration/cutover, compatible peer updates and rollback. Preserve originals and a legacy-to-current mapping. Do not repair or fabricate missing historical coordinates during migration.

This supports the suggestion to create a better structure and then map old data into it, without starting by rebuilding all scientific records or moving to a new database engine.

## Code evidence

- Logical run hierarchy and authoritative output selection: [logical_runs.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/logical_runs.py:65).
- Previous campaign-oriented summary/table: [complex_analysis_dashboard.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/complex_analysis_dashboard.py:60).
- Current group-based entry point: [complex_analysis_collection_ui.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/complex_analysis_collection_ui.py:1420).
- Existing relational schema and representation identity: [molecular_database.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/molecular_database.py:159).
- Automatic-ingestion policy: [collection_policy.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/collection_policy.py:50).
- Applicable acceptance handling: [molecular_database_ingest.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/molecular_database_ingest.py:495).
- Run-group provenance and prediction retention: [molecular_database_ingest.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/molecular_database_ingest.py:1718).
- Explicit saved selection/group creation: [molecular_database_results.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/molecular_database_results.py:612).
- Current one-run acceptance inventory: [molecular_database_results.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/molecular_database_results.py:193).
- Current collection action: [molecular_database_results_ui.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/molecular_database_results_ui.py:249).
- Existing Evaluation projection: [evaluations.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/evaluations.py:179).
- Peer metadata and file boundary: [molecular_database_peer.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/molecular_database_peer.py:1).

## Audit outcome

**The current organization needs a focused relational and workflow revision.** The prior collection-list visibility fixes help find records, but do not solve campaign membership, scientific counting or design-to-model navigation. No further picker or manual collector should be the foundation of that revision.

No application code, production records, services, replication state or historical artifacts were changed for this audit. No migration or new database was created. No new tests or calculations were required for this read-only assessment.
