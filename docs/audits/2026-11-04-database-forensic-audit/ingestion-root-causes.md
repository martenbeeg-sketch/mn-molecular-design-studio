# Ingestion root causes and confirmed failure mechanisms

This report identifies causes supported by current code or production metadata. It does not imply that the corresponding records should be deleted. Source Results remain the authoritative execution history.

## 1. Completed jobs are indexed by two paths, and browsing can trigger a scan

`record_completed_run()` records Protein Library inputs and force-ingests successful jobs from the completion hook ([jobs.py](../../../mn_molecular_design_studio/engine/workflow_engine/core/jobs.py:665), [molecular_database_ingest.py](../../../mn_molecular_design_studio/molecular_database_ingest.py:1621)). Separately, `sync_workspace()` enumerates accessible jobs and calls `ingest_run()` ([molecular_database_ingest.py](../../../mn_molecular_design_studio/molecular_database_ingest.py:1554)). The Database API invokes that scan from ordinary GET `/groups` and `/metrics` requests ([molecular_database_api.py](../../../mn_molecular_design_studio/molecular_database_api.py:81)). A signature cache makes unchanged ingests mostly idempotent, but changed candidate/structure signatures or `force=True` can replay indexing.

**Observed effect:** 1,707 ingestion signatures, 920 job-origin groups, and 15,794 input memberships. Browsing does not always create duplicates, but it is an indexing trigger rather than a read-only view. A clean rebuild followed by the same automatic policy will repopulate job-derived groups.

**Smallest correction:** Make GET queries read-only and run controlled background indexing only under an explicit per-workflow inclusion contract. Keep Results completion independent of Library promotion.

## 2. Inputs from large screens become persistent Database identities

`_find_or_create_input()` creates an `inputs` group for new identities and attaches the source to it ([molecular_database_ingest.py](../../../mn_molecular_design_studio/molecular_database_ingest.py:523)). The compound-set path applies the same behavior to each input compound ([molecular_database_ingest.py](../../../mn_molecular_design_studio/molecular_database_ingest.py:1304)). This is useful for reuse, but the policy has no scale threshold or automatic “screen Results only” default.

**Observed effect:** 15,794 input memberships and 15,782 distinct input identities. Three UniDock Pro input groups contain 5,000, 5,001, and 5,002 members, approximately 15,000 screened compounds plus three receptor entries. The database stores these as input identities and observations; this does not mean it stores 15,000 independent selected poses.

**Smallest correction:** For screening, keep the complete ligand library and all candidate rows in Results; promote only an explicit hit set (including exact selection snapshot and target context). The collection policy exists but the default remains `automatic` when no per-run value is set ([collection_policy.py](../../../mn_molecular_design_studio/collection_policy.py:10)).

## 3. Empty job groups survive after no useful membership is stored

`output_group()` constructs per-run output groups for generated/intermediate results and a shared evaluation group ([molecular_database_ingest.py](../../../mn_molecular_design_studio/molecular_database_ingest.py:958)). A run may create an output group on the first qualifying candidate, then the candidate may fail downstream parsing, be pruned, or have its last stale membership removed. The `groups()` query only omits empty groups from ordinary display; it does not delete the rows ([molecular_database.py](../../../mn_molecular_design_studio/molecular_database.py:303)).

**Observed effect:** 377 persisted empty groups: 286 `generated` and 91 legacy/fallback `intermediate`. They are hidden from the normal group list because empty non-user groups are filtered, but they remain in SQLite and replication metadata.

**Smallest correction:** Delay group creation until a valid identity/representation will be inserted, and remove a just-created group inside the same transaction if the ingest leaves it empty. Existing empty rows remain an archival/review choice.

## 4. Acceptance is not a single flag, and “no status” is not an acceptance decision

`_candidate_rejected()` rejects explicit failed/rejected/invalid outcomes and required campaigns whose validation is not accepted ([molecular_database_ingest.py](../../../mn_molecular_design_studio/molecular_database_ingest.py:338)). It does not convert missing acceptance metadata into an explicit unresolved state. `ingest_run()` can place valid generated sequences in the generated group even without a validation decision ([molecular_database_ingest.py](../../../mn_molecular_design_studio/molecular_database_ingest.py:1031)).

**Observed effect:** Generated memberships include 1,231 accepted, 31 unresolved, and 517 with no explicit `validation_acceptance.status`; the latter are 468 distinct identities. The current Database has zero memberships with an explicit rejected status after the November 3 audit's 975 rejection removals. This does not erase rejection history from Results. Acceptance-policy identity remains significant: the PD-L1 saved-policy set has 162 accepted records even though a separate native PPI diagnostic split 118 pass / 44 fail; those are distinct criteria.

**Smallest correction:** Require policy-specific accepted status for reusable generated-library placement; store unresolved and unreviewed candidates in Results or an explicitly non-library review view. Never merge status labels across policy digests.

## 5. Model classification and group provenance are partly implicit

`_run_mode()` classifies jobs by substring matches in task/tool/workflow labels ([molecular_database_ingest.py](../../../mn_molecular_design_studio/molecular_database_ingest.py:356)). Legacy group classification is also inferred at read time: old job groups become `intermediate`, while non-job groups default to `curated` ([molecular_database.py](../../../mn_molecular_design_studio/molecular_database.py:273)).

**Observed effect:** 248 current groups lack an explicit `collection_kind`; the code fallback accounts for 244 job groups as intermediate and four non-job groups as curated. This is a compatibility behavior, not proof that every legacy classification matches scientific intent. Current display counts are therefore computed using the code fallback, not a persisted normalized enum.

**Smallest correction:** Persist explicit workflow purpose and collection kind at creation; retain the original `origin` payload as provenance, but avoid repeated heuristic classification in reports/UI.

## 6. Structure role and representation evidence are not uniformly complete

Identity fingerprints include a sorted multiset of component type and sequence/canonical isomeric SMILES; chain ID and role are retained as record metadata but excluded from the fingerprint ([molecule_imports.py](../../../mn_molecular_design_studio/molecule_imports.py:25)). This appropriately lets chain roles vary across runs while keeping composition as identity. Experimental, supplied, prediction, complex, and pose coordinates remain separate structure rows with context. Complex Analysis discovery selects an exact accessible molecule, structure, group, digest, target, and chain association ([complex_analysis_sources.py](../../../mn_molecular_design_studio/complex_analysis_sources.py:208)).

**Observed effect:** The independent audit verified 537 accepted R3/Tau structures and exact target/chain links. A broader simple metadata check found 2,230 structure rows without `representation_kind`, `structure_type`, `source_type`, or `tool`; the current file paths existed, but these rows need classification before selective migration. No confirmed global target/binder reversal was found; 537/537 audited accepted structures used target chain A and binder chain B. That cohort result must not be generalized to all structures.

**Smallest correction:** Keep the identity/representation split; require representation type, exact source artifact/run, chain map, and digest when attaching a structure. Preserve original PDB/mmCIF separately from derived predictions and record chain correspondence plus missing-residue mapping where known.

## 7. Evaluation rows and metric values are projections, not independent evidence

The raw `molecule_observations` rows hold namespace, run ID, raw metrics, and context. `molecule_evaluations` and `molecule_metric_values` are typed query projections built from observations by [evaluations.py](../../../mn_molecular_design_studio/evaluations.py:87). All current 41,198 observations have one projection at `evaluation.v5:studio.output.v1`; 746,644 metric-index rows are derived values. Method-specific namespaces and native units must remain separate.

**Observed effect:** No orphan observations/projections were found in the local database. The 3 November audit verified the exact saved structures for the audited cohort and recorded ipSAE as unavailable (not failed) for 412 RFD4 structures without compatible PAE. Projection values are rebuildable; source observations and exact input identity/structure context are the records to preserve.

## 8. Replication tombstones are large, while the current live snapshot is not fully reconciled

The peer catalog uses author/counter records and durable deletion tombstones ([replication.py](../../../mn_molecular_design_studio/cluster/replication.py:218)). The local catalog has 653,406 molecular-entry keys, of which 631,100 are tombstones (96.6%). The snapshot is correct to retain tombstones: discarding them without a coordinated epoch permits deleted records to replay.

**Observed effect:** The peer export views for PC004032 and PC004033 agree with each other on current source-owned counts. Their aggregate has 22,303 memberships, 21,534 structure rows, 41,070 observations, and 33,476 source links, versus the local direct database's 22,306 / 21,537 / 41,198 / 33,604. The difference is three memberships and structure rows plus 128 observations and source links. The export alone does not establish which specific relation rows account for the 128-row difference. The replication catalogs also contain one conflict and the job catalogs differ in live record counts.

**Smallest correction:** Build a source-author keyed reconciliation report and repair delivery/replay semantics only after matching each relation to its owner. Do not use general compaction or tombstone deletion as a substitute for an epoch protocol.

## 9. Evidence that a rebuild alone would not fix

1. `automatic` remains the default collection policy; completion and scan hooks would repopulate generated, input, evaluation, and intermediate records after a rebuild.
2. 15,000 screen inputs would again be inserted unless screening selection is explicit.
3. Per-job and shared evaluation group creation can recreate empty or generic collections.
4. Missing acceptance status can still be mistaken for eligible generation output.
5. Legacy classification and incomplete representation context remain until explicit metadata is written.
6. Peer sync needs a coordinated snapshot/epoch. Replacing only the primary SQLite file cannot reconcile old tombstones or all remote source-owned records.

## Historical and working-tree context

The October 20 behavior recheck documented automatic scanning and empty group creation in an earlier dataset. The current implementation has since added signature-based idempotency, policy inheritance, exact saved Results selections, Protein Library source input retention, and exact Database analysis-source validation. The present audit confirmed the current records rather than assuming earlier findings still applied.

The inspected repository baseline is `main` at `d5df4d6`. The working tree already contained uncommitted modifications in the analysis dashboard, `molecular_database_ingest.py`, and two tests. The ingest diff changes `validation_acceptance` mode classification and defers accepted binders without a structure. These edits were preserved; they were not treated as deployed on the peers. Current peer runtime fingerprints differed from the local checkout.
