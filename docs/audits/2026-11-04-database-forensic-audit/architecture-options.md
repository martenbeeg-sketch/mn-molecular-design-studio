# Remediation options: evidence and trade-offs

## Recommendation

The evidence supports **Option B after first applying the policy changes required by Option A**: keep the existing molecule / structure / observation / collection concepts, stop indiscriminate job indexing, then reconstruct a clean logical index from authoritative Results plus explicitly curated, imported, reference, and policy-accepted records. Preserve the current SQLite database and replication catalogs as read-only rollback evidence until all three workstations converge.

The current schema's main concepts fit the desired library. The clearest problems are default-automatic ingestion, broad input persistence, groups created as a side effect of indexing, incomplete acceptance classification, legacy metadata fallbacks, and peer-state divergence. No observed failure requires a brand-new schema. Option C would add a second representation and migration surface before the current logical model has been tested under a corrected ingestion boundary.

This is a recommendation from the audit, not permission to mutate production. No migration or cleanup was performed.

## Comparison

| Measure | A — Repair in place | B — Rebuild existing index | C — New schema and controlled migration |
|---|---|---|---|
| Data model | Keep tables; fix policy and UI. | Keep identity/structure/observation concepts; repopulate new index from selected authoritative sources. | Define new version and migrate selected information into it. |
| Complexity | Medium. Multiple ingestion entry points and peer semantics still need correction. | High. Requires full source inventory, selection manifest, crosswalk, and 3-node cutover. | Very high. Adds schema translation, compatibility/dual-read or dual-write, and cluster protocol work. |
| Retention estimate | 100% of current rows if no deletion; clutter remains unless views are restricted. | Initial reusable identity target is 3,848 distinct IDs (3,080 curated + 43 reference + 725 explicitly accepted generated, disjoint in the current query), plus selected structures and scientifically valuable observations. Additional screen hits can be explicitly promoted. | Same policy-defined retained content as B; schema novelty does not determine which science is worth keeping. |
| Risk of losing unique scientific records | Low if records remain untouched; risk of hiding useful data in views. | Medium until imports, selections, source-only proteins, acceptance evidence, observations, and remote files have checksummed migration coverage. | High until new schema encodes every retained identity, structure, policy, evaluation, and ownership edge and full round-trip is proven. |
| Cluster effect | Patch sync in place; old tombstone/history burden remains. | Must reconcile source-author data and install the same new index/epoch on pc002975, PC004032, and PC004033. | Same, plus all replication record kinds and snapshots need schema-version compatibility. |
| Results and scientific history | Unchanged. | Keep Results and artifacts intact; rebuild only the reusable index. | Keep Results and artifacts intact; new DB must continue linking exact jobs, inputs, and source structures. |
| Future maintainability | Improves policy but retains accumulated rows and implicit legacy classifications. | Good if there is one canonical ingestion contract and rebuildable projections; old schema still needs explicit versioning. | Potentially good only if it removes real model limitations; adds versioning and maintenance burden immediately. |
| Rebuild-only clutter risk | Not applicable. | High unless automatic policy, screening selection, empty-group creation, and unknown acceptance are corrected first. | High for the same reason; new schema would receive the same noisy inputs. |
| Rollback | Preserve DB snapshot and revert application release. | Keep old SQLite/files and replication catalogs untouched; switch back to old index/epoch if validation fails. | Keep both old and new stores with explicit version routing until full rollback window expires. |

## Required changes before Option B

1. Make run inclusion explicit by workflow. The current `automatic` fallback is incompatible with the requested library philosophy.
2. Keep binders in Results until a policy-specific explicit acceptance decision and required structure are present. Preserve unreviewed state as unknown.
3. Default prediction, analysis, affinity, docking, and screening to Evaluation-on-existing-input or Results-only; never create a fresh per-run reusable collection for an analysis-only workflow.
4. Keep protein/ligand input source capture separate from output promotion. For large screens, save selected hit snapshots only.
5. Make group creation transactional with the first surviving membership; prevent empty groups by construction.
6. Enforce exact identity reuse; structures remain attached to the canonical subject and carry source digest, model, chain mapping, and run provenance.
7. Keep raw observations; rebuild `molecule_evaluations` and `molecule_metric_values` as projections. Do not merge metrics from different engines into one universal score.
8. Define a cluster snapshot/epoch cutover. Preserve tombstones until a known-good new epoch is installed on all three peers.

## Safe staged sequence if B is later authorized

1. **Freeze the specification, not production:** define keep rules by workflow and acceptance policy. Do not infer missing status from metrics or names.
2. **Inventory source availability:** for each identity/structure/observation, record authoritative workstation, Results run, artifact path, SHA-256, and whether the owner is reachable. Use `unknown` for missing evidence.
3. **Create verified offline backups:** SQLite snapshots using SQLite's backup API, file-store manifests and checksums, raw peer catalogs, import-source files, and Results inventories. Keep at least one untouched copy per node.
4. **Build a non-production candidate index:** consume only imported/reference source facts, saved selections, accepted policy evidence, and selected Results. Keep a molecule/structure/evaluation/provenance crosswalk to the old database.
5. **Compare before cutover:** identity fingerprints, membership counts, target/binder chain links, PDB reference correspondence, exact accepted-policy digests, structure hashes, raw observations, metric projections, imported row audits, and selected compound sets. Test each peer's source-owned rows independently.
6. **Exercise rollback offline:** prove that old DB and replication catalog restore without touching source Results. Confirm a candidate failure leaves old content intact.
7. **Cut over one isolated node only after approval:** pin the new epoch, keep the prior DB and replica catalog mounted read-only, and verify API/UI/report results before any other node changes.
8. **Converge the three nodes from one owner manifest:** compare source-author keys and tombstones, then verify all structures by digest with owner authorization. Do not let peers independently rebuild different subsets.
9. **Retain both versions through a review window:** delete or compact nothing until the user accepts reconciliation evidence and the new epoch has proven stable through a normal sync cycle.

## User-visible concepts

- **My Molecules** should be the general identity view.
- **My Proteins** and **My Compounds** are type filters over the same identity index, not separate identity stores.
- **Structures** is a representation view; it distinguishes supplied/experimental coordinates, predictions, complexes, and selected poses.
- **Reference Structures** is a source/representation filter with accession, source file, chain correspondence, and sequence relationship.
- **Collections** are explicit reusable sets with provenance and ownership.
- **Protein Library** should be the reference/source collection view within the general molecular library, sharing canonical identities with generated designs and evaluations.
- **My Results** remains the full run/artifact history, including failures, every candidate, intermediates, logs, and unselected screen poses.

This puts the recent Protein Library work in the general library without creating a second identity system. It also keeps collection membership separate from identity and structure representation.
