# Molecular Database Phase 1 verification

**Checkpoint date:** 2026-11-05

**Specification:** [Approved Phase 1 contract](database-phase1-contract.md)

**Source instructions:** [Original request with user-recovered continuation and subsequent clarifications](database-phase1-source-instructions.md)

**Status:** C1–C9 and F1–F3 corrective regressions pass in isolated tests. The final F2 annotation-compatibility correction below also passes its focused regression and surrounding verification. Phase 1 remains in progress and awaits one narrow independent Sol review. This report is not production verification or Phase 1 approval.

## Scope and safety

This checkpoint covers the corrective implementation, isolated test databases and fixtures, and read-only R3 audit exports. It did not modify or rebuild production Database content, run a bulk import, contact or synchronize production peers, restart a service, deploy, or push changes. The three-machine test uses simulated peer nodes and temporary isolated stores.

## Test environment and result

Tests were run with the app’s configured Conda interpreter, `/home/user/mambaforge/envs/mn-molecular-design-studio/bin/python` (Python 3.11, Streamlit 1.65.0, Plotly 6.9.0, pytest 9.1.1). This is the interpreter configured by the Studio UI service. No package was installed or changed in that environment.

After the final C1–C9 corrections, the expanded relevant regression selection passed **1,000 tests, 15 skipped, 0 failed**, in 276.77 seconds across the same 56 test files. The 15 skips require explicit QC-container or GPU-test opt-in: 10 QC-container checks and 5 GPU-inference checks. Pytest emitted 84 warnings (80 Biopython alignment deprecation warnings and 4 Pydssp runtime warnings); they did not affect the result. The selection includes the previously affected Complex Analysis and ligand-preparation tests. The first expanded attempt exposed four reference-placement regressions caused by applying Results binder/target roles to a Database source component selection. After limiting the allowance to an exact Database source-to-calculation chain/sequence map, the four affected workflows passed and the full expanded selection was rerun successfully.

The focused eleven-module C1–C9 and preserved H/M regression selection passed **254 tests in 40.27 seconds** after the final code correction. Tests ran with `conda run -n mn-molecular-design-studio`; the interpreter is `/home/user/mambaforge/envs/mn-molecular-design-studio/bin/python` (Python 3.11). No production Database, Results, projection, replication catalog, service, or peer was changed or contacted.

After the specific F1–F3 corrections, the three-module focused selection (`test_molecular_database_policy.py`, `test_complex_analysis.py`, `test_bound_docking.py`) passed **121 tests in 13.64 seconds**. The six parameterized F1–F3 cases passed; the prepared-PLIP preservation and Results-source peer controls passed in the same eight-item focused run. The expanded relevant regression selection covering Database policy, peers, Evaluation projections/comparisons, historical import, PAE normalization, Complex Analysis, bound docking, and placement workflows passed **343 tests, 10 skipped, 0 failed in 78.81 seconds**. The skips are the opt-in QC-container tests; 50 Biopython alignment deprecation warnings were reported. The expanded selection includes the frozen-Evaluation, peer replication/delivery-order, projection-upgrade, historical-import and PAE regressions. The isolated actual-legacy Database reproduction also passed against a temporary SQLite backup, as detailed below.

The selection was:

```bash
/home/user/mambaforge/envs/mn-molecular-design-studio/bin/python -m pytest -q \
  tests/test_accepted_design_import.py \
  tests/test_api_extended.py \
  tests/test_automatic_batches.py \
  tests/test_boltzgen_candidate_normalization.py \
  tests/test_bulk_screening_collection.py \
  tests/test_candidate_evidence.py \
  tests/test_distributed_campaign.py \
  tests/test_legacy_migration.py \
  tests/test_mimosa_integrity_filter.py \
  tests/test_mimosa_workflow_wiring.py \
  tests/test_native_collection_execution.py \
  tests/test_prediction_measurements.py \
  tests/test_quality_tools.py \
  tests/test_rfd4_proteina.py \
  tests/test_scientific_gpu.py \
  tests/test_two_stage_mpnn.py \
  tests/test_validation_acceptance.py \
  tests/test_cofolding_msa.py \
  tests/test_docking_preparations.py \
  tests/test_md_database.py \
  tests/test_historical_results_import.py \
  tests/test_output_vocabulary.py \
  tests/test_evaluations_and_comparisons.py \
  tests/test_molecular_database.py \
  tests/test_molecular_database_background.py \
  tests/test_molecular_database_removal.py \
  tests/test_molecular_database_run_cleanup.py \
  tests/test_molecular_database_policy.py \
  tests/test_molecular_database_integration.py \
  tests/test_molecular_database_peer.py \
  tests/test_molecular_database_results_ui.py \
  tests/test_protein_library.py \
  tests/test_library_preparation.py \
  tests/test_api_integrations.py \
  tests/test_bound_docking.py \
  tests/test_docking_results.py \
  tests/test_cofolding_results.py \
  tests/test_complex_analysis*.py \
  tests/test_pipeline_adoption.py \
  tests/test_pipeline_canvas_audit.py \
  tests/test_pipeline_graph.py \
  tests/test_pipeline_job_hierarchy.py \
  tests/test_sequence_pipeline_handoff.py
```

## Distributed ingestion and replication correctness

`tests/test_molecular_database_peer.py::test_three_workstations_ingest_offline_converge_retrieve_evaluate_and_tombstone` passes using isolated simulated nodes named `pc002975`, `PC004032`, and `PC004033`. It verifies that each node ingests locally, one node remains offline while the other nodes replicate, reconnect/replay converges records without duplicate canonical identities, and the independently authored structures, source memberships, observations, and Evaluation records survive.

The same fixture verifies owner provenance on replicated structures, authorized retrieval and byte equality, denial to an unauthorized user, and continued access after the owning node goes offline once the content has been materialized. It then analyzes a molecule/structure authored by another node: the exact molecule ID, structure ID, digest, source group, owner, and source node remain in Evaluation context, while the new Evaluation membership is local to the analyzing node. A source-author deletion creates a tombstone scoped to that author’s group and does not remove other authors’ memberships or Evaluation records.

This establishes correctness for the simulated transport and isolated database implementation. It does not establish that production peers currently converge or reconcile any production discrepancy.

## Evaluation attachment coverage

The tested paths include:

- Results and pipeline outputs attaching to the exact Results molecule and structure, with original digest and source kind retained (`test_results_and_pipeline_evaluations_attach_to_exact_source_structure`).
- Complex Analysis outputs attaching to the exact selected Database representation, including ligand analysis whose evaluator-prepared structure has a different digest (`test_database_analysis_ingestion_attaches_evaluation_to_exact_existing_representation`, `test_database_ligand_analysis_ingestion_keeps_original_digest_after_evaluator_preparation`).
- A registered standalone analysis workflow retaining exact source structure and Evaluation context (`test_registered_workflows_and_database_adapter_preserve_existing_structure_and_ingest_evaluation`).
- A local Evaluation attached to a remotely authored molecule and structure in the three-peer test above.
- Wrong source digests, inaccessible representations, or mismatched models being rejected in Complex Analysis source and dashboard tests.

This covers the shared ingestion/adaptation paths represented by these fixtures. It is not an exhaustive test of every registered standalone tool, every pipeline adapter, or every historic output schema.

## Structure retention and reconstruction

Passing isolated tests establish these specific paths:

- Curated protein designs retain distinct binder and target identities and keep binder, target, and target–binder complex representations with portable content and checksums (`test_curated_design_keeps_binder_complex_and_supplied_target_as_distinct_representations`). Sequence-only target records remain valid without fabricated coordinates.
- Prediction/cofolding indexing retains only the selected best structure for an input and preserves the other evaluation observations (`test_prediction_indexes_only_the_best_structure_per_input`). This fixture does not establish every engine’s target-structure association behavior.
- Docking associates the canonical ligand identity with one selected representative pose, keeps an exact receptor structure separately, and does not create compound identities per pose (`test_docking_indexes_only_the_top_pose_for_each_compound`, `test_compound_and_target_remain_separate_and_share_portable_receptor`).
- Mapped PDBQT reconstruction produces coordinate-bearing SDF in the saved ligand atom order and validates receptor/pose coordinate-frame identity. Ambiguous maps and ligand-identity mismatches fail (`test_md_staging_keeps_mapped_coordinates_atom_order_and_database_identities`, `test_md_rejects_ambiguous_mapping_instead_of_inventing_coordinates`, `test_md_rejects_pose_chemistry_that_does_not_match_saved_ligand_identity`).
- Exported bound complexes retain exact receptor, pose, digest, ligand identity, and atom mapping provenance; registration is idempotent and the complex can reopen in Database analysis (`test_bound_complex_is_attached_to_same_record_with_chemical_identity_and_provenance`, `test_bound_complex_registration_is_idempotent_and_respects_group_write_access`).
- Authorized remote structure retrieval is checked by the three-peer fixture and the existing portable structure peer test.

The tests do not prove that every workflow retains every supplied and prepared ligand SDF/PDBQT, or that all source/experimental/predicted/reconstructed labels and historical chain/atom mappings are complete. Those remain open as described in the contract.

## Historical BindCraft and R3 reconciliation

Synthetic importer tests pass for BindCraft1/BindCraft2 layouts, distinct trajectory/refold/ranked stages, exact hash-and-sequence joins, preservation of unknown acceptance, explicit selection, source-content identity across relocation/renamed ZIP, unsafe path rejection, and read-only/ACL-limited identity lookup.

For the R3 campaigns, the original run directories reported by the audit are not present on this host. The available read-only export contains native table rows and recorded original file checksums, but not the original source CSV bytes or structural files. A read-only preview was run in isolated temporary folders using those exported native rows. It reconciled:

| BindCraft2 run | Trajectories | Refold attempts | Native passing refolds | Ranked designs | Exact ranked hash+sequence joins |
|---|---:|---:|---:|---:|---:|
| `20261019-181522-21f9a96c` | 5,000 | 9,475 | 143 | 65 | 65 |
| `20261021-184403-b0c849c0` | 3,665 | 6,614 | 89 | 35 | 35 |
| **Production total** | **8,665** | **16,089** | **232** | **100** | **100** |

All 100 ranked rows joined to passing refold rows in the captured native tables. The importer preview reported 9,475 and 6,614 distinct hash/sequence-associated candidates for the two campaigns. The separate one-design smoke run is not included in this production total.

This is reconciliation of the audit’s native-row snapshot, **not direct verification against the original campaign files**. Original CSV checksums could not be recomputed from missing source bytes, and no R3 structure files were available to match or checksum. No historical Database import or commit was performed. Thus direct real-file R3 importer/structure validation remains unsupported pending access to those source files.

## Ordinary Database reads

UI/API fixtures pass showing ordinary Database reads, background refresh paths, Results listing, and peer metadata reads do not trigger ingestion. Explicit refresh remains the ingestion trigger. Tests include `test_database_ui_read_does_not_start_background_ingestion`, `test_database_reads_do_not_start_ingestion_but_explicit_refresh_does`, `test_database_and_results_reads_do_not_trigger_ingestion`, and `test_results_listing_is_a_read_and_leaves_candidates_unindexed`.

## Already implemented, corrected, and open

### Already implemented before this Phase 1 work

- Shared canonical molecule identities with separate structures, source records, memberships, observations, and Evaluations.
- Existing owner-authoritative replication, ACL checks, content-addressed structure storage, remote retrieval, and deletion tombstones.
- Complex Analysis exact-source authorization and bound-docking reconstruction/export with atom-map validation.
- Results retained separately from curated Database records.

### Corrected by this implementation

- Policy-less or unclassified historical runs no longer imply automatic promotion. Rejected, unknown, intermediate, and bulk screen inputs remain Results-only unless explicitly selected/accepted under the established policy.
- Source capture and output promotion are handled separately; reindexing is additive and does not reinterpret earlier saved acceptance decisions.
- Database UI/API reads and background refresh do not initiate ingestion; explicit refresh/import paths remain available.
- RosettaLigand's registered docking workflow is classified as docking while registered Rosetta analysis remains exact-source analysis. Structural analyses require source-to-calculation molecular correspondence; unresolved or mismatched calculations do not create identities or Evaluations. Sequence-only analysis does not invent structure references.
- Attached and historical structural evidence is checked for finite coordinates, declared composition, chains, sequence and selected model. Binder-only representations and coordinate-free sequence targets remain valid in their declared contexts. Changed coordinates or calculations create immutable structure/Evaluation revisions; exact replay is idempotent and frozen comparisons retain their original evidence.
- Identical coordinate evidence deterministically reuses legacy structure IDs; different coordinate/model/candidate evidence remains distinct. Evaluation identity retains nested model, chain, ligand-residue selections and exact historical scientific parameters. Missing historical settings do not inherit current defaults.
- Database-source and Results-source Evaluations preserve owner-qualified structure references. Pending associations reconcile as owner records arrive, and local Evaluations on readable remote-only molecules do not modify the remote owner's collection.
- Evaluation projection version v6 upgrades stale rows in bounded batches to query readiness while preserving raw observations. Public queries beyond 500 legacy observations return current filters, selections, sorting and model/instance identity.
- Historical preview validates actual molecular content, supports native BC1 model/rank naming with explicit compatibility limits, and rejects symlinked/unsafe/duplicate archive paths. Both supported RFD4 sequence-carrier contents remain nonstructural after rename. No historical commit path is exposed.
- BindCraft2 generator/validation PAE and RFD4 bidirectional PAE are separate normalized measurements with stage and aggregation provenance.
- Both RFD4 sequence-carrier producers are classified as artificial sequence carriers; they remain valid sequence inputs and are not recorded as native structures.
- Input receptor structures attached during ingestion are made portable in the existing Database store.

## Independent review resolution matrix

| Finding | Original defect | Correction | Regression coverage | Result | Remaining limitation |
|---|---|---|---|---|---|
| H1 | Real analysis jobs could fall through to identity creation when source resolution failed; a broad `rosetta` marker also suppressed RosettaLigand docking. | Analysis-only classification is purpose-specific; structural analysis requires exact sources, and registered RosettaLigand stays on the docking path. | `test_real_analysis_workflows_without_exact_sources_create_no_library_records`; `test_results_and_pipeline_evaluations_attach_to_exact_source_structure`; `test_registered_rosettaligand_docking_keeps_selected_pose_and_evaluation`. | Passed in final focused and expanded selections. | Not every registered adapter or historic output schema has a dedicated fixture; sequence-only matching is limited to supported sequence workflows. |
| H2 | Secondary structures or calculations could be associated despite wrong molecule, chain, sequence or model. | Validate calculation composition against the exact selected representation, including its full protein chain/sequence map and selected ligand evidence; allow reciprocal Database complex roles only when full source/calculation chain→sequence mappings are identical. | `test_results_analysis_rejects_calculation_sequence_mismatch_with_exact_source`; `test_database_analysis_ingestion_attaches_evaluation_to_exact_existing_representation`; `test_binder_owned_complex_analysis_matches_selected_representation_composition`; `test_ligand_owned_bound_complex_plip_checks_selected_representation`; `test_secondary_candidate_structures_must_match_declared_chains_and_sequences`. | Passed after the F1/F2 correction in the final focused and expanded selections. | Unsupported chain/residue changes and ligand preparation that alters heavy-atom geometry without an existing mapping remain unresolved. |
| H3 | Coordinate or score replacement at a reused path/key could change old Evaluation evidence and frozen comparisons. | Changed calculation inputs create distinct immutable observation/Evaluation revisions; exact evidence replay is idempotent; old structure IDs and digests are retained. | `test_changed_coordinates_at_a_portable_path_create_a_new_measurement_revision`; `test_changed_calculation_evidence_preserves_frozen_evaluation_and_is_idempotent`; `test_identical_job_refresh_is_idempotent_and_changed_scores_create_immutable_revision`. | Passed; frozen score 0.8 retains its original Evaluation and structure after score 0.4 arrives. | Evaluator-prepared calculation files remain in Results unless an existing supported retention path attaches them; their digest stays in Evaluation evidence. |
| H4 | A remote-source Evaluation could replicate with dangling structure IDs or fail to resolve when source records arrived later. | Preserve owner workstation/group/structure separately from local replica IDs for Database and Results sources; reconcile pending references when owner records arrive. | `test_results_source_evaluation_preserves_owner_reference_across_three_peers_and_delivery_orders`; `test_database_source_evaluation_reconciles_when_owner_records_arrive_late`; three-workstation offline/replay/tombstone test. | Passed both delivery orders, duplicate replay and A→B→C→A in isolated peers. | Simulated transport is verified; production peer convergence and permanent owner-file loss are not exercised. |
| H5 | Nested model, ligand residue and component selections were omitted from Evaluation identity. | Build a canonical nested selection digest and include it in Evaluation subject/instance projections. | `test_nested_model_and_ligand_residue_selections_are_part_of_evaluation_identity`. | Passed in the expanded suite. | Unrecognized selection schemas are retained as native metadata but are not normalized into a common selection vocabulary. |
| H6 | Missing historical settings were treated as today's pH/cutoff/protocol/image defaults during reuse. | Missing pH/cutoffs remain unknown; reuse also requires exact stored protocol revision and installed image digest. | `test_qc_reuse_requires_saved_parameters_protocol_and_image_identity`. | Passed in the expanded suite. | Reuse stays conservative when an old run lacks parameter, protocol or image provenance. |
| H7 | Local docking analysis of a remote-only compound could create another identity or attempt a write in its owner's collection. | Resolve the accessible canonical identity across readable memberships and create the Observation in the local Evaluation group. | `test_docking_evaluates_remote_only_compounds_in_a_local_evaluation_group`. | Passed in the expanded suite. | Requires the remote molecule to be readable and replicated locally; no identity is fabricated when exact chemistry does not resolve. |
| M1 | A bounded single-batch v6 upgrade left public queries stale or incomplete beyond the first 500 observations. | Public query readiness drains bounded projection batches before returning data; filters, selections, sorting and nested model/instance identity use current projections while raw observations remain unchanged. | `test_public_queries_do_not_return_stale_or_incomplete_rows_during_bounded_upgrade`; `test_legacy_evaluation_projection_upgrade_is_bounded_and_preserves_raw_observations`. | Passed with 501 stale legacy observations in an isolated fixture. | No production projection rebuild or migration/readiness operation was run. |
| M2 | Historical structure checks did not require finite coordinates or declared target–binder composition. | Parse actual structure content, reject NaN/Inf, and enforce target-chain correspondence for complex evidence while allowing declared binder-only structures and sequence-only targets. | `test_historical_complex_evidence_requires_finite_coordinates_and_full_composition`; `test_declared_binder_only_structure_and_sequence_only_target_remain_valid`; `test_historical_preview_requires_actual_matching_coordinates_for_structures`. | Passed in synthetic fixtures; invalid cases produce no proposed Evaluations. | R3 original structures and raw source files are unavailable here, so no real-campaign coordinate or checksum claims are made. |
| M3 | Native BC1 `<Design>_modelN.pdb` and rank-prefix files were missed; unsupported metrics could imply full compatibility. | Recognize exact design/model suffix and ranked prefixes; preserve native metrics and report unsupported mappings/layouts explicitly. | `test_bindcraft1_reports_unrecognized_native_metric_definitions_and_layouts`; `test_bindcraft1_native_tables_are_detected_without_inventing_acceptance`. | Passed in synthetic fixtures. | No verified normalized mapping for averaged native iPTM and no user-visible commit/reconciliation workflow. |
| M4 | Symlinked source tables could alter consumed evidence without changing historical import identity. | Reject root/nested symlinks, unsafe archive paths and duplicate normalized ZIP entries; source identity covers the accepted file tree. | `test_historical_source_rejects_symlinked_files_and_root_symlinks`; `test_historical_zip_rejects_duplicate_normalized_paths`; unsafe ZIP path test. | Passed in synthetic fixtures. | Historical import stays preview-only; source trees outside the tested directory/ZIP formats are unsupported. |
| M5 | BindCraft2 validation PAE projected as generator PAE, and BC2 directional/RFD4 bidirectional values shared a metric identity. | Split generator/validation and directional/bidirectional normalized fields; apply stage restrictions and retain method/aggregation provenance. | `test_generator_confidence_has_declared_units_scope_and_stage`; `test_rfd4_oracle_aliases_preserve_native_scale_and_require_rfd4_generation`; `test_bindcraft_and_rfd4_interface_pae_share_31a_normalization_but_keep_provenance`. | Passed in the expanded suite. | No universal comparison across these scientifically different aggregations is provided. |
| M6 | RFD4 sequence-carrier classification depended on producer filename, so rename to `target.cif` promoted artificial coordinates as a structural reference. | Recognize supported artificial-carrier content markers regardless of basename; sequence identity remains valid without coordinates. | `test_rfd4_idr_sequence_carrier_is_indexed_as_sequence_not_structure`; `test_rfd4_sequence_input_carrier_is_never_registered_as_native_reference`; `test_rfd4_artificial_carrier_remains_nonstructural_after_copy_to_target_cif`. | Passed for both supported producers and renamed files. | Unknown third-party artificial carriers without supported content/provenance markers remain unclassified. |

## Final corrective resolution matrix: C1–C9

Permanent regressions for C1–C9 were exercised against the pre-correction behavior before implementation. The failing pre-correction behaviors are recorded below; C6 dependency replay and C8 finite-coordinate/composition checks were also reconfirmed by disabling only the relevant correction in isolated fixtures after validating the fixtures. Final evidence is from temporary Database instances and simulated peers. Production data and replication state were not used.

| Finding | Original failure and root cause | Code correction | Regression and pre-correction → final result | Remaining limitation |
|---|---|---|---|---|
| C1 | `ligand.docking.rosetta_ligand` yielded no docking Evaluation because a generic `rosetta` marker classified it as analysis. | Removed the generic token; `analysis.rosetta` stays classified by its registered analysis purpose/name. | `test_registered_rosettaligand_docking_keeps_selected_pose_and_evaluation`: zero observations before → selected pose and Evaluation after. Positive Rosetta analysis exact-source controls pass. | Unregistered producer aliases require explicit workflow purpose or specific classification. |
| C2 | Valid KLM calculation coordinates could be attached to selected ACD source, and a valid ligand-owned bound complex was rejected because canonical identity components were treated as the complete selected system. | Validate the exact selected representation's parsed protein/ligand composition, sequences, chain/residue mappings and supported preparation evidence against the calculation. Preserve the canonical owner identity and fail closed when correspondence is unknown. | `test_binder_owned_complex_analysis_matches_selected_representation_composition`; `test_ligand_owned_bound_complex_plip_checks_selected_representation`; existing full-complex reciprocal-role, prepared-PLIP and receptor-substitution controls. Before → after: EFG source/KLM calculation **1 → 0**; unchanged receptor–ethanol PLIP **0 → 1**; receptor substitution remains **0**. | Exact chain maps are required unless a complete explicit map is present; selected ligand residue and heavy-atom geometry must remain identifiable through preparation. |
| C3 | Changed score/coordinates reused an observation key, replacing its projected Evaluation and changing the frozen comparison. | Store deterministic immutable observation revisions keyed by scientific evidence; reuse exact evidence on replay; keep coordinate digests in structure/Evaluation identity. | `test_changed_calculation_evidence_preserves_frozen_evaluation_and_is_idempotent`; `test_changed_coordinates_at_a_portable_path_create_a_new_measurement_revision`; background revision test: score 0.8/old structure persists, score 0.4 is distinct, exact replay adds no revision. | Unknown historic settings remain unknown; no missing provenance is fabricated. |
| C4 | Reingesting unchanged coordinates generated a competing structure ID when the actual legacy row had no `model_id`, making exact Results source resolution ambiguous. | Reuse missing model metadata only when exact digest/provenance match and parsing proves a single coordinate model whose actual model ID matches every supplied modern selector. Missing metadata remains unknown for multi-model or unresolved evidence. | `test_legacy_structure_reference_is_reused_for_identical_reingestion[None]`: missing legacy field creates duplicate before → original ID reused after; the Results source resolves to that same ID before and after. A different selector remains distinct. Actual legacy backup reproduction: **2→2** rows; source resolution **true→true**. | Missing/ambiguous multi-model legacy provenance remains distinct; no historical field is fabricated. |
| C5 | Results-source replication carried a local structure ID but not A's original owner-qualified reference. | Persist source owner node/group/structure for Results-source matches separately from B's Evaluation authorship and local replica IDs. | `test_results_source_evaluation_preserves_owner_reference_across_three_peers_and_delivery_orders`: missing owner link/dangling C reference before → C resolves A's exact structure while author remains B. | Transport simulation does not establish current production peer state. |
| C6 | Evaluation-before-source replay stayed unresolved after the source arrived because no dependency reconciliation was triggered. | Reconcile owner-qualified pending observations after owner group/entry application; resolution is exact and idempotent. | `test_database_source_evaluation_reconciles_when_owner_records_arrive_late`, plus both delivery orders, duplicates and return path in the C5 Results test: unresolved after arrival when reconciliation disabled → resolved automatically after fix. | Reconciliation follows supported peer record application; manual SQL edits are not watched. |
| C7 | A single bounded v6 batch left >500-observation public queries stale or missing metrics without readiness indication. | Drain bounded versioned batches before public projection-backed queries; raw observation rows are retained. | `test_public_queries_do_not_return_stale_or_incomplete_rows_during_bounded_upgrade`: 501 legacy observations with filters, Evaluation selection, metric sorting and nested model identity; incomplete before → complete current public rows after. | Upgrade behavior was tested only on isolated legacy fixtures; production projections were not rebuilt. |
| C8 | Historical preview marked NaN/Inf coordinates and binder-only coordinates as a target–binder complex available because names/checksums did not validate composition. | Shared parser rejects nonfinite selected-model coordinates; historical complex validation requires matching target and binder; declared binder-only and sequence-only-target cases remain supported. | `test_historical_complex_evidence_requires_finite_coordinates_and_full_composition`: NaN, Inf and incomplete complex were falsely available before → rejected with no Evaluation proposals after. `test_declared_binder_only_structure_and_sequence_only_target_remain_valid` passes. | This validator asserts protein chain correspondence; unsupported modified/nonprotein component semantics are not claimed. |
| C9 | Artificial CIF classification changed after renaming because it relied on producer filenames. | Read supported carrier content markers independent of basename; keep basename fallback only when source content is unavailable. | `test_rfd4_artificial_carrier_remains_nonstructural_after_copy_to_target_cif`: both producer outputs were treated as structures after rename before → remain sequence carriers after. Ordinary structure and sequence-only controls pass. | Unknown third-party artificial carrier content without registered markers/provenance remains unclassified. |

## Final F1–F3 corrective verification matrix

These three cases were independently reported as the remaining Phase 1 blockers. The permanent regressions were run before and after correction in isolated Database fixtures. The regression fixtures use the real registered workflow names and preserve the source identity while varying only the calculation or legacy metadata under test.

| Finding | Original failure and root cause | Correction | Regression test and before → after evidence | Remaining limitation |
|---|---|---|---|---|
| **F1 — binder-owned complex target substitution** | The source identity contained binder ACD only, while its authorized representation contained binder ACD plus target EFG. The matcher read only canonical identity components, so EFG→KLM was not checked and the KLM calculation was ingested. | Resolve and parse the exact authorized representation/model; apply its recorded binder/target chain roles; confirm the identity components are present in that representation; compare the complete source and calculation protein chain→sequence maps before Evaluation creation. | `test_binder_owned_complex_analysis_matches_selected_representation_composition` uses `analysis.biotite_interface`. Before: unchanged EFG calculation ingested (1), and substituted KLM calculation also ingested (1). After: EFG remains accepted (1), KLM is rejected (0); the canonical ACD identity and molecule count remain unchanged. | A changed chain label requires a complete explicit source-to-calculation map. |
| **F2 — ligand-owned bound-complex PLIP annotation compatibility** | After the representation-aware F2 correction, `_record()` replaced the selected coordinate residue's CCD with optional SMILES. The exact-source matcher still required the selected CCD, so an unchanged valid calculation with `ligand_smiles=["CCO"]` was rejected (0), while the same coordinates without that annotation were accepted (1). | Keep coordinate-derived CCD/residue identity alongside the optional SMILES. When a calculated ligand SMILES is present, canonicalize it and require chemical equivalence to the selected source ligand before continuing to the existing element-inventory and heavy-atom-geometry checks. Database canonicalization still reduces the owner identity to ethanol `CCO`. | `test_ligand_owned_bound_complex_plip_checks_selected_representation` uses the real `analysis.plip` Database source adapter and exact selected ligand `Z:1/LIG/CCO`. Before: the correct-annotation case failed (`ingest_run` returned 0, expected 1) on the pre-correction implementation. After: unannotated and correct `CCO` cases both ingest (1); conflicting `COC`, receptor substitution, ligand element substitution, altered heavy-atom geometry and wrong selected residue all reject (0). Positive cases assert identical source/calculation coordinate text; the correct case is replayed and remains one Evaluation. Exactly one canonical ethanol identity remains. Seven parameterized PLIP cases pass. | Unsupported ligand renumbering or geometry changes without an authorized mapping remain unresolved. |
| **F3 — legacy single-model representation reuse** | The legacy row omitted `model_id`; current ingestion supplied `"0"`. Direct `None != "0"` comparison created a second structure row even for unchanged single-model coordinates. | Compare exact digest and all recorded scientific fields. When model metadata is absent, parse the coordinate bytes with required single-model selection and reuse only if the actual model ID agrees with every supplied selector. | `test_legacy_structure_reference_is_reused_for_identical_reingestion[None]` constructs the legacy row with no `model_id` from the start. Before: a second same-digest ID was created and exact Results lookup became ambiguous. After: the old ID is reused, one structure row remains, and exact Results lookup returns that ID; model selector `"1"` remains distinct. Sol's actual-legacy reproduction was rerun from read-only source `/tmp/phase1-independent-probes-50vbhvcg/mismatched-binder/database/molecular_library.sqlite3` using a temporary SQLite backup at `/tmp/phase1-final-boundaries-y2draz47/actual-legacy-copy/database`: rows **2→2**, source resolution **true→true**, original ID retained. | A multi-model file, unknown actual model selector or non-coordinate representation without enough evidence is not treated as equivalent. |

The F1–F3 acceptance and preservation run used the configured `mn-molecular-design-studio` Conda interpreter. The expanded command covered 16 relevant test modules and completed with **343 passed, 10 opt-in QC-container skips, 0 failed**. It included the frozen Evaluation score/structure comparison, three-peer owner-qualified references and arrival-order reconciliation, 501-row projection readiness, historical source safety/composition, and distinct generator/validation/directional PAE normalization. Peer retrieval in the Results-source test was mocked to return the exact owner-published bytes; no production peer was contacted.

The final expanded regression command was:

```bash
/home/user/mambaforge/envs/mn-molecular-design-studio/bin/python -m pytest -q -p no:cacheprovider \
  tests/test_molecular_database_policy.py tests/test_molecular_database_peer.py \
  tests/test_evaluations_and_comparisons.py tests/test_historical_results_import.py \
  tests/test_output_vocabulary.py tests/test_complex_analysis_dashboard.py \
  tests/test_molecular_database_background.py tests/test_protein_library.py \
  tests/test_molecular_database_integration.py tests/test_complex_analysis.py \
  tests/test_validation_acceptance.py tests/test_complex_analysis_tag_fusion.py \
  tests/test_complex_analysis_placement_workflow.py tests/test_complex_analysis_liabilities.py \
  tests/test_bound_docking.py tests/test_quality_tools.py
```

## Final F2 annotation-compatibility correction

Sol's 2026-11-05 closure report reproduced one remaining F2 failure. I extended the permanent registered `analysis.plip` regression first and ran it against the unchanged implementation. The exact `correct-smiles` case failed as reported: the unchanged selected complex with `raw_metadata.ligand_smiles=["CCO"]` returned **0 ingestions** where **1** was expected; the logged rejection was that the analyzed chains did not correspond to the exact selected source molecule. The no-annotation case accepted the same coordinate construction.

The fix in [molecular_database_ingest.py](/home/user/programs/git-projects/mn-molecular-design-studio/mn_molecular_design_studio/molecular_database_ingest.py) retains the parsed CCD while adding SMILES. The correspondence matcher checks any calculated ligand SMILES against the canonical SMILES selected from the authorized source, then continues the existing CCD/residue, heavy-atom inventory and geometry checks. Canonical molecular identity handling remains unchanged and drops the redundant CCD when normalizing the identity. No Evaluation is admitted when the supplied chemistry conflicts.

The permanent regression uses the actual Database source adapter for `analysis.plip` and source ligand `chain=Z, residue_key=Z:1, CCD=LIG, SMILES=CCO`. The two positive cases build identical source and calculation coordinate text; only the optional raw SMILES annotation differs. Correct `CCO` annotation is accepted; conflicting `COC`, receptor sequence substitution, ligand element substitution, heavy-atom geometry change, and selecting `Z:2` reject. The test also verifies that ethanol remains one canonical molecule and that replaying accepted evidence leaves one PLIP Evaluation.

Pre-correction focused evidence: `test_ligand_owned_bound_complex_plip_checks_selected_representation[correct-smiles-1]` failed with `0 == 1`. Post-correction focused evidence: all seven parameterized cases passed (**7 passed**). The surrounding focused selection below passed **134 tests, 0 failed** in 39.02 seconds. It includes all bound-docking/PLIP tests, Complex Analysis correspondence tests, Database policy including `test_changed_calculation_evidence_preserves_frozen_evaluation_and_is_idempotent`, and Database integration.

```bash
/home/user/mambaforge/envs/mn-molecular-design-studio/bin/python -m pytest -q -p no:cacheprovider \
  tests/test_bound_docking.py tests/test_complex_analysis.py \
  tests/test_molecular_database_policy.py tests/test_molecular_database_integration.py
```

Code and permanent regression are committed as `3ce5a0c5f39f6629734c6f45de0b0818f813d8d2`. The Phase 1 verification update is a separate documentation commit. The three independent Sol review reports were left unchanged; their SHA-256 values after this correction remain `c4ef36aa4c8519ecefff27f9dba4e83832ab60a7d25b5de2d0ac8c9930083d69` (original), `17e2773e16d68d4565e0e2455118b401213f47c134abfbef9583b49d3bb917cd` (corrective), and `8a9bd0a699e7d71093bc857e8b2a450484975044b7d8bbfef3c0e9a5afe4c117` (final independent verification). The closure report was also not modified.

### Unsupported or still open

- Direct checksum and structural validation of the R3 source files is blocked because the original campaign directories and structures are missing from this host; only the read-only audit-row export was available.
- Historical BindCraft import/commit and cross-peer replay of committed imports are not implemented or authorized in Phase 1.
- Full preservation of every original/prepared SDF/PDBQT and complete structural mapping across all workflow producers remains unverified. Unsupported or ambiguous coordinate-model and chain mappings are rejected rather than inferred.
- Evaluation attachment was verified for representative standalone, Database, Complex Analysis, Results, pipeline, and peer-local paths, not every registered adapter or legacy schema.
- Historical parameter/protocol/image provenance remains unavailable for older QC runs; they conservatively do not qualify for exact reuse.
- Artificial sequence carriers from the two known RFD4 producers are classified; other third-party carrier layouts need explicit provenance and remain unknown.
- Existing forensic differences remain unresolved: 3 memberships, 3 structures, and 128 observations/source links across peer audit snapshots. The audit’s 631,100 tombstones remain untouched; 2,230 representation rows remain unclassified; possible content-store orphans were not deleted.
- Production peer convergence, migration/rebuild, and any cutover require separate authorization and are outside this Phase 1 checkpoint.

## Local corrective commits

The earlier verified Phase 1 implementation remains in these local code/test commits:

- `86ec13b` — exact Evaluation source, coordinate revision and ownership provenance corrections.
- `92517ec` — historical campaign coordinate, naming and source-path validation.
- `b2d0908` — separate generator, validation and directional PAE definitions.

The final corrective code and regression tests are committed in `30032c29d5cf1e3ba638d4811a9a383958b0b3a7`. The contract and this verification report are updated in the following documentation commit recorded in the handover. The earlier Phase 1 commits remain in history. No commit was pushed.

The F1/F2 representation-aware correspondence correction and regression cases are committed in `0fabd6e9a482ab764c8995d27d5fb153fc1b8df4`. The F3 single-model legacy reuse correction and missing-`model_id` regression are committed in `fdc315ef047916196ad86618161928621829eec5`. This F1–F3 verification update is committed separately; no commit was pushed.
