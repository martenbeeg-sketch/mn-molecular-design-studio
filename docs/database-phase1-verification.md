# Molecular Database Phase 1 verification

**Checkpoint date:** 2026-11-04

**Specification:** [Approved Phase 1 contract](database-phase1-contract.md)

**Source instructions:** [Original request with user-recovered continuation and subsequent clarifications](database-phase1-source-instructions.md)

**Status:** Corrective regression selection passes. Phase 1 remains in progress and awaits another independent Sol review. This report is not production verification or Phase 1 approval.

## Scope and safety

This checkpoint covers the corrective implementation, isolated test databases and fixtures, and read-only R3 audit exports. It did not modify or rebuild production Database content, run a bulk import, contact or synchronize production peers, restart a service, deploy, or push changes. The three-machine test uses simulated peer nodes and temporary isolated stores.

## Test environment and result

Tests were run with the app’s configured Conda interpreter, `/home/user/mambaforge/envs/mn-molecular-design-studio/bin/python` (Python 3.11, Streamlit 1.65.0, Plotly 6.9.0, pytest 9.1.1). This is the interpreter configured by the Studio UI service. No package was installed or changed in that environment.

After the corrective changes, the expanded relevant regression selection passed **987 tests, 15 skipped, 0 failed**, in 284.01 seconds across 56 test files. The 15 skips require explicit QC-container or GPU-test opt-in: 10 QC-container checks and 5 GPU-inference checks. Pytest emitted 84 warnings (80 Biopython alignment deprecation warnings and 4 Pydssp runtime warnings); they did not affect the result. The selection includes the previously affected Complex Analysis and ligand-preparation tests.

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
- Every real analysis workflow resolves to an exact source or sequence-defined subject; unresolved workflows do not create Database identities. Sequence-only analysis does not invent a structure association.
- Attached coordinate files are checked against molecular composition, chains, sequence, and selected coordinate model. Coordinate content participates in immutable structure identity; Evaluations preserve both the source representation digest and the validated calculation-coordinate digest.
- Evaluation identity includes nested model, chain, ligand-residue selections, and scientific parameters. Missing historical settings do not inherit current defaults for compatibility or reuse.
- Local Evaluations on readable remote-only molecules preserve owner-qualified source structure references and create local Evaluation membership without modifying the remote owner's collection.
- Evaluation projection version v6 upgrades stale rows in bounded batches and preserves raw observations.
- Historical preview validates actual coordinate content and molecular correspondence, supports native BC1 model/rank naming with explicit compatibility limits, and rejects symlinked/unsafe/duplicate archive paths. No commit path is exposed.
- BindCraft2 generator/validation PAE and RFD4 bidirectional PAE are separate normalized measurements with stage and aggregation provenance.
- Both RFD4 sequence-carrier producers are classified as artificial sequence carriers; they remain valid sequence inputs and are not recorded as native structures.
- Input receptor structures attached during ingestion are made portable in the existing Database store.

## Independent review resolution matrix

| Finding | Original defect | Correction | Regression coverage | Result | Remaining limitation |
|---|---|---|---|---|---|
| H1 | Real `analysis.*` jobs could fall through to identity creation when Results/Database source resolution failed. | Analysis purpose is classified explicitly; structural analyses require exact source run/candidate/workstation and digest; sequence-only QC requires an exact explicit or source molecule and sequence. Unresolved analysis remains Results-only. | `test_real_analysis_workflows_without_exact_sources_create_no_library_records`; `test_results_and_pipeline_evaluations_attach_to_exact_source_structure`; `test_protein_qc_reuses_the_exact_frozen_database_identity`. | Passed in the expanded suite. | Not every registered adapter or historic output schema has a dedicated fixture; sequence-only matching is limited to supported sequence workflows. |
| H2 | Parseable secondary binder/target/complex/validation files could be attached despite wrong molecule, chain, sequence or model. | Validate each attached representation's coordinates and protein composition against declared chains/sequences and primary structure; require an actual selected model for multi-model coordinates; validate target input coordinates separately. | `test_secondary_candidate_structures_must_match_declared_chains_and_sequences`; `test_candidate_structure_model_selection_is_applied_before_attachment`; `test_mismatched_input_target_coordinates_do_not_attach_to_sequence_identity`. | Passed in focused and expanded tests. | Unsupported or ambiguous chain mappings and multi-model files without a numeric coordinate-model selector remain unresolved and are not attached. |
| H3 | Replacing coordinates at a reused path could leave the old portable structure attached to new metrics. | Structure IDs include the bytes' SHA-256; attachment checks the validated digest again; Evaluation context and identity distinguish source-structure and calculation-coordinate digests. | `test_changed_coordinates_at_a_portable_path_create_a_new_measurement_revision`; `test_database_ligand_analysis_ingestion_keeps_original_digest_after_evaluator_preparation`. | Passed in focused and expanded tests. | Calculation artifacts remain in Results unless an existing supported path stores them as Database representations; their validated digest is preserved in the Evaluation. |
| H4 | A remote-source Evaluation could replicate with a structure ID that did not resolve on a third workstation. | Preserve source owner workstation, group and structure IDs separately from local replica IDs through capture/apply and onward replication. | `test_three_workstations_ingest_offline_converge_retrieve_evaluate_and_tombstone` now covers owner A → Evaluation on B → C → return to A. | Passed in isolated three-node simulation. | This verifies the simulated transport and temporary stores, not production peer convergence or every network-failure mode. |
| H5 | Nested model, ligand residue and component selections were omitted from Evaluation identity. | Build a canonical nested selection digest and include it in Evaluation subject/instance projections. | `test_nested_model_and_ligand_residue_selections_are_part_of_evaluation_identity`. | Passed in the expanded suite. | Unrecognized selection schemas are retained as native metadata but are not normalized into a common selection vocabulary. |
| H6 | Missing historical settings were treated as today's pH/cutoff/protocol/image defaults during reuse. | Missing pH/cutoffs remain unknown; reuse also requires exact stored protocol revision and installed image digest. | `test_qc_reuse_requires_saved_parameters_protocol_and_image_identity`. | Passed in the expanded suite. | Reuse stays conservative when an old run lacks parameter, protocol or image provenance. |
| H7 | Local docking analysis of a remote-only compound could create another identity or attempt a write in its owner's collection. | Resolve the accessible canonical identity across readable memberships and create the Observation in the local Evaluation group. | `test_docking_evaluates_remote_only_compounds_in_a_local_evaluation_group`. | Passed in the expanded suite. | Requires the remote molecule to be readable and replicated locally; no identity is fabricated when exact chemistry does not resolve. |
| M1 | Evaluation projection version changes had no safe incremental upgrade path; old metric rows could be destructively cleared. | Versioned `evaluation.v6` projection rebuild processes bounded batches and replaces metrics per observation without deleting unrelated or raw observations. | `test_legacy_evaluation_projection_upgrade_is_bounded_and_preserves_raw_observations`. | Passed on an isolated legacy fixture. | Production projections were not rebuilt; deployment-time continuation remains outside this task. |
| M2 | Historical structure checks relied on filenames/checksums without verifying valid coordinates and molecular correspondence. | Preview parses actual structures and checks binder/target chain sequences; ranked and refolded structures remain distinct. | `test_historical_preview_requires_actual_matching_coordinates_for_structures`; BC2 stage-separation tests. | Passed in synthetic fixtures. | R3 original structures and raw source files are unavailable here, so no real-campaign coordinate or checksum claims are made. |
| M3 | Native BC1 `<Design>_modelN.pdb` and rank-prefix files were missed; unsupported metrics could imply full compatibility. | Recognize exact design/model suffix and ranked prefixes; preserve native metrics and report unsupported mappings/layouts explicitly. | `test_bindcraft1_reports_unrecognized_native_metric_definitions_and_layouts`; `test_bindcraft1_native_tables_are_detected_without_inventing_acceptance`. | Passed in synthetic fixtures. | No verified normalized mapping for averaged native iPTM and no user-visible commit/reconciliation workflow. |
| M4 | Symlinked source tables could alter consumed evidence without changing historical import identity. | Reject root/nested symlinks, unsafe archive paths and duplicate normalized ZIP entries; source identity covers the accepted file tree. | `test_historical_source_rejects_symlinked_files_and_root_symlinks`; `test_historical_zip_rejects_duplicate_normalized_paths`; unsafe ZIP path test. | Passed in synthetic fixtures. | Historical import stays preview-only; source trees outside the tested directory/ZIP formats are unsupported. |
| M5 | BindCraft2 validation PAE projected as generator PAE, and BC2 directional/RFD4 bidirectional values shared a metric identity. | Split generator/validation and directional/bidirectional normalized fields; apply stage restrictions and retain method/aggregation provenance. | `test_generator_confidence_has_declared_units_scope_and_stage`; `test_rfd4_oracle_aliases_preserve_native_scale_and_require_rfd4_generation`; `test_bindcraft_and_rfd4_interface_pae_share_31a_normalization_but_keep_provenance`. | Passed in the expanded suite. | No universal comparison across these scientifically different aggregations is provided. |
| M6 | `sequence_input.cif` placeholder coordinates were treated as native structural evidence. | Classify both RFD4 artificial sequence-carrier outputs as sequence-only conditioning artifacts across generic source capture. | `test_rfd4_idr_sequence_carrier_is_indexed_as_sequence_not_structure`; `test_rfd4_sequence_input_carrier_is_never_registered_as_native_reference`. | Passed in the expanded suite. | Other third-party artificial carriers require explicit producer/representation provenance; unknown layouts remain unclassified. |

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

The verified implementation is recorded in three local code/test commits:

- `86ec13b` — exact Evaluation source, coordinate revision and ownership provenance corrections.
- `92517ec` — historical campaign coordinate, naming and source-path validation.
- `b2d0908` — separate generator, validation and directional PAE definitions.

The contract and this verification report are updated in the documentation commit recorded in the handover. The earlier Phase 1 commits remain in history. No commit was pushed.
