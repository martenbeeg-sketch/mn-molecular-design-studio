# P0 Database Simplification Implementation

Implemented the bounded changes from the 2026-11-05 simplification plan using the existing Database, Results, collection-policy, and replication architecture. No Phase 2B reconstruction was started.

## Before and after

| Behavior | Before | After |
|---|---|---|
| New active-learning final shortlist | Graph schema default, controller fallback, and API capability reported a 1,000-compound import limit. | Each defaults to 0. A saved historical graph with an explicit 1,000 value still migrates with that value unchanged. |
| Ordinary bulk docking | Docking import checkbox and direct/queued docking function defaults enabled Database ingestion. | Checkbox, direct-call defaults, and omitted queued-call defaults are Results-only. Explicit `import_to_database=True` remains supported. |
| Historical SpiroSpace shortlists | An isolated fixture showed both audited shortlist IDs in curated Library and All indexed. | The two exact owner/run collections are excluded from Library and remain in All indexed and Results. Independent curated membership keeps a shared molecular identity in Library. |
| Development visibility | Name or test/development-account flags could hide unrelated groups. | Only the eight audited owner/run pairs are hidden by the existing toggle. Name-only, ambiguous, R3, Tau, and PD-L1 fixture groups remain visible. |
| Explicit screening selection | Final import offered a default top-1,000 action. | Exact IDs or a positive reviewed top-N import require a user action. Exact imports use a separate collection and preserve the earlier automatic Evaluation record. |

The historical collection IDs are `screening-final-a5ba730c962dc449fb40bbab7702a25b` and `screening-final-7483c3d8846ffad1263b3b8d78078d3d`. The implementation checks the audited owner and run ID before excluding either collection from the curated view. It does not globally blacklist molecule IDs.

## Isolated evidence

- Baseline reproduction used an isolated fixture: import limit was `1000`; both shortlist IDs appeared in Library and All indexed.
- New graph/controller/capability defaults are asserted as `0`; the historical explicit `1000` graph value is asserted to survive graph migration.
- The shortlist fixture retains both memberships, the original observation, structure representation, and direct ID/history access. A molecule shared with an intentional import remains visible in Library; the shortlist-only identity does not.
- Exact selected import adds only the requested IDs, reuses an existing canonical ID, replays idempotently, leaves the previous automatic observation unchanged, and does not modify the Parquet results file.
- Results UI fixture starts with no selected IDs, top-N `0`, and disabled import buttons. It ingests only after an explicit selection and button click. Historical shortlist links open All indexed.
- Peer regression coverage exercises existing ACL, owner-authoritative behavior, convergence, and tombstone handling on isolated test state.

## Verification

All tests were run with `/home/user/mambaforge/envs/mn-molecular-design-studio/bin/python`, the app's Conda environment. The final combined Database, Results, screening, import, UI, and replication regression run passed **337 tests in 65.96 seconds**. It covered `test_database_simplification_p0`, `test_bulk_screening_collection`, `test_molecular_database_results_ui`, `test_active_learning_graph`, `test_active_learning_stopping`, `test_api_integrations`, `test_ligand_tools`, `test_molecule_imports`, `test_historical_results_import`, `test_accepted_design_import`, `test_docking_results`, `test_molecular_database_peer`, `test_molecular_database_policy`, `test_molecular_database`, `test_molecular_database_integration`, `test_screening_convergence`, `test_peer_workstations`, and `test_candidate_peer_transfer`. `git diff --check` was clean.

## Scope and remaining P1 work

The 53 identity-integrity cases and 128 local-calculation capture issue remain unchanged. Separate P1 work remains for the exact obsolete docking-reference review, missing original R3 generator/coordinate evidence, broader native Results inspection, and the local-analysis peer-capture defect. The 217 ambiguous smoke/name evidence units remain visible; only the eight documented executions are covered by the development visibility toggle. The four verified intentional compound imports, including HY-L126, were not reimported or modified.

Only source code, isolated fixtures/tests, and this report were changed. No production Database/Results data, peer state, service, deployment, synchronization, migration, or historical run configuration was modified. No push occurred.
