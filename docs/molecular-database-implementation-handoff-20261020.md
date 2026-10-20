# Molecular Database implementation handoff — 20 October 2026

## Status

The Database implementation and three-workstation metadata sync are live on `pc002975`, `pc004032`, and `pc004033`. No agents were used for this finalization. No scientific jobs were submitted.

Each machine now reports the same synchronized database totals: **159 molecules, 347 groups, and 295 observations**. The peer catalog contains 735 molecular records (347 group records and 388 molecule-entry records), with no replication conflicts. The local Database page was refreshed: Library shows 108 molecules, and All indexed shows 159. Library remains the default view; automatic 15-second rerendering is removed, with a manual refresh control available.

The API exposes the molecular database routes on all three machines, and the API and peer-agent services are active. Backfill completed with zero indexing errors on both remote workstations. A remote PDB structure was fetched from `pc002975` onto `pc004033` through the authenticated peer API and passed its SHA-256 check. Structure files stay on their source workstation until requested, then are cached for that user on the receiving workstation.

The synchronization rollout does not make the scientific Python environments identical. `pc002975` currently advertises code revision `4554a30a…`; `pc004032` and `pc004033` advertise `b0eea858…`. This mismatch already existed before this Database rollout because `pc002975` has additional uncommitted Python changes from other work. Database synchronization is verified, but an automatic scientific submission that requires matching `code_revision` may still reject the remote nodes until those broader code differences are addressed. No placement or scheduler policy was changed here.

## Tamarind model and Studio policy

See `molecular-database-review-20261020.md` and `molecular-database-behavior-recheck-20261020.md` for the investigation and sources. The account's empty Tamarind Database did not establish a manual-only policy; its documentation describes automatic molecular registration.

Studio stores canonical molecular identities once. Groups collect related molecules; structures, provenance, and scores stay attached to the relevant molecule and group context. Successful final design outputs enter Library groups. Inputs, evaluation results, and intermediate outputs remain discoverable in All indexed without creating a new Library collection for every calculation. Failed or empty runs are not automatically collected. Results pages also offer an explicit save of a chosen molecular subset to a selected or new group.

## Implementation

- `molecular_database.py`: SQLite/WAL identities, groups, memberships, structures, score observations, access checks, filters, and exact compound snapshots. Library and All indexed are separate views.
- `molecular_database_ingest.py`: successful-run indexing, legacy group migration, generation/evaluation classification, idempotent score and structure ingestion, and refresh throttling.
- `molecular_database_results.py` and `molecular_database_results_ui.py`: result-page molecule listing, exact-subset save-to-Database, and navigation back to source groups.
- `molecular_database_peer.py` and `cluster/replication.py`: peer metadata snapshots, reconciliation, read-only replica groups, access checks, and authenticated source-file retrieval.
- `molecular_database_api.py` and `api.py`: group/import/search/entity/schema/snapshot/structure routes and run-result save routes.
- `database_ui.py`, `molecule_imports.py`, `molecule_exports.py`, `molecule_schema.py`, `molecule_table.py`, `molecule_sketcher.py`, `database_selections.py`, and `docking_preparations.py`: Tamarind-style database/import UI, parsing/export, schema validation, selection, and tool handoff.
- Existing integrations also touch `run_app.py`, `vina_ui.py`, `cofolding_input.py`, `workspace.py`, pipeline modules, and workflow job completion. These files include unrelated changes from other chats; do not blanket-reset or commit the worktree.

The complete narrow rollout package is preserved under `/tmp/studio-database-rollout-20261020-8qaks_1r/`. It was based on clean revision `c95a6100a877c397d5f6459f38ea6b20feb133e8`, checksum `e956b41b315e3e97314e4c89fe4bd80fac8592c30bd9dc11f627f67abc18fbce`, and applied only to the clean remote checkouts. Local and remote source/SQLite backups are under each Studio home’s `backups/database-rollout-*` directory. No rollback was needed. Nothing was committed.

## Verification

The app Conda environment `mn-molecular-design-studio` ran the final focused database, peer, results UI, cluster, and pipeline suites: **121 passed, 1 warning** (an installed Starlette/httpx TestClient deprecation warning). These were focused suites, not the entire repository suite.

```bash
conda run --no-capture-output -n mn-molecular-design-studio python -m pytest -q \
  tests/test_molecular_database.py tests/test_molecular_database_integration.py \
  tests/test_molecular_database_policy.py tests/test_molecular_database_results_ui.py \
  tests/test_molecular_database_peer.py tests/test_cluster.py tests/test_cluster_integration.py \
  tests/test_pipeline_graph.py tests/test_sequence_pipeline_handoff.py
# 121 passed, 1 warning
```

The three live APIs advertise the new routes; agents on all machines are online; each replica converged to the same molecular-record counts with zero conflicts. The remote structure retrieval check passed. The Library and All indexed counts were confirmed in the live Studio UI.

## Remaining limitations

- Full-repository testing was not run. Import and export options were visually inspected, but no new import/export file was committed as part of this Database rollout.
- Existing formula/schema review leads (for example formula dependency cycles and imported-property sorting) remain separate follow-up work; no verified rollout blocker was found in this pass.
- The scheduler code-revision mismatch described above remains. No scientific task was dispatched as part of this rollout.
- Synchronized metadata is available on all three workstations, but source structures still require the owner workstation to be reachable for first retrieval. Cached structures remain available locally afterward.

Do not store account passwords, peer tokens, or browser authentication state in handoff files or logs.
