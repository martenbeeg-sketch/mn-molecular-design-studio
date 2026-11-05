# P0 F4/F5 correction report

Date: 2026-11-05
Scope: the remaining F4 and F5 findings from the independent review of `82ef6fd591443071a47213faafa320c7dec013bc`.

## F4 — exact selected Evaluation identity

The prior duplicate check joined an observation to any source row sharing its group and canonical molecule. That allowed an unrelated measurement in the same group, or a foreign-authoritative replica retaining colliding hostname/source labels, to suppress the selected local Evaluation.

The check now derives the expected observation and source IDs from the historical publisher's exact source key. It accepts an existing local observation only when the group, observation ID, source key, source node, run and compound match. For a replica, it also verifies the derived replica group, authoritative replica node, original observation ID, replica source ID, source key, and replica markers. An unmatched or ambiguous record does not suppress the selected measurement; historical records are left intact.

Permanent tests cover wrong observations in the right group, same-source replay, local-authority replica reuse, a foreign-authority replica with colliding labels, and the private foreign-owner ACL boundary. The selected score and source association are checked, and replay remains idempotent.

## F5 — per-item Vina import through batch children

Non-coalesced batch children lacked the trusted Studio owner/project marker required by completion ingestion. Their Vina CSVs could complete successfully while the explicit Database import was skipped. Importing the aggregate candidates instead would also risk broadening selection or duplicating child ingestion.

The batch coordinator now reads authorized owner/project lineage from its persisted Studio marker and copies only those trusted values to each child marker and metadata. User-supplied item fields cannot grant this authorization. Each child retains its serialized effective import policy. For non-coalesced batches, the parent remains Results-only so only explicitly automatic children ingest; a coalesced all-true job continues through its existing automatic completion path.

The permanent real-path regression covers omitted, false, true, mixed, and reversed choices, with compatible and differing Vina settings. It exercises request submission, persisted worker-request restoration, batch execution, the Vina adapter, output staging, completion ingestion, and retry. Process spawning is run synchronously in the test; Docker execution alone is replaced with deterministic score and pose files. Assertions cover exact observation counts and scores, canonical ID reuse, source/run/settings/owner/project attribution, retained CSV and pose bytes, and retry stability.

## Verification

Before correction, the independent review reproduced F4 as **2 failures** and F5 as **9 failures, 5 passes**. After correction, all **58 independent review cases pass**. The focused permanent F4/F5 selection passes **17 tests**; the additional pose-retention assertion passes all **14 Vina batch matrix cases**.

The full 25-module P0 regression selection passes **586 tests** using the application's Conda interpreter:

```text
/home/user/mambaforge/envs/mn-molecular-design-studio/bin/python
586 passed in 131.87s
```

All scientific writes occurred in temporary fixture databases and run directories. No production database, datasets, peer state, services, or deployed processes were accessed or changed. No migration, synchronization, restart, deployment, push, or Phase 2B work was performed. Docker execution and live inference quality are outside these deterministic ingestion tests.
