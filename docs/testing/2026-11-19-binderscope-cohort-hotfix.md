# Single-job BinderScope submission and browser cache hotfix

## Incident and correction

The 505-design submission timed out while preparing inputs synchronously, before
the job was queued. Validation and submission each prepared the full cohort.
Preparation also scanned Results for SQL-only selections and repeatedly hydrated
the molecular catalog. Increasing calculation CPUs does not fix this HTTP path.

Saved SQL cohorts now receive inexpensive, batched identity and access checks at
admission. The API claims the existing idempotency record before queue creation.
One normal background job prepares and calculates the entire cohort, retaining
the requested CPU allocation (including 32 CPUs), project and owner. Its status
reports preparation progress. Structural integrity, chain roles, target context,
account availability, access and cancellation are checked again by that worker.
The worker adopts the original queued folder: it does not create per-design jobs.
Attached analyses use automatic collection so they attach to their original SQL
designs and representations without admitting new designs.

The browser now uses a logical SQL change counter rather than invalidating its
entire index for unrelated database bookkeeping. Triggers include real changes
from peer replay and exclude no-op updates. Fresh ACLs remain part of every cache
lookup. Two bounded metadata snapshots can be reused across sessions with the
same user, role, filters, database identity, revision and accessible groups;
per-session navigation dictionaries are not shared.

## Verification

- The isolated API, browser, SQL identity/policy/peer, source adapter and
  scheduler-lifecycle regression run passes: 599 passed, eight skipped and six
  warnings in 134.28 seconds. The eight skips are covered by the enabled image run.
- The 505-selector API regression queues exactly one 32-CPU/no-GPU job, avoids
  structure preparation and full catalog hydration in the HTTP request, and
  returns the same job on an exact idempotent retry. A conflicting retry fails.
- Worker tests cover same-folder preparation/calculation, changed selection
  content, disabled accounts, revoked collection access and cancellation.
- With the offline image enabled, all 57 BinderScope submission/core/hotspot
  tests pass, including a real queued-worker calculation and ingestion into the
  original SQL design, as well as replicated-design and target-context cases.
- Cache tests cover actual edits, unchanged writes, unrelated bookkeeping,
  access revocation, database replacement and cross-session isolation.
- On an isolated copy of the real catalog (17,716 indexed molecular records),
  metadata index lookup measured 5.27 seconds cold, 0.023 seconds warm, 0.017
  seconds in a new session and 0.016 seconds after an unrelated SQL write.
  These are metadata timings, not complete browser-rendering measurements.

## Limits and activation

The real 505-design scientific cohort is not calculated as a smoke test. The
regression cohort uses synthetic saved structures; the real container test uses
one synthetic complex and one CPU. No claim is made about total calculation time
or full application/GPU pipeline coverage. Cold loading still has to read and
organize the catalog, and actual scientific changes correctly invalidate it.

Install one immutable committed source archive on the three existing machines,
preserving production authentication, all scientific data, ports and node
identity. Recheck active work before restarting services. Retain the prior full
recovery copies and add fresh consistent database/configuration snapshots for
this additive bookkeeping migration. Keep the main UI on port 8502; do not revive
the retired read-only preview or automatically queue research analyses.
