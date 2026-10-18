# Prediction queue audit — 16 October 2026

The audit followed the complete submission path: the molecular editor's submit
token, the Streamlit API client, REST validation/submission, workflow registration,
queue creation, central and local dispatch, worker job adoption, and both the
filesystem and SQLite result listings.

## Findings and fixes

- The shared workflow registration supplies `optimized_prediction` as its tool.
  Queue creation previously saved that generic value, although `params.engine`
  already identified the selected engine. The worker replaced it only when it
  adopted the run. Job creation now resolves `optimized_prediction_<engine>` and
  the user job name before dispatch. Submission responses include this identity.
- Prediction input preparation explicitly set an already dispatched job back to
  `queued`. It now records `preparing`, then `running` when prediction starts.
- API validation checked parameter types but did not apply the worker's molecule,
  mode, option, sample, and resource checks. The API and worker now share the same
  input preparation checks. Invalid API submissions fail before creating a run.
  Validation does not fetch MSAs, launch containers, or generate random seeds.
- The indexed listing could interpret a resumed `preparing` or `holding` job
  with an older failed result as failed. It now preserves those active/held states.
- Job JSON files were written in place. Writes now replace the file atomically,
  so concurrent readers cannot observe a partially written metadata document.

The latest warm smoke batch created 14 requested runs and the workers reused those
same run directories. The scheduler's existing idempotency check also prevents a
second dispatch for an already queued/running scheduler job. The editor consumes
each submit token once; rendering the same component snapshot does not resubmit.

## Existing records

The 13 remaining generic historical records were 12 cancelled runs and one failed
run. Their terminal states were confirmed against MN Compute Scheduler. Their
engine labels and job names were repaired, and stale scheduler status labels were
synchronized. Original JSON records were saved under
`/tmp/mds-queue-audit-20261016/legacy-records`.

The API was reloaded. Live checks rejected empty requests through both validation
and submission with HTTP 422, accepted a valid preflight request, and kept the job
count at 108. No additional prediction was submitted during these checks. The
final listing contained no generic `optimized_prediction` rows and no active jobs.

## Pytest verification

`python -m pytest -q`: **146 passed**, with one upstream Starlette/httpx
deprecation warning. Coverage includes all 16 registered prediction profiles,
rejected submissions creating no jobs, successful and failed preparation reusing
the original job, local/central dispatch, scheduler idempotency, legacy listing
compatibility, editor submit-token reuse, and atomic metadata writes. Prediction
containers and dispatch services are replaced with test doubles in these checks.

The initial suite had one stale ESMFold assertion: it expected seed and optimization
controls even though that engine exposes neither. The contract test now checks
each engine's declared capabilities.
