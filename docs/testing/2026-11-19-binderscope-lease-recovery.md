# BinderScope lease-loss recovery

## Incident evidence

Production job 9A53E (`20261119-102125-9a53e92f`) selected 505 saved PD-L1
RFD4 validation complexes in one 32-CPU job. The scheduler marked its first
attempt failed at 09:33 UTC after the Studio bridge reported repeated failed
API heartbeats. The local worker was terminated, but its Docker container
continued. The complete original report subsequently contained all 505 results.
The original worker never finalized its candidates or attached those results to
SQL. Studio metadata incorrectly continued to show `running`.

The bridge discarded the underlying heartbeat exceptions, so the exact original
transport/service failure cannot be reconstructed. Scheduler health/status calls
were responsive during investigation and there were no remaining reservations.
It would be unjustified to claim that CPU count alone caused the incident.

Confirmed implementation defects were: optional worker-registry heartbeats
preceded and could prevent job-lease renewal; a short four-second API timeout;
process-group termination did not stop daemon-owned Docker containers; failed
bridge executions did not get terminal Studio result files; and a completed
report had no verified finalization-only recovery path.

## Corrections

- Renew the job lease first, independently of optional worker registration.
  Allow a bounded 12-second response window, retain the existing failure grace,
  and initiate shutdown earlier when the recorded lease approaches expiry.
- Save heartbeat exception type/detail, failure count and timestamps, with the
  lease token redacted. A successful renewal clears the active failure detail.
- On heartbeat-loss shutdown, terminate the worker and stop only Docker
  containers with an exact bind mount of that run. Substring matches are not
  ownership. Failed cleanup is explicit and does not voluntarily release the
  scheduler claim. Infrastructure failure can still prevent confirmed cleanup;
  the operator must investigate that condition, rather than assume success.
- Write a failed Studio result rather than leave a phantom running job. A
  durable successful result wins over a late launcher/heartbeat failure.
- An ordinary retry of the same saved BinderScope job detects an existing
  complete report and finalizes it without Docker or input re-preparation.
  Incomplete or mismatched reports fail closed rather than silently recompute.

Recovery rechecks the current account, project and collection access, immutable
selection hash, exact candidate inventory (including duplicate detection),
structure/reference digests, per-design settings, complete report source fields,
upstream manifest and runtime adapter hashes. Persisted portable-path encoding
is compared in its original serialized form. Original calculation image and
protocol identities are retained even when the finalization wrapper changes.
A separate receipt records recovery, report hash, subject count, and that the
original container exit status was not observed; no success exit is fabricated.

The original job folder, report and consistent molecular/authentication database
snapshots are retained in the private production recovery directory. Production
authentication and design admission decisions are not changed.

## Verification

- 242 focused isolated regressions pass across scheduler/dependencies, workflow
  lifecycle, artifact validation, collection/admission, measurements and SQL
  attachment, including 45 focused recovery/bridge checks.
- All 67 BinderScope tests pass with the offline image enabled, including real
  queued execution followed by same-folder report-only recovery and idempotent
  attachment to the original SQL design.
- Recovery rejects missing, extra or duplicate subjects, replaced SQL selections,
  changed structures, changed settings and mismatched runtime provenance.
- Bridge regressions cover failed optional worker heartbeat with a healthy job
  lease, lease-loss container cleanup and terminal status, token redaction,
  exact-mount isolation and a successful result followed by a heartbeat error.
- The original 505-result production report independently matches all original
  request fields, authorized source selections and input-file digests.

These checks do not simulate an extended production network partition, daemon
failure or an OS kill of the entire bridge. They do not guarantee that future
infrastructure failures are impossible; verified recovery and diagnostics make
such failures visible and avoid unnecessary repeated scientific calculations.
