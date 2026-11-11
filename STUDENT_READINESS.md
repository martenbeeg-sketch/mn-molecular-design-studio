# Student readiness — supervised pilot, not production rollout

Date: 2026-11-11. Application candidate:
`f0e768235eff799d7abbe06a07c4970ca5d05e78`.
Tracked source archive SHA-256:
`c4f90bddd91be0d6e51320cb5a96843fd1ab19199b6fe58b4bd2fc4304bdcf37`.

## Current decision

**GO for instructor-supervised evaluation in the isolated teaching pilot.**
**NO-GO for general student access to the research production system.**

The user authorized preparing the recovery, validation, account and teaching
checks. Production services, production authentication, scientific records and
the 436 cleanup candidates have not been changed. The old research deployment
remains `p11-25e003c/source` on all three hosts. Only the independent teaching
container, its loopback relay and its separate institute-LAN HTTPS entry point
are newly running. The agreed class size is three students and one instructor;
they connect through browsers on the institute network. See `TEACHING_ACCESS.md`
for the URL, one-time certificate setup and account handoff.

## Open the pilot

- UI: <http://127.0.0.1:18502/>
- Institute-LAN UI: <https://pc002975.marionegri.it:18502/>, after trusting the
  separate teaching certificate as described in `TEACHING_ACCESS.md`.
- API: <http://127.0.0.1:18766/api/v1/health>
- Current release/state directory:
  `/home/user/programs/deployments/mn-molecular-design-studio/teaching-pilot-current`
- Private demonstration credentials: `teaching-state/credentials.json` inside
  that directory, mode `0600`, inside private `0700` release/state directories.
  Do not publish, commit, email or reuse this file for real students.

The accounts are `teaching_instructor`, `teaching_student_one`,
`teaching_student_two` and `teaching_student_three`. They are four distinct
teaching seats, not verified real-person identities. Acceptance exercised the
students' mandatory temporary-password change and session revocation. All three
student passwords were then reset to new individual temporary passwords in the
private file, with a fresh mandatory change for first human sign-in. Assign one
seat per student and record that mapping; do not share demonstration logins.

The pilot contains two instructional small molecules and a public structure
already shipped with the repository. It contains no copied production library,
campaign, student personal information or important research data. The selected
chain's sequence-liability analysis is a real native CPU calculation, not a
new structure prediction or evidence of affinity or experimental success.

The browser sign-in page was visibly verified. Authenticated journeys were
exercised over real HTTP, plus maintained Streamlit UI regressions; this is not
a completed manual review of every authenticated browser page or teaching lesson.

## Isolation and resource boundaries

- Container: `mds-teaching-pilot-20261111`, immutable image
  `sha256:b3366c6f1e14ffe3d0962938440503a3126685645d7f9a3fee0e81db6b57691a`.
- The host Docker daemon is rootless. Container UID 0 maps to the ordinary host
  user, not host root; all capabilities are dropped and no-new-privileges is set.
- No container network, GPUs, Docker socket, production Database, production
  Results directory or production cluster state is mounted into the pilot.
- Candidate source and the installed Python runtime are read-only. Only separate
  teaching state and temporary storage are writable.
- Actual cgroup limits were checked: 8 GiB RAM, four CPU cores and 512 processes.
  These are **container-wide limits, not per-student quotas or fair scheduling**.
- Persistent teaching-state storage has no filesystem or per-student disk quota.
  Supervise small examples, monitor free space, and do not enable large uploads
  or campaigns. A classroom storage budget remains a production gate.
- A bounded host-side relay forwards only loopback ports 18502 and 18766 into
  the network-disabled container. A separate Caddy HTTPS gateway now binds only
  `172.18.24.170:18502`, allows the connected institute subnet `172.18.24.0/22`
  and forwards only the teaching UI. The API remains loopback-only. No firewall
  changes, public-Internet listener or production proxy changes were made.
- No production service definitions were changed. The pilot is not configured
  for automatic startup after reboot.
- Pilot Streamlit usage telemetry is disabled and the toolbar uses viewer mode.

Use only the prescribed browsing/import/native-CPU/review/export exercises.
GPU prediction, design campaigns and tools needing the host Docker daemon are
not available in this container. Do not interpret an unavailable tool or missing
measurement as a failed scientific result or a measured zero.

Do not give students the shared `compute`/`user` workstation login or Docker
access as a substitute for application accounts. That would bypass filesystem
isolation and expose private operator state. Class access is through the
separate institute-LAN TLS front end, for the confirmed three students and
instructor. A real-person seat mapping and first authenticated-browser lesson
are still required. Students in other institute VLANs may need a separately
approved source-subnet rule; network membership is not a replacement for login.

## Defects discovered and corrected during acceptance

### Immutable deployment request staging

Real tests against a read-only source mount found that default admin API requests
created temporary request JSON in the deployed source directory. Both ordinary
users and admins now default to their own private application request directory.
Explicit server-side `base_dir` remains admin-only. Existing admin clients using
relative server paths should specify `base_dir` explicitly or use owned uploaded
file references; a caller-selected directory still must support request staging.

### Post-ingestion scientific admission identity

The first pilot exposed repeat work after ingestion enriched descriptive
component roles. The same explicit structure/model/chains acquired a different
admission fingerprint; this was observed in a restart/replay exercise.

Admission v2 excludes only descriptive roles inside
`database_identity_components`. Actual binder/target chain selections, sequences,
reference roles, scientific settings and protocol versions remain significant.
New reservations use versioned `v2:` keys. Existing v1 reservations are checked
against their retained input manifest and original request fingerprint; completed
matches reuse the earliest recorded run deterministically without deleting old
records. Pending/uncertain matches remain blocked. Unverifiable legacy inputs
fail closed and require operator reconciliation, not an automatic retry.

The old experimental state and its two duplicate trial jobs were preserved.
An in-place v1-to-v2 test reused the original run without adding a third job.
The final fresh teaching state has one native assessment across two student
sessions, ingestion and a real service restart. This trial-state replacement
did not clean or migrate production data. Production currently has zero
`scientific` reservation rows on all three hosts.

Admission-v1 writers do not understand v2 keys. Downgrade requires keeping
submissions closed; do not delete v2 rows to make an old binary submit again.
Upgrade every writer before reopening admission.

## Validation evidence

| Check | Result |
|---|---|
| Final full local suite at the candidate | 3,099 passed, 37 skipped, 103 warnings, 899.33 seconds |
| PC004032 isolated release/API/peer/review/student integration | 367 passed, 4 warnings, 53.44 seconds |
| PC004033 isolated release/API/peer/review/student integration | 367 passed, 4 warnings, 52.86 seconds |
| Focused native role-enrichment and legacy admission regressions | 61 passed |
| Focused cluster-account and teaching-tool checks | 35 passed |
| Live authenticated pilot HTTP journey | 30 checks passed |
| Native saved Evaluation and review/export | Exact method version, run, structure digest and Evaluation ID retained; stale review rejected |
| Fresh pilot complete offline data restore | 61 regular files/links and four SQLite Databases verified; Evaluation identity preserved |
| Pilot restart/replay | Original run reused, one analysis job |
| Institute LAN follow-up | 73 focused isolated tests passed; four distinct sign-in sessions rendered over verified TLS WebSockets from each of three PCs |
| Four teaching seats | Four concurrent authenticated reads; students denied admin; temporary passwords reset for first human sign-in |
| Production Database-only snapshots/restores | Six per host, 18 total; integrity and foreign-key checks passed |

Remote tests use the installed workstation Python runtime in a read-only,
network-disabled container, not a host-security-policy workaround. The inspected
test image is `sha256:9b9845b375dc5589a69102b574cc01e73b3df3b96543f313c77ff4c2a1566d12`.
No production storage or Docker socket is mounted. Bubblewrap remains unavailable
on the peers; it was not weakened or bypassed. Local broad tests use Bubblewrap.
The local Python environment's pytest is installed in user site-packages, outside
its narrowly mounted Conda runtime, so an optional local container-test attempt
could not import pytest; the successful local suite uses the existing isolated
launcher instead. No runtime packages were installed or upgraded.

Real CPU containers for Sequence QC, Biotite Interface and Interface Clustering
each completed two retained-fixture items on the primary and PC004033. Their QC
image IDs differ; this does not certify image equivalence. PC004032 lacks
`mn-protein-qc:1.0.0`; no image was downloaded, rebuilt or silently substituted.
Live GPU inference and every scientific engine remain outside this acceptance.

Earlier full runs and failed probes remain diagnostic history: one broad run
exposed five mocked account IDs (`owner`) that violate the real account-ID
contract; the cluster fixtures now use UUID identities, not relaxed validation.
A subsequent pre-admission-v2 broad run passed 3,094 tests with 37 skips. Neither
that result nor earlier partial probes replaces the final candidate run.

Final JUnit/logs, sanitized backup summaries, immutable-source checksums, live
pilot/recovery receipts and CPU-container evidence are packaged in
`STUDENT_READINESS_EVIDENCE.tar.gz`, SHA-256
`7c6b46c32a331fc8d15fdd8c0a119a9a350cd9196999c49d9ee1adfab4d9e155`.
It excludes all scientific Database files, private credentials and complete
teaching-state/recovery directories.

## Recovery boundaries and remaining production gates

The new production snapshot sets are private under each host's own
`~/mn-mds-recovery/20261111T183255Z-student-readiness-online`. They were created
with read-only-source SQLite backup connections and verified in separate restore
copies. They are **online per-Database snapshots**, not a coordinated quiet-point
recovery set for all related Databases, scientific files, queues and replicas.
They were not substituted between workstations or included in public evidence.

The teaching pilot separately has a complete quiet-point data copy, consistent
SQLite snapshots and an offline read-only restore verification using the pinned
runtime. This is recovery evidence for the disposable pilot only. It does not
clear research-production recovery gates or constitute a clean-machine restore.

Before production activation or unsupervised class use:

1. The class size (three students plus instructor) and institute-browser access
   are confirmed. Finish the roster/seat mapping, certificate trust on the actual
   laptops, teaching workflows and resource budgets; validate authenticated
   browser journeys and scheduling at that intended scale.
2. Review the production account/ownership transition. All three production
   account stores currently have zero real users and explicitly use development
   authentication. Preserve existing owner/project IDs; do not just flip the
   authentication environment or assign all historical data to a new UUID.
3. Reconcile the retained legacy queued/running/collecting records without
   treating age, smoke names or single-member groups as deletion authorization.
4. Complete a coordinated per-host scientific-file/queue/replica recovery set and
   restore rehearsal in the documented maintenance window. Existing large
   historical archives are retained; no terabyte-scale file recovery is claimed
   by the new Database-only snapshots.
5. Establish pinned runtime/image availability for each required workflow and
   execute representative real GPU/design jobs through the intended scheduler.
6. Deploy all writers together, verify runtime provenance and derived index
   reconciliation, then perform an instructor-led acceptance lesson before
   opening student admission. Follow `docs/database-simplification-p0-deployment.md`.

AUD-008's pre-pagination filtering, unknown historical links and other limitations
in `INDEPENDENT_SOL_REMEDIATION.md` remain. The original independent baseline audit
and its evidence have not been rewritten as a certification of these later fixes.

## Stop only the teaching pilot

`docker stop mds-teaching-pilot-20261111` stops the pilot without stopping any
research service. Keep its state and recovery directories. The separate loopback
relay's PID file is `teaching-relay.pid` in the current teaching release directory;
verify that PID's command line before stopping it. A restart uses the same retained
state; do not rerun bootstrap against the production account store. To close
only LAN access, run `systemctl --user stop mds-teaching-lan.service`.
