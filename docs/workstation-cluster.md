# Studio on multiple workstations

Studio can share its job catalog and route new workflows among workstations,
while each workstation retains its scientific files and local compute scheduler.
Workstation sharing is disabled by default. Installing an updated Studio does
not change existing scheduler jobs, allocations, paused campaigns, or Results folders.

## Components and ownership

- One Studio API also acts as the coordinator. Its private SQLite database is
  accessed through the API; it is never mounted as a shared database file.
- Every workstation, including the coordinator, runs a Studio background agent.
  Agents connect outward to the coordinator, publish resources and job metadata,
  and claim commands addressed to their permanent workstation UUID.
- The existing MN Compute Scheduler on **each workstation** remains responsible
  for CPU, GPU, RAM and scratch reservations across local MN applications.
- Jobs have a global UUID and an owning workstation. Local GPU numbers are scoped
  to that workstation. Identical run folders on two machines remain distinct.
- Existing result references remain local. Federated references use
  `cluster://WORKSTATION_UUID/TASK_GROUP/RUN_ID/relative/file` and refuse to
  resolve on another workstation.

New target files, sequences, settings and selected MSAs can be staged to the
chosen workstation. Existing Results files and model reference directories are
not bundled. Inputs use SHA-256 checks and a bounded, traversal-safe archive.
The staging limit is 512 MB per submission.

Analysis and follow-up sequence design/folding requests carry the source job's
machine affinity. A source on an offline workstation waits for that workstation.
The coordinator never moves an accepted scientific run to another machine.
Whole visual pipelines execute on one workstation, preserving all intermediate
files there. Requests combining existing source files from different workstations
are rejected; run each branch on its source workstation and compare summary metrics.

The shared catalog contains user/project IDs, job status, scores and a file
inventory. Full structures and logs stay on the owning workstation. A user can
request individual files for viewing/download; these requested coordinator copies
expire after 24 hours. Input bundles waiting for a workstation remain available.
Accepted input files are retained on the execution workstation for queued jobs,
restarts and provenance.

## Installation

Install the same Studio code and `mn-compute-scheduler` package on all machines.
Use each machine's own app workspace, Results folder and reference folder. Do not
copy a configured Studio workspace or its workstation UUID onto another machine.
Start and configure the local MN Compute Scheduler on each machine as usual.

The coordinator owns sign-in, users, project IDs and access grants. Connected
installations use those accounts and cache the project catalog. Create the
coordinator admin in production mode before exposing its API to a private LAN:

```bash
conda activate mn-molecular-design-studio
mn-molecular-design-studio create-admin
mn-molecular-design-studio cluster init --name current-workstation
```

The API listens on loopback by default. Existing SSH access can connect the
workstations without opening a new network-facing API port. On the coordinator:

```bash
mn-molecular-design-studio cluster install-services --enable
mn-molecular-design-studio cluster enroll --name 4032 --output studio-4032.json
mn-molecular-design-studio cluster enroll --name 4033 --output studio-4033.json
```

Transfer each enrollment file privately to its intended workstation. Each file
contains a credential for one node. On workstation 4032, replace the SSH destination
below with the actual coordinator SSH account and hostname:

```bash
mn-molecular-design-studio cluster install-tunnel --ssh-host USER@COORDINATOR --enable
mn-molecular-design-studio cluster join studio-4032.json --coordinator-url http://127.0.0.1:18766
mn-molecular-design-studio cluster install-services --enable
```

Repeat on 4033 using its own enrollment file. Port 18766 is the worker's SSH
tunnel to the coordinator. Port 8766 remains that worker's own Studio API.
Keep the SSH host-key checks and key-based access normally used on these machines.

For an authenticated private-LAN API instead of SSH, start the coordinator API
with `api --host PRIVATE_LISTEN_ADDRESS` and configure the workers with that
coordinator URL. Development authentication cannot bind to a non-loopback address.

The generated user systemd services use the Python executable from the active
Conda environment. They require the existing local scheduler service. User services
normally run during the user's session; use the workstation's usual systemd
linger policy if agents must stay available after logout.

Inspect services without starting scientific work:

```bash
mn-molecular-design-studio cluster status
systemctl --user status studio-api studio-agent studio-coordinator-tunnel
journalctl --user -u studio-agent -n 50
```

Use a single API process on each configured port. When switching from a manually
started Studio API to the generated service, stop that API process before enabling
the service. The local compute scheduler and its scientific workers keep running.

`cluster agent --once` publishes one snapshot and processes at most one assigned
command. It can execute an already queued command; use `cluster status` for a
strictly read-only check.

## User interface

Settings → Workstations shows connection status, free resources, the local
accept-new-jobs setting, and coordinator enrollment controls.

When sharing is enabled, the sidebar offers **Automatic**, **This workstation**,
and the other enrolled workstations, plus a resource summary for all signed-in
users. Automatic prefers the origin when compatible
resources are free; otherwise it chooses another ready workstation. It checks
adapter code revision, engine image identity where known, reference availability,
GPU memory, CPU, RAM, disk and reservations from the local scheduler.

My results becomes a paginated shared catalog with workstation/project filters and
separate job status and data availability. Local results retain their normal result view.
Remote result pages show synchronized scores and offer on-demand result inspection,
structure viewing, file downloads, cancellation, and follow-up analysis/redesign/
Boltz-2 prediction on the owning workstation.

Every API workflow and visual-pipeline submission accepts an optional execution
envelope, for example:

```json
{
  "schema_version": 1,
  "workflow": "analysis.ipsae",
  "project_id": "SHARED_PROJECT_ID",
  "source_job": "cluster:GLOBAL_JOB_UUID",
  "parameters": {"pae_cutoff": 10.0, "dist_cutoff": 10.0},
  "execution": {"placement": "auto", "idempotency_key": "client-generated-unique-submission-id"}
}
```

Use `placement: "local"` for ordinary local execution, including while the
coordinator is unavailable. An explicit remote source always enforces its own
workstation affinity. A repeated submission ID with the same request returns the
same global job; using it for different inputs is rejected.

Coordinator API routes are documented under `/docs`, including `/cluster/nodes`,
`/cluster/jobs`, `/cluster/submissions`, `/cluster/jobs/{id}/actions`, and command
status/file downloads. Node endpoints use separate workstation credentials.

## Offline and failure behavior

- Heartbeats normally arrive every 15 seconds. A node becomes unreachable after
  45 seconds and offline after 120 seconds. Completion status remains unchanged.
- Already running local jobs continue through coordinator/network outages. The
  agent publishes their eventual state after reconnection.
- The interface keeps a per-user cached catalog. Previously verified worker
  sign-in sessions may operate locally until the original session expiry; no new
  central identity or project permission is created while offline.
- An acknowledged job is never automatically retried on another workstation.
  Command lease recovery addresses the same workstation, whose receipt journal
  returns the original local submission result after a lost reply.
- If the agent crashes in the narrow interval between beginning local submission
  and persisting its result, the operation is marked **submission uncertain**.
  Check the local run carrying that `cluster_submission_id` before making another
  submission. This deliberately avoids silently duplicating an expensive run.
- A storage sentinel detects an unavailable Results mount and disables new work
  without creating a replacement storage identity. Missing files remain cataloged
  as unavailable; records are not erased because a workstation or disk is absent.
- Changing **Accept new jobs** drains the workstation from new automatic
  placement. Existing scheduler jobs continue.

The reference-file inventory verifies presence, not full content hashes of every
multi-gigabyte checkpoint. Container identities and code revisions are compared;
keep matching model distributions on the eligible machines. Existing node-local
user/project databases are not automatically merged into coordinator accounts.
Import their projects and map ownership explicitly when migrating older installs.

## Validation

Cluster tests run in the Studio Conda environment and simulate independent node
identities, busy resources, offline hosts, path separation, protected transfers,
shared accounts and lost replies. Integration tests exercise the real coordinator
API and workstation agent while replacing the scientific adapter with a small
local test job. They do not launch a model, SSH session or scheduler GPU job.
Real deployment still requires the selected coordinator and the actual SSH or
private-network addresses of workstations 4032 and 4033.
