# Independent Studio workstations

Studio supports independent installations that exchange account, project, resource
and job metadata. Each installation has its own API, account database, local compute
scheduler, job catalog and Results storage. In **Independent workstations** mode
(`mode=peer`), no particular workstation must remain switched on.

## Normal operation and outages

- Sign in, create projects, submit local jobs and analyze local results even when
  every other workstation is offline. Sessions are issued and checked locally.
- Each background agent publishes its local catalog into its local shared store,
  then synchronizes with every configured reachable peer. Synchronization has
  durable cursors and resumes after an outage or process restart.
- Automatic placement prefers this workstation when compatible resources are
  available. Otherwise it selects a compatible peer using its latest resource
  report. Explicit **This workstation** submissions use the local queue.
- Each originating installation keeps its own durable submission queue. Execution
  machines poll their peers directly for assigned commands. A third workstation
  does not relay or authorize this communication.
- A queued automatic command can be reassigned if its selected machine disappears
  before claiming it. Claiming and reassignment use one database transaction.
  Once claimed, it stays on that execution machine, including after a lost reply.
  Durable execution receipts return the original submission result on retries.
- A submitted request that still needs inputs or acceptance from an offline origin
  waits for that origin to return. It is not duplicated elsewhere. Retry the same
  submission ID against the same originating installation.
- Direct liveness pulses normally run every 15 seconds, independently of job-catalog
  reconciliation and metadata replication. After 45 seconds without direct contact
  a peer is unreachable; after 120 seconds it is offline. Synchronizing an old
  catalog entry through another peer does not make the unavailable machine online.
- Every catalog retains the machine name, original path and last known job state.
  An offline machine's data remain listed and are marked unavailable.

## Scientific data and scheduler ownership

Full results, structures, logs, checkpoints and alignment databases stay local.
Metadata replication never opens another machine's filesystem paths. Source files
use `cluster://WORKSTATION_UUID/TASK_GROUP/RUN_ID/relative/file` references.

Analysis and follow-up redesign/folding requests run on the workstation that owns
those results. Whole visual pipelines execute there too. A disconnected source
waits for its owner. Inputs from different machines cannot silently be combined
into one local run. Requested individual files may be temporarily downloaded for
viewing; those copies expire after 24 hours.

New input bundles are bounded to 512 MB and checked by SHA-256. Pending submission
inputs are retained until accepted. A storage sentinel disables new execution when
a Results disk is unavailable; it never creates a replacement storage identity.

Each local MN Compute Scheduler continues to own GPU/CPU/RAM reservations. Studio
peer migration does not resume paused campaigns, change scheduler allocations or
cancel existing scientific jobs. Matching engine code, container identities and
required reference files are checked before remote execution.

## Accounts, projects and trust

Only explicitly enrolled, trusted workstations participate. Their private machine
credentials authenticate peer APIs over SSH or a trusted authenticated network.
User API credentials cannot read private synchronization records.

Account IDs, password hashes, disabled-account flags, project IDs and memberships
are replicated. Passwords and session tokens are not copied. A user can sign in
locally with the same account after its account record has first synchronized.
A password or disabled-account update invalidates that user's existing sessions
on a peer when the update reaches it.

Changes use Lamport versions and a stable workstation ID tie breaker, rather than
wall clocks, so different system clock settings do not choose a stale update.
Concurrent edits of one record converge deterministically. Membership removals
and project removals are retained as tombstones and cannot be resurrected by old
cached records. Different users created with the same username on disconnected
machines remain separate identities: the collision is retained for administrator
review under Settings → Workstations. Jobs cannot change execution-machine or
user/project ownership through synchronization.

Account and permission revocations cannot reach disconnected machines until they
reconnect. These machines are trusted copies of the organization's account data;
protect their disks and machine credentials accordingly. Application sessions
remain local. No SQLite database is mounted or continuously copied between hosts.

## Initial setup and migration

Install Studio and the local MN Compute Scheduler on each machine. Keep permanent
workstation UUIDs, local Results directories and reference paths distinct. Create
and initialize an independent first installation with the Studio Conda Python:

```bash
conda activate mn-molecular-design-studio
mn-molecular-design-studio create-admin
mn-molecular-design-studio cluster init --name workstation-name --independent
mn-molecular-design-studio cluster install-services --enable
```

For enrolled installations, a private, administrator-prepared manifest migrates
existing node identities and credentials:

```bash
mn-molecular-design-studio cluster enable-peers /private/workstation-peers.json
systemctl --user restart studio-api studio-agent
mn-molecular-design-studio cluster status
```

The manifest has schema version 1, the receiving `node_id`, its own `api_url`, a
`members` array of enrolled IDs/names/credential hashes/disabled flags, and a
`peers` array of IDs/names/API URLs as reachable **from that machine**. Optional
`identity_records` seed existing accounts/projects using the replication format.
Transfer this manifest only through an already authenticated private channel.
The receiving node's existing private credential must match its membership hash.

Migration uses consistent SQLite backups in `cluster/before-peers-*` and copies the
previous configuration and organization metadata. It retains the workstation UUID,
credential, project/job identities, Results path and scheduler. Review synchronization
conflicts before removing those backups.

## Direct SSH connectivity

Each pair of workstations needs a connection that works while the third is off.
With loopback APIs on port 8766, a direct SSH forward on 4032 can expose 4033 as:

```bash
ssh -N -o BatchMode=yes -o ExitOnForwardFailure=yes \
  -o StrictHostKeyChecking=yes -o ServerAliveInterval=30 \
  -L 127.0.0.1:18763:127.0.0.1:8766 compute@pc004033.marionegri.it
```

Configure 4033's peer URL as `http://127.0.0.1:18763` **on 4032**. Set up the
opposite direction too. Keep these tunnels in restartable systemd user services.
Dedicated SSH keys can be restricted to forwarding to `127.0.0.1:8766`, with no
shell, PTY, agent or X11 forwarding. Validate host keys using an existing trusted
SSH connection; do not disable host-key verification.

On pc002975, the existing tunnels to 4032/4033 can retain reverse forwards so the
workers can contact pc002975 without exposing its API. Those links are optional
when pc002975 is off: the direct 4032↔4033 links remain available.

Inspect without starting a scientific job:

```bash
mn-molecular-design-studio cluster status
systemctl --user status studio-api studio-agent
journalctl --user -u studio-agent -n 50
```

User services require the usual login/linger policy to remain available after
logout. Keep API listeners on loopback when using development authentication.

## Legacy deployments

`coordinator` and `worker` modes remain supported for compatibility. In those modes
shared sign-in, projects and dispatch require one coordinator; cached sessions only
work until their original expiry. `cluster init` without `--independent`, `enroll`,
`join` and `install-tunnel` retain their original behavior. Migrate enrolled nodes
with `enable-peers` to remove that dependency.

## Verification

The peer tests exercise separate real SQLite account/catalog databases, HTTP peer
APIs, local authentication with no network, offline project creation, metadata
catch-up, concurrent edits, permission tombstones, ownership checks, private peer
credentials, durable cursors and duplicate submission receipts. Model execution is
replaced with a recorded submission, so these checks do not use a GPU or start a
scientific workflow.
