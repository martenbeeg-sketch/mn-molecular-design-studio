# Projects and cleanup review

Projects is the default landing page. All tools retains its original separate
position at the top of the sidebar. Below its divider, the order is Projects,
Pipelines, My results, Database, Complex Analysis, then account/administrative
controls. Existing job-detail routes and tool back links keep their original
destinations.

## Project workspace

Search and order projects, then select one project. The workspace displays:

- Logical runs, with separate counts of finished executions and active/unknown
  executions. A pipeline or distributed run can contain several executions.
- A paginated run overview, including execution counts and workstations.
- Separately launched refolding and analysis executions join their recorded
  upstream Run, rather than appearing as unrelated single-job runs. The collapsed
  **Linked executions & cancellation** panel shows each execution and its source.
- Saved/imported collections in a collapsed section. SQL collection counts are
  memberships, not a count of admitted designs or passing candidates.
- Links to project-scoped Results and Database. Opening Results clears previous
  selections, confirmations and search filters before selecting the project.
  Opening Database resets stale browse/result filters and uses the normal
  admitted-results Library view with all accessible collections in that project.

The inventory reads indexed job metadata and SQL group counts, not coordinates
or workdir discovery. Local and remote rows are qualified by workstation and
run reference; a local copy preserves the linked cluster ID. Hidden executions
remain in review counts. The session cache refreshes after 60 seconds or a manual
Refresh, and invalidates when the visible project IDs or user/role changes.
Current permission and run-state checks still occur in existing purge endpoints;
the cached overview is never an authorization or deletion receipt.

## Cleanup boundary

Projects first creates a **review**, never an automatic deletion cascade.
Its downloadable JSON identifies the project, exact finished execution references,
blocked executions, protected collections, generation time and catalog completeness.
No candidate is automatically selected because its name contains smoke/testing,
because it has one entry, or because old acceptance evidence is absent.

Projects now provides three explicit actions (see [Project lifecycle](project-lifecycle.md)):

- Archive / Restore changes project visibility without removing files or SQL.
- Clean selected runs removes a reviewed, dependency-complete scope of finished
  executions and automatic SQL memberships, retaining the project.
- Delete whole project includes all registered project executions and requires
  explicit inclusion of saved/imported collections before completion.
- Data removal is administrator-only and requires the exact project name and
  acknowledgement of the automatic verified recovery backup. Every trusted peer
  must implement `project-lifecycle.v1` and confirm the same workstation roster.
- Shared molecular identities and global structure files are retained. Shared
  cohorts, incomplete lineage, missing peers and active jobs block cleanup.

Deleting a folder manually does **not** reliably remove SQL records, saved
structures, catalog entries or replicated copies. Back up SQLite and referenced
files before any approved purge; validate restore separately. Never use a folder
name, project name or raw SQL deletion as a substitute for provenance checks.

Recovery files remain in `.project-recovery` beside the managed runs directory;
cleanup does not automatically reclaim their disk space. A deleted project keeps
a hidden identity tombstone. This is not secure erasure, automatic garbage
collection, or a one-click production restore workflow. Unassigned legacy data
is never attributed to a project by name alone.

## Dependency-aware cancellation

**Cancel selected job + downstream** stops the selected execution and its
transitive recorded descendants, including queued or paused work. It traverses
completed intermediate steps without changing their saved results. Cancelling a
downstream analysis does not cancel its upstream design or an independent branch.
Connections are qualified by project, owner, workstation and execution identity;
display names, identical sequences and shared logical Run labels are not grounds
for cancelling a job.

Local workers check durable upstream cancellation at dispatch and immediately
before execution. Both scheduler bridges also enforce the durable cancellation
receipt if it races worker startup/PID publication. The per-run
`cancellation_request.json` survives an older worker rewriting ordinary metadata;
only explicit resume clears it. Shared catalogs retain cancellation requests across restart,
apply them to later submissions/heartbeats and queue idempotent stop commands on
owning workstations. A remote request is not proof that a process has stopped:
offline computers can only acknowledge cancellation after reconnecting. All
participants must run compatible updated application code for late-child guards.
Legacy records without an unambiguous source/parent link cannot safely be
included automatically; their scope must be reviewed rather than guessed.
Dependency reconciliation uses compact metadata and scoped cancellation receipts;
an idle dispatch does not rescan old cancellation histories or load structures.

Cancellation is separate from cleanup: it does not delete run files, SQL entries,
accepted designs or completed evaluation results. To deliberately restart local
cancelled work, use the existing explicit resume workflow; other cancelled
executions remain cancelled until separately resumed.

## Development review

The focused development preview uses a private SQL copy and a read-only export
of scalar local job metadata plus the last-known workstation catalog. Those
snapshots are not transactionally simultaneous and are marked incomplete.
Projects and Database use the application components; the preview Results list
is read-only. Other sidebar destinations are explicitly unavailable there rather
than connecting to production services. Normal sign-in is retained; production
credentials and authentication settings are not changed.

The research snapshot remains read-only. Select **Disposable cleanup sandbox**
in the sidebar to exercise Archive, Clean and Delete with simulated data in a
separate private SQLite database and run tree. No real research cleanup,
production deployment, credential change or peer operation occurs during review.
