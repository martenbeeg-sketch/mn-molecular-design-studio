# Project lifecycle and recovery

## User workflow

Open **Projects → Cleanup review**. The project overview is a cached navigation
aid; Prepare cleanup review independently checks authoritative SQL and run files.

1. **Archive project** hides it from the default Projects list and new-job
   selectors. Show archived projects exposes Restore. Existing Results and
   Database records remain available. Only the creator or an administrator can
   change visibility; an active cleanup cannot be unarchived to clear its lock.
2. **Clean selected runs** requires an explicit run selection. Selecting a logical
   run includes its recorded tool executions, including downstream analyses.
   Imported and explicitly saved collections are excluded unless their separate
   checkbox is selected. That checkbox includes the project's saved/imported
   collections, not just collections matching a name or timestamp.
3. **Delete whole project** reviews all registered executions. Excluded retained
   collections, active/waiting jobs, inconsistent ownership, unindexed owned run
   folders, uncovered workstations and retained raw-file dependencies block it.

Only an active administrator can prepare/confirm data removal. Confirmation
requires the exact reviewed project name and acknowledgement of recovery storage.
The action does not run during ordinary rendering or metadata refresh.

## Coordinated execution

`ProjectCoordinator` persists the immutable reviewed scope and per-node stage
receipts under `organization/project-lifecycle`. Each source saves its own plan
and operation journal. All trusted peers must answer the authenticated node
endpoint with `project-lifecycle.v1`, the expected node identity and the same
participant roster. Offline/older peers fail closed; legacy coordinator/worker
mode does not enable project-wide cleanup.

Begin revalidates file manifests, project metadata, SQL rows and Database reviews,
then establishes project/SQLite locks on every participant before any execution.
Current API submission entry points serialize acceptance with project changes
and reject archived, busy or deleted projects. SQL triggers reject stale member,
structure, source, observation and review writes, including older direct SQLite
clients. Ordinary per-run purge cannot bypass the active project lock.

Each source creates and integrity-checks a private SQLite recovery snapshot and
copies the selected retained structures with checksums. Existing recorded
structure hashes are checked, not replaced with the hash of changed bytes.
Referenced catalog databases are also snapshotted when available. Finished run
directories are atomically moved to the recovery store. Hard-linked regular
files are detached there; run checksums and structure receipts are verified
before SQL removal. A rename across incompatible filesystems fails rather than
silently deleting originals. External symlink targets are not recursively
deleted or copied as if they were owned by the project.

Only reviewed group memberships and their scoped scores/structures/sources are
removed. Canonical molecule identities and global CAS files remain, including
those used by other collections/projects. Explicit membership tombstones use the
existing source-authoritative removal outbox. A partial selection sharing one
automatic SQL cohort with unselected runs is blocked: the implementation does
not guess which unqualified legacy scores belong to which execution.

Every source must report data removal before finalization begins. Remote sources
finalize before the initiating workstation. Only then does whole-project deletion
write a hidden, persistent project tombstone. Peer replay cannot clear a busy
operation or resurrect a deleted identity. SQLite/path locks also prevent stale
indexers from recreating removed project data or new references to removed runs.

Interrupted operations retain their scope, ownership, locks and recovery receipts.
Use Continue confirmed cleanup; completed stages are not repeated. A changed
preparation can be aborted before data movement and must be reviewed again.
An operation remains visible in the pending-operations panel even if a process
stopped after writing its final project tombstone.

## Recovery and limits

- Storage is `<runs-parent>/.project-recovery/<operation-uuid>/`: SQLite snapshots,
  selected structure copies, moved run folders and a checksum receipt. Source
  projects/scientific identities remain qualified by their original IDs.
- Recovery copies still consume disk space. There is no automatic backup expiry,
  CAS garbage collection, database VACUUM, secure erasure or permanent trash purge.
- Restore checks use isolated fixtures. Restoring a deleted project's data in a
  live installation still requires an administrator-led, scoped procedure; never
  overwrite a running database with an old full snapshot, which would discard
  unrelated later edits. Archive has its own immediate Restore button.
- Shared SQL cohorts must be selected completely. Unknown legacy score lineage
  blocks partial cleanup. Data lacking a recorded project assignment is left
  unassigned, not matched by smoke/test names or one-entry run heuristics.
- Dependency discovery covers recorded upstream job identities, managed input /
  pipeline references and retained SQL structure paths. Arbitrary external scripts
  or unregistered references cannot be proven safe by this workflow.
- Distributed interruption, protocol mismatch and replay tests use isolated
  workstation simulations. Production peers must be upgraded and separately
  verified before enabling these actions on real research data.

The development preview keeps research copies read-only and offers an explicitly
disposable sandbox for exercising the same components. No scheduler or production
peer is connected to that sandbox.
