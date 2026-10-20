# Database behavior recheck — 20 October 2026

## What was checked again

The authenticated Tamarind account shows three completed jobs in My Results: an AutoDock Vina example, a Boltz-2 test, and a BindCraft2 batch. The Database remained empty after selecting Organization, all projects, all users, and all four import categories. The Vina and Boltz result pages expose structures, score tables, files, downloads, and downstream tools. No explicit Save to Database control was visible on the result views inspected. No jobs were submitted or account records imported during this investigation.

Initially, the Database filter selected only Tool inputs and User uploads, hiding Tool output. Explicitly enabling all four categories still returned no groups. The empty view therefore cannot be explained solely by those initial filters.

This observation does **not** establish that Tamarind intentionally excludes completed results or requires manual promotion. Its current documentation explicitly describes automatic molecular registration. The reason these particular account jobs are absent could not be established from the visible interface. Account/tool integration and historical backfill behavior remain unverified.

Sources:

- [Database documentation](https://docs.tamarind.bio/tamarind/database)
- [Molecules API concepts](https://app.tamarind.bio/api-docs/molecules-concepts)
- Authenticated [Database](https://app.tamarind.bio/database) and [My Results](https://app.tamarind.bio/app/results)

## Tamarind's stated purpose

The database keeps molecular identities and their accumulated scientific data available for reuse. A molecule may belong to several groups. Reimporting the same composition retains its existing identity and scores. Structures are representations associated with molecules and labeled by source.

The documented distinction is:

| Operation | Documented database behavior |
| --- | --- |
| Import molecules | Create or reuse molecular identities and place them in a group. |
| Generate/design new molecules | Automatically save molecules and create an output group. |
| Predict structures or calculate scores | Attach representations and scores to the input molecules/groups; create an input group if needed. |
| Run a job | Keep execution history and downloadable output on the job page. Completion alone does not mean that every artifact is a separate molecule. |

For example, repeated predictions for one protein should produce multiple source-labeled structures on that protein record. Docking one compound against a target should preserve the receptor, pocket, pose, settings, and source run as score context. It should not create a chemically new compound for each pose or every time a scoring tool runs.

Export is a snapshot for use outside the application. A persistent database additionally supports identity, accumulated scores, comparison across tools, reusable selections, provenance, and shared access.

## Studio's current difference

`molecular_database_ingest.ingest_run()` creates an output group before checking whether a run contains supported molecules. `sync_workspace()` scans accessible jobs, and completion hooks also invoke ingestion. Prediction/evaluation runs can therefore create their own collections, and runs without usable molecular output leave empty groups.

The local library inspected during this recheck contained **152 molecular identities, 306 groups, 77 empty groups, and 282 observations**. These counts are a local implementation finding, not Tamarind behavior. Existing identity deduplication already avoids some repeated molecule rows; collection creation and visibility still make the page resemble the job catalog.

The 15-second timer around the Database renderer was removed. Page entry, ordinary widget interactions, and manual refresh still update its view. The isolated import-wizard UI regression passed in the application's Conda environment.

## Recommended Studio policy

Keep a complete execution history in My Results. Use the Database primarily for reusable molecular collections and comparisons. Implement the following behavior without discarding existing indexed records:

1. Uploaded and explicitly imported molecules belong in the library immediately.
2. Generative tools publish final molecular output groups, scoped to their campaign/project. Keep acceptance/rejection state and allow rejected or intermediate candidates to be hidden by default.
3. Prediction, docking, rescoring, and affinity tools attach structures and contextual scores to their input molecular records. Reuse a source group when one exists, rather than creating a fresh output group for every calculation.
4. Provide **Save selected to Database** on result pages for choosing a subset and its destination group. Provide **Open in Database** when results are already associated. This makes selective reuse clear without requiring export and reimport.
5. Separate curated groups from automatically indexed collections in the default view. Keep an explicit all-indexed view for discovery and migration. Never delete prior records merely to simplify that view.
6. Keep logs, checkpoints, temporary preparation files, failed executions, and large trajectories in job/artifact storage. Store retrievable references and scientifically useful observations where appropriate.
7. Avoid creating empty molecular groups. Make automatic registration an explicit per-tool contract, not a blanket side effect of scanning every job folder.

An optional stricter preference could require promotion for ad hoc jobs whose inputs were never in the library. That would be a Studio product choice; Tamarind's current documentation instead says it creates missing input groups automatically. This distinction should stay explicit.

No ingestion or library-inclusion policy was changed during this research recheck.

## Cross-workstation status at the initial recheck

The local Studio is configured as peer `pc002975` with peers `pc004032` and `pc004033`, and the local agent is active. The added molecular adapter supports metadata replication plus authenticated structure retrieval from the owning workstation. At inspection time, the live replica had no molecular records and the local library had no replicated molecular groups.

Luna's read-only process and API checks established why configuration and liveness are insufficient here:

- Local `studio-agent.service` started on 19 October, before the molecular replication adapter and integration were added on 20 October. It keeps the previously imported replication module loaded and has no automatic code reload.
- Local `studio-api.service` also started before the new routes were added. Its live OpenAPI document contains the older peer endpoints but no molecular database or remote molecular structure route.
- Authenticated read-only peer queries succeeded on both `pc004032` and `pc004033`. Both exposed 253 job and 12 project records, with no molecular records. Their live OpenAPI documents also lacked the molecular database and structure retrieval routes.

At that point, the molecular API/replication rollout was not active across the three running workstations. No services were restarted and no remote code or records were changed during that diagnosis.

## Follow-up activation — 20 October 2026

The narrow database implementation was subsequently applied to both clean remote checkouts, and the API and peer agent services were reloaded on all three machines. Each peer now has the same 159 molecules, 347 groups, and 295 score observations; the replication catalog has 735 molecular records with zero conflicts. A remote PDB structure was fetched successfully from `pc002975` onto `pc004033` using the authenticated API and passed its checksum check. The full rollout and the existing scheduler code-revision mismatch are recorded in `molecular-database-implementation-handoff-20261020.md`.

Remote files remain on their source workstation until fetched. Metadata synchronization is not a full offline mirror of all structure files or trajectories.
