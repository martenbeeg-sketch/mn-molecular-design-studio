# Runs and result outputs

The scientific hierarchy is **Project → Run → named result outputs**. A Run can
contain native tool jobs, distributed batches, retries and later analyses.
Jobs remain the execution records used by the existing scheduler and ledgers.
Their directories, native parameters, filtering decisions and provenance stay
intact. Results groups jobs before pagination; individual jobs are available
under Execution details. Native singleton viewers and docking pipeline viewers
retain their specialized Results screens.

Analysis, acceptance filtering, sequence design and compatible cofolding inputs
select a complete Run output. One finite operation freezes its source list and
submits native collection jobs through the existing API/queue on each source
workstation. Processing produces one named output under the source Run. No
prediction or acceptance filter is added automatically. Whole-Run selection
does not apply an unqualified candidate ID or per-batch limit across sources.
An explicit saved subset should be used for that purpose.

The Run results graph input uses the same source list and the graph controller's
existing candidate-batch ledger. Compatible tool connections are checked against
their saved-collection signatures. Each native source collection stays intact;
the owning workstation is recorded in its immutable batch request. Subsequent
cohort ranking/merging follows an analysis or prediction step. No additional
scheduler or long-running controller is introduced.
Clustering retains its original whole-cohort protocol; it is never silently
split into independent native source collections.

Run identity derives from the project, owner and durable job/parent references.
Explicit Run links take precedence. Historical batch-name grouping retains the
previous project/owner/tool/date boundaries. Ambiguous or inaccessible parent
links are not guessed. Ownerless historical analyses are not silently assigned
to another user's project.
Historical distributed attempts can be linked by their exact campaign-ledger
job receipts. The small `run_membership.json` projection records the evidence
and matching owner/project, so Results need not repeatedly parse old ledgers.

Candidate IDs remain native. `run_candidate_id` qualifies them by their source
job; exports retain both IDs, workstation and native metrics. A parent's combined
collection is authoritative when it covers its children, avoiding double counts.
Results read normalized candidate metadata without re-running scientific scores.
Source structures and PAE matrices remain on their owning workstations. Requested
remote result metadata is loaded asynchronously; offline or failed sources are
reported explicitly. Complete exports are disabled while sources are unavailable.

## API

- `GET /api/v1/runs?project_id=…` lists logical Runs and named outputs.
- `GET /api/v1/runs/{reference}?project_id=…` includes execution/source jobs.
- `GET /api/v1/runs/{reference}/results?project_id=…` returns combined candidates
  and pending/error/excluded sources. `export=csv` or `export=jsonl` exports a
  complete output. A reference can include `~output_id`.
- Existing workflow validate/submit endpoints accept the output reference as
  `source_job`. Existing job references and explicit paths remain supported.
- `POST /api/v1/run-operations/{operation_id}/retry` resumes dispatch failures
  from the frozen operation manifest. Existing receipts and child idempotency
  keys prevent duplicate work. Native job retries remain with the existing queue.

Operation manifests live under the Studio workspace's `run_operations` directory.
They store the original source snapshot, request settings, dispatch receipts and
errors. Collection selection does not move scientific files or duplicate them
into the molecular database.
