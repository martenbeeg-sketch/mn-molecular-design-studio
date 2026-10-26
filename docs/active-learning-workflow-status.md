# Active-learning docking workflow

## Bulk results and Database import

Campaign, docking-batch and model jobs persist a manual Database collection policy
before dispatch, including SSH workers. Intermediate compound inputs, scores and
poses are never automatically indexed as molecular records. Job progress,
original library references and output locations remain in Studio's catalog.

The full ranking stays in `artifacts/active_learning/docking_scores.parquet`
(original ID, name, SMILES, score, rank, rounds and worker provenance). Failures,
round history, models and retained poses also stay external. The Final results
node independently sets a Database shortlist of 0–1,000 final compounds; 0
disables import. Physical pose retention is a separate setting. Results pages
preview external data without importing it and offer bounded final promotion.

Standalone docking tiles expose **Automatically import compounds and results into
Database**, enabled for ordinary small jobs. Turn it off for bulk screens.
Restrictions survive job completion, refresh and forced reindexing.

`bulk_screening_maintenance.remove_automatic_bulk_records` removes accidental
past imports using explicit job references. Back up the Database and archive extra
portable poses first; this adapter never deletes files. Curated selections and
unrelated job relations survive. Publish replica tombstones before resuming sync.

## Current implementation

### Campaign results update — 2026-10-26

My results groups scheduled work across computers under its pipeline, using the
organizer's persistent node ID and parent run ID. Active-learning campaigns show
one row per iteration in the Steps dialog; opening an iteration reveals its
individual docking batches and model jobs. Counts reflect jobs created so far,
so the final iteration count does not have to be known at submission.

The **Final results** action opens the combined external ranking and convergence
plots. These show the best docking score, mean scores for the top 100, 1,000 and
10,000, and the largest rank-wise score change for those groups. A dashed line
shows the configured tolerance. The comparison table distinguishes score changes
from changes in compound order. Plotting additional groups does not change the
campaign's configured stopping criterion.

New runs save compact summaries after every round and one previous top-score
vector for comparison. Existing completed runs can reconstruct these summaries
from their external table only after checking that no score changed after its
first round, cumulative counts match, and saved top vectors match. The resulting
small `score_convergence.json` sidecar avoids scanning the compound table on
every page visit. No additional compounds enter the molecular Database.

The completed E91DA SpiroSpace campaign has 8 iterations, 160 docking batch jobs
and 7 model jobs across three computers. It stopped with three stable top-100
comparisons at iterations 6–8. The top-1,000 mean improved through iteration 7;
the top-10,000 mean still improved slightly in iteration 8.

### Controller design

The active-learning canvas now submits a scheduled campaign job through Studio's normal workflow queue. The controller runs the editable graph as a round-based workflow and writes node/worker progress to `graph_execution.json` for the pipeline dialog.

- The centroid library and full screening library are independent registered sources. A random, reproducible centroid sample carries its ID, name, SMILES, 1024-bit fingerprint and source-library ID into the first docking round. Fingerprint settings must match.
- The initial docking adapter is UniDock Pro in its normal `balance` mode. A round is shuffled and split evenly across selected, online Studio computers; every shard is queued through that computer's scheduler. Remote docking uses the Studio Conda environment and SSH user configured on the docking node. Score tables return to the organizer after each shard so cumulative ranking and model training can continue.
- Completed scores are accumulated by `(library_id, compound_key)`, retaining the best score for each compound. After each round, the cumulative list is sorted and its top-N compound identities and scores are compared with the previous round. A comparison is stable only when the same compounds remain in the same order and every rank's score change is within the configured tolerance. Consecutive stable rounds, screening-library exhaustion, or the optional iteration limit can stop the loop; `max_iterations=0` means no limit.
- Binder labels are recalculated from all accumulated docking scores before each training run using the configured percentage schedule. The controller streams fingerprints from the external screening dataset for prediction and selects the best predicted compounds plus random exploration compounds, excluding compounds already docked from that screening library.
- The final organizer-side output is a compressed Parquet score table, round history, model checkpoint and up to the configured number of best poses. Large working databases, full-library predictions and unselected poses are not retained. Per-ligand preparation files are pruned after the selected poses have been copied.
- The docking computer selector now defaults each selected machine to its available scheduler CPU slots, reserving one organizer slot for the campaign coordinator. Run details show the computers, CPU/GPU allocation and scheduler job link for docking and model steps.
- Other docking engines are not connected to this controller yet; the current active-learning graph only offers UniDock Pro so an unsupported choice cannot be submitted accidentally.

The OpenVS-style network uses PyTorch with three 3,000-node hidden layers. Studio environments use PyTorch 2.14.1 with CUDA 13.0 so the organizer can reuse its selected GPU for training and inference after each docking round. RDKit fingerprinting remains CPU-based. UniDock Pro continues using the GPUs selected for docking.

## Initial verification (historical)

The focused active-learning, pipeline, library-preparation, scheduler, and workstation test set passed: 94 passed, with one existing Starlette/httpx deprecation warning. The CUDA PyTorch install and a small CUDA operation still need verification on the three Studio machines. No multi-computer docking smoke run has been launched yet.

## Registered test data

The prepared library in `Spirochem library test random 1m seed 42` contains:

- 1,000,000 screening fingerprints: `/mnt/db/reference_files/SpiroChem/Spirochem library test random 1m seed 42/fingerprints/manifest.json`.
- 177,430 centroid representatives: `/mnt/db/reference_files/SpiroChem/Spirochem library test random 1m seed 42/centroid-preparation/centroids/manifest.json`.

These are manifests for externally stored files; the full compound table is not imported into Studio's molecule database.

## Initial smoke plan (historical; completed)

No docking campaign has been launched from this controller yet. The first smoke run should use a small initial sample and a low maximum iteration count. Confirm the receptor and docking box, which computers are available and their UniDock Pro GPU IDs, and the sample/iteration limits. The draft defaults retain the requested 10,000 compounds per docking round and the configured three-computer selection.

The smoke run should be observed in the Studio job table and pipeline dialog, confirming scheduler dispatch to each selected computer, returned scores, cumulative reranking, dynamic binder labels, model training and the stop decision. Do not use the 1,000,000-compound screening dataset for a full campaign until that small run has been reviewed.
