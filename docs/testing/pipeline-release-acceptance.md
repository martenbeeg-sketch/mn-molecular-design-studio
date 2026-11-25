# Studio pipeline release acceptance

Run the quick pass after a pipeline fix and the full campaign before a release or worker image update. Use the **Studio Release Acceptance** project and fixtures reserved for software testing. Start with one coordinator and schedule work across the available GPUs through Studio. Allow about six hours for the full pass and record any case still queued or running at the deadline.

## Start a run

1. Use Studio MCP to identify the signed-in account, release project and workstations. Record the Studio revision, each worker's available resources, image digest and compatibility, and the fixture checksums. Set limits for queue occupancy, runtime and test artifact storage.
2. Snapshot the Capabilities inventory. Key each case by capability, engine, mode, target context and output format. Link each case to its input fixture, expected output and required nodes. Carry forward the unresolved cases from the previous run and add rows for new modes.
3. For the quick pass, run the short pipelines and inspect their Results. Start from the saved **Release smoke · synthetic sequence to structure and sequence QC** pipeline in the release project. Run at least one case with a regular user in Streamlit and one through MCP.
4. For the full pass, continue through the windows below. Validate and save the graphs before launching them. Submit each case with a stable idempotency key built from the release, case ID and fixture revision. Follow the returned run IDs to a terminal state or the test budget.
5. For every terminal run, inspect the output and attempt the next compatible handoff. Open the corresponding UI Results. Compare the expected project records with the actual database records. Record the exact observed result before moving a case to passed.
6. End with a coverage summary by capability. Include counts for live executions, UI verification, passing controls, fresh passes, blocked cases and unfinished cases. Use the open case IDs on the next run.

### MCP checkpoints

Use `validate_pipeline` to confirm readiness, then `save_pipeline_draft` and `run_pipeline` to record the graph revision and run ID. Check `list_jobs` for bounded execution status and `list_runs` for the consolidated Results entry. Read a small output with `run_results`. If it reports that inspection is required, use `studio_api_request` to send `POST` to the returned `inspection_url` path with the project ID as a query parameter. Poll any returned command IDs with `command_status`, then reread Results. For an individual worker job, use `request_job_inspection` and read its completed command response. Record execution status and output availability separately. Use a bounded `result_presentation` query to check that a representative row can be rendered. A large run can exceed the MCP response limit for `get_job`, so use paged Results for candidate checks.

For a short app regression check before using the GPU queue, run the focused suites from the checkout.

```bash
python -m pytest -q tests/test_pipeline_graph.py tests/test_pipeline_dataflow_contracts.py tests/test_candidate_execution.py tests/test_candidate_sets.py tests/test_artifact_collection.py tests/test_candidate_peer_transfer.py
python -m pytest -q tests/test_mcp_server.py tests/test_api_integrations.py tests/test_run_results.py tests/test_remote_results_evaluations.py tests/test_result_presentation_identity.py tests/test_molecular_database_results_ui.py
python -m pytest -q tests/test_quality_tools.py
```

On a test workstation with the QC images installed, run `MDS_QC_DOCKER_TESTS=1 python -m pytest -q tests/test_quality_tools.py` to exercise the native QC containers.

### Fixture rules

- Use fixed small-molecule decks, inert synthetic sequences and archived nonfunctional structures with checksums. Use fixture artifacts and contract tests for capabilities that lack a suitable live test input. Record those cases as contract coverage.
- Include a valid handoff, a missing or incompatible artifact, and a fixture with multiple chain or model IDs. Give accepted, rejected and incomplete fixtures distinct stable IDs. Test archive reads separately from newly generated test outputs.
- Put any intentionally saved records in the release project. Record their IDs so the next campaign can distinguish old and new assertions. Keep test library artifacts in a test-scoped catalog, or use temporary fixtures for the preparation tests if the deployment only offers a shared catalog.

## Coverage and execution

Build the case inventory from the live Capabilities tables and MCP `tool_connections`, `workflow_setup` and `list_workflows`. The previous [tool execution inventory](tool-execution-coverage.md) supplies starting fixtures and identifies engines that needed a patched replay. Use the live catalogue to catch newly added modes.

### Capabilities coverage and student entry point

1. Make a row for every actionable table cell and every workflow variant it exposes. Record the capability, task, target input, format, engine or checkpoint, and workflow variant together. For protein design, a miniprotein pass leaves linear peptide, macrocyclic peptide, VHH and any other advertised formats as separate cases for that engine.
2. As a regular user, open the pipeline builder action from the cell. Inspect the setup options, open the canvas, and verify that the expected nodes, typed connections and selected settings appear. Save and reopen the pipeline. Compare the saved graph with the equivalent MCP preset. Record browser setup, graph validation and job execution separately.
3. Check every advertised producer-to-consumer connection with fixture artifacts. Cover each alternative shown in the setup dialog at least once, including optional steps on and off. Test representative multi-step combinations, including a valid connection with missing input and a deliberately incompatible connection. Record the exact composition in the evidence row so coverage of one combination stays specific to that combination.
4. For a case with a suitable live fixture, submit the graph built from the table and follow it to Results. Record the run ID, observed output and first downstream handoff. Use archived artifacts for graph and handoff checks where live execution is unavailable. Keep catalogue-only, contract-tested, browser-verified, live-passed, failed, blocked and untested cases visible in the coverage report.
5. Before sharing a preset with students, check the evidence for that exact preset and its advertised variants. Record open cases at the six-hour cutoff so the next session starts with a known list.

Run the preset contract suite from the release checkout. Its named cases form the baseline checklist for every table route, refolding choice and supported starting input. Save the JUnit output with the acceptance evidence. A contract pass records the assembled graph and its serialized round trip. Record the actual browser steps and live run IDs separately.

```bash
python -m pytest -q tests/test_capabilities_acceptance.py --junitxml=capabilities-contracts.xml
```

For an isolated browser check, launch both the UI and API from the same checkout. Use a fresh local workspace with an explicit Results directory and a regular fixture account. The default Results directory can be shared even when the workspace directory is temporary. Use a dedicated scheduler or local execution for this fixture so production workers cannot claim the test jobs. For a GPU smoke test with local execution, set an explicit GPU index. Repeat the distributed queue checks in the approved Studio acceptance project. Verify that the browser and workers load the intended revision before comparing their results.

Count visible table buttons separately from the workflow variants inside each dialog. Sweep the buttons as a regular user and check that the setup can be opened by keyboard and pointer after horizontal scrolling. Use the widget tests to verify each binder variant through draft save and reopen. Include two native handoffs in the browser pass. For library preparation, connect fingerprints to centroids and confirm the rejected row and shared library ID. For database input, confirm that a frozen selection reaches Results with its source molecule reference and that a snapshot from another project is rejected.

| Window | Pipeline cases | Required evidence |
| --- | --- | --- |
| First 30 minutes | Short editable pipelines for sequence input, structure prediction, sequence QC, library preparation, and a pose analysis. Open one from a Capability action, save and reopen it in Streamlit, then submit an equivalent draft through MCP. | Required inputs are visible and incompatible connections are rejected. The same settings reach the same workflow, Results opens, and expected structures, tables, downloads or manifests are available. |
| Hours 1–3 | Protein workflow contracts by capability. Cover structured, sequence-only, IDR and multidomain input, plus each advertised binder format. Include new folds, fixed-backbone redesign, motif scaffolding and nucleic-acid-conditioned nodes in the connection and fixture checks. Run end-to-end cases on the designated inert synthetic fixtures. | Verify that the modular generator, sequence-design and refolding nodes route their expected artifacts. Sequence-free backbones select the sequence-design handoff. Co-designed sequences retain their identity. Verify connected-template and server/cache MSA provenance where applicable. Verify filter receipts, any eligible database save and the BinderScope link from the saved row. |
| Hours 2–4 | Small-molecule generation, docking, affinity prediction, pose analysis and a short multi-round library screen. Use small fixed decks and a bounded shortlist. Exercise the import, prepare, screen and inspect paths. | Compare input and output IDs, validate chemical records, open poses and scores, verify unique shortlist retention across rounds, and check external library storage. Small-molecule design acceptance is undefined. Report output validity and any selected analysis metrics. |
| Hours 3–5 | Structure preparation, short OpenMM and GROMACS simulations, conformational sampling, saved run reuse and native ligand analysis on trajectory frames. | Keep topology and trajectory paired. Verify frame counts, preparation provenance, model or pose identity and downstream Results. |
| Final hour | Run a long graph with branching and fan-in. Include a short final batch, an empty accepted set, a reloaded upstream run and a recoverable worker failure on test jobs. Rerun at least one short graph from Streamlit as a regular user. Exercise project draft editing and shared-template permissions using the test fixtures. | One Results row per pipeline run, linked worker batches, accurate partial results, preserved provenance, reuse of completed work after a retry, correct project and template access, and exactly one eligible database record per candidate. |

Submit independent cases together where resources permit. Keep each case's candidate count, batch size, retry limit and runtime budget in the request. Give the quick check one compatible GPU. In the distributed check, compare actual candidate counts and native process launches with the planned batch sizes and recheck the same counts after a retry. Use fresh executions for acceptance evidence. Reuse a verified result explicitly when testing downstream handoffs or recovery. Keep expensive engine cases separate so a blocked engine does not delay unrelated checks.

## Protein acceptance

Record two independent outcomes for each protein test case. The **workflow** outcome records whether the expected artifacts reached Results. The **filter** outcome records how many new samples passed the configured filter on the designated inert fixture. Use a bounded attempt count and retain rejection reasons. A completed run with zero passes is an execution success if its expected artifacts are complete. Use a known passing control to exercise the exact-structure save path, and rejected and incomplete controls to exercise exclusion. Count controls separately from new samples.

Run the same configured acceptance rules on both the control and the fresh case. Keep the rule set and schema version in the run. Check each accepted row against the exact structure it saves in the test project. Review rejected and incomplete rows through Results. Small-molecule cases can validate chemical records and pose geometry while acceptance criteria are defined for individual studies.

## Evidence per case

Record case ID, capability, workflow and mode, pipeline ID and revision, project, fixture, expected output and acceptance policy, worker job IDs, duration, native invocation count, terminal status, generated and scored counts, fresh passes, controls, Results links and observed database writes. Append the release code revision and test conclusion to the evidence record. Use **passed**, **failed**, **blocked**, and **unfinished** for execution. Track **fresh pass**, **zero fresh passes**, and **acceptance undefined** separately from execution. Label coverage as **catalogue**, **contract**, **live**, or **live and UI**. A live case needs an inspected output and downstream handoff.

For the same saved run, record the time to first usable Results table, the wait to open a worker step, and whether scrolling stays responsive during table and structure browsing. Capture the workload size and distinguish the first load from a repeat load. Compare with the previous release on the same fixture and workstation.

Use MCP to validate, save, run and inspect each case. Use the Streamlit run details to verify representative Results for each family and any reported UI failure. Retest a failed case after a fix on the same fixture before changing its status. Keep the failed run ID in the evidence.

### Run record

Copy one row per case into the release notes. Record the review date, Studio revision, fixture revision, and tested workstations above the table.

| Case | Fixture and graph revision | Run and worker references | Expected → observed output | Coverage | Execution | Filter outcome | Database IDs and Results link |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `SMOKE-01` Sequence to structure and QC | Synthetic sequence · saved smoke pipeline | Record run and worker IDs | Predicted structure and a sequence composition row → record actual files and rows | Live and UI after checking both | Pending | Acceptance undefined | Record actual writes and link |

### First smoke execution

On 25 November 2026, Studio MCP submitted the saved `SMOKE-01` pipeline. Runs `CBED6` and `6441F` exposed problems carrying standalone prediction records and their recorded sequences to Sequence QC. Replay `589C2` reused the saved ESMFold prediction and published one structure and one sequence composition row for the 20-residue fixture. The regular-user Streamlit follow-up `9F0D8` submitted an equivalent ESMFold2 Full pipeline from Capabilities and verified the Results view. [Smoke and browser receipts](evidence/capabilities-acceptance-2026-11-25.md)

Distributed replay `F92AC` in the BBF-14 benchmark project collected 106 refolds across three workstations. All 106 had ipSAE measurements and filter decisions, including 11 passes. This replay covered matching worker images, PAE transfers from the original refold jobs and paged MCP Results. Parent-level cancellation of the earlier mixed-workstation run `05688` returned an artifact-location error. Its queued worker jobs were cancelled individually. Keep mixed-workstation parent cancellation in the next replay.

The 26 November acceptance pass covered all three TMol modes, the Ligand SDF input, all three BioEmu-2 modes and 11 native QC containers. It verified saved-sequence reuse through TMol and stable Results links during job submission after two fixes. [Run and worker receipts](evidence/capabilities-native-2026-11-26.md)

## Open acceptance cases

These cases carry forward the gaps in the [distributed run](evidence/release-acceptance-2026-11-25.md), the [regular-user browser sweep](evidence/capabilities-follow-up-2026-11-25.md) and the [native mode replay](evidence/capabilities-native-2026-11-26.md). A case closes when the required evidence below is linked to the tested release and fixture.

| Case | Coverage to date | Evidence needed |
| --- | --- | --- |
| `LIB-02` Multi-round screening | The browser run covered fingerprint preparation and centroids. Retention and round selection have regression coverage. | Run at least two docking rounds on the same small deck across two workstations. Check first-round centroids, later-round selection and exclusion of previously docked IDs. Compare the final unique shortlist and retained poses against an independent sort of all round scores. Verify that cleanup preserves the shortlist and failure records. |
| `RUN-02` Distributed retry and cancellation | One mixed-workstation cancellation and successful fan-in were recorded in `921A0` and `FB731`. | Exercise retry and cancellation with workers active on two hosts at once. Verify that completed work is reused, prerequisite failures block dependents, every cancelled worker stops, and the parent has one durable Results row. |
| `UI-02` Non-binder browser presets | All visible setup buttons were opened as a regular user. All 82 binder presets were exercised in Streamlit widget tests through save and reopen. | Save and reopen the remaining non-binder combinations and compare the settings and typed edges with their API equivalents. Use the per-case checklist to record each browser result. |
| `UI-03` Pipeline and template permissions | Project and template lifecycle have regression coverage. | Verify distinct regular-user and instructor accounts. The regular user can open and copy a shared template and save a project pipeline. The instructor can update the shared template, and existing run snapshots still open at their recorded version. |
| `DB-02` Admission controls | Browser and API checks exercised frozen database input and access checks. | Run accepted, rejected and incomplete archived controls through admission on an isolated database. Verify the exact saved structure and scores, one publication across a retry, and analysis attachment to the accepted record. |
| `FMT-02` Generator, format and target-context combinations | Preset contracts and browser setup cover the advertised formats. Native evidence covers specific recorded fixtures. | Record one execution and downstream handoff per feasible format, engine, and target context using the designated benign fixture. Record contract-only and unavailable cases explicitly. Verify the optional refolder, sequence design, and analysis settings through the graph they create. |
| `AN-02` BioEmu propagation from Results | Propagation ran successfully using the coordinate file from a saved prediction. | Check the completed-design input route end to end. Direct use of a named Run output requires support for an ensemble-producing consumer in the saved-run contract. |
| `PERF-01` Results response time | Worker links and Results rendering were checked on small fixtures. | Measure first usable table, worker drill-down and scroll response on a fixed 500-row result and a larger archived run. Record first and repeat loads on the same workstation and compare with the prior release. |
