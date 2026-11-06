# Remote Results Evaluation completion — 2026-11-06

Status: **implemented and verified locally; not deployed**. Production remains
`259cd17` on pc002975, PC004032 and PC004033. No production service was restarted,
no scientific calculation was repeated, and no scientific/curated Database row
was inserted, promoted, rewritten or resynchronized during this task.

## Confirmed failure and correction

The existing production calculation is valid:

- PC004033 source: `design/20261102-081701-97f7d42e`.
- Original candidate: `rfd4_000050_protein_binder_rep9`.
- Analysis: `analysis/20261106-180651-cee86ec8`, Biotite Interface.
- Source SHA-256: `bdbc6172e1d7ea156434791d0c1d515ddc1bd364ea94cab090aa9a71e6928438`.
- Chains B/A, 5 Å; 430 atom contacts, 63 residue pairs, 1 salt bridge.
- Hydrogen-bond measurements are unavailable because explicit hydrogens are required.

The Results Board used `_accessible_result_jobs()` and local `Path(run_dir)`
candidate scans. It consequently excluded peer-owned saved analyses even though
My results could read them. It also returned before rendering saved measurements
when the exploratory acceptance route selected no candidates.

The Board now uses the existing `logical_runs.catalog` and `run_results.collect`
transport for authorized source collections. Only selected structure artifacts
are retrieved, through the existing job-scoped artifact action. Their original
bytes are staged in private user inputs and checked against recorded digests.
The Board links analyses by producer/upstream node and job, candidate, exact
structure digest, model and chain roles. It preserves missing/pending/error
states separately from numeric zero. Saved measurements remain visible even
when a new analysis route has no accepted candidates.

Inspector → Results Evaluation Board and Board → Inspector/source Results retain
the exact remote source and selected candidate. Local canvas submission remains
available; peer-owned Board sources use the existing contextual Inspector
submission path. This task does not add remote canvas execution capabilities.
The remote Inspector uses the same saved-measurement lookup for its Biotite
panel, recognizes an exact completed result, and does not offer a duplicate
Biotite submission. Pending/denied history is shown as unavailable and disables
submission until the lookup is resolved.

An inspect receipt is now explicitly a cache, not an access grant: the source
job is reauthorized before every cached remote read. Revoked project access and
changed source identity cannot expose cached candidate measurements.

## Results and Database Evaluations

The scientific evidence is the persisted normalized candidate result. A
Results-associated Evaluation is a read-only descriptor of that saved evidence,
using the existing `evaluations.evaluation_descriptor` schema and identifier
helpers. Its `results-evaluation-*` ID is deterministic from the origin node,
producing run/candidate/method, recorded model, source structure, chain/ligand
selection and scientific parameters. Human-readable workstation aliases, local
cache paths and worker counts do not change its identity. Unknown model identity
stays unknown; a materialized filename is not a model ID.
The descriptor is derived on read from persistent Results records. The ID is not
written back into the original candidate, and no separately stored Evaluation
record is created. The same durable evidence reconstructs the same descriptor
and ID after restart or on another workstation.

A Database-associated Evaluation is the existing projection of a canonical
`molecule_observations` record associated with an exact existing molecule and
representation. Its existing ID is retained when recorded. A curated structure
representation requires the existing explicit collection/selection operation.
Reading Results creates neither representation nor observation nor molecule.

The production run has no Database observation on any of the three hosts. Its
Results source was not an exact collected Database subject, so the guarded
Database projection did not create one. A null raw `evaluation_id` therefore
does not imply missing scientific measurements. It is left unchanged.

The isolated saved-evidence Board check produces:

`results-evaluation-45998d32b97b255dfbce08eb259bcb08`

This identifies the Results-associated Evaluation; it is not a claimed Database
observation ID. No competing Evaluation database or additional observation log
was introduced. Repeated reads use the same saved record and same identity.

## Verification evidence

Evidence root (ignored, retained locally):
`.tmp/final-production-integration-20261106/`.

- `saved-biotite.json`: bounded read of the original source/analysis candidate,
  metadata and analysis manifest from PC004033.
- `source-structure.cif`: only the selected original CIF; SHA-256 verified above.
- `live-catalog-evidence.json`: the two existing producer-qualified catalog rows.
- `saved-board-evidence.json`: normal Results Board rendered by Streamlit
  AppTest using those exact saved candidates and structure bytes, with isolated
  replay of authorized metadata transport responses. It shows the original
  candidate, source job, parameters, structure and saved measurements. This is
  code/UI verification, not a claim that the undeployed production Board changed.
- Focused tests cover actual federated collection transport, stable identities,
  missing versus zero, pending results, model zero, wrong node/run rejection,
  external-path rejection, cached-read access revocation, actual coordinator
  project denial, Board rendering and Inspector/source navigation.
- Production API mode tests reject missing sessions and the development bearer.
  Existing node-token tests reject unauthenticated/incorrect tokens and publishing
  another node's jobs. All tests use isolated stores.
- Final verification: **35 passed, 64 deselected in 2.44 seconds**. Command:

  ```text
  python -m pytest -q tests/test_remote_results_evaluations.py tests/test_run_results.py tests/test_results_board_identity.py tests/test_central_scheduler_worker.py tests/test_complex_analysis.py tests/test_cluster_ui.py tests/test_cluster.py::test_api_enrollment_catalog_and_agent_commands_are_authenticated tests/test_cluster.py::test_node_token_cannot_publish_another_nodes_jobs -k 'remote or results_board or scheduler_bridge or local_preview or peer_results or model_zero or complex_analysis_initialization or biotite_receipt or production_api or api_enrollment_catalog or node_token_cannot_publish'
  ```

## Actual production authentication

Production is **http://127.0.0.1:8502/**. The previous 8510 URL was an isolated
preview address and is not the production endpoint.

The existing base UI/API service definitions explicitly set
`MN_MOLECULAR_DESIGN_STUDIO_AUTH_MODE=development`. The release drop-ins changed
source/PYTHONPATH only and inherited that flag. This is explicit configuration,
not an automatic authentication fallback. The old development-named service
configuration explains its provenance; no human intention to use it permanently
in production is assumed.

On all three hosts the actual account database contains **zero users, zero
enabled administrators and zero project memberships**. An unauthenticated
`GET /api/v1/auth/me` on the primary returns the synthetic `developer` administrator
(`00000000000000000000000000000001`). Consequently user/project ACL functions
run with admin privileges in this deployed mode; named-user authorization has
not been established. Existing enrolled-node tokens and peer permission checks
remain enforced independently. No identity, membership or enrollment was invented.

### Separately approved authentication change — prepared, not applied

Current unit definitions are backed up locally with mode 0600 under
`auth-change/{pc002975,pc004032,pc004033}-units.before.txt` in the evidence root.
The prepared `production-auth.conf` contains only the supported production-mode
environment override; it has not been installed in any production service.

1. Have the operator choose the legitimate initial administrator and enter their
   password using the existing interactive `python -m mn_molecular_design_studio.cli
   create-admin` command in the primary's actual Studio home/environment. This
   task does not create that account or supply credentials.
2. Allow existing trusted identity replication to deliver that same account;
   verify the actual ID, enabled status and required project memberships on all
   intended sign-in hosts. Do not create separate placeholder identities or
   rewrite historical owners. An administrator can initially inspect historical
   records under the existing admin policy.
3. Validate real sign-in and denial cases using an isolated authentication canary
   with the supported production setting, disabled cluster and local scheduler.
   The focused code test already confirms that production mode rejects absent
   and development-mode sessions. Actual credentials remain an operator gate.
4. In a separately approved quiet window, install the production override for
   both UI and API, preserve each host's home/node/project settings, and restart
   only the relevant Studio services. Confirm the effective environment and real
   `/auth/me` account on each workstation, then check project-member/outsider access.
5. Retain unit backups and the prior code release. On failed sign-in, keep the
   interface protected while correcting configuration; never silently fall back
   to development admin access. Rollback that would restore development access
   requires an explicit operator decision and loopback-only exposure.

## Preview isolation — applied only to the isolated preview

The old preview API PID 3410265 had no active preview jobs, scheduler claims or
scientific children. Its cluster config was still `peer` despite `accept_jobs:
false`; that flag alone does not stop central bridge registration.

Only that API was gracefully stopped and restarted as PID 4117915. It uses its
original preview home/Results, `SCHEDULER_MODE=local`, an isolated cluster home
with `mode: disabled`, no peers, and an isolated local scheduler state directory.
The original preview configuration is backed up with mode 0600. Preview files,
saved Results and code were preserved. Both central startup paths now refuse
registration in local mode, including direct bridge CLI invocation. The saved
preview launcher retains these settings for subsequent restarts.

Production API/UI/agent PIDs remained 4005534/4005536/4005535. No production
worker or MN Compute Scheduler service was stopped or restarted. The old preview
bridge's six scheduler registrations are offline; they have not been deleted.

## Bounded production integrity snapshot

All SQLite checks used `mode=ro`, without app schema initializers, forced
replication, migrations, reimports, promotion, or cleanup.

| Host | Molecular identities | All indexed membership identities | Library identities | Groups | Observations | Replica conflicts |
|---|---:|---:|---:|---:|---:|---:|
| pc002975 | 21,473 | 20,560 | 2,341 | 1,050 | 41,198 | 1 |
| PC004032 | 21,470 | 20,557 | 2,341 | 1,071 | 41,070 | 37 |
| PC004033 | 21,470 | 20,557 | 2,341 | 1,071 | 41,070 | 37 |

Primary Database:
`/home/user/.local/share/mn-molecular-design-studio/database/molecular_library.sqlite3`.
Both remotes:
`/home/compute/.local/share/mn-molecular-design-studio/database/molecular_library.sqlite3`.
Primary Results are under `/mnt/data/RESULTS/mn-molecular-design-studio/workdir/runs`;
remote Results are under `/data/RESULTS/mn-molecular-design-studio/workdir/runs`.

Library uses the existing curated/generated/reference membership policy and
existing shortlist exclusion. All indexed additionally exposes other indexed
memberships. Ungrouped identities are included in the total identity count, not
the membership count. Historical PD-L1 groups and molecular memberships remain
present. No mass loss is indicated by these counts.

The unfiltered production browser shows **1,771 Runs**, All projects, All
workstations, blank search and no status filter. Primary local catalog: 2,840
jobs; its shared catalog: 4,388 jobs. Peer local catalogs: 884 (4032), 830 (4033);
their shared catalogs: 5,050 and 4,428 respectively. Jobs and grouped Runs are
different counts. The earlier eight-run view was project-filtered.

Existing replication cursors cover the observed peer sequence tips, but conflict
rows remain and the primary has three extra identities / 128 extra observations.
The historical Library identity count agrees. Full replica equality is **not**
claimed; those bounded discrepancies are preserved for the separate Database
task rather than being repaired here.

## Minimal later release and acceptance

1. Approve the focused local commit and the existing coordinated Studio release
   procedure. No dependency, transport protocol or Database schema change is
   required. Keep each host's paths, node enrollment and data intact; take the
   established quiet-window code/config rollback point. Publish no Git changes
   automatically. Keep authentication approval/configuration separate as above.
2. Install the pinned code through the existing release mechanism on the Studio
   installations that must display the Board. The old remote inspect/artifact
   responses are compatible. Use the same verified revision/fingerprint across
   the participating Studio processes before subsequent scientific submissions.
3. On primary 8502, open Complex Analysis → Results Evaluation Board. Select
   PC004033 `design/20261102-081701-97f7d42e` and the exact candidate above. Use
   Refresh saved Results while metadata/artifact requests are pending.
4. Verify the existing Biotite row and Results Evaluation ID above, 5 Å, B/A,
   430/63/1, and unavailable H-bonds. Open in Scientific Inspector and return to
   source Results; the same candidate and SHA must remain selected. Do not press
   Run analysis: the completed calculation is reused.
5. Under real sign-in, confirm an authorized project member can inspect it and
   an unrelated user cannot read its job, cached receipt or artifact. Confirm
   the saved Evaluation remains available after a UI/API restart in that approved
   window. Verify no new curated membership or scientific observation was created.

Remaining gates: deploy this code separately; establish real production accounts
and sign-in separately; perform the saved-run acceptance in that environment.
Remote canvas execution from this Board and replica-conflict cleanup are not
expanded by this task.
