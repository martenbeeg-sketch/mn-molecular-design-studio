# MN Molecular Design Studio

For authenticated REST access, external library/results queries, file exports
and the optional MCP adapter, see [REST and MCP integration](docs/api-and-mcp.md).

Molecular Design Studio is a separate project with its own code, workspace,
job catalog, Results folder, database, and organization data. Workflow modules
needed by Studio are kept inside this repository; the previous app is a
reference only and its job results are not imported.

The [LocalFold integration](docs/localfold.md) adds the current sokrypton JAX
fork with 19 model choices, including Mol-1 and BioEmu2. Its model files and
shared caches live under the configured `reference_files/localfold` directory.

The Results folder starts empty, defaults to the same data volume as the
previous app, and can be changed in Settings. The Database includes only
molecule sets imported or created through Molecular Design Studio.

### Browse connected Database results

The Database keeps its three-pane layout: choose and filter on the left,
select a connected result in the center, and inspect its structure and scores
on the right.

- **Molecules** defaults to **Targets & ligands**. A recorded binding target or
  receptor is distinct from a designed binder; **Designed molecules** and
  **All molecular entities** are separate choices. Recorded target labels take
  precedence over carrier filenames or a complex's binder alias. **Whole molecule**
  keeps complete complexes, partners and stoichiometry separate.
- **Runs** shows the individual results of a recorded execution. Campaigns
  combine only explicitly linked executions, across their source collections.
- **Collections** are imported libraries or folders of selected molecules,
  not new calculations. Existing selection, editing and export remain here.

Each evaluation/model retains its own scores and representations. A structure
without a recorded evaluation is shown without borrowing another job's scores.
Repeated target carriers and other input-only history are hidden from outputs
unless **Include input / reference history** is enabled. Structures/evaluations
appear first; this order is not a scientific score ranking.

Links use accessible exact contents, explicit target IDs, or an unambiguous
origin-qualified run's recorded target input. Ambiguous or conflicting legacy
relationships are not guessed from names. Session-local metadata caching checks
fresh collection access and all SQLite commits, including peer edits and deletion;
the first load after new data can still take several seconds. Browsing does not
merge, rename or delete original records.

### Experimental Data browser

Project and target/ligand names now have [canonical SQL labels](docs/canonical-entity-names.md),
contextual Properties dialogs, searchable previous names and stale-edit protection.
Scientific IDs, structures and filesystem paths are unchanged. Shared label
editing stays on its source workstation; all peers need the naming-capable
release for consistent target labels.

The separate **Data browser** page adds a folder tree without replacing Database:

- **Projects** is the default: project → target or compound → activity → experiment.
  Single-activity levels are omitted. Imported libraries open into their compounds;
  reusable reference libraries expose their targets directly, not as storage folders.
  **Start at** also offers targets, compounds and deliberate libraries/datasets as
  shortcuts to the same SQL records. Automatic per-job registration batches are
  not datasets in the navigation tree.
- Select an activity or experiment to list its saved outcomes or individual
  structures. Search and optional refinements apply to the selected level;
  pagination covers every matching row. Project paths retain project context;
  a molecule opened from a collection can show its other accessible calculations.
- Experiments group only recorded logical/campaign lineage within the same owner
  and project, never matching names alone. **Source jobs & provenance** retains
  execution batches and workstation identities. Results without a verified subject
  remain accessible under **Unassigned results** rather than being assigned a target.
- Docking exposes the recorded engine (for example Vina, GNINA or Uni-Dock Pro),
  with separate experiments underneath when several exist. Engine selections keep
  only that engine's evaluations and poses, even when a saved outcome spans tools.
- Unambiguous exact sequence crops can be browsed under their registered reference
  family in the same owner/project context. This is navigation only: exact target
  IDs, constructs and evaluations remain separate; **Target construct** refines the
  family listing. Ambiguous containment and different PDB contexts are not folded.
- The right pane keeps original measurements, parameters, provenance and attached
  post-analysis together. Design-internal refolding is not another design row.
  Missing analysis is reported as missing, never as a zero or an inferred pass.
  Protein–ligand complexes use protein cartoons without side-chain sticks and keep
  native mmCIF atom names, coordinates and connectivity. Missing legacy pose links
  are distinguished from score-only results. Recovery requires a unique exact file
  checksum with matching source context; conflicting references are never guessed.
  Historical MN-Ligand docking links can also be reconciled with
  `scripts/reconcile_legacy_docking.py`: `plan` compares the original attempt,
  engine, compound, replicate, numeric scores, pose path and exact receptor frame;
  `apply` requires that unchanged plan, an administrator and a new SQLite recovery
  backup. It stores portable copies of the verified coordinates and an audit of
  the replaced links, without importing molecules, modifying scores or deleting
  historical representations. Apply on the owning workstation, not a peer copy;
  normal replication distributes the corrected source records. Source conflicts
  and concurrent SQL edits fail closed. Rehearse against a separate database first.
  Molecular-entry replication snapshots allow up to 8 MiB for accumulated
  scientific evidence; other metadata remains limited to 2 MiB. Synchronization
  pages are byte-bounded to 16 MiB and retain their continuation cursor.
- **Molecular dynamics** opens a study → replica tree, with saved RMSD, RMSF,
  gyration, contact and thermodynamic plots in the middle pane and the exact
  source pose / recorded topology snapshot on the right. Endpoint energies and
  trajectory analyses remain separate, source-qualified calculations. See
  [MD browser semantics and limits](docs/md-data-browser.md).
- Recorded PDB references have separate target folders, with exact protein
  constructs underneath when several exist. An identical sequence does not
  merge different PDB studies or their receptor preparations. Method folders
  select that method's own saved evaluations and coordinates. Nesso affinity-only
  calculations appear under **Binding affinity**, not as invented structures.
  Calculations without an explicit PDB reference remain accessible under a
  **PDB unspecified** target context instead of being guessed into a PDB study.
- Native result ingestion and later Results selections use the same exact-target
  registration contract. Target accessions come from explicit input metadata,
  an accessible selected/byte-identical SQL representation, or the coordinate
  file's accession field, never from a display name or a sequence-only guess.
  `studio.target-reference.v1` preserves recorded, unspecified and conflicting
  references separately from canonical molecule identity. Conflicts withhold the
  PDB assignment without merging studies. Missing reference metadata is not a pass.
  Original structures, scores and attached analysis carry the same reference;
  post-analysis inherits the exact selected source rather than its output name.
  Registered artifact-contract workflow IDs also supply purpose when a native
  adapter does not repeat its workflow ID in the job's input or metadata.
  Successful coordinate-free affinity calculations are recorded calculations,
  not passed designs. Unreviewed screening inputs remain outside saved outcomes.
  Database-bound OpenMM/GROMACS submissions retain owner-qualified pose lineage;
  ingestion checks the saved pose, receptor frame and staged input-complex digest
  before attaching MD measurements to the existing ligand. Changed/missing source
  inputs are withheld, not promoted as another molecule. These contracts apply to
  new registrations; existing data is not bulk-reingested by opening the browser.
- **Folder properties** and **Result / molecule properties** open shared, typed
  editors for project, dataset, compound/target, saved outcome and explicit dataset
  membership annotations. The left tree stays navigation-only. Dialogs identify
  the affected scope; Save is explicit and Cancel/Close discard drafts. Project
  properties use the same editor on the Projects page. Scientific fields remain
  read-only; project moves, cleanup and access management stay separate actions.
- The page opens SQLite read-only, including when derived Evaluation indexes are
  missing. It never indexes workdir jobs, runs an analysis or admits candidates
  merely because somebody browses them.
  Explicit metadata saves use a separate writable transaction, not this reader.

Automatically named designs use `Design_<generating-computer>_<creation-time-UTC>`.
The time comes from the individual native output, never the run start, SQL
registration, import or transfer time. RFD4 retains its embedded CIF creation
stamp; original native structure-write times are explicitly distinguished when
an engine did not record a separate sequence-generation timestamp. Missing
creation evidence is not replaced by a fabricated time. Short stable collision
suffixes are used only when needed. There is no distributed counter or
score-dependent renumbering. Target, engine, run, replicate and native names remain in
`metadata.studio_naming` and searchable aliases; molecular fingerprints,
candidate IDs and Evaluation IDs do not depend on display names. Meaningful
user-supplied names are retained. Unresolved target names are not guessed.

Target folders in both browser pages use the shared SQL-backed component name
resolver. It prefers a descriptive registered target name, then an accessible
named input/reference with the exact same molecular contents. Designed-binder
names and whole-complex/compound names are not used as protein-target names.
Existing aliases stay searchable. If SQL has no descriptive name, the browser
retains an explicit recorded label or shows a type/length fallback rather than
guessing a biological identity. No target-specific naming rules or online lookup
are needed while browsing.

`scripts/curate_library_names.py` provides a separate administrator-only
plan/apply workflow for historical naming and exact-target registration. It
edits only source-owned collections, requires a new integrity-checked SQLite
backup, rejects stale reviewed rows atomically, and preserves original target
references alongside canonical links. Optional reviewed reference constructs
can name an exact target or uniquely matching crop; names alone never link data.
Apply on owning workstations and synchronize through normal metadata replication.
`scripts/reconcile_legacy_targets.py` separately plans and backs up missing
MN-Ligand target links. Docking uses the actual copied receptor, not a reused
artifact ID or an unverified reference checksum. Cofolding/affinity uses the
literal protein chains in the exact compound's preserved model input. The plan
retains source-file hashes and original target references, rejects conflicts,
checks source snapshots again before applying, and never imports more candidates
or changes scores, coordinates or admission decisions. Browsing does not run this
maintenance operation. Historical repairs must run on their owning workstation;
mixed-version readers need the new browser code to expose all activity folders.
See the [real-data verification and limitations](docs/audits/2026-11-14-readable-library/verification.md)
for the initial preview and historical naming maintenance.
The [creation-time follow-up](docs/design-creation-naming.md) explains source
evidence, retention during downstream processing, and the backed-up recovery tool.

Studio also supports an optional [workstation cluster](docs/workstation-cluster.md):
shared job visibility, machine availability, remote submission through local
schedulers, and analysis on the workstation that owns its input data. Enable it
in Settings → Workstations or with the `cluster` CLI commands. It is disabled by
default and does not move existing results or scheduler jobs.

## Create and activate the environment

From this directory:

```bash
conda env create -f environment.yml
conda activate mn-molecular-design-studio
```

The environment installs Molecular Design Studio, the sibling
`mn-compute-scheduler` client package, and CUDA-enabled PyTorch for active-
learning training and inference on NVIDIA GPUs. The same PyTorch build also
supports `pydssp`. It does not install or edit the previous protein design app.
The scheduler API is a separate service and must already be running before
submitted Studio jobs can be claimed.

The scheduler checkout is expected next to this repository at
`../mn-compute-scheduler`. For a pip-only setup, install the optional scheduler
extra with `pip install -e '.[scheduler]'` instead. The scheduler client reads
`MN_COMPUTE_SCHEDULER_API_URL` (default `http://127.0.0.1:8765`) and the optional
`MN_COMPUTE_SCHEDULER_API_TOKEN` from the environment. Set
`MN_MOLECULAR_DESIGN_STUDIO_SCHEDULER_MODE=local` only to explicitly use the
Studio's local worker queue.

## Start the local API

```bash
mn-molecular-design-studio api
```

Authentication defaults to production mode. For local development and UI
testing, both the API and Streamlit process can be run with authentication
disabled and a synthetic local admin account. Set this in each terminal before
starting that process:

```bash
export MN_MOLECULAR_DESIGN_STUDIO_AUTH_MODE=development
mn-molecular-design-studio api
```

In a second terminal, export the same setting before running Streamlit. The
development identity is not stored as a real user, and account creation and
password management are disabled in this mode. Keep development mode bound to
localhost. Unset the variable or set it to `production` before testing real
user accounts; changing modes requires restarting the API and Streamlit.

```bash
export MN_MOLECULAR_DESIGN_STUDIO_AUTH_MODE=development
streamlit run mn_molecular_design_studio/run_app.py
```

The API binds to `127.0.0.1:8766`. Interactive API documentation is available
at `http://127.0.0.1:8766/docs` while the server is running. Admin API requests
can use absolute paths or managed `app:///`, `reference:///`, and `runs:///`
paths; relative paths resolve from the server's working directory, or an admin
may pass `?base_dir=/path/to/inputs`. Regular users are limited to their private
uploads, shared reference files, and files inside jobs they can access; they
cannot choose a server-side `base_dir`.

Initial endpoints cover workflow discovery and schemas, individual workflow
and campaign validation/submission, and job listing, inspection, result
retrieval, artifact listing, and cancellation. `GET /api/v1/jobs` accepts
search, status, workflow, and task-group filters plus `limit`, `offset`, or a
cursor, and returns `items`, the filtered `total`, and `next_cursor`. Job facets
are available from `/api/v1/jobs/facets`.

## Start the Streamlit catalog

```bash
streamlit run mn_molecular_design_studio/run_app.py
```

The Streamlit app includes the workflow catalog, job results, molecule
database/import, and multi-engine pipelines.

## Run tests

Install the test extra in the Studio environment, then run the full offline
suite:

```bash
pip install -e '.[test]'
python -m pytest
```

The tests cover Boltz and BoltzGen input contracts, BindCraft and BindCraft 2
settings and target validation, workflow validation, job metadata, and the
scheduler submission boundary using a fake scheduler client. They do not start
containers, use GPUs, or submit jobs to the live scheduler.

## Accounts and project access

In production mode, create the first administrator account once, from the
computer running the Studio:

```bash
mn-molecular-design-studio create-admin
```

The command prompts for the admin username and password. Once signed in, use
**Users & access** to create user accounts with temporary passwords and grant
access to selected projects. New users must change their temporary password at
first sign-in. Passwords are stored as scrypt hashes; session tokens are stored
as hashes and expire after 12 hours.

If you forget the admin password, recover it from the Studio host with
`mn-molecular-design-studio reset-admin-password`.

Users can see their own jobs and jobs in projects they have been granted.
Admins can see all jobs and projects, including existing runs that predate
accounts. Job listing, result details, artifact downloads, cancellation, and
the workflow API all check access on the server. The API is still bound to
`127.0.0.1`; API clients sign in with `POST /api/v1/auth/login` and send the
returned bearer token on later requests.

By default, Molecular Design Studio stores its data in
`~/.local/share/mn-molecular-design-studio`. Set
`MN_MOLECULAR_DESIGN_STUDIO_HOME` before starting the API and Streamlit app to
choose another workspace. Both processes must use the same value. Configure
the job results directory from the app's Settings page; both processes pick up
the saved directory without a restart.

The Settings page stores the Results directory as `results_directory` in the
workspace's `settings.json`. Studio loads that value into
`MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR` for API requests and scheduler workers. Do
not set that variable to a different path in only one terminal when the app is
using its saved Settings value. See [the scheduler worker runbook](docs/scheduler-worker.md)
for path checks and safe recovery steps.

Job records are indexed in SQLite at
`<workspace>/database/job_catalog.sqlite3`; the catalog stores searchable
metadata and relative links to the files under the configured results folder.
Inputs, results, logs, structures, and other run artifacts remain in their
per-job folders. The catalog is rebuilt from those folders when needed, so it
does not become the only copy of job data.

Projects are created in the Studio-owned organization and are required when
launching a tool, pipeline, or molecule import. The local API rejects job
submissions without a valid project. Rosetta and ipSAE show completed complex
results from the selected project only, and save each analysis job to that same
project. Their project selectors can also create a project in place. Existing
jobs without a project can be assigned one from the job details page.

Each job keeps its files in the Results folder and stores its project link in
the Studio marker and SQLite job catalog, so Results can search and filter runs
from different tools together.

New workflow submissions are registered with `mn-compute-scheduler` under the
`mn-molecular-design-studio` app ID. The Studio worker bridge starts only after
a user submits a workflow, claims only Studio workflow jobs, heartbeats active
leases, and passes scheduler allocations to the existing app-owned worker.
Scheduler records keep paths to the result and metadata files; structures,
logs, inputs, and result data remain in each run folder. Studio jobs wait for
the central API if it is unavailable; they are not silently sent to the local
queue, which could bypass shared resource accounting.
