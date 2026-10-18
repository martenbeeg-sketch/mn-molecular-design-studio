# MN Molecular Design Studio

Molecular Design Studio is a separate project with its own code, workspace,
job catalog, Results folder, database, and organization data. Workflow modules
needed by Studio are kept inside this repository; the previous app is a
reference only and its job results are not imported.

The Results folder starts empty, defaults to the same data volume as the
previous app, and can be changed in Settings. The Database includes only
molecule sets imported or created through Molecular Design Studio.

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

The environment installs Molecular Design Studio and the sibling
`mn-compute-scheduler` client package. It does not install or edit the previous
protein design app. The scheduler API is a separate service and must already be
running before submitted Studio jobs can be claimed.

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
