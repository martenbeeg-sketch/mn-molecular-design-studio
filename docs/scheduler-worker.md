# Scheduler worker paths and recovery

The MN Compute Scheduler is a separate, shared service. Molecular Design
Studio starts its own worker bridge when a job is submitted; the bridge claims
only Studio jobs. Do not restart or edit the shared scheduler to fix a Studio
worker-path problem.

## Which path setting controls a job?

Use **Settings → Results directory** in Molecular Design Studio. The saved
value is `results_directory` in:

```text
~/.local/share/mn-molecular-design-studio/settings.json
```

If `MN_MOLECULAR_DESIGN_STUDIO_HOME` is set, substitute that workspace for
`~/.local/share/mn-molecular-design-studio`. Both the Streamlit app and API
must use the same `MN_MOLECULAR_DESIGN_STUDIO_HOME` value. At runtime Studio
loads the saved Results directory into `MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR`;
the API and worker bridge use that effective path. Let Studio set this variable
from Settings rather than exporting a different value in one shell.

The worker's singleton lock and log live under the stable Studio workspace,
not under Results, so changing Results does not create a second bridge or leave
the active lock behind in the old Results directory. The bridge reloads saved
workspace paths at startup and before each new scheduler claim.

## Read-only diagnosis

Check Studio jobs and scheduler jobs before taking recovery action:

```bash
curl -fsS 'http://127.0.0.1:8766/api/v1/jobs?statuses=running&limit=100'
curl -fsS 'http://127.0.0.1:8765/v1/jobs?app_id=mn-molecular-design-studio&status=running&limit=100'
curl -fsS 'http://127.0.0.1:8765/v1/workers?app_id=mn-molecular-design-studio'
```

Compare the worker's `labels.results_root` with the Results directory shown in
Studio Settings. New worker registrations include both `app_home` and
`results_root`; older registrations may only show `app_home`. For an older
bridge, inspect its process environment for
`MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR` and compare it with the saved setting.

## Safe recovery for an old, misconfigured bridge

If an old bridge still points at a different Results directory:

1. Wait for Studio jobs already running on that bridge to finish. Confirm there
   are no Studio jobs in `running` state in either read-only query above.
2. Verify the PID belongs to
   `mn_molecular_design_studio.engine.workflow_engine.core.central_scheduler_worker
   --service` and that its environment has the stale Results path.
3. Send `SIGTERM` to that worker bridge PID only. Its handler stops new claims;
   it does not stop or restart MN Compute Scheduler. Do not use broad `pkill`
   patterns or signal the scheduler API process.
4. Submit the next Studio job normally. Studio starts the worker bridge with
   the saved Results directory and records the new root on the worker.

If a Studio job is still running, leave the bridge alone and wait. Never force
stop the scheduler or delete the worker lock while a Studio run is active.
