# Production API access for preserved data owners

## Defect and correction

`StudioClient` treated the former development-owner account ID as proof that
the API was running in development mode. The supported production bootstrap
can deliberately preserve that ID so existing research data retain their owner.
Consequently, a correctly authenticated production instructor could use the
Streamlit application and REST API but be rejected by the generic REST/MCP client.

The client now requires the explicit boolean `development_mode` returned by
`GET /api/v1/auth/me`. Account IDs no longer determine authentication mode.
Development mode remains rejected for HTTP/MCP clients, including when local
development opt-in is enabled. Missing, malformed or non-boolean declarations
fail closed; older endpoints without this declaration must be upgraded rather
than silently trusted. The current API already supplies the declaration.

The same ID-only inference was present in distributed-job metadata. Job tagging
now also requires the execution process to be explicitly configured for
development mode. Production jobs retain their real owner without acquiring a
development marker. No command schema, account, password, access grant, scientific
record or authentication configuration is changed by this correction.

## Verification

- Twenty-five isolated client cases cover preserved production identities,
  ordinary identities, both transport modes, actual development-mode rejection,
  explicit local development opt-in and invalid identity/mode responses.
- Authenticated HTTP tests bootstrap a production instructor with the preserved
  ID, sign in through the actual login endpoint and reach the normal REST routes
  through `StudioClient` and ASGI HTTP transport, without dependency overrides
  or direct submission-handler calls.
- A 505-selector request validates and queues exactly one 32-CPU/no-GPU job.
  HTTP preparation remains deferred to the worker; an identical retry returns
  that job and a conflicting retry returns 409.
- Missing/invalid bearer tokens, unauthorized project access and unavailable
  structure selectors are rejected without queueing work.
- The real MCP HTTP transport accepts the preserved production instructor.
  Existing MCP tests retain cross-user identity and development-mode rejection.
- With `MDS_TEST_BINDERSCOPE_IMAGE=1`, an authenticated HTTP submission runs the
  actual offline BinderScope image on one synthetic saved complex and one CPU.
  Normal ingestion attaches one observation to the original molecule, source
  group, representation and checksum. Repeated ingestion/submission does not
  duplicate it; original membership and coordinate records remain unchanged.
- The API/client/BinderScope submission regression passes 123 tests with four
  numerical warnings from existing synthetic DSSP fixtures. Cluster regression
  passes 73 tests, including production/development job tagging for the preserved
  owner ID.
- The combined isolated release check passes 275 tests with the same four
  warnings, including enabled offline-image execution, preserved-owner bootstrap,
  analysis-owner identity and scheduler-heartbeat/recovery regressions.

All scientific test inputs and accounts are isolated fixtures. These checks do
not submit new research analyses, reset production credentials or exercise an
actual Streamlit mouse-click/browser session. They verify the authenticated
HTTP-to-worker-to-SQL path, not every application tool or GPU pipeline.
