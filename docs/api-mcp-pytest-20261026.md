# Full pytest verification · 26 October 2026

## Result

**1,334 passed, 15 skipped, zero failures/errors**, 17 warnings, in 107.88 seconds.
All 1,349 collected cases were included in the standard full-suite invocation.
The suite ran on the organizer with the Studio conda environment's Python:

```text
/home/user/mambaforge/envs/mn-molecular-design-studio/bin/python -m pytest -q --tb=short --junitxml=<temporary-workspace>/pytest.xml
```

Application, results, reference and cluster paths used private temporary storage.
The scientific workflow and scheduler contracts use small fixtures and mocked
external execution. The default optional live-container checks remain opt-in.

Detailed local artifacts:

- Log: `/tmp/studio-pytest-full-acGwmp/pytest.log`
- JUnit XML: `/tmp/studio-pytest-full-acGwmp/pytest.xml`

## Added coverage

`tests/test_api_extended.py` adds 25 cases covering:

- Kejia IDR and Peptide Trap setup schemas, CPU defaults and parameter constraints.
- Template browsing, coordinate downloads, private staging, examples and previews.
- Reference-library symlink protection and private input ownership.
- Every canvas starter through API draft creation, with project access checks.
- Completed-candidate continuation with candidate and source-project validation.
- Structure inventories and private saved camera/representation settings.
- Local job identity, results, hiding/restoring and resumable controls.
- Standalone resource discovery and distributed inspection without GET side effects.
- Connected-worker forwarding under the caller's account.
- Candidate viewers and private CSV, SVG and HTML exports.
- Convergence reports with top 100, 1,000 and 10,000 curves and recorded engine units.
- Separate bound-target/ligand layers, including legacy target links without copied files.
- Workstation administration, credential redaction and collection schema permissions.
- Actual MCP tool calls against the REST API using an isolated student identity.
- Authenticated binary download descriptors retaining their query parameters.

## Fixes found while running tests

- Local API resolution now includes the authorized catalogue identity and current
  hidden flag. This fixes source-project checks for continuation and job restoration.
- The shared Boltz-2 prediction adapter retains the generated binder's ligand CCD
  definition. Custom fragments remain rejected, and AlphaFold2 rejects ligand
  complexes. Ligand-only handoffs do not request a target MSA.
- Standalone fingerprint jobs can select their local scheduler without enrolling
  in workstation sharing. Unregistered remote workers remain rejected; live local
  scheduler availability and slot limits remain enforced.
- Kejia IDR and Peptide Trap are included in the shared generator output vocabulary.
- Test reference files/settings are explicitly isolated. The library workflow test
  uses real local fingerprint products behind a mocked distributed execution boundary.
- Distributed result tests use the existing `command_id` response field.

Additional regression checks cover Boltz-2 ligand handoffs with and without source
run provenance, and standalone fingerprint worker discovery and scheduler failure.

## Skips and warnings

- Three MD-analysis cases require optional `mdtraj`, which is absent from this
  environment.
- Ten quality-container cases require `MDS_QC_DOCKER_TESTS=1`.
- Two native RFD4 container cases require `MDS_RFD4_DOCKER_TESTS=1`.
- Fifteen DuckDB warnings concern deprecated `fetch_record_batch()` usage.
- Two pydssp numerical warnings arise from the synthetic structure fixture.

This run covers the complete standard pytest suite; optional live-container
checks and the three unavailable MDTraj cases remain unexecuted.
