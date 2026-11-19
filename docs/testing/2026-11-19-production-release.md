# Scientific browser and attached-analysis release gate

Date: 2026-11-19. Production activation was explicitly requested for pc002975,
PC004032 and PC004033. Remote login remains `compute`. No analysis is to be
submitted automatically, and no research records or credentials are to be reset.

## Verification

- 1,162 regression tests passed, 10 skipped, 12 warnings in 139.25 seconds.
  This covers every staged Python test file plus central-worker, workflow-contract,
  post-Phase-2 and complex-analysis integration regression tests. Execution used
  the network-disabled isolated release harness and separate writable data roots.
- The seven explicit offline BinderScope image cases were then enabled in a
  separate isolated-data run: all 49 BinderScope/core-hotspot tests passed,
  including actual container calculation, original-design SQL attachment and
  replicated-design ingestion. No research analysis was queued.
- Five JavaScript inspector tests passed: linked residue selection, local tabs,
  exact geometry overlays, confidence handling and untrusted labels.
- The wheel builds successfully, and the vendored 50-file source manifest and
  container adapter hashes were checked. The integrated runtime image is
  `mn-binderscope:c61ef9c`, local image ID
  `sha256:2ab83b3bbc723b2ebc7b170dcf4f7ad8e4df80dd8c8920630fe6ffdf65a535d4`.
- The application diff passes whitespace checks; imported upstream whitespace
  is preserved to retain the source inventory. Python compilation succeeds.

The central-worker dependency check exposed a name collision between a local
cancellation boolean and an imported predicate. The predicate now has a distinct
name; a successful monitored child is covered by a permanent regression test.
The MD bound-docking fixture now isolates its cluster home alongside its SQL root.
Browser tests explicitly request optional columns before asserting those columns.

## Coverage limits

This is not a claim that the complete application suite or every GPU/tool pipeline
was rerun. An earlier whole-suite attempt was stopped after 608 passing tests,
13 skips and two MD fixture-isolation failures; the fixture correction and its
dependent tests pass in the final regression run. Three MD trajectory cases still
skip because the host environment lacks optional `mdtraj`; the seven offline
BinderScope skips in the primary run are covered by the separate enabled run.
No real scientific job is started as a deployment smoke test.

## Activation requirements

Use one immutable committed source archive on all three machines, preserving
each machine's own database, node identity, authentication and configured paths.
Capture and verify per-host recovery copies before activation. Preserve the
existing primary UI on port 8502 and keep previously inactive compute-node UIs
inactive. Install and verify the same integrated BinderScope image on all nodes.
Do not connect production to the read-only preview database on port 8514.
Record full release SHA, process provenance, preservation checks and private
recovery locations outside the repository in the deployment receipts.
