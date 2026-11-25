# Capabilities acceptance, 25 November 2026

Tested the Studio `c8d9506` checkout with the acceptance patches. The live Studio MCP catalogue and the checkout both exposed 82 binder routes across 13 formats and three target contexts.

| Coverage | Verified |
| --- | --- |
| Preset graph contracts | 346 passed. These include 82 default binder route/context cases, 114 binder option variations, 144 other operation presets and starting inputs, five guided compositions, and one prediction input label check. |
| Browser setup under a regular user | BindCraft2 linear peptide, macrocyclic peptide and VHH reached the canvas. The linear peptide draft saved and reopened. |
| Browser submission and Results under a regular user | ESMFold2 Full with the optional Sequence Composition step, run `9F0D8`. The browser showed both worker steps as completed. The run produced one predicted structure and one analysis result for the same 20-residue synthetic sequence. |
| Preset and browser graph comparison | The saved linear peptide and monomer prediction graphs matched their corresponding preset node kinds and typed connections. |
| Related regressions | 669 passed and one was skipped in the host test run. The skipped MDTraj ligand protonation checks passed separately inside the GROMACS image. |

The [case checklist](capabilities-acceptance-2026-11-25.json) records browser and live evidence alongside each contract result. Per-format binder execution and most regular-user browser combinations remain open.

## Findings

- The live Capabilities page failed during import because the UI process still held an older module after an in-place update. The current files imported successfully in a fresh process. The primary Studio UI was restarted and its health endpoint responded. The regular-user browser replay used an isolated Studio instance.
- The Molecular inputs editor labelled the stored monomer selection as **Protein–protein**. The editor now uses the prediction goal labels. The saved monomer reopened as **Protein monomer**.
- The compatibility catalogue omitted antibody and peptide canvas nodes when classifying their existing workflows. The live MCP compatibility response now lists their graph routes.
- The first browser run, `3DBB6`, used the older editable installation in a worker subprocess and failed schema validation. The worker now launches from the checkout root. That first fixture left a failed test record in the shared Results directory. The replay used explicit isolated Results storage.
- A shared scheduler worker also claimed an isolated prediction and rejected its work directory. The test pipeline was cancelled and replayed against the local scheduler with one explicit GPU. This keeps the local browser case separate from production queue coverage. The local attempt using the **all** GPU selector failed before model execution.
- The generic Results header counted the structure and its Sequence Composition row as two final candidates. The step table and artifacts identify them as two records for one input.

The three source corrections were copied to the workstation deployments. The primary UI was restarted. The test suite and source changes remain in the isolated worktree awaiting a commit.

The [follow-up acceptance pass](capabilities-follow-up-2026-11-25.md) adds the full table sweep, per-variant widget tests and native library and database-selection handoffs.
