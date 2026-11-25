# Capabilities follow-up, 25 November 2026

The browser sweep used a signed-in regular user in an isolated Studio workspace. All 160 visible table actions opened a setup dialog with an enabled **Create pipeline** action. The numbers count table buttons. A button can expose more than one workflow variant.

| Table group | Setup dialogs checked |
| --- | ---: |
| Protein binders with a structure, multidomain target or disordered sequence | 68 |
| Other protein design tasks | 24 |
| Structure prediction and cofolding | 10 |
| Docking and affinity | 7 |
| Structure preparation | 2 |
| Molecular dynamics | 2 |
| Library preparation | 2 |
| Conformational sampling | 3 |
| Small-molecule design | 12 |
| Analysis and validation | 30 |

The folded or unknown sequence context correctly showed an empty table. The Streamlit widget suite separately built, saved and reopened all 82 binder variants, including the variants sharing a button. These tests exercised the real setup code with a simulated table event. Native binder generation was outside this run.

## Native pipeline runs

| Case | Run | Evidence |
| --- | --- | --- |
| Prepare library → Create centroids | `8D088` | Built from Library preparation in the student browser. The fingerprint step retained four of five input rows and recorded the invalid SMILES as a rejection. The centroid step used the connected fingerprint manifest and the same library ID, producing four centroids. Both steps completed. |
| Saved prediction → Protein-Sol | `B2C70` | The named ESMFold2 output carried one 20-residue sequence into analysis. The worker recorded Protein-Sol’s minimum-length rejection and the browser displayed the failed run. |
| Database selection → Protein-Sol → Results | `2ED48` | The same frozen 40-residue database fixture passed readiness after the lookup repair and completed after the evidence export repair. The regular user saw one result record. The final evidence CSV includes the source molecule reference. |

## Fixes

- The sticky engine column could cover a button when the table was scrolled horizontally. The table reserves the column’s width when scrolling to focused controls. The previously covered button was visible and clickable after moving focus back from the far-right column.
- The Prepare library editor accepts one file path per line and shows the configured output root. The centroid checkbox was removed from this editor because centroid preparation has its own connected step.
- The pipeline API validates a frozen database selection in its recorded task group and still checks access, the exact directory and project membership. A cross-project snapshot failed the regression test.
- Results and CSV export format structured molecule parent references from existing snapshots. The Results header counts result records, which can include analyses and predictions from the same molecule.

The [case checklist](capabilities-acceptance-2026-11-25.json) distinguishes graph coverage, browser setup and native execution. Native binder-format runs, live execution of every optional branch, and a full browser save-and-reopen for every non-binder combination remain open.

The final regression run passed 625 tests. Six source files were synchronized to all three Studio workstation installations. The primary API and UI restarted and both health checks passed. The changes remain uncommitted in the isolated worktree.
