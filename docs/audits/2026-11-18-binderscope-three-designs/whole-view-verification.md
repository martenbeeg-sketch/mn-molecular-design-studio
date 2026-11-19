# Whole-view table and visible analysis action

19 November 2026 — development preview `browser-tree-v66-whole-view-binderscope`.
No production deployment, authentication changes, data edits or analysis submission.

## Changes

- BinderScope scope selection and Analyze button are above the table in all
  workspaces, not inside the results/options popover.
- All designs in this view includes the current folder, search and refinements
  across every page, independently of the display-window size. The dialog checks
  exact complexes/settings and exclusions before any submission.
- Rows shown defaults to All for views up to 1,000 entries. Users can restore
  paging at 50, 100, 250 or 500 entries. Larger views retain paging. Changing the
  row window resets selection rather than transferring row positions incorrectly.
- Table and row-action identities include the display window. Only chosen
  scientific columns and row actions are sent to the browser; full evidence is
  retained for inspection, optional columns and complete CSV exports.

## Read-only benchmark: real PD-L1 project

506 admitted design entries, production SQLite opened read-only. Timings are
local observations, not guarantees for every client or folder.

| Measurement | Seconds |
| --- | ---: |
| Initial accessible SQL index | 5.912 |
| Initial PD-L1 folder inventory | 3.056 |
| First 506-row table projection with existing schema | 2.823 |
| Repeated projection using cached cohort metadata | 0.119 |
| First projection in repeat run | 2.604 |
| Repeated projection in repeat run | 0.135 |
| Default-column table serialization, repeat run | 0.013 / 0.009 |

Serializing every field, including hidden report/provenance metadata, transferred
3,814,696 Arrow bytes. The default 16 scientific columns transferred 136,024 bytes
before adding the small row-action column: approximately 96.4% less. Normal UI
reruns additionally reuse the prepared page/window table. These figures do not
include network/rendering latency or the separate per-complex analysis preflight.
The native grid scrolls through the metadata; it does not instantiate 506 3D viewers.

SQLite SHA256 stayed unchanged before/after both benchmark runs:
`9fd3917c01cb1bc2818d9eb8117af544431e54ae1b10b280569379bf89e9f3e4`.

## Validation

173 targeted Python tests passed, including the offline SQL/peer/hotspot tests.
Coverage includes 506 entries in All mode, restoring page 3 at 50 rows, stable
display-window identities, the 1,000-row guard, analysis selection across 506
matching rows, the visible analysis action outside popovers, source-bound row
inspection and keeping hidden report JSON out of the default browser table while
retaining full CSV evidence. Formatting/compilation checks passed.

The authenticated preview shows All 506 with no Page control, the explicit
all-pages analysis scope and an enabled Analyze · BinderScope (506) button. Opening
that button starts input checks only; final job submission remains disabled in
the read-only preview.
The live check opened the all-designs input-validation dialog but did not submit
an analysis or establish eligibility of all 506 complexes. Per-complex input
validation is slower than displaying the table and is not included in the table
benchmark. Screenshot: `binderscope-whole-view-analysis.png` in the chat
visualization directory.

## Follow-up: clearer selection and binder format

19 November 2026 — development preview `browser-tree-v68-binder-format-labels`.
The toolbar is one row: result/window control, analysis scope and the visible
Analyze button. Scope explanations and the all-pages/checked-subset choice live
in the scope popover rather than a permanently displayed radio group and captions.
Scope-widget identities are bound to the checked design IDs so stale state from
a previous subset cannot silently broaden the next analysis.

Preflight defaults to unassessed complexes only. Completed-report counting checks
the saved representation, input/report digests, current group access, project and
chain selection, plus the explicit target association when supplied. Failed,
pending, other-model/chain, conflicting/missing-report and other-project evidence
does not suppress analysis. Each source's collection identity remains separate.
The editable Analyze checkbox, Include all, Unassessed only and Clear selection
allow explicit overrides without launching work. Viewing saved reports remains a
separate row action with no checkbox requirement.

The settings column is Binder format, with Protein, Linear peptide, Cyclic peptide
(head-to-tail) and Not recorded labels. The old internal parameter strings remain
compatible; this presentation correction neither reclassifies saved molecules
nor infers format from sequence length.

184 targeted tests passed, including real temporary-SQL completed-assessment
defaults, explicit inclusion/restoration, project/access/checksum exclusions,
checked-subset/all-pages scope changes, UI labels and all previous table/inspector
regressions. No production deployment, authentication or scientific-data changes;
no real analysis was submitted.

Authenticated real-data UI check: PD-L1 design
`Design_pc002975_20261101T221447Z` showed one completed assessment, an unchecked
Analyze box and zero selected by default. Include all changed the selection to
one; Unassessed only restored zero. The displayed format was Protein. Submission
stayed disabled in the read-only preview. Screenshot:
`binderscope-protein-unchecked.png` in the chat visualization directory.
