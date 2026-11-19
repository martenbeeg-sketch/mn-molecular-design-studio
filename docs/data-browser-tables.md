# Data browser tables and CSV

The tree shows scientific methods, recorded dates and stable experiment numbers.
Import history remains in **History** and original collection/run identifiers;
“Imported study” is not used as a navigation label. These presentation changes do
not rename SQL identities, merge studies, or assert a filter pass.

PDB target folders use one activity-based route: target → docking → tool →
experiment, or target → molecular dynamics → study → replica. Protein constructs
are not parallel folders containing the same results. **Refine this level → Target
construct** selects an exact recorded protein target when multiple constructs are
present. The default **All constructs** also retains results whose exact construct
was not recorded; selecting a construct does not guess such a link. Exact target
identifiers remain distinct, with **Target construct** and **Target sequence**
available in the table and CSV. The existing Database page is unchanged.
Construct labels use the calculation's recorded receptor/PDB context, so a shared
sequence named under another PDB does not relabel the selected experiment.

## Navigation performance

Project/root outlines are precreated from the current accessible metadata, without
reading result inventories. Their arrows expand/collapse in the browser without a
server rerun. Result inventories are prepared on first selection and reused when
returning to a folder; activity/tool branches and counts are also reused for the
same snapshot. The session-local inventory cache holds at most eight scopes and
32,000 raw evidence rows. Oversized folders still open but are not retained there.

Every server interaction still checks fresh group access and the SQLite data version.
Local/peer edits, deletion, permission changes, database replacement and reopening
the database watch invalidate cached views. Nothing is shared between users or
persisted as a second scientific catalog. Very large folders and the first index
load can still take time; **Refresh SQL view** discards the session cache explicitly.
PDB-specific observation queries exclude other receptor contexts in SQL before
decoding rows, while keeping structure/source evidence needed for legacy links.

## Flexible workspace

**Workspace → Split** retains the original folders/table/details layout.
**Table** hides the details pane and gives the table the available width;
**Structure** hides the table and expands the structure/details pane, with a
taller 3D viewer. Focus modes automatically unpin folders. **Pin folders** can
keep the tree beside either focused pane, or unpin it in Split to widen both
remaining panes. When unpinned, **Folders** opens the same navigation and
refinements as an overlay, without leaving the current workspace.

Folder, search, page, table row selection and chosen structure are retained
across layout changes. Structure focus offers **Result on this page** to inspect
another result without restoring the table. Batch checkboxes remain a separate
selection for analysis; inspecting another result does not add it to that batch.
Search, column choices and CSV controls remain in **Results, filters & analysis**
in focused workspaces. A compact toolbar above the table shows the result count,
analysis scope and **Analyze · BinderScope**. Open the result-count control for
**Rows shown** and paging. MD studies retain their analysis/production selectors
under **MD analysis & snapshot selection**. Hidden details are not rendered in
Table mode. These are session-local presentation choices, not SQL changes.

The middle table always offers project, method and linked target, even when all
rows share those values. **Table columns & CSV export** adds native result fields,
recorded scientific input parameters, replicas, attempts and source identifiers.
Sequence/SMILES, linked target identifiers and recorded PDB references are also
available as columns and in the CSV.

### Whole-view tables and batch analysis

**Result count → Rows shown → All** displays the complete matching view without paging for up to
1,000 entries and is the initial choice for folders within that limit. The grid
scrolls rather than creating a separate structure viewer for each row. Switch to
50, 100, 250 or 500 rows to restore paging. Larger folders retain paging to limit
transfer and memory. Changing the row window resets checkbox selection; changing
the workspace alone retains it. Search and refinements still define the same view.

**Whole view · N designs → Analyze · BinderScope** reviews every matching admitted
design, not just loaded rows. Checking rows at the far left of the table switches
the scope to **N checked designs**; clearing them restores the whole view. The scope
control also lets you explicitly choose all pages while keeping a checked subset.
Table checkboxes choose a batch, not pass/fail or completed-analysis status. Opening
a saved report uses the row's **View BinderScope** action and needs no checkbox.
The dialog checks
exact complexes and settings and lists exclusions before submission. Only eligible
protein-design complexes are analyzed; docking ligands and MD records are not
silently included. Opening the dialog does not launch jobs. Read-only development
preview allows inspection but keeps final submission disabled.

In that dialog, **Analyze** checkboxes are the final queue selection. Designs with
completed, verified BinderScope reports for the exact saved complex and chain
selection in the same project start unchecked. Failed/incomplete assessments,
other models or component selections, and missing/conflicting report digests do
not suppress a new assessment. This is not automatic protocol/version reuse.
Individually tick an assessed design, choose **Include all**, or restore
**Unassessed only**; none of these actions queues work until **Queue BinderScope**
is confirmed. Existing results and original admission decisions stay unchanged.

The settings column is labelled **Binder format**, not topology: **Protein**,
**Linear peptide**, **Cyclic peptide (head-to-tail)**, or **Not recorded**.
Linear refers to backbone connectivity, not an unfolded or
disordered conformation. Disulfide/side-chain loops do not imply head-to-tail
cyclization. Historical parameter values remain unchanged for compatibility;
renaming these UI labels does not reclassify molecules or rewrite prior results.

The table sends only chosen columns to the browser. Hidden provenance/report JSON
remains server-side until requested through a column choice, inspection or full
CSV export. SQL summaries are batched and cached for the displayed window; cached
metadata is still invalidated by current access and database changes. Initial SQL
index/folder preparation and analysis-input validation are distinct from the cost
of displaying or scrolling 500 table rows.

Automatic column presets follow the selected activities: designs prioritize
admission and generator confidence; docking prioritizes native affinity/pose
scores; prediction/cofolding prioritizes confidence; binding-affinity calculations
show native affinity outputs; MD prioritizes recorded duration and RMSD summaries.
The preset can be overridden, including **Compact** and **All columns**, without
changing the tree or workspace layout. Missing values remain blank; zero
and false are preserved. Metrics absent from SQL are not calculated from files.

An MD study is labelled by its recorded simulation engine, not the namespace of
an arbitrary analysis observation. Production is the representative when present;
MM/GBSA, MM/PBSA and MD-analysis records remain attached to the same study. The MD
table shows whether production is recorded, replica and distinct production-segment
counts, and saved analysis types. A displayed segment duration is not a cumulative
study duration: continuation segments are not summed without verified lineage.
Endpoint-only imports explicitly say **Production: Not recorded**. The expanded
evaluation CSV retains each calculation's own method and original observation ID.

Each table row uses its displayed original Evaluation. Scores from other models,
tools are not merged or ranked into it. Native scales are
not automatically converted between tools. Actual recorded per-row seeds take
precedence for display, with a differing configured seed retained as
`Input · requested_seed`. Admission-only measurements are separately labelled
**Admission evidence**, never reconstructed as a new Evaluation.

## Attached analysis

Post-analysis is loaded in one batched, accessible molecular-cohort query for the
visible page, using the same design/target relationship contract as the details
pane. Each tool has its own status and native score columns. Numeric cells require
one completed evaluation on the selected structure, with no conflicting target or
declared component selection. Analysis on another model, a derived structure,
sequence-only evidence, unknown status and unresolved links is labelled but does
not supply a same-structure numeric score. Multiple completed evaluations are not
collapsed, even if their scores agree. Design-internal refolding and filter
decisions are not presented as post-analysis.

MD analysis additionally requires an exact recorded production-segment link,
replica and study in the same accessible owner/project/source-node scope;
sharing a compound or a static structure is insufficient.
Its status says **Recorded MD segment**, rather than asserting a static-structure
measurement. Different analysis windows remain separate evaluations.

The absence label is **Not recorded**, not “not analysed” or “passed.” Every linked
evaluation, its version, parameters, source reference, relationship and original
metrics remain in the per-tool JSON **Records** CSV column. A dedicated scored
evaluation identifier identifies the source of each populated numeric column.
Reference/input library rows say **Not applicable** rather than borrowing the
post-analysis of a design or pose elsewhere in the same molecular cohort.

### BinderScope from the table

Click **View BinderScope** directly beside a design's name in the table. The
name and viewer column stay visible while scrolling horizontally. No checkbox
selection is required, and inspection does not change the analysis batch.
Rows with multiple assessments offer **Choose assessment (N)**; rows without a
verifiable saved report have no action. The same popup opens in Split or the
full-width Table workspace without a separate toolbar or details pane. Structure
workspace retains its focused-result shortcut. These are UI actions, not extra
scientific columns in the prepared CSV export.
Batch selection for **Analyze · BinderScope** remains a separate action.
Opening a saved assessment never queues a calculation or changes
admission. Multiple completed same-model assessments require an explicit choice;
neither the newest nor the best is silently selected. Access, current completion
and exact model relationship are checked again when opening.

Automatic/activity presets include BinderScope status, buried interface area
(Å²), H-bonds, salt bridges, buried unpaired polar atoms and contacted/requested
hotspots when recorded. **Table columns & CSV export → Show optional BinderScope
details** adds assessed model, analysis version, mean binder pLDDT, N/C tag
accessibility, terminal pLDDT/confidence and target-reference Cα RMSD. Individual
columns can still be chosen manually. Compact remains compact; All columns and
CSV include all available details regardless of the toggle.

These are separate post-analysis columns, never substitutes for generator
confidence or filter evidence. Missing assessments say **Not assessed**; missing
numbers remain blank, while zero contacts and zero unpaired atoms remain zero.
Unspecified/unavailable hotspots are labelled, not converted to 0/0. Confidence
requires an explicit recorded pLDDT interpretation; no B-factor guessing occurs.
Tag accessibility retains family-specific differences; “free” is a geometric
prediction, not an experimental verdict. Cyclic termini are not applicable.

The page reads a small batched summary from saved SQL reports, without loading
coordinates or full reports per row. A conflicting report/SQL checksum does not
supply scored cells or an inspector shortcut. Legacy recorded metrics without a
complete report remain visible but cannot open a guessed report. No schema or
production-data migration is needed. Full report/coordinates are read only on
explicit inspection and retain the existing access and checksum checks.

Use **Prepare CSV**, then **Download CSV**:

- **Current listing** exports every matching row across all pages, including all
  available columns, regardless of which columns are displayed.
- **All original evaluations** expands those same selected results into separate
  original evaluation rows for their recorded models/replicates. It does not add
  post-analysis entries or unrelated calculations outside the selected folder.

CSV uses UTF-8 with a BOM, native field names and unrounded numeric values.
Arrays/objects occupy JSON-valued cells. Strings that could execute spreadsheet
formulas receive a leading apostrophe; negative numeric measurements remain
numbers. The complete CSV is prepared only on request. Normal navigation batches
SQL reads for the visible page rather than hydrating the whole database.
