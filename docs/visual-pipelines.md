# Visual pipelines

Open **Pipelines → Node canvas**, select a project, then choose a starter or a blank canvas.
The Foundry starter connects Target → RFdiffusion3 → Foundry MPNN → RosettaFold3 → Results.

## Working on the canvas

- Drag a block to any position. Drag between its sockets to connect steps.
- Branch an output to compare downstream engines; Results can collect several branches.
- Click a block to open its input form. Main controls and **Optional settings** start with the saved node values. **Reset to tool defaults** loads defaults from the same workflow schema as its standalone tool page.
- **Apply to pipeline** saves the form to the draft; **Cancel** keeps the existing values. Applied engine settings can be undone on the canvas.
- Click a Target block to use the familiar structure upload/fetch, chain, trim, and hotspot editor. The Successful project designs block opens the project’s eligible design selection. **Target & inputs** above the canvas also provides shortcuts to these forms.
- Connected upstream structures and sequences appear as inherited inputs in engine forms; editing the engine settings preserves those connections.
- Scroll to zoom; drag the background to pan; press **F** to fit.
- Shift-select several blocks. **Ctrl+D** duplicates; **Delete** removes selected blocks or wires.
- **Ctrl+Z** undoes; **Ctrl+Shift+Z** redoes. Double-click a wire to disconnect it.
- **Focus canvas** expands the working area. **Map** toggles a small overview.

The library distinguishes a design engine (generation only) from a complete design pipeline.
Sequence design, prediction, filtering, and collection are explicit downstream steps.
Connections are checked for compatible inputs, structure formats, and cycles.
Not every studio engine has a canvas adapter yet; the library lists the supported steps.

## Drafts, versions, and results

Drafts autosave within the selected project. **Save version** keeps an immutable copy;
restoring a version replaces the draft without changing that version or an earlier run.
**Run pipeline** submits a CPU coordinator and ordinary child jobs through the existing scheduler.
It freezes the graph and target input for that run and preserves original design artifacts.

With **Reuse unchanged completed steps** enabled, a new run can reuse successful outputs
from this user and project when settings, upstream inputs, target contents, and adapter
implementation match. Layout and labels do not affect reuse. A setting change invalidates
that step and its descendants. Disable reuse to request fresh work.

Only successful design jobs with explicit filter-passing designs are offered as existing
design inputs. No surviving candidates skips downstream computation. A failed step blocks
its descendants while independent branches can finish. Canceling a pipeline also cancels
its active owned child jobs. Previous runs remain available in the results section.

The Results block collects candidate structures and metrics for the studio results UI;
it does not run a separate new RMSD calculation by itself.

Selections of compounds, designs, sequences and candidates remain collections at
every connected calculation step. Cofolding Inputs can reference a Database
compound snapshot or a FASTA sequence collection; all selected records are passed
to the engine's native batch runner in one job. See [automatic collection execution](automatic-collections.md).

## Verification

The graph/API/queue tests use temporary files and simulated child jobs; they do not launch
GPU predictions. They cover dependency validation, branches, failure propagation, target
locking during sequence design, scoped inputs, immutable versions, reuse, cancellation,
queue adoption, and Streamlit save/restore handling. The live browser was checked for
block dragging, socket wiring, parameter edits, undo, saving, and version restoration.
