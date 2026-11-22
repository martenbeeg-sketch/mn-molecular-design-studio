# Visual pipelines

**Pipelines** opens one canvas for project drafts. Open **Connections**, select a
step, and choose an output to see compatible tools and their input sockets. Click
a match to add and connect that tool, or drag a wire between existing sockets.
**Saved-result handoffs** lists tools that read a saved job or manifest. The app,
API and MCP use the same declared compatibility matrix. In MCP, use
`tool_connections` with `view="connections"`, `operation`, `direction` and `port`.
Run validation also checks the graph's inputs and scientific constraints.

Shift-select steps and click **Group** to name, collapse or duplicate them. A
collapsed group exposes its boundary sockets. It keeps the original runnable
steps and wires in the draft. **Open group** shows the contents on the canvas.

In the large-library screen, **Screening round** opens a block outline inside the
canvas. **Repeat docking rounds** frames the real stages. The If/Otherwise
decision either exits to Final results or prepares the next batch. Select a step
in the outline to edit that stage's settings. The screening controller manages
the recipe's wires.

The screening stack shows the selected centroid sample for round one, the
stopping criteria, and how the next batch is built. The selection block uses the
controller's sizing and rounding rules. Its counts assume enough unattempted
entries remain. The model-ranked and uniform random portions exclude previous
attempts from the screening library, and the random portion also excludes the
compounds already selected by the model for that round.

The screen uses two prepared datasets. **Predict** and **Select next docking
batch** refer to the same screening fingerprint library. **Centroid library**
draws the first-round sample from a completed, registered **Create centroids**
output. The centroid dataset may come from a different library with compatible
fingerprint settings. **Prepare library** and **Create centroids** are reusable
tools. The stages managed by the active-learning controller are available within
the complete **Large-library screen** recipe.

Open **Pipelines** and select a project. In the **Pipeline browser**, expand **Starter pipelines** and choose
binder design, peptide/IDR design, structure prediction and cofolding, docking
and library screening, molecular dynamics, or analysis of saved results.
Pick a method and click **Create editable pipeline** to open it on the canvas.
Under **Saved drafts**, select an existing pipeline to review its steps, then click **Open in canvas**.
The binder settings offer complete native, comparative and modular workflows.
Specialist recipes, the blank canvas and legacy templates are under **Advanced**.
The MCP `pipeline_templates` tool identifies main starters, variants with their
parent IDs, and advanced starters. Existing template IDs remain callable.

The default binder recipe connects Target → Binder Design Context → RFdiffusion3 →
initial ProteinMPNN → Boltz-2 → Results. Backbone methods require initial sequence
design; co-design methods retain their sequences. Native BindCraft2, BoltzGen,
Complexa, Genie3 and PXDesign pipelines keep their internal stages and do not
receive redundant external redesign or prediction by default.

**Compare selected methods independently** creates separate branches, including
sequence, prediction, optional selection and Results nodes. Shared attempt counts
are not equal runtime or native acceptance budgets. Successful branches remain
available when another method fails; partial comparisons are labeled explicitly.
Failed measurement rows remain inspectable and do not count as successful branches.
Raw model ranking scores are not calibrated across engines. Workload summaries
report encoded attempt/sample counts, without inventing runtime estimates.

Molecular Structure Prediction has separate monomer, protein-complex, ligand,
DNA, RNA and cyclic-partner modes. Method choices respect molecular capabilities,
ligand notation and target-template requirements. It defaults to one engine;
samples, seeds, comparison and advanced engine settings remain configurable.

Optional SolubleMPNN retains parents alongside redesigned candidates and protects
the mapped interface, fixed motifs, antibody framework and cysteines using the
Phase 4 implementation. It is disabled by default and does not guarantee improved
experimental solubility. Optional evidence-based selection uses the existing
filter and Results infrastructure. Additional validation can run after an explicit
shortlist. See [Phase 5 implementation and validation](testing/task-catalogue-consolidation.md).

## Working on the canvas

Under **Molecular dynamics**, choose **Analyze a saved MD run** or **Prepare a
complex and run OpenMM**. Both starters connect the exact trajectory and
topology to the CPU analysis and Results. The simulation starter includes the
existing complex preparation tool. Connect its
MD measurements output to Results; branch the saved source to compare sampling
windows/settings. Source files stay read-only on the pipeline workstation.
This route calculates RMSD, Cα/ligand RMSF, mass-weighted Rg and geometric
contacts; it does not automatically perform GB/PB endpoint calculations or
chemical hydrogen-bond classification. See [MD analysis contracts and limits](md-data-browser.md#executable-md-analysis-tile-and-pipelines).

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

## Legacy RFD4-Proteina starter workflows

The following IDs retain their historical settings for compatibility. New task
recipes start without example-specific atom or motif selectors; supply constraints
for the actual input through task settings or the existing RFD4 node form.

The starter menu includes three RFD4-Proteina paths:

- **Protein binder → Boltz-2 → interface analysis** connects a Target and Binder Design Context to RFD4, then passes its co-designed sequence and structure through Boltz-2, ipSAE, and Rosetta metrics.
- **Small-molecule pocket design → Boltz-2** starts from a protein–ligand mmCIF, runs RFD4's small-molecule binder task, then predicts the generated complex independently with Boltz-2. Its ligand restraints are prefilled for the tutorial's 5ENI/5QA example.
- **Enzyme motif scaffolding → Boltz-2** starts from the tutorial's 5ENI/5QA example. Its motif residues, partial atom anchors, ligand chlorine restraints, and hydrogen-bond roles are prefilled; edit them when using another structure.

RFD4 generates sequence and structure together, so these starters do not add ProteinMPNN. The task-specific `options_json` field is multiline. For atom-level ligand constraints, use options such as `ligand_atoms_to_hold`, `ligand_interactions`, `ligand_hbond_donors`, and `ligand_hbond_acceptors`; enzyme motif workflows use `motif_residues`, `motif_atoms`, and `active_site_atoms`. The standalone RFD4 form exposes these controls directly.

The [Python API tutorial](https://github.com/RosettaCommons/RFD4-Proteina/blob/main/docs/examples/python_api_tutorial.ipynb) demonstrates building an annotated CIF with pocket residues, fixed and tip atoms, ligand hydrogen-bond roles, and coordinate variants. Its enzyme starter uses those same example selections and expects the 5ENI/5QA complex; the small-molecule starter uses the ligand atom restraints from the tutorial. Change chain, residue, and atom names if you use another input. You can also pass your own annotated CIFs to the RFD4 node by setting **Input method** to **annotated_cif** and listing their paths in **Input structure files**. The [Hugging Face ZeroGPU Space](https://huggingface.co/spaces/nvidia/rfd4-proteina-zerogpu) uses AtomWorks Studio to edit annotated CIF conditions before calling the native sampler; these starters borrow its edit-to-generate pattern and add Studio's own downstream prediction and analysis.

## Drafts, versions, and results

Drafts autosave within the selected project. **Save version** keeps an immutable copy;
restoring a version replaces the draft without changing that version or an earlier run.
Use **Archive draft** to take a recipe out of the active list. The **Archived drafts**
panel restores it with its versions intact. In MCP, call `list_saved_pipelines` with
`include_archived=true` and use `archive_pipeline_draft` with the current revision
and `archived=false` to restore it. Saved versions and Results can still be read
while a draft is archived. Restore it before running it again.

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

## Binder and IDR outcomes

Open a completed design or pipeline in **My Results** and expand **Binder design
outcomes**. **Native filter** reads the generator receipt. **Final acceptance**
reads the recorded validation or campaign policy. **Database** checks the
candidate's run lineage, project and stored structure digest. A missing decision
or a different saved representation is shown explicitly. Select a row to see
its recorded decision details and open a matching BinderScope assessment.
The link requires the saved representation, input digest and chain roles to agree.
A saved refold with matching lineage is the default design model in Data Browser.
When its BinderScope assessment is available, Results opens it first as
**Refolded model**. Other matching assessments can be opened on their own models.
A remote job shows this information after loading full Results from its workstation.
BindCraft2's ranked refolded complex is named in the viewer. Historical records
use this label when the producing result identifies the same saved structure,
candidate, run, and checksum as the ranked output.

The local MCP `result_presentation` tool returns the same per-candidate
`binder_outcome` evidence, including the preferred assessment's observation ID.
Follow `next_offset` until it is null. MCP retries with a smaller page when a
response exceeds its size limit.

## Follow-up analysis from My Results

- Choose a saved **Docked pose** and open **Analyze this pose**. The review
  page checks the score row, retained coordinates, receptor, and consumed
  ligand record before offering PLIP, PandaMap, or PoseBusters.
- Choose a completed **Production result** and open **Analyze trajectory**.
  The MD Analysis form opens on that segment and shows the matching
  trajectory and topology. Remote inputs stay on their owning workstation.
- Open **Analyze selected structure** next to a prediction or cofolding
  model, or select a binder in **Binder design outcomes**. The binder
  model picker prefers a verified saved refold. The inspector checks
  compatible analyses and reports when the selected model lacks PAE.

The review page opens before submission. After queueing a pose or MD analysis,
use **Open analysis Results**. A retry from the same review session uses the
original request key.

## Verification

The graph/API/queue tests use temporary files and simulated child jobs; they do not launch
GPU predictions. They cover dependency validation, branches, failure propagation, target
locking during sequence design, scoped inputs, immutable versions, reuse, cancellation,
queue adoption, and Streamlit save/restore handling. The live browser was checked for
block dragging, socket wiring, parameter edits, undo, saving, and version restoration.
