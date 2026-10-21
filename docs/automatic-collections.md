# Automatic collection execution

Selecting several independent compounds, designs, sequences or candidate structures
means one batch for that tool/engine and pipeline step. There is no batch-mode flag.
Chains inside one molecular complex remain one calculation; they are not split by
their chain count. Seeds and samples remain settings of each calculation.

## Inputs and execution

- Prediction tiles expose input sources directly: new molecular inputs, compound
  collections, sequence collections, Database molecules, or selected designs.
  Target-template-only engines retain their existing template/design selectors.
- FASTA records and plain sequence lines describe independent calculations. A
  colon keeps multiple chains together within a FASTA record.
- A compound collection uses the shared protein target for every compound.
  Boltz-2 collection predictions default to potentials and affinity prediction.
- A compound selection runs in full by default. The optional compound limit is
  `0` for all selected compounds; an explicit positive limit selects a subset.
  A selection exceeding the tool's supported run limit is rejected, not truncated.
- Pipeline cofolding inputs accept `compound_set_file` and `sequence_file` and
  expand them through the same native prediction adapter. Database protein
  selections can be predicted from their sequences without requiring a PDB.
- Existing design validation, MPNN, inverse folding and analysis workflows keep
  their candidate collections intact when passed between pipeline steps.

Native collection interfaces are preferred: the prediction images receive one
`mn-predict batch` manifest (ESMFold receives one combined FASTA), Nesso receives
an input directory, Uni-Dock Pro receives a ligand index, AutoDock-GPU receives a
file list, RosettaLigand receives ligand chunks, and Vina/GNINA run bounded workers
inside one container. Existing images already provide these interfaces; no image
rebuild is needed for the shared submission change.

## Shared API contract

Normal workflow requests still accept a parameters object containing a native
collection such as `tasks`, `validation_sources`, `candidates_jsonl`, `batch_inputs`,
`compound_set_file` or `sequence_file`.

For queued calculation workflows, `parameters` can also be an ordered list of
input objects. This automatically creates one scheduled job:

```json
{
  "schema_version": 1,
  "workflow": "refolding.optimized_prediction",
  "project_id": "PROJECT_ID",
  "parameters": [
    {"engine": "boltz2", "protein_sequences": "ACDEFGHIK", "ligand_smiles": "CCO"},
    {"engine": "boltz2", "protein_sequences": "ACDEFGHIK", "ligand_smiles": "CCN"}
  ]
}
```

Matching prediction settings collapse into a single native prediction batch.
Compatible compound selections with the same receptor, pocket and engine settings
are merged before invoking the docking/affinity engine once. Other queued scalar
adapters execute their items in order inside one parent job and one allocation;
they can still require a container invocation per item. Grouping alone does not
give an engine a native in-process batch interface that it does not implement.
Orchestrators that already schedule multiple stages retain their native collection
fields rather than accepting arrays of entire orchestrator requests.

The scheduler reserves the maximum resources needed by an item, not the sum of
resources for every selected item. Native per-compound workers share that allocation.
An ordered fallback batch uses that allocation serially. Each item is validated
before any job is created, and normal input-path/project authorization still applies.

## Results, failures and cancellation

One Results entry owns the selection. `batch_manifest.json` records the collection.
Scalar fallback artifacts live below `artifacts/batch_items`, so purging the parent
removes the whole batch. Per-input statuses and errors are shown together. Failed
items do not prevent later inputs from running; a partly successful batch retains
its successful candidates for downstream steps and reports the failed count.
Cancellation prevents remaining inputs from starting and retains completed outputs.

The collection manifest is independent of the machine's current resource report;
it describes logical inputs rather than additional scheduler reservations.

## Verification

Regression tests use the application's Conda environment and simulated execution.
They verify one API/queue job for every installed prediction profile, full compound
selection defaults, native compound merging, pipeline collection expansion, chain
semantics, per-item failure isolation and resource accounting. No scientific GPU
job is launched by these tests.
