# MPNN implementation handoff · 2026-10-17

Implementation, pytest and the requested MPNN smoke runs are complete. The user
authorized the smoke after switching models. Keep the V7 campaign unchanged.

## What is implemented

- A single **ProteinMPNN** tile replaces the two catalogue entries. Legacy
  `sequence.ligandmpnn` and `sequence.foundry_mpnn` APIs remain available.
- The new public entry point is `sequence.mpnn`, implemented in
  `engine/workflow_engine/workflows/mpnn.py`.
- `mpnn_ui.render_inputs` is the common input form for the tool and both saved
  MPNN node kinds. Foundry, Proteina-Complexa, PXDesign, RFdiffusion, Genie3 and
  Protpardelle starter templates all use it. The redundant Foundry MPNN node is
  hidden from the palette but retained in saved graphs and the Foundry template.
- Tamarind-inspired green controls, PDB/CIF selector, fetch/upload, designed
  chains and residues, linked structure/residue viewer without the right tools
  panel, PDB trimming, sequence count, model pills, optional settings, project
  selection and existing-results batch selection.
- Model choices: ProteinMPNN, LigandMPNN, SolubleMPNN, HyperMPNN, AbMPNN and both
  global/per-residue membrane models. Author residue numbering, insertion codes
  and long CIF chain IDs are mapped to the native runner input.
- Settings include checkpoint training noise, sampling temperature, optional
  seed, omitted amino acids (`empty` omits none), global/per-residue biases and
  omissions, fixed residues, homooligomer tying, ligand side-chain context,
  membrane labels and optional side-chain packing. Homooligomer copies must
  already exist; the code ties designed chains only and does not fabricate
  copies or redesign fixed targets.
- Automatic backend chooses Foundry for compatible CIF inputs with the core
  models, and shared LigandMPNN for other families and packing. Both receive
  the same applicable constraints; explicit backend selection is available.
- Optional Verify Sequences queues the maintained ESMFold prediction adapter
  after successful sequence design, one task per designed chain. Additional
  prediction engines remain configurable as downstream pipeline nodes.
- Existing-result sources are limited to completed successful jobs in the
  selected project; the HTTP API checks the source project too.
- The 9F6O example matches Tamarind: chain B, residues B26–32, B52–56 and B95–102
  (20 residues), two sequences. The second example is 1BC8, protein C conditioned
  on DNA A/B. Examples are copied into the user's input folder when loaded.

## Images and weights

- Existing `mn-ligandmpnn:latest`, source revision
  `26ec57ac976ade5379920dbd43c7f97a91cf82de`; no new image was required.
- Existing `mn-foundry:cu128`, Foundry 0.1.9.
- Core and membrane checkpoints already exist under
  `/mnt/db/reference_files/protpardelle-1c/model_params/LigandMPNN`.
- Added to `/mnt/db/reference_files/mpnn`:
  - `v48_002_epoch240_hyper.pt`
  - `v48_010_epoch300_hyper.pt`
  - `v48_020_epoch300_hyper.pt`
  - `v48_030_epoch300_hyper.pt`
  - `abmpnn.pt`
- `checkpoints.json` alongside them records URLs, source checksums and SHA256.
  Installer: `scripts/download_mpnn_checkpoints.py --reference-dir /mnt/db/reference_files`.
  HyperMPNN source is pinned to commit
  `1877ea181d23fe84fb33147748e1367e5a02b6a9`. AbMPNN comes from Zenodo 8164693
  (published MD5 `911df5d74f2955373286c3bb3df69409`).
- Each selected checkpoint is bound read-only into its container. No model is
  selected merely by a display label without a matching checkpoint.

## Verification already completed

- Full pytest: **418 passed**, one existing Starlette/httpx deprecation warning,
  15.53 seconds. Log `/tmp/mpnn-final-pytest.log`.
- New coverage includes all seven model requests, native checkpoint flags,
  target and motif protection, per-backbone constraint files, CIF chain and
  residue mapping, nucleotide context, shared page/node forms, all starter
  templates, project boundaries and one queued run reused by the worker.
- All five downloaded checkpoints loaded strictly into the installed
  ProteinMPNN model on CPU. This checked tensor/key compatibility only; no
  prediction forward pass or sequence sampling was run.
- Pipeline frontend rebuilt successfully.
- Live UI checked with the 9F6O example, without submitting a job.
- Preview saved as `docs/mpnn-input-preview.png`; selected chain tags, model
  buttons and optional controls use the green styling.
- Both local APIs serve `/api/v1/workflows/sequence.mpnn/schema` (28 parameters).

## Smoke run results

All jobs were submitted through the API and existing scheduler into the
**Generator smoke tests** project (1a4f7ce287d348d9af6a79bc80cbfeff). Each
successful run generated two sequences and was allocated both GPUs (0,1).

| Model/backend | Example | Run ID | Result |
| --- | --- | --- | --- |
| ProteinMPNN / shared | 9F6O, nanobody chain B, 20 CDR residues | 20261017-225504-3162bdce | 2 sequences; unchanged residues and target verified |
| ProteinMPNN / Foundry | 9F6O, same selection | 20261017-225546-f776fe2b | 2 sequences; unchanged residues and target verified |
| SolubleMPNN / shared | 9F6O, same selection | 20261017-225606-18ccd3f8 | 2 sequences; unchanged residues and target verified |
| HyperMPNN / shared | 9F6O, same selection | 20261017-225618-e4e3a496 | 2 sequences; unchanged residues and target verified |
| AbMPNN / shared | 9F6O, same selection | 20261017-225631-29f07514 | 2 sequences; unchanged residues and target verified |
| LigandMPNN / shared | 1BC8, protein C with DNA A/B context | 20261017-225643-fc5217af | 2 sequences; both DNA chains retained |
| LigandMPNN / Foundry | 1BC8, protein C with DNA A/B context | 20261017-225956-634e4489 | 2 sequences; both DNA chains retained and roles mapped correctly |

Two runtime issues surfaced and were fixed before the final reruns:

- Docker needs quotes around a multi-GPU selector. The shared GPU argument
  builder now passes the quoted device=0,1 selector so Docker no longer
  interprets the comma as a conflicting GPU count. One first ProteinMPNN
  submission failed from the old argument; the successful retry is listed above.
  Regression coverage was added.
- Foundry rewrites mmCIF auth chain IDs to label chain IDs. The result adapter
  now maps source IDs through the staged input CIF. A first Foundry LigandMPNN
  run produced sequences but mislabeled one DNA chain as a binder; the corrected
  run is listed above. Its output has DNA chains A/C and protein binder F, each
  with the expected residue counts.

Full pytest after these fixes: **421 passed**, with one existing
Starlette/httpx deprecation warning. The global/per-residue membrane modes were
not smoke-tested because the chosen examples are not membrane structures.
Optional side-chain packing and ESMFold sequence verification were not enabled.
No campaign state was changed; the V7 campaign was untouched. No jobs remain
queued or running.

## Services

- Repo `/home/user/programs/git-projects/mn-molecular-design-studio`.
- Python `/home/user/mambaforge/envs/mn-molecular-design-studio/bin/python`.
- Studio stays on port 8502.
- API port 8766 refreshed to PID 3103921, log `/tmp/studio-api-8766-mpnn.log`.
- API port 8767 refreshed to PID 3103990, log `/tmp/studio-api-8767-mpnn.log`.
- Scheduler PID 3525059 was not restarted; campaign states were not changed.
- The repository is entirely untracked in Git. Do not reset or remove prior
  work. Before-change copies of affected original files are in
  `/tmp/mpnn-change-before`.

## Follow-up tests, if requested

Use an appropriate membrane structure to smoke-test the global and per-residue
MembraneMPNN modes. Side-chain packing and downstream ESMFold verification were
left as separate optional checks.

Public request shape (use a managed uploaded/example path and a real project):

```json
{
  "schema_version": 1,
  "workflow": "sequence.mpnn",
  "project_id": "SELECT-A-PROJECT",
  "parameters": {
    "job_name": "mpnn-ab-example",
    "structure_path": "/path/to/managed/9F6O.pdb",
    "design_chains": "B",
    "designed_residues": "B26-32 B52-56 B95-102",
    "model_type": "ab_mpnn",
    "num_seq_per_target": 2,
    "sampling_temp": 0.1,
    "omit_aas": "C",
    "seed": 123,
    "backend": "auto",
    "verify_sequences": false
  },
  "resources": {"gpu_device": "all"}
}
```
