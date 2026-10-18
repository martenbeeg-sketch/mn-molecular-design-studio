# Inverse-folding implementation handoff · 2026-10-17

The user requested independent ESM-IF1, Caliby and BoltzGen inverse-folding
tools, with Tamarind-inspired input forms and interchangeable pipeline nodes.
The user explicitly requested pytest and **no smoke tests** in this turn.
No model forward pass, sampling job, design job or GPU smoke was started.
Leave the paused V7 campaign and its scheduler state unchanged.

## Implemented

- Tiles/API workflows: `sequence.esm_if1`, `sequence.caliby`,
  `sequence.boltzgen_inverse`.
- `inverse_folding_ui.render_inputs` is shared by standalone pages and canvas
  node dialogs. Green controls, PDB/CIF selector, fetch/upload, chain and linked
  Mol* residue selection, examples, project selection, existing-result selection,
  multi-file/ZIP batch input, optional settings and scheduler submission.
- Tamarind ESM-IF1 uses a single designed chain, selected/fixed positions and a
  sequence count. The Studio adds temperature, seed and multichain backbone
  context. It uses upstream sampling so partial-sequence constraints work;
  the kit's CLI-only batched fast route is not exposed for this adapter.
- Caliby matches the two task choices (sequence design/scoring), five model
  pills, generated/provided ensembles, sequence/conformer counts, omissions and
  position constraints CSV. Distilled models use a single structure. Scoring
  accepts a single or supplied ensemble. Author IDs in CSV constraints are
  translated to native numbering. Native Potts energies are saved as
  `caliby_energy`; these are model scores, not binding affinity or confidence.
- Caliby optional acceleration: upstream; fast for single structure; exact for
  the generated 32-conformer protocol. Scoring runs upstream.
- BoltzGen uses the installed package's `inverse_fold_only.yaml` with a file
  entity and an explicit residue design mask. Its Predict task is instantiated
  directly. It loads the IF checkpoint and does not run backbone generation,
  refolding, affinity, analysis or filtering. The input may originate from any
  generator or prediction model. Temperature and residue exclusions are wired
  into `override.inverse_fold_args`.
- All adapters accept a direct PDB/mmCIF, structure batches or selected normalized
  candidates. Nodes accept Target, Backbones, Sequences or Predictions and emit
  Sequences. Connect any compatible folding node afterwards. Existing starter
  templates still use ProteinMPNN; replace or branch that node as desired.
- Chains/author residue IDs/insertion codes are mapped into job-owned native
  input files and restored in result structures. Connected targets and verified
  fixed motifs cannot be redesigned. No source artifact is modified.
- Outputs: FASTA, score CSV, native input/settings/output files, and normalized
  sequence-design candidates. Result structures retain input backbone geometry
  and nonprotein context. Redesigned side chains are removed; they are not new
  structure predictions. Previous sequence's prediction confidence is discarded.
- Project boundaries and source-collection containment are checked. API queue
  adoption reuses the scheduled run and does not create a generic extra job.

## Images and parameters

- Built `mn-esm-if1-opt:2.0.1-2b369911` from the local optimization kit, with
  fair-esm 2.0.1, PyTorch 2.4.0+cu121 and PyG. Image ID
  `sha256:2ad84d925267af57351282158df2d72b224b1fb3fe09ec563a8df22297271251`.
- Built `mn-caliby-opt:0.1-41d31560` from the pinned kit: Caliby commit
  `41d31560c3c73d7980d94f40f3c852b90bfab5c0`, atomworks-caliby 1.0.0,
  Protpardelle 1.3.2, Python 3.12 and PyTorch 2.6+cu124. Image ID
  `sha256:6b89991ea850dedec2afaf6c7fdd6ec2c0dd1068f40efc15e69ff95e3b73b809`.
  The kit's optional `_jitcache` COPY referred to a missing directory; the saved
  Dockerfile omits that optional cache import. Inference libraries are unchanged.
- Reuse `mn-boltzgen-opt:a07fea1-cu130-batch3`, code BoltzGen 0.3.2. Its normal
  batch entrypoint is overridden with Python for this independent adapter.
- Existing ESM-IF1 weights:
  `/mnt/db/reference_files/esm_models/esm_if1_gvp4_t16_142M_UR50.pt`.
- Downloaded Caliby assets to `/mnt/db/reference_files/caliby`, including all five
  sequence-design checkpoints, the Protpardelle cc95 ensemble model/config and
  ProteinMPNN files required by its ensemble environment. The shared HF snapshot
  also contains packing/large-model checkpoints; they were fetched but are not
  extra UI modes. SHA256/size/source revision are recorded in `checkpoints.json`.
  HF revision `a51f011f4ec7ffe2daad3ab9b2bfb67ff628a096`.
- BoltzGen IF weights and canonical-molecule archive already exist in
  `/mnt/db/reference_files/boltzgen-cache`. No duplicate model download.
- Download installer: `scripts/download_inverse_folding_checkpoints.py`.
  Image recipes: `scripts/inverse-folding/`. Their build script takes the local
  `uplifting-biomolecular-modeling` checkout as its build-context argument.
- Every inference container binds reference files read-only. Temporary caches
  and native outputs are in the new job's artifacts.

## Verification and services

Image import checks passed for ESM-IF1 and Caliby, and the BoltzGen Predict API and
inverse-only configuration loaded. These checks did not load a model for a
forward pass or run sequence design.

Pytest log: `/tmp/studio-inverse-full-pytest.log`. Focused checks:
`/tmp/studio-inverse-pytest.log`. New tests use fake native model APIs and cover
native masks/configuration, fixed targets/motifs, numbering/insertion codes,
nonprotein context, malformed outputs, Caliby model combinations, batches,
project boundaries, queue adoption, graph connections, forms and tile routing.
The final test count is recorded at the bottom of this handoff.

API ports 8766 and 8767 were refreshed to PIDs 3139489 and 3139553, with the same
environment/arguments as before. Logs:
`/tmp/studio-api-8766-inverse-folding.log` and
`/tmp/studio-api-8767-inverse-folding.log`. Both serve all three new schemas.
Scheduler PID 3525059 was not restarted. No campaign state was changed.

## Smoke test results · 2026-10-18

The user subsequently authorized running the smoke tests. Jobs used the
development API and the central scheduler, in the existing project **Generator
smoke tests · 2026-10-17**. GPU jobs were assigned to GPU 0. The V7 campaign
`overath-target-batch-api-current-plus-july-20261012-v7` was left paused and
untouched.

### Inverse-folding engines

- ESM-IF1 used the 3HTN example, chain A, residues 62–70, multichain context,
  and produced two sequences: `20261017-234506-4c665845`.
- Caliby single-structure design completed for all five checkpoints, two
  sequences each: `20261017-234538-2db8e46e` (Caliby),
  `20261017-234619-b4459e32` (SolubleCaliby),
  `20261017-234620-e13da633` (SolubleCaliby v1),
  `20261017-234620-59014e3c` (Caliby-distill), and
  `20261017-234620-b11e3f3e` (SolubleCaliby-distill).
- Caliby generated-ensemble design (8 conformers), native-sequence scoring,
  provided-ensemble design from all eight saved conformers, and fast mode all
  completed: `20261017-234831-963c20ad`, `20261017-234942-fdd540af`,
  `20261017-235139-10b8e928`, and `20261017-235217-6887716d`.
- Caliby's exact accelerated 32-conformer run matched the upstream 32-conformer
  run exactly for the sampled sequences and energy scores:
  `20261017-235300-8ecaa42d` (upstream) and
  `20261017-235433-1d5705ee` (exact).
- BoltzGen inverse folding on 1CRN and on a Protpardelle-generated backbone
  completed with two sequences per run: `20261017-235653-e2988b9e` and
  `20261017-235726-0b37ce92`.

### Full canvas handoff

Run `20261018-000031-f7d7be18` completed in 4m27s with 8 nodes and no failed or
cached steps. It used one selected 128-residue Protpardelle backbone and
branched it through ESM-IF1, Caliby (single-structure mode), and BoltzGen
inverse folding. Each generated sequence then flowed through ESMFold into a
single Results node. All three branches produced one normalized candidate with
a PDB and sequence. ESMFold pLDDT scores were 71.7, 81.6, and 79.4; all three
passed the configured minimum of 70. The corresponding ESMFold RMSDs to their
inverse-folding input structures were 1.89 Å, 1.11 Å, and 1.62 Å.

This is a workflow/connectivity smoke on a generated single-chain backbone. It
does not validate a binder-target interface or binding affinity. The optional
secondary-structure metric reports `pydssp unavailable`; it did not affect
sequence generation, folding, or the result candidates.

The first canvas attempt found a path-base mismatch: the request file lived in
the pipeline's `requests/` directory while managed paths had been stored
relative to the pipeline run directory. `_submit` now writes operational child
requests without rewriting those paths, so the workflow validator resolves
the selected upstream candidate files correctly. The failed first attempt
created no child model jobs; the successful rerun confirmed all branches,
folding handoffs, and Results collection.

Earlier standalone retries also exposed and fixed three adapter mismatches:
queued-job adoption now JSON-normalizes `Path` inputs; the BoltzGen image uses
`pytorch_lightning` rather than the absent `lightning` module; and Caliby can
reuse topology-matched backbone-only conformers by restoring the primary
sequence labels on staged copies (real sequence/topology mismatches still
fail validation). The failed attempts remain in history; all retry jobs passed.

No queued or running jobs remain in the smoke-test project. The one-shot failed
pipeline attempt is retained in job history for diagnosis. No generic
`optimized_prediction` job was created.

## Verification after the smoke fixes

Full suite: **471 passed**, one existing Starlette/httpx deprecation warning,
18.62 seconds. The first `pytest` executable invocation used the base Python
and could not import FastAPI; rerunning with the configured
`mn-molecular-design-studio` environment passed. Earlier implementation
verification was 470 passed before the additional pipeline request-path fix.

Existing ESM-IF1 checkpoint SHA256 was checked against the kit pin and matches.
There are 49 feature-specific tests. Live browser checks confirmed all three
custom input pages, example loading, Caliby's five model choices and
design/scoring switching, and reference-file readiness. A routing error found
during visual inspection was fixed and covered by three dispatcher tests.
Previews: `docs/esm-if1-input-preview.jpg`, `docs/caliby-input-preview.jpg`, and
`docs/boltzgen-inverse-input-preview.jpg`.
