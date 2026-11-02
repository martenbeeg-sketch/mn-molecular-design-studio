# LocalFold

Studio's **LocalFold** tile and pipeline node run the `sokrypton/alphafold3`
fork at `e68d628df1d815cccb48bb1118d8763a594bfa23`. The fresh image is
`mn-localfold:e68d628-cu13`, with Python 3.12, JAX 0.10.2 and CUDA 13.

Build it from the neighboring container repository:

```bash
../mn-tool-containers/build.sh localfold
```

LocalFold is an engine of `refolding.optimized_prediction`, so it uses the
existing project access checks, scheduler, molecule editor, collections,
pipeline inputs, candidate identities, result viewer, and artifact exports.
Select `engine: localfold` and set `engine_options.model` in API requests.
The earlier ColabFold preview workflow and historical AF3 image profiles keep
their existing versioned execution contracts.

## REST API

LocalFold uses the existing authenticated workflow endpoints:

- `GET /api/v1/workflows/refolding.optimized_prediction/schema` discovers the
  engine and its contract at `engine_contracts.localfold`. The contract lists
  every model, model-specific defaults, option types/ranges/choices, supported
  inputs, PAE availability, and terms requirements. `website_model` distinguishes
  the website's choices from the additional JAX models.
- `GET /api/v1/workflows/refolding.optimized_prediction/setup` includes the same
  schema and the LocalFold pipeline node.
- `POST /api/v1/workflows/validate` validates inputs without running inference.
- `POST /api/v1/workflows` queues the prediction through the existing scheduler.

Example request body for validation or submission:

```json
{
  "schema_version": 1,
  "project_id": "<existing-project-id>",
  "workflow": "refolding.optimized_prediction",
  "parameters": {
    "engine": "localfold",
    "job_name": "localfold-mol1",
    "protein_sequences": "ACDEFGHIKLMNPQRSTVWY",
    "samples": 1,
    "seed": 42,
    "engine_options": {
      "model": "mol1",
      "msa_mode": "automatic",
      "weights_precision": "int8"
    }
  }
}
```

Use `Authorization: Bearer <Studio session token>`. The top-level `project_id`
must identify a project the caller can access. Seeds, samples and recycles are
workflow parameters; model settings belong in `parameters.engine_options`.
For multiple independent complexes, pass `parameters.tasks` instead of a single
`protein_sequences` input. Modifications, bonds, crosslinks and conditioning use
the same task representation and model-specific validation as the Studio form.
AF3 requires `parameters.accept_alphafold3_terms: true`; HelixFold3 requires its
own acceptance in `engine_options`. Submit only after reviewing those terms.
The MCP workflow schema, validation and submission tools use this same API.

The Studio API keeps explicit LocalFold terms acknowledgments per Studio user
in `<Studio home>/localfold_terms_acceptance.json`. After a user accepts the
current AlphaFold 3 or HelixFold3 terms once, later API requests for that model
reuse the acknowledgment automatically. Records are tied to fingerprints of
the reviewed source text; changing the pinned terms manifest requires a fresh
acknowledgment. The model terms and their use restrictions still apply to each
prediction.

## Models and controls

The selector keeps localfold.org's 12 families together: AlphaFold 3, OpenBind-0,
OpenDDE, Boltz-2, Protenix-v2, Mol-1, IntelliFold-2, RoseTTAFold3, AF2-mono,
AF2-multi, ESMFold2 and BioEmu2. ESMFold2 has four variants. A separate
**Additional JAX models** group retains Chai-1, HelixFold3, OpenFold3 preview 2
and Protenix V1, all documented in the downloaded JAX fork. These are not
entries on the public website; together the two groups offer 19 parameter sets.

The entity editor follows the LocalFold per-chain layout. **Add entity** opens
the type picker above the rows, with separate CCD/SMILES ligands, modifications,
bonds, and crosslinks where supported. Modification rows target a polymer
residue; bonds target two named atoms; crosslinks target two protein residues
and use the fork’s crosslinker definitions. These rows do not consume chain IDs.
Protein, DNA and RNA modifications, covalent bonds, and crosslinks are written
to native AF3 JSON and validated before execution.

Menus follow the selected model and variant. AF2 accepts unmodified proteins.
ESMFold2 supports protein, DNA, RNA, ligand and modification inputs; the small
600M/300M LocalFold paths omit bonds, and every ESMFold2 variant omits crosslinks,
cyclic inputs and templates. BioEmu2 supports proteins, DNA, RNA, ligands,
modifications, bonds and contact conditioning. Contacts are also available for
Boltz-2 and Chai-1; pocket conditioning is available for those two models.
Unsupported rows retained during a model switch block submission with an error.

The selected model’s controls are available in the browser immediately. Form
changes rerun the prediction panel; schema discovery inspects only the chosen
workflow. The **My molecules** picker keeps the existing library selection scope
for refolding. It queries the database only while open; the closed picker does
no database work during form creation or model changes. Profiling the running
Studio form found 6.8 seconds spent in the closed picker out of a 7.3-second
render; deferring those queries reduced that server render to 0.044 seconds.
These numbers measure server execution, excluding browser and network time.
Advanced settings exclude the primary controls and options unused by
the selected model.

Each protein can use
the job alignment or a single sequence, and can opt out of uploaded templates.
The model selector and that model's active controls are shown together. The
controls select sample count, seed, recycles, diffusion steps, alignment source,
and MSA depth; additional model-specific controls remain under **Advanced
settings**.

**Automatic** alignment first reuses an exact-sequence A3M from
`<reference_files>/msa_cache` or the compatible legacy MSA repository. This
shared cache is populated by earlier predictions. On a cache miss, LocalFold
requests an alignment from the MMseqs2 server and saves the full result back to
the shared cache. **Reference files only** disables the server fallback.
Single-sequence mode skips MSA lookup and search.

The selected model sets the initial recycle, diffusion, and MSA-depth values.
These starting values and the preset choices follow LocalFold's web interface;
the Studio passes the selected values explicitly to the pinned JAX runner:

| Model | Recycles | Diffusion steps | MSA depth |
| --- | ---: | ---: | ---: |
| AlphaFold 3 | 3 | 25 | 128 |
| OpenBind-0 | 3 | 25 | 128 |
| Boltz-2 | 3 | 25 | 128 |
| Protenix V2 | 3 | 25 | 128 |
| Mol-1 | 3 | 25 | 128 |
| RoseTTAFold3 | 3 | 25 | 128 |
| IntelliFold 2 | 3 | 25 | 128 |
| OpenDDE | 3 | 16 | 128 |
| ESMFold2 · ESM-C 6B | 3 | 15 | 16384 |
| ESMFold2 Fast · ESM-C 6B | 3 | 15 | 16384 |
| ESMFold2 · ESM-C 600M | 3 | 15 | — |
| ESMFold2 · ESM-C 300M | 3 | 15 | — |
| BioEmu2 | 3 | 200 | 4092 |
| AlphaFold 2 pTM | 3 | — | 128:256 |
| AlphaFold 2 Multimer v3 | 3 | — | 128:256 |

An em dash means the control is not applicable for that model. AF2's sample
count is fixed at one; choose how many parameter sets to run (1–4 or all 5).
Its MSA depth is shown as clustered rows:extra rows.
Stored weight precision is `int8`, `fp16`, or `fp32`; it selects a weight blob,
not inference arithmetic. Official AF3 and AF2 use their original parameters.
The source's `auto` attention selection is the default; colabfold-kernels 0.4.0
provides its Pallas path. The image retains the documented Tokamax GPU capability
guard needed by the RTX 4090. The previous image's
`exact`, `fast`, and two-GPU `big` patches are not part of this new image.

BioEmu2 exposes structure and equilibrium modes, EDM and steered samplers,
Contact entity rows, potential YAML, and an optional initial mmCIF for propagation.
The API also accepts legacy contact JSON; it cannot be combined with Contact rows.
Each structure retains its physicality report. BioEmu2 has no confidence head:
pLDDT, PAE, pTM and ipTM are unavailable, including when a structure contains
zero B-factors. AF2 is protein-only and samples parameter sets rather than
diffusion samples.

The pinned fork writes a full confidence JSON beside each supported structure.
Studio attaches that native `*_confidences.json` as the candidate's PAE artifact,
so the existing result viewer can display its PAE matrix without copying the
matrix into a second file. BioEmu2 has no PAE artifact.

## Shared reference files

All LocalFold persistent data are in **`<reference_files>/localfold`**. This
machine uses `/mnt/db/reference_files/localfold`:

```text
localfold/
  weights/         converted blobs, separate model/precision directories; LM towers
  ccd/             sparse chemical component dictionaries
  runtime/<rev>/   JAX compilation and executable caches
  bioemu2/         steering reference tables
  home/ xdg/ huggingface/
```

Studio mounts the parent reference folder read-only at `/reference`, then
mounts its `localfold` child writable at `/reference/localfold`. This preserves
relative links on the host and in the container. Existing compatible converted
blobs and LM towers are linked from `mn-cofolding/cache/alphafold3/weights`;
original AF3 weights are linked from `alphafold3`. Original AF2 parameters are
linked from `alphafold_models` into `localfold/weights/alphafold2`. New downloads and unpacked data go into the new
LocalFold folder. The earlier files are retained in place.

The cached older `esmfold2_lm300m` and `esmfold2_lm600m` predictor blobs predate
the trained confidence heads in this source revision. They are deliberately
excluded from reuse; their matching LM towers are reusable. A first run fetches
the current converted weights into the new LocalFold folder.

Each job stores its input manifest, source/image identity, selected model,
precision, seeds, native command, structures, confidence sidecars where
available, logs, and physicality reports in its result folder. The manifest is
`artifacts/raw/optimized_prediction/mn_predict_manifest.json`; the native receipt
is `artifacts/raw/optimized_prediction/output/mn-batch-result.json`.

## Separate installations and input formats

LocalFold has its own image, Dockerfile and reference directory. Existing engine
tiles retain their own images and native adapters. Shared editor changes are
enabled by the LocalFold capability flag, and shared batch dispatch adds an
explicit `localfold` branch.

Studio first normalizes molecules into its internal task representation. For
LocalFold, it writes AF3-style JSON containing sequences, modifications, atom
bonds and supported conditioning constraints. The pinned fork then selects the
model-specific featurizer and runner. This common input envelope does not imply
that every model supports every field. Separate installations continue to use
their native YAML, JSON or sequence inputs through their own adapters.

## Validation on this workstation

The image was built and exercised on the RTX 4090. Short synthetic inputs passed
through the real Studio worker, container, output adapter, and candidate store:

| Model / mode | Inputs | Saved structures |
| --- | ---: | ---: |
| Mol-1, int8 | 2 | 2 |
| BioEmu2 structure, int8 | 2 | 2 |
| OpenFold3, existing linked int8 weights | 2 | 2 |
| BioEmu2 propagation, two steps | 1 | 2 |
| AF2 pTM, existing linked original parameters | 1 | 1 |

These installation checks used one recycle and two diffusion steps where
applicable. They check integration and artifact identity, not prediction
accuracy or capacity. All 19 models are checked at the Studio validation
boundary; the remaining models have not had a GPU inference run in this change.
Focused tests also cover precision selection, cache links, preflight errors,
confidence availability, batch identity, schema discovery, and pipeline wiring.
