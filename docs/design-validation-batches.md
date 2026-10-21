# Refolding selected designs

Open a cofolding engine tile, select a project, then choose **Input source → Validate
existing designs**. Select design campaigns and the individual designs to predict.
Submission is manual and uses the existing API and compute scheduler.
The full selection is automatically one native batch job; no batch-mode switch is required.

## Eligibility and project scope

- The project is required before campaigns can be listed or submitted.
- Only `design` jobs with status `completed`/`succeeded` and `result.success == true`
  are eligible. A successful container exit alone is insufficient.
- A campaign appears only if at least one normalized candidate has an affirmative
  overall `native_pass_filters` or `pass_filters` flag. An explicit false flag vetoes
  a pass. Boolean values, `true`/`false` strings and 1/0 are recognized.
- Legacy BindCraft candidates explicitly recorded in its native `Accepted` pool
  are also accepted. Generated/prefilter pools, failed filters, and unknown filter
  status are excluded. The app does not invent quality thresholds for tools that
  did not record a filter decision.
- Only eligible designs appear in the design multiselect. Eligibility is checked
  again during API preflight and when the scheduled worker prepares its inputs.
- The API also rejects campaigns outside the selected project, including campaigns
  the user could otherwise access. Source job access checks remain in force.

## Inputs and prediction modes

**Complex refolding** uses the designed target and binder protein sequences. Target
MSAs are optional: existing reference alignments are reused, missing alignments are
fetched through the shared ColabFold path, and the job retains its own A3M files.
The binder is always a single sequence, with no MSA search or reference MSA reuse.
If selected and supported, a template contains only target coordinates.

**Binder alone** predicts the binder sequence without target coordinates or an MSA.
The same engine settings and inference profiles remain available.

AlphaFold2 Initial Guess requires a target structure and currently supports one
designed binder chain per complex. BoltzGen Fold requires one protein target and
one binder; neither offers binder-alone mode. Other engines accept multiple protein
target/binder chains where their normal input adapter supports them. Complex
comparison currently requires a protein target; ligand-only target designs can use
binder-alone prediction. Identical target/binder sequences cannot have different
MSA policies in adapters that share homomer alignments (Boltz, ColabFold and Chai);
disable target MSA for that case.

## Optional structure templates

The shared editor places **Templates** in optional settings beside the MSA controls.
PDB/mmCIF uploads are matched to protein chains by sequence, or restricted to a
selected input chain. Template files and native converted inputs stay with the job.

| Engine | Installed input route |
| --- | --- |
| AlphaFold3 JAX / ColabFold 2.0 preview | Inline mmCIF with residue index mappings |
| AlphaFold3 JAX with OpenFold3 P2 parameters | Inline mmCIF with residue index mappings |
| AlphaFold3 PyTorch with OpenFold3 P2 parameters | Inline mmCIF with residue index mappings |
| Boltz-2 | CIF template and explicit query/template chain mapping |
| ColabFold 1.6.1 | Per-item custom template PDB |
| Protenix V2 | `templatesPath` pointing to mapped inline-CIF JSON |
| OpenFold3 OpenBind-0 | `template_cif_paths` and template chain IDs |
| RosettaFold3 | Structure components and `template_selection`; exact sequence match required |
| Chai-1 | Local CIF files and explicit template-hit table |
| BoltzGen Fold | Required target PDB input |
| AlphaFold2 Initial Guess | Required target PDB input, target-only template scope |

The pinned Protenix V1 and OpenDDE adapters have HHR/A3M template-search routes,
but no verified direct PDB/mmCIF upload serializer. They retain their existing
search controls without offering the new structure upload. ESMFold, ESMFold2
Fast/Full and AtlasFold do not expose supplied structure templates in their installed
adapters.

Chai's installed batch adapter needs a job-local bridge for local template files;
the AF3 preview batch adapter needs its existing runner model families exposed.
Guarded scripts apply these changes only inside each disposable job container.
They do not change the base image, model weights, or prediction parameters. If a
pinned adapter changes, the bridge fails explicitly instead of silently ignoring
the requested input.

## Comparison in results

Every predicted sample is compared against a job-local copy of its selected source
design. Protein chains and Cα atoms are matched by sequence alignment.

- **Complex:** fit all matched target Cα atoms, apply that transform to the whole
  prediction, then calculate binder Cα RMSD without fitting the binder again.
- **Binder alone:** fit matched binder Cα atoms and report fitted binder RMSD.

The results show the RMSD in Å, compared atom count, coordinate coverage, target fit
quality where applicable, and a grey-reference/teal-prediction overlay. At least
three non-collinear matched atoms are required for fitting. Missing or ambiguous
chain mappings produce an explicit comparison failure. Prediction files are
retained if comparison fails; the job records the failed analysis.

Download the aligned prediction or the batch-wide `artifacts/validation/comparisons.csv`
report. Source run/candidate IDs are retained for every comparison. These are result
analyses, not prediction filters or additional setup controls.

## API example

POST `/api/v1/workflows/validate` before submitting the same payload to
`/api/v1/workflows`:

```json
{
  "schema_version": 1,
  "workflow": "refolding.optimized_prediction",
  "project_id": "PROJECT_ID",
  "parameters": {
    "engine": "boltz2",
    "job_name": "selected-design-refolding",
    "optimization_mode": "off",
    "samples": 5,
    "validation_mode": "complex",
    "target_use_msa": true,
    "use_target_template": true,
    "validation_sources": [
      {
        "source_run_dir": "/managed/results/runs/design/DESIGN_RUN_ID",
        "candidate_ids": ["bindcraft2_00001"]
      }
    ]
  }
}
```

For AlphaFold2 Initial Guess, use `refolding.af2_design_validation` with the same
`validation_sources` plus `num_recycles`, `model_count`, `multimer`,
`target_use_msa`, `msa_max_depth`, `cpu_cores` and `gpu_device`.

Verification: **216 pytest tests passed**, including the new eligibility, project,
MSA, template serialization, native CIF and RMSD/overlay cases. Both Boltz-2 GPU
smoke modes completed successfully using accepted hPDL1 design A75EF. See the
[smoke report](design-validation-smoke-20261016.md) for inputs, job IDs and results.
