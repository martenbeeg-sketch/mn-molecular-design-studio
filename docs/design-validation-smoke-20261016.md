# Selected-design validation smoke · 16 October 2026

## Outcome

**216 tests passed**, with one dependency deprecation warning from Starlette's
httpx test client. The two Boltz-2 smoke predictions and their structure comparisons
completed successfully through the REST API and MN Compute Scheduler.

| Mode | Job | GPU | Total elapsed¹ | Binder Cα RMSD | Binder coverage |
| --- | --- | --- | --- | --- | --- |
| Complex with target template/MSA | 9257E | 0 | 154 s | 1.314 Å after target alignment | 94/94 |
| Binder alone, no template/MSA | 3E724 | 1 | 104 s | 1.251 Å after binder alignment | 94/94 |

¹Creation to completion, including scheduling, preparation and result analysis;
these are not inference-only timings. The complex target fit used all 115 target
Cα atoms and had RMSD 1.117 Å. The binder was not refitted in the complex comparison.

## Source selection

- Project: **BindCraft2 modality smoke · 2026-10-15**
- Project ID: `721f5cac3fac4a7f900ac874eb85a442`
- Source campaign: `hPDL1-miniprotein-serious-2026-10-15-pae-retry`
- Source job: **A75EF**, `20261015-143814-a75ef696`
- Candidate: `bindcraft2_00001`
- Recorded native filter result: `native_pass_filters: true`
- Source chains: target A, 115 residues; binder B, 94 residues.

Both smoke results were submitted to this same project. Preflight and submission
used `/api/v1/workflows/validate` and `/api/v1/workflows`, respectively. The engine
was `boltz2`, mode `off`, seed 20261016, one sample, three recycles and 200 diffusion
steps. No other engines were launched for this GPU smoke. Scheduler campaign state
was not changed.

## Input checks

The complex job staged one target-only CIF template containing chain A with 115
residues. Boltz's processed output includes a template NPZ for that chain. No binder
coordinates were supplied as a template.

The target MSA was fetched through the common ColabFold path: 5,317 rows were saved
in the shared `reference_files/msa_cache`, while the job used its requested 1,024
rows. The binder input contained its query sequence only. The binder-alone job
contained one 94-residue protein, no template and no MSA search.

## Issues found and fixed

1. **Native BindCraft2 CIF parsing:** the accepted CIF omits the optional occupancy
   column. Preflight initially rejected it, before creating any jobs. Design
   validation now uses the existing tolerant candidate parser. Both CIF and gzipped
   CIF cases are covered by regression tests.
2. **Reference absent from overlay:** the viewer's direct CIF loader skipped this
   native reference. Comparison overlays now use Cα coordinates from that same
   tolerant parser for both structures, preserving the already fitted coordinate
   frame. Both objects load in the browser after the fix.

Additional tests cover positive overall filter evidence, failed/unknown filter
exclusion, required project selection, cross-project rejection without queueing,
single job creation, target-only templates, binder query-only MSAs across MSA-capable
adapters, each supported native template serializer, and target-fit versus binder-fit
RMSD behavior. The other engine serializers were tested without launching their
GPU predictors.

## Saved results

- [Complex result folder](/mnt/data/RESULTS/mn-molecular-design-studio/workdir/runs/refolding-validation/20261016-234215-9257efd9)
- [Complex comparison CSV](/mnt/data/RESULTS/mn-molecular-design-studio/workdir/runs/refolding-validation/20261016-234215-9257efd9/artifacts/validation/comparisons.csv)
- [Binder result folder](/mnt/data/RESULTS/mn-molecular-design-studio/workdir/runs/refolding-validation/20261016-234216-3e724091)
- [Binder comparison CSV](/mnt/data/RESULTS/mn-molecular-design-studio/workdir/runs/refolding-validation/20261016-234216-3e724091/artifacts/validation/comparisons.csv)
- [Pytest output](/tmp/mds-design-validation-20261016/smoke/pytest-final.txt)
- [Exact API requests and smoke summaries](/tmp/mds-design-validation-20261016/smoke)
