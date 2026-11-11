# Production tools and pipeline acceptance — 2026-11-11

## Decision

**NO-GO for an unrestricted, all-tools production rollout.**

The concern about testing the actual tools was justified. The CPU suite passes,
and several real scientific workflows work, but fresh native executions expose
three concrete blockers: Boltz-2 output finalization, ESM-1v's container
dependency, and Genie3's incompatible installed image on PC004033. There is also
an observed pipeline-export limitation and substantial remaining live coverage.

This is a computational integration assessment, **not a certification of
scientific predictive accuracy, useful binder affinity or experimental success**.
It does not approve switching research production or importing important data
into the teaching pilot.

## Baseline and safety

- Tested committed source: `3daf0bb94193dc8e55d784d39d10f945b6c4a519`.
- Separate detached local worktree:
  `/home/user/.codex/worktrees/production-tool-acceptance/mn-molecular-design-studio`.
- Both `compute` workstations have separate detached test worktrees at the same
  commit. Their development checkouts remain on `d5df4d6616b624c69776645020d40307588ce6e3`.
- No application code, production authentication, existing scientific records,
  campaign membership or old-job status was changed. No image, model or Python
  package was downloaded or installed. The teaching instance was not repurposed.
- Live primary-host workers ran with a read-only host filesystem, separate
  writable test state, offline containers, installed images pinned to their
  observed immutable IDs and read-only reference mounts. Production model/cache
  mounts were not made writable. An ESM cache volume was replaced with a private
  test bind mount, not the existing shared Docker volume.
- GPU tests reserved only GPU 0 through a separate
  `mn-studio-production-acceptance` scheduler application. These **new diagnostic
  scheduler records** are the intentional shared-state additions; no existing
  research job was cancelled, modified or relabelled. All diagnostic leases were
  released. CPU-only binary probes allocated no GPU.
- Research services remain on `p11-25e003c/source` on all three machines. Their
  service PIDs and working directories were checked unchanged. Production
  authentication remains development mode, with zero real production accounts;
  a reviewed account/ownership transition is still required.

The test results concern the committed candidate, **not an independent
recertification of the older research deployment currently running**.

## What was checked

| Area | Fresh result | What it establishes |
|---|---|---|
| Workflow registry | 122/122 entrypoints import and are callable on all three hosts | Routing/import availability, not engine success |
| Canvas catalogue | 94 registered nodes inventoried | Includes input/control/collection nodes, not 94 independent engines |
| Pipeline factories | 47/48 validate topology | The deliberately unfinished blank canvas requires connections; not a defect |
| Optimized predictors | 18 profiles inventoried | Only Boltz-2 received fresh optimized-predictor inference here |
| Ligand manifest | 32 entries inventoried | Existing `validated` labels were not accepted as new evidence |
| MPNN checkpoints | All 19 registered noise/model checkpoint paths exist on each host | Presence, not all-variant inference or cross-host byte equivalence |
| Full isolated CPU suite | **3,099 passed, 55 skipped, 103 warnings**, 905.22 s | Broad contracts, mocks, UI regressions and lifecycle tests |
| Installed image-side adapter fixtures | **70 passed**, 80.41 s | Actual installed adapter source with controlled model/process fixtures; no new inference |
| Initial real GPU batch | **4 passed, 8 failed, 1 skipped**, 1,186.04 s | Seven failures share the Boltz defect; one was a harness cache-volume refusal |
| ESMFold2 VHH follow-up | Native generation passes; downstream Boltz fails | The harness refusal was corrected without touching the shared cache |
| All 12 QC/interaction tools | **11 return valid measurements; ESM-1v fails** | One or two public-fixture items per tool, not a validation benchmark |
| Vina → PLIP → PandaMap → PoseBusters | All four real workers succeed | One prepared 1IEP/imatinib example, native Vina exit 0, one scored compound |
| Actual CPU graph controller | Completes QC → cohort clustering → Results | Two candidates, retained method-qualified measurements, one shared cluster and representative |
| Parent graph ZIP | Builds, but contains no coordinate files | Not accepted as a self-contained scientific pipeline archive |

The CPU sweep's 18 additional skips relative to the earlier readiness report
come from image-side source fixtures that are absent beside a fresh worktree.
The installed Boltz/Protenix adapter sources were then copied from newly created,
never-started inspection containers, and the 70-test fixture file passed,
including those 18 cases. Do not add 70 to 3,099 as if they were all unique tests.
Remaining skips include explicitly opt-in inference/container checks and optional
dependencies; they are not passes.

`PRODUCTION_TOOL_MATRIX.csv` contains **273 rows**: 122 workflows, 94 canvas
nodes, 18 optimized-predictor variants, 32 ligand-manifest tools and seven MPNN
models. These families overlap. It is not a claim of 273 distinct tools tested.
Unexecuted methods and unresolved image mappings remain explicitly labelled.

## Release-blocking findings

### PROD-001 — Boltz inference succeeds, then the application rejects its output

**Severity:** P1. **Class:** reproduced application/adapter integration defect.

Eight fresh Boltz jobs produced **11 native prediction outputs**. Every native
batch reports `execution.return_code = 0` and completed outputs. Nevertheless,
application finalization marks the job failed because the native batch adapter
leaves an input-staging symlink such as:

```text
artifacts/raw/optimized_prediction/output/.mn-predict-boltz2-inputs/<item>.yaml
  -> /work/artifacts/raw/optimized_prediction/inputs/<item>.yaml
```

The canonical input YAML is retained inside the run. `/work` is a container
namespace, so the alias is dangling on the host. The strict collector classifies
anything under an `output` directory as scientific output and rejects this
input alias. The required PDB/PAE outputs are present; suppressing integrity
checks or marking missing structures successful would be the wrong fix.

**Concrete first run:** `20261111-222049-d41a89d7`.
Its native batch completed one prediction. The application error is:
`Required scientific output is not accessible in persistent Results: Dangling scientific output is not declared or materialized as a canonical Results artifact`.

**Affected fresh paths:** sampling/recycling controls; reference prediction to
QC/selection; prediction to ipSAE/Rosetta; real RFD4 generation to Boltz; mixed
protein/ATP/zinc prediction; RFD4 two-stage sequence variants to Boltz; and native
ESMFold2 VHH generation to Boltz. Downstream analyses blocked before invocation
are **not** counted as passing analyses.

**Independent isolation:** in a separate diagnostic copy, materializing only
the exact canonical input alias changes collector validation from incomplete to
complete. All 26 original regular files remain byte-identical, including the
actual prediction. The original run remains failed and was not repaired or
promoted. See `boltz-finalization-reproduction.json` in the evidence archive.

**Source:** `engine/workflow_engine/core/artifact_collection.py:402` and `:418`;
`engine/workflow_engine/workflows/optimized_prediction.py:3242`, under the
`mn_molecular_design_studio/` package. Relevant lines have history in
`be6bcc49` and `2f20f9b1`; this test alone does not establish a Phase-2 regression.

**Recommendation:** normalize/materialize the adapter's verified native input
aliases before finalization, or declare an exact verified execution-only alias.
Retain strict checks for genuinely missing structures, PAE and scientific
outputs. Add the native directory layout to CPU finalization regressions and
rerun fresh Boltz plus the dependent pipelines. A completed GPU process is not a
completed Studio pipeline.

### PROD-002 — ESM-1v scoring cannot run in its registered image

**Severity:** P1 for ESM-1v-dependent workflows. **Class:** reproduced runtime
dependency mismatch, not a missing model download.

Both fresh ESM-1v items fail with `No module named 'Bio'`; no probability
measurements are returned. All five registered ESM-1v checkpoint files are
present. The registered image is `mn-esm-if1-opt:2.0.1-2b369911`, observed ID:

`sha256:2ad84d925267af57351282158df2d72b224b1fb3fe09ec563a8df22297271251`.

A separate CPU-only dependency probe confirms `Bio_available: false` in that
exact image. The ESM-1v branch invokes `sequence_metrics()` for sequence
validation, which depends on Biopython. The same image successfully produces
real ESM-IF1 likelihoods, demonstrating that image presence and a neighbouring
tool's success do not establish ESM-1v compatibility.

**Source:** `engine/workflow_engine/workflows/quality_runner.py:345`, including
the `sequence_metrics()` call, and `quality_tools.py:11` in the package.

**Recommendation:** either ship the required dependency in a newly pinned image
or separate pure sequence validation from composition calculations. Preserve the
amino-acid and 1,022-residue protocol checks, five-model averaging and failed-item
semantics. Installing Biopython only in the host Conda environment will not fix
this container. Rerun real five-model scoring and verify finite probabilities.

### PROD-003 — PC004033's Genie3 image cannot execute on its RTX 5090

**Severity:** P1 for dispatch to that engine/host. **Class:** deployment defect.

PC004033 has an RTX 5090, compute capability 12.0. Its `mn-genie3:cu128` image is:

`sha256:23ff4ee8da51d25b6d5df97e4017a9134c7d660a981fd3a113621fe6f63f1bf9`.

The actual PyTorch runtime reports CUDA **12.6**, with targets only through
`sm_90`, despite the image's `cu128` tag. A scheduler-leased, offline CUDA tensor
operation fails with **`no kernel image is available for execution on the device`**.
No model inference or checkpoint download was needed to reproduce the failure.
The compiled-target probe independently rejects this image/GPU combination.

**Recommendation:** provision and pin a compatible image, prove native kernel
and model execution on the target host, and gate eligibility by the actual
target-host runtime. Do not repair this by changing the host driver, trusting a
tag, relabelling the job, or copying an arbitrary 4090-only build. Other
same-named images differ between hosts; different IDs alone are not proof of a
defect or of scientific non-equivalence.

## Other operational observations

### Fleet coverage is not uniform

| Host | GPU | Available source/catalog-referenced image tags |
|---|---|---|
| PC002975 | Two RTX 4090, 24 GiB each | 86/86 |
| PC004032 | RTX 5090, 32 GiB | 42/86 |
| PC004033 | RTX 5090, 32 GiB | 44/86 |

The union includes optional/compatibility aliases, so not all 86 images are
requirements for every workflow. The matrix lists exact declared requirements
where known. A concrete retained gap is `mn-protein-qc:1.0.0` on PC004032; that
host cannot execute the associated QC worker as configured. This extends the
previously documented missing-image risk, rather than treating it as a newly
discovered source defect. No images were silently substituted or synchronized.

Compiled-runtime probes inspect 39 referenced images. Some use JAX, nonstandard
Python entrypoints or CPU-only binaries; missing Torch inventories are **unknown**,
not compatibility passes or proof that those engines are broken. In particular,
the generic BindCraft2 probe cannot find `python` at its default entrypoint.
Fresh BindCraft2 generation was not performed here.

### Parent run ZIP is not a complete pipeline export

The real CPU graph preserves upstream Sequence Composition measurements through
method-qualified lookup, even after clustering changes the immediate flat
metric map. Its two identical inputs share one cluster and one representative;
the graph and both native child jobs complete successfully.

Coordinates remain accessible from the retained sibling child run. However,
`build_run_zip(parent)` exports the parent's local artifacts, not those sibling
coordinates: no PDB/CIF is present in the generated ZIP. Treat this as an
**observed export/portability limitation**, not evidence that the live coordinates
were lost or that the graph failed computationally. The complete-pipeline
archive acceptance check is not cleared. Keep the related child jobs and the
complete Results recovery set; do not use this ZIP alone as a research backup.

### Receptor preparation matters

The first Vina probe supplied the unprepared crystal complex. Meeko attempted
an external CCD download for retained `STI`, which the offline sandbox correctly
blocked. That setup failure is not counted as a Vina engine defect. With an
explicit chain-A, protein-only receptor and its selected ligand-centred box,
real docking succeeds: one scored compound, zero failures, followed by valid
PLIP, PandaMap and PoseBusters measurements. The deliberately small docking
search is an integration check; its score is not experimental affinity evidence.

## Coverage limits and next gates

1. **Fix PROD-001 and PROD-002; replace/quarantine PROD-003's incompatible image.**
   Preserve the failing evidence. Do not manufacture successful historical jobs
   or weaken canonical-output validation to turn the suite green.
2. Repeat one fresh Boltz run and all affected end-to-end chains, then run them
   through the authenticated application API and representative browser forms,
   ingestion, reuse, review, inspection and downloads. Most GPU tests here call
   the real workflow workers directly, not a full authenticated UI journey.
3. Finish the explicit matrix. Seventeen other optimized predictor profiles,
   most native generators, nondefault MPNN/noise variants, MD/free-energy/QM,
   ADMET and many ligand generators have **not** received fresh model/native
   execution. Present images, checkpoint paths and old reports are insufficient.
4. Test approved cross-workstation dispatch, actual target-host image/weight
   provenance, cancellation, worker restart and recovery at the intended scale.
   This run does not certify fresh distributed multi-stage inference, sustained
   throughput, long campaigns or numerical equivalence across 4090/5090 builds.
5. Complete the previously documented coordinated scientific-file/database/
   queue/replica recovery and production account/ownership transition, then
   activate a pinned release on all writers in a maintenance window. Existing
   old-job and cleanup candidates remain preserved; no cleanup gate is cleared
   by these new disposable tests.

Input/reference fixtures are deliberately small and public. No wet-lab,
prospective predictive-validity, chemically exhaustive or production-load
benchmark was performed. A valid negative scientific assessment is not an engine
failure; conversely, a scheduler/controller's exit zero is not proof that every
item produced a valid measurement.

## Evidence

- `PRODUCTION_TOOL_MATRIX.csv`: complete inventory with explicit live/unknown
  status and per-host exact-image availability where resolved.
- `PRODUCTION_TOOL_ACCEPTANCE_EVIDENCE.tar.gz`: sanitized inventories, CUDA
  profiles, JUnit/logs, tool/graph receipts and reproducible diagnostic scripts.
  It excludes databases, credentials, lease tokens, full image environment
  dumps, model weights and Docker images.
  SHA-256: `0cfd8fa18577b46f1a515d69616921842af4fa655c3373f1170998707ca56f63`.
- Complete raw scientific-run evidence is retained privately under
  `/home/user/programs/deployments/mn-molecular-design-studio/production-acceptance-3daf0bb94193/`.
  It is a test-evidence snapshot, not a production database or restorable live
  deployment. Original scratch runs remain under
  `/tmp/mds-production-acceptance-20261111/` without repaired failure statuses.
  The private `complete-scientific-runs.tar.gz` SHA-256 is
  `69ef8de8fe02f2d7036d3a02044b0de684bdb456944806fcdce495b693bed286`.

The original independent audit and prior readiness evidence remain unchanged.
This report does not turn the earlier narrow teaching-pilot decision into an
all-tools production certification.
