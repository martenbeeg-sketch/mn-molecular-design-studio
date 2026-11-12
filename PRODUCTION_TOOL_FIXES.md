# Production tool corrections — 2026-11-12

Subsequent deployment status: see `PRODUCTION_RELEASE.md`. The historical native
test results below are retained; the corrected research release is now active.

## Scope and release decision

This follow-up addresses the three concrete release blockers and the final-result
export limitation reproduced in `PRODUCTION_TOOL_ACCEPTANCE.md`. That report and
its failure evidence remain historical baseline evidence, not overwritten passes.
The original independent audit worktree and the committed acceptance worktree
remain unchanged. Corrections start from `fd036f81fca0a0ea7a24ac706edd688aff823b16`.
Application correction commit: `b333f3a2b36f5083d81ed53b001300e391efe320`.
Export hardening follow-up / final application candidate:
`86804ff29f50d9846ba285d8d5b7f502c6d6dbdd`.

**A corrected candidate is not blanket approval for all-tools production.**
Unexecuted methods in `PRODUCTION_TOOL_MATRIX.csv`, the remaining review-filter
limitation, research account/ownership migration and coordinated recovery rehearsal
remain separate gates. No prediction result establishes experimental efficacy.

## Corrections

### PROD-001: retain Boltz's exact native input alias

The optimized Boltz adapter now materializes only the exact input aliases created
by the installed batch adapter, after native execution and before collection.
Each alias must match that item's submitted `/work/.../inputs/<id>.yaml` path.
The canonical input must remain a regular file inside the same run and match its
pre-execution SHA-256 identity. Redirected directories, traversal, changed inputs
and unexpected targets fail closed. Already materialized matching inputs are
idempotent. Unrecognized aliases and genuinely missing scientific outputs remain
subject to the unchanged strict collector.

This does not repair or promote any historical failed job. Native coordinates and
PAE files are not altered. CPU regressions reproduce the installed adapter layout
and exercise negative integrity cases in `tests/test_boltz_input_aliases.py`.

### PROD-002: remove composition-library coupling from sequence validation

Pure standard-amino-acid validation no longer imports Biopython. ESM-1v uses this
validation directly, normalizes case consistently and retains its 1,022-residue
limit, five separate checkpoints and unmasked wild-type marginal protocol.
Pre-existing item errors and invalid/nonfinite native probabilities cannot become
successful scores. Sequence-composition metrics still use Biopython; no fake
composition or replacement probability is introduced. Protein-Sol's validation
also no longer calculates unused Biopython properties.

`tests/test_esm1v_protocol.py` blocks Biopython imports and exercises all five
checkpoint calls, averaging, invalid sequences, item errors and invalid native
probabilities with controlled model fixtures. Those tests are not inference.

Fresh real scoring in the previously failing installed ESM image completes both
public 7EOW items using all five installed ESM-1v checkpoints. Each has native
probability `0.8726456580079687`; both ESM-IF1 controls also remain valid. No image,
host package or checkpoint modification was required for this correction.

### PROD-003: separate the compatible Genie3 runtime from ambiguous old tags

The application catalogue and both Genie3 manifests select
`mn-genie3:cu128-sm120-7549829b`, instead of the divergent `mn-genie3:cu128` tag.
Its installed source image from PC004032 is:

`sha256:7549829b69548c7b45ecb3c251879f1595fe75fdb4ecc83dca42958c495fc1df`.

The primary's containerd image store reports the imported image as
`sha256:208571a8d0814d894b7248156b8557b2f9bb7daa0dc89213e69a9cc43f6599e6`.
All three hosts have identical executable-content fingerprints:
`sha256:bc2e4ae37d2aa4038bc20b4a72a704da5bb37833af5d2ead16632f3ed15a40df`.
The different store-specific IDs are not a different runtime.

Its actual PyTorch is `2.7.1+cu128`, CUDA 12.8, with compiled targets
`sm_75 sm_80 sm_86 sm_90 sm_100 sm_120 compute_120`. The older PC004033 image
remains available under its original tag for reproducible history and rollback;
the host driver is unchanged. The new image is transferred by immutable ID and
tagged separately, never by blindly overwriting a shared production tag.

The installed checkpoint and configuration also match byte-for-byte on all three
hosts: checkpoint SHA-256
`57039afcd4e5a66e58c52ba1ece256bc311a4294a31c37fca7301dd3a4c4479e`,
configuration SHA-256
`a01d7ad8b908e4380f45a11c8d691ccf27af3f50a9dfb209823ee99989a343e5`.

The catalogue's historical Genie3 source-version label is not independent proof
of an upstream Git commit inside a stripped image. Installed image content and
checkpoint identity, not that label, identify these executions.

### Final-candidate scientific ZIP

The UI's job ZIP now reauthorizes the parent run and each referenced producing
job, requires catalogue membership for linked artifacts and enforces recorded
artifact hashes. An additional `scientific_bundle/` contains final candidate
coordinates, preserved measurements, deduplicated SHA-256-identified structures
and archive-root-relative canonical references. Source files are not rewritten.

Graph snapshots/execution receipts and artifact manifests are also retained.
Arbitrary outward, directory or dangling symlinks are not followed; their names
are recorded in `export_manifest.json`. The same linked-artifact authorization
helper is shared with inspection, rather than importing UI code into the engine.

Adversarial review also reproduced an export confidentiality gap with an aliased
top-level `artifacts` directory. The final follow-up skips that directory rather
than walking into arbitrary host files. A missing/dangling candidate collection
cannot be advertised as a successful scientific bundle. These cases use only
private sentinel files in regressions; no production/private file was exported.

The default programmatic `build_run_zip(run_dir)` remains a local run archive;
cross-run inclusion requires the explicit authenticated `user` argument, which
the UI supplies. The scientific bundle covers final candidates, **not every raw
intermediate child-run file, trajectory or model checkpoint**. Original run JSON
is retained as provenance and may still contain its original managed references;
portable final records live in `scientific_bundle/candidates.jsonl`.

Fresh actual CPU orchestration of reference input → SequenceQC → cohort
SequenceClustering → Results completes. The authenticated export retains two
candidate records, their method-qualified upstream measurements and one deduplicated
coordinate file with a verified SHA-256 identity. All original runs remain intact.

## Validation ledger

- Final focused Boltz/ESM/export/quality/collector/input suite: **112 passed,
  10 opt-in container skips**, 17.55 seconds. An initial test assertion used the wrong collector status
  key; corrected before this successful run, without changing collector behavior.
- Export/controller/UI regression sweep: **89 passed, one optional `chaperone`
  dependency skip**, 955.85 seconds.
- Installed image-side native adapter fixtures: **70 passed**, 91.52 seconds.
  The first broad sweep exposed 12 old Boltz fixture failures because that harness
  translated the container's input paths into host paths but retained its resulting
  host-target symlinks. The fixture now translates those exact aliases back to
  `/work/...`, matching actual container output. The application does not whitelist
  arbitrary host aliases. A clean full sweep is rerunning with this correction.
- Full final-source CPU suite in two disjoint test-file shards: **3,154 passed,
  37 skipped, 103 warnings**. Core shard: 3,110 passed / 37 skipped, 593.59 seconds;
  workflow-mode shard: 44 passed, 559.06 seconds. The union matches all **3,189
  collected test cases**, plus two module-level optional-dependency skips.
  No test IDs are duplicated between the shards. These are completed shards,
  not a claimed successful monolithic run. Earlier monolithic attempts were
  interrupted or time-limited, including a 30-minute GPU batch. An extended CPU
  attempt reached 3,053 reported passes and 37 skips
  before interruption inside OpenMM. Partial JUnit records can include an
  unfinished final case and are not used to fabricate passing totals.
  The primary checkout's default fixture discovery finds the neighbouring
  `mn-tool-containers/shared/batch_api` source; the older frozen worktree did not.
  This accounts for 18 fewer skips than the initial acceptance baseline. That
  neighbouring source is not assumed to be the installed image: the separate
  70-case run explicitly uses adapter source copied from the installed images.
  Those overlapping cases are not added again to the unique full-suite total.
- Final export-hardening/inspection/review/quality regression slice: **195 passed,
  10 opt-in skips**, 60.63 seconds; committed final ZIP adversarial file rerun:
  **11 passed**, 0.12 seconds. The later ZIP-only follow-up is not misrepresented
  as having been present at the start of the baseline broad run.
- Fresh GPU prediction/generation/downstream cases: **12 unique cases pass**,
  one retained-artifact replay test is intentionally skipped. The first bounded
  batch reports nine completed passes before timeout interrupts a later case;
  the separate frozen-final-source tail batch completes **four passes**, 645.32
  seconds, including one repeated case. Counts are deduplicated, not 13 passes.
  The union covers both AF2 conditioning modes, Boltz sampling/recycling,
  prediction → QC/selection, ipSAE/Rosetta, real RFD4 generation, mixed ATP/zinc
  inputs, both two-stage sequence-design routes, AbMPNN framework preservation
  and native ESMFold2 VHH → Boltz. The timed-out attempt remains failed; reruns
  create new private jobs rather than relabelling it.
- PC004033 compatible Genie3: CUDA tensor probe **passes**; one fresh native
  unconditional 40-residue sample **passes**, native exit 0. The retained PDB has
  40 C-alpha coordinates, all finite, SHA-256
  `f1d63fbf8fcfeff3b17437c4704997d24d4133971b6ac0ef5e4cbf7725b074f7`.
  This tests a backbone generator, not sequence design or binder efficacy.
- Primary-host final-source Genie3 application worker: **passes** one fresh native
  generation through normalization and strict finalization. It retains one
  40-residue C-alpha backbone, with `binder_sequence: null`, not a fabricated
  designed sequence. All three required artifacts validate. A generator result
  does not supply an assessment `scientific_status` or imply binding efficacy.
- PC004032 final-source isolated CPU release/peer/review/inspection slice:
  **212 passed, 10 opt-in container skips**, 5.69 seconds.
- PC004033 same final-source slice: **212 passed, 10 skips**, 5.38 seconds.
- All seven CPU QC tools rerun on frozen final source: **all seven return valid
  measurements** for both public-fixture items, including changed Protein-Sol
  validation. Together with the two real ESM tools, nine protein QC methods have
  fresh successful correction-run evidence. The three ligand QC workers retain
  their earlier baseline acceptance; they are not falsely counted as rerun here.

The isolated real QC helper records item-level outcomes; its process exit alone
is not a scientific success criterion. The export harness initially inspected a
method-qualified measurement as a flat metric key; the corrected harness uses the
application's `metric_value` contract, and its fresh graph/export run passes.

Early GPU/QC checks exercise the correction working tree. The final GPU tail,
Genie3 worker and graph checks use the separate extracted `86804ff` source;
the final CPU shards run the committed checkout with application/test sources
verified unchanged before and after execution. Both compute slices mount the
staged final source read-only. These are distinguished rather than claiming
every development-stage attempt used one immutable commit.

Genie3's first sandbox probes exposed harness permissions rather than another
CUDA failure: capability-dropped root on a rootful daemon does not bypass a
normal user's directory permissions, and native logging also needs a writable
location. The successful target-host check runs as the owning UID/GID, with only
its private output and two private log directories writable. An earlier slow
startup probe was intentionally stopped; its failed/interrupted evidence remains
retained. No production permissions, image contents or reference caches were
relaxed to obtain the successful native run.

## Preserved boundaries and remaining work

The final committed candidate is staged, not activated, at
`~/programs/deployments/mn-molecular-design-studio/tools-86804ff29f50/source`
on the primary and both `compute` PCs. The compressed tracked-source archive is
`8db889f998e93cbe906e56713374f92cd50a10da9a7284bcf66b31940d185741`
(SHA-256), verified after transfer. Peer development branches remain unchanged;
the commit is imported into their object databases/FETCH_HEAD only.

Evidence is retained in `PRODUCTION_TOOL_FIXES_EVIDENCE.tar.gz` (66 files,
495,138 bytes), SHA-256
`06f025805058d3add50cf8b5d8ce41f03eb704648a262a16df4b551dbbf8e7c2`.
It includes exact-case CPU coverage reconciliation, JUnit results, native leases,
failed/limited attempt logs, image/checkpoint identities, harnesses and the actual
authorized final-candidate export. No authentication databases, credentials,
private CA keys, model weights or complete image-inspection environment dumps
are included.

Raw application runs and the successful target-host native run are retained
separately at `tools-86804ff29f50/acceptance/native-runs.tar.gz` under the deployment
directory (90,045,662 bytes), SHA-256
`58efeb71388a1908f4ba0fac675be68f196cefc950c4da65be60274b08cdbe59`.
This is test-evidence retention, not a coordinated research-production recovery
backup. Earlier target-host sandbox attempts also remain in their separate
`compute` scratch directories.

- Research services and authentication are not switched, and the isolated teaching
  pilot is not repurposed. Staging source/images is not service activation.
  The running research deployment remains `p11-25e003c/source`, not this corrected
  candidate. Users must not assume its existing UI already runs these fixes.
- No important campaigns/imports, old job statuses, production databases or model
  caches are cleaned or rewritten. Smoke/testing names and single-entry campaigns
  are not safe deletion criteria.
- New versioned image copies, normal agent image-inventory updates and our own
  diagnostic scheduler records are intentional operational additions. No existing
  research job is cancelled. All diagnostic GPU leases are released, no labelled
  test containers remain, and both research and teaching API/UI health checks pass.
- Generic pipeline/campaign admission is not automatically covered by exact
  Database-assessment idempotency. Mixed-version writers remain unsafe for CAS.
- AUD-008's cohort/search/measurement filters are still not all applied before
  model pagination. Do not interpret visible-page counts as whole-cohort totals.
- Missing optional fleet images and unexecuted scientific methods remain explicit
  coverage limitations. A compatible Genie3 image is not proof of other engines.
- Research activation still needs a reviewed legacy-owner/account mapping and a
  coordinated scientific-file, database, queue and replica restore rehearsal.
