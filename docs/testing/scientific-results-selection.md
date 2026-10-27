# Scientific results and candidate selection — Phase 2

Implemented on `fix/scientific-pipeline-validation`, starting at `d0d77ca`.
No catalogue reorganization, new prediction engines, model downloads or GitHub push.

## Outcomes and scheduler compatibility

`result.json` distinguishes `execution_status` from `scientific_status` (`valid`,
`partial`, `failed`). Execution completion describes the batch adapter, not the
validity of its requested measurements. Original container logs remain authoritative
for computational errors. Candidates carry `raw_metadata.measurement_status`, with
the validity of each required field and its original error. QC retains its native
`qc_status` as well. Physical output domains are validity checks, not binding-quality
thresholds.

A batch with some usable candidates keeps scheduler `success=true` and reports
`scientific_status=partial`. A batch with zero usable required measurements has
`success=false`. Missing structures are retained as per-item failures. Missing PAE
now produces a failed analysis collection with candidate causes instead of dropping
the candidates or aborting before writing results. The graph can expose such failed
analysis outputs to Results while retaining the failed node/job status.

ipSAE requires its actual `ipsae` measurement; an auxiliary finite PAE cannot mask
a missing ipSAE result. Rosetta interface analysis requires `rosetta_interface_dG`.
QC requires the tool's declared measurement fields. A finite old score plus a new
error is still failed evidence. Native tool acceptance remains separate.

## Prediction/QC evidence

Current-stage native keys and `common.` keys retain their current-stage meaning.
Use an explicit selector to read a preceding stage:

```text
evidence/complex_refolding/optimized_prediction_boltz2/iptm
evidence/complex_refolding/optimized_prediction_boltz2/common.confidence.validation.complex.iptm
evidence/analysis/ipsae/ipsae
evidence/analysis/sequence_qc/qc_sequence_gravy
```

The tool is the exact normalized `source_tool`. Repeated stage/tool observations
are ambiguous unless the candidate ID is supplied before the metric:
`evidence/<stage>/<tool>/<candidate_id>/<metric>`.

The selector traverses saved `source_candidate` lineage without copying ancestral
metrics into the current metrics dictionary. Analysis children contain their own
measurements; prediction and graph handoffs retain source snapshots and parent IDs.
Optimized prediction tools now participate in the existing common metric vocabulary.
`final_evidence.csv` and the Results evidence table identify source candidate, stage,
tool, status, native/common field, units, direction and an exact selector. Effective
settings, structures, PAE paths and full ancestry remain in the candidate snapshots
and production job artifacts.

Different-sequence ancestors, including parental antibodies, remain visible as
comparators with their exact sequence and structure. Their measurements cannot be
borrowed for the redesigned sequence, including through an intermediate sequence
change. This preserves an existing parental comparator; it does not automatically
run a matched parental prediction/QC experiment or invent parental measurements.

## Optional selection

Use the existing Filter candidates node. No starter receives a mandatory filter.
The existing single threshold, native-pass flag, direction and limit still work;
additional JSON lists support required fields, multiple rules and manual IDs.

The reusable Python contract is:

```python
from mn_molecular_design_studio.candidate_selection import evaluate_selection, selected_candidates

outcomes = evaluate_selection(candidates, {
    "required_metrics": required_measurement_selectors,
    "rules": task_specific_rules,  # [{"field": selector, "operator": ">=", "value": configured_value}]
    "ranking_metric": chosen_measurement_selector,
    "direction": "higher",  # or "lower", according to that measurement
    "limit": experimental_capacity_per_protocol,
    # Optional, using existing sequence/interface/Foldseek clustering output:
    # "diversity_field": "evidence/analysis/sequence_clustering/diversity_cluster",
    # "max_per_cluster": 1,
})
shortlist = selected_candidates(outcomes)
```

Every required/ranking/cluster field must be finite, valid and from an unambiguous
usable source. Rules are an AND conjunction. Status is `pass`, `reject`, `missing`,
`failed` or `not_evaluated`; `selected` is a separate decision. Passing candidates
omitted by shortlist/diversity quotas remain `pass` with an omission reason.
Manual IDs permit explicit manual selection; configured requirements still apply.
A node with no policy forwards candidates without asserting acceptance.

Ranking uses one measurement within a recorded protocol, with candidate ID as the
deterministic tie breaker. Common measurements enforce their declared direction.
Native fields require the configured direction; no unit or direction is inferred
from numerical magnitude. Cohorts include tool/stage, effective protocol and
prediction settings, target identity and image information. Different PAE filenames
and worker counts do not split an otherwise equivalent analysis protocol. Unknown
protocols are conservatively separated by candidate identity. The limit is **per
comparable protocol**, so it does not promise a global limit across incompatible or
unknown protocols. No universal ipTM/ipSAE acceptance thresholds or weighted score
are supplied.

Each filter retains `artifacts/pipeline_nodes/<node>/artifacts/selection/`:

- `all_candidates.jsonl`: every decision, including failed and rejected candidates.
- `selection_decisions.csv` and `final_evidence.csv`: inspectable decisions/evidence.
- `selected_candidates.jsonl`: exact sequences, full provenance and exported structures.
- `selected_sequences.fasta`: one record per designed chain; JSON preserves the exact
  colon-separated chain sequence, chain IDs and original candidate IDs.
- `structures/`: byte-identical selected structures, with SHA-256 in each export record.
- `selection_manifest.json`: policy, counts and intended downstream evaluation purpose.

Only selected candidates proceed to computational child nodes. A direct Results
consumer also reads the filter's complete review collection, so failed/rejected rows
are visible beside selected rows. Earlier filter archives remain available when
further prediction/analysis nodes intervene. Selected exports require every declared
structure file to exist; sequence-only candidates remain explicitly sequence-only.

The legacy ranking workflow now preserves all normalized candidates and records a
shortlist instead of deleting rows outside its top-N window. Its implicit universal
threshold defaults were removed. `analysis_score` remains a null compatibility
column; no composite is calculated. Callers depending on old defaults or inherited
top-level analysis metrics must configure their criteria and use qualified evidence.
The candidate schema version, parent IDs and scheduler terminal states are unchanged.

## Reproductions and verification

Before repairs, `test_measurement_semantics.py` and `test_candidate_evidence.py`
produced **9 failures and 1 pass**. They reproduced completed QC with no finite
required measurements, auxiliary-score masking, dropped missing-structure inputs,
and the inaccessible prediction score after QC. Explicit outcome/partial-batch
regressions also failed before implementation. The previous phase's all-failed
ipSAE test already passed at the starting revision; it was retained, not claimed as
a newly reproduced defect.

CPU coverage includes successful/partial/all-failed measurement batches, missing and
nonfinite values, physical domains, stale errors, stage-qualified access, ambiguous
stages, parental/child lineage, predictor vocabulary, no-policy/manual selection,
threshold direction, protocol separation, deterministic ties, diversity quotas,
all-candidate retention and byte-exact sequence/structure export. A real in-process
graph source → filter → Results test verifies both selected and rejected review rows.

The real-container tests add a fresh experimental-reference candidate → Boltz-2 →
ipSAE → sequence QC → shortlist/export chain and a deterministic replay of a genuine
previous GPU prediction. Each includes a missing-structure/invalid-sequence control.
The replay uses real analysis containers but performs **no new GPU inference**.
Fresh RFD4 generation → Boltz-2, native ESMFold2 VHH generation → Boltz-2, and AbMPNN
framework/cysteine protection use the existing integration infrastructure.

Commands (existing Conda environment):

```bash
python -m pytest -q
MDS_GPU_TESTS=1 MDS_TEST_GPU=0 python -m pytest tests/test_scientific_gpu.py -k real_prediction_to_qc -v -s
MDS_GPU_TESTS=1 MDS_TEST_GPU=0 MDS_REPLAY_PREDICTION=/absolute/path/to/retained/prediction/run \
  python -m pytest tests/test_scientific_gpu.py -k retained_prediction -v -s
MDS_GPU_TESTS=1 MDS_TEST_GPU=0 python -m pytest tests/test_scientific_gpu.py tests/test_vhh_gpu.py \
  -k 'real_rfd4_generation_boltz2 or real_esmfold2_vhh_generation_boltz2 or real_abmpnn_preserves' -v -s
```

Execution results, runtime/image information, artifact paths and source digests are
recorded in `scientific-results-selection-2026-10-27.json`. Logs are retained under
`.tmp/phase2/logs/`, real engine artifacts under `.tmp/gpu-tests/`.

Final execution results on 27 October 2026:

| Verification | Result |
|---|---|
| Full CPU suite | **1,396 passed, 22 skipped**, 115.92 s; optional MD/QC/preparation/GPU prerequisites explain skips |
| Final focused status/evidence/selection/graph/comparison regressions | **92 passed**, 5.25 s, after the final validity/protocol changes |
| Fresh Boltz-2 → real ipSAE → real sequence QC → selection/export | **PASS**, including failed controls; experimental 7EOW reference encoded as a candidate, not an RFD4-generated design |
| Genuine retained prediction → real ipSAE/QC → selection | **PASS**, with no new inference; rerun again on the fresh Phase 2 prediction after final selection changes |
| Fresh RFD4 generation → Boltz-2 | **PASS**, exact target/binder identities and parent handoff retained |
| Fresh AbMPNN CDR3 redesign | **PASS**, fixed frameworks/cysteines/unselected CDRs retained |
| Fresh native ESMFold2 VHH generation → Boltz-2 | **PASS**, exact VHH sequence and parent handoff retained |

Five distinct real integration tests passed, with one additional final replay
execution: six test executions and twelve retained production jobs. Four test paths
performed new GPU inference; the replay path only reran CPU analyses. All used the
existing RTX 4090 device 0 and installed caches/images. No inference was mocked.
Fresh inference retained seed 20261027 and the existing minimal settings: Boltz
one sample/one recycle/20 diffusion steps without target MSA/template; RFD4 one
30-residue design, 60 steps; native VHH one design/one optimization step. These
settings test integration, not useful binder quality. The original multichain/CCD
GPU test was not rerun in this phase; its CPU identity contracts remain covered.

## Commits and changed files

- `af2956981f85b037e4b4b57c52db624a75023941`: measurement/status correction.
- `a492b7892b736250ebea779ead59245636568042`: prediction/QC provenance integration.
- The commit adding this guide: candidate selection, Results/UI integration, tests
  and execution evidence. Obtain its hash with `git log -1 --oneline`.

Production files in `mn_molecular_design_studio/`: `measurement_status.py`,
`candidate_selection.py`, `output_vocabulary.py`, `design_validation.py`,
`pipeline_graph.py`, `pipeline_node_ui.py`, `result_presentation.py`, `run_app.py`,
and engine workflows `analysis.py`, `quality.py`, `optimized_prediction.py`,
`graph_pipeline.py`. Tests: `test_measurement_semantics.py`,
`test_candidate_evidence.py`, `test_candidate_selection.py`, `test_quality_tools.py`,
`test_evaluations_and_comparisons.py`, `test_scientific_gpu.py`. Documentation: this
guide, the execution JSON and the link in `scientific-pipelines.md`.

## Scientific boundaries and next phase

Selection is a transparent experimental prioritization decision, not evidence of
binding, affinity, activity, solubility or clinical utility. Confidence measures
remain dependent on their prediction methods and settings. QC composition metrics
are descriptive heuristics. Existing Rosetta/ipSAE chain-support boundaries remain;
no engine capabilities were broadened. No prospective wet-lab calibration, matched
parental benchmarking or broad native-generator campaign was performed. Old records
with incomplete ancestry/protocols cannot acquire missing evidence retroactively.
No automatic join across independent sibling analyses is guessed; connect a serial
evidence chain or select explicitly sourced records. Archive portability still
depends on retaining the referenced upstream job artifacts; the selected structures
and sequence JSON/FASTA are copied into the export bundle.

Recommended next implementation phase: repair F4/F6 (accurate AF2 conditioning labels
and effective Boltz controls), then establish small matched task-specific validation
sets with parental/negative controls to calibrate configurable acceptance policies.
Use that evidence before the later task-oriented catalogue redesign. Full hashing of
external dependencies and broader cache/protocol identity remain the separate F9 work.
