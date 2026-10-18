# BindCraft2 and BoltzGen regression check — 16 October 2026

## Scope

Review of current input forms, API workflow schemas, staging, scheduler dispatch,
container execution, candidate normalization, and confidence artifact readers.
Jobs use the installed images and common reference weights, through the live API
and scheduler on GPUs 0 and 1. The paused V7 campaign is outside this check.

These are functional smoke runs. BindCraft2 uses one trajectory, one step per
design stage, one redesign, one model/recycle, and relaxed confidence filters.
BoltzGen uses one 60-residue design, two generation sampling steps, and a budget
of one. These settings do not assess binder quality or scientific performance.

## Fixes

- The shared path reader interpreted inline `design_yaml` and `target_msa_a3m`
  content as relative filenames. A queued YAML specification could gain a job
  directory prefix before its `entities` key. Inline content now stays literal
  during job storage, worker restoration, and non-admin API path validation.
  Actual file fields retain their path resolution and access checks.
- A ligand-only BoltzGen run wrote `str(None)` as its target structure in a later
  input snapshot. This now remains an empty value, consistent with initial staging.
  The prediction itself did not require this nonexistent target.
- The old BoltzGen YAML normalizer remains for compatibility with older requests.
- BoltzGen's native refolding archive contains five samples. The result reader
  previously left its three-dimensional PAE unreadable for a ranked complex.
  It now selects the same sample as the pinned engine writer and analysis code:
  `argmax(0.8 * design_to_target_iptm + 0.2 * design_ptm)`. Ranked filenames are
  matched to the archive, and native zero-based sample filenames are distinguished
  from the cofolding adapter's one-based exports. The Streamlit confidence cache
  function was updated so cached results are recomputed with this reader.
- BindCraft2 sequence targets may export a sampled crop instead of the entire
  FASTA sequence. The candidate normalizer now recognizes target substrings and
  uses the native binder sequence for joined crops, avoiding reversed target and
  binder chain labels. Tests cover contiguous and joined target windows.

## Automated checks

`python -m pytest -q`: **163 passed**, one existing Starlette/httpx deprecation
warning. Added cases cover inline YAML/A3M/text preservation, actual Path worker
arguments, non-admin path checks, queued BoltzGen YAML staging, and ligand-only
input snapshots. Further tests check native sample selection and prevent use of
an unrelated archive's selected sample.

## Live jobs

| Path | Job code | Result |
| --- | --- | --- |
| BindCraft2 structure target | 5036C | Passed; one exported ranked complex, pLDDT and 175 × 175 PAE readable |
| BindCraft2 sequence target | 3BFFA | Engine completed, exit 0; one refolded candidate rejected by native Target_pLDDT filter, no accepted final design |
| BoltzGen protein, full pipeline | 1C220 | Passed; ranked complex, metrics/report, pLDDT and selected 253 × 253 PAE readable |
| BoltzGen ligand, full pipeline | CED47 | Passed; ranked complex, finite affinity outputs, pLDDT and selected 76 × 76 PAE readable |
| BoltzGen custom YAML, generation | 29AA2 | In progress; inline and staged YAML retain the correct entities root |

Exact request payloads and run paths are recorded locally in
`/tmp/mds-design-regression-20261016/jobs.json`.
