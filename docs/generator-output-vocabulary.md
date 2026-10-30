# Shared generator outputs

Studio uses `output_vocabulary.py` as a versioned dictionary for generator
outputs. It is applied by the shared candidate reader and writer, rather than
being implemented separately in each tool page. Both tile jobs and pipeline
steps therefore emit the same `normalized_output` contract.

The installed co-fold prediction score and PAE sidecar catalog is described in
[prediction-measurements.md](prediction-measurements.md) and is returned by the
output-vocabulary API.

## Generator coverage

All 15 design workflows currently exposed as tiles are covered. The coverage
test checks the actual tool catalog so adding a design tile requires reviewing
this list.

| Generator family | Output represented | Treatment |
| --- | --- | --- |
| RFdiffusion, including beta pairing | Generated backbone | Preserve unknown sequence; common binder geometry and secondary structure |
| RFdiffusion3 | Backbone; sequence or refolding stages in native pipelines | Keep stages distinct; retain native whole-structure statistics |
| BindCraft2 and BindCraft2 Generator | Gradient design trajectory or validated design | Separate generator confidence from validation confidence; declare confidence scales and regions |
| BoltzGen and BoltzGen Generator | Generated structure with sequence; native validation | Same entity contract; native pTM/iPTM enter scoped confidence families |
| PXDesign | Backbone; native pipeline designs | Shared structure measurements; native evaluator scores remain identifiable |
| Genie3 | Backbone; optional native pipeline evaluation | Same backbone contract; exact native pTM/iPTM aliases for evaluated designs |
| Protpardelle-1C | Backbone; optional native sequence/evaluation stages | Same shared structural fields, without inventing confidence for backbones |
| Proteina-Complexa and its Generator | Generated structure with sequence; native pipeline results | Same entity, chain-role and geometry contract; native scores retained with their original meaning |
| RFantibody | Evaluated antibody complex | Preserve heavy/light chains and target context; common structure measurements |
| RFpeptides | Generated peptide backbone | Preserve peptide chain and unknown designed sequence; common structure measurements |
| ESMFold2 binder design | Designed/evaluated complex | Same structure and sequence contract; distogram interface proxies stay distinct from predicted iPTM |

Sequence tools also pass through this layer. They are downstream consumers of
generator results, not the main scope of this change.

## What is common

Each candidate retains its original `source_tool`, `metrics`, chain lists,
structure paths and parents. `normalized_output` adds:

- Canonical engine family, stage and output kind (`backbone_proposal`,
  `designed_structure`, `sequence_design`, `refolded_complex`, etc.).
- Designed sequence and ordered designed/target chains, plus declared backbone,
  sequence and complex-context capabilities. Unknown backbone sequences stay
  unknown. A requested length range is never treated as a measured length.
- Common measurements with a definition, scope, unit, method and direction.
- Original field/value, conversion, source tool, stage and lineage information.
- Conflicting aliases and invalid ranges, which do not produce a common value.

Examples:

| Common field | Meaning |
| --- | --- |
| `common.design.sequence_length` | Total length of the actual designed sequence, preserving its chain context |
| `common.geometry.binder.residue_count` | Residues with binder coordinates, distinct from sequence length |
| `common.geometry.binder.ca_rg` | Studio Cα radius of gyration, Å |
| `common.geometry.binder.heavy_atom_rg` | Studio unweighted heavy-atom radius of gyration, Å; requires the shared geometry-method marker |
| `common.geometry.target.alignment_rmsd` | Retained target Cα alignment RMSD, Å |
| `common.secondary_structure.binder.dssp.helix_fraction` | Studio pydssp C3 helix fraction |
| `common.secondary_structure.binder.ca_trace.helix_like_fraction` | Cα-trace heuristic; separate from DSSP |
| `common.confidence.generator.complex.iptm` | Native interface TM prediction during generation |
| `common.confidence.validation.complex.iptm` | Native interface TM prediction in the validation stage |

The same predicted quantity may have different model estimators. Those methods
and source runs remain visible; a common column does not assert identical
calibration between models.

## Audit findings that affect interpretation

- RFdiffusion3 native `num_residues` can include the retained target, whereas
  `generated_residue_count` counts generated residues. They are not aliases.
- RFdiffusion3 native `radius_of_gyration` and `helix_fraction` do not explicitly
  identify the same region/method as Studio's binder measurements. They remain
  native fields.
- Older hotspot calculations used `binder_radius_of_gyration` for a Cα
  measurement. The common heavy-atom field requires the current geometry marker;
  it never guesses the method from the value.
- Native BindCraft2 source was inspected in the installed
  `mn-bindcraft2:v1.0.4-cu13-pae-v1` image. Its `pLDDT` is binder confidence in
  the complex, `Target_pLDDT` is target confidence, and
  `Unbound_Binder_pLDDT` refers to the binder alone. The declared 0–1 readings
  become separate 0–100 fields. Native `i_pAE` is divided by 31; the common
  interface mean PAE multiplies it by 31 to report Å. RFD4 native
  `oracle_ppi.ipae_coi` is a chain-of-interest/partner mean PAE in Å averaged
  across both directions; COI is the chain of interest, ordinarily the designed
  chain. Its configured oracle engine is retained on the candidate (the current
  Studio default is OpenFold3 JAX). `common.confidence.generator.interface.pae_norm31`
  divides the raw value by 31 to expose a shared numeric scale, not a calibrated
  score equivalence.
  Aggregation definitions, raw scores and producing predictors remain distinct.
- A generic `plddt`, a whole-complex pLDDT, a binder pLDDT and a composite
  confidence score are not interchangeable. Minimum, mean and directional PAE
  also remain distinct. RMSDs with different atoms, alignment or references do
  not become one `rmsd` field.
- ESMFold2's distogram iPTM proxies are not renamed to predicted interface TM.
- Several existing MPNN outputs only stored `overall_confidence` and `seq_rec`
  inside FASTA headers. These are recovered in the common projection; new runs
  also collect them as native metrics. Exact header parsing prevents `score`
  from accidentally matching `global_score`. Sequence score families retain
  their model names; Caliby energy and ESM-IF1 likelihood remain separate.

For the sequence-score definitions, see the upstream
[ProteinMPNN README](https://github.com/dauparas/ProteinMPNN) and
[LigandMPNN README](https://github.com/dauparas/LigandMPNN). Generator alias
rules use the adapter implementations and installed native source, rather than
guessing units from numeric magnitude.

## Database and downstream use

Database keeps native observations unchanged and projects common fields when
reading them. Common fields appear first in the column chooser, have readable
unit labels, and support sorting, filtering, plotting and export. The sheet now
selects a coherent evaluation, preferring the most recent representative model
or pose. Target/evaluator/run selectors restrict that choice. A missing metric
is never borrowed from another evaluation. `score_sources` identifies the exact
observation and run. Details retain every run and its original scores.

The observation context stores the producing adapter, stage, chain roles,
sequence and compact model/checkpoint settings. This also prevents an upstream
generator's name in inherited metadata from replacing the sequence tool that
actually produced the result. Sequence-only designs are not registered under
the molecular identity of a target-only PDB.

No raw metric rename is required. An additive, rebuildable SQL metric index
accelerates selection and filtering; raw observations remain authoritative. Existing local runs
are reindexed on the next ordinary Database refresh to add missing context;
observations already carrying enough context can be projected immediately.
Synchronized observations use the same dictionary on the receiving machine.
All machines need this code version for the new common view.

The shared comparison view operates on the entire filtered cohort, rather than
only the current table page. It separates native ranks from declared comparison
ranks and records missing coverage. A saved comparison freezes evaluation IDs,
metric directions, weights, thresholds and rank method in a curated group.
These definitions use the existing peer synchronization protocol; the derived
SQL index is rebuilt on each computer. See [Evaluation and QC guide](evaluations-and-qc.md).

Pipeline filters accept the same common IDs, for example
`common.geometry.binder.ca_rg`, across candidates from several generators.
The existing candidate paths, chain roles and parents remain the handoff to
MPNN, folding and evaluation steps, including existing pipeline fan-in.

Authenticated API discovery:

```text
GET /api/v1/molecules/output-vocabulary
GET /api/v1/molecules/metrics?project_id=...
```

The dictionary endpoint reports all canonical IDs, alias/conversion rules,
generator workflows, stage restrictions, units, scopes and methods. The metrics
endpoint only reports fields belonging to molecules the caller can access.
