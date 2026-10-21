# Evaluation comparisons and batch QC

## Evaluation selection

Molecular identity is reusable. Its evaluations are separate calculations with
their own target, target revision, representation, chain roles, method, run,
settings, status and provenance. A repaired receptor is a different structure
revision even if its sequence is unchanged. Historical unknown context remains
explicitly unknown; it is not reconstructed by guessing.

The Database's default sheet selects one evaluation per molecule, preferring
the most recent representative model/pose, then the most recent calculation.
Use **Evaluation selection** to restrict the target, evaluator or run, and the
detail panel's **Evaluation** selector to inspect a specific model. The displayed
structure and scores refer to that selection. Raw history is retained. Only the
best prediction/pose remains the default portable Database representation;
other models and bulky artifacts remain accessible through Results.

The local `molecule_evaluations` and `molecule_metric_values` tables are derived
indexes. They rebuild from the existing observations, invalidate on updates and
deletion, and do not introduce a second replicated scientific data source.
Canonical evaluation IDs survive peer remapping of physical SQL row IDs.

## Compare a collection

In Database, use **View → Compare evaluations**. Results exposes the same view
under **Compare evaluations**. The entire filtered cohort is included, independent
of sheet pagination. The view offers:

- One selected evaluation per molecule, or explicitly every matching evaluation.
- Candidate name/compound ID, producing method, target, run, native rank,
  comparison rank, metric coverage and missing measurements.
- Explicit favorable directions, positive weights and threshold rules.
- Mean, worst or product of metric ranks; Pareto fronts for up to 5,000 rows.
- Scatter box/lasso selection, distributions, individual/combined passing counts,
  failure reasons, pairwise Spearman correlations and CSV export.
- A reusable selection and **Save comparison**, which freezes exact evaluation
  IDs and rules in a curated group that follows the ordinary peer sync protocol.

To compare prediction confidence together with subsequent QC, enable **Show every
matching evaluation**, then **Contributing evaluations → Combine evaluations of
the same structure**. This joins only the same molecule, representation digest,
target revision and chain roles. The evaluator/run selector chooses which
calculations contribute. Conflicting repeated values or protocols are left
missing, with the conflicting metric names shown; a later calculation does not
silently replace an earlier one. CSV exports and saved comparisons retain all
contributing evaluation IDs. Reopening a saved comparison restores these choices.

Missing or failed measurements never receive an implicit zero. A row must have
every selected ranking metric and pass the thresholds to receive a comparison
rank. Different targets, methods/scopes, or scientific parameter signatures
block pooled binding-score ranking. Execution placement, CPU counts and source
paths do not change the scientific signature. Native ranks are never overwritten.

The legacy analysis adapter also stopped suffix alias matching, whole-complex
confidence substitution for binder pLDDT, numeric-magnitude scale guesses, and
sums with inconsistent metric coverage. Its derived score uses only measurements
available with a consistent method throughout that cohort.

## Analysis tiles and pipeline nodes

Select completed project results on a tile's input page. Alternatively select
protein rows in Database and choose **Use in tool**, or connect a Database
selection / generator / prediction collection to an analysis node in Pipelines.
Every collection stays one visible job. CPU workers are capped by the scheduler
grant; per-candidate errors remain in the batch and do not disappear.

| Tile | Protocol / required input |
| --- | --- |
| Sequence Composition | Standard amino acids; length, fractions, pI, declared-pH charge, entropy, GRAVY, aromaticity and instability index. Unknown residues cause an explicit item failure. |
| Protein-Sol | Original sequence predictor, at least 21 amino acids; raw, scaled and population solubility. |
| Biotite Interface | Declared binder and target chain groups; heavy-atom and residue contacts, residue maps, 4 Å charged-side-chain salt bridges. Hydrogen bonds require explicit H atoms; their absence is reported as unavailable. |
| Sequence Clustering | Global alignment identity divided by longer sequence length; single-linkage threshold, at most 2,000 items. |
| Interface Clustering | Jaccard similarity of contacted target residues; one identical target revision/numbering, at most 2,000 items. |
| Foldseek | Version 10-941cd33, TM alignment and coverage threshold; one selected protein chain per item. |
| PEP-Patch | Selected protein chains; Jain-scale hydrophobic patches and electrostatic patches at declared pH, areas and surface fractions, per-residue artifacts. |
| ESM-1v Scoring | Five-model ensemble mean native amino-acid probability with unmasked wild-type marginals. This is **not** masked pseudo-likelihood. Limit 1,022 residues; no implicit windowing. |
| ESM-IF1 Scoring | Native-sequence conditional log likelihood with complex-backbone context; declared sequence must match the selected structure chains. No sequence sampling is performed. |
| PLIP | Protein–ligand interaction profiles on best saved docking poses or declared cofolded ligands. |
| PandaMap | Protein–ligand interaction maps and its explicitly labelled empirical estimates; source protein residue numbering is retained. |
| PoseBusters | Pose plausibility checks using a coordinate-bearing ligand whose connectivity is matched to the source SMILES, plus its receptor. A failed plausibility check is a completed evaluation with a negative result. |

The ligand analysis tiles accept an existing docking job even when it predates
the normalized-candidate contract. They select the best available pose per
compound using the source ranks/scores. Cofolding complexes need declared ligand
SMILES and an unambiguous organic ligand residue. Native atom order is not used
as a substitute for a chemical graph match. Ambiguous ligands and failed topology
assignment produce explicit per-item failures. Combined analyzed complexes remain
linked to the compound's Database identity.

Clustering preserves all rows. Cluster representatives follow source order;
cluster numbers are identifiers, not quality scores. The ordinary downstream
filter node can select `common.diversity.representative >= 1` when one item per cluster is
desired. Multichain Foldseek clustering is not silently approximated by treating
independent chain clusters as a complex cluster.

The Results page displays batch progress, item failures, a batch metric table
and source-numbered interface contact maps. Native contact/patch maps and
residue probabilities remain downloadable artifacts. Numerical scalar results
are ingested with the same evaluation/vocabulary contract as existing tools.

## Containers and reproducibility

Recipes live in the sibling `mn-tool-containers` repository:

```bash
./build.sh protein-qc
./build.sh proteinsol
./build.sh foldseek-qc
./build.sh pep-patch
```

The resulting images are `mn-protein-qc:1.0.0`, `mn-proteinsol:1.0.0`,
`mn-foldseek:10-941cd33-qc1`, and `mn-pep-patch:65da35b-qc1`. All four have been
built on PC002975. Protein-Sol's source archive is checksum-pinned. Foldseek,
PEP-Patch and scientific Python dependencies have explicit versions/revisions.
The ESM scoring tiles reuse `mn-esm-if1-opt:2.0.1-2b369911` and require the
declared `reference_files/esm_models` checkpoint files before execution.
PLIP, PandaMap and PoseBusters reuse their installed MN Ligand analysis images
and batch adapters; no dependency on the external MN Ligand application is added.

Each job records the actual image digest, checkpoint SHA-256 values where used,
settings and a hash of the adapter/configuration code. It copies its standalone
runner and helper code into `artifacts/runtime` before launching the container.
Later source edits therefore do not change those saved runtime files. Containers
run without network access. Rootless Docker is supported.

OVO helper adaptations retain their MIT license and upstream revision under
`workflows/qc_ovo`. Studio does not depend on an OVO checkout, OVO's application
database, Nextflow, or the separate MN Protein Design / MN Ligand applications.

## API

All routes use existing Studio authentication and group access controls:

```text
GET  /api/v1/molecules/output-vocabulary
GET  /api/v1/molecules/evaluations/options
POST /api/v1/molecules/evaluations/search
POST /api/v1/molecules/comparisons/preview
POST /api/v1/molecules/comparisons
GET  /api/v1/workflows/analysis.sequence_qc/schema
POST /api/v1/workflows/validate
POST /api/v1/workflows
```

Example preview (IDs must belong to the authenticated user's accessible groups):

```json
{
  "selection": {"group_id": "GROUP", "evaluation_filters": {"run_id": "RUN"}},
  "settings": {
    "fields": ["common.interface.ipsae.minimum"],
    "directions": {"common.interface.ipsae.minimum": "higher"},
    "thresholds": [{"field": "common.interface.ipsae.minimum", "operator": ">=", "value": 0.5}],
    "method": "mean"
  }
}
```

Submitting an analysis uses the usual workflow payload with `source_job` or
explicit managed `source_run_dir`/`candidates_jsonl` paths. A multisequence
collection is already a batch; callers should not submit one workflow per row.

## Checks and remaining optional extensions

Run pytest with the Studio Conda Python. The independent Docker checks can be
enabled with `MDS_QC_DOCKER_TESTS=1 pytest tests/test_quality_tools.py`. They
exercise two-item batches in all ten CPU analysis protocols. Each item has a
private working directory, including native tools that use fixed intermediate
filenames. ESM adapters have
schema/graph/placement checks; their full model inference needs a separate smoke
run. Peer tests cover frozen comparison selection across synchronization.

Promb remains optional and unregistered pending an explicit versioned reference
proteome. No scores are labelled as immunogenicity, experimental affinity, or
validated solubility.
