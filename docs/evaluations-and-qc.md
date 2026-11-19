# Evaluation comparisons and batch QC

## Evaluation selection

Molecular identity is reusable. Its evaluations are separate calculations with
their own target, target revision, representation, chain roles, method, run,
settings, status and provenance. A repaired receptor is a different structure
revision even if its sequence is unchanged. Historical unknown context remains
explicitly unknown; it is not reconstructed by guessing.

The Database's **Collections** sheet selects one evaluation per molecule, preferring
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

## Scientific outcomes in Molecules and Runs

The three-column browser defaults to **Scientific outcomes**, not a list of
calculations. Choose a target/ligand or a run on the left; select a retained
outcome in the centre; inspect its structure and scores on the right. Left-hand
filters include molecular type and whether the reference itself has recorded
coordinates or is sequence-only. A design complex does not make a sequence-only
target a structure-backed reference.

- **Design** groups exact designed chemical contents for one recorded target.
  Only retained generated/curated sources contribute primary rows; intermediate
  sequence/backbone stages and analyses do not become additional designs.
  Chaperone is a reference for the design-with-attached-history concept, not a
  replacement database or a requirement for this view. Its optional bridge is
  unchanged.
- **Structure prediction** describes independent monomer predictions.
  **Cofolding** describes independent complex predictions, including ligand
  complexes. **Binding/docking** keeps docking cases distinct from cofolding.
- Explicit design-pipeline refolding is attached to the design. Merely using a
  design sequence as input to a separate prediction does not make that prediction
  pipeline refolding. New ingestion records workflow purpose in provenance;
  ambiguous historical purposes stay **Unclassified / Needs review**. This does
  not force reingestion of unchanged historical jobs or a database migration.
- **Admitted results** is the default. Designs appear only after a recorded
  run-filter pass or an explicit later selection. Historical SQL membership alone
  is not admission. Unreviewed legacy records stay preserved in diagnostic
  history, not mixed into the admitted design list. A later explicit selection
  can retain a non-passing design without changing its original filter decision.
- **Passed / recorded** restricts the view to explicitly passing retained designs and prediction
  or binding cases with recorded coordinates. **Filter evidence** identifies the
  recorded acceptance stage: native filtering is not labelled final acceptance.
  Completed execution and a successful analysis are never acceptance evidence.
  A required campaign gate takes precedence over native diagnostics.
- Explicit manual selections without a pass decision are **Manually selected**.
  Unknown or conflicting legacy acceptance remains **Needs review**. Rejected
  evidence remains in **Calculation history** unless the design is explicitly
  retained later; that selection is not relabelled a filter pass.

One navigational design can have original evaluations from several runs. Every
source retains its own measurements, model, target revision, version and settings;
values are neither averaged nor replaced with arbitrary latest-tool scores.
Binding/cofolding cases preserve target revisions; missing revision evidence
keeps executions separate. Unknown chemical identities/targets are not guessed.

**Calculation history** restores individual evaluation/structure/record rows,
including optional input/reference history. **Collections** remains the folder
and imported-library view. These are read-only projections: no historical record
is deleted, scientifically relabelled, or automatically reimported. Policy-less
historical runs remain Results-only during ingestion.

## Connected analyses for one design

In the Database's **Molecules** and **Runs** views, select a result in the middle
column. In Scientific outcomes, **This molecule** shows **Original result scores**
and **Downstream evaluations** together beneath the existing structure viewer.
The middle table contains admitted designs, not separate post-analysis entries
or post-analysis score columns. Linked ipSAE, Biotite, PEP-Patch and other analysis
data are accessed through the selected design's structure viewer on the right.
Its status table also counts failed attempts,
unknown execution status and results on other models. A dash means no verified
linked SQL result, not proof an analysis never ran. Detailed scores and technical
identifiers remain under **Analysis scores & provenance** and **Models & calculation
history**. Names shown in the overview never replace scientific identity keys.
These include ipSAE, PEP-Patch and Biotite Interface when recorded, and explicit
design refolding. Calculation history also offers **Design analyses (N)**.
Each analysis keeps its own native measurements, run, model/pose,
source representation, calculation digest, engine version and scientific settings.
The **Evaluation** selector can open any verified linked calculation with its
own scores and recorded structure. Changing the displayed evaluation does not
hide the original producing-tool scores in Scientific outcomes. **Scientific
context** exposes representation, target revision, chains and setting digests;
**Evaluation details** retains the selected calculation's full recorded context.

Browsing reads saved SQL evidence, not arbitrary work-directory outputs. The
metadata and derived target lists are cached per session and access scope; fresh
group permissions and SQLite commits invalidate them, including peer edits and
deletions. Version-valid evaluation descriptors avoid repeatedly parsing full
observation contexts. Missing or stale projections fall back to the original
observations without requiring browser-driven migration or rewriting history.
The first load still scans authorized metadata; subsequent selections reuse it.

Binder-only identities and complexes containing the exact same designed chains
can share this navigation. Changed sequences, ligand stereoisomers and different
binder stoichiometries remain separate. Explicitly different binding targets are
excluded. Analyses of another structure of the same design, derived structures
and sequence-only checks are labelled separately. Names and candidate aliases
alone never establish a relationship. Unresolved source/target links are counted
but are not automatically attached; their evidence remains in result history.

The three-column layout and left-hand filters are preserved. If a design has
more than 20 post-evaluations, use **Analysis scores page** on the left to page
through the score panels. Access changes, peer updates and deletions invalidate
the session-local metadata cache; opening an individual evaluation rechecks
access and loads only that evaluation's scores and representations.

This is read-only navigation, not a new scientific acceptance policy or a repair
of historical relationships. The 2FLU records without acceptance evidence need
deliberate review before they can appear in the passed shortlist; their absence
from that shortlist is not data loss. Cold index loading still scans the existing
catalogue; session caches avoid repeating the projection on unchanged reruns.

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
