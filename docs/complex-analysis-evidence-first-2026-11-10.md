# Complex Analysis: evidence-first saved work

The Database-only workspace now opens source collections as saved evidence,
without requiring a structure or fetching coordinates. No new scientific
identity, schema, scheduler, replication mechanism or Results scan was added.

- Saved work → Evidence → Inspect & Visualize → Shortlist & Export.
- Discovery includes structureless saved observations and curated molecular
  libraries. Input-only memberships remain distinct; no promotion is performed.
- Existing Evaluation target references supplement missing collection labels.
  Multiple target/preparation references remain separate; names do not establish
  correspondence or merge campaigns.
- Evidence has summary cards, method-filtered pages (25 observations), exact
  metrics, original parameters/provenance and the existing metric vocabulary.
  Existing Evaluation descriptors are reused; absent projections use the normal
  canonical descriptor. These are read-only views, not newly stored Evaluations.
- GB and PB, native confidence, affinity predictions and docking scores retain
  their namespaces, scope, missing values and measured zero values.
- Model metadata is paged. Explicit model loading reuses authorized Database
  retrieval, CandidateSet acceptance evidence and the existing Mol* viewer,
  target alignment/contact maps, Inspector and shared shortlist/export controls.
- Protein assessment choices follow molecular inputs. Protein–ligand assessment
  reuses the existing PLIP/PandaMap/PoseBusters Database adapter for docking or
  cofolding complexes; the source engine does not select a separate UI.
- No calculations start on opening evidence, choosing a method or loading a
  structure. New assessments require explicit submission.

## Verification

73 focused regressions passed across saved evidence, collection state/acceptance,
collection labels, authorized catalog and Database hydration. New checks cover
structureless MD evidence, paging, ACL rejection/revalidation, measured zero,
GB/PB separation, original Evaluation/model identity, target references,
metadata-only opening and selection → evidence navigation.

On the isolated production Database snapshot, summary lookup for 673 accessible
groups took 0.098 s; the first 25 observations of the THRalpha MD reference group
(79 total) took 0.069 s. These are query timings, not whole-page browser timings.

Browser: the Piccirillo collection opened 16 saved Vina/GNINA/Boltz-2 observations
without coordinates. Explicit loading rendered its original Vina ligand pose
file/ensemble in Mol*. No receptor was invented or substituted for this
ligand-only representation. Screenshot: evidence-first-viewer.jpg in the
isolated preview directory. The final stable browser session opened the THRalpha MD collection without a
structure and displayed all 79 observations / four methods. Selecting MM/GBSA
showed its 23 records separately from PB, with exact target/preparation IDs.
Screenshots: evidence-first-md-dashboard.jpg and evidence-first-md-gb.jpg in the
same isolated directory. No assessment submission or scientific calculation was
performed. The viewer preserved the source representation as a file/ensemble;
it did not invent individual pose/model identities.

Preview data/configuration: /tmp/mds-scientific-work-review-20261110, port 8510.
Production data/services are unchanged. The preview file watcher was disabled
for final testing after module hot-reload errors; no application service change
or deployment was performed.

## Boundaries

The evidence view can review MD summaries, GB/PB endpoints and their saved source
references without a 3D representation. A full trajectory player or new MD
submission adapter is not introduced. Sequence/chemical identities without
coordinates are inspectable; the existing cohort batch adapter still requires
selected structural inputs. Remote coordinates use existing retrieval and its
availability/authorization requirements. No missing measurements, chain roles,
chemical identities or acceptance decisions are inferred.

## Relevant implementation

- complex_analysis_catalog.py: bounded, authorized observation and target-reference queries.
- complex_analysis_evidence_ui.py: saved-work navigation, evidence dashboard and explicit model selection.
- complex_analysis_collection_ui.py: general discovery, input-appropriate assessment actions and inline Mol* inspection.
- complex_analysis_batch.py: existing evaluator choices by molecular input.
- complex_analysis_ui.py: workspace description; existing Results deep links are unchanged.
- test_complex_analysis_saved_evidence.py and the existing collection-label tests.

## Summary wording correction

The initial screenshot showed only the selected MD reference collection, not a
Database-wide inventory. Its 79 observations were a technical record count,
not 79 designs or simulations. The summary now identifies the selected project
and work type, displays MD evidence by method (18 trajectory-analysis, 23 GB,
23 PB and 15 production records), and keeps the aggregate record count under
Sources & provenance. Other work uses molecular, structural and evidence-type
counts. Historical/superseded entries are not relabelled as independent runs.
Focused follow-up: 34 saved-evidence and collection-label tests passed.
