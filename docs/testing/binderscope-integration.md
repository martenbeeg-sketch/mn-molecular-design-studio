# Integrated BinderScope verification — 2026-11-17

## Source and packaging

The actual analysis library and offline review resources from BinderScope
revision `c61ef9ca3dfc64dbde018ae7329d7538d8c583f1` are vendored with their MIT
licence. `SOURCE_MANIFEST.json` records SHA256 values for all 50 imported files.
Only the final newline of the bundled 3Dmol JavaScript was normalized; Python
analysis source is unchanged. The runtime verifies that inventory before work.

The built wheel contains every manifest-listed file, including validation
schemas/scripts, the licence and offline viewer assets. The CPU image contains
both the engine and the Studio adapter. Its build does not clone BinderScope.
Its runtime does not mount a separate source checkout or fetch code. Adapter
hashes, engine source inventory hash and the actual image ID are retained.

## Executed checks

184 focused tests passed with `MDS_TEST_BINDERSCOPE_IMAGE=1` across:

```
tests/test_binderscope.py
tests/test_data_browser.py
tests/test_data_browser_analyses.py
tests/test_data_browser_table.py
tests/test_database_browser_cache.py
tests/test_binder_interface_chemistry.py
tests/test_complex_analysis_comparison.py
tests/test_workflow_contracts.py
tests/test_post_phase2_release.py
```

The new tests exercise exact-chain role mapping without chain-order fallback,
canonical-to-author mapping and insertion codes, unchanged source coordinates,
legacy explicit confidence/topology requirements, per-design recorded settings,
unknown-confidence geometry-only execution, single-model versus ambiguous-model
selection, replica authority, ACL denial, checksum rejection, API preflight,
explicit-only dialog submission without navigating away, read-only preview
submission blocking, table-subset versus all-pages selection, project-boundary
partitioning, source inventory and the standalone offline image.
Confirmation-table Include checkboxes preserve exact design identities and edits
while excluding unchecked designs. Clear-selection and select-all callbacks are
tested, including disabled submission for an empty selection. A live browser
check with isolated synthetic records confirmed selection in the main browser
table and deselection inside the confirmation dialog; these records were never
added to the research database.
The real container-to-SQL test attaches an Evaluation to the original molecular
and representation IDs, creates no extra molecule, retains the complete
per-design result, and denies its inspection to another user.
Both legacy settings and the new per-design geometry-only path were exercised
through actual offline containers and SQL ingestion. A mixed-topology cohort
test confirms that a cyclic binder has no terminal assessment, unknown B-factors
do not appear as pLDDT, geometry tracks remain present, and target B-factors stay
uninterpreted even when binder confidence is declared.

An additional manual native-campaign check used eight synthetic upstream toy
designs and a supplied full-reference structure. The integrated image ran the
native pipeline, validated the export and produced all offline review pages.
It was also rerun with only the input/output data folder mounted: neither the
external checkout nor the Studio runtime source was mounted. Original native
CSV counts passed the upstream validator. The fixture generator was used only
to create disposable test inputs, not as an application dependency.

The earlier signed-in read-only 8514 preview was checked without resetting credentials.
A 301-row Tau R3 selection spans multiple generators and all seven table pages;
preflight prepared 201 exact structures and explicitly excluded 100 historical
rows instead of guessing roles. The initial source-workstation restriction was
subsequently removed: verified peer structures now use the existing local
Evaluation collection and owner-qualified source-reference synchronization path.
That earlier subset reached the standalone tile; it is now superseded by the
browser-first confirmation dialog. The updated UI is additionally checked using
isolated, synthetic admitted designs from two generators, without bypassing the
research preview's sign-in. Submission remains disabled in previews; no research-
data calculations or production writes were made.

## Limits

These are focused integration checks, not the complete application suite or a
validation of all upstream modalities. Specificity with supplied offline species
data, optional antibody backends, every scaffold format, licensed ConSurf data,
and distributed launch/cancellation were not exercised end-to-end here. Native
campaign/version support remains bounded by the upstream readers. Saved SQL
subsets cannot reconstruct rejected attempts, native funnel/stopping statistics
or missing scaffold/reference associations; the tile and documentation state
this distinction instead of claiming input-independent parity.

No production service deployment, credential change, data cleanup, commit or
push was performed for this integration.

## Replicated-design correction — 2026-11-18

The source-workstation-only preflight was incorrect for this workflow. The
existing shared-database ingestion already supports new locally authored
Evaluations on exact peer structures without modifying the owner's collection.
BinderScope now uses that path. The method contract is
`studio.binderscope.saved-cohort.v3`.

The reported real design, `Design_pc004032_20261019T163649.411014Z`, exists in
both research SQL and the preview snapshot. Its verified saved complex contains
31 target residues and 100 binder residues. The legacy SQL identity describes
only the binder, with an explicit member-specific SQL target link to Tau R3.
Both sequences now resolve unique chains against the saved coordinates. Role
evidence is revalidated during staging and ingestion; deleted, inaccessible,
conflicting or nonmatching target associations cannot authorize that mapping.
The calculated candidate describes the complete coordinate complex, while its
Evaluation still attaches to the original binder identity and saved structure.

That native mmCIF lacks the optional occupancy column. The standalone runtime
now supplies occupancy only in memory before normalization. Missing or unknown
B-factors stay unknown rather than becoming Gemmi's default value; compressed
CIF input is also covered. Source coordinates and author/insertion mappings
remain unchanged.

Executed verification:

- **212 tests passed**, including all files listed above and
  `tests/test_molecular_database_peer.py`, with the real offline image enabled.
- **146 additional ingestion, identity-integrity, policy and bound-docking
  tests passed**: `test_molecular_database.py`,
  `test_molecular_database_integration.py`,
  `test_molecular_database_identity_integrity.py`,
  `test_molecular_database_policy.py` and `test_bound_docking.py`.
  Total: **358 passing regression tests**, with no skips in these two runs.
- An offline container test starts from a replicated binder-only SQL identity
  and native CIF, stages and calculates BinderScope, ingests one attached
  Evaluation without creating a new molecule, and replays it idempotently to
  the structure owner and a third database. Source-owned rows remain unchanged.
- A read-only legacy-preflight test forbids broad `list_molecules` hydration,
  preventing accidental SQL projection writes from the development snapshot.
- In the signed-in live 8514 preview, the exact reported selection changed to
  **1 eligible / 0 excluded**; its full 65-design cohort reports
  **65 eligible / 0 excluded**, including both table pages.
- An additional offline runtime check on an isolated copy of the actual
  reported complex completed successfully: buried area 2204.56 Å² and 52
  contacting residue pairs. Original bytes were unchanged and no research SQL
  observation was inserted. Unknown topology/confidence still correctly leaves
  confidence-dependent terminal assessment unavailable.

The rebuilt local image is
`sha256:458d15808e12a13459d96ba3a965893cb1d7e3415cbe5ec2aa7cd06f986bfd9e`.
The live development preview is refreshed; it remains read-only, so queueing is
deliberately disabled there. Production services and credentials are unchanged.

## Hotspots and original input targets — 2026-11-18

The saved-cohort contract is now `studio.binderscope.saved-cohort.v4`.
`binderscope_sources.py` obtains source-specific native hotspot declarations and
verified original target references through SQL access checks. The source
snapshot is rechecked at staging and Evaluation ingestion. Reference files are
copied into the authorized input bundle and checked again before execution.
The independent runtime helper is packaged in the image, hashed in requests and
reports, and included in the protocol revision.

The normalization map now composes original-target residue identity with saved
target numbering and canonical BinderScope numbering. Tests cover changed chain
names and arbitrary new residue numbers, author insertion codes, unique cropped
segments, ambiguous repeated segments, changed sequences and unavailable
references. Uncertain mappings remain unmapped. Every design/run retains its
own intended hotspot set; the upstream interface and binding-mode calculations
receive the correctly translated residues rather than one cohort-wide union.

Verification:

- **368 integration/regression tests passed** across the previous 15 test files
  plus `tests/test_binderscope_hotspots.py`, with offline image tests enabled.
  An additional **1 result-viewer test passed**, verifying original and
  re-predicted residue labels and the sequence-only fixed-fold warning.
- Actual container-to-SQL tests with both local and replicated original targets
  attach hotspot contacts, distances and original-target RMSD to the original
  design Evaluation. Reference checksum changes and stale hotspot settings are
  rejected. Admission decisions remain unchanged.
- The live read-only 8514 preview shows the real PD-L1 design
  `Design_pc004033_20261101T222355Z` as eligible with **3 requested hotspots**
  and an **available original target reference**. No research job was queued.
- An isolated offline check of that actual design maps A40, A99 and A107
  unambiguously and reports **0 contacted, 3 missed, 0 unmapped** at 4.5 Å.
  Minimum heavy-atom distances are 7.994, 5.723 and 6.870 Å respectively.
  The 115-residue target Cα RMSD to its verified original input is 0.170 Å.
  This new contact assessment is not a reinterpretation of the original filter
  pass: historical admission may use different criteria. No research SQL result
  was inserted by this isolated verification.

The updated local image is
`sha256:2ab83b3bbc723b2ebc7b170dcf4f7ad8e4df80dd8c8920630fe6ffdf65a535d4`.
Production services, original structures, existing scores and credentials remain
unchanged. Reference specificity, noncontiguous/mutated-target alignment and
full-length scaffold/partner placement are not inferred by this new SQL path.
