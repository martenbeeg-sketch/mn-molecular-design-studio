# Molecular-dynamics studies in Data browser

## Navigation and storage

Open any accessible project, target, ligand or collection, then **Molecular
dynamics → study → replica**. Selecting a study opens its analyses in the middle
pane without replacing the right-hand structure viewer. Study labels use the
compound, receptor, PDB accession, originating pose method and MD engine, where
recorded. A replica's production segments remain separately selectable.

This is a read-only projection of existing SQL observations and their Results
artifacts. It does not admit candidates, rerun simulations, or insert trajectories,
checkpoints, solvated systems or numerical arrays into SQL. Intermediate MD
evaluation collections are discoverable only through already accessible library
molecules. This exception does not make intermediate design outputs library
entries.

MD grouping includes source workstation identity, owner, project, exact molecule,
study ID, selected complex and PDB context. Production, analysis and endpoint
imports belonging to that same study appear together. A shared compound or import
job alone is insufficient to merge studies. Unknown study IDs stay run-scoped.

## Analyses

- **Overview · all replicas:** a compact MN-Ligand-style plot matrix from one
  saved aggregate analysis, with replicas 1–3 (or other recorded repeats) side by
  side. The default shows interaction maps, protein RMSD, ligand RMSD and contact
  occupancies. Plot options expose the other supported saved curves and summaries,
  separate/overlay/recorded-consensus display and replica selection. Earlier and
  later analysis windows stay separately selectable; the overview does not merge
  production segments. Original tables are collapsed under the plots.
  **Expand overview** opens the selected matrix in a wide dialog while retaining
  the ordinary right-hand structure viewer on the main page. Individual panels
  can also be enlarged or saved as SVG. The structure selectors do not restrict
  the multi-replica overview.
- **Dynamics:** stored backbone/ligand RMSD, protein and ligand RMSF, radius of
  gyration and secondary-structure fractions.
- **Contacts:** stored contact occupancies, distance series, reference-site
  retention, water/salt bridges and interaction-network tables, where available.
- **Thermodynamics:** stored temperature, density, potential energy and volume.
  Long series are sampled to at most 2,000 display points without rewriting data.
- **Endpoint energy:** separately selected MM/GBSA or MM/PBSA calculations, with
  methods, units, source production and supersession metadata retained.
- **SQL summary / Provenance:** saved scalar measurements and the context of each
  calculation, even when a source artifact is offline.

Attached trajectory-analysis reports open on explicit request, avoiding eager
loading of every report. Original HTML/JavaScript dashboards are not executed;
charts are reconstructed from the saved numerical reports, without recalculation.
Endpoint and attached analysis selectors require the selected production-segment
link when it is recorded; a matching replica alone does not establish that link.

**Contacts → Ligand atom–residue interactions** places the linked ligand's
2D chemical identity in the centre and recorded contacting residues around it.
Colours identify stored hydrogen-bond, hydrophobic, water-bridge, salt-bridge or
general geometric-contact occupancy. Display options select an interaction type,
minimum occupancy and up to 16 residues; the full numerical tables remain below.
The collapsed **Saved interaction occupancies (%) · all recorded residues** table in both
Data browser and MD Results shows every valid saved residue, including interactions
below the diagram threshold or outside its residue limit. It presents each
recorded interaction type separately as a percentage to two decimal places;
missing/invalid values stay blank and measured zeros remain zero. Changing the
selected production report or replica changes this table, without merging reports
or recalculating the historical analysis. Geometric-only analyses do not acquire
unmeasured classified interaction columns.
The diagram can be enlarged or downloaded as SVG with its source segment and
recorded analysis-window metadata. It is a trajectory-averaged schematic, not a
single-frame interaction map, a spatial binding-site layout or a binding-energy
estimate. Dashed lines now end on the named ligand atom recorded in each residue's
`top_ligand_atom`, rather than a central panel or an arbitrary atom. The ligand
block embedded in that same verified production report provides the atom names;
its heavy-atom graph and defined stereochemistry must match the linked SMILES.
Continuation reports without an embedded ligand block can use a digest-verified
topology only when its manifest pairing identifies the selected production, or
the report records the exact same prepared-system identity and filename. Merely
sharing a compound, replica or basename is insufficient. A contradictory pairing
or ambiguous ligand residue fails closed. Source atom names survive chemical-bond assignment and 2D layout; generated
SMILES atom indices are never substituted for names. Missing, duplicated,
unmatched or chemically inconsistent names yield an explicit missing-endpoint
notice and no atom link. The assignment table retains the recorded protein atom
as a separate statistic, not an asserted jointly observed atom pair.

**Historical-data limitation:** these reports preserve a representative atom and
residue-level occupancies, not per-interaction ligand/protein atom-pair occupancy.
The representative can aggregate nearest-contact and other interaction counts.
It must not be relabelled as the exact hydrogen-bonding atom for every frame.
The diagram, export and UI label this distinction explicitly. Exact pair-specific
interactions would require a separate trajectory analysis; no such analysis is
silently run or reconstructed from a snapshot.
The overview's compact interaction maps use the original report's saved named
2D ligand depiction and source hotspot ordering, with separate percentage labels
on each coloured interaction line. Replica-specific values come from the exact
replica column in `contact_matrix`, never from another replica or the consensus
as a fallback. Saved consensus values are explicitly labelled and not recomputed.
Missing/ambiguous atom names omit the corresponding link; the original app's
arbitrary nearest-atom fallback is intentionally not reproduced. Unlike the
single-segment diagram's optional threshold, this overview includes every positive
recorded interaction for the selected hotspot residues (including a 9% bridge).
General contacts are shown only where no classified interaction meets the display
threshold. Missing/invalid occupancies are not converted to zero or percentages;
missing ligand SMILES yields an explicit labelled placeholder, not an inferred
chemical structure. Native MD Analysis results also show their recorded geometric
contacts without labelling them hydrogen bonds; where they did not save atom-name
associations, no atom endpoints are invented.

Each curve uses its own recorded time axis. Stored `per_phase_local.production`
RMSD can be selected independently of the original recorded timeline. Continuation
segments are not silently concatenated, and phase-local reference coordinates are
not inferred from a topology file. Missing values remain missing; a missing ligand
series does not suppress a recorded backbone curve. Failed production is not
presented as successful analytics.

## Executable MD Analysis tile and pipelines

The existing **MD Analysis** tool also has a typed visual-pipeline adapter:
**Saved MD run → MD Analysis → Results**. Create this route from the Pipelines
page or the standalone MD Analysis tool. The source picker selects one completed
project job, one exact production segment/replica and a recorded matching
topology. Static predictions cannot connect to the trajectory input. Branch the
source into multiple analysis nodes to compare explicit sampling windows or
settings without implicitly concatenating segments.

The pipeline uses the ordinary scheduler, cancellation and child-job lineage.
The CPU analysis mounts the source job read-only, checks source hashes before
and after computation, and produces numerical artifacts and native plots, not
just a report attachment. Analysis outputs are collected separately from design
candidates; running this route does not admit new designs to the molecular
library. Imported studies retain their recorded source-study link for SQL
attachment. A native run without a verified study link remains available in
Results rather than being assigned to an inferred SQL molecule.

MD Analysis streams the selected trajectory window and reports RMSD, RMSF,
radius of gyration, and frame-by-frame protein–ligand interactions. It writes
per-frame residue counts, per-residue occupancies, and an atom-pair CSV with
protein and ligand atom names, frame counts, and occupancies. A residue is
counted only once per class per frame, even if several of its atoms interact.
Units, sampling, timing, imaging policy, and input hashes are retained.

In **Pipelines → Capabilities → Analysis & validation → Protein–ligand pose
analysis**, **Native Ligand Interactions** applies the same interaction classes
to one docked or cofolded pose. It is also available from **Analyze this pose**
or **Analyze selected structure** in My results. Use an exact ligand selection
when a structure has multiple ligands. Pose results show atom pairs and measured
distances. With typed donor and acceptor sites but no hydrogen coordinates,
hydrogen-bond and water-bridge candidates are marked as distance-only evidence.

The classifier uses MN-Ligand's geometric rules for hydrophobic and salt-bridge
contacts. MD hydrogen bonds and water bridges use bonded hydrogens with
Wernet–Nilsson geometry. Waters contribute only when their coordinates are
present. Salt bridges require ligand formal charges from the selected topology,
a verified match to the linked study ligand, or an explicit atom-name charge
map. The report records which method was available for each interaction class.
The visual pipeline runs where the
trajectory already resides; remote trajectory transfer or automatic relocation
of the graph coordinator is not implemented. Standalone MD Analysis retains its
existing authenticated owner-workstation submission route.

## Structure identity

For legacy imports, the original selected complex is resolved by the recorded
selected-complex ID and preparation-source SHA-256. A trajectory's topology
snapshot uses its manifest path and SHA-256. The source manifest must match the SQL
study, pose, compound, accession and preparation provenance. New MD imports retain
the manifest SHA-256 in SQL as well as each report's SHA-256.
This additive metadata does not force completed historical jobs to re-ingest.

For native bound-pose MD, the original SQL structure and retained pose digest are
used. New ingestion records the MD summary path/digest when that artifact exists.

**A paired topology is not necessarily the production-start coordinates.** In
particular, a continued production segment can start from a checkpoint even when
it retains the original prepared-system topology. The UI labels these roles
explicitly. The historical imports do not uniformly archive a separately identified
production-start / RMSD-reference frame; no frame is invented or substituted.

## Artifact access and limitations

Artifact access rechecks current SQL visibility and that the selected observation
still exists with unchanged context. Local runs additionally require matching
catalog ownership/project and job access. Remote requests use the normal
authenticated, source-qualified cluster artifact API. Reports and snapshots are
opened on demand; Refresh checks pending remote requests without duplicating them.

Legacy MD observations with a recorded PDB accession but no canonical target ID
remain discoverable at that PDB's top level. They are not assigned to an exact
sequence/construct folder without evidence, and this fallback never overrides an
existing canonical target link.

An optional database-local `md-browser-cache` holds offline copies of small
artifacts, keyed by source node, task group, run and relative path. Each copy has
source metadata and a SHA-256 checked on read; a retained SQL digest is also checked
when present. Cache files are not an authorization mechanism: deletion or revoked
SQL visibility still prevents access. Old unversioned source manifests have no
retroactively invented SQL digest.

Only JSON, CSV, PDB, CIF and SDF artifacts up to 32 MiB are accepted. Absolute paths,
traversal, escaping symlinks, changed hashes and trajectory formats are rejected.
An unavailable artifact leaves SQL measurements accessible and shows an explicit
notice; a different pose or analysis is never silently substituted.

## Verification (2026-11-14)

- Regression command: `python -m pytest -q tests/test_md*.py tests/test_database*.py tests/test_data_browser.py tests/test_molecular_database*.py tests/test_bound_docking.py tests/test_target_registration.py tests/test_legacy_md*.py tests/test_entity*.py --tb=short`.
- Result: **484 passed, 3 skipped**. The skips require the optional MDTraj runtime,
  which is not installed in the host test environment. No MD simulation was rerun.
- The isolated preview's real SQL data contains **47 MD studies**. Read-only
  retrieval from PC004033 produced **771 verified small artifacts**: 675 JSON
  reports/manifests and 96 structure snapshots. All 665 distinct SQL summary
  references were recovered without a digest mismatch or source-export skip;
  279 production reports contain the stored rich structural-dynamics fields.
- Real-data browser checks cover study/replica navigation, contact curves,
  endpoint provenance, and original selected-pose viewing. Both the original
  selected complex and paired topology were independently digest-checked.
- Preview verification uses the isolated, digest-checked offline cache. The new
  UI's authenticated cross-workstation transfer path is covered by the existing
  artifact client contract/tests, but was not exercised end-to-end against the
  production API during this preview. The production deployment, SQL records and
  trajectory storage were not changed.

### Executable pipeline checks (2026-11-14)

- Extended regression selection: **571 passed, 3 skipped, 2 deselected**. The two
  previously failing candidate-handoff tests (`test_execution_reuses_only_unchanged_steps_and_preserves_prior_results`
  and `test_parallel_branches_collect_distinct_candidate_ids`) were excluded;
  this is not a claim that the complete application suite passes.
- Final MD pipeline, source-picker, CPU-submission and adapter-signature checks:
  **19 passed**. Scheduler-child tests cover successful calculation collection,
  failure blocking, input types, project/access checks, exact topology/segment
  identity and tampered numerical outputs. The scheduler is mocked in these tests.
- The existing `mn-md:cu128` container ran all **12 MD runtime tests**, including
  the three skipped on the host. Synthetic trajectories check Å/ns units,
  alignment, population RMSF, mass-weighted Rg, DCD timing, PBC handling and
  changed-source hash rejection. No retained research simulation was rerun.
- Pipeline frontend build, Python compilation and `git diff --check` succeeded.
  The three-step draft is visible in the isolated development preview on port
  8514. Production services, SQL data and remote application installations were
  not changed; owner-affinity execution of a real retained trajectory remains
  a deployment/integration check, not an exercised production workflow.

### Ligand–residue diagram checks (2026-11-14)

- `python -m pytest -q tests/test_md*.py tests/test_data_browser.py tests/test_database_browser.py --tb=short`:
  **183 passed, 3 skipped** (optional host MDTraj runtime).
- Diagram tests cover chemical depiction, recorded interaction types, exact
  occupancy thresholds, deterministic truncation, chain-qualified residues,
  invalid values, escaped labels, missing chemistry and geometric-only analyses.
  App tests check display options, the enlarged dialog and that switching replicas
  changes the diagram without mixing source reports.
- The isolated preview renders the saved ZINCtv000001S7sB / TRα / 4LNW / replica 2
  production analysis from verified imported artifacts. No simulation, interaction
  recalculation, SQL migration or production deployment is performed.

### Named ligand-atom endpoints (2026-11-14)

- Extended MD/browser regression selection: **193 passed, 3 skipped** (host MDTraj
  unavailable). Atom tests cover source-name preservation under atom reordering,
  exact SVG endpoint coordinates, chemistry/stereochemistry mismatches, duplicate
  names, missing names, chain-qualified residues, exact replica selection and
  production/preparation topology-pairing contradictions. No nearest/arbitrary
  depicted atom is substituted for unavailable evidence.
- A read-only check of one longest stored production per retained MD study found
  usable named-atom mappings for **46 of 47 studies**, with **653 recorded atom
  names matched and none unmatched**. This is a study-representative check, not
  a claim to have inspected every production segment. The remaining 4LNW reference
  control lacks a selected ligand residue name in its sampled report; it stays
  explicitly unmapped rather than being assigned a guessed ligand.
- Original source code was checked: legacy `top_ligand_atom` can aggregate nearest
  contacts and other interaction counts; typed edge weights remain residue-level.
  Neither the revised diagram nor its assignment table asserts exact per-type
  atom-pair occupancies. Those would require new trajectory calculations.

### Compact multi-replica overview (2026-11-14)

- Regression selection `python -m pytest -q tests/test_md*.py tests/test_data_browser.py tests/test_database_browser.py --tb=short`:
  **218 passed, 3 skipped** (the optional host MDTraj runtime). Python compilation
  and `git diff --check` also succeeded. Tests cover all three maps, exact replica
  matrices, missing/invalid occupancies, report-source validation, failed endpoint
  exclusion, the wide overview and per-panel enlargement.
- The reference was MN-Ligand's `_render_md_plot_matrix`, including parallel
  replica panels and percentage-labelled interaction maps, not just its raw
  contact tables. The new shared viewer reconstructs that compact arrangement
  from imported aggregate reports in both MD Results and Data browser.
- Up to 19 plot kinds are supported from the saved report: interaction maps,
  recorded RMSD/displacement/distance/Rg/secondary-structure curves, per-type
  occupancies, RMSF, contact matrices and scalar endpoint/retention/throughput
  panels. This does **not** claim complete parity with the original 29-panel
  implementation: additional BB/SC composition, interaction-composition and
  persistence/interface-network presentations, mean/SD summary styling, layout
  ordering and combined PNG ZIP/workbook exports remain outside this change.
- Data browser requires every overview replica's source production ID and
  replica number to match an accessible production in the selected study. Results
  checks recorded source IDs against its study manifest. Missing matrix values
  are not filled from consensus or another replica, failed replicas/endpoint
  values are not presented as successful plots, and explicitly inapplicable
  salt-bridge analyses are not offered as measured occupancy panels.
- The verified ZINCtv000001S7sB / TRα / 4LNW report spans 0–98.97 ns and contains
  three completed source replicas. The compact maps preserve their different
  SER277 hydrogen-bond occupancies (38.33%, 50.24%, 65.64%) and water-bridge
  occupancies (17.64%, 9%, 1.85%). All three are shown together in the isolated
  development preview, including percentages on the actual coloured lines.
- No trajectories were loaded, analyses recalculated, SQL observations rewritten
  or production services restarted. The preview SQL file's size and modification
  timestamp remained unchanged. Its authenticated session was preserved.

### Attachment-aware interaction layout (2026-11-14)

- Both the compact maps and single-segment diagram use the attachment-aware
  `distribute_2d_labels` geometry from MN-Ligand's static interaction viewer:
  ligand-atom angular ordering, equal **arc-length** spacing around an ellipse,
  and deterministic phase selection. Hotspot score still selects which residues
  to show; it no longer determines their positions around the ligand.
- A bounded slot-assignment refinement reduces straight leader crossings without
  moving ligand atoms or changing perimeter spacing. Compact maps additionally
  fan different interaction types into separate routes and place their percentage
  labels against actual rendered bounds, avoiding other labels, residue bubbles
  and the central molecule where possible. Positions remain consistent across
  replicas. The layout is a schematic, not a reconstructed three-dimensional site.
- The displayed ZINCtv000001S7sB / TRα / 4LNW example previously had 25 straight
  residue-leader crossings in score order. The source spacing routine reduces
  that to 3; the refinement reduces it to 0. Sampling each rendered curve at
  65 points finds no proper crossings between different-atom/different-residue
  routes in any of its three replica maps (16, 17 and 14 interaction lines), and
  no percentage-label bounding-box overlaps. Shared atom endpoints are intentional
  and excluded from that crossing count. This does not guarantee a crossing-free
  embedding for every ligand, graph or residue limit.
- Regression selection: **235 passed, 3 skipped** (optional host MDTraj runtime).
  New tests cover perimeter spacing, deterministic routing, shared atoms, stable
  replica positions, percentage-label collisions, invalid geometry and the
  representative 25→0 crossing reduction. Compilation and `git diff --check`
  succeed. Only presentation geometry changes; saved values, source atom names,
  SQL data and production services are unchanged.
