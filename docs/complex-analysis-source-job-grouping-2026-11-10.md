# Complex Analysis: source jobs, campaigns and scientific counts

## Correction

The Database selection view groups recorded campaign memberships first, then
reuses Results' historical batch grouping as explicitly labelled presentation.
It does not merge independent runs merely because their target is the same.
Saved registration cohorts remain cohorts, not invented parent campaigns.
Project/user boundaries and source-qualified job identities are preserved.

The R3 `10-design batch` entries dated 2026-10-30 comprise 18 source jobs across
pc002975, PC004032 and PC004033. Checked source metadata/catalog records have
no common parent run ID. They appear in one historical batch group with the
individual source jobs accessible in the source breakdown; this does not
create or repair a scientific parent relationship. The native target ID and
Database target-molecule ID are separate references: their equivalence is not
inferred from labels.

## Counts and input targets

The main table separates source attempts, native acceptance, recorded final
acceptance, Database-retained designs and retained models. Missing totals are
unknown, not zero. Parent aggregate counts are not added to worker totals.
Evaluation coverage is no longer presented as a campaign pass rate.

Read-only source evidence for PC004032 `design/20261021-184403-b0c849c0`:

- Native `.campaign_state.json`: 3,665 trajectories; 6,614 redesigned candidates
  scored; 6,525 redesigned candidates rejected; 35 accepted/ranked designs.
- These are different populations; no combined pass-rate denominator is inferred.
- Its older normalized `result.json` reports zero trajectory counts; both source
  facts remain preserved, with the native state explicitly authoritative here.
- Original input FASTA: `Tau_2N4R_R3_residues_306-336`,
  `VQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQ` (31 residues).
- Native states for the other inspected Tau source jobs report 5,000/65 and
  1,638/8 trajectories/ranked designs on PC004032, and 3,874/17 on PC004033.

Completion/explicit indexing captures the existing native campaign-state
reader's typed counts and the exact in-job target FASTA, with source references
and SHA-256 evidence. Input sequences remain distinct from model-specific crops;
FASTA does not establish coordinate chain assignments. Saved target context is
also read from the exact producing Evaluation, not unrelated downstream analyses.

## Verification and scope

The 8510 preview uses a 2026-11-10 snapshot of the actual pc002975 Database,
with dispatch and peer synchronization disabled. Only that copy received the
four explicit Tau source summaries. No production Database, Results, scientific
outputs or protected campaign-settings files were modified.

Browser checks confirmed the single 18-job batch row and membership summary,
and the Tau 35-design selection with 3,665 attempts and 35 retained models.
Coordinates are not loaded during selection. The source-context metadata query
for 439 authorized Database groups took 0.309 seconds in this staging copy.

Other historical runs require explicit source summary publication to show
unindexed totals/target correspondence; ordinary navigation never scans Results
or contacts remote workstations. No deployment or push was performed.

Focused regressions: 90 passed across collection labels, authorized catalog,
scientific lineage, Database hydration and collection evidence tests. These
cover explicit versus historical grouping, stage/count units, missing versus
zero values, source-qualified lookup and ACL rejection, original versus cropped
target context, idempotent capture and rejection of out-of-job FASTA paths.

## Source-count repair and legacy scientific collections

The follow-up recovered 299 RFD4 source summaries into the isolated preview,
using an explicit manifest of Database-linked Michele_tau jobs. Reads used the
recorded local roots and SSH as `compute` on PC004032/PC004033. Only metadata,
normalized candidate records and declared metric CSVs were inspected; no
coordinates, predictions or scientific calculations were generated. Full
rejected-candidate records remain excluded from replicated publications.

The existing `molecular_library.sqlite3` is the correct relational home:

- `scientific_jobs`: source-qualified job identity, scope, typed counters and
  source references/digests in `facts`; publication remains owner-authored.
- `scientific_execution_edges` / `scientific_campaign_jobs`: explicit execution
  relationships and campaign membership; no name-inferred parent identity.
- `scientific_campaigns`: derived campaign summaries.
- Existing molecular, structure and Evaluation tables: retained scientific
  records, unchanged by source-count repair. No new Database/schema was needed.

RFD4 capture now counts unique native per-example generation rows, independently
of requested designs or retained structures. A row missing complete native PPI
metrics is unresolved, even if `ppi.pass` is zero. Duplicate aggregate/part rows
are counted once; conflicting decisions remain unresolved. Counts carry their
method/stage, population units, declared CSV references and SHA-256 evidence.
These are **scored generations**, excluding launches/sampling failures that
produced no metric row; they are not presented as exhaustive attempts.

Distributed-controller counters retain their original units and campaign scope.
A recorded parent total replaces its workers' totals only with explicit shared
campaign membership. Final pass counts require the saved `validated` acceptance
level. Requested budgets, native acceptance and final decisions remain distinct.
The validated local controller `design/20261101-104153-71432d5e` records 6,700
started design slots / returned candidates and 100 final validated candidates.
The earlier group display name's requested 200 passes is not used as an outcome.

Real preview evidence:

- R3 18-job historical batch group: **180 scored generations**, 174 retained
  models, 180 incomplete native PPI decisions, no recorded final decision.
- Accepted-registration cohort: 72 generations recovered from 72 of 73 source
  jobs; the unavailable historical source remains explicit rather than zero.
- Four source references remain unresolved: historical source directory
  `20261030-014030-a12be5a3` is unavailable, and the referenced source metadata
  for `20261101-231254-aa0b147a`, `20261101-234700-34442405` and
  `20261102-081647-0ff0921a` differs in owner/project from the requesting Database
  groups. No ownership reassignment or automatic cross-scope mapping was made.
- Capturing both a retained design and a supplemental retained target exposed a
  publication validation bug. The summary now distinguishes published retained
  links from candidate/model counts; targets never become generated designs.
- Completion capture records the local configured workstation label alongside
  its unchanged source node ID. Historical manifests must supply the verified
  source label when indexing an offline copy.

The user's additional data were present in the same snapshot. The selector was
also excluding imported/evaluation groups with real saved evidence. It now reads
all authorized group categories and includes saved imported/evaluated structural
collections without inventing source runs or computational acceptance. Import
collections open in `All designs`, not `Collected passing designs`.

| Saved data | Actual Database location / limitation |
| --- | --- |
| SpiroSpace | Project **Spirochem library test**, two final shortlists of 1,000 members each, 1,000 and 1,001 structure records. Recognized as docking from its explicit screening-selection flag. |
| Piccirillo | Project **Piccirillo FABP Ligand Selectivity**, collection **Piccirillo FABP reference ligands and evidence**: 3 members, 16 structure records, 16 observations; now selectable. |
| HY-L126 | Project **Thyroid Hormon Recptor interaction HY-L126**: library 751 members, no structures; **Evaluated structures** has 748 members / 2,174 structure records and is now selectable. Legacy **Tool inputs** also contains 7,118 structure records / 10,506 observations and remains in Database; input groups are not automatically treated as curated analysis outputs. Counts overlap and are not unique compound/model totals. |
| THRalpha MD | Project **Thyroid Hormon Recptor interaction MD references**, **MD evaluations**: 79 saved observations, no structure attached to that group. Evidence is available through Database; it cannot be offered as a retained 3D model without an explicitly linked representation. |

Browser verification confirmed 180 recorded generations on the selected batch
and restored Piccirillo discovery. Metadata-only source-context lookup remained
under 0.4 seconds in the measured staging query. Canonical molecule/group/
structure counts remain 21,473 / 1,050 / 21,537. Production data/services were
unchanged; no deployment, backfill publication to peers or push was performed.
Focused regression suites: 99 passed. Reproducible read-only capture manifest,
source hashes and repair report are retained under
`/tmp/mds-scientific-work-review-20261110/rfd4-count-repair/`.
