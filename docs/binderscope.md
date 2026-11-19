# Integrated BinderScope

## Live analysis workspace versus preview

A read-only preview (`MDS_READ_ONLY_PREVIEW=1`, or the intentionally disconnected
API origin `http://127.0.0.1:1`) cannot queue calculations. Removing only the
disabled-button condition does not make an isolated database copy executable.
Use a pinned, tested source release for an execution-enabled UI and API, with
production authentication, the real Studio home and the normal scheduler.
`MDS_ANALYSIS_WORKSPACE=1` labels that UI as live and initially opens Data Browser;
it does not bypass authentication, grant access, change data roots, or start jobs.
Keep the old preview snapshot for recovery, and never copy its accounts or SQL
contents over live data. Saved-cohort BinderScope runs use the authenticated local
queue, retain exact source/project/checksum validation, and need the matching CPU
image. Final submission remains an explicit user action.

The BinderScope tile runs the actual MIT-licensed analysis source in
`mn_molecular_design_studio/vendor/binderscope`, not a service or code link to a
separate checkout. `PROVENANCE.md` records the imported revision. Studio's
existing interface/liability/tag tools remain separate methods; their numbers
are not relabelled as BinderScope measurements.

Build the CPU runtime from the repository root:

```
docker build -f containers/binderscope/Dockerfile -t mn-binderscope:c61ef9c .
```

Only the base image and ordinary Python dependencies need downloading at build
time. BinderScope source and the Studio runtime adapter are copied from this repository, not mounted from a separate checkout. Analysis containers run
with networking disabled. No BinderScope launcher, server, external checkout or
BindCraft installation is used. No BindCraft engine code is redistributed.

## Launch and attachments

In Data Browser, select any tree level and use the table's row checkboxes, then
**Analyze · BinderScope**. The default scope is the selected rows on the current
page; **All matching designs** includes every page and the current folder's
refinements. The existing tree, table and right-hand viewer stay in place.
An in-browser confirmation dialog checks coordinates, checksums and recorded
roles only when explicitly opened, never during normal tree navigation.
Missing or ambiguous roles exclude a design; alphabetical chain order is not
used. A scope spanning projects queues a separate authorized job per project,
not a mixed-access collection. The dialog lists exclusions and existing exact-
structure assessments. It never queues work merely because a folder was opened.
The confirmation table also has an explicit **Include** checkbox per design,
plus **Select all** and **Clear selection**. Unchecked designs are not submitted;
clearing the selection disables submission. Changing analysis settings does not
change the chosen design identities or their admission status.
Access and checksums are checked again during API preparation and execution.
Replicated designs can be analyzed on a workstation with an accessible, verified
copy of their exact saved complex. The result is authored in a new local
Evaluation collection, not in the source workstation's collection. Owner-qualified
structure references preserve the attachment when evaluations are synchronized.
This does not readmit the design or alter its original filter decision.

For older binder-only SQL identities, an explicit accessible SQL target link can
supply the target sequence. Both binder and target must match unique coordinate
chains by exact sequence; missing, conflicting or ambiguous evidence is excluded.
The API rechecks this evidence before staging. No source workdir or alphabetical
chain ordering is needed. Native mmCIF files without optional occupancy or
B-factor columns are normalized on a copy; absent confidence remains unknown.

Saved-mode analysis requires one exact protein binder chain and one target chain.
RFD4's generated design and native-oracle refold are separate representations of
one saved design. The adapter verifies the exact native filename, molecular and
chain identities, and candidate-specific confidence receipt before choosing the
refolded complex. Both coordinate files are copied into the portable Database
store; the refolded representation is preferred for evaluation and default
BinderScope selection. Historical recovery is additive and preserves the original
observations and admission decisions. A generated-only legacy record is not
silently substituted when refolding evidence exists: default analysis explains
the missing validation coordinates. Explicit generated-model selection remains
available for design-stage geometry, without borrowing refolded confidence.
Database-owned coordinates are content-addressed files indexed by SQL, not links
to the generator's workdir. Portable peer copies retain their checksums and remain
usable after the source workdir is unavailable. A generated source model and the
selected refold are labelled separately in score history; an ancestry link does
not mean their coordinates or measurements are interchangeable.
Topology defaults come from recorded binder format/task settings, not sequence
length or the target's modality. Confidence defaults require an explicit encoding
declaration on the exact representation or its linked Evaluation for the same
binder chain; an average
pLDDT, generator name or different model is insufficient. Conflicts remain unknown.
The dialog provides per-design overrides and optional bulk declarations for
unknown values; mixed datasets do not share one global topology/encoding.
Settings and their evidence are part of the authorized collection identity.
Unknown confidence leaves geometry, contacts, patches, liabilities, structural
comparisons and binding modes available. Its B-factors are not exposed as pLDDT,
and confidence-dependent terminal assessments are unavailable. Unknown topology
also disables terminal assessment; cyclic binders explicitly have no free termini.
Declared pLDDT outside 0–100 is downgraded to uninterpreted with a recorded reason.
These declarations concern the binder's Cα B-factors. Target B-factors remain
uninterpreted, including when the binder has confirmed pLDDT.
Different exact target sequences are evaluated in separate cohorts; topology
subgroups retain their own native format rules within a cohort. Source
coordinates are never overwritten; normalized coordinates retain a
canonical-to-author chain/residue/insertion-code map and their own checksum.

Each saved result attaches to its existing authorized SQL molecular identity and
exact source representation using the normal Evaluation ingestion contract.
It is post-analysis, not an additional admitted design. Original generator
scores, filter pass and human decisions are unchanged. The job retains the full
cohort report and downloadable archive; source-bound per-design details remain
available in Data Browser. Shared cohort matrices are not duplicated into every
SQL Evaluation.

## Native Studio inspector

Select a saved design in Data Browser, open **Scores & analyses → BinderScope**,
then **Open BinderScope inspector · 3D & interface**. This opens a wide inspection
dialog without replacing the browser tree or changing the selected design. The
Results job view uses the same inspector. The presentation is Studio-owned code,
not an embedded upstream report, external BinderScope server or link to another
checkout. The linked view uses a locally bundled JavaScript/3Dmol frontend inside
an isolated Streamlit component. Plotly/Mol* remain the fallback and detailed
measurement views. Scientific measurements still come from the bundled,
checksum-verified upstream analysis engine.

- **Interface · 2D & 3D:** a residue-pair contact map and a crossing-reduced
  bipartite residue network. Interaction types retain overlapping H-bond/salt
  assignments, marker sizes show atom-pair counts, and mapped intended hotspots
  are marked. Filters change the display, not the stored measurements. All pairs
  remain available in the contact table and CSV export.
- The linked inspector places cartoon 3D, sequence tracks and a crossing-reduced
  amino-acid-labelled contact diagram together. Hover and click link these views
  locally, without a Streamlit rerun. Its **Findings**, **Confidence** and
  **Footprint** tabs are local views of the same saved assessment. The footprint
  colors the target surface by recorded burial for **one complex**, not campaign-
  wide contact frequency. Per-residue burial is not an energy or alanine scan.
- Hydrogen-bond overlays use the recorded donor/acceptor atom coordinates;
  salt bridges use the upstream charged-side-chain group centroids. Missing
  endpoints are omitted and counted, never replaced with Cα coordinates.
  An orange Cα connector is only a selection guide, not a chemical bond.
  Protein chains remain cartoons, including when residues are highlighted.
  Confirmed confidence coloring is restricted to the binder. The component
  contains its own assets, blocks network connections, and exports a structure
  image and the contact diagram. The original assessed coordinate file and full
  contact CSV remain downloadable outside the component.
- Clicking a contact marker, or choosing a pair, highlights the corresponding
  residues in the verified saved complex. Protein chains remain cartoons.
  Chain-role, binder pLDDT, hydrophobicity and secondary-structure coloring are
  available; pLDDT coloring requires confirmed binder encoding for that model.
  Target B-factors are not relabelled as confidence.
- **Residue profile:** binder burial and confirmed confidence tracks, secondary
  structure, exposure, interface membership and liabilities with CSV export.
- **Chemistry & liabilities:** atom-labelled hydrogen bonds, salt bridges,
  buried unpaired polar atoms, exposed hydrophobic patches and sequence flags.
- **Hotspots, target & tags:** original-to-saved hotspot mapping and distances,
  supported target-reference comparisons, and terminal tag/fusion probes.
  Unavailable evidence remains unavailable, not a successful check.
- **Parameters & evidence:** exact settings, immutable original source inventory,
  residue mapping and downloadable per-design JSON.

Plot labels and 3D selections use saved author chains, residue numbers and
insertion codes, not BinderScope's internal renumbering. SQL observation bindings
resolve the local representation after peer replay; the original embedded source
IDs remain in the evidence. Coordinates are authorized again and checksum-bound
to the analysis before display. A deleted/unresolved association or changed model
does not fall back to another structure. The 3D inspector is only loaded on an
explicit inspection action, not during ordinary tree expansion.

## Parity boundary

### Original targets and intended hotspots

The browser recovers per-design hotspot declarations from member-specific native
input settings and evaluations bound to the exact source representation. It
does not borrow hotspot settings from another design, model or pipeline branch.
Conflicting declarations are retained as unassessed evidence. Original reference
coordinates must be a readable, checksum-verified SQL target reference within
the design's project boundary; arbitrary workdir paths are not opened as proof.

Re-prediction can rename chains or renumber residues. Original hotspot residue
identities are therefore mapped by the target sequence to the exact saved
complex, then translated to BinderScope's internal numbering. Whole-chain exact
matches and unique contiguous crops are supported, including author insertion
codes. Repeated/ambiguous segments, mutations, noncontiguous crops and absent
residues are left **unmapped**, not guessed from a numeric offset. Explicit
`hotspot_numbering=saved_target_author` can instead identify already mapped
saved-structure declarations.

Each attached result shows requested, contacted, missed and unmapped residues,
both original and saved labels, and minimum heavy-atom distances using the
upstream 4.5 Å contact rule. Full intended-hotspot coverage is only reported
when every requested residue maps. An absent declaration is not a pass. This
post-analysis never changes the original generator's filter/admission decision;
those filters may use different definitions or may not require hotspot contact.

Where the original target can be mapped, its Cα conformation is compared with
the saved target after a proper rigid-body fit, reporting the mapped scope and
RMSD. For sequence-only/disordered targets, fixed-fold preservation is explicitly
not applicable; hotspot contacts can still describe the saved conformation when
their residue mapping is supported. Full-length reference placement, partners,
specificity and scaffold checks remain separate native-mode capabilities.

The implementation uses the same upstream calculations, not approximations
claimed to be equivalent. Available evidence nevertheless depends on input:

| Analysis | Saved SQL designs | Native BindCraft campaigns |
| --- | --- | --- |
| Contacts, SASA, H-bond geometry, salt bridges, buried polar flags, patches, secondary structure and liabilities | Actual upstream engine | Actual upstream engine |
| Binding modes and structural comparisons | Selected exact-target cohort only | Upstream campaign analysis |
| Intended hotspots | Source-specific declarations, exact residue mapping and upstream contact distances; otherwise explicitly unassessed | Original run hotspot definitions |
| Original input target conformation | Verified SQL target reference, exact sequence/unique crop, mapped Cα RMSD; not a fixed-fold verdict for IDRs | Supplied input target and upstream placement |
| Native scores and filtering | Studio originals retained separately; admission is not inferred from BinderScope | Upstream readers retain original campaign scores and failures |
| Tag/fusion assessment | Recorded/declared single-chain topology and confirmed confidence encoding; otherwise explicitly unavailable | Upstream campaign format rules |
| Full-reference placement/context | Not automatically joined from SQL; use explicit native mode for now | Supplied reference, chain and partner selections |
| Specificity | Not inferred from saved designs | Requires supplied profile/accession and offline reference data |
| Scaffold/paratope | Unavailable without a native scaffold edit specification | Upstream campaign edit specification |
| Rejected designs, near misses, funnel and stopping advice | Cannot be established from an admitted-only SQL subset | Upstream readers, with their documented version/modality limitations |

Thus this is integrated **engine parity**, not a claim that an admitted-only SQL
subset contains all information of a native campaign. Unsupported/missing inputs
must not be presented as successful checks, zero findings or experimental
validation. Native mode returns the upstream validated export and offline review
pages in the archive; it does not automatically promote native designs to SQL.
Original campaign inputs are available under the tile's **Advanced** expander,
not as a competing source selector in the normal browser flow.
ConSurf is not enabled. There is no new universal ranking or automatic filtering.

The focused tests include an explicit Docker smoke test:

```
MDS_TEST_BINDERSCOPE_IMAGE=1 python -m pytest -q tests/test_binderscope.py
```

Payload/identity/security regressions are in `tests/test_binderscope_view.py`.
Client-side DOM interaction tests are separate from actual WebGL/browser testing:

```
npm --prefix tests/frontend install --ignore-scripts
npm --prefix tests/frontend test
```

They test linked selection, local tabs, overlays, confidence scope and safe text
rendering with a stubbed 3D renderer. They do not establish browser WebGL support.
