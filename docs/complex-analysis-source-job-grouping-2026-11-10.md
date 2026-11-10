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
