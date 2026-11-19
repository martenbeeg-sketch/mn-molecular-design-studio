# Design creation times and display names

Names follow `Design_<generating-workstation>_<UTC-creation-time>`.
For example: `Design_pc004033_20261101T094421Z`.
Native file-write times can retain microseconds, for example
`Design_pc004032_20261020T040618.170398Z`.
Names with indistinguishable computer/time combinations receive a stable short
collision suffix. Immutable molecular, candidate and Evaluation IDs remain the
scientific identities; display names never join results.

## Creation is not launch or registration

`metadata.design_creation` records the source event using
`studio.design-creation.v1`. It includes the UTC instant, precision, actual
generating workstation/node, original artifact path and SHA-256, binder sequence
digest, and the recovery basis. The matching `metadata.studio_naming` policy is
`studio.scientific-names.v2`.

Two kinds of retained evidence are distinguished:

- `native_cif_creation`: the original RFD4 generated CIF's `_entry.date` and
  `_entry.time`. The pinned native writer uses its container clock; the reviewed
  image uses UTC. Non-UTC overrides/localtime mounts are not assumed to be UTC.
  Future RFD4 adapter launches explicitly set `TZ=UTC`. The source stamp must
  agree with its recorded execution window and exact designed chain sequence.
- `native_structure_written`: a sequence-matched original native output on the
  generating workstation. For BindCraft2 this is the unique matching
  `2_Refolded/Complexes` output identified through `!_Refolded.csv`, not the later
  rank/promotion copy. This is the original design structure's write time. The
  exact instant ProteinMPNN generated its sequence was not separately recorded,
  and the evidence explicitly retains that limitation.

File birth/ctime, copied file timestamps, aggregate CSV timestamps, run folder
names, run start/completion times and SQLite row creation times are not substituted
for per-design creation. A coordinator can read an embedded native stamp from a
copied file, but cannot assert that copy's mtime is the original creation time.
The recorded execution node takes precedence over the machine registering data.

## Future processing

The normalized-candidate writer captures creation evidence before artifact
materialization can copy files. The information is retained in candidate JSONL,
the candidate's raw metadata, SQL membership metadata and observation provenance.
Re-normalizing an existing candidate preserves its first recorded creation event.
SQL re-ingestion also preserves the existing creation event and display name.

Post-analysis/refolding can retain a source design's event only when its designed
sequence matches. A new sequence cannot inherit another design's creation time.
An adapter with no trustworthy per-design output time remains unresolved rather
than assigning the current time during import. New adapters should provide native
creation evidence when producing their output.

Target labels, source engines, candidate numbers, replicates, batches, run IDs and
original names stay in provenance and History. Creation metadata is not part of
Evaluation identity, scientific settings, filter acceptance or score calculation.

The current library correction names 807 RFD4 design events from embedded native
stamps and 125 BindCraft2 design events from original native structure-write times.
All 932 admitted designs are retained. Twenty-six RFD4 names require collision
suffixes because their embedded timestamps have only second precision.

## Historical recovery

`scripts/curate_design_creation.py` plans and applies source-owned naming updates.
It verifies the designed sequence against the original coordinates, rejects
ambiguous matches and unknown generating hosts, and never changes original
workdir files or scientific observations. Apply requires a new integrity-checked
SQLite backup, an unchanged reviewed plan and unchanged original artifacts.
Names/provenance are the only changed scientific-library fields.

```bash
python scripts/curate_design_creation.py plan \
  --root /path/to/database --runs-root /path/to/workdir/runs \
  --plan /private/new-plan.json --admin-id EXISTING_ENABLED_ADMIN_ID
python scripts/curate_design_creation.py apply \
  --root /path/to/database --runs-root /path/to/workdir/runs \
  --plan /private/new-plan.json --admin-id EXISTING_ENABLED_ADMIN_ID \
  --backup /private/new-recovery-backup.sqlite3
```

Apply on each owning workstation, then use normal source-authored replication.
Never alter a peer-owned membership locally to make labels appear synchronized.
The browser also disambiguates separately authored colliding labels after peer
replay without changing or merging their scientific identities.

Recovery backups from this correction are retained on each workstation at
`~/.local/share/mn-molecular-design-studio/backups/20261114-design-creation/library-before.sqlite3`.

Runtime activation is a separate deployment step: the future capture/naming hooks
are implemented and tested in the development source. This metadata maintenance
does not restart production services or deploy the experimental browser.
