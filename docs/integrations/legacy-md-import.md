# MN-Ligand MD import — 23 October 2026

## Selected data

The import publishes **completed 100 ns studies only**, including their required
preparation, initial production segments and continuation history. Standalone
short runs are excluded. A requested duration alone is insufficient: completed
production durations or recorded steps and timestep must establish the total.

| Project suffix | Studies | Completed 100 ns replicas | Retrospective jobs |
| --- | ---: | ---: | ---: |
| Medchem | 24 | 72 | 3 |
| Large library screen | 20 | 60 | 3 |
| MD references | 3 | 9 | 4 |
| Total | 47 | 141 | 10 |

The full project names retain the existing `Thyroid Hormon Recptor interaction`
prefix. Each job groups one dataset and tool; targets and replicas remain
separate studies and records inside that job. The reference project includes
both OpenMM and GROMACS. Compound studies derive from the original Medchem or
large-library selection; overlap with HY-L126 does not change their source
dataset.

The reviewed plan contains 279 production records, 292 MMGBSA results,
94 aggregate analyses and 759 source directories. It preserves 1 + 99 ns and
50 + 50 ns histories without concatenating trajectories or double-counting
an appended file.

## Storage and execution

Data owner: **PC004033**.

- Source: `/data/mn-ligand/workdir/runs`
- Destination: `/data/results/mn-molecular-design-studio/workdir/runs`
- Normal job groups: `ligand-md`, `ligand-mmgbsa`, `md-analysis`
- Exact source copies: `<group>/<job>/artifacts/legacy/<original-run-id>/`
- Normalized study manifest: `<group>/<job>/artifacts/md/study.json`

The source is copied on PC004033 into independently stored files. Each copied
file receives SHA-256 verification and destination readback. Conflicts fail;
retries verify existing copies. No symlink or hardlink is used for the import.
The bulk trajectories remain on their owning workstation under Studio's storage
policy; other workstations receive catalog and molecular-database metadata.
Small selected Results summaries can be fetched through authenticated artifact
requests. A `cluster://` reference identifies the **copied Studio artifact**,
not the old source directory.

Copying is staged under `_mn_ligand_md_migration/staging`. The reviewed import is
approximately 852 GiB. Jobs are published after the copy finishes and the
registration validation gate matches the plan and code hashes. Publication uses
the normal catalog, preserves original source dates, and backs up the owner
molecular database before adding observations.

The independent user service is `studio-legacy-md-import-20261023`. Its working
files are under:

`~/.local/share/mn-molecular-design-studio/validation/md-import-20261023/`

Read status on PC004033:

```bash
systemctl --user status studio-legacy-md-import-20261023
cat /data/results/mn-molecular-design-studio/workdir/runs/_mn_ligand_md_migration/progress.json
journalctl --user -u studio-legacy-md-import-20261023 -n 30
```

The status changes from `copy` to `copy_complete`, then `registered`.
`registration-report.json` records observation counts, withheld reasons and
trajectory snapshot warnings. The controller stops on verification or
registration errors. The service runs independently of the initiating SSH
connection; a workstation reboot needs a new service start.

## Scientific provenance

- Preserve the exact selected complex and pose; rank is not assumed to mean
  “best”. Manual selections remain manual selections.
- Preserve prepared-target identifiers and preparation compatibility contracts.
- Preserve modeled ligand chemistry separately from source-dataset identity.
- Reuse existing canonical molecule identities. New writable MD evaluation
  groups hold observations; peer-owned source dataset groups are not edited.
- Keep GB and PB observations separate. Retain original units, method and
  sampling settings; convert supported energy terms once to kcal/mol.
- Compare endpoint-recorded trajectory hashes with copied production hashes.
  Unknown or mismatched snapshots remain explicit provenance warnings.
- Retain raw summaries, interaction maps, checkpoints and analysis artifacts.
  Large residue matrices stay in Results rather than becoming database columns.

## MD Analysis tile

`ligand.md.analysis` analyzes one selected completed production trajectory with
its matching topology. It runs on the data owner and mounts source artifacts
read-only. Outputs include protein-aligned backbone RMSD, ligand RMSD,
geometric heavy-atom residue contact occupancy, CSV tables and a representative
PDB. The analysis records source hashes, selections, time window and PBC method.

DCD timing is read from a validated CHARMM CORD header, including the first-step
offset. Native XTC timestamps are retained. Time windows therefore refer to the
actual stored time axis; they are not reset to zero. PBC imaging requires bond
connectivity and fails explicitly when it cannot be applied safely.

The Results viewer also displays the original recorded RMSF, gyration,
secondary-structure, retention, contact, bridge and network data. Recalculating
the richer historical hydrogen-bond/water-bridge/salt-bridge analysis is not yet
implemented. Geometric contacts are labeled accordingly. Separate segments are
not silently stitched into a new 100 ns trajectory, and the linked historical
stages are not a new executable visual pipeline.

Cross-workstation submission still enforces Studio's code-revision compatibility
check. These workstations also have unrelated development differences; owner-local
submission on PC004033 avoids that mismatch. The MD importer does not bypass
the compatibility check or synchronize unrelated source changes.

## Validation

- Exact path preflight: **279/279** trajectory/topology pairs resolved.
- Header preflight: **273** valid DCDs and **6** XTCs.
- Scientific parser preflight: **292/292** endpoints support both GB and PB;
  **94/94** aggregate analyses have three successful replicas. No parser
  withholding or recognized unit ambiguity.
- Isolated registration smoke on PC004033: existing molecular identity reused,
  separate GB/PB observations, 1 + 99 ns lineage preserved, repeat ingestion
  adds zero observations.
- Actual bundled MD container fixtures: **12 passed**, including DCD timing,
  chunk/stride windows and periodic imaging.
- Focused migration, storage, database, Results and authorization checks:
  **70 passed** locally. Three MDTraj cases skipped in that environment were
  exercised successfully in the actual MD container above.

Full scientific trajectories have not been recalculated as part of migration.
