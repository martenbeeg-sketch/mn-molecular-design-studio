# Complex Analysis coverage: accepted R3 and Tau IDR designs

Audit date: 2026-11-03. This covers accepted design records registered in the Database, not rejected or unresolved candidates.

## Result

All **537 accepted structures** now have a completed Biotite interface-geometry analysis on the exact saved pose. The analysis context resolves to the expected target identifier and stored collection/model for **537/537** records, with target chain A and binder chain B. No acceptance decisions changed.

| Target | Engine | Accepted structures | Biotite | ipSAE |
|---|---|---:|---:|---|
| R3 | BindCraft2 | 100 | 100 completed | 100 measured |
| R3 | RFD4 | 306 | 306 completed | unavailable: no compatible pairwise PAE |
| Tau IDR | BindCraft2 | 25 | 25 completed | 25 measured |
| Tau IDR | RFD4 | 106 | 106 completed | unavailable: no compatible pairwise PAE |
| **Total** |  | **537** | **537 completed** | **125 measured; 412 unavailable** |

That is 406 R3 and 131 Tau IDR accepted structures. The inventory contains 537 distinct saved binder sequences. Every Biotite observation is attached to the exact Database molecule, accepted collection, source model, target ID, and chain assignment.

## Measurements

Biotite used a 5.0 Å heavy-atom contact cutoff. Across the 537 poses, heavy-atom contact pairs ranged from 366 to 1,744, residue-pair contacts from 47 to 217, binder interface residues from 20 to 73, and target interface residues from 15 to 60. Every pose had a nonzero interface-contact count. These describe geometry in one predicted pose; they are not binding-affinity predictions or experimental evidence.

ipSAE is available for all 125 BindCraft2 structures with compatible saved PAE: 100 existing R3 results, one existing IDR result, and 24 newly run IDR analyses. All 125 have numeric scores and completed observations attached to the exact structure, target ID, and chain assignment. Scores span 0.127386–0.466779. The 24 new jobs completed on the source workstations: 8 on PC004032 and 16 on PC004033.

The 412 accepted RFD4 structures do not have compatible pairwise PAE in the saved evidence. ipSAE is therefore unavailable for these records, not failed and not an acceptance rejection. Biotite results are available for all 412.

## Workstation visibility

The primary Database contains all 537 Biotite observations and all 125 BindCraft2 ipSAE measurements. The sharing flag had been disabled on the primary and both source workstations, which prevented the IDR observations from appearing centrally. Sharing is now enabled on all three, and the 25 IDR BindCraft2 results have replicated into the primary. Each score was verified against its exact structure digest, source collection/run, analysis run, and target ID. One accepted design has two matching ipSAE observations from separate analysis runs; it remains one distinct design. The agents continue synchronizing the remaining historical Database metadata.

No BindCraft2 or RFD4 filter was rerun, and no pass/fail state was changed. Biotite and ipSAE are additional measurements, not new acceptance criteria.

`accepted-complex-analysis.csv` contains one row per accepted structure with the source collection/model IDs, exact structure digest, target/chain link, Biotite measurements, and ipSAE result or availability status.
