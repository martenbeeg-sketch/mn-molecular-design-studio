# Generator route evidence from archived runs

The 26 November review used Studio `267a0bc` to read logical Runs through MCP and inspect retained output manifests. The [Capabilities checklist](capabilities-acceptance-2026-11-25.json) covers 82 combinations of route and target context. This review matched 17 combinations to an archived native generator and its next pipeline step. Twelve use the structured-target miniprotein routes and five use structured-target VHH routes.

A match requires a completed generator worker, a retained structure for every counted generator record, and a recorded candidate identity in the next completed step. The check follows identity aliases in the collection metadata. It uses each step's ledger and retained files because a pipeline with successful branches can have a Completed status alongside failed or blocked branches.

## Verified pipeline handoffs

| Route ID | Pipeline / generator worker | Generator records matched in next step | Last verified stage in the archived chain |
| --- | --- | ---: | --- |
| `rfd4-external` | `66D12` / `F53D8` | 10 | Mol-1, ipSAE and acceptance decisions for all 10 |
| `rfd4-internal` | `CE342` / `34442` and linked batches | 194 | Native refold metrics, collection and recorded filter decisions for all 194 |
| `rfd3` | `66D12` / `81CDA` | 10 | Sequence design produced a variant for each backbone |
| `rfdiffusion` | `66D12` / `5DCAE` | 10 | Sequence design produced a variant for each backbone |
| `genie3` | `66D12` / `03FF2` | 10 | Sequence design produced a variant for each backbone |
| `pxdesign` | `4983A` / `273F5` | 10 | Sequence design, Mol-1, ipSAE and acceptance decisions for all 10 |
| `protpardelle-binder` | `66D12` / `CF6B8` | 10 | Sequence design produced a variant for each backbone |
| `bindcraft2-generator` | `43B42` / `3B921` | 6 | Mol-1, ipSAE and acceptance decisions for all 6 |
| `boltzgen-generator` | `74BE2` / `1B502` | 10 | Existing sequences carried through sequence design, Mol-1, ipSAE and acceptance |
| `complexa-generator` | `2D9A1` / `0E17E` | 10 | Mol-1, ipSAE and acceptance decisions for all 10 |
| `esmfold2-protein` | `66D12` / `E3C06` | 10 | Mol-1, ipSAE and acceptance decisions for all 10 |
| `rfdiffusion-beta-pairing` | `66D12` / `6D604` | 10 | Sequence design produced a variant for each backbone |
| `rfantibody-vhh` | `9A6DA` / `0E35A` | 1 | Saved design through SAbR, AbMPNN, Boltz-2 and Sequence QC |
| `rfantibody-mimosa-vhh` | `50698` / `62129` | 1 | Saved design through SAbR, AbMPNN, Boltz-2 and Sequence QC |
| `boltzgen-vhh` | `7D7A2` / `90114` | 1 | Recorded generator-only chain through AbMPNN and Boltz-2 |
| `boltzgen-vhh-native` | `9809A` / `729E5` | 1 | Saved design through SAbR, AbMPNN, Boltz-2 and Sequence QC |
| `boltzgen-mimosa-vhh` | `C55C2` / `EF804` | 1 | Saved design through SAbR, AbMPNN, Boltz-2 and Sequence QC |

The 12 miniprotein chains contain 300 generator records with matched downstream identities. The VHH checks trace one selected record per mode. A separate RFD4 external check found all 18 generator records in the three Boltz-2 branches of `0D091`. Where the archive used Mol-1 or a different sequence-design option, the evidence applies to the recorded composition. For five branches in `66D12`, sequence design completed and the subsequent Mol-1 step failed. The later ColabFold assessment `F92AC` has separate lineage and coverage evidence in the playbook.

## Standalone native evidence

| Route ID | Context | Worker | Retained result |
| --- | --- | --- | --- |
| `rfpeptides` | Structure | `2BDC3` | One cyclic backbone record |
| `complexa` | Structure | `C98CA` | One design record and one internally refolded complex |
| `boltzgen-miniprotein` | Structure | `B0857` | One design record and one internally refolded complex |
| `peptide-trap-idr` | Disordered sequence | `88B99` | One generator-only design record |
| `kejia-idr` | Disordered sequence | `BBFC6` | Twelve generator-only backbone records |

The standalone evidence confirms those native modes produced structures. Their handoff from the Capabilities pipeline remains open. Across the 82 combinations, 65 still lack an independently verified pipeline handoff in this review, including the five above. Historical handoffs document their original configurations. Fresh release acceptance and optional workflow variants retain their own evidence requirements.

On Studio `267a0bc`, 726 Capabilities, generator input, generator capability and pipeline dataflow regression tests passed.
