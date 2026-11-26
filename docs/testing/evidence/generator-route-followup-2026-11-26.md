# Generator coverage follow-up

The 26 November follow-up examined all 65 combinations that remained after the [first archive review](generator-route-archive-2026-11-26.md). It checked accessible logical Runs through MCP and scanned 1,025 archived design-job receipts across PC002975, PC004032 and PC004033. A route match requires the recorded engine mode, output format and target context. In particular, an archived sequence-only run does not establish the structured-target route.

## Verified native handoffs

Two additional disordered-sequence routes have native outputs linked to completed downstream assessments through the Database.

| Route | Native receipts | Assessment | Exact links checked |
| --- | --- | --- | ---: |
| `rfd4-idr` | `71432` | `B84A5` | 100 |
| `bindcraft2-idr-miniprotein` | `4C669` and `7EA80` | `72CF9` | 25 |

For each row, the read-only check matched the native run and candidate ID, the SHA-256 of the native structure bytes, the saved Database record and a completed BinderScope observation. Where the assessment used a registered derivative, it matched the derivative digest and followed the saved source-model link to the original. These are historical Database handoffs. The same-run Capabilities graphs for these two routes still need live verification.

Together with the earlier 17 pipeline chains, this establishes historical native downstream evidence for 19 of the 82 combinations. The other 63 still lack such evidence. Five of those have standalone native structures. The peer BindCraft2 R3 archive declares sequence input, and the archived ESMFold2 VHH job also declares sequence-only input. The structured-target routes therefore remain open.

## Results regression

The fixed-fixture replay exercised the real graph controller and Results writer for every advertised combination. It checked that candidate lineage, sequence and exact structure bytes survived the handoff. It also found two software failures.

- A generator-to-Results graph could fail during retained replay because the target identity was only computed when a collection-processing node was present. Retained replay now derives the identity when the upstream output lacks it, including cached outputs.
- Fresh-run input validation could block retained RFdiffusion β-pairing, Kejia or Peptide Trap results. The graph requests native generation settings when executing a fresh worker. The preflight and replay still validate the retained worker identity, project and target.

All 82 fixed-fixture handoffs passed. This suite covers the generator-to-Results boundary before route-specific downstream scientific stages.

## Evidence for the 65 combinations reviewed

| Route | Target context | Native evidence |
| --- | --- | --- |
| `rfd4-generator` | Structure | No matching receipt located |
| `rfd4-generator` | Multidomain | No matching receipt located |
| `rfd4-external` | Multidomain | No matching receipt located |
| `rfd4-internal` | Multidomain | No matching receipt located |
| `bindcraft2-generator` | Multidomain | No matching receipt located |
| `bindcraft2` | Structure | No matching receipt located |
| `bindcraft2` | Multidomain | No matching receipt located |
| `rfd3` | Multidomain | No matching receipt located |
| `bindcraft2-peptide` | Structure | No matching receipt located |
| `bindcraft2-peptide` | Multidomain | No matching receipt located |
| `bindcraft2-cyclic` | Structure | No matching receipt located |
| `bindcraft2-cyclic` | Multidomain | No matching receipt located |
| `bindcraft2-vhh` | Structure | No matching receipt located |
| `bindcraft2-vhh` | Multidomain | No matching receipt located |
| `boltzgen-vhh` | Multidomain | No matching receipt located |
| `esmfold2-vhh` | Structure | No matching receipt located |
| `rfd4-idr` | Disordered sequence | Completed archived Database handoff, `71432` → `B84A5` |
| `genie3` | Multidomain | No matching receipt located |
| `pxdesign` | Multidomain | No matching receipt located |
| `complexa-generator` | Multidomain | No matching receipt located |
| `complexa` | Structure | Standalone output `C98CA` |
| `complexa` | Multidomain | No matching receipt located |
| `bindcraft2-large-protein` | Structure | No matching receipt located |
| `bindcraft2-large-protein` | Multidomain | No matching receipt located |
| `bindcraft2-vhh-native` | Structure | No matching receipt located |
| `bindcraft2-vhh-native` | Multidomain | No matching receipt located |
| `bindcraft2-scfv` | Structure | No matching receipt located |
| `bindcraft2-scfv` | Multidomain | No matching receipt located |
| `bindcraft2-fab` | Structure | No matching receipt located |
| `bindcraft2-fab` | Multidomain | No matching receipt located |
| `bindcraft2-arp` | Structure | No matching receipt located |
| `bindcraft2-arp` | Multidomain | No matching receipt located |
| `bindcraft2-multidomain` | Structure | No matching receipt located |
| `bindcraft2-multidomain` | Multidomain | No matching receipt located |
| `bindcraft2-homo-oligomer` | Structure | No matching receipt located |
| `bindcraft2-homo-oligomer` | Multidomain | No matching receipt located |
| `boltzgen-generator` | Multidomain | No matching receipt located |
| `boltzgen-miniprotein` | Structure | Standalone output `B0857` |
| `boltzgen-miniprotein` | Multidomain | No matching receipt located |
| `boltzgen-peptide` | Structure | No matching receipt located |
| `boltzgen-peptide` | Multidomain | No matching receipt located |
| `boltzgen-cyclotide` | Structure | No matching receipt located |
| `boltzgen-cyclotide` | Multidomain | No matching receipt located |
| `boltzgen-vhh-native` | Multidomain | No matching receipt located |
| `esmfold2-scfv` | Structure | No matching receipt located |
| `rfd4-vhh` | Structure | No matching receipt located |
| `rfpeptides` | Structure | Standalone output `2BDC3` |
| `rfd3-cyclic` | Structure | No matching receipt located |
| `rfd4-cyclic` | Structure | No matching receipt located |
| `rfantibody-scfv` | Structure | No matching receipt located |
| `bindcraft2-idr-miniprotein` | Disordered sequence | Completed archived Database handoff, `4C669` and `7EA80` → `72CF9` |
| `bindcraft2-idr-large-protein` | Disordered sequence | No matching receipt located |
| `bindcraft2-idr-peptide` | Disordered sequence | No matching receipt located |
| `bindcraft2-idr-cyclic-peptide` | Disordered sequence | No matching receipt located |
| `bindcraft2-idr-vhh` | Disordered sequence | No matching receipt located |
| `bindcraft2-idr-scfv` | Disordered sequence | No matching receipt located |
| `bindcraft2-idr-fab` | Disordered sequence | No matching receipt located |
| `bindcraft2-idr-arp` | Disordered sequence | No matching receipt located |
| `bindcraft2-idr-multidomain` | Disordered sequence | No matching receipt located |
| `bindcraft2-idr-homo-oligomer` | Disordered sequence | No matching receipt located |
| `bindcraft-classic` | Structure | No matching receipt located |
| `peptide-trap-idr` | Disordered sequence | Standalone output `88B99` |
| `kejia-idr` | Disordered sequence | Standalone output `BBFC6` |
| `rfd4-antibody` | Structure | No matching receipt located |
| `rfd4-tcr` | Structure | No matching receipt located |
