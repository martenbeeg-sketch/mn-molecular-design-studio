# Native pipeline acceptance, 26 November 2026

Studio worker `pc002975`, RTX 4090. The coordinator and workers came from the acceptance checkout based on `c8d9506`. A regular fixture account used the authenticated Studio API to validate, save, submit and inspect each pipeline in the isolated acceptance project. Input structures came from a saved 20-residue synthetic prediction and a generated ethanol SDF. The checks made no research database writes.

| Case | Pipeline / worker | Observed result |
| --- | --- | --- |
| Saved prediction → TMol → Sequence Composition | `76D25` / `19265` | TMol Score, Cartesian minimization and FastRelax produced structures and finite energy values. Sequence Composition failed because the TMol output dropped the recorded protein-chain sequences. |
| Same inputs after the fix | `C5C2F` / `049BA`, `4652F`, `F94C2`, `67B5E` | Three TMol modes completed and Sequence Composition produced one valid 20-residue result. The output retains the parent identity and the original sequence record. |
| BioEmu-2 equilibrium sampling and structure prediction | `F458A` / `BE426`, `8ABCD` | Each mode published an XTC with a matching mmCIF topology and two frames that passed the configured physicality checks. The sequence modes reused an exact cached A3M. |
| BioEmu-2 propagation from an explicit mmCIF input | `89D2B` / `9CDDA` | Propagation used the saved prediction's coordinate file and published an XTC and matching topology with two frames that passed the physicality checks. |
| Ligand SDF → TMol Score | `2083D` / `AB4D6` | The native TMol worker parameterized the ethanol input and published one scored structure with a finite energy. |

Streamlit component tests with the two real TMol run snapshots verified the failed-step warning, the recorded error and navigation to the Sequence Composition result. The replay showed a completed step and the corresponding result link. The Results API grouped the terminal pipelines with their workers, and the BioEmu trajectory preview resolved every XTC/topology pair.

While the standalone ligand graph was running, Results briefly listed its TMol worker under a separate logical Run ID. That reference returned 404 after the worker was grouped into its actual pipeline. The submission code now writes the pipeline and step identities when creating the worker job, before dispatch and Results indexing. On replay `73134`, 71 polls from submission through completion found one stable Run reference and worker `11E9A` under it.

The separate BioEmu Propagate mode still accepts an explicit structure file or a completed design collection. Passing a named Run output directly to Propagate remains outside the current saved-run collector contract. BioEmu's outputs are trajectories, while that collector merges normalized candidate records. The acceptance run used the explicit structure input.

The focused Capabilities and QC regression suite returned 506 passes with the native QC tests disabled. All 11 CPU QC container tests passed separately, including PLIP, PandaMap and PoseBusters. The QC test assertions cover all three tools' optional source-provenance arguments. The pipeline graph and Results regression suite passed all 120 tests after the submission fix.

The two tested runtime files were synchronized and checksum-verified on PC002975, PC004032 and PC004033. The live checks ran on PC002975. The source changes remain in the acceptance checkout.
