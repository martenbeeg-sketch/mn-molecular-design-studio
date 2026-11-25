# Studio release acceptance, 25 November 2026

The live campaign ran from 21:17 to 22:31 Berlin time in **Studio Release Acceptance**. It used synthetic protein inputs, the bundled imatinib example, and isolated library and database fixtures. Live runs were submitted through Studio MCP across PC002975, PC004032 and PC004033, using the Studio `c8d9506` deployment with the fixes recorded below. All scheduled acceptance runs reached a terminal state.

| Case | Initial run | Verified result or replay |
| --- | --- | --- |
| Synthetic sequence smoke | `10124` | ESMFold and Sequence QC completed. MCP returned the 20-residue result. |
| Vina | `F3BC6` | Completed. The saved imatinib pose was available as a typed pipeline input. |
| Boltz-2 affinity | `17ADB` | Completed. MCP returned one result with affinity fields. |
| Linked pose checks | `38456` | Replay `67F17` returned one result each from PLIP, PandaMap and PoseBusters on the same saved pose. |
| Native ligand interaction classifier | `C49E3` | Completed on the same saved pose with one measured Results row and 16 contact residues. |
| Batched monomer prediction | `58835` | Replay `D1FD9` completed three sequence predictions in one native launch, with three rows in Results. |
| OpenMM and saved trajectory reuse | `9CE14` | Replay `3256A` completed with one production frame. Reuse run `E752E` and analysis job `2D047` reread the exact saved DCD and bonded topology and selected one frame. |
| GROMACS | `453A0` | Completed equilibration, production and frame analysis in children `F5131` and `15EA5`. The result exposed a separate GROMACS sampling bug. Replaying the production stage from the saved equilibration with the fix produced two trajectory frames at 0 and 5 ps. |
| Mixed-workstation fan-in | `53C7B` | Replay `FB731` collected seven Results rows. |
| Mixed-workstation cancellation | `49C54` | Replay `921A0` stayed cancelled and its active remote job was cancelled. The other worker job had yet to be submitted. |
| Candidate filter | `EF9A0` | Replay `ADC00` selected one row at the lower cutoff and zero at the higher cutoff. The rejected candidate remained visible for review. |
| Saved MSA and connected template | `0D154` | Replay `10FE7` completed with Mol-1. Chain A came from the recorded reference alignment cache, and the connected template was staged through its portable reference. |
| Unconditional fold to sequence design and refold | `FCA14`, `51EC3` | Replay `B88E0` delivered one ESMFold structure after Protpardelle, ProteinMPNN and Sequence QC. Exact parent identity resolved a shared native alias, and the refolder accepted the QC record with the verified original sequence. |

## Prediction engine coverage

All 15 configured synthetic monomer cases completed and each produced one Results row with a 20-residue Sequence QC measurement.

| Engine or mode | Verified run |
| --- | --- |
| Boltz-2 | `0E4B8` |
| ESMFold2 Fast | `F2058` |
| Mol-1 | `C4F84` |
| AF3 PyTorch with OpenFold3 P2 weights | `4338A` |
| Chai-1 | `45DF0` |
| Protenix V2 | `05B5D` |
| OpenFold3 P2 JAX | `1212B` |
| LocalFold with OpenFold3 weights | `4FAB6` |
| ESMFold2 Full | `EA6FC` |
| ColabFold | `C7D6A` |
| OpenFold3 OB0 | `A7B9D` |
| Protenix V1 | `BC7A7` |
| AtlasFold | `FC267` |
| OpenDDE native | `F0AAD` on PC004032, `D381A` on PC004033, `9A8C4` on PC002975 |
| RF3 | `E6E75` |

## Deployment and remaining coverage

The test exposed worker-slot starvation when pipeline controllers were counted as execution workers. Studio now starts a separate controller pool through the same scheduler. Test jobs interrupted by the shared-pool stall were requeued under the new routing. The GROMACS test also exposed a missing Docker CPU quota on a cgroup v2 host. Its running container was given the four CPUs assigned by Studio, and the quota detection now recognises the cgroup v2 root controller list.

The short MD runs now have a reporting interval capped at the production length. The GROMACS runtime writes that interval into the production MDP and uses the same interval for analysis. The original GROMACS result remains available with its original single frame. Mixed-workstation cancellation resolves scheduler job IDs directly. Saved-pose handoffs preserve the job identifier and the selected model provenance. The live Runs API lists the trajectory analysis replay independently and reports the successful GROMACS pipeline as completed. The cancelled GROMACS attempt remains in its execution history.

ColabFold had failed on PC004033 because the multimer parameter cache was missing. The five weight assets and completion marker were copied from the existing PC002975 cache and verified on both peers. OpenDDE failed on PC004032 due to a missing NVML Python binding. PC004033 needed the RTX 5090 profile. Both passed with the same verified image and pinned NVML overlay. PC002975 passed with its existing RTX 4090 image. LocalFold's official AF3 checkpoint requires its separate access terms. The replay used OpenFold3 weights. [NVML package version](https://pypi.org/project/nvidia-ml-py/12.570.86/)

The regular-user Streamlit check reached sign-in and remains unverified. The library, approval and database-admission checks used isolated fixtures and contract tests. Live passing-design publication and the linked BinderScope workflow were outside the fixture coverage. A focused MCP, Results, database-policy and project-lifecycle suite passed 161 tests. The library preparation and active-learning suite passed 123 tests. The final GROMACS and logical-run regression group passed 97 tests.

This run has no per-cell Capabilities coverage record. The 15 prediction cases cover monomer prediction and Sequence QC. Per-engine binder-format and workflow-variant executions, regular-user pipeline creation starting from the table, and the resulting combinations of optional steps remain open in the [acceptance playbook](../pipeline-release-acceptance.md#capabilities-coverage-and-student-entry-point).

The runtime fixes are deployed to the three workstations. The code changes are in the isolated `codex/bbf14-pipeline-handoffs` worktree and await a commit.
