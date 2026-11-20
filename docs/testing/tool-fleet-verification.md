# Studio fleet verification · 20 November 2026

The initial verified build was `b41a988`. The canvas release is `d8d0479d0ff6`. The three production workstations are pc002975, PC004032 and PC004033. The regression work used separate Git worktrees.

## Production workflow retests

The cases ran through the paired Studio MCP connection under **Tool Validation**. Each verified row has a successful native Results record. Manual collection kept the calculation results outside the reusable molecule library.

14/14 successful named production retests
| Case | Workstation | Native run | Status |
| --- | --- | --- | --- |
| BindCraft2 | PC004032 | `20261120-163141-724e2f50` | Results verified |
| BindCraft2 Generator | PC004033 | `20261120-163144-7c4383cd` | Results verified |
| Interface Clustering | pc002975 | `20261120-163147-3142c0f0` | Results verified |
| AutoDock-GPU | pc002975 | `20261120-162745-339a5866` | Results verified |
| Surf2Spot | pc002975 | `20261120-162838-82106095` | Results verified |
| Protein Surface Analysis | pc002975 | `20261120-162840-3611986e` | Results verified |
| Mol-1 | pc002975 | `20261120-162837-be3699e7` | Results verified |
| Complexa Generator | pc002975 | `20261120-162838-d4f1c290` | Results verified |
| Complexa | pc002975 | `20261120-162836-c98cae97` | Results verified |
| OpenMM | pc002975 | `20261120-162839-4cae4de1` | Results verified |
| GROMACS | pc002975 | `20261120-162840-823b7536` | Results verified |
| BoltzGen Fold | pc002975 | `20261120-165704-ea179877` | Results verified |
| OMTRA | PC004033 | `20261120-163546-28d729cf` | Results verified |
| MD Analysis | pc002975 | `20261120-163835-d3769a35` | Results verified |

The additional BoltzGen Fold transfer to PC004032 returned both models on the final release. The replay used portable `reference:///` paths after a manually reconstructed cohort with bare relative paths failed file resolution. The successful run was `20261120-170241-5f73e221`.

## Canvas and Results

- All 95 tiles and their modes map to registered canvas operations. The live catalogue contains 129 operations and 4,866 typed port matches. A port match identifies a potential connection. File requirements, permissions and scientific constraints still apply. [Connection matrix](tool-connection-matrix.json)
- Prepare library → Create centroids → Results completed with three input compounds and three centroids. The run was `20261120-163017-b6bdb653`.
- LDDM → Vina → Results completed. Vina scored all four generated compounds from the saved compound set. The run was `20261120-164037-b1580e1c`.
- Complex Preparation → OpenMM → MD Analysis → Results completed. The analysis used the exact recorded topology and trajectory, with file hashes checked at the handoff. The run was `20261120-163310-d804d568`.
- The released production canvas also completed Saved MD run → MD Analysis → Results through paired MCP. It linked one analysis in `20261120-170012-33ea61a1` and MCP returned all 95 tile mappings.
- Headless Streamlit opened 15 Results views in isolated workspaces using copied native test artifacts.

## Database admission and deployment

- BindCraft2 executed fresh test runs on both compute hosts. Its native zero-pass control registered zero design observations under automatic collection. The separate positive control used one genuinely accepted historical output. It registered one design with a byte-identical structure, and repeating collection added no duplicate.
- The OMTRA protein-conditioned checkpoint and container content were verified on all three hosts. OMTRA executed on the coordinator and on PC004033. The pocket-only test generated disconnected molecules. The Vina preparation stage rejected these as expected and preserved the generator outputs for inspection.
- The coordinator had a separate scheduler bridge that survived application service restarts. Refreshing the bridge on the compute hosts resolved the initial BindCraft2 failures. The final activation verified the bridge's active source directory on every workstation. [Bridge refresh](../../scripts/refresh_studio_scheduler_bridge.py)
- AlphaFold3 accepted the existing versioned terms acknowledgment for the paired account. Changing the upstream terms version requires a new acknowledgment.
- The rollout compared account and membership tables and scientific row counts at activation. All SQLite recovery copies passed integrity checks. The deployment preserved accounts, memberships and existing research records.
- The three live SQL copies passed quick integrity checks. Checking 21 recorded validation run IDs found zero result observations on every host.

## Regression scope

The canvas release passed 135 focused regression checks. The full regression run on `b41a988` reported 4027 passed, 64 skipped, 103 warnings in 2058.05s (0:34:18). The selected native cases verify input handling, execution, artifact integrity and Results presentation. Coverage is limited to the inputs and handoffs listed here. The fresh BindCraft2 run had zero passing candidates. Its positive admission evidence comes from the historical replay.
