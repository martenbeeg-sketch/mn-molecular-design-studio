# RFdiffusion3 task smoke report — 2026-10-17

All task types were submitted through the Studio API and shared scheduler using the existing `mn-foundry:cu128` image. Every job exited successfully, produced one native result, and normalized to one `generation.backbone` candidate.

Project: `RFdiffusion3 task smoke · 2026-10-17` (`bf69d1b3c2bf4ce69e8052feb3a03f38`).

| Task | Run ID | Elapsed | Output check |
|---|---|---:|---|
| Protein binder (`50E23`) | `20261017-000358-50e2305b` | 22m 10s | Binder chain A, 100 aa; target chain B, 115 aa. Target alignment maps 115 CA atoms, RMSD 0.048 Å, status `ok`. |
| DNA binder (`96B9E`) | `20261017-000358-96b9e13c` | 13m 22s | DNA chains A and B, 12 residues each; protein binder chain C, 100 aa. Result viewer loaded the generated structure. |
| Ligand binder (`19871`) | `20261017-000359-1987191a` | 1m 14s | Protein binder chain A, 100 aa; ligand IAI retained in chain B. Reported minimum ligand distance 4.06 Å and zero ligand clashes. |
| Enzyme, atom-conditioned (`06035`) | `20261017-000359-06035874` | 2m 24s | Enzyme backbone chain A, 100 aa; NAI and ACT retained. Source motif residues A108, A139, A152, A156 map to output A11, A18, A90, A44. |
| Custom JSON (`C16F3`) | `20261017-000359-c16f32e2` | 38s | One 64-aa protein backbone candidate. |
| Enzyme, residue-conditioned (`E7FCC`) | `20261017-000359-e7fcc679` | 2m 55s | Enzyme backbone chain A, 104 aa, with nonprotein ligands retained. Motif residues A108, A139, A152, A156 map to output A21, A42, A63, A84. |

## Checkpoint-mode caveat

These smokes were submitted with `low_memory_mode=True`. In the pinned image, the runner logged missing chunked pairwise embedding weights and reinitialized those modules. Runs completed, but this is not evidence that low-memory outputs match the intended checkpoint quality. The workflow and UI now default to `low_memory_mode=False`, matching the native engine default. Low-memory mode remains available as an explicitly experimental option. No extra GPU runs were started after making that change.

One optional secondary-structure metric noted that `pydssp` was unavailable; geometry checks and candidate normalization succeeded. All results are diffusion backbones only: they still require sequence design and cofolding/refolding, and no filter-passing or designability claim was made.

## Verification

- Full pytest suite: **253 passed**, one existing Starlette/httpx deprecation warning.
- Manual UI check: RFdiffusion3 task form, green styling, and result viewer were inspected in the running local app.
- The target-chain ID remapping in the protein candidate was fixed in normalization using the original input target-chain metadata, then refreshed from the saved native output (no prediction rerun).
- The paused V7 campaign was not modified.
