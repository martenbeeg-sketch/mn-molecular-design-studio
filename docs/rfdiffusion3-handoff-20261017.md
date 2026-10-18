# RFdiffusion3 implementation handoff — 2026-10-17

## Status

Implementation and requested verification are complete. The RFdiffusion3 tool is available in Molecular Design Studio using the existing `mn-foundry:cu128` image. The input flow follows the Tamarind task structure while retaining Studio projects, scheduler, results, and sequence-design integration. The interface uses the Studio's green styling.

See [the smoke report](rfdiffusion3-smoke-report-20261017.md) for all six completed task runs, run IDs, timings, output checks, and limitations. The full Python test suite passed: **253 passed**, with one existing Starlette/httpx deprecation warning.

## Important validation limitation

All six task smokes completed successfully through the scheduler and produced one normalized backbone candidate each. They were submitted before the UI/workflow default was changed to checkpoint-matched `low_memory_mode=False`; those smoke jobs used experimental low-memory mode and the pinned image logged missing chunked pairwise embedding keys that it reinitialized. The smoke runs establish task execution, output normalization, and result rendering, not comparative design quality. No additional predictions were launched after the default change. New jobs default to `low_memory_mode=False`; the optional low-memory checkbox is marked experimental.

The outputs are generated backbones. They have not passed sequence design, refolding, or designability filters. The paused V7 campaign was left untouched.

## Runtime state

- Repository: `/home/user/programs/git-projects/mn-molecular-design-studio`
- Python: `/home/user/mambaforge/envs/mn-molecular-design-studio/bin/python`
- Results: `/mnt/data/RESULTS/mn-molecular-design-studio/workdir/runs`
- References: `/mnt/db/reference_files`
- API: `http://127.0.0.1:8767` (PID 2648663 during final verification)
- UI: `http://127.0.0.1:8502`
- Shared scheduler: PID 3525059; it was not restarted
- Six smoke jobs remain in project `RFdiffusion3 task smoke · 2026-10-17` (ID `bf69d1b3c2bf4ce69e8052feb3a03f38`).

The repository contents are untracked in this checkout; preserve them and do not reset the working tree.
