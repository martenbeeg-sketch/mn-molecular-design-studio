# Generator engines handoff — 2026-10-17

## User's latest constraint

Run pytest; **do not run smoke tests**. The user will do those later with another
model. Do not launch jobs automatically when resuming this handoff.

## Completed

PXDesign, RFdiffusion, Genie3 and Protpardelle-1C have generation-only tiles,
shared green input pages/node forms, 109 bundled examples, native runners,
normalized result displays and editable upstream-inspired pipeline templates.
See [implementation notes](generator-engines.md) for tasks, limitations and
references. All native generator workflow IDs end in `_generate`.

Full pytest result: **382 passed**, one existing Starlette/httpx deprecation
warning, 14.62 seconds. Log: `/tmp/generator-full-pytest.log`.
The frontend and wheel builds passed. CPU tests never run Docker inference.

Two API-discovered bugs were fixed: named RFdiffusion checkpoints are not
interpreted as file paths; target chain identity matching excludes source binder
chains. Motif ranges are preserved when examples pass through the shared picker.

## Runtime

- Repo: `/home/user/programs/git-projects/mn-molecular-design-studio`
- Python: `/home/user/mambaforge/envs/mn-molecular-design-studio/bin/python`
- Streamlit: port 8502, existing process, unchanged.
- Studio APIs refreshed: port 8766 PID 3008140; port 8767 PID 3008212.
  Logs: `/tmp/studio-api-8766-native-generators.log`,
  `/tmp/studio-api-8767-native-generators.log`.
- Scheduler PID 3525059 was not restarted. Preserve paused V7 campaign state.
- Weights: configured shared reference folder, currently
  `/mnt/db/reference_files`. Two missing RFdiffusion task checkpoints were
  downloaded there; provenance is recorded alongside them.
- Before-change backups: `/tmp/generator-pages-before/`.
- Whole repository is untracked in Git. Preserve all prior changes.

## Later smoke tests, only when authorized

Use the existing scheduler/API and the provided task examples. First validate
generation alone for each task, then explicitly choose a small example for each
pipeline template. Do not assume mocked command tests prove GPU inference.
Readiness checks confirm expected paths exist; they do not load model weights.

Genie3's template uses Studio's AF2 Initial Guess adapter rather than the exact
upstream ColabFold evaluation wrapper. RFdiffusion omits optional FastRelax.
Protpardelle's template starts with an unconditional monomer; replace downstream
prediction steps as appropriate when switching to complex tasks. These differences
are disclosed in each template's About section and the implementation notes.
