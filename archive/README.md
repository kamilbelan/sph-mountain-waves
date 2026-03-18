# archive/

This directory contains historical code kept for *reference* only. **None of it is called by any active code path and never should**

## Contents

- `legacy_formulations/` — each script contains the domain setup, physics, and I/O. These were superseded by a modular approach in `src/formulations/` + `src/core/`.
  - `obsolete/` — even earlier standalone scripts, before the legacy modules.
- `legacy_scripts/` — the old per-formulation run scripts from, superseded by the unified `scripts/run_sim.jl`.

## Do not add active code here

All new code should stay in `src/formulations/` and `src/core/`. This directory is intentionally outside `src/` so that `src/` contains only active code.
