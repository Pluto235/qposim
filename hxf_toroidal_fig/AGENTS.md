# hxf_toroidal_fig instructions

These instructions supplement `../AGENTS.md`. Paths below are relative to
`hxf_toroidal_fig/` unless stated otherwise.

## Source and case layout

- `MHD1-2D/`, `HD/`, `MHD1/`, `MHD2/`, `MHD3/`, and `MHD2w/` are the six
  simulation cases. Use the selected case's `definitions.hpp`, `init.cpp`,
  makefile, runtime input, and launch script.
- `common/gpluto.mk` supplies shared build logic; `common/rhd3d_init.cpp` and
  `common/rmhd3d_init.cpp` supply shared 3D initialization. Inspect include
  relationships before editing: a shared change can affect multiple cases.
- MHD1-2D now owns its solver, plotting and validation implementation;
  the old top-level forwarding arrangement has been removed.
- `postprocess/` contains diagnostics and regression tests; `postprocess/visualization/`
  contains visualization tools, while `viz/` retains historical outputs.
  `postprocess/radiation/kmonty/` is the maintained kmonty source path in
  each of the two qposim worktrees described below.
  `postprocess/radiation/kmonty_*/` directories contain only tests, configs,
  validation manifests, reports, and frozen evidence. The external
  `tools/RaptorP/kmonty` tree is the user-designated original/stable source
  baseline, now imported into qposim; preserve it as historical provenance,
  not as the current κ development source.
  `postprocess/radiation/` owns radiation source, configs, tests and reports;
  old run-local source paths are compatibility links. See `docs/layout.md`.
- Use the README files for scientific context. Verify grids, solver settings,
  output cadence, and restart paths in the actual inputs and launch scripts;
  historical descriptions are not authoritative runtime configuration.

Non-default case inputs live in `configs/`; helper and continuation scripts
live in `scripts/`. Historical names may be compatibility links, not duplicate
sources. Keep radiation simulation input paths explicit in `paths.json`.

## kmonty version and experiment management

- Use one maintenance repository (`/mnt/mydisk/qposim`) with two independent
  worktrees. Both versions are Git-tracked at the source path above:
  - Stable: branch `codex/kmonty-stable`, worktree
    `/mnt/mydisk/qposim-kmonty-stable`, baseline tag `kmonty-stable-20261006`.
    This preserves the user's chosen `tools/RaptorP/kmonty` version, including
    its historical modifications; "stable" is not a new physics certification.
  - Active κ development: branch `codex/kmonty-pump-fix`, worktree
    `/mnt/mydisk/qposim`, baseline tag `kmonty-kappa-checkpoint-20261006`.
    `codex/kmonty-g0-g1` remains a historical branch, not the latest source.
- Keep baseline tags fixed. Make new changes as scoped commits; do not switch
  the active development directory to stable, restore historical G3 trees,
  or automatically merge κ changes into stable. Review and validate any
  explicitly requested transfer. Never stage unrelated changes with `git add .`.
- Save experiments separately under
  `/mnt/mydisk/qposim-kmonty-experiments/stable/<experiment>/` or
  `/mnt/mydisk/qposim-kmonty-experiments/kappa/<experiment>/`.
  Keep validation reports/evidence in the relevant `kmonty_*` directory.
  Copy the selected executable and inputs into each experiment directory;
  do not link to a development binary that may be rebuilt during a run.
- Keep build directories, outputs and writable numerical caches independent.
  Read-only simulation inputs may be shared. Record source commit and dirty
  diff, build command/macros/toolchain, binary SHA256, input paths, seed and
  MPI/OMP layout for each experiment. Do not commit large generated outputs.
- In the relocated stable worktree, explicitly set
  `GSL_PREFIX=/mnt/mydisk/qposim/tools/raptorp-deps/prefix/usr`; do not rely on
  the old makefile's relative dependency path. Use the verified MPI toolchain
  at `/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/comm_libs/mpi/bin` and separate
  build directories. Follow the existing copy-before-clean rules.
- See `docs/kmonty_repair/version_management.md` and the branch-specific
  `version_baseline.json` for provenance and source hashes. These worktrees
  share local Git history; they are not independent backups. No qposim remote
  was configured at setup; verify current remotes before claiming remote backup.

## Build and execution

- Build from the selected case using its makefile and the inherited
  `PLUTO_DIR` default. MPI builds use `ARCH=linux.mpicxx.acc.defs`.
- Inspect the case's `submit.slurm` and the script it invokes. Shared 3D
  launch logic is in `common/run_3v100.sh`; MHD3's main submission invokes
  `MHD3/scripts/run_staged.sh` with staged runtime inputs.
- Existing executables may be reused by launch scripts. After changing
  sources, compile-time options, or architecture, explicitly rebuild the
  affected executable before testing or submitting it.
- Apply the parent's copy-before-clean rule in the selected case directory.
  `common/gpluto.mk` removes only build-directory `.o` files; do not
  broaden cleanup into run data or unrelated builds.
- Use a fresh run directory for smoke tests. Preserve the executable, input,
  and logs needed to identify what ran. For continuation work, verify the
  intended restart index, snapshot, and matching inputs before launching.

## Validation entry points

- `MHD1-2D/`, `HD/`, and `MHD1/` provide `configs/pluto_smoke.ini`. Copy the selected
  input into a fresh test directory as `pluto.ini` and bound the run with
  `-maxsteps 10`, using the appropriate GPU/MPI allocation.
- MHD2 and MHD2w need test inputs retaining their own runtime field parameters;
  sharing MHD1's compiled physics does not make its runtime input equivalent.
- MHD3 has no top-level `pluto_smoke.ini`. Its
  `experiments/axial_solver_smoke/README.md` documents a full-grid solver test;
  `experiments/hxf_ref_smoke/README.md` documents a separate HXF reference compatibility
  test. Check which source and input each test exercises before using it as
  validation for a change.
- Validate MHD1-2D initial output with
  `python3 MHD1-2D/validate_initial_state.py <run-directory>`.
- Run post-processing regression tests from this project directory with
  `python3 -m unittest discover -s postprocess/tests -v`.
  Follow `postprocess/README.md` for affected diagnostics and figure generation.
- For radiation changes, use the tests and instructions of the specific
  `postprocess/radiation/` pipeline. For kmonty changes, build from
  `postprocess/radiation/kmonty/` and keep test/report artifacts in the
  relevant `kmonty_*` directory. Do not substitute an unrelated solver smoke
  test for radiation validation.

## Results protection

Treat simulation snapshots and restart metadata as inputs to post-processing;
keep analysis read-only with respect to those inputs. Preserve run-local source,
configuration, reports, and existing results when changing analysis code.
Do not infer that a directory is disposable from its job ID or location under
`runs/`.
