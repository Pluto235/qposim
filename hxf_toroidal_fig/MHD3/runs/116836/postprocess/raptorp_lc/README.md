# RaptorP LC Stage 1

This worktree is an isolated Hu et al. (2025) Figure 5 methodology
reproduction on the available mapped-grid MHD3 run, `t=300..1000`.  It is not
a strict numerical reproduction: the available grid is mapped `320 x 320 x
450`, whereas the paper uses uniform `500 x 500 x 600`; the available series
contains 351 samples rather than the paper's 401.

`source/` is copied from the validated Figure 3 worktree without its Git
worktree pointer, build directory, or executable.  The provenance chain is
recorded in `provenance/` as upstream commit plus validated Figure 3 patch plus
the LC-only patch.  No Figure 3 source file is modified by this workflow.

Stage 1 creates four targets from one physics source: `RAPTOR_LC4_IMAGE`,
`RAPTOR_LC4_SPEC`, `RAPTOR_SED50_IMAGE`, and `RAPTOR_SED50_SPEC`.  The LC
frequency file is generated from `scripts/frequency_contract.py`; the optical
and X-ray values are not duplicated rounded literals.

The volume control is separate from the released fractional tracer control.
`RAPTORP_VOLUME_MASK=kink-only` masks the original physical coordinate with
the strict open interval `40 < z/R_j < 60`, before fluid coefficients are
evaluated.  `whole` is retained only for the specified pilot control.  The
production runner rejects any other mode, and this Stage 1 worktree contains no
production submit script.

Run local gates with:

```bash
python3 tests/test_stage1_contracts.py
python3 scripts/preflight.py
```

Pilot arrays are submitted only through `scripts/submit_pilot.py`, which caps
concurrency and appends the returned Slurm ID to `pilot/job_registry.json`.
