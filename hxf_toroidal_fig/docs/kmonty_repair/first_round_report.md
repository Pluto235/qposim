# kmonty G0/G1 first-round report

Branch: `codex/kmonty-g0-g1`; base commit `1c25e1e1`. Existing worktree modifications were preserved.

Implemented in `postprocess/radiation/kmonty_readonly_audit_20260916/source/kmonty`:

- `src/compton.c`: Compton electron sampling dispatches through `sample_beta_distr`, so thermal sampling does not call the κ width solver. Rejection exhaustion records a diagnostic and exits non-zero instead of changing `Thetae`.
- `src/kappa_sampler_num.c`: removed the fixed cold-electron (`gamma=1.000001`) fallback.
- `src/electron_distribution.c`: added adaptive log-gamma moment integration, complete Maxwell–Jüttner target moments, the explicit `kappa_exp_v1` smooth exponential cutoff (`gamma_cut=1e5` as a physical scale), log-space bracketing/bisection, and thread-local cache identity including κ and cutoff.
- `tests/model_identity.json`: records the model identity and legacy reference name.

The reproducible build command was:

    make -C postprocess/radiation/kmonty_readonly_audit_20260916/source/kmonty sphere -j2

It stopped before compilation because the active MPI compiler environment lacks `gsl/gsl_integration.h`. Consequently all numerical acceptance rows are `NOT_RUN`; raw logs and manifests are under `validation/20260916_g0g1/G0` and `validation/20260916_g0g1/G1`. No production spectra or G2–G5 work was run.

Thermal dispatch, failure handling, and cold fallback removal are bug fixes. The κ PDF definition is an explicit model change (`kappa_exp_v1`). Remaining work is a formal GSL/MPI development environment, independent reference CDF/moment validation, and complete build/cache manifests before claiming PASS.
