# Stage 1 Pilot Report

## Recommendation

**GO for a separately authorized Stage 2 production array only.** Do not submit it from Stage 1.
The pilot selects FOV `10x16 R_j` and resolution `320x512`. The `6x18 R_j` FOV fails the formal convergence gate; `10x16 R_j` passes against `12x18 R_j` for all times and frequencies.

## Resources

The selected 320x512 image pilot averages 19.94 minutes/epoch. The largest measured RSS across all pilots is 12.33 GiB. Recommend `--cpus-per-task=32`, `--mem=24G`, `--time=01:00:00`, and array concurrency `%2`. At two concurrent tasks, 351 selected-resolution epochs project to 58.3 hours (70.0 hours with 20% margin), based on image pilots; spectrum-only production is expected to be no slower.

## Scope and Evidence

All 23 unique Stage 1 pilot cases have valid `COMPLETE.json`, finite four-frequency spectra, physical Stokes checks, and per-case resource records. The mask audit confirms that whole-domain/kink-only and tracer-off/tracer-on controls produce distinct integrated values, while the strict `40 < z/R_j < 60` mask remains independent of tracer selection. Job IDs and final queue evidence are in `pilot/job_registry.json` and `pilot/analysis/latest_scheduler_status.json`; `sacct` is unavailable on this cluster. It is a Hu et al. (2025) Figure 5 methodology reproduction on the available mapped `320x320x450` run, `t=300..1000`, with adopted camera/azimuth/temperature/distance values documented separately from paper-published values.

Stage 2 remains unauthorized: no 351-epoch production array and no 50-frequency campaign has been submitted.
