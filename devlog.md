# qposim devlog

Lightweight running log of Claude-mediated modifications.
Newest at top. Format: `- **<version>** · YYYY-MM-DD HH:MM · [<sha>] <summary>`.
Detailed design rationale -> `docs/`. Authoritative history -> `git log`.

---

- **qposim** · 2026-10-06 21:30 · [uncommitted] Version kmonty stable and kappa sources in separate qposim worktrees; preserve source bytes, freeze baseline tags, and separate experiment directories.

- **qposim** · 2026-09-11  · [uncommitted] Add initial energy-matched κ distribution and finite-torus dust seed contract with source/scatter diagnostics; baseline and κ/dust builds plus standalone tests pass.

- **qposim** · 2026-09-10 17:36 · [uncommitted] Launch isolated PKS exploratory fast-light production after hash-matched seed fixes, 24 transport/interface assertions and three successful pilots; add momentum-selected 0.8–1.6 degree clustered scattering diagnostics and fail-closed checks; submit three Ns=1e7 seeds as Slurm 118308, retain unchanged MHD scaling and flag unconverged pilot GeV weights.
- **qposim** · 2026-09-10 17:20 · [uncommitted] Repair PLUTO kmonty seed support/interpolation, rejection checks, linear roulette reduction and count quadrature; trace the first NaN against fixed upstream, pass 24 prior checks and final 12 seed/6 transport runs, validate zero emission, floor invariance and 1/4-thread statistics; preserve 44 prior artifacts and archive provenance without PKS production.
- **qposim** · 2026-09-10 14:55 · [uncommitted] Repair initial kmonty boundary/energy-guard/shared-state/logging defects and add isolated PKS 1510-089 preflight; pass 20 original assertions and four transport checks, but reproduce 33,709 invalid weights in 167,929 generated seeds and 33,727 invalid escape records; stop production and SED plotting at the failed validation gate, preserve source/build/config/log provenance.
- **qposim** · 2026-09-10 10:35 · [uncommitted] Add isolated PLUTO→kmonty audit with six Slurm seed/sample runs, scatter-order and clustered-weight diagnostics, tracer/origin attribution and snapshot checks; pass 18 analytic assertions and reproduce two baseline defects (z-face lookup and high-energy photon rejection), identify logged luminosity 4π error and unresolved high-energy convergence; preserve original transport and results.
- **qposim** · 2026-09-07 21:16 · [uncommitted] Complete PKS 1510-089 MJD 58200–58850 LAT likelihood SED (7 detections and 1 upper limit), calibrated UVOT photometry (53 accepted exposures), and audited historical Herschel SPIRE reference; assemble 40 traceable SED rows and fit 22 event detections with agnpy/CXC Sherpa from three starts. Export component/residual PNG/PDF, provenance and scripts; flag p1 boundary, super-Eddington thin-disk assumption and initialization degeneracy rather than finalize physical parameters. Pass 31 stage-two tests and numerical model checks; preserve prior XRT and simulation products.
- **qposim** · 2026-09-07 17:13 · [uncommitted] Convert 3038 PKS 1510-089 event XRT counts into nine response-folded spectra and five SED constraints; fit shared Γ=1.296 with independent epoch normalizations, grouping/NH/curvature checks and 200 native-channel count simulations. Add separately labeled Chandra 2001/2010 pileup-model references and Suzaku XIS0 2006/2009 spectra with regenerated attitude-aware responses and eight historical band constraints; export traceable CSV/JSON/PNG/PDF and scripts, preserve earlier products, and pass 21 stage-two plus 10 existing regression tests.
- **qposim** · 2026-09-07 16:35 · [uncommitted] Add explicit slow-light populated-cell audit and rebuild all six maintained RaptorP targets while the isolated query-event fast/slow production rerun continues.
- **qposim** · 2026-09-07 16:03 · [uncommitted] Build PKS 1510-089 event-anchored SED preview from ALMA/ASAS-SN/NEOWISE with separate AllWISE/LAT historical references and explicitly non-fitting Swift model guide; preserve raw provenance, export CSV/PNG/PDF, pass 13 stage-two and 10 existing tests; event LAT/XRT/UVOT spectral extraction remains pending.
- **qposim** · 2026-09-07 15:02 · [uncommitted] Add PKS 1510-089 stage-two event definition and verified multiwavelength coverage inventory with cached archive provenance, Swift anchor candidates, reproducible plots, eight new tests and ten existing tests passing; defer SED extraction and fitting.
- **qposim** · 2026-09-07 14:13 · [uncommitted] Add isolated job 116836 query-event primitive reconstruction and local temperature conversion; pass analytic tests and all nine pilot cases, then launch job 117621 for 339 fast plus 339 slow 521 nm cases with transactional audits and preserved historical products (production in progress).
- **qposim** · 2026-09-06 01:13 · [uncommitted] Add PKS 1510-089 cadence/index/brightest-cycle robustness analysis with official LCR data, eight gap-matched AR(1) tests, WWZ/GLS and phase comparison plots, and ten passing regression tests.
- **qposim** · 2026-09-05 23:35 · [uncommitted] Extend cycle-band and phase-folding plots to all ten Fermi-LAT event windows with matched cycle colors, phase tables, and a complete figure gallery.
- **qposim** · 2026-09-05 17:10 · [uncommitted] Reorganize hxf cases and radiation sources with compatible legacy paths, preserve run data, remove 10 unused backups, and validate builds, tests and figure equivalence.
- **qposim** · 2026-09-05 16:51 · [uncommitted] Add an A1 cycle-band and phase-folding demonstration with cycle colors, expected maxima, and weighted phase-bin means.
- **qposim** · 2026-09-05 16:48 · [uncommitted] Overlay uncertainty-weighted sinusoid fits at each Fermi-LAT window's highest WWZ power period, label boundary maxima, and save fit parameters.
- **qposim** · 2026-09-05 16:41 · [uncommitted] Reuse Foster WWZ for ten Fermi-LAT QPO windows and add Morlet CWT panels, 1000 gap-matched AR(1) simulations per window, reproducible spectra, and significance summaries.
- **qposim** · 2026-09-04 18:40 · [uncommitted] Add reproducible Fermi-LAT weekly light-curve downloads, MJD event tables, provenance, and figures for five A-level short-QPO sources.
- **qposim** · 2026-09-02 12:41 · [uncommitted] Add a shared-scale two-row fast/slow 521 nm WWZ map with dashed highest-global-peak period annotations.
- **qposim** · 2026-09-02 12:21 · [uncommitted] Add matched 521 nm fast/slow Foster WWZ analysis, peak annotations, shared-scale time-frequency maps, and periodicity comparison report.
- **qposim** · 2026-09-02 12:12 · [uncommitted] Add labeled gray dashed fast- and slow-light mean flux-density lines to the full 521 nm light-curve figure.
- **qposim** · 2026-08-25 17:03 · [uncommitted] Authorize and launch the 339-case 521 nm slow-light production array with persistent retry, audit, and final-report controller.
- **qposim** · 2026-08-25 16:33 · [uncommitted] Add and validate isolated 521 nm fast/slow SPEC light-curve targets, six-case pilot, and gated 339-case slow production workflow for job 116836.
- **qposim** · 2026-08-25 14:10 · [uncommitted] Add and validate an isolated paired fast/slow-light RaptorP workflow for job 116836, including a six-case scientific pilot and comparison report.
- **qposim** · 2026-08-25 12:16 · [uncommitted] Remove failed t=300 HLLD branch and PPM/PLM smoke artifacts, reclaiming 24.8 GiB while retaining the TVDLF baseline and axial memory-validation case.
- **qposim** · 2026-08-25 12:14 · [uncommitted] Conclude t=300 HLLD tests: CFL 0.1 becomes catastrophic by t=308.85, while CFL 0.05 with PPM or PLM still triggers energy-solve and negative-pressure fallback near t=300.3; retain TVDLF baseline.
- **qposim** · 2026-08-25 11:13 · [uncommitted] Add valid HLLD+PLM restart smoke: it completes 300 steps but fallback starts at t=300.314, so LINEAR reconstruction does not resolve the instability.
- **qposim** · 2026-08-25 10:54 · [uncommitted] Record HLLD CFL 0.05 restart smoke: it completes 300 steps but recurrent fallback begins at t=300.337, so lower CFL does not resolve instability.
- **qposim** · 2026-08-25 10:52 · [uncommitted] Add devlog.md with recent Git history backfill.
- **qposim** · 2026-08-04 22:24 · [1550789] Complete RaptorP LC Stage 1 pilot validation
- **qposim** · 2026-08-04 18:24 · [831eab1] Add pilot convergence analysis tooling
- **qposim** · 2026-08-04 18:21 · [c666f4d] Add isolated RaptorP light-curve pilot workflow
- **qposim** · 2026-07-11 23:19 · [18ab54c] Port RMHD jet setup to gPLUTO GPU
