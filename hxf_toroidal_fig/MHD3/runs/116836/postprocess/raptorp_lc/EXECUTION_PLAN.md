# RaptorP Light-Curve Execution Plan

## 1. Objective

Implement and validate a Hu et al. (2025) Figure 5 methodology reproduction on
the available MHD3 mapped-grid run, using the existing snapshots only:

- simulation interval: `300 <= t/(R_j/c) <= 1000`;
- cadence: `Delta t = 2 R_j/c`;
- epochs: 351;
- observer inclination: 10 degrees;
- observer azimuth: 0 degrees;
- electron distribution: PIC-informed variable kappa distribution;
- emitting/propagating volume: `40 < z/R_j < 60`;
- four exact representative frequencies;
- integrated polarized Stokes light curves, polarization diagnostics, PSDs,
  correlations, and Q-U evolution.

This is a **mapped-grid, truncated-interval methodology reproduction**, not a
strict numerical reproduction of the paper. The paper uses a uniform
`500 x 500 x 600` simulation and `t=300..1100`; the available run uses a mapped
`320 x 320 x 450` grid and stops at `t=1000`.

## 2. Authorization Boundary

The work is divided into two separately authorized stages.

### Stage 1: authorized after the user says to begin

Stage 1 may:

1. create the isolated LC worktree and directory structure;
2. implement the kink-volume mask and multi-frequency/spectrum-only support;
3. build and test all required binaries;
4. generate and audit the 351-epoch snapshot manifest;
5. submit, monitor, retry, and complete the pilot matrix;
6. generate the pilot convergence and resource report.

Stage 1 must stop after the pilot report. It must not submit the 351-epoch
production array.

### Stage 2: requires a new explicit user approval

Only after the user reviews the pilot report may Stage 2 submit the production
array:

```text
351 epochs x 4 frequencies, 32 CPU cores per epoch, maximum 2 concurrent tasks
```

The separate 50-frequency SED campaign is also outside Stage 1 and requires
later approval.

## 3. Evidence and Provenance

### Paper source

- Local paper text: `/home/luoji/qposim/documents/2511.03140.md`
- Local paper PDF: `/home/luoji/qposim/documents/2511.03140.pdf`
- Relevant locations:
  - Section 2.2: fast-light polarized transfer in Minkowski spacetime;
  - Section 4: `t=300..1100`, cadence 2, 50 logarithmic frequencies, 10-degree
    inclination;
  - Section 4.2: kink region and four representative frequencies;
  - Section 4.3: total polarization degree and intensity-PD correlations;
  - Appendix E: scaled-chi-square QPO significance calculation.

### Validated RaptorP baseline

- Baseline worktree:
  `/home/luoji/qposim/hxf_toroidal_fig/MHD3/runs/116705/postprocess/raptorp_t400/source`
- Upstream commit:
  `8e71ba667cff419522fb61ef2fc1662aa59a1902`
- The baseline contains validated local changes for safe absolute paths,
  argument handling, and runtime tracer selection.
- Do not modify this baseline in place.
- Do not use `/home/luoji/qposim/tools/RaptorP` as the LC source baseline; that
  checkout has additional unrelated local modifications.

Before changing the isolated LC source, save:

1. the upstream commit;
2. the complete baseline diff;
3. SHA256 for all source files and the existing executable;
4. a second patch containing only the new LC changes.

The final provenance model is:

```text
upstream commit + validated Figure 3 patch + isolated LC patch
```

## 4. Frozen Scientific Decisions

| Item | Frozen value |
|---|---|
| Main deliverable | Figure 5 four-frequency LC/QPO methodology reproduction |
| Time interval | `t=300..1000` |
| Sampling | `Delta t=2`, 351 epochs |
| Inclination | 10 degrees |
| Azimuth | 0 degrees |
| Spacetime/time treatment | Minkowski, fast-light |
| Electron distribution | Variable kappa, turbulence PIC prescription |
| Temperature closure | `R_LOW=1`, `R_HIGH=5` |
| Distance | 16.8 Mpc |
| Effective RaptorP code mass | `6.5e9 M_sun` |
| Physical paper mass | `6.5e7 M_sun` |
| Length/time unit | `R_j=100 R_g`, `R_j/c=32023.8 s` |
| Tracer | Off; include jet, mixing layer, and environment inside the slab |
| Transfer volume | Strict open interval `40 < z/R_j < 60` |
| Outside slab | Emission, absorption, rotation, and conversion coefficients all zero |
| Production image output | Disabled; integrated spectrum only |
| Pilot image output | Enabled for visual and numerical QA |
| Production parallelism | 32 CPU cores/task, at most 2 tasks concurrently, no GPU |

The paper does not explicitly publish the LC azimuth, camera, distance, or
`R_LOW/R_HIGH`. The adopted values preserve continuity with the validated
Figure 3 reproduction and must be labelled as adopted code/reproduction
assumptions rather than explicit paper values.

## 5. Frequency Contract

### Four-frequency LC target

Use the exact named frequencies, not nearest points from the 50-point SED grid:

| Label | Frequency |
|---|---:|
| 1.5 GHz | `1.500000000e9 Hz` |
| 219 GHz | `2.190000000e11 Hz` |
| 521 nm | compute and record `c / (521 nm)`, approximately `5.75e14 Hz` |
| 1.2 keV | compute and record `(1.2 keV) / h`, approximately `2.90e17 Hz` |

Generate the optical and X-ray frequencies from one constants module. Do not
duplicate rounded literals across shell, C, and Python code. Record both the
exact generated frequency and the display label.

### Future 50-frequency SED target

Generate exactly:

```text
numpy.logspace(9, 18, 50)
```

The four-frequency LC and the 50-frequency SED are independent frequency
products. Do not silently substitute nearest SED-grid points into Figure 5.

## 6. Input Snapshot Contract

### Source run

- Run directory:
  `/home/luoji/qposim/hxf_toroidal_fig/MHD3/runs/116836`
- Grid:
  `320 x 320 x 450`, mapped Cartesian, `x1,x2=[-5,5]`, `x3=[0,60]`
- Primitive order:
  `rho vx1 vx2 vx3 Bx1 Bx2 Bx3 prs tr1`
- Expected snapshot size:
  `3,317,760,000 bytes`

### Epoch selection

Read `dbl.out`; do not infer times from filenames alone. Select exactly the
even indices:

```text
30, 32, 34, ..., 730
```

These correspond to the intended times:

```text
300, 302, 304, ..., 1000
```

The manifest must contain exactly 351 rows and at least:

```text
epoch_ordinal
snapshot_index
actual_time_Rj_c
target_time_Rj_c
time_error
snapshot_path
resolved_path
is_symlink
size_bytes
primitive_order
validation_status
```

Manifest acceptance criteria:

- exactly 351 records;
- first/last index `30/730`;
- first/last intended time `300/1000`;
- absolute time error below `0.01 R_j/c`;
- every symlink resolves;
- every snapshot has the expected size;
- all `dbl.out` records have identical byte order, file mode, and primitive
  order;
- no destructive changes to any snapshot or symlink.

## 7. Work Directory Contract

Use:

```text
/home/luoji/qposim/hxf_toroidal_fig/MHD3/runs/116836/postprocess/raptorp_lc/
```

Required layout:

```text
raptorp_lc/
|-- EXECUTION_PLAN.md
|-- README.md
|-- assumptions.md
|-- config/
|-- source/
|-- provenance/
|-- manifests/
|-- scripts/
|-- tests/
|-- pilot/
|   |-- raw/
|   |-- analysis/
|   |-- figures/
|   `-- logs/
|-- production/
|   |-- epochs/
|   |-- logs/
|   |-- status/
|   `-- failed/
|-- aggregate/
|-- analysis/
|   |-- psd/
|   |-- correlation/
|   `-- q_u/
`-- figures/
```

Do not place generated outputs in the validated Figure 3 directories.

## 8. Isolated Source Design

Create one isolated source tree and build these targets from the same physics
source:

| Target | Frequencies | Images | Purpose |
|---|---:|---:|---|
| `RAPTOR_LC4_IMAGE` | 4 | on | pilot image and spectrum QA |
| `RAPTOR_LC4_SPEC` | 4 | off | 351-epoch production |
| `RAPTOR_SED50_IMAGE` | 50 | on | optional one-frame test only |
| `RAPTOR_SED50_SPEC` | 50 | off | future SED production |

Requirements:

1. Make frequency count a build parameter; do not maintain divergent source
   copies for 4 and 50 frequencies.
2. Read explicit `frequencies.txt` files.
3. Record frequency count, image-output state, git/provenance hashes, compiler,
   and build command in machine-readable metadata.
4. The image and spectrum-only builds must produce identical integrated
   spectra for the same pilot case within floating-point/text-output tolerance.
5. Keep tracer control independent from the kink volume mask.
6. Implement the volume mask in the fluid/transfer path so that outside
   `40 < z < 60` all plasma coefficients are zero.
7. A whole-domain mode may exist for pilot controls, but production must refuse
   to run unless `kink-only` is explicit.
8. Use structured parsing and explicit validation; do not use brittle string
   replacement to generate scientific configuration files.

## 9. Output and Recovery Contract

Each epoch directory must be deterministic, for example:

```text
production/epochs/t0300_i0030/
```

Each completed epoch must contain:

```text
spectrum.dat
run_metadata.json
stdout/stderr or linked log references
executable.sha256
source/provenance reference
COMPLETE.json
```

Rules:

- write outputs into a temporary epoch directory;
- validate four frequencies and all Stokes values;
- atomically rename into the final epoch directory;
- write `COMPLETE.json` last;
- if a valid completion record exists, skip the epoch;
- never overwrite a valid completed result;
- move an invalid/incomplete prior directory under `production/failed/` with a
  timestamp before retrying;
- aggregation writes a temporary output and atomically publishes only after
  full validation;
- final figures read only the validated aggregate, never ad hoc log scans.

## 10. Pilot Matrix

Use `t=400,600,800` and calculate all four frequencies in every case.

### 10.1 Mask semantics at t=400

Run:

1. whole domain, tracer off;
2. kink-only, tracer off;
3. kink-only, tracer on.

Production physics is case 2. Cases 1 and 3 prove that the volume mask and
tracer mask are independent.

### 10.2 FOV convergence

At `t=400,600,800`, compare:

| FOV | Role |
|---|---|
| `6 x 18 R_j` | Figure 3 display-crop control |
| `10 x 16 R_j` | full-slab candidate |
| `12 x 18 R_j` | margin/reference candidate |

Choose the smallest FOV converged against the larger reference. The formal LC
must target integrated slab flux rather than visually matching the Figure 3
crop.

### 10.3 Resolution convergence

At the selected FOV and `t=400,600,800`, compare:

```text
160 x 256
320 x 512
480 x 768
```

Use the highest resolution as the reference. Select `320 x 512` if it passes
all thresholds; otherwise select `480 x 768`.

### 10.4 Convergence thresholds

For every pilot time and frequency:

```text
relative integrated-I difference < 2 percent
absolute linear-polarization-degree difference < 0.003
sqrt((Delta Q)^2 + (Delta U)^2) / I_reference < 0.003
EVPA circular difference < 3 degrees when m_L >= 0.005
```

Do not enforce EVPA convergence where `m_L < 0.005`, because EVPA is unstable
when integrated linear polarization is nearly zero.

### 10.5 Additional pilot checks

- all image arrays and spectra finite;
- integrated `I > 0`;
- `sqrt(Q^2+U^2) <= I` and `abs(V) <= I`, allowing documented floating-point
  tolerance;
- inspect intensity at all camera boundaries to detect clipped flux;
- image-build spectrum equals spectrum-only-build spectrum;
- repeat one epoch twice to verify deterministic numerical output;
- interrupt/restart a test epoch to verify recovery behavior;
- record wall time, CPU utilization, peak RSS, and filesystem I/O;
- ensure peak RSS is below 80 percent of requested Slurm memory;
- ensure epoch wall time is below 70 percent of the Slurm time limit.

The pilot is expected to contain roughly 18 non-duplicate cases after reusing
overlapping configurations.

## 11. Production Gate

Do not submit production until all of the following are true and the user has
explicitly approved the pilot report:

1. the 351-row manifest passes;
2. all frozen physical values are present in immutable config metadata;
3. the kink and tracer controls pass;
4. FOV and resolution pass the accepted convergence thresholds;
5. all pilot Stokes/image values pass physical and finite-value checks;
6. image/spectrum builds agree;
7. deterministic rerun and interrupted recovery tests pass;
8. measured memory and wall time fit the Slurm request with the required
   margins;
9. the pilot report states the selected FOV, resolution, expected total wall
   time, and residual numerical uncertainty.

The planned production resources, subject to pilot measurement, are:

```bash
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=24G-or-32G
#SBATCH --time=04:00:00
#SBATCH --array=0-350%2
```

Do not request a GPU. RaptorP is an OpenMP CPU calculation in this workflow.
The node has 72 logical CPUs; two 32-core tasks leave capacity for the system
and light post-processing.

## 12. Production Completion Gate

Production is complete only when:

- 351/351 epochs have valid `COMPLETE.json` records;
- 1404/1404 epoch-frequency records exist;
- all frequencies match the four-frequency contract;
- all Stokes values pass finite/physical checks;
- no epoch directory is silently duplicated or overwritten;
- aggregate time ordering and cadence pass;
- audit JSON records zero missing or invalid rows.

Do not generate files labelled `final` before this gate passes.

## 13. Aggregate Data Products

The authoritative table stores raw Stokes fluxes in Jy:

```text
epoch_ordinal
snapshot_index
time_Rj_c
time_seconds
frequency_label
frequency_hz
I_Jy
Q_Jy
U_Jy
V_Jy
run_status
```

Derived values are calculated only after spatial integration:

```text
P_Jy = sqrt(Q_Jy^2 + U_Jy^2)
m_L = P_Jy / I_Jy
EVPA_rad = 0.5 * atan2(U_Jy, Q_Jy)
```

Save both EVPA modulo 180 degrees and a separately unwrapped EVPA time series.
Never average pixel polarization degrees to obtain total polarization degree.

Raw Jy is the scientific source of truth. For Figure 5-style display, compute
per-frequency z-scores over the complete available 351-point interval:

```text
I_norm = (I - mean(I_300:1000)) / std(I_300:1000, ddof=0)
```

Name and label all normalized products `t300_1000_dt2`. Do not imply that their
normalization is identical to the paper's unpublished formula.

## 14. Time-Series Analysis

### 14.1 Main Figure 5-style panels

Generate:

1. four-frequency normalized intensity LC;
2. intensity PSD for 219 GHz, 521 nm, and 1.2 keV;
3. four-frequency total linear polarization degree;
4. polarization-degree PSD for 219 GHz, 521 nm, and 1.2 keV;
5. 521 nm intensity and PD for `t=428..456`, with `t=430..454` shaded;
6. 521 nm Q-U trajectory for `t=430..454`, colored by time.

Panel 5 uses two physical y axes:

- left: raw 521 nm intensity in Jy;
- right: linear polarization degree.

Do not invent the paper's unpublished arbitrary intensity rescaling. A
z-score overlay may be supplied as a supplemental figure.

### 14.2 Supplemental figures

- four-frequency raw Jy light curves;
- wrapped and unwrapped four-frequency EVPA time series;
- all four raw and processed PSDs;
- pilot FOV/resolution convergence figures.

### 14.3 Correlation analysis

For all four frequencies, using all 351 unsmoothed samples, calculate:

- Pearson `r` and two-sided P-value;
- Spearman rho as a supplemental robustness statistic;
- bootstrap confidence intervals with a fixed and recorded seed/method.

Use 219 GHz, 521 nm, and 1.2 keV as the main paper-comparison table. Retain
1.5 GHz in machine-readable and supplemental results.

### 14.4 PSD analysis

Paper-approximation track:

- intensity: Hann window, `nperseg=228`;
- polarization degree: Hann window, `nperseg=178`;
- `noverlap=50 percent`;
- constant detrending per segment;
- `scipy.signal.welch(..., scaling="density")`;
- separately save unnormalized PSD and the paper-style PSD divided by its
  lowest nonzero-frequency power.

Audit track:

- raw unsmoothed periodogram;
- Welch PSD before paper-style normalization;
- overlap/window sensitivity checks;
- complete algorithm and version metadata.

Always report frequency first in `c/R_j`. Convert with
`R_j/c = 32023.8 s` into Hz and periods in days. Do not repeat the paper's
internally inconsistent literal `mHz` label without explaining the unit issue.

### 14.5 QPO significance

Follow Appendix E while making missing implementation choices explicit:

1. fit `P proportional to f^(-alpha)` over `0.01 < f < 0.1 c/R_j`;
2. record all excluded candidate-QPO bins before calculating significance;
3. estimate the scaled-chi-square parameters from the smoothed PSD moments;
4. calculate local P-values using the regularized upper incomplete gamma
   expression in Appendix E;
5. also report a global/look-elsewhere significance across searched bins;
6. retain blind/no-exclusion and paper-candidate fits as separate outputs.

Because the available series has 351 rather than 401 samples, do not claim
one-to-one equality of PSD amplitude, lowest-frequency normalization, or
P-values. Compare QPO peak locations and qualitative behavior first.

## 15. Final Deliverables

Stage 1 deliverables:

- isolated, provenance-recorded LC source;
- four build targets and build/test logs;
- 351-row validated snapshot manifest;
- immutable scientific configs and frequency files;
- pilot raw images/spectra/logs;
- mask-semantics audit;
- FOV convergence table;
- resolution convergence table;
- image/spec equivalence and deterministic/recovery tests;
- selected production FOV and resolution;
- measured memory, per-epoch runtime, and projected campaign duration;
- pilot report with a clear GO/NO-GO recommendation;
- no production array submission.

Stage 2 deliverables, only after new approval:

- validated 351-epoch raw Stokes archive;
- aggregate and derived LC CSV/JSON products;
- six-panel main figure and supplemental figures;
- correlation, PSD, QPO, and Q-U products;
- final audit and reproduction documentation.

Future SED deliverables, separately approved:

- 50-frequency per-epoch spectra;
- every-five-frequencies binning into 10 characteristic bins;
- whole-region and kink-region time-averaged SED comparison.

## 16. Scientific Claim Boundary

Use this result name:

```text
Hu et al. (2025) Figure 5 methodology reproduction on the available
mapped-grid MHD3 run, t=300..1000
```

Do not use `strict reproduction`, `exact reproduction`, or `complete Figure 5
reproduction`.

Every report must state:

- mapped `320 x 320 x 450` versus paper uniform `500 x 500 x 600`;
- 351 samples over `300..1000` versus paper 401 over `300..1100`;
- adopted unpublished camera/azimuth/temperature/distance assumptions;
- fast-light limitation;
- whether QPO agreement is qualitative or statistically comparable.

## 17. Safety and Operational Rules

- Read repository `AGENTS.md` files before implementation.
- Preserve the dirty worktree and all user changes.
- Do not delete, move, rewrite, compress, or relink any MHD snapshot.
- Do not clean large run directories to create space.
- Do not modify the validated Figure 3 source/worktree in place.
- Do not overwrite successful pilot or production results.
- Do not submit the production array without a new explicit user message.
- Do not launch the 50-frequency campaign in Stage 1.
- Keep Slurm jobs under active supervision until they reach terminal states.
- Inspect failed logs and retry only after identifying and recording the cause.
- Avoid GPU requests for RaptorP LC jobs.
- Record all job IDs, terminal states, exit codes, elapsed time, and resource
  metrics in the pilot report.

## 18. Stage 1 Definition of Done

Stage 1 is done when, and only when:

1. the isolated implementation and all tests are complete;
2. all pilot jobs have terminal successful states, or failures are fully
   resolved and rerun;
3. the pilot convergence gates have an objective result;
4. a production FOV/resolution and Slurm resource request are recommended;
5. the pilot report and machine-readable audits are present;
6. the agent has stopped before production submission and asked the user for
   the separate Stage 2 authorization.
