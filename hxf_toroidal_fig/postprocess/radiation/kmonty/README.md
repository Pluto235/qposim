# kmonty

`kmonty` is a Monte Carlo radiation transport program derived from the public
kmonty code. It generates synchrotron superphotons, propagates them through a
simulation domain, and includes Compton upscattering.

When using this code, cite Dolence et al. (2009),
<https://arxiv.org/abs/0909.0708>, and the relevant kmonty/RaptorP application
paper.

## Build

The code requires an MPI C compiler and GSL. On PDC, the makefile automatically
uses the GSL installation under `../../raptorp-deps/prefix/usr`.

```bash
make pluto
```

For a synchrotron-only control built from the same PLUTO reader and electron
distribution:

```bash
make pluto-no-compton
```

This writes `kmonty_no_compton` without changing the default `kmonty`
executable.

Other retained model targets are:

```bash
make sphere
make bhac
```

The checked-in targets build from `src/`. The `src_var_kappa/` directory is an
upstream work-in-progress copy and is not selected by the makefile.

The main build remains the validated thermal/synchrotron baseline. Defining
`KAPPA=1 THERMAL=0` enables the energy-matched relativistic-κ branches. The
width `w` is solved numerically per electron temperature by matching the
truncated κ mean kinetic energy to the truncated Maxwell-Jüttner mean; an
unrepresentable cell fails closed. `kappa_synch` may be set to 4 or 8 at
compile time and must be greater than 3.

Defining `EXTERNAL_DUST=1` enables a finite-torus dust seed stream. The torus
has `R=6.5e18 cm`, cross-section radius `0.1 R`, centre offset `3e18 cm`,
`T=1000 K`, and `L=0.6*6.7e45 erg/s`. Each external packet is sampled from a
Planck spectrum and a torus surface point; its energy weight represents
`L_DT*T_unit` and is divided by the MPI-wide packet count. The first hook uses
a fixed `DUST_PACKET_FRACTION=0.5` allocation within the existing emitting
cells. This is an explicitly documented volume-injection closure, not yet a
boundary ray-trace or a full local energy-density solution, and must not be
interpreted as a production EC model without the external-field validation
report.

When external dust is enabled, `components_<index>.dat` records source and
scattering-order contributions (synchrotron/dust × nscatt 0/1/2+). Unscattered
dust packets are kept separate and are not labelled EC.

## PLUTO Run

Run from this directory because the PLUTO reader expects `grid.out` in the
current working directory:

```bash
ln -s /path/to/run/grid.out grid.out
export OMP_NUM_THREADS=48
mpirun -np 1 ./kmonty model_pluto.in /path/to/run/data.0096.dbl 96
```

The output is `spec_96.dat`. The final argument is only the output index; the
number of generated superphotons and physical scaling are read from
`model_pluto.in`.

Set `KMONTY_SEED` to an integer for repeatable diagnostic runs. Keep the same
OpenMP/MPI layout when comparing runs because each worker has its own random
number stream.

`model_smoke.in` is a bounded 100-superphoton configuration for deployment
checks; it is not intended for production spectra.

For a Slurm run:

```bash
sbatch job.slurm /path/to/run/data.0096.dbl 96
```

## Configuration

Runtime parameters are in `model_pluto.in`. Compile-time physics and detector
settings are in `decs.h`. The supplied build uses:

- Monte Carlo Compton scattering enabled;
- Cartesian flat spacetime;
- thermal electrons (`THERMAL=1`);
- `N_THBINS=9` observer-angle bins;
- seed frequencies from `1e8` to `1e20` Hz;
- no z trimming (`TRIM=0`).

When `TRIM=0`, `TRIM_ZI` and `TRIM_ZF` in `model_pluto.in` are not read. The
PLUTO reader consumes the first eight primitive arrays and does not currently
apply a tracer mask.

The `hotcross.dat` lookup table is generated on first use when it is absent.
It and all `spec_*.dat` files are ignored build/run artifacts.
Delete `hotcross.dat` before running a build with a different electron
distribution or modified Compton cross-section implementation so that the
table is regenerated for the active physics.
