# gPLUTO GPU port

This directory ports the original `initial_setup/` PLUTO problem to the
OpenACC GPU build of gPLUTO.

## Build

```bash
export PLUTO_DIR=/home/luoji/opt/gPLUTO-870aa807
make -j
```

The makefile selects `linux.nvc++.acc.defs`, which uses NVIDIA `nvc++` and
OpenACC managed memory.

## Smoke test

```bash
./pluto -i pluto_smoke.ini -maxsteps 2 -no-write
```

## Production input

```bash
./pluto -i pluto.ini
```

The production input preserves the original `500 x 500 x 600` base grid.
gPLUTO currently uses a static grid, so the old Chombo AMR blocks have been
removed. Scaling from the measured smoke-test allocation, this grid needs
roughly 76 GiB for the fluid arrays alone. It does not fit on one 16 GB V100,
and four such GPUs do not leave enough aggregate memory after overhead.

For a practical single-V100 run, use the reduced `200 x 200 x 240` input:

```bash
./pluto -i pluto_v100.ini
```

This keeps the original domain and parameters while reducing all three grid
dimensions by a factor of 2.5. A one-step V100 test allocated 5237.92 MB.

`MACH` is retained for input compatibility with the original setup, although
the original `init.c` does not use it in the jet profile.
