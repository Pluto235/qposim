# κ + dust-seed implementation checkpoint

This directory records the first implementation checkpoint for the kmonty
κ-electron and external-dust extension. It is not a PKS 1510−089 fit.

Implemented in the nested kmonty source:

- `KAPPA=1 THERMAL=0`: one numerically energy-matched κ width per `Thetae`,
  shared by synchrotron, absorption, Compton cross sections and electron
  sampling. `kappa_synch=4` and `8` both build.
- `EXTERNAL_DUST=1`: finite torus geometry, Planck sampling, AGN-frame ray
  direction, fluid-tetrad transform, energy-normalised packet weight, and
  `seed_source` propagation through scattering.
- `components_<index>.dat`: source × scattering-order diagnostic output.

Validation run from `tools/RaptorP/kmonty`:

```text
make pluto -j4                                      PASS
make pluto ... -DTHERMAL=0 -DKAPPA=1 -Dkappa_synch=4 PASS
make pluto ... -DTHERMAL=0 -DKAPPA=1 -Dkappa_synch=8 PASS
tests/test_external_dust.c (EXTERNAL_DUST=0/1)       PASS
tests/test_kappa_matching.c (KAPPA=1)               PASS
```

The dust hook currently allocates half of the existing cell packet stream to
external seeds. This gives a reproducible luminosity-weighted packet
normalisation, but it is not a substitute for a validated boundary injection
or a full local anisotropic radiation-field integral. Snapshot production is
therefore intentionally not started at this checkpoint.

Remaining hard gates: analytic κ moment comparison over the full temperature
range, independent dust-ring energy-density and angular-distribution tests,
thermal/κ × dust-off/on regression on a transparent medium, and serial versus
OpenMP statistical agreement.
