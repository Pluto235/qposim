# Assumptions and Claim Boundary

The paper explicitly supplies fast-light polarized transfer in Minkowski
spacetime, a 10-degree viewing inclination, the `t=300..1100` cadence-two
series, a kink region `40 < z/R_j < 60`, and the four representative
frequencies.  It also describes the PIC-informed variable-kappa prescription.

The following are adopted reproduction/code assumptions inherited from the
validated Figure 3 setup, not values publicly tabulated by the paper: azimuth
0 degrees, distance 16.8 Mpc, `R_LOW=1`, `R_HIGH=5`, and the effective RaptorP
code mass `6.5e9 M_sun`.  The physical paper mass remains `6.5e7 M_sun`; the
effective code mass represents one PLUTO coordinate as `R_j=100 R_g` and gives
`R_j/c=32023.8 s`.

The available data ends at `t=1000`.  Stage 1 therefore audits 351 epochs
from `t=300..1000` but does not start an MHD continuation or any 50-frequency
campaign.  Any later statistical comparison must describe QPO agreement as
qualitative unless the shorter mapped-grid interval permits a justified
statistical comparison.
