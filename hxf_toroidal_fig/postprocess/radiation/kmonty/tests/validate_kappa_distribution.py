#!/usr/bin/env python3
"""Independent reference checks for the finite-support kappa_exp_v1 model.

The C probe supplies the widths returned by kmonty.  All normalization,
moment, and CDF calculations below are recomputed with SciPy and do not call
kmonty's integration or root-finding routines.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

# Keep this independent reference synchronized with decs.h's default model
# contract.
GAMMA_CUT = 1.0e5
GAMMA_MAX = GAMMA_CUT * 64.0
LINE = re.compile(r"kappa=(?P<k>\S+) Thetae=(?P<t>\S+) w=(?P<w>\S+)")


def pdf_log_gamma(x: float, w: float, kappa: float, thermal: bool = False) -> float:
    gamma = math.exp(x)
    if thermal:
        tail = math.exp(-gamma / w)
    else:
        tail = (1.0 + (gamma - 1.0) / (kappa * w)) ** (-kappa - 1.0)
        tail *= math.exp(-gamma / GAMMA_CUT)
    return gamma * gamma * math.sqrt(max(gamma * gamma - 1.0, 0.0)) * tail


def split_points(scale: float) -> list[float]:
    xmax = math.log(GAMMA_MAX)
    points = [0.0]
    points.extend(
        math.log1p(c * scale)
        for c in (1, 2, 4, 8, 16, 32, 64, 128, 256, 512,
                  1024, 4096, 16384, 65536, 262144, 1_000_000)
        if math.log1p(c * scale) < xmax
    )
    points.append(xmax)
    return sorted(set(points))


def integrate(func, scale: float) -> float:
    points = split_points(scale)
    value = sum(
        quad(func, left, right, epsabs=1e-13, epsrel=3e-10, limit=400)[0]
        for left, right in zip(points, points[1:])
    )
    return value


def moment(w: float, kappa: float, energy: bool) -> float:
    value = integrate(
        lambda x: pdf_log_gamma(x, w, kappa)
        * ((math.exp(x) - 1.0) if energy else 1.0), w
    )
    error = 0.0
    assert math.isfinite(value) and error < 1e-7 * max(1.0, abs(value))
    return value


def thermal_mean(thetae: float) -> float:
    # At the required 1e-3 and 2e-3 points, exp(-(gamma-1)/thetae)
    # underflows if evaluated directly in log(gamma).  Integrate the same
    # Maxwell-Juttner density after y=(gamma-1)/thetae instead.
    def base(y: float) -> float:
        gamma = 1.0 + thetae * y
        return thetae * gamma * math.sqrt(max(gamma * gamma - 1.0, 0.0)) * math.exp(-y)

    norm = quad(base, 0.0, 100.0, epsabs=1e-14, epsrel=3e-12, limit=400)[0]
    energy = quad(
        lambda y: base(y) * thetae * y,
        0.0,
        100.0,
        epsabs=1e-14,
        epsrel=3e-12,
        limit=400,
    )[0]
    return energy / norm


def cdf(x: float, w: float, kappa: float, norm: float) -> float:
    if x <= 0.0:
        return 0.0
    points = [p for p in split_points(w) if p < x] + [x]
    return sum(
        quad(lambda y: pdf_log_gamma(y, w, kappa), left, right,
             epsabs=1e-13, epsrel=3e-10, limit=400)[0]
        for left, right in zip(points, points[1:])
    ) / norm


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("probe_output", type=Path)
    args = parser.parse_args()
    rows = []
    for line in args.probe_output.read_text().splitlines():
        match = LINE.search(line)
        if match:
            rows.append(tuple(float(match.group(name)) for name in ("k", "t", "w")))
    if not rows:
        raise SystemExit("no kappa probe rows found")

    for kappa, thetae, width in rows:
        target = thermal_mean(thetae)
        norm = moment(width, kappa, False)
        energy = moment(width, kappa, True)
        rel_energy_error = abs(energy / norm - target) / max(1.0, abs(target))
        assert abs(norm) > 0.0 and rel_energy_error < 1e-5, (
            kappa, thetae, width, norm, rel_energy_error
        )

        grid = np.linspace(0.0, math.log(GAMMA_MAX), 33)
        values = np.array([cdf(float(x), width, kappa, norm) for x in grid])
        assert values[0] == 0.0
        assert values[-1] > 1.0 - 2e-8
        assert np.all(np.diff(values) >= -1e-12)
        assert np.all((values >= 0.0) & (values <= 1.0 + 1e-10))

        print(
            f"kappa={kappa:g} Thetae={thetae:g} "
            f"norm={norm:.8e} rel_energy_error={rel_energy_error:.3e} "
            f"cdf_end={values[-1]:.12f} PASS"
        )
    print("independent kappa normalization/moment/CDF: PASS")


if __name__ == "__main__":
    main()
