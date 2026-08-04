#!/usr/bin/env python3
"""Validate image/spec equivalence, determinism, and interrupted recovery."""

from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "pilot" / "raw"


def spectrum(case_id: str) -> list[list[float]]:
    path = RAW / case_id / "spectrum.dat"
    return [[float(value) for value in line.split()] for line in path.read_text(encoding="ascii").splitlines()
            if line and not line.startswith("#")]


def compare(left: str, right: str) -> dict[str, object]:
    left_rows, right_rows = spectrum(left), spectrum(right)
    if len(left_rows) != len(right_rows):
        raise SystemExit(f"row count mismatch: {left} {right}")
    maximum = 0.0
    for first, second in zip(left_rows, right_rows):
        for a, b in zip(first, second):
            maximum = max(maximum, abs(a - b) / max(abs(a), abs(b), 1.0))
            if not math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-300):
                raise SystemExit(f"spectrum mismatch: {left} {right}")
    return {"left": left, "right": right, "maximum_relative_difference": maximum, "passed": True}


def main() -> None:
    checks = [
        compare("smoke_kink_t0400_2x4_32x64", "validation_equivalence_spec_t0400_2x4_32x64"),
        compare("validation_deterministic_a_t0400_2x4_32x64", "validation_deterministic_b_t0400_2x4_32x64"),
    ]
    recovery = "validation_recovery_interrupt_t0400_10x16_320x512"
    failed = list((ROOT / "pilot" / "failed").glob(f"{recovery}_*"))
    if not failed or not (RAW / recovery / "COMPLETE.json").is_file():
        raise SystemExit("interrupted recovery evidence is incomplete")
    payload = {"checks": checks, "recovery_failed_attempts": [str(path) for path in failed], "passed": True}
    output = ROOT / "pilot" / "analysis" / "validation_checks.json"
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="ascii")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
