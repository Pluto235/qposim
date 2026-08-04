#!/usr/bin/env python3
"""Generate immutable Stage 1 frequency/configuration metadata."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from frequency_contract import lc4_frequencies, sed50_frequencies, write_frequency_file


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=ROOT / "config")
    args = parser.parse_args()
    output = args.output
    output.mkdir(parents=True, exist_ok=True)

    lc4 = lc4_frequencies()
    sed50 = sed50_frequencies()
    write_frequency_file(output / "frequencies_lc4.txt", [value for _, value in lc4])
    write_frequency_file(output / "frequencies_sed50.txt", sed50)
    metadata = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "lc4": [{"label": label, "frequency_hz": f"{value:.30E}"}
                for label, value in lc4],
        "sed50_definition": "numpy.logspace(9, 18, 50)",
        "sed50_count": len(sed50),
    }
    (output / "frequency_contract.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="ascii")


if __name__ == "__main__":
    main()
