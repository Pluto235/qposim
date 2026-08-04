"""Fast local checks for immutable Stage 1 contracts."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_frequency_contract() -> None:
    subprocess.run([sys.executable, str(ROOT / "scripts" / "generate_configs.py")], check=True)
    contract = json.loads((ROOT / "config" / "frequency_contract.json").read_text())
    labels = [item["label"] for item in contract["lc4"]]
    assert labels == ["1.5 GHz", "219 GHz", "521 nm", "1.2 keV"]
    assert len((ROOT / "config" / "frequencies_lc4.txt").read_text().splitlines()) == 4
    assert len((ROOT / "config" / "frequencies_sed50.txt").read_text().splitlines()) == 50


def test_manifest() -> None:
    subprocess.run([sys.executable, str(ROOT / "scripts" / "make_manifest.py")], check=True)
    with (ROOT / "manifests" / "snapshot_manifest.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 351
    assert rows[0]["snapshot_index"] == "30"
    assert rows[-1]["snapshot_index"] == "730"
    assert all(row["validation_status"] == "valid" for row in rows)


def test_strict_mask_is_in_transfer_path() -> None:
    source = (ROOT / "source" / "model.c").read_text()
    assert "RAPTORP_VOLUME_MASK" in source
    assert "40. &&" in source and "< 60." in source
    assert "get_fluid_params" in source


if __name__ == "__main__":
    test_frequency_contract()
    test_manifest()
    test_strict_mask_is_in_transfer_path()
    print("Stage 1 contract tests passed")
