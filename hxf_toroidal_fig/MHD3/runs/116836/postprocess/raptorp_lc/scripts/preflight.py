#!/usr/bin/env python3
"""Fail closed on every Stage 1 production-gate prerequisite available pre-pilot."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    manifest = json.loads((ROOT / "manifests" / "snapshot_manifest_audit.json").read_text())
    assert manifest["passed"], manifest
    assert manifest["actual_records"] == 351, manifest
    physics = json.loads((ROOT / "config" / "physics.json").read_text())
    assert physics["frozen_physics"]["transfer_volume"] == "strict open interval 40 < z/R_j < 60"
    build = json.loads((ROOT / "provenance" / "build_metadata.json").read_text())
    assert all(item["returncode"] == 0 and item["binary_exists"]
               for item in build["targets"].values())
    cases = json.loads((ROOT / "config" / "pilot_phase1_cases.json").read_text())
    assert all(case["target"] == "RAPTOR_LC4_IMAGE" for case in cases)
    assert all(case["volume_mask"] in {"whole", "kink-only"} for case in cases)
    assert not list((ROOT / "production" / "epochs").glob("*/COMPLETE.json"))
    print("Stage 1 preflight passed: manifest, immutable physics, builds, and pilot-only scope")


if __name__ == "__main__":
    main()
