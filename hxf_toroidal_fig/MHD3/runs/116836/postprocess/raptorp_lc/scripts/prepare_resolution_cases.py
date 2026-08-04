#!/usr/bin/env python3
"""Create the required resolution cases only after FOV convergence is known."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    selection = json.loads((ROOT / "pilot" / "analysis" / "fov_selection.json").read_text())
    if not selection["gate_passed"]:
        raise SystemExit("FOV convergence did not pass")
    fov = selection["selected_fov"]
    cases = []
    # The selected-FOV 320 x 512 cases were already completed as the FOV
    # matrix.  Only submit the two new resolutions, then reuse those records
    # in analyze_pilots.py.
    for resolution in ((160, 256), (480, 768)):
        for time, snapshot_index in ((400, 130), (600, 330), (800, 530)):
            width, height = resolution
            cases.append({
                "case_id": f"resolution_{width}x{height}_t{time:04d}",
                "time": time, "snapshot_index": snapshot_index, "volume_mask": "kink-only",
                "tracer_mask": "off", "fov": fov, "resolution": list(resolution),
                "target": "RAPTOR_LC4_IMAGE", "roles": ["resolution"],
            })
    output = ROOT / "config" / "pilot_resolution_cases.json"
    output.write_text(json.dumps(cases, indent=2, sort_keys=True) + "\n", encoding="ascii")
    print(output)


if __name__ == "__main__":
    main()
