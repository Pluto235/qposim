#!/usr/bin/env python3
"""Audit completed pilot outputs and evaluate FOV/resolution convergence."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
THRESHOLDS = {"relative_I": 0.02, "absolute_mL": 0.003, "qu_over_I": 0.003, "evpa_deg": 3.0}


def metadata(case_id: str) -> dict[str, Any]:
    path = ROOT / "pilot" / "raw" / case_id / "run_metadata.json"
    return json.loads(path.read_text(encoding="ascii"))


def stokes(record: dict[str, Any]) -> dict[float, dict[str, float]]:
    return {row["frequency_hz"]: row for row in record["spectrum"]}


def comparison(candidate: dict[str, Any], reference: dict[str, Any]) -> list[dict[str, Any]]:
    result = []
    for frequency, ref in stokes(reference).items():
        value = stokes(candidate)[frequency]
        ref_i = ref["I_Jy"]
        m_l = math.hypot(value["Q_Jy"], value["U_Jy"]) / value["I_Jy"]
        ref_m_l = math.hypot(ref["Q_Jy"], ref["U_Jy"]) / ref_i
        evpa = 0.5 * math.degrees(math.atan2(value["U_Jy"], value["Q_Jy"]))
        ref_evpa = 0.5 * math.degrees(math.atan2(ref["U_Jy"], ref["Q_Jy"]))
        evpa_difference = abs((evpa - ref_evpa + 90.0) % 180.0 - 90.0)
        row = {
            "frequency_hz": frequency,
            "relative_I": abs(value["I_Jy"] - ref_i) / ref_i,
            "absolute_mL": abs(m_l - ref_m_l),
            "qu_over_I": math.hypot(value["Q_Jy"] - ref["Q_Jy"], value["U_Jy"] - ref["U_Jy"]) / ref_i,
            "evpa_deg": evpa_difference,
            "reference_mL": ref_m_l,
        }
        row["passes"] = (row["relative_I"] < THRESHOLDS["relative_I"] and
                         row["absolute_mL"] < THRESHOLDS["absolute_mL"] and
                         row["qu_over_I"] < THRESHOLDS["qu_over_I"] and
                         (ref_m_l < 0.005 or evpa_difference < THRESHOLDS["evpa_deg"]))
        result.append(row)
    return result


def write_json(name: str, payload: Any) -> None:
    path = ROOT / "pilot" / "analysis" / name
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="ascii")


def fov_analysis(cases_path: Path) -> None:
    cases = json.loads(cases_path.read_text(encoding="ascii"))
    fov_cases = [case for case in cases if "fov" in case["roles"]]
    by_key = {(tuple(case["fov"]), case["time"]): case for case in fov_cases}
    reference_fov = (12.0, 18.0)
    output = []
    fov_passes: dict[tuple[float, float], bool] = {}
    for fov in ((6.0, 18.0), (10.0, 16.0), (12.0, 18.0)):
        passed = True
        for time in (400, 600, 800):
            candidate = by_key[(fov, time)]
            reference = by_key[(reference_fov, time)]
            rows = comparison(metadata(candidate["case_id"]), metadata(reference["case_id"]))
            passed &= all(row["passes"] for row in rows)
            for row in rows:
                output.append({"fov": list(fov), "time": time, "case_id": candidate["case_id"],
                               "reference_case_id": reference["case_id"], **row})
        fov_passes[fov] = passed
    choices = [fov for fov in ((6.0, 18.0), (10.0, 16.0), (12.0, 18.0)) if fov_passes[fov]]
    selection = {"thresholds": THRESHOLDS, "passes": {f"{x}x{y}": passed for (x, y), passed in fov_passes.items()},
                 "selected_fov": list(choices[0]) if choices else None,
                 "gate_passed": bool(choices)}
    write_json("fov_convergence.json", output)
    write_json("fov_selection.json", selection)
    with (ROOT / "pilot" / "analysis" / "fov_convergence.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)
    if not choices:
        raise SystemExit("no FOV passes convergence thresholds")


def resolution_analysis(cases_path: Path) -> None:
    selection = json.loads((ROOT / "pilot" / "analysis" / "fov_selection.json").read_text())
    if not selection["gate_passed"]:
        raise SystemExit("cannot analyze resolution before FOV gate")
    selected_fov = tuple(selection["selected_fov"])
    cases = json.loads(cases_path.read_text(encoding="ascii"))
    by_key = {(tuple(case["resolution"]), case["time"]): case for case in cases}
    output = []
    passes_320 = True
    for resolution in ((160, 256), (320, 512), (480, 768)):
        for time in (400, 600, 800):
            candidate = by_key[(resolution, time)]
            reference = by_key[((480, 768), time)]
            rows = comparison(metadata(candidate["case_id"]), metadata(reference["case_id"]))
            if resolution == (320, 512):
                passes_320 &= all(row["passes"] for row in rows)
            for row in rows:
                output.append({"fov": list(selected_fov), "resolution": list(resolution), "time": time,
                               "case_id": candidate["case_id"], "reference_case_id": reference["case_id"], **row})
    selection_out = {"selected_fov": list(selected_fov), "selected_resolution": [320, 512] if passes_320 else [480, 768],
                     "thresholds": THRESHOLDS, "mid_resolution_passes": passes_320, "gate_passed": True}
    write_json("resolution_convergence.json", output)
    write_json("resolution_selection.json", selection_out)
    with (ROOT / "pilot" / "analysis" / "resolution_convergence.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("fov", "resolution"))
    parser.add_argument("--cases", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "fov":
        fov_analysis(args.cases)
    else:
        resolution_analysis(args.cases)


if __name__ == "__main__":
    main()
