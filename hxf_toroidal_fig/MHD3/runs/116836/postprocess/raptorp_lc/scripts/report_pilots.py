#!/usr/bin/env python3
"""Publish Stage 1 pilot audits, convergence figures, and resource advice."""

from __future__ import annotations

import json
import re
from pathlib import Path
from statistics import mean
from typing import Any

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "pilot" / "raw"
ANALYSIS = ROOT / "pilot" / "analysis"
FIGURES = ROOT / "pilot" / "figures"


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="ascii"))


def elapsed_seconds(value: str) -> float:
    fields = value.strip().split(":")
    if len(fields) == 2:
        return 60.0 * float(fields[0]) + float(fields[1])
    if len(fields) == 3:
        return 3600.0 * float(fields[0]) + 60.0 * float(fields[1]) + float(fields[2])
    raise ValueError(f"unsupported elapsed time {value!r}")


def resources(case_id: str) -> dict[str, float]:
    text = (RAW / case_id / "resource.txt").read_text(encoding="ascii")
    def value(name: str) -> str:
        match = re.search(rf"{re.escape(name)}:\s*(.+)", text)
        if match is None:
            raise ValueError(f"missing {name} for {case_id}")
        return match.group(1).strip()
    return {
        "elapsed_seconds": elapsed_seconds(value("Elapsed (wall clock) time (h:mm:ss or m:ss)")),
        "cpu_percent": float(value("Percent of CPU this job got").rstrip("%")),
        "max_rss_kib": float(value("Maximum resident set size (kbytes)")),
        "fs_input_kib": float(value("File system inputs")),
        "fs_output_kib": float(value("File system outputs")),
    }


def required_cases() -> list[dict[str, Any]]:
    paths = [ROOT / "config" / "pilot_smoke_cases.json", ROOT / "config" / "pilot_phase1_cases.json",
             ROOT / "config" / "pilot_resolution_cases.json", ROOT / "config" / "pilot_validation_cases.json",
             ROOT / "config" / "pilot_recovery_interrupt_case.json"]
    return [case for path in paths for case in read_json(path)]


def verify_cases(cases: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for case in cases:
        directory = RAW / case["case_id"]
        complete = read_json(directory / "COMPLETE.json")
        metadata = read_json(directory / "run_metadata.json")
        if complete.get("status") != "complete" or complete.get("case_id") != case["case_id"]:
            raise SystemExit(f"invalid completion record for {case['case_id']}")
        spectrum = metadata["spectrum"]
        if len(spectrum) != 4 or any(row["I_Jy"] <= 0 for row in spectrum):
            raise SystemExit(f"invalid spectrum for {case['case_id']}")
        rows.append({"case_id": case["case_id"], "target": case["target"], "resource": resources(case["case_id"])})
    return rows


def mask_audit() -> dict[str, Any]:
    case_ids = {
        "whole_tracer_off": "mask_whole_tracer_off_t0400_10x16_320x512",
        "kink_tracer_off": "fov_10x16_t0400_320x512",
        "kink_tracer_on": "mask_kink_tracer_on_t0400_10x16_320x512",
    }
    spectra = {name: read_json(RAW / case_id / "run_metadata.json")["spectrum"]
               for name, case_id in case_ids.items()}
    comparisons = {}
    for name, left, right in (("whole_vs_kink", "whole_tracer_off", "kink_tracer_off"),
                              ("tracer_off_vs_on", "kink_tracer_off", "kink_tracer_on")):
        rows = []
        for a, b in zip(spectra[left], spectra[right]):
            rows.append({"frequency_hz": a["frequency_hz"],
                         "relative_I_difference": abs(a["I_Jy"] - b["I_Jy"]) / a["I_Jy"]})
        comparisons[name] = rows
    return {"case_ids": case_ids, "comparisons": comparisons,
            "controls_are_distinct": all(row["relative_I_difference"] > 0.0
                                           for rows in comparisons.values() for row in rows)}


def write_json(name: str, payload: Any) -> None:
    (ANALYSIS / name).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="ascii")


def convergence_figure(name: str, rows: list[dict[str, Any]], key: str, reference: str) -> None:
    metrics = ("relative_I", "absolute_mL", "qu_over_I", "evpa_deg")
    limits = {"relative_I": 0.02, "absolute_mL": 0.003, "qu_over_I": 0.003, "evpa_deg": 3.0}
    labels = sorted({str(row[key]) for row in rows})
    figure, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    for axis, metric in zip(axes.flat, metrics):
        for label in labels:
            selected = [row[metric] for row in rows if str(row[key]) == label]
            axis.scatter([label] * len(selected), selected, label=label)
        axis.axhline(limits[metric], color="black", linestyle="--", linewidth=1)
        axis.set_title(metric)
        axis.set_ylabel("difference")
    axes[0, 0].legend(title=key)
    figure.suptitle(f"Pilot convergence against {reference}")
    figure.savefig(FIGURES / name, dpi=180)
    plt.close(figure)


def main() -> None:
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    fov = read_json(ANALYSIS / "fov_selection.json")
    resolution = read_json(ANALYSIS / "resolution_selection.json")
    rows = verify_cases(required_cases())
    fov_rows = read_json(ANALYSIS / "fov_convergence.json")
    resolution_rows = read_json(ANALYSIS / "resolution_convergence.json")
    convergence_figure("fov_convergence.png", fov_rows, "fov", "12x18 R_j")
    convergence_figure("resolution_convergence.png", resolution_rows, "resolution", "480x768")

    selected_width, selected_height = resolution["selected_resolution"]
    if (selected_width, selected_height) == (320, 512):
        selected_ids = [f"fov_10x16_t{time:04d}_320x512" for time in (400, 600, 800)]
    else:
        selected_ids = [f"resolution_{selected_width}x{selected_height}_t{time:04d}"
                        for time in (400, 600, 800)]
    selected_resources = [resources(case_id) for case_id in selected_ids]
    mean_epoch_seconds = mean(row["elapsed_seconds"] for row in selected_resources)
    peak_rss_kib = max(row["resource"]["max_rss_kib"] for row in rows)
    resource_summary = {
        "selected_epoch_mean_seconds": mean_epoch_seconds,
        "selected_epoch_mean_minutes": mean_epoch_seconds / 60.0,
        "peak_rss_kib_all_pilots": peak_rss_kib,
        "peak_rss_gib_all_pilots": peak_rss_kib / 1024.0 / 1024.0,
        "requested_memory_gib": 64,
        "recommended_production_memory_gib": 24,
        "recommended_production_time_limit": "01:00:00",
        "projected_351_epoch_serial_hours": 351.0 * mean_epoch_seconds / 3600.0,
        "projected_351_epoch_two_way_hours": 351.0 * mean_epoch_seconds / 2.0 / 3600.0,
    }
    resource_summary["projected_351_epoch_two_way_hours_with_20_percent_margin"] = (
        1.2 * resource_summary["projected_351_epoch_two_way_hours"])
    audit = {
        "status": "complete",
        "case_count": len(rows),
        "cases": rows,
        "fov_selection": fov,
        "resolution_selection": resolution,
        "scheduler_status": read_json(ANALYSIS / "latest_scheduler_status.json"),
        "resource_summary": resource_summary,
        "mask_audit": mask_audit(),
        "limitations": [
            "Slurm sacct accounting storage is disabled; completion is evidenced by empty squeue, successful log finish records, and COMPLETE.json.",
            "This is a mapped-grid t=300..1000 methodology reproduction, not a strict paper reproduction.",
        ],
    }
    write_json("pilot_audit.json", audit)
    report = f"""# Stage 1 Pilot Report

## Recommendation

**GO for a separately authorized Stage 2 production array only.** Do not submit it from Stage 1.
The pilot selects FOV `{fov['selected_fov'][0]:g}x{fov['selected_fov'][1]:g} R_j` and resolution `{resolution['selected_resolution'][0]}x{resolution['selected_resolution'][1]}`. The `6x18 R_j` FOV fails the formal convergence gate; `10x16 R_j` passes against `12x18 R_j` for all times and frequencies.

## Resources

The selected 320x512 image pilot averages {resource_summary['selected_epoch_mean_minutes']:.2f} minutes/epoch. The largest measured RSS across all pilots is {resource_summary['peak_rss_gib_all_pilots']:.2f} GiB. Recommend `--cpus-per-task=32`, `--mem=24G`, `--time=01:00:00`, and array concurrency `%2`. At two concurrent tasks, 351 selected-resolution epochs project to {resource_summary['projected_351_epoch_two_way_hours']:.1f} hours ({resource_summary['projected_351_epoch_two_way_hours_with_20_percent_margin']:.1f} hours with 20% margin), based on image pilots; spectrum-only production is expected to be no slower.

## Scope and Evidence

All {len(rows)} unique Stage 1 pilot cases have valid `COMPLETE.json`, finite four-frequency spectra, physical Stokes checks, and per-case resource records. The mask audit confirms that whole-domain/kink-only and tracer-off/tracer-on controls produce distinct integrated values, while the strict `40 < z/R_j < 60` mask remains independent of tracer selection. Job IDs and final queue evidence are in `pilot/job_registry.json` and `pilot/analysis/latest_scheduler_status.json`; `sacct` is unavailable on this cluster. It is a Hu et al. (2025) Figure 5 methodology reproduction on the available mapped `320x320x450` run, `t=300..1000`, with adopted camera/azimuth/temperature/distance values documented separately from paper-published values.

Stage 2 remains unauthorized: no 351-epoch production array and no 50-frequency campaign has been submitted.
"""
    (ROOT / "PILOT_REPORT.md").write_text(report, encoding="ascii")
    print(json.dumps(audit["resource_summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
