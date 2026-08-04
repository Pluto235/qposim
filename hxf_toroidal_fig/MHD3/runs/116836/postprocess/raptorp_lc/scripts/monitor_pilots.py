#!/usr/bin/env python3
"""Capture queue and accounting state for every submitted pilot array."""

from __future__ import annotations

import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run(command: list[str]) -> str:
    return subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False).stdout


def main() -> None:
    registry = json.loads((ROOT / "pilot" / "job_registry.json").read_text(encoding="ascii"))
    observations = []
    for job in registry["jobs"]:
        job_id = job["job_id"]
        observations.append({
            "job_id": job_id,
            "squeue": run(["squeue", "-h", "-j", job_id, "-o", "%A_%a|%T|%M|%l|%R"]),
            "sacct": run(["sacct", "-X", "-n", "-P", "-j", job_id,
                           "--format=JobIDRaw,State,ExitCode,Elapsed,TotalCPU,MaxRSS,MaxDiskRead,MaxDiskWrite"]),
        })
    record = {"observed_utc": datetime.now(timezone.utc).isoformat(), "jobs": observations}
    status_dir = ROOT / "pilot" / "analysis"
    status_dir.mkdir(parents=True, exist_ok=True)
    (status_dir / "latest_scheduler_status.json").write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="ascii")
    print(json.dumps(record, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
