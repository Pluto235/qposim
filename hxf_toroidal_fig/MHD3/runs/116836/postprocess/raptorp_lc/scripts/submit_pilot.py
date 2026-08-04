#!/usr/bin/env python3
"""Submit a bounded Stage 1 pilot array and append its Slurm provenance."""

from __future__ import annotations

import argparse
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", type=Path, required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--max-concurrent", type=int, default=2)
    args = parser.parse_args()
    cases_path = args.cases.resolve()
    cases = json.loads(cases_path.read_text(encoding="ascii"))
    if not cases:
        raise SystemExit("cannot submit an empty case list")
    case_file = ROOT / "pilot" / f"{args.label}_case_ids.txt"
    case_file.write_text("".join(f"{case['case_id']}\n" for case in cases), encoding="ascii")
    command = [
        "sbatch", f"--array=0-{len(cases) - 1}%{args.max_concurrent}",
        f"--export=ALL,CASE_FILE={cases_path},CASE_IDS_FILE={case_file}",
        str(ROOT / "pilot" / "run_pilot.slurm"),
    ]
    output = subprocess.check_output(command, text=True).strip()
    job_id = output.split()[-1]
    registry_path = ROOT / "pilot" / "job_registry.json"
    registry = json.loads(registry_path.read_text()) if registry_path.exists() else {"jobs": []}
    registry["jobs"].append({
        "submitted_utc": datetime.now(timezone.utc).isoformat(),
        "label": args.label,
        "job_id": job_id,
        "array_size": len(cases),
        "max_concurrent": args.max_concurrent,
        "case_file": str(cases_path),
        "case_ids_file": str(case_file),
        "sbatch_command": command,
        "sbatch_output": output,
    })
    registry_path.write_text(json.dumps(registry, indent=2, sort_keys=True) + "\n", encoding="ascii")
    print(job_id)


if __name__ == "__main__":
    main()
