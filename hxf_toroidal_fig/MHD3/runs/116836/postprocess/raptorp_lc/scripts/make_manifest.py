#!/usr/bin/env python3
"""Create and validate the read-only 351-epoch MHD3 snapshot manifest."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT.parents[1]
EXPECTED_SIZE = 3_317_760_000
PRIMITIVES = "rho vx1 vx2 vx3 Bx1 Bx2 Bx3 prs tr1"


def parse_dbl_out(path: Path) -> dict[int, dict[str, object]]:
    records: dict[int, dict[str, object]] = {}
    for line in path.read_text(encoding="ascii").splitlines():
        fields = line.split()
        if len(fields) < 10:
            raise ValueError(f"invalid dbl.out record: {line!r}")
        index = int(fields[0])
        records[index] = {
            "actual_time_Rj_c": float(fields[1]),
            "file_mode": fields[4],
            "byte_order": fields[5],
            "primitive_order": " ".join(fields[6:]),
        }
    return records


def main() -> None:
    records = parse_dbl_out(RUN / "dbl.out")
    selected = list(range(30, 731, 2))
    rows: list[dict[str, object]] = []
    failures: list[str] = []
    modes: set[str] = set()
    orders: set[str] = set()
    byte_orders: set[str] = set()
    for ordinal, index in enumerate(selected):
        record = records.get(index)
        if record is None:
            failures.append(f"missing dbl.out index {index}")
            continue
        path = RUN / f"data.{index:04d}.dbl"
        target_time = 300.0 + 2.0 * ordinal
        try:
            resolved = path.resolve(strict=True)
            size = resolved.stat().st_size
            resolves = True
        except FileNotFoundError:
            resolved = path
            size = -1
            resolves = False
        time_error = float(record["actual_time_Rj_c"]) - target_time
        status = "valid"
        if not resolves:
            status = "broken_symlink"
        elif size != EXPECTED_SIZE:
            status = "wrong_size"
        elif abs(time_error) >= 0.01:
            status = "time_error"
        if status != "valid":
            failures.append(f"index {index}: {status}")
        modes.add(str(record["file_mode"]))
        orders.add(str(record["primitive_order"]))
        byte_orders.add(str(record["byte_order"]))
        rows.append({
            "epoch_ordinal": ordinal,
            "snapshot_index": index,
            "actual_time_Rj_c": f"{record['actual_time_Rj_c']:.9f}",
            "target_time_Rj_c": f"{target_time:.9f}",
            "time_error": f"{time_error:.9f}",
            "snapshot_path": str(path),
            "resolved_path": str(resolved),
            "is_symlink": path.is_symlink(),
            "size_bytes": size,
            "primitive_order": record["primitive_order"],
            "validation_status": status,
        })
    if len(rows) != 351 or modes != {"single_file"} or orders != {PRIMITIVES} or byte_orders != {"little"}:
        failures.append("dbl.out metadata consistency failure")
    fieldnames = list(rows[0]) if rows else []
    csv_path = ROOT / "manifests" / "snapshot_manifest.csv"
    with csv_path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    audit = {
        "expected_records": 351,
        "actual_records": len(rows),
        "first_snapshot_index": selected[0],
        "last_snapshot_index": selected[-1],
        "expected_size_bytes": EXPECTED_SIZE,
        "file_modes": sorted(modes),
        "byte_orders": sorted(byte_orders),
        "primitive_orders": sorted(orders),
        "failures": failures,
        "passed": not failures,
    }
    (ROOT / "manifests" / "snapshot_manifest_audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if failures:
        raise SystemExit("manifest validation failed: " + "; ".join(failures))


if __name__ == "__main__":
    main()
