#!/usr/bin/env python3
"""Run one pilot/production epoch transactionally with strict output checks."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT.parents[1]
SOURCE = ROOT / "source"


def load_case(path: Path, case_id: str) -> dict[str, Any]:
    for case in json.loads(path.read_text(encoding="ascii")):
        if case["case_id"] == case_id:
            return case
    raise ValueError(f"unknown case_id {case_id}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def model_input(case: dict[str, Any]) -> str:
    physics = json.loads((ROOT / "config" / "physics.json").read_text(encoding="ascii"))
    adopted = physics["adopted_reproduction_assumptions"]
    width, height = case["resolution"]
    fov_x, fov_y = case["fov"]
    return "\n".join([
        f"MBH (Msun) {adopted['effective_raptor_code_mass_msun']:.9g}",
        f"DISTANCE (kpc) {adopted['distance_kpc']:.9g}",
        "M_UNIT (g) 1.e26",
        f"R_LOW (-) {adopted['r_low']:.9g}",
        f"R_HIGH (-) {adopted['r_high']:.9g}",
        "INCLINATION (deg) 10",
        "azmuth (deg) 0",
        f"IMG_WIDTH (pixels) {width}",
        f"IMG_HEIGHT (pixels) {height}",
        "SHIFT (Rg) 0",
        f"CAM_SIZE_X (Rg) {fov_x:.9g}",
        f"CAM_SIZE_Y (Rg) {fov_y:.9g}",
        "FREQS_PER_DEC (-) 1",
        "FREQ (Hz) 1.5e9",
        "STEPSIZE (-) 0.02",
        "STEPSIZE_MAX (-) 0.02",
        "STEPSIZE_MIN (-) 0.01",
        "",
    ])


def validate_spectrum(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for line in path.read_text(encoding="ascii").splitlines():
        if not line or line.startswith("#"):
            continue
        values = [float(value) for value in line.split()]
        if len(values) < 5:
            raise ValueError(f"malformed spectrum line: {line!r}")
        frequency, intensity, q, u, v = values[:5]
        if not all(math.isfinite(value) for value in values):
            raise ValueError("non-finite spectrum value")
        if intensity <= 0:
            raise ValueError("non-positive integrated Stokes I")
        if math.hypot(q, u) > intensity * (1.0 + 1e-8) or abs(v) > intensity * (1.0 + 1e-8):
            raise ValueError("unphysical integrated Stokes values")
        rows.append({"frequency_hz": frequency, "I_Jy": intensity, "Q_Jy": q, "U_Jy": u, "V_Jy": v})
    if len(rows) != 4:
        raise ValueError(f"expected 4 frequency rows, found {len(rows)}")
    return rows


def validate_image(path: Path, rows: list[dict[str, float]]) -> None:
    with h5py.File(path, "r") as handle:
        for row in rows:
            suffix = f"{row['frequency_hz']:.6e}"
            for stokes in "IQUV":
                dataset = handle[f"{stokes}{suffix}"]
                values = dataset[()]
                if not np.isfinite(values).all():
                    raise ValueError(f"non-finite {dataset.name}")
            image_i = handle[f"I{suffix}"][()]
            edge = np.concatenate((image_i[0], image_i[-1], image_i[:, 0], image_i[:, -1]))
            if not np.isfinite(edge).all():
                raise ValueError("non-finite camera-boundary intensity")


def complete_is_valid(directory: Path) -> bool:
    try:
        complete = json.loads((directory / "COMPLETE.json").read_text(encoding="ascii"))
        return complete["status"] == "complete" and (directory / "spectrum.dat").is_file()
    except (FileNotFoundError, json.JSONDecodeError, KeyError):
        return False


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", type=Path, required=True)
    parser.add_argument("--case-id", required=True)
    parser.add_argument("--kind", choices=("pilot", "production"), required=True)
    args = parser.parse_args()
    case = load_case(args.cases, args.case_id)
    if args.kind == "production" and case["volume_mask"] != "kink-only":
        raise SystemExit("production refuses any volume mode other than explicit kink-only")
    base = ROOT / args.kind / ("raw" if args.kind == "pilot" else "epochs")
    base.mkdir(parents=True, exist_ok=True)
    grid_link = base / "grid.out"
    if not grid_link.exists() and not grid_link.is_symlink():
        grid_link.symlink_to(RUN / "grid.out")
    final = base / case["case_id"]
    if complete_is_valid(final):
        print(f"already complete: {final}")
        return
    if final.exists():
        failed = ROOT / args.kind / "failed"
        failed.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        os.replace(final, failed / f"{case['case_id']}_{stamp}")
    temporary = base / f".{case['case_id']}.tmp-{uuid.uuid4().hex}"
    temporary.mkdir()
    try:
        (temporary / "model.in").write_text(model_input(case), encoding="ascii")
        shutil.copyfile(ROOT / "config" / "frequencies_lc4.txt", temporary / "frequencies.txt")
        (temporary / "grid.out").symlink_to(RUN / "grid.out")
        snapshot = RUN / f"data.{case['snapshot_index']:04d}.dbl"
        binary = SOURCE / "bin" / case["target"]
        if not binary.is_file():
            raise FileNotFoundError(f"missing built target {binary}")
        env = os.environ.copy()
        env["RAPTORP_VOLUME_MASK"] = case["volume_mask"]
        env["RAPTORP_TRACER_MASK"] = "1" if case["tracer_mask"] == "on" else "0"
        env.setdefault("OMP_NUM_THREADS", env.get("SLURM_CPUS_PER_TASK", "32"))
        command = [str(binary), str(temporary / "model.in"), str(snapshot), str(case["time"])]
        completed = subprocess.run(
            ["/usr/bin/time", "-v", "-o", str(temporary / "resource.txt"), *command],
            cwd=temporary, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            check=False)
        (temporary / "stdout_stderr.log").write_text(completed.stdout, encoding="utf-8")
        if completed.returncode:
            raise RuntimeError(f"RaptorP exited {completed.returncode}")
        output = temporary / "output"
        spectrum = next(output.glob("spectrum_*.dat"), None)
        if spectrum is None:
            raise FileNotFoundError("RaptorP did not write spectrum")
        shutil.copyfile(spectrum, temporary / "spectrum.dat")
        rows = validate_spectrum(temporary / "spectrum.dat")
        image = next(output.glob("img_data_*.h5"), None)
        if case["target"].endswith("IMAGE"):
            if image is None:
                raise FileNotFoundError("image target did not write HDF5 output")
            validate_image(image, rows)
        metadata = {
            "case": case,
            "command": command,
            "completed_utc": datetime.now(timezone.utc).isoformat(),
            "executable_sha256": sha256(binary),
            "snapshot_resolved_path": str(snapshot.resolve(strict=True)),
            "spectrum": rows,
        }
        (temporary / "run_metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="ascii")
        (temporary / "executable.sha256").write_text(sha256(binary) + "\n", encoding="ascii")
        (temporary / "COMPLETE.json").write_text(json.dumps({"status": "complete", "case_id": case["case_id"]}, indent=2) + "\n", encoding="ascii")
        os.replace(temporary, final)
    except BaseException:
        if temporary.exists():
            failed = ROOT / args.kind / "failed"
            failed.mkdir(parents=True, exist_ok=True)
            stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
            os.replace(temporary, failed / f"{case['case_id']}_{stamp}")
        raise


if __name__ == "__main__":
    main()
