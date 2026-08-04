#!/usr/bin/env python3
"""Build all isolated RaptorP LC targets and record reproducibility metadata."""

from __future__ import annotations

import json
import os
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source"
TARGETS = ("RAPTOR_LC4_IMAGE", "RAPTOR_LC4_SPEC", "RAPTOR_SED50_IMAGE", "RAPTOR_SED50_SPEC")


def main() -> None:
    environment = os.environ.copy()
    deps = Path("/home/luoji/qposim/tools/raptorp-deps")
    environment["PATH"] = f"{deps / 'bin'}:{environment['PATH']}"
    environment["LD_LIBRARY_PATH"] = (
        f"{deps / 'prefix/usr/lib/x86_64-linux-gnu'}:"
        f"{environment.get('LD_LIBRARY_PATH', '')}"
    )
    results: dict[str, object] = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "compiler": subprocess.check_output(["gcc", "--version"], text=True, env=environment).splitlines()[0],
        "h5cc_wrapper_sha256": __import__("hashlib").sha256(
            (deps / "bin" / "h5cc").read_bytes()).hexdigest(),
        "platform": platform.platform(),
        "targets": {},
    }
    for target in TARGETS:
        completed = subprocess.run(
            ["make", "-C", str(SOURCE), target, "-j8"],
            text=True, capture_output=True, env=environment, check=False)
        log_path = ROOT / "tests" / f"build_{target}.log"
        log_path.write_text(completed.stdout + completed.stderr, encoding="utf-8")
        binary = SOURCE / "bin" / target
        results["targets"][target] = {
            "command": ["make", "-C", str(SOURCE), target, "-j8"],
            "returncode": completed.returncode,
            "binary": str(binary),
            "binary_exists": binary.is_file(),
        }
        if completed.returncode:
            raise SystemExit(f"{target} build failed; see {log_path}")
    (ROOT / "provenance" / "build_metadata.json").write_text(
        json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="ascii")


if __name__ == "__main__":
    main()
