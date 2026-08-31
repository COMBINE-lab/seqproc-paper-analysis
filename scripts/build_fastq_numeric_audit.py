#!/usr/bin/env python3
"""Build the dependency-free publication FASTQ auditor and record a receipt."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path

from benchmark_harness import sha256_file


ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "tools" / "bin" / "fastq-numeric-audit",
    )
    args = parser.parse_args()
    source = ROOT / "tools" / "fastq_numeric_audit.rs"
    rustc = shutil.which("rustc")
    if rustc is None:
        parser.error("rustc is not available; load an appropriate compiler module")
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    command = [
        rustc,
        "-C",
        "opt-level=3",
        "-C",
        "target-cpu=x86-64-v3",
        "-o",
        str(output),
        str(source),
    ]
    subprocess.run(command, check=True)
    receipt = {
        "schema_version": "1.0.0",
        "command": command,
        "rustc_version": subprocess.run(
            [rustc, "--version"], check=True, capture_output=True, text=True
        ).stdout.strip(),
        "source": str(source),
        "source_sha256": sha256_file(source),
        "output": str(output),
        "output_sha256": sha256_file(output),
    }
    receipt_path = output.with_suffix(".build.json")
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
