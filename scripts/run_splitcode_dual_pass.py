#!/usr/bin/env python3
"""Run splitcode on forward and precomputed reverse-complement FASTQs.

This wrapper deliberately performs no union, duplicate reconciliation, or other
post-processing.  Its report records the two pass runtimes and their sum so the
paper can present splitcode's best-faith dual-orientation cost without claiming
native dual-orientation support.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from collections.abc import Sequence
from pathlib import Path


def splitcode_command(
    binary: Path,
    config: Path,
    threads: int,
    input_path: Path,
    output_path: Path,
    mapping_path: Path,
) -> list[str]:
    return [
        str(binary),
        "--config",
        str(config),
        "--assign",
        "--mapping",
        str(mapping_path),
        "--nFastqs",
        "1",
        "--threads",
        str(threads),
        "--output",
        str(output_path),
        str(input_path),
    ]


def run_pass(
    label: str,
    binary: Path,
    config: Path,
    threads: int,
    input_path: Path,
    output_path: Path,
    mapping_path: Path,
) -> dict[str, object]:
    command = splitcode_command(
        binary, config, threads, input_path, output_path, mapping_path
    )
    started = time.monotonic()
    completed = subprocess.run(command, check=False)
    elapsed = time.monotonic() - started
    if completed.returncode != 0:
        raise RuntimeError(
            f"splitcode {label} pass failed with exit code {completed.returncode}"
        )
    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise RuntimeError(f"splitcode {label} pass did not produce {output_path}")
    return {
        "label": label,
        "command": command,
        "wall_seconds": elapsed,
        "output": str(output_path.resolve()),
        "output_bytes": output_path.stat().st_size,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--forward-input", type=Path, required=True)
    parser.add_argument("--reverse-input", type=Path, required=True)
    parser.add_argument("--forward-output", type=Path, required=True)
    parser.add_argument("--reverse-output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args(argv)

    if args.threads <= 0:
        parser.error("--threads must be positive")
    for path in (args.binary, args.config, args.forward_input, args.reverse_input):
        if not path.is_file():
            parser.error(f"required input does not exist: {path}")
    for path in (args.forward_output, args.reverse_output, args.report):
        if path.exists():
            parser.error(f"refusing to overwrite {path}")
        path.parent.mkdir(parents=True, exist_ok=True)

    mapping_forward = args.report.with_name("splitcode-forward.mapping.txt")
    mapping_reverse = args.report.with_name("splitcode-reverse.mapping.txt")
    passes = []
    try:
        passes.append(
            run_pass(
                "forward",
                args.binary,
                args.config,
                args.threads,
                args.forward_input,
                args.forward_output,
                mapping_forward,
            )
        )
        passes.append(
            run_pass(
                "reverse-complement",
                args.binary,
                args.config,
                args.threads,
                args.reverse_input,
                args.reverse_output,
                mapping_reverse,
            )
        )
    finally:
        mapping_forward.unlink(missing_ok=True)
        mapping_reverse.unlink(missing_ok=True)

    payload = {
        "schema_version": "1.0.0",
        "measurement": "two sequential splitcode passes",
        "threads": args.threads,
        "passes": passes,
        "summed_pass_wall_seconds": sum(float(item["wall_seconds"]) for item in passes),
        "reverse_complement_precomputed_outside_measurement": True,
        "duplicate_reconciliation_performed": False,
        "duplicate_reconciliation_included_in_timing": False,
    }
    args.report.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
