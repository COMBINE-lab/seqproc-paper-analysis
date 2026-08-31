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
    x_only: bool = False,
) -> list[str]:
    command = [
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
    ]
    if x_only:
        command.append("--x-only")
    else:
        command.extend(["--output", str(output_path)])
    command.append(str(input_path))
    return command


def run_pass(
    label: str,
    binary: Path,
    config: Path,
    threads: int,
    input_path: Path,
    output_path: Path,
    mapping_path: Path,
    discard_output: bool = False,
    working_directory: Path | None = None,
    x_only: bool = False,
) -> dict[str, object]:
    command = splitcode_command(
        binary, config, threads, input_path, output_path, mapping_path, x_only
    )
    started = time.monotonic()
    completed = subprocess.run(command, check=False, cwd=working_directory)
    elapsed = time.monotonic() - started
    if completed.returncode != 0:
        raise RuntimeError(
            f"splitcode {label} pass failed with exit code {completed.returncode}"
        )
    if not discard_output and not x_only and (
        not output_path.is_file() or output_path.stat().st_size == 0
    ):
        raise RuntimeError(f"splitcode {label} pass did not produce {output_path}")
    return {
        "label": label,
        "command": command,
        "wall_seconds": elapsed,
        "output": None if x_only else str(output_path.resolve()),
        "output_bytes": (
            None if discard_output or x_only else output_path.stat().st_size
        ),
        "output_discarded": discard_output or x_only,
        "working_directory": (
            str(working_directory.resolve()) if working_directory is not None else None
        ),
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
    parser.add_argument(
        "--discard-output",
        action="store_true",
        help="require both sequence outputs to be /dev/null and skip file validation",
    )
    parser.add_argument(
        "--mapping-sink",
        type=Path,
        help="write both per-pass mappings to this sink (timing runs use /dev/null)",
    )
    parser.add_argument(
        "--extraction-output",
        action="append",
        default=[],
        metavar="NAME.fastq",
        help="config-derived extraction filename to isolate for each pass",
    )
    parser.add_argument(
        "--x-only",
        action="store_true",
        help="emit only config-derived extraction FASTQs, not retained full reads",
    )
    args = parser.parse_args(argv)

    if args.threads <= 0:
        parser.error("--threads must be positive")
    for path in (args.binary, args.config, args.forward_input, args.reverse_input):
        if not path.is_file():
            parser.error(f"required input does not exist: {path}")
    sequence_paths = () if args.x_only else (args.forward_output, args.reverse_output)
    for path in (*sequence_paths, args.report):
        if path.exists() and path != Path("/dev/null"):
            parser.error(f"refusing to overwrite {path}")
        if path != Path("/dev/null"):
            path.parent.mkdir(parents=True, exist_ok=True)

    if args.discard_output and (
        args.forward_output != Path("/dev/null")
        or args.reverse_output != Path("/dev/null")
    ):
        parser.error("--discard-output requires both outputs to be /dev/null")
    for name in args.extraction_output:
        if Path(name).name != name or name in ("", ".", ".."):
            parser.error(f"invalid extraction output name: {name!r}")

    mapping_forward = args.mapping_sink or args.report.with_name(
        "splitcode-forward.mapping.txt"
    )
    mapping_reverse = args.mapping_sink or args.report.with_name(
        "splitcode-reverse.mapping.txt"
    )
    passes = []
    try:
        forward_work = args.report.with_name("splitcode-forward-work")
        reverse_work = args.report.with_name("splitcode-reverse-work")
        for work in (forward_work, reverse_work):
            work.mkdir(parents=True, exist_ok=False)
            if args.discard_output:
                for name in args.extraction_output:
                    (work / name).symlink_to("/dev/null")
        passes.append(
            run_pass(
                "forward",
                args.binary,
                args.config,
                args.threads,
                args.forward_input,
                args.forward_output,
                mapping_forward,
                args.discard_output,
                forward_work,
                args.x_only,
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
                args.discard_output,
                reverse_work,
                args.x_only,
            )
        )
    finally:
        if args.mapping_sink is None:
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
        "sequence_output_sink": "/dev/null" if args.discard_output else None,
        "mapping_sink": str(args.mapping_sink) if args.mapping_sink else None,
        "x_only": args.x_only,
    }
    args.report.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
