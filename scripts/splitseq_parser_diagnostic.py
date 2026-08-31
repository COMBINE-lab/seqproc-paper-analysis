#!/usr/bin/env python3
"""Decompose needletail/paraseq performance on paired SPLiT-seq input."""

from __future__ import annotations

import argparse
import json
import os
import random
import statistics
import subprocess
import time
from pathlib import Path

from real_seqproc_io_benchmark import (
    CONFIGS,
    filesystem_metadata,
    parse_time_report,
    sha256_file,
    stage_gzip,
)


ROOT = Path(__file__).resolve().parents[1]
MODES = (
    ("identity_plain", "identity", False, False),
    ("identity_gzip", "identity", True, False),
    ("identity_gzip_rapid", "identity", True, True),
    ("full_plain", "full", False, False),
    ("full_gzip", "full", True, False),
    ("full_gzip_rapid", "full", True, True),
)


def labeled_path(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("expected LABEL=PATH")
    label, raw = value.split("=", 1)
    return label, Path(raw)


def summarize(runs: list[dict[str, object]]) -> list[dict[str, object]]:
    groups: dict[tuple[str, str], list[dict[str, object]]] = {}
    for run in runs:
        groups.setdefault((str(run["parser"]), str(run["mode"])), []).append(run)
    result = []
    for (parser, mode), values in sorted(groups.items()):
        result.append(
            {
                "parser": parser,
                "mode": mode,
                "replicates": len(values),
                "median_wall_seconds": statistics.median(
                    float(value["wall_seconds"]) for value in values
                ),
                "median_fragments_per_second": statistics.median(
                    float(value["fragments_per_second"]) for value in values
                ),
                "median_peak_rss_kib": statistics.median(
                    float(value["peak_rss_kib"]) for value in values
                ),
                "median_user_seconds": statistics.median(
                    float(value["user_seconds"]) for value in values
                ),
                "median_system_seconds": statistics.median(
                    float(value["system_seconds"]) for value in values
                ),
            }
        )
    return result


def comparisons(summary: list[dict[str, object]]) -> list[dict[str, object]]:
    indexed = {(str(row["parser"]), str(row["mode"])): row for row in summary}
    result = []
    for mode, _, _, _ in MODES:
        needletail = indexed[("needletail", mode)]
        paraseq = indexed[("paraseq", mode)]
        base = float(needletail["median_wall_seconds"])
        candidate = float(paraseq["median_wall_seconds"])
        result.append(
            {
                "mode": mode,
                "paraseq_speedup": base / candidate,
                "paraseq_wall_change_percent": 100.0 * (candidate / base - 1.0),
                "paraseq_rss_change_percent": 100.0
                * (
                    float(paraseq["median_peak_rss_kib"])
                    / float(needletail["median_peak_rss_kib"])
                    - 1.0
                ),
                "needletail_wall_seconds": base,
                "paraseq_wall_seconds": candidate,
            }
        )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", action="append", type=labeled_path, required=True)
    parser.add_argument("--r1", type=Path, required=True)
    parser.add_argument("--r2", type=Path, required=True)
    parser.add_argument("--identity-geometry", type=Path, required=True)
    parser.add_argument(
        "--full-geometry", type=Path, default=CONFIGS / "splitseq_filter_edit.geom"
    )
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--reads", type=int, default=1_000_000)
    parser.add_argument("--replicates", type=int, default=9)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--affinity", default="0-7")
    args = parser.parse_args()

    binaries = {label: path.resolve() for label, path in args.binary}
    if set(binaries) != {"needletail", "paraseq"}:
        parser.error("supply exactly needletail=PATH and paraseq=PATH")
    for label, path in binaries.items():
        if not path.is_file() or not os.access(path, os.X_OK):
            parser.error(f"binary is not executable: {label}={path}")
    if min(args.reads, args.replicates, args.threads) <= 0:
        parser.error("reads, replicates, and threads must be positive")

    args.workdir.mkdir(parents=True, exist_ok=True)
    args.outdir.mkdir(parents=True, exist_ok=True)
    gzip_r1 = args.workdir / f"splitseq-{args.reads}-r1.fastq.gz"
    gzip_r2 = args.workdir / f"splitseq-{args.reads}-r2.fastq.gz"
    if not gzip_r1.exists():
        stage_gzip(args.r1, gzip_r1)
    if not gzip_r2.exists():
        stage_gzip(args.r2, gzip_r2)

    geometry = {
        "identity": args.identity_geometry.resolve(),
        "full": args.full_geometry.resolve(),
    }
    schedule = [
        (label, mode, geometry_kind, gzip_input, rapid, replicate)
        for label in binaries
        for mode, geometry_kind, gzip_input, rapid in MODES
        for replicate in range(1, args.replicates + 1)
    ]
    random.Random(0xFA57_2026).shuffle(schedule)
    runs = []
    for index, (label, mode, geometry_kind, gzip_input, rapid, replicate) in enumerate(
        schedule, 1
    ):
        input1, input2 = (gzip_r1, gzip_r2) if gzip_input else (args.r1, args.r2)
        time_report = args.workdir / f"{label}-{mode}-{replicate}.time"
        command = ["/usr/bin/time", "-v", "-o", str(time_report)]
        if args.affinity:
            command.extend(("taskset", "-c", args.affinity))
        command.extend(
            (
                str(binaries[label]),
                "run",
                "--geom",
                str(geometry[geometry_kind]),
                "--file1",
                str(input1),
                "--file2",
                str(input2),
                "--out1",
                "/dev/null",
                "--out2",
                "/dev/null",
                "--threads",
                str(args.threads),
                "--preserve-order",
            )
        )
        if rapid:
            command.extend(
                (
                    "--accelerated-gzip-input",
                    "--gzip-input-threads",
                    "1",
                    "--gzip-input-chunk-size",
                    str(256 * 1024),
                )
            )
        print(f"[{index}/{len(schedule)}] {label} {mode} rep={replicate}", flush=True)
        started = time.perf_counter()
        completed = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
        elapsed = time.perf_counter() - started
        if completed.returncode != 0:
            raise RuntimeError(f"run failed: {command}\n{completed.stderr}")
        runs.append(
            {
                "parser": label,
                "mode": mode,
                "replicate": replicate,
                "wall_seconds": elapsed,
                "fragments_per_second": args.reads / elapsed,
                "command": command,
                **parse_time_report(time_report),
            }
        )
        time_report.unlink()

    summary = summarize(runs)
    artifact = {
        "schema_version": "1.0.0",
        "scope": "paired SPLiT-seq parser/decompressor decomposition",
        "binaries": {
            label: {"path": str(path), "sha256": sha256_file(path)}
            for label, path in binaries.items()
        },
        "inputs": {
            "r1": {"path": str(args.r1.resolve()), "sha256": sha256_file(args.r1)},
            "r2": {"path": str(args.r2.resolve()), "sha256": sha256_file(args.r2)},
            "gzip_r1": {"path": str(gzip_r1), "sha256": sha256_file(gzip_r1)},
            "gzip_r2": {"path": str(gzip_r2), "sha256": sha256_file(gzip_r2)},
        },
        "geometries": {
            label: {"path": str(path), "sha256": sha256_file(path)}
            for label, path in geometry.items()
        },
        "host": {
            "hostname": os.uname().nodename,
            "machine": os.uname().machine,
            "cpu_affinity": args.affinity,
        },
        "filesystem": filesystem_metadata(args.workdir),
        "parameters": {
            "reads": args.reads,
            "replicates": args.replicates,
            "threads": args.threads,
            "preserve_order": True,
            "randomized_schedule": True,
        },
        "runs": runs,
        "summary": summary,
        "comparisons": comparisons(summary),
    }
    destination = args.outdir / "splitseq-parser-diagnostic-results.json"
    destination.write_text(json.dumps(artifact, indent=2) + "\n")
    print(destination)


if __name__ == "__main__":
    main()
