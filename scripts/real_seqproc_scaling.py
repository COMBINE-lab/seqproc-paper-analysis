#!/usr/bin/env python3
"""Run a randomized seqproc scaling matrix on the staged paper FASTQs.

This matrix measures parsing, transformations, FASTQ serialization, and writes
to `/dev/null`; it is a compute/I/O-read benchmark, not a storage-output test.
One statistics-enabled validation run per dataset records accepted/rejected
counts and match-distance summaries. Performance runs leave statistics off.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import statistics
import subprocess
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONFIGS = ROOT / "configs" / "seqproc"


DATASETS = (
    {
        "name": "splitseq_pe",
        "geometry": "publication_splitseq_pe.geom",
        "r1": "SRR6750041_10M_R1.fastq",
        "r2": "SRR6750041_10M_R2.fastq",
        "input_reads": 10_000_000,
        "additional": (),
    },
    {
        "name": "lr_splitseq",
        "geometry": "publication_lr_splitseq_dual_core.geom",
        "r1": "SRR13948564_1M.fastq",
        "r2": None,
        "input_reads": 1_000_000,
        "additional": (),
    },
    {
        "name": "10x_short",
        "geometry": "10x_v2.geom",
        "r1": "10x_short/SRR8315379_1M_R1.fastq",
        "r2": "10x_short/SRR8315379_1M_R2.fastq",
        "input_reads": 1_000_000,
        "additional": (),
    },
    {
        "name": "sciseq",
        "geometry": "sciseq3_edit.geom",
        "r1": "SRR7827254_1M_1.fastq",
        "r2": "SRR7827254_1M_2.fastq",
        "input_reads": 1_000_000,
        "additional": (),
    },
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def parse_time_report(path: Path) -> dict[str, float]:
    values: dict[str, float] = {}
    mapping = {
        "User time (seconds)": "user_seconds",
        "System time (seconds)": "system_seconds",
        "Maximum resident set size (kbytes)": "peak_rss_kib",
        "File system inputs": "filesystem_inputs",
        "File system outputs": "filesystem_outputs",
    }
    for line in path.read_text().splitlines():
        stripped = line.strip()
        for prefix, key in mapping.items():
            if stripped.startswith(prefix + ":"):
                values[key] = float(stripped.rsplit(":", 1)[1].strip())
    return values


def command_for(
    binary: Path,
    data_dir: Path,
    dataset: dict[str, object],
    threads: int,
    mode: str,
    affinity: str | None,
    time_report: Path,
    summary: Path | None = None,
) -> list[str]:
    command = ["/usr/bin/time", "-v", "-o", str(time_report)]
    if affinity:
        command.extend(("taskset", "-c", affinity))
    command.extend(
        (
            str(binary),
            "run",
            "--geom",
            str(CONFIGS / str(dataset["geometry"])),
            "--file1",
            str(data_dir / str(dataset["r1"])),
            "--out1",
            "/dev/null",
            "--threads",
            str(threads),
        )
    )
    if dataset["r2"] is not None:
        command.extend(
            (
                "--file2",
                str(data_dir / str(dataset["r2"])),
                "--out2",
                "/dev/null",
            )
        )
    if mode == "staged":
        command.append("--staged-pipeline")
    elif mode == "ordered":
        command.append("--preserve-order")
    elif mode != "legacy":
        raise ValueError(f"unknown mode: {mode}")
    for additional in dataset["additional"]:
        command.extend(("--additional", str(CONFIGS / str(additional))))
    if summary is not None:
        command.extend(("--summary", str(summary)))
    return command


def execute(command: list[str], time_report: Path) -> tuple[float, dict[str, float]]:
    started = time.perf_counter()
    completed = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    elapsed = time.perf_counter() - started
    if completed.returncode != 0:
        raise RuntimeError(
            f"command failed ({completed.returncode}): {command}\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )
    timing = parse_time_report(time_report)
    time_report.unlink()
    return elapsed, timing


def summarize(runs: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, int], list[dict[str, object]]] = {}
    for run in runs:
        key = (str(run["dataset"]), str(run["mode"]), int(run["threads"]))
        grouped.setdefault(key, []).append(run)
    result = []
    for (dataset, mode, threads), values in sorted(grouped.items()):
        rates = [float(value["input_reads_per_second"]) for value in values]
        result.append(
            {
                "dataset": dataset,
                "mode": mode,
                "threads": threads,
                "replicates": len(values),
                "median_input_reads_per_second": statistics.median(rates),
                "min_input_reads_per_second": min(rates),
                "max_input_reads_per_second": max(rates),
                "median_wall_seconds": statistics.median(
                    float(value["wall_seconds"]) for value in values
                ),
                "median_peak_rss_kib": statistics.median(
                    float(value["peak_rss_kib"]) for value in values
                ),
            }
        )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 2, 4, 8])
    parser.add_argument("--affinity", default="0-7", help="taskset CPU list; empty disables")
    parser.add_argument("--only", action="append", help="dataset name; repeat to select several")
    parser.add_argument(
        "--validation-only",
        action="store_true",
        help="record statistics-enabled validation reports without performance replicates",
    )
    args = parser.parse_args()
    if args.replicates <= 0 or any(value <= 0 for value in args.threads):
        parser.error("replicates and thread counts must be positive")
    binary = args.binary.resolve()
    data_dir = args.data_dir.resolve()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    selected = set(args.only or [])
    known = {str(dataset["name"]) for dataset in DATASETS}
    if selected - known:
        parser.error(f"unknown datasets: {', '.join(sorted(selected - known))}")
    datasets = [
        dataset for dataset in DATASETS if not selected or str(dataset["name"]) in selected
    ]
    for dataset in datasets:
        for key in ("r1", "r2"):
            relative = dataset[key]
            if relative is not None and not (data_dir / str(relative)).is_file():
                parser.error(f"missing input: {data_dir / str(relative)}")

    provenance = json.loads((data_dir / "paper_data_provenance.json").read_text())
    input_records = {
        item["output"]: item for item in provenance["outputs"]
    }
    top_threads = max(args.threads)
    validations = []
    for dataset in datasets:
        label = str(dataset["name"])
        report_path = outdir / f"{label}.validation.json"
        time_path = outdir / f"{label}.validation.time"
        command = command_for(
            binary,
            data_dir,
            dataset,
            top_threads,
            "staged",
            args.affinity or None,
            time_path,
            report_path,
        )
        print(f"[validation] {label} t={top_threads}", flush=True)
        elapsed, timing = execute(command, time_path)
        report = json.loads(report_path.read_text())
        report_path.unlink()
        validations.append(
            {
                "dataset": label,
                "wall_seconds": elapsed,
                **timing,
                "report": report,
            }
        )

    modes = [("legacy", 1)]
    modes.extend(("staged", threads) for threads in sorted(set(args.threads)))
    modes.append(("ordered", top_threads))
    schedule = [] if args.validation_only else [
        (dataset, mode, threads, replicate)
        for replicate in range(1, args.replicates + 1)
        for dataset in datasets
        for mode, threads in modes
    ]
    random.Random("seqproc-real-scaling-v1").shuffle(schedule)
    runs = []
    for index, (dataset, mode, threads, replicate) in enumerate(schedule, 1):
        label = str(dataset["name"])
        print(
            f"[{index}/{len(schedule)}] {label} {mode} t={threads} replicate={replicate}",
            flush=True,
        )
        time_path = outdir / f"{label}-{mode}-t{threads}-r{replicate}.time"
        command = command_for(
            binary,
            data_dir,
            dataset,
            threads,
            mode,
            args.affinity or None,
            time_path,
        )
        elapsed, timing = execute(command, time_path)
        input_reads = int(dataset["input_reads"])
        runs.append(
            {
                "dataset": label,
                "mode": mode,
                "threads": threads,
                "replicate": replicate,
                "wall_seconds": elapsed,
                "input_reads": input_reads,
                "input_reads_per_second": input_reads / elapsed,
                **timing,
            }
        )

    inputs = {}
    for dataset in datasets:
        files = []
        for key in ("r1", "r2"):
            relative = dataset[key]
            if relative is None:
                continue
            record = input_records[str(relative)]
            files.append(
                {
                    "path": str((data_dir / str(relative)).resolve()),
                    "bytes": record["bytes"],
                    "sha256": record["sha256"],
                    "records": record["records"],
                }
            )
        inputs[str(dataset["name"])] = {
            "files": files,
            "geometry": str((CONFIGS / str(dataset["geometry"])).resolve()),
            "geometry_sha256": sha256_file(CONFIGS / str(dataset["geometry"])),
        }
    artifact = {
        "schema_version": "1.0.0",
        "scope": "real paper FASTQ compute/read-I/O scaling; output sink is /dev/null",
        "binary": {"path": str(binary), "sha256": sha256_file(binary)},
        "host": {
            "hostname": os.uname().nodename,
            "machine": os.uname().machine,
            "cpu_affinity": args.affinity or None,
        },
        "parameters": {
            "replicates": args.replicates,
            "thread_counts": sorted(set(args.threads)),
            "randomized_schedule": True,
            "statistics_in_performance_runs": False,
            "output_sink": "/dev/null",
        },
        "inputs": inputs,
        "validations": validations,
        "runs": runs,
        "summary": summarize(runs),
    }
    artifact_path = outdir / (
        "real-seqproc-validation-results.json"
        if args.validation_only
        else "real-seqproc-scaling-results.json"
    )
    artifact_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {artifact_path}")


if __name__ == "__main__":
    main()
