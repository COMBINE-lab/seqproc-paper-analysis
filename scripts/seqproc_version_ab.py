#!/usr/bin/env python3
"""A/B benchmark two or more seqproc binaries on identical paper workloads.

The primary ``paper`` mode reproduces the manuscript's uncompressed-input,
uncompressed-output execution shape. ``sink`` keeps parsing, transformation,
FASTQ serialization, and all graph work, but sends output to /dev/null to
separate compute from persistent-output cost. Runs are randomized and pinned;
single-thread control runs require byte-identical outputs across versions.
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

from real_seqproc_io_benchmark import (
    CONFIGS,
    DATASETS,
    copy_fastq_prefix,
    count_fastq_records,
    filesystem_metadata,
    parse_time_report,
    sha256_file,
)


ROOT = Path(__file__).resolve().parents[1]


def labeled_path(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("expected LABEL=PATH")
    label, raw_path = value.split("=", 1)
    if not label or not raw_path:
        raise argparse.ArgumentTypeError("expected nonempty LABEL=PATH")
    return label, Path(raw_path)


def fastq_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def command_for(
    binary: Path,
    dataset: dict[str, object],
    staged: dict[str, Path | None],
    output1: Path,
    output2: Path | None,
    threads: int,
    affinity: str | None,
    time_report: Path,
) -> list[str]:
    input1 = staged["r1"]
    input2 = staged["r2"]
    assert input1 is not None
    command = ["/usr/bin/time", "-v", "-o", str(time_report)]
    if affinity:
        command.extend(("taskset", "-c", affinity))
    # Use the legacy flag-only interface because it is common to the preprint
    # binary and the current backward-compatible CLI.
    command.extend(
        (
            str(binary),
            "--geom",
            str(dataset.get("geometry_path", CONFIGS / str(dataset["geometry"]))),
            "--file1",
            str(input1),
            "--out1",
            str(output1),
            "--threads",
            str(threads),
        )
    )
    if input2 is not None and output2 is not None:
        command.extend(("--file2", str(input2), "--out2", str(output2)))
    for additional in dataset["additional"]:
        command.extend(("--additional", str(CONFIGS / str(additional))))
    return command


def summarize(runs: list[dict[str, object]]) -> list[dict[str, object]]:
    groups: dict[tuple[str, str, str], list[dict[str, object]]] = {}
    for run in runs:
        key = (str(run["version"]), str(run["dataset"]), str(run["mode"]))
        groups.setdefault(key, []).append(run)
    result = []
    for (version, dataset, mode), values in sorted(groups.items()):
        walls = [float(value["wall_seconds"]) for value in values]
        rss = [float(value["peak_rss_kib"]) for value in values]
        result.append(
            {
                "version": version,
                "dataset": dataset,
                "mode": mode,
                "replicates": len(values),
                "median_wall_seconds": statistics.median(walls),
                "min_wall_seconds": min(walls),
                "max_wall_seconds": max(walls),
                "median_input_fragments_per_second": statistics.median(
                    float(value["input_fragments_per_second"]) for value in values
                ),
                "median_peak_rss_kib": statistics.median(rss),
                "min_peak_rss_kib": min(rss),
                "max_peak_rss_kib": max(rss),
            }
        )
    return result


def comparisons(
    summary: list[dict[str, object]], baseline: str, candidate: str
) -> list[dict[str, object]]:
    indexed = {
        (str(row["version"]), str(row["dataset"]), str(row["mode"])): row
        for row in summary
    }
    result = []
    for key, base in indexed.items():
        version, dataset, mode = key
        if version != baseline:
            continue
        current = indexed[(candidate, dataset, mode)]
        base_wall = float(base["median_wall_seconds"])
        current_wall = float(current["median_wall_seconds"])
        base_rss = float(base["median_peak_rss_kib"])
        current_rss = float(current["median_peak_rss_kib"])
        result.append(
            {
                "dataset": dataset,
                "mode": mode,
                "baseline": baseline,
                "candidate": candidate,
                "speedup": base_wall / current_wall,
                "wall_time_change_percent": 100.0 * (current_wall / base_wall - 1.0),
                "rss_change_percent": 100.0 * (current_rss / base_rss - 1.0),
                "baseline_wall_seconds": base_wall,
                "candidate_wall_seconds": current_wall,
                "baseline_peak_rss_kib": base_rss,
                "candidate_peak_rss_kib": current_rss,
            }
        )
    return sorted(result, key=lambda row: (str(row["dataset"]), str(row["mode"])))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", action="append", type=labeled_path, required=True)
    parser.add_argument("--revision", action="append", type=labeled_path, default=[])
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    parser.add_argument("--input-workdir", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--reads", type=int, default=1_000_000)
    parser.add_argument("--replicates", type=int, default=5)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--affinity", default="0-7")
    parser.add_argument("--mode", action="append", choices=("paper", "sink"))
    parser.add_argument("--only", action="append")
    parser.add_argument(
        "--geometry-override",
        action="append",
        type=labeled_path,
        default=[],
        help="repeatable DATASET=PATH override, applied identically to every binary",
    )
    args = parser.parse_args()

    if min(args.reads, args.replicates, args.threads) <= 0:
        parser.error("reads, replicates, and threads must be positive")
    binaries = dict(args.binary)
    if len(binaries) != len(args.binary):
        parser.error("binary labels must be unique")
    if args.baseline not in binaries or args.candidate not in binaries:
        parser.error("baseline and candidate must name supplied binaries")
    for label, path in list(binaries.items()):
        resolved = path.resolve()
        if not resolved.is_file() or not os.access(resolved, os.X_OK):
            parser.error(f"binary is not executable: {label}={resolved}")
        binaries[label] = resolved

    selected = set(args.only or [])
    known = {str(dataset["name"]) for dataset in DATASETS}
    if selected - known:
        parser.error(f"unknown datasets: {', '.join(sorted(selected - known))}")
    geometry_overrides = {label: path.resolve() for label, path in args.geometry_override}
    if geometry_overrides.keys() - known:
        parser.error(
            "unknown geometry override datasets: "
            + ", ".join(sorted(geometry_overrides.keys() - known))
        )
    datasets = []
    for source_dataset in DATASETS:
        if selected and str(source_dataset["name"]) not in selected:
            continue
        dataset = dict(source_dataset)
        name = str(dataset["name"])
        if name in geometry_overrides:
            dataset["geometry_path"] = geometry_overrides[name]
        datasets.append(dataset)
    modes = args.mode or ["paper", "sink"]
    revisions = {label: str(path) for label, path in args.revision}

    args.input_workdir.mkdir(parents=True, exist_ok=True)
    args.output_root.mkdir(parents=True, exist_ok=True)
    args.outdir.mkdir(parents=True, exist_ok=True)

    staged_inputs: dict[str, dict[str, Path | None]] = {}
    input_manifest: dict[str, object] = {}
    for dataset in datasets:
        name = str(dataset["name"])
        staged: dict[str, Path | None] = {}
        records = []
        for lane in ("r1", "r2"):
            relative = dataset[lane]
            if relative is None:
                staged[lane] = None
                continue
            source = args.data_dir / str(relative)
            destination = args.input_workdir / f"{name}-{args.reads}-{lane}.fastq"
            if not destination.exists():
                copy_fastq_prefix(source, destination, args.reads)
            staged[lane] = destination
            records.append(
                {
                    "lane": lane,
                    "source": str(source.resolve()),
                    "source_sha256": sha256_file(source),
                    "staged": str(destination.resolve()),
                    "staged_sha256": sha256_file(destination),
                    "bytes": destination.stat().st_size,
                }
            )
        staged_inputs[name] = staged
        geometry = Path(
            dataset.get("geometry_path", CONFIGS / str(dataset["geometry"]))
        )
        input_manifest[name] = {
            "records": records,
            "geometry": str(geometry.resolve()),
            "geometry_sha256": sha256_file(geometry),
        }

    # A deterministic one-thread control establishes semantic equivalence even
    # though unordered multithreaded output is allowed in the timed paper mode.
    correctness: dict[str, dict[str, object]] = {}
    for dataset in datasets:
        name = str(dataset["name"])
        version_results = {}
        for version, binary in binaries.items():
            out1 = args.output_root / f"correctness-{version}-{name}-r1.fastq"
            out2 = (
                args.output_root / f"correctness-{version}-{name}-r2.fastq"
                if staged_inputs[name]["r2"] is not None
                else None
            )
            report = args.output_root / f"correctness-{version}-{name}.time"
            command = command_for(
                binary, dataset, staged_inputs[name], out1, out2, 1, args.affinity, report
            )
            completed = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
            if completed.returncode != 0:
                raise RuntimeError(f"correctness run failed: {command}\n{completed.stderr}")
            outputs = [out1] + ([out2] if out2 is not None else [])
            version_results[version] = {
                "records": [count_fastq_records(path, False) for path in outputs],
                "sha256": [fastq_digest(path) for path in outputs],
            }
            report.unlink()
            for output in outputs:
                output.unlink()
        reference = version_results[args.baseline]
        correctness[name] = {
            "equivalent": all(result == reference for result in version_results.values()),
            "versions": version_results,
        }

    schedule = [
        (dataset, mode, version, replicate)
        for dataset in datasets
        for mode in modes
        for version in binaries
        for replicate in range(1, args.replicates + 1)
    ]
    random.Random(0x5E9A_B2026).shuffle(schedule)
    runs = []
    for index, (dataset, mode, version, replicate) in enumerate(schedule, 1):
        name = str(dataset["name"])
        if mode == "sink":
            out1 = Path("/dev/null")
            out2 = Path("/dev/null") if staged_inputs[name]["r2"] is not None else None
        else:
            out1 = args.output_root / f"{version}-{name}-{mode}-r1.fastq"
            out2 = (
                args.output_root / f"{version}-{name}-{mode}-r2.fastq"
                if staged_inputs[name]["r2"] is not None
                else None
            )
        time_report = args.output_root / f"{version}-{name}-{mode}-{replicate}.time"
        command = command_for(
            binaries[version],
            dataset,
            staged_inputs[name],
            out1,
            out2,
            args.threads,
            args.affinity or None,
            time_report,
        )
        print(
            f"[{index}/{len(schedule)}] {version} {name} {mode} rep={replicate}",
            flush=True,
        )
        started = time.perf_counter()
        completed = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
        elapsed = time.perf_counter() - started
        if completed.returncode != 0:
            raise RuntimeError(f"run failed: {command}\n{completed.stderr}")
        output_paths = [] if mode == "sink" else [out1] + ([out2] if out2 else [])
        output_records = [count_fastq_records(path, False) for path in output_paths]
        runs.append(
            {
                "version": version,
                "dataset": name,
                "mode": mode,
                "replicate": replicate,
                "threads": args.threads,
                "input_fragments": args.reads,
                "output_fragments": output_records,
                "wall_seconds": elapsed,
                "input_fragments_per_second": args.reads / elapsed,
                "command": command,
                **parse_time_report(time_report),
            }
        )
        time_report.unlink()
        for output in output_paths:
            output.unlink()

    summary = summarize(runs)
    artifact = {
        "schema_version": "1.0.0",
        "scope": "preprint-versus-current seqproc A/B",
        "binaries": {
            label: {
                "path": str(path),
                "sha256": sha256_file(path),
                "revision": revisions.get(label),
            }
            for label, path in binaries.items()
        },
        "host": {
            "hostname": os.uname().nodename,
            "machine": os.uname().machine,
            "cpu_affinity": args.affinity or None,
        },
        "filesystem": filesystem_metadata(args.output_root),
        "parameters": {
            "reads": args.reads,
            "replicates": args.replicates,
            "threads": args.threads,
            "modes": modes,
            "randomized_schedule": True,
        },
        "inputs": input_manifest,
        "correctness": correctness,
        "runs": runs,
        "summary": summary,
        "comparisons": comparisons(summary, args.baseline, args.candidate),
    }
    destination = args.outdir / "seqproc-version-ab-results.json"
    destination.write_text(json.dumps(artifact, indent=2) + "\n")
    print(destination)


if __name__ == "__main__":
    main()
