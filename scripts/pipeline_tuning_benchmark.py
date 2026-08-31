#!/usr/bin/env python3
"""Measure staged-pipeline batch/window tradeoffs on manuscript geometries."""

from __future__ import annotations

import argparse
import json
import random
import shutil
import statistics
from pathlib import Path

from protocol_integration_benchmark import (
    CONFIGS,
    PROTOCOLS,
    generate_protocol,
    run_once,
    sha256_file,
)


SCHEMA_VERSION = "1.0.0"


def configurations(workers: int) -> list[dict[str, object]]:
    result = []
    for batch_size in (256, 512, 1024, 2048, 4096):
        result.append(
            {
                "name": f"tight_b{batch_size}",
                "threads": workers,
                "batch_size": batch_size,
                "queue_capacity": 2,
                "max_in_flight_batches": workers,
                "flags": [
                    "--staged-pipeline",
                    "--batch-size",
                    batch_size,
                    "--queue-capacity",
                    2,
                    "--max-in-flight-batches",
                    workers,
                ],
            }
        )
    result.extend(
        (
            {
                "name": "balanced_b1024",
                "threads": workers,
                "batch_size": 1024,
                "queue_capacity": max(2, workers // 2),
                "max_in_flight_batches": workers + max(2, workers // 2),
                "flags": [
                    "--staged-pipeline",
                    "--batch-size",
                    1024,
                    "--queue-capacity",
                    max(2, workers // 2),
                    "--max-in-flight-batches",
                    workers + max(2, workers // 2),
                ],
            },
            {
                "name": "current_default",
                "threads": workers,
                "batch_size": 2048,
                "queue_capacity": workers * 2,
                "max_in_flight_batches": workers * 3,
                "flags": ["--staged-pipeline"],
            },
        )
    )
    return result


def summarize(runs: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = {}
    for run in runs:
        key = (str(run["protocol"]), str(run["mode"]))
        grouped.setdefault(key, []).append(run)
    result = []
    for (protocol, mode), values in sorted(grouped.items()):
        throughputs = [float(value["reads_per_second"]) for value in values]
        item = {
            "protocol": protocol,
            "mode": mode,
            "batch_size": values[0]["batch_size"],
            "queue_capacity": values[0]["queue_capacity"],
            "max_in_flight_batches": values[0]["max_in_flight_batches"],
            "replicates": len(values),
            "median_reads_per_second": statistics.median(throughputs),
            "min_reads_per_second": min(throughputs),
            "max_reads_per_second": max(throughputs),
            "median_peak_rss_kib": statistics.median(
                float(value["peak_rss_kib"]) for value in values
            ),
        }
        result.append(item)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--reads", type=int, default=100_000)
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    if args.reads <= 0 or args.replicates <= 0 or args.workers <= 0:
        parser.error("reads, replicates, and workers must be positive")
    binary = args.binary.resolve()
    if not binary.is_file():
        parser.error(f"binary does not exist: {binary}")

    outdir = args.outdir.resolve()
    data_dir = outdir / "inputs"
    work_dir = outdir / "work"
    data_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    generated = {
        protocol.name: generate_protocol(protocol, args.reads, data_dir)
        for protocol in PROTOCOLS
    }
    configs = configurations(args.workers)
    schedule = [
        (protocol, config, replicate)
        for replicate in range(1, args.replicates + 1)
        for protocol in PROTOCOLS
        for config in configs
    ]
    random.Random("seqproc-pipeline-tuning-v1").shuffle(schedule)

    runs = []
    reference_digests: dict[str, str] = {}
    for index, (protocol, config, replicate) in enumerate(schedule, 1):
        print(
            f"[{index}/{len(schedule)}] {protocol.name} {config['name']} "
            f"replicate={replicate}",
            flush=True,
        )
        r1, r2 = generated[protocol.name]
        run = run_once(binary, protocol, r1, r2, config, replicate, work_dir)
        for key in ("batch_size", "queue_capacity", "max_in_flight_batches"):
            run[key] = config[key]
        previous = reference_digests.setdefault(
            protocol.name, str(run["output_multiset_sha256"])
        )
        if run["output_multiset_sha256"] != previous:
            raise RuntimeError(f"output divergence for {protocol.name}")
        if run["output_reads"] != args.reads:
            raise RuntimeError(f"unexpected output count for {protocol.name}")
        runs.append(run)

    artifact = {
        "schema_version": SCHEMA_VERSION,
        "scope": "protocol-shaped synthetic staged-pipeline tuning benchmark",
        "binary": {"path": str(binary), "sha256": sha256_file(binary)},
        "parameters": {
            "reads_per_protocol": args.reads,
            "replicates": args.replicates,
            "workers": args.workers,
            "randomized_schedule": True,
        },
        "inputs": {
            protocol.name: {
                "geometry": str(CONFIGS / protocol.geometry),
                "geometry_sha256": sha256_file(CONFIGS / protocol.geometry),
                "r1_sha256": sha256_file(generated[protocol.name][0]),
                "r2_sha256": (
                    sha256_file(generated[protocol.name][1])
                    if generated[protocol.name][1]
                    else None
                ),
            }
            for protocol in PROTOCOLS
        },
        "reference_output_multiset_sha256": reference_digests,
        "runs": runs,
        "summary": summarize(runs),
    }
    artifact_path = outdir / "pipeline-tuning-results.json"
    artifact_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    shutil.rmtree(work_dir)
    print(f"wrote {artifact_path}")


if __name__ == "__main__":
    main()
