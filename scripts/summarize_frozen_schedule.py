#!/usr/bin/env python3
"""Aggregate valid attempts from a frozen benchmark schedule."""

from __future__ import annotations

import argparse
import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any, Mapping, Sequence

from benchmark_harness import sha256_file
from run_frozen_schedule import _load_mapping, build_schedule, load_verified_schedule


AGGREGATE_SCHEMA_VERSION = "1.0.0"


def successful_attempts(run_root: Path) -> list[tuple[Path, dict[str, Any]]]:
    attempts = []
    for path in sorted(run_root.glob("attempt-*/run.json")):
        record = json.loads(path.read_text())
        if record.get("success") is True:
            attempts.append((path.parent, record))
    return attempts


def collect_rows(
    manifest: Mapping[str, Any], schedule: Mapping[str, Any], runs_root: Path
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    specs = {str(run["id"]): run for run in manifest["runs"]}
    rows: list[dict[str, Any]] = []
    exclusions: list[dict[str, Any]] = []
    for entry in schedule["entries"]:
        condition_id = str(entry["condition_id"])
        attempts = successful_attempts(runs_root / str(entry["run_id"]))
        if len(attempts) != 1:
            exclusions.append(
                {
                    "condition_id": condition_id,
                    "run_id": entry["run_id"],
                    "reason": "missing_success" if not attempts else "multiple_successful_attempts",
                    "successful_attempts": [str(path) for path, _ in attempts],
                }
            )
            continue
        attempt_path, record = attempts[0]
        metadata = dict(specs[condition_id]["spec"].get("metadata", {}))
        output_counts = [int(item["records"]) for item in record["output_counts"]]
        if output_counts and len(set(output_counts)) != 1:
            exclusions.append(
                {
                    "condition_id": condition_id,
                    "run_id": entry["run_id"],
                    "reason": "paired_output_count_mismatch",
                    "output_counts": output_counts,
                }
            )
            continue
        normalized = [item.get("normalized_sha256") for item in record["outputs"]]
        if not normalized or any(value is None for value in normalized):
            exclusions.append(
                {
                    "condition_id": condition_id,
                    "run_id": entry["run_id"],
                    "reason": "missing_normalized_output_digest",
                }
            )
            continue
        input_records = int(metadata["input_records"])
        rows.append(
            {
                "condition_id": condition_id,
                "run_id": entry["run_id"],
                "attempt": int(record["attempt"]),
                "attempt_dir": str(attempt_path),
                "dataset": metadata["dataset"],
                "tool": metadata["tool"],
                "execution_mode": metadata["execution_mode"],
                "threads": int(metadata["threads"]),
                "replicate": int(metadata["replicate"]),
                "input_records": input_records,
                "emitted_records": output_counts[0] if output_counts else 0,
                "wall_seconds": float(record["wall_seconds"]),
                "user_cpu_seconds": float(record.get("user_cpu_seconds", 0.0)),
                "system_cpu_seconds": float(record.get("system_cpu_seconds", 0.0)),
                "peak_rss_kib": int(record["peak_rss_kib"]),
                "input_records_per_second": input_records / float(record["wall_seconds"]),
                "output_bytes": sum(int(item["bytes"]) for item in record["outputs"]),
                "normalized_output_sha256": ";".join(str(value) for value in normalized),
            }
        )
    return rows, exclusions


def summarize(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        key = (row["dataset"], row["tool"], row["execution_mode"], row["threads"])
        groups[key].append(row)
    summaries = []
    for key, values in sorted(groups.items()):
        wall = [float(item["wall_seconds"]) for item in values]
        throughput = [float(item["input_records_per_second"]) for item in values]
        rss = [int(item["peak_rss_kib"]) for item in values]
        summaries.append(
            {
                "dataset": key[0],
                "tool": key[1],
                "execution_mode": key[2],
                "threads": key[3],
                "replicates": len(values),
                "wall_seconds_mean": statistics.fmean(wall),
                "wall_seconds_median": statistics.median(wall),
                "wall_seconds_sd": statistics.stdev(wall) if len(wall) > 1 else 0.0,
                "input_records_per_second_mean": statistics.fmean(throughput),
                "input_records_per_second_median": statistics.median(throughput),
                "peak_rss_kib_mean": statistics.fmean(rss),
                "peak_rss_kib_max": max(rss),
            }
        )
    return summaries


def correctness(rows: list[dict[str, Any]]) -> dict[str, Any]:
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[str(row["dataset"])].append(row)
    datasets = []
    for dataset, values in sorted(groups.items()):
        digests = sorted({str(item["normalized_output_sha256"]) for item in values})
        emitted_counts = sorted({int(item["emitted_records"]) for item in values})
        datasets.append(
            {
                "dataset": dataset,
                "conditions": len(values),
                "normalized_digest_count": len(digests),
                "normalized_output_sha256": digests,
                "emitted_record_counts": emitted_counts,
                "identical_across_modes_threads_replicates": len(digests) == 1,
            }
        )
    return {
        "all_identical": bool(datasets)
        and all(item["identical_across_modes_threads_replicates"] for item in datasets),
        "datasets": datasets,
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--schedule", type=Path, required=True)
    parser.add_argument("--runs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)

    manifest_path = args.manifest.resolve()
    manifest = _load_mapping(manifest_path)
    expected = build_schedule(manifest, manifest_path)
    schedule = load_verified_schedule(args.schedule, expected)
    rows, exclusions = collect_rows(manifest, schedule, args.runs.resolve())
    summaries = summarize(rows)
    correctness_report = correctness(rows)
    args.output.mkdir(parents=True, exist_ok=True)
    aggregate = {
        "schema_version": AGGREGATE_SCHEMA_VERSION,
        "manifest_sha256": sha256_file(manifest_path),
        "schedule_sha256": sha256_file(args.schedule),
        "scheduled_conditions": len(schedule["entries"]),
        "valid_conditions": len(rows),
        "excluded_conditions": len(exclusions),
        "correctness": correctness_report,
        "summaries": summaries,
        "exclusions": exclusions,
    }
    (args.output / "aggregate.json").write_text(
        json.dumps(aggregate, indent=2, sort_keys=True) + "\n"
    )
    (args.output / "correctness.json").write_text(
        json.dumps(correctness_report, indent=2, sort_keys=True) + "\n"
    )
    write_csv(args.output / "runs.csv", rows)
    write_csv(args.output / "summary.csv", summaries)
    print(json.dumps({"valid": len(rows), "excluded": len(exclusions), **correctness_report}))
    return 0 if not exclusions and correctness_report["all_identical"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
