#!/usr/bin/env python3
"""Aggregate valid attempts from a frozen benchmark schedule."""

from __future__ import annotations

import argparse
import csv
import json
import statistics
from collections import defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from benchmark_harness import sha256_file
from run_frozen_schedule import _load_mapping, build_schedule, load_verified_schedule

AGGREGATE_SCHEMA_VERSION = "1.2.0"


def successful_attempts(run_root: Path) -> list[tuple[Path, dict[str, Any]]]:
    attempts = []
    for path in sorted(run_root.glob("attempt-*/run.json")):
        record = json.loads(path.read_text())
        if record.get("success") is True:
            attempts.append((path.parent, record))
    return attempts


def collect_rows(
    manifest: Mapping[str, Any],
    schedule: Mapping[str, Any],
    runs_root: Path,
    datasets: frozenset[str] | None = None,
    tools: frozenset[str] | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    specs = {str(run["id"]): run for run in manifest["runs"]}
    rows: list[dict[str, Any]] = []
    exclusions: list[dict[str, Any]] = []
    for entry in schedule["entries"]:
        condition_id = str(entry["condition_id"])
        metadata = dict(specs[condition_id]["spec"].get("metadata", {}))
        if datasets and str(metadata.get("dataset", "")) not in datasets:
            continue
        if tools and str(metadata.get("tool", "")) not in tools:
            continue
        attempts = successful_attempts(runs_root / str(entry["run_id"]))
        if len(attempts) != 1:
            exclusions.append(
                {
                    "condition_id": condition_id,
                    "run_id": entry["run_id"],
                    "reason": "missing_success"
                    if not attempts
                    else "multiple_successful_attempts",
                    "successful_attempts": [str(path) for path, _ in attempts],
                }
            )
            continue
        attempt_path, record = attempts[0]
        measurement_track = str(metadata.get("measurement_track", "legacy"))
        output_counts = [int(item["records"]) for item in record["output_counts"]]
        output_length_validity = [
            {
                "path": item["path"],
                "records": int(item["records"]),
                "min_sequence_length": item.get("min_sequence_length"),
                "max_sequence_length": item.get("max_sequence_length"),
                "sequence_length_counts": item.get("sequence_length_counts", {}),
                "nominal_sequence_lengths": item.get(
                    "nominal_sequence_lengths", []
                ),
                "non_nominal_sequence_records": int(
                    item.get("non_nominal_sequence_records", 0)
                ),
            }
            for item in record["output_counts"]
        ]
        if (
            output_counts
            and len(set(output_counts)) != 1
            and not bool(metadata.get("allow_output_count_mismatch", False))
        ):
            exclusions.append(
                {
                    "condition_id": condition_id,
                    "run_id": entry["run_id"],
                    "reason": "paired_output_count_mismatch",
                    "output_counts": output_counts,
                }
            )
            continue
        normalized = [
            item["normalized_sha256"]
            for item in record["outputs"]
            if item.get("normalized_sha256") is not None
        ]
        if measurement_track != "timing" and not normalized:
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
                "measurement_track": measurement_track,
                "sequence_output_policy": metadata.get("sequence_output_policy"),
                "execution_mode": metadata["execution_mode"],
                "threads": int(metadata["threads"]),
                "replicate": int(metadata["replicate"]),
                "input_records": input_records,
                "emitted_records": output_counts[0] if output_counts else None,
                "non_nominal_sequence_records": sum(
                    item["non_nominal_sequence_records"]
                    for item in output_length_validity
                ),
                "output_length_validity_json": json.dumps(
                    output_length_validity, separators=(",", ":"), sort_keys=True
                ),
                "wall_seconds": float(record["wall_seconds"]),
                "user_cpu_seconds": float(record.get("user_cpu_seconds", 0.0)),
                "system_cpu_seconds": float(record.get("system_cpu_seconds", 0.0)),
                "peak_rss_kib": int(record["peak_rss_kib"]),
                "input_records_per_second": input_records
                / float(record["wall_seconds"]),
                "output_bytes": sum(int(item["bytes"]) for item in record["outputs"]),
                "normalized_output_sha256": (
                    ";".join(str(value) for value in normalized)
                    if normalized
                    else None
                ),
            }
        )
    return rows, exclusions


def summarize(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        key = (
            row["measurement_track"],
            row["sequence_output_policy"],
            row["dataset"],
            row["tool"],
            row["execution_mode"],
            row["threads"],
        )
        groups[key].append(row)
    summaries = []
    for key, values in sorted(groups.items()):
        wall = [float(item["wall_seconds"]) for item in values]
        throughput = [float(item["input_records_per_second"]) for item in values]
        rss = [int(item["peak_rss_kib"]) for item in values]
        summaries.append(
            {
                "measurement_track": key[0],
                "sequence_output_policy": key[1],
                "dataset": key[2],
                "tool": key[3],
                "execution_mode": key[4],
                "threads": key[5],
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
        if row.get("normalized_output_sha256") is None:
            continue
        groups[str(row["dataset"])].append(row)
    datasets = []
    for dataset, values in sorted(groups.items()):
        digests = sorted({str(item["normalized_output_sha256"]) for item in values})
        emitted_counts = sorted({int(item["emitted_records"]) for item in values})
        condition_groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(
            list
        )
        for item in values:
            condition_groups[(str(item["tool"]), str(item["execution_mode"]))].append(
                item
            )
        determinism = []
        for (tool, execution_mode), condition_values in sorted(
            condition_groups.items()
        ):
            condition_digests = sorted(
                {str(item["normalized_output_sha256"]) for item in condition_values}
            )
            condition_counts = sorted(
                {int(item["emitted_records"]) for item in condition_values}
            )
            non_nominal_counts = sorted(
                {
                    int(item.get("non_nominal_sequence_records", 0))
                    for item in condition_values
                }
            )
            determinism.append(
                {
                    "tool": tool,
                    "execution_mode": execution_mode,
                    "conditions": len(condition_values),
                    "normalized_digest_count": len(condition_digests),
                    "emitted_record_counts": condition_counts,
                    "non_nominal_sequence_record_counts": non_nominal_counts,
                    "deterministic_across_threads_replicates": len(condition_digests)
                    == 1
                    and len(condition_counts) == 1,
                }
            )
        datasets.append(
            {
                "dataset": dataset,
                "conditions": len(values),
                "normalized_digest_count": len(digests),
                "normalized_output_sha256": digests,
                "emitted_record_counts": emitted_counts,
                "identical_across_modes_threads_replicates": len(digests) == 1,
                "all_tools_modes_deterministic": all(
                    item["deterministic_across_threads_replicates"]
                    for item in determinism
                ),
                "determinism": determinism,
            }
        )
    all_deterministic = bool(datasets) and all(
        item["all_tools_modes_deterministic"] for item in datasets
    )
    return {
        "all_deterministic": all_deterministic,
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
    parser.add_argument(
        "--dataset",
        action="append",
        help="summarize only this dataset/technology block; may be repeated",
    )
    parser.add_argument(
        "--tool",
        action="append",
        help="summarize only this tool; may be repeated",
    )
    args = parser.parse_args(argv)

    manifest_path = args.manifest.resolve()
    manifest = _load_mapping(manifest_path)
    expected = build_schedule(manifest, manifest_path)
    schedule = load_verified_schedule(args.schedule, expected)
    selected_datasets = frozenset(args.dataset or ())
    selected_tools = frozenset(args.tool or ())
    known_datasets = {
        str(run["spec"].get("metadata", {}).get("dataset", ""))
        for run in manifest["runs"]
    }
    unknown = sorted(selected_datasets - known_datasets)
    if unknown:
        parser.error(
            f"unknown dataset block(s): {', '.join(unknown)}; "
            f"available: {', '.join(sorted(known_datasets))}"
        )
    known_tools = {
        str(run["spec"].get("metadata", {}).get("tool", ""))
        for run in manifest["runs"]
    }
    unknown_tools = sorted(selected_tools - known_tools)
    if unknown_tools:
        parser.error(
            f"unknown tool(s): {', '.join(unknown_tools)}; "
            f"available: {', '.join(sorted(known_tools))}"
        )
    rows, exclusions = collect_rows(
        manifest, schedule, args.runs.resolve(), selected_datasets, selected_tools
    )
    summaries = summarize(rows)
    correctness_report = correctness(rows)
    require_cross_tool_identity = bool(
        manifest.get("study", {}).get("require_cross_tool_identity", False)
    )
    args.output.mkdir(parents=True, exist_ok=True)
    aggregate = {
        "schema_version": AGGREGATE_SCHEMA_VERSION,
        "manifest_sha256": sha256_file(manifest_path),
        "schedule_sha256": sha256_file(args.schedule),
        "scheduled_conditions": sum(
            1
            for run in manifest["runs"]
            if not selected_datasets
            or str(run["spec"].get("metadata", {}).get("dataset", ""))
            in selected_datasets
            if not selected_tools
            or str(run["spec"].get("metadata", {}).get("tool", ""))
            in selected_tools
        ),
        "datasets": sorted(selected_datasets or known_datasets),
        "tools": sorted(selected_tools or known_tools),
        "valid_conditions": len(rows),
        "excluded_conditions": len(exclusions),
        "correctness": correctness_report,
        "require_cross_tool_identity": require_cross_tool_identity,
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
    print(
        json.dumps(
            {"valid": len(rows), "excluded": len(exclusions), **correctness_report}
        )
    )
    has_correctness_rows = any(
        row.get("normalized_output_sha256") is not None for row in rows
    )
    correctness_ok = not has_correctness_rows or (
        correctness_report["all_deterministic"]
        and (
            correctness_report["all_identical"]
            or not require_cross_tool_identity
        )
    )
    return 0 if not exclusions and correctness_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
