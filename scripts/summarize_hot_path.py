#!/usr/bin/env python3
"""Aggregate development hot-path runs without copying numbers by hand."""

from __future__ import annotations

import argparse
import json
import statistics
from pathlib import Path
from typing import Any, Iterable, Sequence


def load_run(run_dir: Path) -> dict[str, Any]:
    run = json.loads((run_dir / "run.json").read_text())
    if not run.get("success"):
        raise ValueError(f"run did not succeed: {run_dir}")
    payload = json.loads((run_dir / "stdout.txt").read_text())
    required = {"mode", "statistics", "reads", "threads", "seconds"}
    missing = required.difference(payload)
    if missing:
        raise ValueError(f"{run_dir} payload missing {sorted(missing)}")
    seconds = [float(value) for value in payload["seconds"]]
    if not seconds or any(value <= 0 for value in seconds):
        raise ValueError(f"{run_dir} has invalid timing values")
    return {
        "run_id": run["run_id"],
        "attempt": run["attempt"],
        "run_dir": str(run_dir.resolve()),
        "peak_rss_kib": run.get("peak_rss_kib"),
        "user_cpu_seconds": run.get("user_cpu_seconds"),
        "system_cpu_seconds": run.get("system_cpu_seconds"),
        "mode": payload["mode"],
        "statistics": bool(payload["statistics"]),
        "statistics_level": payload.get(
            "statistics_level", "detailed" if payload["statistics"] else "off"
        ),
        "reads": int(payload["reads"]),
        "threads": int(payload["threads"]),
        "seconds": seconds,
    }


def summarize(run_dirs: Iterable[Path]) -> dict[str, Any]:
    runs = [load_run(path) for path in run_dirs]
    conditions: list[dict[str, Any]] = []
    for run in runs:
        seconds = run["seconds"]
        mean_seconds = statistics.fmean(seconds)
        conditions.append(
            {
                **run,
                "n": len(seconds),
                "mean_seconds": mean_seconds,
                "median_seconds": statistics.median(seconds),
                "sample_standard_deviation_seconds": (
                    statistics.stdev(seconds) if len(seconds) > 1 else 0.0
                ),
                "reads_per_second_from_mean": run["reads"] / mean_seconds,
            }
        )

    effects: list[dict[str, Any]] = []
    modes = sorted({condition["mode"] for condition in conditions})
    for mode in modes:
        matching = [condition for condition in conditions if condition["mode"] == mode]
        off = [condition for condition in matching if not condition["statistics"]]
        on = [condition for condition in matching if condition["statistics"]]
        if len(off) != 1 or not on:
            raise ValueError(
                f"mode {mode!r} requires exactly one statistics-off run and at least one enabled run"
            )
        levels = [condition["statistics_level"] for condition in on]
        if len(levels) != len(set(levels)):
            raise ValueError(f"mode {mode!r} has duplicate enabled statistics levels")
        off_condition = off[0]
        for on_condition in sorted(on, key=lambda condition: condition["statistics_level"]):
            if (off_condition["reads"], off_condition["threads"]) != (
                on_condition["reads"],
                on_condition["threads"],
            ):
                raise ValueError(f"mode {mode!r} conditions are not comparable")
            effects.append(
                {
                    "mode": mode,
                    "statistics_level": on_condition["statistics_level"],
                    "reads": off_condition["reads"],
                    "threads": off_condition["threads"],
                    "statistics_off_run_id": off_condition["run_id"],
                    "statistics_on_run_id": on_condition["run_id"],
                    "mean_time_overhead_pct": 100
                    * (on_condition["mean_seconds"] / off_condition["mean_seconds"] - 1),
                    "median_time_overhead_pct": 100
                    * (on_condition["median_seconds"] / off_condition["median_seconds"] - 1),
                    "mean_time_reduction_pct": 100
                    * (on_condition["mean_seconds"] - off_condition["mean_seconds"])
                    / on_condition["mean_seconds"],
                    "median_time_reduction_pct": 100
                    * (on_condition["median_seconds"] - off_condition["median_seconds"])
                    / on_condition["median_seconds"],
                    "throughput_gain_from_mean_pct": 100
                    * (on_condition["mean_seconds"] / off_condition["mean_seconds"] - 1),
                    "throughput_gain_from_median_pct": 100
                    * (on_condition["median_seconds"] / off_condition["median_seconds"] - 1),
                }
            )
    return {
        "schema_version": "1.1.0",
        "status": "development-exploratory",
        "conditions": conditions,
        "effects": effects,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dirs", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)

    result = summarize(args.run_dirs)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
