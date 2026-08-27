#!/usr/bin/env python3
"""Generate, validate, and execute an immutable randomized benchmark schedule.

The manifest is YAML or JSON and contains an explicit list of run specifications.
Each run belongs to a block (normally dataset/thread/replicate); runs are shuffled
within blocks, and blocks are shuffled, from the manifest's recorded seed.
"""

from __future__ import annotations

import argparse
import contextlib
import fcntl
import json
import os
import random
import re
import sys
import time
from collections import defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import yaml
from benchmark_harness import (
    HarnessError,
    canonical_json,
    execute_spec,
    prepare_spec,
    sha256_bytes,
    sha256_file,
)

SCHEDULE_SCHEMA_VERSION = "1.0.0"
GENERATOR_VERSION = "seqproc-frozen-schedule-1"
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
COMMIT_RE = re.compile(r"^[0-9a-f]{40}$")


class ScheduleError(RuntimeError):
    """The manifest, schedule, or execution state is invalid."""


def _load_mapping(path: Path) -> dict[str, Any]:
    try:
        value = yaml.safe_load(path.read_text())
    except (OSError, yaml.YAMLError) as error:
        raise ScheduleError(f"cannot read {path}: {error}") from error
    if not isinstance(value, dict):
        raise ScheduleError(f"{path} must contain a mapping at its root")
    return value


def _walk(value: Any, location: str = "$"):
    yield location, value
    if isinstance(value, Mapping):
        for key, item in value.items():
            yield from _walk(item, f"{location}.{key}")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            yield from _walk(item, f"{location}[{index}]")


def _path_from_manifest(value: str, manifest_dir: Path) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = manifest_dir / path
    return path.resolve()


def _validate_sha256(value: Any, location: str) -> str:
    digest = str(value)
    if not SHA256_RE.fullmatch(digest):
        raise ScheduleError(f"{location} must be a lowercase SHA-256 digest")
    return digest


def _validate_publication_spec(spec: Mapping[str, Any], location: str) -> None:
    for collection in ("inputs", "configs", "executables"):
        for index, entry in enumerate(spec.get(collection, [])):
            if not isinstance(entry, Mapping):
                raise ScheduleError(
                    f"{location}.{collection}[{index}] must declare path and sha256"
                )
            _validate_sha256(
                entry.get("sha256"), f"{location}.{collection}[{index}].sha256"
            )
    for index, repository in enumerate(spec.get("repositories", [])):
        if not isinstance(repository, Mapping):
            raise ScheduleError(f"{location}.repositories[{index}] must be a mapping")
        commit = str(repository.get("commit", ""))
        if not COMMIT_RE.fullmatch(commit):
            raise ScheduleError(
                f"{location}.repositories[{index}].commit must be a full commit digest"
            )
        if repository.get("allow_dirty"):
            raise ScheduleError(
                f"{location}.repositories[{index}] permits a dirty tree"
            )


def validate_manifest(
    manifest: Mapping[str, Any], manifest_path: Path
) -> list[dict[str, Any]]:
    for location, value in _walk(manifest):
        if isinstance(value, str) and value.strip() == "REQUIRED":
            raise ScheduleError(f"incomplete manifest field at {location}")
        if isinstance(value, str) and ":latest" in value.lower():
            raise ScheduleError(f"floating container tag at {location}: {value}")

    if str(manifest.get("schema_version")) != "1.0.0":
        raise ScheduleError("manifest schema_version must be 1.0.0")
    mode = manifest.get("mode")
    if mode not in ("development", "publication"):
        raise ScheduleError("manifest mode must be 'development' or 'publication'")
    study = manifest.get("study")
    if not isinstance(study, Mapping) or not isinstance(study.get("random_seed"), int):
        raise ScheduleError("study.random_seed must be an integer")

    manifest_dir = manifest_path.resolve().parent
    artifact_digests: dict[str, str] = {}
    for index, artifact in enumerate(manifest.get("artifacts", [])):
        if not isinstance(artifact, Mapping) or "path" not in artifact:
            raise ScheduleError(f"artifacts[{index}] must declare path and sha256")
        expected = _validate_sha256(
            artifact.get("sha256"), f"artifacts[{index}].sha256"
        )
        path = _path_from_manifest(str(artifact["path"]), manifest_dir)
        if not path.is_file():
            raise ScheduleError(f"frozen artifact does not exist: {path}")
        observed = sha256_file(path)
        if observed != expected:
            raise ScheduleError(
                f"frozen artifact digest mismatch for {path}: expected {expected}, observed {observed}"
            )
        path_key = str(path)
        if path_key in artifact_digests:
            raise ScheduleError(f"duplicate frozen artifact path: {path}")
        artifact_digests[path_key] = observed

    runs = manifest.get("runs")
    if not isinstance(runs, list) or not runs:
        raise ScheduleError("manifest runs must be a non-empty list")
    seen_ids: set[str] = set()
    normalized: list[dict[str, Any]] = []
    execution = manifest.get("execution", {})
    if not isinstance(execution, Mapping):
        raise ScheduleError("execution must be a mapping")
    default_timeout = execution.get("timeout_seconds")
    environment_allowlist = execution.get("sanitized_environment_allowlist")

    for index, run in enumerate(runs):
        location = f"runs[{index}]"
        if not isinstance(run, Mapping):
            raise ScheduleError(f"{location} must be a mapping")
        condition_id = str(run.get("id", ""))
        block_id = str(run.get("block", ""))
        spec_value = run.get("spec")
        if not condition_id or condition_id in seen_ids:
            raise ScheduleError(
                f"{location}.id is empty or duplicated: {condition_id!r}"
            )
        if not block_id:
            raise ScheduleError(f"{location}.block is required")
        if not isinstance(spec_value, Mapping):
            raise ScheduleError(f"{location}.spec must be a mapping")
        seen_ids.add(condition_id)
        spec = dict(spec_value)
        if default_timeout is not None and "timeout_seconds" not in spec:
            spec["timeout_seconds"] = default_timeout
        if environment_allowlist is not None:
            spec.setdefault("environment_mode", "allowlist")
            spec.setdefault("environment_allowlist", environment_allowlist)
        metadata = dict(spec.get("metadata", {}))
        metadata.setdefault("condition_id", condition_id)
        metadata.setdefault("block_id", block_id)
        spec["metadata"] = metadata
        if mode == "publication":
            _validate_publication_spec(spec, f"{location}.spec")
        try:
            prepared = prepare_spec(spec)
        except HarnessError as error:
            raise ScheduleError(f"{location}: {error}") from error
        if mode == "publication":
            frozen_files = [
                *prepared["identity"]["executables"],
                *prepared["identity"]["inputs"],
                *prepared["identity"]["configs"],
            ]
            for frozen in frozen_files:
                path_key = str(Path(frozen["path"]).resolve())
                expected = artifact_digests.get(path_key)
                if expected is None:
                    raise ScheduleError(
                        f"{location} uses an executable/input/config absent from artifacts: {path_key}"
                    )
                if expected != frozen["sha256"]:
                    raise ScheduleError(
                        f"{location} artifact digest disagrees with run specification for {path_key}"
                    )
        normalized.append(
            {
                "condition_id": condition_id,
                "block_id": block_id,
                "run_id": prepared["run_id"],
                "spec": spec,
            }
        )
    return normalized


def build_schedule(manifest: Mapping[str, Any], manifest_path: Path) -> dict[str, Any]:
    runs = validate_manifest(manifest, manifest_path)
    seed = int(manifest["study"]["random_seed"])
    rng = random.Random(seed)
    blocks: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for run in runs:
        blocks[run["block_id"]].append(run)
    block_ids = sorted(blocks)
    rng.shuffle(block_ids)
    entries: list[dict[str, Any]] = []
    for block_id in block_ids:
        block = sorted(blocks[block_id], key=lambda item: item["condition_id"])
        rng.shuffle(block)
        for item in block:
            entries.append(
                {
                    "ordinal": len(entries) + 1,
                    "condition_id": item["condition_id"],
                    "block_id": block_id,
                    "run_id": item["run_id"],
                }
            )
    return {
        "schema_version": SCHEDULE_SCHEMA_VERSION,
        "generator_version": GENERATOR_VERSION,
        "manifest_sha256": sha256_file(manifest_path),
        "random_seed": seed,
        "entries": entries,
    }


def _write_new(path: Path, content: bytes) -> None:
    if path.exists():
        raise ScheduleError(f"refusing to overwrite immutable file: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    temporary.write_bytes(content)
    os.replace(temporary, path)


def write_schedule(schedule: Mapping[str, Any], path: Path) -> str:
    content = json.dumps(schedule, indent=2, sort_keys=True).encode() + b"\n"
    digest = sha256_bytes(content)
    _write_new(path, content)
    _write_new(
        path.with_suffix(path.suffix + ".sha256"), f"{digest}  {path.name}\n".encode()
    )
    return digest


def load_verified_schedule(path: Path, expected: Mapping[str, Any]) -> dict[str, Any]:
    sidecar = path.with_suffix(path.suffix + ".sha256")
    if not sidecar.is_file():
        raise ScheduleError(f"missing schedule digest sidecar: {sidecar}")
    fields = sidecar.read_text().strip().split()
    if not fields or not SHA256_RE.fullmatch(fields[0]):
        raise ScheduleError(f"invalid schedule digest sidecar: {sidecar}")
    observed = sha256_file(path)
    if observed != fields[0]:
        raise ScheduleError(
            f"schedule digest mismatch: expected {fields[0]}, observed {observed}"
        )
    try:
        schedule = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ScheduleError(f"cannot parse schedule {path}: {error}") from error
    if canonical_json(schedule) != canonical_json(expected):
        raise ScheduleError("schedule does not match deterministic manifest expansion")
    return schedule


def _successful_attempt(output_root: Path, run_id: str) -> Path | None:
    run_root = output_root.resolve() / run_id
    for attempt in sorted(run_root.glob("attempt-*")):
        record = attempt / "run.json"
        if not record.is_file():
            continue
        try:
            if json.loads(record.read_text()).get("success") is True:
                return attempt
        except (OSError, json.JSONDecodeError):
            continue
    return None


def _append_log(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as handle:
        handle.write(json.dumps(value, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


@contextlib.contextmanager
def execution_lock(output_root: Path):
    """Hold a process-scoped exclusive lock for one benchmark output root."""
    output_root = output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    lock_path = output_root / ".coordinator.lock"
    handle = lock_path.open("a+")
    try:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            handle.seek(0)
            owner = handle.read().strip() or "unknown owner"
            raise ScheduleError(
                f"another coordinator holds {lock_path}: {owner}"
            ) from error
        handle.seek(0)
        handle.truncate()
        handle.write(
            json.dumps(
                {
                    "pid": os.getpid(),
                    "host": os.uname().nodename,
                    "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                },
                sort_keys=True,
            )
            + "\n"
        )
        handle.flush()
        os.fsync(handle.fileno())
        yield
    finally:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
        finally:
            handle.close()


def _execute_schedule_unlocked(
    manifest: Mapping[str, Any],
    manifest_path: Path,
    schedule: Mapping[str, Any],
    output_root: Path,
    max_runs: int | None = None,
    datasets: frozenset[str] | None = None,
    tools: frozenset[str] | None = None,
) -> tuple[int, int, int]:
    runs = validate_manifest(manifest, manifest_path)
    by_id = {item["condition_id"]: item for item in runs}
    mode = manifest["mode"]
    if max_runs is not None and mode == "publication":
        raise ScheduleError("--max-runs is forbidden in publication mode")
    known_datasets = {
        str(item["spec"].get("metadata", {}).get("dataset", "")) for item in runs
    }
    if datasets:
        unknown = sorted(datasets - known_datasets)
        if unknown:
            raise ScheduleError(
                f"unknown dataset block(s): {', '.join(unknown)}; "
                f"available: {', '.join(sorted(known_datasets))}"
            )
    known_tools = {
        str(item["spec"].get("metadata", {}).get("tool", "")) for item in runs
    }
    if tools:
        unknown = sorted(tools - known_tools)
        if unknown:
            raise ScheduleError(
                f"unknown tool(s): {', '.join(unknown)}; "
                f"available: {', '.join(sorted(known_tools))}"
            )
    selected = [
        entry
        for entry in schedule["entries"]
        if (
            not datasets
            or str(
                by_id[entry["condition_id"]]["spec"]
                .get("metadata", {})
                .get("dataset", "")
            )
            in datasets
        )
        and (
            not tools
            or str(
                by_id[entry["condition_id"]]["spec"]
                .get("metadata", {})
                .get("tool", "")
            )
            in tools
        )
    ]
    if max_runs is not None:
        selected = selected[:max_runs]
    completed = skipped = failed = 0
    log_path = output_root.resolve() / "execution-log.jsonl"
    for entry in selected:
        item = by_id[entry["condition_id"]]
        if item["run_id"] != entry["run_id"]:
            raise ScheduleError(f"run identity changed for {entry['condition_id']}")
        previous = _successful_attempt(output_root, item["run_id"])
        if previous is not None:
            skipped += 1
            _append_log(
                log_path,
                {
                    "event": "resume-skip",
                    "condition_id": entry["condition_id"],
                    "run_id": item["run_id"],
                    "attempt_dir": str(previous),
                    "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                },
            )
            continue
        _append_log(
            log_path,
            {
                "event": "start",
                "condition_id": entry["condition_id"],
                "run_id": item["run_id"],
                "ordinal": entry["ordinal"],
                "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            },
        )
        result, run_dir = execute_spec(item["spec"], output_root)
        completed += 1
        if not result["success"]:
            failed += 1
        _append_log(
            log_path,
            {
                "event": "finish",
                "condition_id": entry["condition_id"],
                "run_id": item["run_id"],
                "attempt": result["attempt"],
                "attempt_dir": str(run_dir),
                "success": result["success"],
                "exit_code": result["exit_code"],
                "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            },
        )
    return completed, skipped, failed


def execute_schedule(
    manifest: Mapping[str, Any],
    manifest_path: Path,
    schedule: Mapping[str, Any],
    output_root: Path,
    max_runs: int | None = None,
    datasets: frozenset[str] | None = None,
    tools: frozenset[str] | None = None,
) -> tuple[int, int, int]:
    with execution_lock(output_root):
        return _execute_schedule_unlocked(
            manifest, manifest_path, schedule, output_root, max_runs, datasets, tools
        )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--schedule", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("benchmark_results/frozen"))
    parser.add_argument(
        "--generate", action="store_true", help="write schedule and digest, then exit"
    )
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument(
        "--max-runs", type=int, help="development-mode smoke-test limit"
    )
    parser.add_argument(
        "--dataset",
        action="append",
        help=(
            "execute or validate only this dataset/technology block; may be repeated. "
            "The relative order remains the order frozen in the full schedule"
        ),
    )
    parser.add_argument(
        "--tool",
        action="append",
        help=(
            "execute or validate only this tool; may be repeated. The relative "
            "order remains the order frozen in the full schedule"
        ),
    )
    args = parser.parse_args(argv)
    if args.max_runs is not None and args.max_runs <= 0:
        parser.error("--max-runs must be positive")
    if args.generate and (args.dataset or args.tool):
        parser.error("--dataset/--tool cannot be combined with --generate")

    try:
        manifest_path = args.manifest.resolve()
        manifest = _load_mapping(manifest_path)
        expected_schedule = build_schedule(manifest, manifest_path)
        if args.generate:
            digest = write_schedule(expected_schedule, args.schedule)
            print(
                json.dumps(
                    {"entries": len(expected_schedule["entries"]), "sha256": digest}
                )
            )
            return 0
        schedule = load_verified_schedule(args.schedule, expected_schedule)
        selected_datasets = frozenset(args.dataset or ())
        selected_tools = frozenset(args.tool or ())
        runs_by_id = {run["id"]: run for run in manifest["runs"]}
        selected_entries = [
            entry
            for entry in schedule["entries"]
            if (
                not selected_datasets
                or runs_by_id[entry["condition_id"]]["spec"]["metadata"]["dataset"]
                in selected_datasets
            )
            and (
                not selected_tools
                or runs_by_id[entry["condition_id"]]["spec"]["metadata"]["tool"]
                in selected_tools
            )
        ]
        known_datasets = {
            str(run["spec"].get("metadata", {}).get("dataset", ""))
            for run in manifest["runs"]
        }
        unknown = sorted(selected_datasets - known_datasets)
        if unknown:
            raise ScheduleError(
                f"unknown dataset block(s): {', '.join(unknown)}; "
                f"available: {', '.join(sorted(known_datasets))}"
            )
        known_tools = {
            str(run["spec"].get("metadata", {}).get("tool", ""))
            for run in manifest["runs"]
        }
        unknown_tools = sorted(selected_tools - known_tools)
        if unknown_tools:
            raise ScheduleError(
                f"unknown tool(s): {', '.join(unknown_tools)}; "
                f"available: {', '.join(sorted(known_tools))}"
            )
        if args.validate_only:
            print(
                json.dumps(
                    {
                        "valid": True,
                        "mode": manifest["mode"],
                        "entries": len(selected_entries),
                        "datasets": sorted(selected_datasets or known_datasets),
                        "tools": sorted(selected_tools or known_tools),
                        "manifest_sha256": schedule["manifest_sha256"],
                        "schedule_sha256": sha256_file(args.schedule),
                    },
                    sort_keys=True,
                )
            )
            return 0
        completed, skipped, failed = execute_schedule(
            manifest,
            manifest_path,
            schedule,
            args.output,
            args.max_runs,
            selected_datasets,
            selected_tools,
        )
    except KeyboardInterrupt:
        print("run-frozen-schedule: interrupted", file=sys.stderr)
        return 130
    except (OSError, ScheduleError, HarnessError) as error:
        print(f"run-frozen-schedule: {error}", file=sys.stderr)
        return 2
    print(json.dumps({"completed": completed, "skipped": skipped, "failed": failed}))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
