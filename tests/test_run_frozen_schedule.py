import json
import sys
from pathlib import Path

import pytest
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from benchmark_harness import sha256_file
from run_frozen_schedule import (
    ScheduleError,
    build_schedule,
    execute_schedule,
    execution_lock,
    load_verified_schedule,
    validate_manifest,
    write_schedule,
)
from summarize_frozen_schedule import collect_rows, correctness, summarize


def make_manifest(tmp_path: Path, fail: bool = False):
    input_path = tmp_path / "input.txt"
    input_path.write_text("input\n")
    artifact = tmp_path / "artifact.txt"
    artifact.write_text("frozen\n")
    runs = []
    for tool in ("alpha", "beta"):
        code = (
            "raise SystemExit(7)"
            if fail and tool == "beta"
            else (
                "from pathlib import Path; import os, sys; "
                "assert 'SCHEDULER_TEST_SECRET' not in os.environ; "
                "Path(sys.argv[1], 'out.fastq').write_text('@r\\nA\\n+\\nI\\n')"
            )
        )
        runs.append(
            {
                "id": f"dataset-t1-r1-{tool}",
                "block": "dataset-t1-r1",
                "spec": {
                    "name": tool,
                    "cwd": str(tmp_path),
                    "command": [sys.executable, "-c", code, "{run_dir}"],
                    "inputs": [str(input_path)],
                    "outputs": []
                    if fail and tool == "beta"
                    else [
                        {
                            "path": "{run_dir}/out.fastq",
                            "format": "fastq",
                            "normalize": "fastq_multiset",
                            "mate": 1,
                            "min_bytes": 1,
                        }
                    ],
                    "metadata": {
                        "dataset": "dataset",
                        "tool": tool,
                        "execution_mode": "default",
                        "replicate": 1,
                        "threads": 1,
                        "input_records": 1,
                    },
                },
            }
        )
    manifest = {
        "schema_version": "1.0.0",
        "mode": "development",
        "study": {"name": "test", "random_seed": 741211},
        "artifacts": [{"path": str(artifact), "sha256": sha256_file(artifact)}],
        "execution": {
            "timeout_seconds": 5,
            "sanitized_environment_allowlist": ["PATH"],
        },
        "runs": runs,
    }
    path = tmp_path / "manifest.yaml"
    path.write_text(yaml.safe_dump(manifest, sort_keys=False))
    return manifest, path


def test_schedule_is_deterministic_and_digest_protected(tmp_path):
    manifest, manifest_path = make_manifest(tmp_path)
    first = build_schedule(manifest, manifest_path)
    second = build_schedule(manifest, manifest_path)
    assert first == second
    assert sorted(entry["ordinal"] for entry in first["entries"]) == [1, 2]

    schedule_path = tmp_path / "schedule.json"
    digest = write_schedule(first, schedule_path)
    assert digest == sha256_file(schedule_path)
    assert load_verified_schedule(schedule_path, first) == first

    schedule_path.write_text(schedule_path.read_text() + " ")
    with pytest.raises(ScheduleError, match="digest mismatch"):
        load_verified_schedule(schedule_path, first)


def test_execute_and_idempotent_resume(tmp_path, monkeypatch):
    monkeypatch.setenv("SCHEDULER_TEST_SECRET", "must-not-leak")
    manifest, manifest_path = make_manifest(tmp_path)
    schedule = build_schedule(manifest, manifest_path)
    output = tmp_path / "results"

    assert execute_schedule(manifest, manifest_path, schedule, output) == (2, 0, 0)
    assert execute_schedule(manifest, manifest_path, schedule, output) == (0, 2, 0)
    events = [
        json.loads(line)
        for line in (output / "execution-log.jsonl").read_text().splitlines()
    ]
    assert [event["event"] for event in events].count("resume-skip") == 2

    rows, exclusions = collect_rows(manifest, schedule, output)
    assert len(rows) == 2
    assert exclusions == []
    assert len(summarize(rows)) == 2
    assert correctness(rows)["all_identical"] is True


def test_dataset_block_selection_preserves_frozen_relative_order(tmp_path):
    manifest, manifest_path = make_manifest(tmp_path)
    second_runs = []
    for run in manifest["runs"]:
        copied = json.loads(json.dumps(run))
        copied["id"] = copied["id"].replace("dataset-", "other-")
        copied["block"] = "other-t1-r1"
        copied["spec"]["metadata"]["dataset"] = "other"
        second_runs.append(copied)
    manifest["runs"].extend(second_runs)
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False))
    schedule = build_schedule(manifest, manifest_path)
    output = tmp_path / "results"

    assert execute_schedule(
        manifest,
        manifest_path,
        schedule,
        output,
        datasets=frozenset({"other"}),
    ) == (2, 0, 0)
    events = [
        json.loads(line)
        for line in (output / "execution-log.jsonl").read_text().splitlines()
        if json.loads(line)["event"] == "start"
    ]
    observed = [event["condition_id"] for event in events]
    expected = [
        entry["condition_id"]
        for entry in schedule["entries"]
        if entry["condition_id"].startswith("other-")
    ]
    assert observed == expected
    rows, exclusions = collect_rows(manifest, schedule, output, frozenset({"other"}))
    assert len(rows) == 2
    assert exclusions == []


def test_unknown_dataset_block_is_rejected(tmp_path):
    manifest, manifest_path = make_manifest(tmp_path)
    schedule = build_schedule(manifest, manifest_path)

    with pytest.raises(ScheduleError, match="unknown dataset block"):
        execute_schedule(
            manifest,
            manifest_path,
            schedule,
            tmp_path / "results",
            datasets=frozenset({"missing"}),
        )


def test_failed_attempt_is_preserved_and_retried(tmp_path):
    manifest, manifest_path = make_manifest(tmp_path, fail=True)
    schedule = build_schedule(manifest, manifest_path)
    output = tmp_path / "results"

    assert execute_schedule(manifest, manifest_path, schedule, output) == (2, 0, 1)
    assert execute_schedule(manifest, manifest_path, schedule, output) == (1, 1, 1)
    failed = next(
        entry for entry in schedule["entries"] if entry["condition_id"].endswith("beta")
    )
    attempts = sorted((output / failed["run_id"]).glob("attempt-*"))
    assert len(attempts) == 2
    assert all((attempt / "run.json").is_file() for attempt in attempts)


def test_cross_tool_difference_is_separate_from_within_tool_determinism():
    rows = []
    for tool, digest in (("alpha", "a" * 64), ("beta", "b" * 64)):
        for replicate in (1, 2):
            rows.append(
                {
                    "dataset": "dataset",
                    "tool": tool,
                    "execution_mode": "default",
                    "replicate": replicate,
                    "normalized_output_sha256": digest,
                    "emitted_records": 10,
                }
            )
    report = correctness(rows)

    assert report["all_deterministic"] is True
    assert report["all_identical"] is False
    assert report["datasets"][0]["all_tools_modes_deterministic"] is True


def test_concurrent_coordinator_is_refused(tmp_path):
    manifest, manifest_path = make_manifest(tmp_path)
    schedule = build_schedule(manifest, manifest_path)
    output = tmp_path / "results"

    with (
        execution_lock(output),
        pytest.raises(ScheduleError, match="another coordinator"),
    ):
        execute_schedule(manifest, manifest_path, schedule, output)


def test_incomplete_manifest_is_refused(tmp_path):
    manifest, manifest_path = make_manifest(tmp_path)
    manifest["study"]["name"] = "REQUIRED"
    with pytest.raises(ScheduleError, match="incomplete manifest"):
        validate_manifest(manifest, manifest_path)


def test_publication_manifest_requires_pinned_inputs(tmp_path):
    manifest, manifest_path = make_manifest(tmp_path)
    manifest["mode"] = "publication"
    with pytest.raises(ScheduleError, match="must declare path and sha256"):
        validate_manifest(manifest, manifest_path)


def test_publication_manifest_requires_every_runtime_file_in_artifacts(tmp_path):
    manifest, manifest_path = make_manifest(tmp_path)
    manifest["mode"] = "publication"
    input_path = tmp_path / "input.txt"
    for run in manifest["runs"]:
        run["spec"]["inputs"] = [
            {
                "path": str(input_path),
                "sha256": sha256_file(input_path),
                "verify": False,
            }
        ]
    with pytest.raises(ScheduleError, match="absent from artifacts"):
        validate_manifest(manifest, manifest_path)
