import json
import os
import subprocess
import sys
from pathlib import Path

import pytest


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from benchmark_harness import HarnessError, execute_spec, prepare_spec, sha256_file


def make_spec(tmp_path: Path, code: str, outputs=None):
    input_path = tmp_path / "input.txt"
    input_path.write_text("input\n")
    return {
        "name": "harness-test",
        "cwd": str(tmp_path),
        "command": [sys.executable, "-c", code, "{run_dir}"],
        "inputs": [str(input_path)],
        "outputs": outputs or [],
        "environment_mode": "allowlist",
        "metadata": {"replicate": 1, "threads": 1},
    }


def test_success_records_timing_output_and_digest(tmp_path):
    spec = make_spec(
        tmp_path,
        "from pathlib import Path; import sys; "
        "Path(sys.argv[1], 'result.txt').write_text('result\\n'); print('ok')",
        ["{run_dir}/result.txt"],
    )
    result, run_dir = execute_spec(spec, tmp_path / "runs")

    assert result["success"] is True
    assert result["exit_code"] == 0
    assert result["outputs"][0]["sha256"] == sha256_file(run_dir / "result.txt")
    assert result["wall_seconds"] >= 0
    assert result["peak_rss_kib"] > 0
    assert (run_dir / "command.json").is_file()
    assert (run_dir / "identity.json").is_file()
    assert json.loads((run_dir / "run.json").read_text())["run_id"] == result["run_id"]


def test_nonzero_exit_is_preserved(tmp_path):
    spec = make_spec(tmp_path, "raise SystemExit(7)")
    result, run_dir = execute_spec(spec, tmp_path / "runs")

    assert result["success"] is False
    assert result["exit_code"] == 7
    assert result["time_exit_status"] == 7
    assert (run_dir / "run.json").is_file()


def test_missing_declared_output_fails_run(tmp_path):
    spec = make_spec(tmp_path, "pass", ["{run_dir}/missing.txt"])
    result, _ = execute_spec(spec, tmp_path / "runs")

    assert result["success"] is False
    assert result["exit_code"] == 0
    assert result["missing_outputs"]


def test_identical_spec_has_stable_run_id_and_new_attempt(tmp_path):
    spec = make_spec(tmp_path, "pass")
    first, first_dir = execute_spec(spec, tmp_path / "runs")
    second, second_dir = execute_spec(spec, tmp_path / "runs")

    assert first["run_id"] == second["run_id"]
    assert first["attempt"] == 1
    assert second["attempt"] == 2
    assert first_dir != second_dir


def test_argv_is_not_interpreted_by_a_shell(tmp_path):
    marker = tmp_path / "shell-expanded"
    token = f";touch {marker}"
    spec = make_spec(
        tmp_path,
        "import sys; assert sys.argv[2].startswith(';touch')",
    )
    spec["command"].append(token)
    result, _ = execute_spec(spec, tmp_path / "runs")

    assert result["success"] is True
    assert not marker.exists()


def test_declared_input_digest_can_avoid_rehash(tmp_path, monkeypatch):
    spec = make_spec(tmp_path, "pass")
    spec["inputs"] = [{"path": str(tmp_path / "input.txt"), "sha256": "abc", "verify": False}]
    prepared = prepare_spec(spec)

    assert prepared["identity"]["inputs"][0]["sha256"] == "abc"
    assert prepared["identity"]["inputs"][0]["digest_source"] == "declared"


def test_bad_declared_digest_is_rejected_when_verified(tmp_path):
    spec = make_spec(tmp_path, "pass")
    spec["inputs"] = [{"path": str(tmp_path / "input.txt"), "sha256": "wrong", "verify": True}]

    with pytest.raises(HarnessError, match="SHA-256 mismatch"):
        prepare_spec(spec)


def test_dirty_repository_requires_explicit_opt_in(tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    subprocess.run(["git", "-C", str(repo), "config", "user.email", "test@example.org"], check=True)
    subprocess.run(["git", "-C", str(repo), "config", "user.name", "Test"], check=True)
    tracked = repo / "tracked.txt"
    tracked.write_text("before\n")
    subprocess.run(["git", "-C", str(repo), "add", "tracked.txt"], check=True)
    subprocess.run(["git", "-C", str(repo), "commit", "-qm", "baseline"], check=True)
    tracked.write_text("after\n")

    spec = make_spec(tmp_path, "pass")
    spec["repositories"] = [{"name": "repo", "path": str(repo)}]
    with pytest.raises(HarnessError, match="is dirty"):
        prepare_spec(spec)

    spec["repositories"][0]["allow_dirty"] = True
    result, run_dir = execute_spec(spec, tmp_path / "runs")
    assert result["identity"]["repositories"][0]["dirty"] is True
    assert (run_dir / "repository-repo.diff").is_file()


def test_dirty_repository_archives_untracked_files(tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    subprocess.run(["git", "-C", str(repo), "config", "user.email", "test@example.org"], check=True)
    subprocess.run(["git", "-C", str(repo), "config", "user.name", "Test"], check=True)
    tracked = repo / "tracked.txt"
    tracked.write_text("tracked\n")
    subprocess.run(["git", "-C", str(repo), "add", "tracked.txt"], check=True)
    subprocess.run(["git", "-C", str(repo), "commit", "-qm", "baseline"], check=True)
    untracked = repo / "nested" / "candidate.txt"
    untracked.parent.mkdir()
    untracked.write_text("candidate source\n")

    spec = make_spec(tmp_path, "pass")
    spec["repositories"] = [
        {"name": "repo", "path": str(repo), "allow_dirty": True}
    ]
    result, run_dir = execute_spec(spec, tmp_path / "runs")

    repository = result["identity"]["repositories"][0]
    assert repository["untracked"] == [
        {
            "path": "nested/candidate.txt",
            "bytes": len("candidate source\n"),
            "sha256": sha256_file(untracked),
        }
    ]
    archived = run_dir / "repository-repo-untracked" / "nested" / "candidate.txt"
    assert archived.read_text() == "candidate source\n"
