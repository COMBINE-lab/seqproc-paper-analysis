import gzip
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from benchmark_harness import (
    HarnessError,
    execute_spec,
    inspect_fastq,
    inspect_and_normalize_fastq,
    normalized_fastq_id_multiset_sha256,
    normalized_fastq_multiset_sha256,
    prepare_spec,
    sha256_file,
)


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


def test_launcher_and_tool_executables_are_both_fingerprinted(tmp_path):
    spec = make_spec(tmp_path, "pass")
    spec["executables"] = ["/usr/bin/time"]
    identity = prepare_spec(spec)["identity"]

    assert identity["binary"]["path"] == str(Path(sys.executable).resolve())
    assert [item["path"] for item in identity["executables"]] == [
        str(Path(sys.executable).resolve()),
        str(Path("/usr/bin/time").resolve()),
    ]


def test_validated_output_can_be_removed_after_recording_digests(tmp_path):
    spec = make_spec(
        tmp_path,
        "from pathlib import Path; import sys; "
        "Path(sys.argv[1], 'result.fastq').write_text('@r1\\nAC\\n+\\nII\\n')",
        [{"path": "{run_dir}/result.fastq", "format": "fastq", "min_bytes": 1}],
    )
    spec["retain_outputs"] = False
    result, run_dir = execute_spec(spec, tmp_path / "runs")

    assert result["success"] is True
    assert result["outputs"][0]["retained"] is False
    assert result["outputs"][0]["sha256"]
    assert not (run_dir / "result.fastq").exists()
    assert (run_dir / "outputs.sha256").read_text().strip()


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


def test_preexisting_declared_output_is_refused(tmp_path):
    stale = tmp_path / "stale.fastq"
    stale.write_text("@old\nA\n+\nI\n")
    spec = make_spec(tmp_path, "pass", [str(stale)])

    with pytest.raises(HarnessError, match="pre-existing declared outputs"):
        execute_spec(spec, tmp_path / "runs")


def test_fastq_output_is_counted_and_validated(tmp_path):
    spec = make_spec(
        tmp_path,
        "from pathlib import Path; import sys; "
        "Path(sys.argv[1], 'result.fastq').write_text('@r1\\nACGT\\n+\\nIIII\\n')",
        [
            {
                "path": "{run_dir}/result.fastq",
                "format": "fastq",
                "normalize": "fastq_multiset",
                "mate": 1,
                "min_bytes": 1,
            }
        ],
    )
    result, run_dir = execute_spec(spec, tmp_path / "runs")

    assert result["success"] is True
    assert result["output_counts"][0]["records"] == 1
    assert json.loads((run_dir / "output-counts.json").read_text())[0]["valid"] is True
    assert inspect_fastq(run_dir / "result.fastq")["records"] == 1
    assert result["outputs"][0]["normalized_sha256"]
    assert (run_dir / "outputs.normalized.sha256").read_text().strip()


def test_normalized_fastq_digest_is_order_independent_and_gzip_independent(tmp_path):
    first = b"@a/1 extra\nAC\n+\nII\n@b/1\nGT\n+\nJJ\n"
    reordered = b"@b/1\nGT\n+\nJJ\n@a/1 extra\nAC\n+\nII\n"
    path1 = tmp_path / "first.fastq"
    path2 = tmp_path / "second.fastq"
    gzip_path = tmp_path / "second.fastq.gz"
    path1.write_bytes(first)
    path2.write_bytes(reordered)
    with gzip.open(gzip_path, "wb") as handle:
        handle.write(reordered)

    digest = normalized_fastq_multiset_sha256(path1, mate=1, chunk_bytes=10)
    assert digest == normalized_fastq_multiset_sha256(path2, mate=1, chunk_bytes=10)
    assert digest == normalized_fastq_multiset_sha256(gzip_path, mate=1, chunk_bytes=10)
    assert digest != normalized_fastq_multiset_sha256(path2, mate=2, chunk_bytes=10)


def test_read_id_digest_ignores_sequence_and_order_but_not_mate(tmp_path):
    first = tmp_path / "first.fastq"
    second = tmp_path / "second.fastq"
    first.write_bytes(b"@a/1 extra\nAC\n+\nII\n@b/1\nGT\n+\nJJ\n")
    second.write_bytes(b"@b/1\nTT\n+\nHH\n@a/1 other\nCC\n+\nGG\n")

    digest = normalized_fastq_id_multiset_sha256(first, mate=1, chunk_bytes=5)
    assert digest == normalized_fastq_id_multiset_sha256(second, mate=1, chunk_bytes=5)
    assert digest != normalized_fastq_id_multiset_sha256(second, mate=2, chunk_bytes=5)


def test_fused_fastq_inspection_matches_independent_digests(tmp_path):
    content = b"@b/1\nGT\n+\nJJ\n@a/1 extra\nAC\n+\nII\n"
    for suffix in (".fastq", ".fastq.gz"):
        path = tmp_path / f"reads{suffix}"
        if suffix.endswith(".gz"):
            with gzip.open(path, "wb") as handle:
                handle.write(content)
        else:
            path.write_bytes(content)
        result = inspect_and_normalize_fastq(
            path,
            mate=1,
            normalization="fastq_id_multiset",
            chunk_bytes=5,
            temp_dir=tmp_path,
        )
        assert result["records"] == 2
        assert result["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
        assert result["normalized_sha256"] == normalized_fastq_id_multiset_sha256(
            path, mate=1, chunk_bytes=5, temp_dir=tmp_path
        )


def test_numeric_accession_set_is_order_independent_and_exact(tmp_path):
    first = tmp_path / "first.fastq"
    reordered = tmp_path / "reordered.fastq"
    first.write_bytes(b"@SRR1.1/1 extra\nAC\n+\nII\n@SRR1.3\nGT\n+\nJJ\n")
    reordered.write_bytes(b"@SRR1.3\nTT\n+\nHH\n@SRR1.1 other\nCC\n+\nGG\n")
    left = inspect_and_normalize_fastq(
        first,
        mate=1,
        normalization="fastq_numeric_accession_set",
        numeric_id_max=3,
    )
    right = inspect_and_normalize_fastq(
        reordered,
        mate=1,
        normalization="fastq_numeric_accession_set",
        numeric_id_max=3,
    )
    assert left["normalized_sha256"] == right["normalized_sha256"]
    assert left["sha256"] != right["sha256"]


def test_compiled_numeric_auditor_matches_python(tmp_path):
    executable = Path(__file__).resolve().parents[1] / "tools" / "bin" / "fastq-numeric-audit"
    if not executable.is_file():
        pytest.skip("compiled numeric FASTQ auditor is not available")
    path = tmp_path / "reads.fastq"
    path.write_bytes(b"@SRR1.3\nTT\n+\nHH\n@SRR1.1 other\nCC\n+\nGG\n")
    expected = inspect_and_normalize_fastq(
        path,
        mate=1,
        normalization="fastq_numeric_accession_set",
        numeric_id_max=3,
    )
    observed = inspect_and_normalize_fastq(
        path,
        mate=1,
        normalization="fastq_numeric_accession_set",
        numeric_id_max=3,
        numeric_audit_executable=executable,
        temp_dir=tmp_path,
    )
    assert observed["records"] == expected["records"]
    assert observed["sha256"] == expected["sha256"]
    assert observed["normalized_sha256"] == expected["normalized_sha256"]
    assert observed["validator"] == "fastq-numeric-audit-v1"


@pytest.mark.parametrize(
    "content, message",
    [
        (b"@not-numeric\nAC\n+\nII\n", "numeric accession ID"),
        (b"@SRR1.4\nAC\n+\nII\n", "outside 1..3"),
        (b"@SRR1.1\nAC\n+\nII\n@SRR1.1\nGT\n+\nJJ\n", "duplicate"),
        (b"@SRR1.1\nAC\n+\nII\n@SRR2.2\nGT\n+\nJJ\n", "changes accession prefix"),
    ],
)
def test_numeric_accession_set_rejects_invalid_ids(tmp_path, content, message):
    path = tmp_path / "reads.fastq"
    path.write_bytes(content)
    with pytest.raises(HarnessError, match=message):
        inspect_and_normalize_fastq(
            path,
            normalization="fastq_numeric_accession_set",
            numeric_id_max=3,
        )


def test_malformed_or_too_small_output_fails_run(tmp_path):
    spec = make_spec(
        tmp_path,
        "from pathlib import Path; import sys; "
        "Path(sys.argv[1], 'bad.fastq').write_text('@r1\\nACGT\\n+\\nIII\\n')",
        [{"path": "{run_dir}/bad.fastq", "format": "fastq", "min_bytes": 100}],
    )
    result, _ = execute_spec(spec, tmp_path / "runs")

    assert result["success"] is False
    assert len(result["invalid_outputs"]) == 2


def test_timeout_preserves_failed_attempt(tmp_path):
    spec = make_spec(tmp_path, "import time; time.sleep(30)")
    spec["timeout_seconds"] = 0.05
    result, run_dir = execute_spec(spec, tmp_path / "runs")

    assert result["success"] is False
    assert result["timed_out"] is True
    assert result["termination_signal"] is not None
    assert json.loads((run_dir / "run.json").read_text())["timed_out"] is True


def test_identical_spec_has_stable_run_id_and_new_attempt(tmp_path):
    spec = make_spec(tmp_path, "pass")
    first, first_dir = execute_spec(spec, tmp_path / "runs")
    second, second_dir = execute_spec(spec, tmp_path / "runs")

    assert first["run_id"] == second["run_id"]
    assert first["attempt"] == 1
    assert second["attempt"] == 2
    assert first_dir != second_dir


def test_output_contract_and_timeout_are_part_of_run_identity(tmp_path):
    spec = make_spec(tmp_path, "pass", ["{run_dir}/result.txt"])
    original = prepare_spec(spec)["run_id"]
    spec["outputs"] = [{"path": "{run_dir}/result.txt", "min_bytes": 10}]
    assert prepare_spec(spec)["run_id"] != original
    with_output_contract = prepare_spec(spec)["run_id"]
    spec["timeout_seconds"] = 10
    assert prepare_spec(spec)["run_id"] != with_output_contract


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


def test_repository_commit_must_match_when_declared(tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    subprocess.run(["git", "-C", str(repo), "config", "user.email", "test@example.org"], check=True)
    subprocess.run(["git", "-C", str(repo), "config", "user.name", "Test"], check=True)
    (repo / "tracked.txt").write_text("tracked\n")
    subprocess.run(["git", "-C", str(repo), "add", "tracked.txt"], check=True)
    subprocess.run(["git", "-C", str(repo), "commit", "-qm", "baseline"], check=True)

    spec = make_spec(tmp_path, "pass")
    spec["repositories"] = [
        {"name": "repo", "path": str(repo), "commit": "0" * 40}
    ]
    with pytest.raises(HarnessError, match="expected"):
        prepare_spec(spec)


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
