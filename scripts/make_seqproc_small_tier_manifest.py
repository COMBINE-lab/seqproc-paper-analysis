#!/usr/bin/env python3
"""Freeze a local seqproc-only integration tier into a scheduler manifest."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path
from typing import Any, Sequence

import yaml

from benchmark_harness import inspect_fastq, sha256_file


ROOT = Path(__file__).resolve().parents[1]
CONFIG_ROOT = ROOT / "configs" / "seqproc"
DATASETS = (
    {
        "name": "splitseq_pe",
        "r1": "splitseq_pe_R1.fastq",
        "r2": "splitseq_pe_R2.fastq",
        "geometry": "splitseq_filter_edit.geom",
        "additional": (
            "splitseq_bc3_seq2seq.tsv",
            "splitseq_bc2_seq2seq.tsv",
            "splitseq_bc1_seq2seq.tsv",
        ),
    },
    {
        "name": "lr_splitseq",
        "r1": "splitseq_se_R1.fastq",
        "r2": None,
        "geometry": "splitseq_singleend_edit_ann.geom",
        "additional": (
            "splitseq_bc3_seq2seq.tsv",
            "splitseq_bc2_seq2seq.tsv",
            "splitseq_bc1_seq2seq.tsv",
        ),
    },
    {
        "name": "tenx_v2",
        "r1": "10x_short_R1.fastq",
        "r2": "10x_short_R2.fastq",
        "geometry": "10x_v2.geom",
        "additional": (),
    },
    {
        "name": "scirnaseq3",
        "r1": "sciseq_R1.fastq",
        "r2": "sciseq_R2.fastq",
        "geometry": "sciseq3_edit.geom",
        "additional": (),
    },
)


def git_commit(repository: Path) -> str:
    return subprocess.run(
        ["git", "-C", str(repository), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def frozen_file(path: Path, verify: bool = False) -> dict[str, Any]:
    return {"path": str(path.resolve()), "sha256": sha256_file(path), "verify": verify}


def build_manifest(
    binary: Path,
    input_dir: Path,
    seqproc_repo: Path,
    antisequence_repo: Path,
    threads: list[int],
    replicates: int,
    seed: int,
) -> dict[str, Any]:
    binary = binary.resolve()
    input_dir = input_dir.resolve()
    seqproc_repo = seqproc_repo.resolve()
    antisequence_repo = antisequence_repo.resolve()
    required_paths = [binary]
    for dataset in DATASETS:
        required_paths.append(input_dir / str(dataset["r1"]))
        if dataset["r2"]:
            required_paths.append(input_dir / str(dataset["r2"]))
        required_paths.append(CONFIG_ROOT / str(dataset["geometry"]))
        required_paths.extend(CONFIG_ROOT / name for name in dataset["additional"])
    required_paths.extend(
        (
            Path(__file__).resolve(),
            ROOT / "scripts" / "benchmark_harness.py",
            ROOT / "scripts" / "run_frozen_schedule.py",
            ROOT / "scripts" / "summarize_frozen_schedule.py",
            ROOT / "requirements.txt",
            ROOT / "requirements.lock",
        )
    )
    missing = sorted(str(path) for path in required_paths if not path.is_file())
    if missing:
        raise FileNotFoundError("missing required files:\n" + "\n".join(missing))

    artifacts_by_path = {
        str(path.resolve()): {"path": str(path.resolve()), "sha256": sha256_file(path)}
        for path in required_paths
    }
    repositories = [
        {
            "name": "seqproc",
            "path": str(seqproc_repo),
            "commit": git_commit(seqproc_repo),
        },
        {
            "name": "antisequence",
            "path": str(antisequence_repo),
            "commit": git_commit(antisequence_repo),
        },
    ]
    dataset_records: dict[str, int] = {}
    dataset_manifest: list[dict[str, Any]] = []
    for dataset in DATASETS:
        r1 = input_dir / str(dataset["r1"])
        r2 = input_dir / str(dataset["r2"]) if dataset["r2"] else None
        records = int(inspect_fastq(r1)["records"])
        if r2 is not None:
            mate_records = int(inspect_fastq(r2)["records"])
            if mate_records != records:
                raise ValueError(
                    f"paired input count mismatch for {dataset['name']}: {records} != {mate_records}"
                )
        dataset_records[str(dataset["name"])] = records
        dataset_manifest.append(
            {
                "name": dataset["name"],
                "records": records,
                "r1": frozen_file(r1),
                "r2": frozen_file(r2) if r2 is not None else None,
                "mates_verified": r2 is None or mate_records == records,
            }
        )
    runs: list[dict[str, Any]] = []
    for replicate in range(1, replicates + 1):
        for dataset in DATASETS:
            r1 = input_dir / str(dataset["r1"])
            r2 = input_dir / str(dataset["r2"]) if dataset["r2"] else None
            geometry = CONFIG_ROOT / str(dataset["geometry"])
            additional = [CONFIG_ROOT / name for name in dataset["additional"]]
            for thread_count in threads:
                block = f"{dataset['name']}-t{thread_count}-r{replicate}"
                for execution_mode in ("default", "staged"):
                    command = [
                        str(binary),
                        "run",
                        "--geom",
                        str(geometry),
                        "--file1",
                        str(r1),
                        "--out1",
                        "{run_dir}/R1.fastq",
                        "--threads",
                        str(thread_count),
                    ]
                    outputs: list[dict[str, Any]] = [
                        {
                            "path": "{run_dir}/R1.fastq",
                            "format": "fastq",
                            "normalize": "fastq_multiset",
                            "mate": 1,
                            "min_bytes": 1,
                        }
                    ]
                    inputs = [frozen_file(r1)]
                    if r2 is not None:
                        command.extend(("--file2", str(r2), "--out2", "{run_dir}/R2.fastq"))
                        outputs.append(
                            {
                                "path": "{run_dir}/R2.fastq",
                                "format": "fastq",
                                "normalize": "fastq_multiset",
                                "mate": 2,
                                "min_bytes": 1,
                            }
                        )
                        inputs.append(frozen_file(r2))
                    if execution_mode == "staged":
                        command.append("--staged-pipeline")
                    for path in additional:
                        command.extend(("--additional", str(path)))
                    configs = [frozen_file(geometry), *(frozen_file(path) for path in additional)]
                    condition_id = f"{block}-{execution_mode}"
                    runs.append(
                        {
                            "id": condition_id,
                            "block": block,
                            "spec": {
                                "name": "seqproc-small-tier",
                                "cwd": str(ROOT),
                                "command": command,
                                "inputs": inputs,
                                "configs": configs,
                                "outputs": outputs,
                                "repositories": repositories,
                                "metadata": {
                                    "dataset": dataset["name"],
                                    "tool": "seqproc",
                                    "execution_mode": execution_mode,
                                    "threads": thread_count,
                                    "replicate": replicate,
                                    "input_records": dataset_records[str(dataset["name"])],
                                    "tier": "integration-small",
                                },
                            },
                        }
                    )
    return {
        "schema_version": "1.0.0",
        "mode": "development",
        "study": {
            "name": "seqproc-current-small-tier",
            "random_seed": seed,
            "purpose": "exercise frozen publication harness before the full cross-tool rerun",
        },
        "artifacts": sorted(artifacts_by_path.values(), key=lambda item: item["path"]),
        "datasets": dataset_manifest,
        "execution": {
            "timeout_seconds": 600,
            "sanitized_environment_allowlist": [
                "PATH",
                "LD_LIBRARY_PATH",
                "LANG",
                "TMPDIR",
                "RUST_BACKTRACE",
                "RUST_LOG",
            ],
        },
        "runs": runs,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--seqproc-repo", type=Path, default=ROOT.parent / "seqproc")
    parser.add_argument("--antisequence-repo", type=Path, default=ROOT.parent / "antisequence")
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 8, 32])
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--seed", type=int, default=741211)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.replicates <= 0 or any(value <= 0 for value in args.threads):
        parser.error("replicates and thread counts must be positive")
    manifest = build_manifest(
        args.binary,
        args.input_dir,
        args.seqproc_repo,
        args.antisequence_repo,
        sorted(set(args.threads)),
        args.replicates,
        args.seed,
    )
    if args.output.exists():
        parser.error(f"refusing to overwrite {args.output}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(yaml.safe_dump(manifest, sort_keys=False))
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
