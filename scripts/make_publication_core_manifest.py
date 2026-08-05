#!/usr/bin/env python3
"""Freeze the full-data, cross-tool publication benchmark manifest."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path
from typing import Any, Sequence

import yaml

from benchmark_harness import sha256_file


ROOT = Path(__file__).resolve().parents[1]
ECOSYSTEM = ROOT.parent
TASKSET = Path("/usr/bin/taskset")
DATASETS = (
    {
        "name": "splitseq_pe",
        "r1": "SRR6750041_R1.fastq",
        "r2": "SRR6750041_R2.fastq",
        "seqproc": "publication_splitseq_pe.geom",
        "matchbox": "publication_splitseq_pe.mb",
        "splitcode": "publication_splitseq_pe.config",
        "seqproc_support": (
            "splitseq_bc23_whitelist.txt",
            "splitseq_bc1_whitelist_6bp.txt",
            "splitseq_bc3_seq2seq.tsv",
            "splitseq_bc2_seq2seq.tsv",
            "splitseq_bc1_seq2seq.tsv",
        ),
        "matchbox_support": ("rt_6bp.csv", "r2_r3.txt"),
        "splitcode_assign": True,
        "equivalence": "best-practical; splitcode linker matching is substitution-only",
    },
    {
        "name": "lr_splitseq",
        "r1": "SRR13948564_full.fastq",
        "r2": None,
        "seqproc": "publication_lr_splitseq.geom",
        "matchbox": "publication_lr_splitseq.mb",
        "splitcode": "publication_lr_splitseq.config",
        "seqproc_support": (),
        "matchbox_support": (),
        "splitcode_assign": True,
        "equivalence": "forward-orientation best-practical; splitcode linker matching is substitution-only",
    },
    {
        "name": "tenx_v2",
        "r1": "SRR8315379_R1.fastq",
        "r2": "SRR8315379_R2.fastq",
        "seqproc": "10x_v2.geom",
        "matchbox": "publication_10x_v2.mb",
        "splitcode": "10x_v2.config",
        "seqproc_support": (),
        "matchbox_support": (),
        "splitcode_assign": False,
        "equivalence": "fixed-position exact extraction",
    },
    {
        "name": "scirnaseq3",
        "r1": "SRR7827254_1.fastq",
        "r2": "SRR7827254_2.fastq",
        "seqproc": "sciseq3_edit.geom",
        "matchbox": "publication_sciseq3.mb",
        "splitcode": "sciseq3.config",
        "seqproc_support": (),
        "matchbox_support": (),
        "splitcode_assign": True,
        "equivalence": "best-practical; indel semantics differ by backend",
    },
)


def git_commit(repository: Path) -> str:
    return subprocess.run(
        ["git", "-C", str(repository), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def frozen(path: Path, digest: str | None = None) -> dict[str, Any]:
    path = path.resolve()
    return {
        "path": str(path),
        "sha256": digest if digest is not None else sha256_file(path),
        "verify": False,
    }


def build_manifest(args: argparse.Namespace) -> dict[str, Any]:
    provenance_path = args.data_provenance.resolve()
    provenance = json.loads(provenance_path.read_text())
    if provenance.get("complete") is not True:
        raise ValueError(f"full-data provenance is incomplete: {provenance_path}")
    by_name = {
        Path(item["fastq"]["path"]).name: item["fastq"] for item in provenance["files"]
    }
    binaries = {
        "seqproc": args.seqproc_binary.resolve(),
        "matchbox": args.matchbox_binary.resolve(),
        "splitcode": args.splitcode_binary.resolve(),
    }
    repositories = {
        "analysis": ROOT,
        "seqproc": args.seqproc_repo.resolve(),
        "antisequence": args.antisequence_repo.resolve(),
        "matchbox": args.matchbox_repo.resolve(),
        "splitcode": args.splitcode_repo.resolve(),
    }
    repository_records = {
        name: {"name": name, "path": str(path), "commit": git_commit(path)}
        for name, path in repositories.items()
    }

    artifacts: dict[str, dict[str, str]] = {}

    def add_artifact(path: Path, digest: str | None = None) -> None:
        resolved = path.resolve()
        if not resolved.is_file():
            raise FileNotFoundError(resolved)
        value = digest if digest is not None else sha256_file(resolved)
        artifacts[str(resolved)] = {"path": str(resolved), "sha256": value}

    add_artifact(TASKSET)
    add_artifact(provenance_path)
    for binary in binaries.values():
        add_artifact(binary)
    for path in (
        Path(__file__),
        ROOT / "scripts" / "benchmark_harness.py",
        ROOT / "scripts" / "run_frozen_schedule.py",
        ROOT / "scripts" / "summarize_frozen_schedule.py",
        ROOT / "requirements.txt",
        ROOT / "requirements.lock",
    ):
        add_artifact(path)

    dataset_records: dict[str, int] = {}
    dataset_entries = []
    resolved_datasets: dict[str, dict[str, Any]] = {}
    for dataset in DATASETS:
        r1_record = by_name[str(dataset["r1"])]
        r2_record = by_name[str(dataset["r2"])] if dataset["r2"] else None
        if r2_record and int(r1_record["records"]) != int(r2_record["records"]):
            raise ValueError(f"paired record-count mismatch for {dataset['name']}")
        r1 = Path(r1_record["path"]).resolve()
        r2 = Path(r2_record["path"]).resolve() if r2_record else None
        add_artifact(r1, str(r1_record["sha256"]))
        if r2 is not None:
            add_artifact(r2, str(r2_record["sha256"]))
        dataset_records[str(dataset["name"])] = int(r1_record["records"])
        resolved_datasets[str(dataset["name"])] = {"r1": r1, "r2": r2}
        dataset_entries.append(
            {
                "name": dataset["name"],
                "records": int(r1_record["records"]),
                "r1": frozen(r1, str(r1_record["sha256"])),
                "r2": frozen(r2, str(r2_record["sha256"])) if r2 else None,
                "equivalence_scope": dataset["equivalence"],
            }
        )
        config_paths = (
            ROOT / "configs" / "seqproc" / str(dataset["seqproc"]),
            ROOT / "configs" / "matchbox" / str(dataset["matchbox"]),
            ROOT / "configs" / "splitcode" / str(dataset["splitcode"]),
            *(ROOT / "configs" / "seqproc" / name for name in dataset["seqproc_support"]),
            *(ROOT / "configs" / "matchbox" / name for name in dataset["matchbox_support"]),
        )
        for path in config_paths:
            add_artifact(path)

    runs: list[dict[str, Any]] = []
    for replicate in range(1, args.replicates + 1):
        for dataset in DATASETS:
            dataset_name = str(dataset["name"])
            r1 = resolved_datasets[dataset_name]["r1"]
            r2 = resolved_datasets[dataset_name]["r2"]
            inputs = [frozen(r1, by_name[r1.name]["sha256"])]
            if r2 is not None:
                inputs.append(frozen(r2, by_name[r2.name]["sha256"]))
            for threads in args.threads:
                block = f"{dataset_name}-t{threads}-r{replicate}"
                affinity = "0" if threads == 1 else f"0-{threads - 1}"
                for tool in ("seqproc", "matchbox", "splitcode"):
                    tool_binary = binaries[tool]
                    command = [str(TASKSET), "--cpu-list", affinity, str(tool_binary)]
                    configs: list[Path]
                    outputs: list[dict[str, Any]] = []
                    if tool == "seqproc":
                        geometry = ROOT / "configs" / "seqproc" / str(dataset["seqproc"])
                        support = [
                            ROOT / "configs" / "seqproc" / name
                            for name in dataset["seqproc_support"]
                        ]
                        configs = [geometry, *support]
                        command.extend(
                            [
                                "run",
                                "--geom",
                                str(geometry),
                                "--file1",
                                str(r1),
                                "--out1",
                                "{run_dir}/R1.fastq",
                                "--threads",
                                str(threads),
                                "--staged-pipeline",
                            ]
                        )
                        outputs.append(
                            {
                                "path": "{run_dir}/R1.fastq",
                                "format": "fastq",
                                "normalize": "fastq_id_multiset",
                                "mate": 1,
                                "min_bytes": 1,
                            }
                        )
                        if r2 is not None:
                            command.extend(
                                ["--file2", str(r2), "--out2", "{run_dir}/R2.fastq"]
                            )
                            outputs.append(
                                {
                                    "path": "{run_dir}/R2.fastq",
                                    "format": "fastq",
                                    "normalize": "fastq_id_multiset",
                                    "mate": 2,
                                    "min_bytes": 1,
                                }
                            )
                        for support_path in support:
                            if support_path.name.endswith("_seq2seq.tsv"):
                                command.extend(["--additional", str(support_path)])
                        repos = [repository_records[name] for name in ("analysis", "seqproc", "antisequence")]
                    elif tool == "matchbox":
                        script = ROOT / "configs" / "matchbox" / str(dataset["matchbox"])
                        support = [
                            ROOT / "configs" / "matchbox" / name
                            for name in dataset["matchbox_support"]
                        ]
                        configs = [script, *support]
                        command.extend(
                            ["--error", "0.2", "--match-mode", "one-best", "--threads", str(threads), "--script-file", str(script), str(r1)]
                        )
                        if r2 is not None:
                            command.extend(["--paired-with", str(r2)])
                        outputs.append(
                            {
                                "path": str(ROOT / "mb_r1.fq"),
                                "format": "fastq",
                                "normalize": "fastq_id_multiset",
                                "mate": 1,
                                "min_bytes": 1,
                            }
                        )
                        if r2 is not None:
                            outputs.append(
                                {
                                    "path": str(ROOT / "mb_r2.fq"),
                                    "format": "fastq",
                                    "normalize": "fastq_id_multiset",
                                    "mate": 2,
                                    "min_bytes": 1,
                                }
                            )
                        repos = [repository_records[name] for name in ("analysis", "matchbox")]
                    else:
                        config = ROOT / "configs" / "splitcode" / str(dataset["splitcode"])
                        configs = [config]
                        command.extend(["--config", str(config)])
                        if dataset["splitcode_assign"]:
                            command.extend(["--assign", "--mapping", "{run_dir}/mapping.txt"])
                        command.extend(["--nFastqs", "2" if r2 is not None else "1", "--threads", str(threads)])
                        output_paths = ["{run_dir}/R1.fastq"]
                        if r2 is not None:
                            output_paths.append("{run_dir}/R2.fastq")
                        command.extend(["--output", ",".join(output_paths), str(r1)])
                        if r2 is not None:
                            command.append(str(r2))
                        outputs.extend(
                            {
                                "path": path,
                                "format": "fastq",
                                "normalize": "fastq_id_multiset",
                                "mate": index,
                                "min_bytes": 1,
                            }
                            for index, path in enumerate(output_paths, start=1)
                        )
                        repos = [repository_records[name] for name in ("analysis", "splitcode")]

                    condition_id = f"{block}-{tool}"
                    runs.append(
                        {
                            "id": condition_id,
                            "block": block,
                            "spec": {
                                "name": "seqproc-publication-core",
                                "cwd": str(ROOT),
                                "command": command,
                                "executables": [frozen(tool_binary)],
                                "inputs": inputs,
                                "configs": [frozen(path) for path in configs],
                                "outputs": outputs,
                                "retain_outputs": False,
                                "repositories": repos,
                                "metadata": {
                                    "dataset": dataset_name,
                                    "tool": tool,
                                    "execution_mode": "best-practical-uncompressed",
                                    "threads": threads,
                                    "replicate": replicate,
                                    "input_records": dataset_records[dataset_name],
                                    "cpu_affinity": affinity,
                                    "equivalence_scope": dataset["equivalence"],
                                    "tier": "publication-core-full-data",
                                },
                            },
                        }
                    )

    return {
        "schema_version": "1.0.0",
        "mode": "publication",
        "study": {
            "name": "seqproc-publication-core-full-data",
            "random_seed": args.seed,
            "purpose": "full-data cross-tool scaling and resource benchmark",
            "require_cross_tool_identity": False,
        },
        "artifacts": sorted(artifacts.values(), key=lambda item: item["path"]),
        "data_provenance": frozen(provenance_path),
        "datasets": dataset_entries,
        "execution": {
            "timeout_seconds": args.timeout_seconds,
            "sanitized_environment_allowlist": ["PATH", "LD_LIBRARY_PATH", "LANG", "TMPDIR"],
            "cpu_policy": "requested threads pinned to physical CPUs 0 through N-1; SMT siblings excluded",
        },
        "runs": runs,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-provenance", type=Path, required=True)
    parser.add_argument("--seqproc-binary", type=Path, default=ECOSYSTEM / "seqproc" / "target" / "release" / "seqproc")
    parser.add_argument("--matchbox-binary", type=Path, default=ECOSYSTEM / "competitors" / "matchbox-v0.3.2" / "target" / "release" / "matchbox")
    parser.add_argument("--splitcode-binary", type=Path, default=ECOSYSTEM / "competitors" / "splitcode-v0.31.6" / "build" / "src" / "splitcode")
    parser.add_argument("--seqproc-repo", type=Path, default=ECOSYSTEM / "seqproc")
    parser.add_argument("--antisequence-repo", type=Path, default=ECOSYSTEM / "antisequence")
    parser.add_argument("--matchbox-repo", type=Path, default=ECOSYSTEM / "competitors" / "matchbox-v0.3.2")
    parser.add_argument("--splitcode-repo", type=Path, default=ECOSYSTEM / "competitors" / "splitcode-v0.31.6")
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 2, 4, 8, 16, 32])
    parser.add_argument("--replicates", type=int, default=5)
    parser.add_argument("--seed", type=int, default=741211)
    parser.add_argument("--timeout-seconds", type=int, default=21600)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.replicates <= 0 or args.timeout_seconds <= 0 or any(value <= 0 or value > 64 for value in args.threads):
        parser.error("replicates/timeout must be positive and threads must be in 1..64")
    if args.output.exists():
        parser.error(f"refusing to overwrite {args.output}")
    manifest = build_manifest(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(yaml.safe_dump(manifest, sort_keys=False))
    print(f"wrote {args.output} with {len(manifest['runs'])} conditions")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
