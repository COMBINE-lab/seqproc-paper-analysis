#!/usr/bin/env python3
"""Freeze the full-data, cross-tool publication benchmark manifest."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import yaml
from benchmark_harness import sha256_file

ROOT = Path(__file__).resolve().parents[1]
ECOSYSTEM = ROOT.parent
TASKSET = Path("/usr/bin/taskset")
DEFAULT_FASTQ_AUDITOR = ROOT / "tools" / "bin" / "fastq-numeric-audit"
DEFAULT_LR_REVERSE_COMPLEMENT = Path(
    "/scratch1/seqproc-benchmark-data/full/fastq/SRR13948564_full_RC.fastq"
)
DEFAULT_LR_REVERSE_COMPLEMENT_PROVENANCE = Path(
    "/scratch1/seqproc-benchmark-data/full/fastq/SRR13948564_full_RC.provenance.json"
)
SPLITCODE_DUAL_WRAPPER = ROOT / "scripts" / "run_splitcode_dual_pass.py"
DEFAULT_THREADS = (1, 4, 16, 32)
DEFAULT_REPLICATES = 3
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
        "splitcode_dual_pass": False,
        "analysis_role": "primary",
        "equivalence": "best-practical; splitcode linker matching is substitution-only",
    },
    {
        "name": "lr_splitseq_dual",
        "r1": "SRR13948564_full.fastq",
        "r2": None,
        "r1_reverse_complement": True,
        "seqproc": "splitseq_singleend_edit_ann.geom",
        "matchbox": "publication_lr_splitseq_dual.mb",
        "splitcode": "splitseq_singleend.config",
        "seqproc_support": (),
        "matchbox_support": (),
        "splitcode_assign": True,
        "splitcode_dual_pass": True,
        "analysis_role": "primary",
        "equivalence": (
            "capability-complete best-practical dual orientation; seqproc and matchbox "
            "process both orientations natively; splitcode is measured as two sequential "
            "passes over forward and precomputed reverse-complement inputs without "
            "duplicate reconciliation"
        ),
    },
    {
        "name": "lr_splitseq_forward",
        "r1": "SRR13948564_full.fastq",
        "r2": None,
        "r1_reverse_complement": False,
        "seqproc": "publication_lr_splitseq.geom",
        "matchbox": "publication_lr_splitseq.mb",
        "splitcode": "publication_lr_splitseq.config",
        "seqproc_support": (),
        "matchbox_support": (),
        "splitcode_assign": True,
        "splitcode_dual_pass": False,
        "analysis_role": "supplementary",
        "equivalence": (
            "forward-orientation controlled supplementary comparison; splitcode linker "
            "matching is substitution-only"
        ),
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
        "splitcode_dual_pass": False,
        "analysis_role": "primary",
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
        "splitcode_dual_pass": False,
        "analysis_role": "primary",
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
    lr_source = by_name["SRR13948564_full.fastq"]
    lr_reverse_complement = args.lr_reverse_complement.resolve()
    lr_reverse_provenance_path = args.lr_reverse_complement_provenance.resolve()
    if not lr_reverse_complement.is_file() or not lr_reverse_provenance_path.is_file():
        raise FileNotFoundError(
            "the precomputed LR reverse-complement FASTQ and its provenance are required"
        )
    lr_reverse_provenance = json.loads(lr_reverse_provenance_path.read_text())
    if (
        lr_reverse_provenance.get("transformation")
        != "reverse-complement FASTQ sequence and reverse quality"
        or lr_reverse_provenance.get("header_preserved") is not True
        or lr_reverse_provenance.get("reconciliation_performed") is not False
        or int(lr_reverse_provenance.get("records", -1)) != int(lr_source["records"])
        or Path(lr_reverse_provenance["source"]["path"]).resolve()
        != Path(lr_source["path"]).resolve()
        or str(lr_reverse_provenance["source"]["sha256"]) != str(lr_source["sha256"])
        or Path(lr_reverse_provenance["output"]["path"]).resolve()
        != lr_reverse_complement
    ):
        raise ValueError(
            "LR reverse-complement provenance does not match the frozen source"
        )
    lr_reverse_sha256 = str(lr_reverse_provenance["output"]["sha256"])
    if sha256_file(lr_reverse_complement) != lr_reverse_sha256:
        raise ValueError(
            "LR reverse-complement FASTQ digest does not match its provenance"
        )
    binaries = {
        "seqproc": args.seqproc_binary.resolve(),
        "matchbox": args.matchbox_binary.resolve(),
        "splitcode": args.splitcode_binary.resolve(),
    }
    fastq_auditor = args.fastq_audit_binary.resolve()
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
    add_artifact(fastq_auditor)
    add_artifact(Path(sys.executable).resolve())
    add_artifact(provenance_path)
    add_artifact(lr_reverse_complement, lr_reverse_sha256)
    add_artifact(lr_reverse_provenance_path)
    for binary in binaries.values():
        add_artifact(binary)
    for path in (
        Path(__file__),
        ROOT / "scripts" / "benchmark_harness.py",
        ROOT / "scripts" / "run_frozen_schedule.py",
        ROOT / "scripts" / "summarize_frozen_schedule.py",
        ROOT / "scripts" / "build_fastq_numeric_audit.py",
        ROOT / "scripts" / "run_splitcode_dual_pass.py",
        ROOT / "scripts" / "stage_lr_reverse_complement.py",
        ROOT / "tools" / "fastq_numeric_audit.rs",
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
        r1_reverse = (
            lr_reverse_complement if dataset.get("r1_reverse_complement") else None
        )
        resolved_datasets[str(dataset["name"])] = {
            "r1": r1,
            "r2": r2,
            "r1_reverse": r1_reverse,
        }
        dataset_entry = {
            "name": dataset["name"],
            "analysis_role": dataset["analysis_role"],
            "records": int(r1_record["records"]),
            "r1": frozen(r1, str(r1_record["sha256"])),
            "r2": frozen(r2, str(r2_record["sha256"])) if r2 else None,
            "equivalence_scope": dataset["equivalence"],
        }
        if r1_reverse is not None:
            dataset_entry["derived_inputs"] = {
                "r1_reverse_complement": frozen(r1_reverse, lr_reverse_sha256),
                "provenance": frozen(lr_reverse_provenance_path),
                "staged_outside_measurement": True,
            }
        dataset_entries.append(dataset_entry)
        config_paths = (
            ROOT / "configs" / "seqproc" / str(dataset["seqproc"]),
            ROOT / "configs" / "matchbox" / str(dataset["matchbox"]),
            ROOT / "configs" / "splitcode" / str(dataset["splitcode"]),
            *(
                ROOT / "configs" / "seqproc" / name
                for name in dataset["seqproc_support"]
            ),
            *(
                ROOT / "configs" / "matchbox" / name
                for name in dataset["matchbox_support"]
            ),
        )
        for path in config_paths:
            add_artifact(path)

    runs: list[dict[str, Any]] = []
    for replicate in range(1, args.replicates + 1):
        for dataset in DATASETS:
            dataset_name = str(dataset["name"])
            r1 = resolved_datasets[dataset_name]["r1"]
            r2 = resolved_datasets[dataset_name]["r2"]
            r1_reverse = resolved_datasets[dataset_name]["r1_reverse"]
            base_inputs = [frozen(r1, by_name[r1.name]["sha256"])]
            if r2 is not None:
                base_inputs.append(frozen(r2, by_name[r2.name]["sha256"]))
            for threads in args.threads:
                block = f"{dataset_name}-t{threads}-r{replicate}"
                last_cpu = args.first_cpu + threads - 1
                affinity = (
                    str(args.first_cpu)
                    if threads == 1
                    else f"{args.first_cpu}-{last_cpu}"
                )
                for tool in ("seqproc", "matchbox", "splitcode"):
                    tool_binary = binaries[tool]
                    command = [str(TASKSET), "--cpu-list", affinity, str(tool_binary)]
                    inputs = list(base_inputs)
                    configs: list[Path]
                    outputs: list[dict[str, Any]] = []
                    if tool == "seqproc":
                        geometry = (
                            ROOT / "configs" / "seqproc" / str(dataset["seqproc"])
                        )
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
                                "normalize": "fastq_numeric_accession_set",
                                "numeric_id_max": dataset_records[dataset_name],
                                "numeric_audit_executable": str(fastq_auditor),
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
                                    "normalize": "fastq_numeric_accession_set",
                                    "numeric_id_max": dataset_records[dataset_name],
                                    "numeric_audit_executable": str(fastq_auditor),
                                    "mate": 2,
                                    "min_bytes": 1,
                                }
                            )
                        for support_path in support:
                            if support_path.name.endswith("_seq2seq.tsv"):
                                command.extend(["--additional", str(support_path)])
                        repos = [
                            repository_records[name]
                            for name in ("analysis", "seqproc", "antisequence")
                        ]
                    elif tool == "matchbox":
                        script = (
                            ROOT / "configs" / "matchbox" / str(dataset["matchbox"])
                        )
                        support = [
                            ROOT / "configs" / "matchbox" / name
                            for name in dataset["matchbox_support"]
                        ]
                        configs = [script, *support]
                        command.extend(
                            [
                                "--error",
                                "0.2",
                                "--match-mode",
                                "one-best",
                                "--threads",
                                str(threads),
                                "--script-file",
                                str(script),
                                str(r1),
                            ]
                        )
                        if r2 is not None:
                            command.extend(["--paired-with", str(r2)])
                        outputs.append(
                            {
                                "path": str(ROOT / "mb_r1.fq"),
                                "format": "fastq",
                                "normalize": "fastq_numeric_accession_set",
                                "numeric_id_max": dataset_records[dataset_name],
                                "numeric_audit_executable": str(fastq_auditor),
                                "mate": 1,
                                "min_bytes": 1,
                            }
                        )
                        if r2 is not None:
                            outputs.append(
                                {
                                    "path": str(ROOT / "mb_r2.fq"),
                                    "format": "fastq",
                                    "normalize": "fastq_numeric_accession_set",
                                    "numeric_id_max": dataset_records[dataset_name],
                                    "numeric_audit_executable": str(fastq_auditor),
                                    "mate": 2,
                                    "min_bytes": 1,
                                }
                            )
                        repos = [
                            repository_records[name]
                            for name in ("analysis", "matchbox")
                        ]
                    else:
                        config = (
                            ROOT / "configs" / "splitcode" / str(dataset["splitcode"])
                        )
                        configs = [config]
                        if dataset["splitcode_dual_pass"]:
                            if r1_reverse is None:
                                raise ValueError(
                                    f"{dataset_name} requires a reverse-complement input"
                                )
                            inputs.append(frozen(r1_reverse, lr_reverse_sha256))
                            command = [
                                str(TASKSET),
                                "--cpu-list",
                                affinity,
                                str(Path(sys.executable).resolve()),
                                str(SPLITCODE_DUAL_WRAPPER),
                                "--binary",
                                str(tool_binary),
                                "--config",
                                str(config),
                                "--threads",
                                str(threads),
                                "--forward-input",
                                str(r1),
                                "--reverse-input",
                                str(r1_reverse),
                                "--forward-output",
                                "{run_dir}/R1.forward.fastq",
                                "--reverse-output",
                                "{run_dir}/R1.reverse.fastq",
                                "--report",
                                "{run_dir}/dual-pass.json",
                            ]
                            for path in (
                                "{run_dir}/R1.forward.fastq",
                                "{run_dir}/R1.reverse.fastq",
                            ):
                                outputs.append(
                                    {
                                        "path": path,
                                        "format": "fastq",
                                        "normalize": "fastq_numeric_accession_set",
                                        "numeric_id_max": dataset_records[dataset_name],
                                        "numeric_audit_executable": str(fastq_auditor),
                                        "mate": 1,
                                        "min_bytes": 1,
                                    }
                                )
                            outputs.append(
                                {"path": "{run_dir}/dual-pass.json", "min_bytes": 1}
                            )
                        else:
                            command.extend(["--config", str(config)])
                            if dataset["splitcode_assign"]:
                                command.extend(
                                    ["--assign", "--mapping", "{run_dir}/mapping.txt"]
                                )
                            command.extend(
                                [
                                    "--nFastqs",
                                    "2" if r2 is not None else "1",
                                    "--threads",
                                    str(threads),
                                ]
                            )
                            output_paths = ["{run_dir}/R1.fastq"]
                            if r2 is not None:
                                output_paths.append("{run_dir}/R2.fastq")
                            command.extend(
                                ["--output", ",".join(output_paths), str(r1)]
                            )
                            if r2 is not None:
                                command.append(str(r2))
                            outputs.extend(
                                {
                                    "path": path,
                                    "format": "fastq",
                                    "normalize": "fastq_numeric_accession_set",
                                    "numeric_id_max": dataset_records[dataset_name],
                                    "numeric_audit_executable": str(fastq_auditor),
                                    "mate": index,
                                    "min_bytes": 1,
                                }
                                for index, path in enumerate(output_paths, start=1)
                            )
                        repos = [
                            repository_records[name]
                            for name in ("analysis", "splitcode")
                        ]

                    condition_id = f"{block}-{tool}"
                    executable_records = [frozen(tool_binary), frozen(fastq_auditor)]
                    if tool == "splitcode" and dataset["splitcode_dual_pass"]:
                        executable_records.append(
                            frozen(Path(sys.executable).resolve())
                        )
                    if dataset_name == "lr_splitseq_dual":
                        execution_mode = (
                            "dual-pass-no-reconciliation"
                            if tool == "splitcode"
                            else "native-dual-orientation"
                        )
                    elif dataset_name == "lr_splitseq_forward":
                        execution_mode = "forward-only-supplementary"
                    else:
                        execution_mode = "best-practical-uncompressed"
                    runs.append(
                        {
                            "id": condition_id,
                            "block": block,
                            "spec": {
                                "name": "seqproc-publication-core",
                                "cwd": str(ROOT),
                                "command": command,
                                "executables": executable_records,
                                "inputs": inputs,
                                "configs": [frozen(path) for path in configs],
                                "outputs": outputs,
                                "retain_outputs": False,
                                "repositories": repos,
                                "metadata": {
                                    "dataset": dataset_name,
                                    "tool": tool,
                                    "execution_mode": execution_mode,
                                    "analysis_role": dataset["analysis_role"],
                                    "orientation": (
                                        "dual"
                                        if dataset_name == "lr_splitseq_dual"
                                        else "forward"
                                        if dataset_name == "lr_splitseq_forward"
                                        else "not-applicable"
                                    ),
                                    **(
                                        {"duplicate_reconciliation_measured": False}
                                        if dataset_name == "lr_splitseq_dual"
                                        else {}
                                    ),
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
            "purpose": (
                "full-data cross-tool scaling and resource benchmark with native dual-"
                "orientation LR-SPLiT-seq primary and forward-only supplementary comparison"
            ),
            "require_cross_tool_identity": False,
        },
        "artifacts": sorted(artifacts.values(), key=lambda item: item["path"]),
        "data_provenance": frozen(provenance_path),
        "derived_data": {
            "lr_reverse_complement": frozen(lr_reverse_complement, lr_reverse_sha256),
            "lr_reverse_complement_provenance": frozen(lr_reverse_provenance_path),
            "staged_outside_measurement": True,
        },
        "datasets": dataset_entries,
        "execution": {
            "timeout_seconds": args.timeout_seconds,
            "sanitized_environment_allowlist": [
                "PATH",
                "LD_LIBRARY_PATH",
                "LANG",
                "TMPDIR",
            ],
            "cpu_policy": (
                f"requested threads pinned to physical CPUs {args.first_cpu} through "
                f"{args.first_cpu} + N - 1; SMT siblings excluded"
            ),
        },
        "runs": runs,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-provenance", type=Path, required=True)
    parser.add_argument(
        "--seqproc-binary",
        type=Path,
        default=ECOSYSTEM / "seqproc" / "target" / "release" / "seqproc",
    )
    parser.add_argument(
        "--matchbox-binary",
        type=Path,
        default=ECOSYSTEM
        / "competitors"
        / "matchbox-v0.3.2"
        / "target"
        / "release"
        / "matchbox",
    )
    parser.add_argument(
        "--splitcode-binary",
        type=Path,
        default=ECOSYSTEM
        / "competitors"
        / "splitcode-v0.31.6"
        / "build"
        / "src"
        / "splitcode",
    )
    parser.add_argument(
        "--fastq-audit-binary", type=Path, default=DEFAULT_FASTQ_AUDITOR
    )
    parser.add_argument("--seqproc-repo", type=Path, default=ECOSYSTEM / "seqproc")
    parser.add_argument(
        "--antisequence-repo", type=Path, default=ECOSYSTEM / "antisequence"
    )
    parser.add_argument(
        "--matchbox-repo",
        type=Path,
        default=ECOSYSTEM / "competitors" / "matchbox-v0.3.2",
    )
    parser.add_argument(
        "--splitcode-repo",
        type=Path,
        default=ECOSYSTEM / "competitors" / "splitcode-v0.31.6",
    )
    parser.add_argument("--threads", type=int, nargs="+", default=list(DEFAULT_THREADS))
    parser.add_argument("--replicates", type=int, default=DEFAULT_REPLICATES)
    parser.add_argument(
        "--lr-reverse-complement",
        type=Path,
        default=DEFAULT_LR_REVERSE_COMPLEMENT,
    )
    parser.add_argument(
        "--lr-reverse-complement-provenance",
        type=Path,
        default=DEFAULT_LR_REVERSE_COMPLEMENT_PROVENANCE,
    )
    parser.add_argument(
        "--first-cpu",
        type=int,
        default=1,
        help="first physical CPU used for affinity (default: 1; avoids the noisy CPU 0 on the benchmark node)",
    )
    parser.add_argument("--seed", type=int, default=741211)
    parser.add_argument("--timeout-seconds", type=int, default=21600)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if (
        args.replicates <= 0
        or args.timeout_seconds <= 0
        or args.first_cpu < 0
        or any(value <= 0 or args.first_cpu + value > 64 for value in args.threads)
    ):
        parser.error(
            "replicates/timeout must be positive and requested physical CPU IDs must be in 0..63"
        )
    if args.output.exists():
        parser.error(f"refusing to overwrite {args.output}")
    manifest = build_manifest(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(yaml.safe_dump(manifest, sort_keys=False))
    print(f"wrote {args.output} with {len(manifest['runs'])} conditions")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
