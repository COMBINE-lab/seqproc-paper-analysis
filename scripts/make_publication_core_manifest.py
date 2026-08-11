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
PRIMARY_OUTPUT_LENGTHS = {
    "splitseq_pe": {
        "seqproc": {1: (66, 66), 2: (30, 30)},
        "matchbox": {1: (66, 66), 2: (30, 30)},
        "splitcode": {1: (66, 66)},
    },
    "lr_splitseq_dual": {
        "seqproc": {1: (32, 32)},
        "matchbox": {1: (32, 32)},
    },
    "lr_splitseq_forward": {
        "seqproc": {1: (32, 32)},
        "matchbox": {1: (32, 32)},
    },
    "tenx_v2": {
        tool: {1: (26, 26), 2: (98, 98)}
        for tool in ("seqproc", "matchbox", "splitcode")
    },
    "scirnaseq3": {
        "seqproc": {1: (30, 30), 2: (56, 56)},
        # Fuzzy anchor matching can absorb up to two bases into matchbox's
        # projected BC1 segment.  This is a best-practical semantic difference,
        # not malformed FASTQ, so audit the full observed range.
        "matchbox": {1: (27, 30), 2: (56, 56)},
    },
}
SPLITCODE_EXTRACTION_LENGTHS = {
    "umi_bc3_bc2_bc1.fastq": (30, 30),
    "prefix.fastq": (18, 18),
    "bc2.fastq": (8, 8),
    "bc1.fastq": (6, 6),
    "splitcode_sci_rna_seq3.fastq": (27, 28),
}
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
        "splitcode_select": (0,),
        "splitcode_no_outb": True,
        "splitcode_extract_outputs": ("umi_bc3_bc2_bc1.fastq",),
        "analysis_role": "primary",
        "semantic_workload": "filter, whitelist-match, and project UMI+BC3+BC2+BC1",
        "equivalence": (
            "best-practical; seqproc and splitcode emit canonical barcode values, "
            "matchbox emits matched segments, and splitcode linker matching is "
            "substitution-only"
        ),
    },
    {
        "name": "lr_splitseq_dual",
        "r1": "SRR13948564_full.fastq",
        "r2": None,
        "r1_reverse_complement": True,
        "seqproc": "splitseq_singleend_edit_ann.geom",
        "matchbox": "publication_lr_splitseq_dual.mb",
        "splitcode": "publication_lr_splitseq.config",
        "seqproc_support": (),
        "matchbox_support": (),
        "splitcode_assign": True,
        "splitcode_dual_pass": True,
        "splitcode_x_only": True,
        "splitcode_extract_outputs": ("prefix.fastq", "bc2.fastq", "bc1.fastq"),
        "analysis_role": "primary",
        "semantic_workload": (
            "filter both orientations and project UMI+BC3+BC2+BC1; splitcode "
            "materializes the projection as three component FASTQs"
        ),
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
        "splitcode_x_only": True,
        "splitcode_extract_outputs": ("prefix.fastq", "bc2.fastq", "bc1.fastq"),
        "analysis_role": "supplementary",
        "semantic_workload": (
            "filter the forward orientation and project UMI+BC3+BC2+BC1; "
            "splitcode materializes the projection as three component FASTQs"
        ),
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
        "splitcode": "publication_10x_v2_filter.config",
        "seqproc_support": (),
        "matchbox_support": (),
        "splitcode_assign": False,
        "splitcode_dual_pass": False,
        "splitcode_trim_only": True,
        "splitcode_extract_outputs": (),
        "analysis_role": "primary",
        "semantic_workload": "reject pairs with R1 shorter than 26 nt; otherwise passthrough",
        "equivalence": "length-filter-only task; no splitcode extraction",
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
        "splitcode_select": (1,),
        "splitcode_no_outb": True,
        "splitcode_extract_outputs": ("splitcode_sci_rna_seq3.fastq",),
        "analysis_role": "primary",
        "semantic_workload": "anchor-filter and project BC1+BC2+UMI",
        "equivalence": "best-practical; indel semantics differ by backend",
    },
)


def sequence_length_contract(
    dataset: str, tool: str, mate: int
) -> dict[str, Any]:
    bounds = PRIMARY_OUTPUT_LENGTHS.get(dataset, {}).get(tool, {}).get(mate)
    if bounds is None:
        return {}
    nominal = list(range(bounds[0], bounds[1] + 1))
    if (
        dataset.startswith("lr_splitseq") and tool in ("seqproc", "matchbox")
    ) or (dataset == "splitseq_pe" and tool == "matchbox" and mate == 2):
        return {
            "nominal_sequence_lengths": nominal,
            "enforce_sequence_lengths": False,
        }
    return {
        "min_sequence_length": bounds[0],
        "max_sequence_length": bounds[1],
        "nominal_sequence_lengths": nominal,
    }


def extraction_length_contract(name: str) -> dict[str, Any]:
    bounds = SPLITCODE_EXTRACTION_LENGTHS[name]
    nominal = list(range(bounds[0], bounds[1] + 1))
    if name in ("bc2.fastq", "bc1.fastq"):
        return {
            "nominal_sequence_lengths": nominal,
            "enforce_sequence_lengths": False,
        }
    return {
        "min_sequence_length": bounds[0],
        "max_sequence_length": bounds[1],
        "nominal_sequence_lengths": nominal,
    }


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
            "semantic_workload": dataset["semantic_workload"],
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

    measurement_track = str(args.measurement_track)
    benchmark_threads = (
        tuple(args.threads)
        if measurement_track == "timing"
        else (int(args.correctness_threads),)
    )
    benchmark_replicates = args.replicates if measurement_track == "timing" else 1
    runs: list[dict[str, Any]] = []
    for replicate in range(1, benchmark_replicates + 1):
        for dataset in DATASETS:
            dataset_name = str(dataset["name"])
            r1 = resolved_datasets[dataset_name]["r1"]
            r2 = resolved_datasets[dataset_name]["r2"]
            r1_reverse = resolved_datasets[dataset_name]["r1_reverse"]
            base_inputs = [frozen(r1, by_name[r1.name]["sha256"])]
            if r2 is not None:
                base_inputs.append(frozen(r2, by_name[r2.name]["sha256"]))
            for threads in benchmark_threads:
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
                    execution_cwd: str | None = None
                    runtime_symlinks: list[dict[str, str]] = []
                    if tool == "seqproc":
                        geometry = (
                            ROOT / "configs" / "seqproc" / str(dataset["seqproc"])
                        )
                        support = [
                            ROOT / "configs" / "seqproc" / name
                            for name in dataset["seqproc_support"]
                        ]
                        configs = [geometry, *support]
                        seqproc_out1 = (
                            "/dev/null"
                            if measurement_track == "timing"
                            else "{run_dir}/R1.fastq"
                        )
                        command.extend(
                            [
                                "run",
                                "--geom",
                                str(geometry),
                                "--file1",
                                str(r1),
                                "--out1",
                                seqproc_out1,
                                "--threads",
                                str(threads),
                                "--staged-pipeline",
                            ]
                        )
                        if measurement_track == "correctness":
                            outputs.append(
                                {
                                    "path": "{run_dir}/R1.fastq",
                                    "format": "fastq",
                                    "normalize": "fastq_numeric_accession_set",
                                    "numeric_id_max": dataset_records[dataset_name],
                                    "numeric_audit_executable": str(fastq_auditor),
                                    "mate": 1,
                                    "min_bytes": 1,
                                    "logical_product": "primary output 1",
                                    **sequence_length_contract(
                                        dataset_name, "seqproc", 1
                                    ),
                                }
                            )
                        if r2 is not None:
                            seqproc_out2 = (
                                "/dev/null"
                                if measurement_track == "timing"
                                else "{run_dir}/R2.fastq"
                            )
                            command.extend(["--file2", str(r2), "--out2", seqproc_out2])
                            if measurement_track == "correctness":
                                outputs.append(
                                    {
                                        "path": "{run_dir}/R2.fastq",
                                        "format": "fastq",
                                        "normalize": "fastq_numeric_accession_set",
                                        "numeric_id_max": dataset_records[dataset_name],
                                        "numeric_audit_executable": str(fastq_auditor),
                                        "mate": 2,
                                        "min_bytes": 1,
                                        "logical_product": "primary output 2",
                                        **sequence_length_contract(
                                            dataset_name, "seqproc", 2
                                        ),
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
                        execution_cwd = "{run_dir}/tool-work"
                        runtime_symlinks.append(
                            {
                                "path": "{run_dir}/tool-work/configs",
                                "target": str(ROOT / "configs"),
                            }
                        )
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
                        if measurement_track == "timing":
                            runtime_symlinks.append(
                                {
                                    "path": "{run_dir}/tool-work/mb_r1.fq",
                                    "target": "/dev/null",
                                }
                            )
                        else:
                            outputs.append(
                                {
                                    "path": "{run_dir}/tool-work/mb_r1.fq",
                                    "format": "fastq",
                                    "normalize": "fastq_numeric_accession_set",
                                    "numeric_id_max": dataset_records[dataset_name],
                                    "numeric_audit_executable": str(fastq_auditor),
                                    "mate": 1,
                                    "min_bytes": 1,
                                    "logical_product": "primary output 1",
                                    **sequence_length_contract(
                                        dataset_name, "matchbox", 1
                                    ),
                                }
                            )
                        if r2 is not None:
                            if measurement_track == "timing":
                                runtime_symlinks.append(
                                    {
                                        "path": "{run_dir}/tool-work/mb_r2.fq",
                                        "target": "/dev/null",
                                    }
                                )
                            else:
                                outputs.append(
                                    {
                                        "path": "{run_dir}/tool-work/mb_r2.fq",
                                        "format": "fastq",
                                        "normalize": "fastq_numeric_accession_set",
                                        "numeric_id_max": dataset_records[dataset_name],
                                        "numeric_audit_executable": str(fastq_auditor),
                                        "mate": 2,
                                        "min_bytes": 1,
                                        "logical_product": "primary output 2",
                                        **sequence_length_contract(
                                            dataset_name, "matchbox", 2
                                        ),
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
                        execution_cwd = "{run_dir}/tool-work"
                        extract_names = tuple(dataset["splitcode_extract_outputs"])
                        if (
                            measurement_track == "timing"
                            and not dataset["splitcode_dual_pass"]
                        ):
                            runtime_symlinks.extend(
                                {
                                    "path": f"{{run_dir}}/tool-work/{name}",
                                    "target": "/dev/null",
                                }
                                for name in extract_names
                            )
                        if dataset["splitcode_dual_pass"]:
                            if r1_reverse is None:
                                raise ValueError(
                                    f"{dataset_name} requires a reverse-complement input"
                                )
                            inputs.append(frozen(r1_reverse, lr_reverse_sha256))
                            forward_output = (
                                "/dev/null"
                                if measurement_track == "timing"
                                else "{run_dir}/R1.forward.fastq"
                            )
                            reverse_output = (
                                "/dev/null"
                                if measurement_track == "timing"
                                else "{run_dir}/R1.reverse.fastq"
                            )
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
                                forward_output,
                                "--reverse-output",
                                reverse_output,
                                "--report",
                                "{run_dir}/dual-pass.json",
                            ]
                            for name in extract_names:
                                command.extend(["--extraction-output", name])
                            if dataset.get("splitcode_x_only"):
                                command.append("--x-only")
                            if measurement_track == "timing":
                                command.extend(
                                    [
                                        "--discard-output",
                                        "--mapping-sink",
                                        "/dev/null",
                                    ]
                                )
                            elif not dataset.get("splitcode_x_only"):
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
                                            "logical_product": "accepted full read",
                                            **sequence_length_contract(
                                                dataset_name, "splitcode", 1
                                            ),
                                        }
                                    )
                            outputs.append(
                                {"path": "{run_dir}/dual-pass.json", "min_bytes": 1}
                            )
                        else:
                            command.extend(["--config", str(config)])
                            if dataset["splitcode_assign"]:
                                command.extend(
                                    [
                                        "--assign",
                                        "--mapping",
                                        "/dev/null"
                                        if measurement_track == "timing"
                                        else "{run_dir}/mapping.txt",
                                    ]
                                )
                            if dataset.get("splitcode_trim_only"):
                                command.append("--trim-only")
                            if dataset.get("splitcode_no_outb"):
                                command.append("--no-outb")
                            selected_outputs = tuple(dataset.get("splitcode_select", ()))
                            if selected_outputs:
                                command.extend(
                                    [
                                        "--select",
                                        ",".join(str(value) for value in selected_outputs),
                                    ]
                                )
                            command.extend(
                                [
                                    "--nFastqs",
                                    "2" if r2 is not None else "1",
                                    "--threads",
                                    str(threads),
                                ]
                            )
                            n_fastqs = 2 if r2 is not None else 1
                            output_paths = [
                                (
                                    f"{{run_dir}}/R{index + 1}.fastq"
                                    if measurement_track == "correctness"
                                    and not dataset.get("splitcode_no_outb")
                                    and (
                                        not selected_outputs
                                        or index in selected_outputs
                                    )
                                    else "/dev/null"
                                )
                                for index in range(n_fastqs)
                            ]
                            if dataset.get("splitcode_x_only"):
                                command.append("--x-only")
                            else:
                                # splitcode 0.31.6 validates --output against nFastqs
                                # even when --select narrows the emitted files.
                                command.extend(["--output", ",".join(output_paths)])
                            command.append(str(r1))
                            if r2 is not None:
                                command.append(str(r2))
                            if (
                                measurement_track == "correctness"
                                and not dataset.get("splitcode_x_only")
                                and not dataset.get("splitcode_no_outb")
                            ):
                                outputs.extend(
                                    {
                                        "path": path,
                                        "format": "fastq",
                                        "normalize": "fastq_numeric_accession_set",
                                        "numeric_id_max": dataset_records[dataset_name],
                                        "numeric_audit_executable": str(fastq_auditor),
                                        "mate": index,
                                        "min_bytes": 1,
                                        "logical_product": "accepted full read",
                                        **sequence_length_contract(
                                            dataset_name, "splitcode", index
                                        ),
                                    }
                                    for index, path in enumerate(output_paths, start=1)
                                    if path != "/dev/null"
                                )
                        if measurement_track == "correctness":
                            extraction_products = (
                                (
                                    (f"splitcode-{orientation}-work/{name}", orientation)
                                    for orientation in ("forward", "reverse")
                                    for name in extract_names
                                )
                                if dataset["splitcode_dual_pass"]
                                else (
                                    (f"tool-work/{name}", "forward")
                                    for name in extract_names
                                )
                            )
                            outputs.extend(
                                {
                                    "path": f"{{run_dir}}/{relative}",
                                    "format": "fastq",
                                    "normalize": "fastq_numeric_accession_set",
                                    "numeric_id_max": dataset_records[dataset_name],
                                    "numeric_audit_executable": str(fastq_auditor),
                                    "mate": 10 + index,
                                    "min_bytes": 1,
                                    "logical_product": f"{orientation} extracted component",
                                    **extraction_length_contract(
                                        Path(relative).name
                                    ),
                                }
                                for index, (relative, orientation) in enumerate(
                                    extraction_products, start=1
                                )
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
                                **(
                                    {
                                        "execution_cwd": execution_cwd,
                                        "create_execution_cwd": True,
                                    }
                                    if execution_cwd is not None
                                    else {}
                                ),
                                **(
                                    {"runtime_symlinks": runtime_symlinks}
                                    if runtime_symlinks
                                    else {}
                                ),
                                "command": command,
                                "executables": executable_records,
                                "inputs": inputs,
                                "configs": [frozen(path) for path in configs],
                                "outputs": outputs,
                                "retain_outputs": measurement_track == "correctness",
                                "repositories": repos,
                                "metadata": {
                                    "dataset": dataset_name,
                                    "tool": tool,
                                    "measurement_track": measurement_track,
                                    "sequence_output_policy": (
                                        "dev-null"
                                        if measurement_track == "timing"
                                        else "materialized-and-validated"
                                    ),
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
                                    **(
                                        {"allow_output_count_mismatch": True}
                                        if dataset_name == "lr_splitseq_dual"
                                        and tool == "splitcode"
                                        and measurement_track == "correctness"
                                        else {}
                                    ),
                                    "threads": threads,
                                    "replicate": replicate,
                                    "input_records": dataset_records[dataset_name],
                                    "cpu_affinity": affinity,
                                    "equivalence_scope": dataset["equivalence"],
                                    "semantic_workload": dataset["semantic_workload"],
                                    "tier": f"publication-core-full-data-{measurement_track}",
                                },
                            },
                        }
                    )

    return {
        "schema_version": "1.0.0",
        "mode": "publication",
        "study": {
            "name": f"seqproc-publication-core-full-data-{measurement_track}",
            "random_seed": args.seed,
            "purpose": (
                "full-data cross-tool benchmark with qualitatively aligned workloads; "
                + (
                    "primary runtime/RSS measurement with biological sequence outputs "
                    "directed to /dev/null"
                    if measurement_track == "timing"
                    else "materialized FASTQ correctness validation outside timing claims"
                )
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
        "--measurement-track",
        choices=("timing", "correctness"),
        default="timing",
        help=(
            "timing directs biological outputs to /dev/null; correctness runs one "
            "materialized full-data replicate"
        ),
    )
    parser.add_argument(
        "--correctness-threads",
        type=int,
        default=16,
        help="thread count for the single materialized correctness replicate",
    )
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
        or args.correctness_threads <= 0
        or args.first_cpu + args.correctness_threads > 64
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
