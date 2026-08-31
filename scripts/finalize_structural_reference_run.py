#!/usr/bin/env python3
"""Finalize a structural-reference run into one machine-verifiable artifact.

The validator streams very large ID lists, so it deliberately does not spend a
second pass hashing them.  This small finalization step hashes and counts every
ID list, checks the count against its summary, attaches the already-verified
input FASTQ checksums, and records the exact reconstruction command.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_PROVENANCE = Path(
    "/scratch1/seqproc-benchmark-data/full/full_data_provenance.json"
)


def digest_and_lines(path: Path) -> tuple[str, int, int]:
    digest = hashlib.sha256()
    lines = 0
    size = 0
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
            lines += chunk.count(b"\n")
            size += len(chunk)
    return digest.hexdigest(), lines, size


def sha256(path: Path) -> str:
    return digest_and_lines(path)[0]


def git_head() -> str:
    return subprocess.run(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def command_for(summary: dict, output: Path, summary_path: Path) -> list[str]:
    command = [
        ".venv/bin/python",
        "scripts/edit_tolerant_validity.py",
        summary["fastq"],
    ]
    if summary.get("mate_fastq"):
        command.extend(["--r2", summary["mate_fastq"]])
    command.extend(["--chem", summary["chem"], "--threads", str(summary["threads"])])
    criteria = summary["criteria"]
    if summary["chem"] == "lr":
        command.extend(
            [
                "--max-linker1-edit",
                str(criteria["max_linker1_edit"]),
                "--max-linker2-edit",
                str(criteria["max_linker2_edit"]),
            ]
        )
    command.extend(
        [
            "--out",
            str(output),
            "--summary-json",
            str(summary_path),
        ]
    )
    return command


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_directory", type=Path)
    parser.add_argument(
        "--data-provenance", type=Path, default=DEFAULT_DATA_PROVENANCE
    )
    parser.add_argument("--output", default="provenance.json")
    args = parser.parse_args(argv)

    run_directory = args.run_directory.resolve()
    source = json.loads(args.data_provenance.read_text())
    inputs = {
        str(Path(item["fastq"]["path"]).resolve()): item["fastq"]
        for item in source["files"]
    }
    validator = ROOT / "scripts" / "edit_tolerant_validity.py"
    artifacts = []
    for summary_path in sorted(run_directory.glob("*.summary.json")):
        summary = json.loads(summary_path.read_text())
        id_path = summary_path.with_name(
            summary_path.name.replace(".summary.json", ".valid_ids.txt")
        )
        digest, lines, size = digest_and_lines(id_path)
        if lines != summary["valid"]:
            raise ValueError(
                f"{id_path}: {lines:,} ID lines != summary valid={summary['valid']:,}"
            )
        input_records = []
        for key in ("fastq", "mate_fastq"):
            path = summary.get(key)
            if not path:
                continue
            resolved = str(Path(path).resolve())
            if resolved not in inputs:
                raise KeyError(f"input absent from data provenance: {resolved}")
            item = inputs[resolved]
            input_records.append(
                {
                    "path": resolved,
                    "sha256": item["sha256"],
                    "records": item.get("records"),
                }
            )
        artifacts.append(
            {
                "name": summary_path.name.removesuffix(".summary.json"),
                "summary": summary,
                "summary_path": str(summary_path.relative_to(ROOT)),
                "valid_ids": {
                    "path": str(id_path.relative_to(ROOT)),
                    "sha256": digest,
                    "bytes": size,
                    "records": lines,
                    "ordering": "input_order",
                },
                "inputs": input_records,
                "command": command_for(summary, id_path, summary_path),
            }
        )

    payload = {
        "schema_version": "1.0.0",
        "repository_base_commit_at_run": git_head(),
        "provenance_note": (
            "Outputs were generated from a working tree based on this commit. The "
            "recorded validator, whitelist, input, summary, and output content hashes "
            "are the authoritative identities; the changes are committed together "
            "with this provenance artifact."
        ),
        "validator": {
            "path": str(validator.relative_to(ROOT)),
            "sha256": sha256(validator),
        },
        "whitelists": [
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": sha256(path),
            }
            for path in (
                ROOT / "configs/seqproc/splitseq_bc8_whitelist.txt",
                ROOT / "configs/seqproc/splitseq_bc1_whitelist_6bp.txt",
            )
        ],
        "data_provenance": {
            "path": str(args.data_provenance),
            "sha256": sha256(args.data_provenance),
        },
        "artifacts": artifacts,
    }
    output = run_directory / args.output
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
