#!/usr/bin/env python3
"""Summarize completed publication correctness runs against structural references.

The frozen harness deliberately removes its temporary accession bitmaps after
validating each FASTQ.  This script recreates compact persistent bitmaps using
the pinned Rust auditor, checks that multiple products from one tool carry the
same read IDs, unions splitcode's two LR orientations, and reports agreement
with the conservative structural reference plus pairwise tool concordance.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_AUDITOR = ROOT / "tools" / "bin" / "fastq-numeric-audit"
DEFAULT_REFERENCES = {
    "splitseq_pe": ROOT
    / "publication_results/structural_reference_2026-08-17/splitseq_pe_core.valid_ids.txt",
    "lr_splitseq_dual": ROOT
    / "publication_results/structural_reference_2026-08-17/lr_splitseq_core_e3_e3.valid_ids.txt",
    "scirnaseq3": ROOT
    / "publication_results/structural_reference_2026-08-17/scirnaseq3_core.valid_ids.txt",
    "tenx_v2": ROOT
    / "publication_results/structural_reference_2026-08-17/tenx_v2_core.valid_ids.txt",
}
TOOLS = ("seqproc", "matchbox", "splitcode")
MAGIC = b"fastq_numeric_accession_set_v1"
POPCOUNT = tuple(bin(value).count("1") for value in range(256))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def successful_attempts(run_root: Path) -> dict[tuple[str, str], tuple[Path, dict, dict]]:
    found = {}
    for identity_path in run_root.glob("*/attempt-*/identity.json"):
        attempt = identity_path.parent
        run_path = attempt / "run.json"
        if not run_path.is_file():
            continue
        run = json.loads(run_path.read_text())
        if run.get("success") is not True:
            continue
        identity = json.loads(identity_path.read_text())
        metadata = identity["metadata"]
        key = (str(metadata["dataset"]), str(metadata["tool"]))
        if key in found:
            raise RuntimeError(f"multiple successful attempts found for {key}")
        found[key] = (attempt, identity, run)
    return found


def canonical_bitmap(path: Path) -> tuple[dict[str, object], bytes]:
    fields = path.read_bytes().split(b"\0", 4)
    if len(fields) != 5 or fields[0] != MAGIC:
        raise RuntimeError(f"unexpected accession bitmap format: {path}")
    return (
        {
            "mate": int(fields[1]),
            "numeric_id_max": int(fields[2]),
            "accession_prefix": fields[3].decode(),
        },
        fields[4],
    )


def audit_fastq(
    fastq: Path,
    bitmap_path: Path,
    numeric_id_max: int,
    auditor: Path,
) -> tuple[bytes, dict]:
    bitmap_path.parent.mkdir(parents=True, exist_ok=True)
    bitmap_path.unlink(missing_ok=True)
    completed = subprocess.run(
        [str(auditor), str(fastq), "1", str(numeric_id_max), str(bitmap_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    audit = json.loads(completed.stdout)
    metadata, bitmap = canonical_bitmap(bitmap_path)
    if metadata["numeric_id_max"] != numeric_id_max:
        raise RuntimeError(f"numeric ID limit differs for {fastq}")
    return bitmap, {
        "fastq": str(fastq),
        "fastq_sha256": sha256(fastq),
        "bitmap": str(bitmap_path),
        "bitmap_sha256": sha256(bitmap_path),
        "audit": audit,
    }


def emitted_bitmap(
    key: tuple[str, str],
    attempt: Path,
    identity: dict,
    run: dict,
    bitmap_root: Path,
    auditor: Path,
) -> tuple[bytes, dict]:
    dataset, tool = key
    numeric_id_max = int(identity["metadata"]["input_records"])
    fastqs = [
        Path(item["path"])
        for item in run.get("outputs", [])
        if item.get("normalization") == "fastq_numeric_accession_set_v1"
    ]
    if not fastqs:
        raise RuntimeError(f"no FASTQ output recorded for {dataset}/{tool}")

    groups: dict[str, list[Path]] = {"accepted": fastqs}
    if dataset == "lr_splitseq_dual" and tool == "splitcode":
        groups = {
            "forward": [path for path in fastqs if "splitcode-forward-work" in str(path)],
            "reverse": [path for path in fastqs if "splitcode-reverse-work" in str(path)],
        }
    union = bytearray((numeric_id_max + 7) // 8)
    sources = []
    group_counts = {}
    for group_name, paths in groups.items():
        if not paths:
            raise RuntimeError(f"empty {group_name} output group for {dataset}/{tool}")
        group_bitmap = None
        for index, fastq in enumerate(paths, 1):
            bitmap_path = bitmap_root / f"{dataset}.{tool}.{group_name}.{index}.bitmap"
            bitmap, source = audit_fastq(fastq, bitmap_path, numeric_id_max, auditor)
            sources.append(source)
            if group_bitmap is None:
                group_bitmap = bitmap
            elif bitmap != group_bitmap:
                raise RuntimeError(
                    f"extracted component ID sets differ for {dataset}/{tool}/{group_name}"
                )
        assert group_bitmap is not None
        group_counts[group_name] = popcount(group_bitmap)
        for offset, value in enumerate(group_bitmap):
            union[offset] |= value
    return bytes(union), {
        "attempt": str(attempt),
        "condition_id": identity["metadata"]["condition_id"],
        "group_records": group_counts,
        "sources": sources,
        "run": {
            "wall_seconds": run.get("wall_seconds"),
            "user_cpu_seconds": run.get("user_cpu_seconds"),
            "system_cpu_seconds": run.get("system_cpu_seconds"),
            "peak_rss_kib": run.get("peak_rss_kib"),
        },
        "repositories": identity.get("repositories", []),
        "executables": identity.get("executables", []),
        "configs": identity.get("configs", []),
    }


def popcount(bitmap: bytes) -> int:
    return sum(POPCOUNT[value] for value in bitmap)


def intersection_count(left: bytes, right: bytes) -> int:
    if len(left) != len(right):
        raise RuntimeError("bitmap lengths differ")
    return sum(POPCOUNT[a & b] for a, b in zip(left, right))


def union_count(left: bytes, right: bytes) -> int:
    if len(left) != len(right):
        raise RuntimeError("bitmap lengths differ")
    return sum(POPCOUNT[a | b] for a, b in zip(left, right))


def reference_bitmap(path: Path, numeric_id_max: int, output: Path) -> bytes:
    raw = bytearray((numeric_id_max + 7) // 8)
    records = 0
    with path.open("rb") as handle:
        for line_number, line in enumerate(handle, 1):
            token = line.strip().split(None, 1)[0]
            if not token:
                continue
            try:
                numeric_id = int(token.rsplit(b".", 1)[1])
            except (IndexError, ValueError) as error:
                raise RuntimeError(f"bad reference ID at {path}:{line_number}") from error
            if not 1 <= numeric_id <= numeric_id_max:
                raise RuntimeError(f"reference ID {numeric_id} outside 1..{numeric_id_max}")
            byte, bit = divmod(numeric_id - 1, 8)
            mask = 1 << bit
            if raw[byte] & mask:
                raise RuntimeError(f"duplicate reference ID {numeric_id} in {path}")
            raw[byte] |= mask
            records += 1
    if records == 0:
        raise RuntimeError(f"empty structural reference: {path}")
    output.write_bytes(raw)
    return bytes(raw)


def reference_argument(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("expected DATASET=PATH")
    dataset, path = value.split("=", 1)
    if not dataset or not path:
        raise argparse.ArgumentTypeError("dataset and path are required")
    return dataset, Path(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--dataset", action="append")
    parser.add_argument("--reference", action="append", type=reference_argument)
    parser.add_argument("--auditor", type=Path, default=DEFAULT_AUDITOR)
    args = parser.parse_args()

    campaign = args.campaign.resolve()
    output = campaign / "aggregates" / "accuracy"
    bitmap_root = output / "bitmaps"
    output.mkdir(parents=True, exist_ok=True)
    bitmap_root.mkdir(parents=True, exist_ok=True)
    attempts = successful_attempts(campaign / "runs")
    available = sorted({dataset for dataset, _ in attempts})
    datasets = sorted(set(args.dataset or available))
    references = dict(DEFAULT_REFERENCES)
    references.update(dict(args.reference or []))

    missing = [
        (dataset, tool)
        for dataset in datasets
        for tool in TOOLS
        if (dataset, tool) not in attempts
    ]
    if missing:
        raise RuntimeError(f"missing successful correctness runs: {missing}")

    emitted = {}
    provenance = {}
    for dataset in datasets:
        for tool in TOOLS:
            key = (dataset, tool)
            emitted[key], provenance[f"{dataset}/{tool}"] = emitted_bitmap(
                key, *attempts[key], bitmap_root, args.auditor.resolve()
            )

    metrics = []
    pairwise = []
    reference_records = {}
    for dataset in datasets:
        numeric_id_max = int(attempts[(dataset, TOOLS[0])][1]["metadata"]["input_records"])
        reference_path = references.get(dataset)
        if reference_path is None or not reference_path.is_file():
            raise RuntimeError(f"missing structural reference for {dataset}: {reference_path}")
        reference = reference_bitmap(
            reference_path.resolve(), numeric_id_max, bitmap_root / f"{dataset}.reference.raw"
        )
        reference_count = popcount(reference)
        reference_records[dataset] = {
            "path": str(reference_path.resolve()),
            "sha256": sha256(reference_path.resolve()),
            "records": reference_count,
            "interpretation": "conservative structural reference, not biological ground truth",
        }
        for tool in TOOLS:
            tool_set = emitted[(dataset, tool)]
            emitted_count = popcount(tool_set)
            intersection = intersection_count(tool_set, reference)
            precision = intersection / emitted_count if emitted_count else 0.0
            recall = intersection / reference_count if reference_count else 0.0
            f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
            metrics.append(
                {
                    "dataset": dataset,
                    "tool": tool,
                    "input_records": numeric_id_max,
                    "reference_records": reference_count,
                    "emitted_records": emitted_count,
                    "intersection_records": intersection,
                    "emitted_fraction": emitted_count / numeric_id_max,
                    "precision": precision,
                    "recall": recall,
                    "f1": f1,
                }
            )
        for left, right in combinations(TOOLS, 2):
            left_set, right_set = emitted[(dataset, left)], emitted[(dataset, right)]
            intersection = intersection_count(left_set, right_set)
            union = union_count(left_set, right_set)
            pairwise.append(
                {
                    "dataset": dataset,
                    "left": left,
                    "right": right,
                    "intersection_records": intersection,
                    "union_records": union,
                    "left_only_records": popcount(left_set) - intersection,
                    "right_only_records": popcount(right_set) - intersection,
                    "jaccard": intersection / union if union else 1.0,
                }
            )

    payload = {
        "schema_version": "1.0.0",
        "campaign": str(campaign),
        "datasets": datasets,
        "summarizer": {
            "path": str(Path(__file__).resolve()),
            "sha256": sha256(Path(__file__).resolve()),
        },
        "auditor": {"path": str(args.auditor.resolve()), "sha256": sha256(args.auditor.resolve())},
        "manifest": {
            "path": str(campaign / "specs/publication-core-correctness-t32.yaml"),
            "sha256": sha256(campaign / "specs/publication-core-correctness-t32.yaml"),
        },
        "schedule": {
            "path": str(campaign / "specs/publication-core-correctness-t32.schedule.json"),
            "sha256": sha256(
                campaign / "specs/publication-core-correctness-t32.schedule.json"
            ),
        },
        "references": reference_records,
        "metrics": metrics,
        "pairwise": pairwise,
        "provenance": provenance,
    }
    result_path = output / "accuracy_metrics.json"
    result_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    with (output / "accuracy_metrics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metrics[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(metrics)
    with (output / "pairwise_metrics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(pairwise[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(pairwise)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
