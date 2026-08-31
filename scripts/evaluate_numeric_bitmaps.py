#!/usr/bin/env python3
"""Evaluate one or more fastq-numeric-audit bitmap unions against a reference."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


MAGIC = b"fastq_numeric_accession_set_v1"
POPCOUNT = tuple(bin(value).count("1") for value in range(256))


def canonical_bitmap(path: Path) -> bytes:
    fields = path.read_bytes().split(b"\0", 4)
    if len(fields) != 5 or fields[0] != MAGIC:
        raise ValueError(f"not a fastq-numeric-audit bitmap: {path}")
    return fields[4]


def bitmap_union(paths: list[Path]) -> bytes:
    bitmaps = [canonical_bitmap(path) for path in paths]
    if not bitmaps:
        raise ValueError("a tool must provide at least one bitmap")
    if any(len(bitmap) != len(bitmaps[0]) for bitmap in bitmaps):
        raise ValueError("bitmap lengths differ")
    result = bytearray(bitmaps[0])
    for bitmap in bitmaps[1:]:
        for index, value in enumerate(bitmap):
            result[index] |= value
    return bytes(result)


def popcount(bitmap: bytes) -> int:
    return sum(POPCOUNT[value] for value in bitmap)


def intersection_count(left: bytes, right: bytes) -> int:
    if len(left) != len(right):
        raise ValueError("bitmap lengths differ")
    return sum(POPCOUNT[a & b] for a, b in zip(left, right))


def union_count(left: bytes, right: bytes) -> int:
    if len(left) != len(right):
        raise ValueError("bitmap lengths differ")
    return sum(POPCOUNT[a | b] for a, b in zip(left, right))


def tool_argument(value: str) -> tuple[str, list[Path]]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("expected LABEL=BITMAP[,BITMAP...]")
    label, path_list = value.split("=", 1)
    paths = [Path(path) for path in path_list.split(",") if path]
    if not label or not paths:
        raise argparse.ArgumentTypeError("label and at least one bitmap are required")
    return label, paths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--input-records", type=int, required=True)
    parser.add_argument("--tool", action="append", type=tool_argument, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    reference = args.reference.read_bytes()
    reference_count = popcount(reference)
    tools: dict[str, bytes] = {}
    sources: dict[str, list[str]] = {}
    for label, paths in args.tool:
        if label in tools:
            parser.error(f"duplicate tool label: {label}")
        tools[label] = bitmap_union(paths)
        sources[label] = [str(path.resolve()) for path in paths]
    if any(len(bitmap) != len(reference) for bitmap in tools.values()):
        raise ValueError("tool and reference bitmap lengths differ")

    metrics = {}
    for label, bitmap in tools.items():
        emitted = popcount(bitmap)
        intersection = intersection_count(bitmap, reference)
        precision = intersection / emitted if emitted else 0.0
        recall = intersection / reference_count if reference_count else 0.0
        f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        metrics[label] = {
            "sources": sources[label],
            "emitted_records": emitted,
            "intersection_records": intersection,
            "emitted_fraction": emitted / args.input_records,
            "precision": precision,
            "recall": recall,
            "f1": f1,
        }

    pairwise = {}
    labels = list(tools)
    for left_index, left in enumerate(labels):
        for right in labels[left_index + 1 :]:
            intersection = intersection_count(tools[left], tools[right])
            union = union_count(tools[left], tools[right])
            pairwise[f"{left}__{right}"] = {
                "intersection_records": intersection,
                "union_records": union,
                "jaccard": intersection / union if union else 0.0,
                f"{left}_only": popcount(tools[left]) - intersection,
                f"{right}_only": popcount(tools[right]) - intersection,
            }

    payload = {
        "schema_version": "1.0.0",
        "input_records": args.input_records,
        "reference": {
            "path": str(args.reference.resolve()),
            "records": reference_count,
            "interpretation": "conservative structural reference, not biological ground truth",
        },
        "tools": metrics,
        "pairwise": pairwise,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
