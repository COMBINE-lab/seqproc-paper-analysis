#!/usr/bin/env python3
"""Reason-code LR-SPLiT-seq calls unique to one or more seqproc modes.

The input bitmaps are produced by ``fastq-numeric-audit``. Competitor bitmaps
are unioned, and each seqproc cohort is defined as seqproc minus that union.
Target reads are then rescanned once and evaluated with the independent
structural validator used to construct V_total.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from edit_tolerant_validity import (
    LINKERS,
    ham1_set,
    revcomp,
    validate_orientation,
)


POPCOUNT = tuple(bin(value).count("1") for value in range(256))


def popcount(value: bytes | bytearray) -> int:
    return sum(POPCOUNT[byte] for byte in value)


def canonical_bitmap(path: Path) -> bytes:
    fields = path.read_bytes().split(b"\0", 4)
    if len(fields) != 5 or fields[0] != b"fastq_numeric_accession_set_v1":
        raise ValueError(f"not a fastq-numeric-audit bitmap: {path}")
    return fields[4]


def raw_bitmap(path: Path) -> bytes:
    return path.read_bytes()


def union_bitmaps(bitmaps: list[bytes]) -> bytearray:
    if not bitmaps:
        return bytearray()
    size = len(bitmaps[0])
    if any(len(bitmap) != size for bitmap in bitmaps):
        raise ValueError("bitmap lengths differ")
    result = bytearray(size)
    for bitmap in bitmaps:
        for index, byte in enumerate(bitmap):
            result[index] |= byte
    return result


def difference(left: bytes, right: bytes | bytearray) -> bytearray:
    if len(left) != len(right):
        raise ValueError("bitmap lengths differ")
    return bytearray(a & ~b for a, b in zip(left, right))


def intersection_count(left: bytes | bytearray, right: bytes | bytearray) -> int:
    if len(left) != len(right):
        raise ValueError("bitmap lengths differ")
    return sum(POPCOUNT[a & b] for a, b in zip(left, right))


def contains(bitmap: bytes | bytearray, numeric_id: int) -> bool:
    bit = numeric_id - 1
    return bool(bitmap[bit // 8] & (1 << (bit % 8)))


def named_bitmap(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("expected LABEL=BITMAP")
    label, path = value.split("=", 1)
    if not label:
        raise argparse.ArgumentTypeError("bitmap label cannot be empty")
    return label, Path(path)


def counter_dict(counter: Counter) -> dict:
    return {str(key): value for key, value in sorted(counter.items())}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fastq", type=Path)
    parser.add_argument("--seqproc", action="append", type=named_bitmap, required=True)
    parser.add_argument("--competitor", action="append", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--bc23-whitelist", type=Path, required=True)
    parser.add_argument("--bc1-whitelist", type=Path, required=True)
    parser.add_argument("--max-linker1-edit", type=int, default=6)
    parser.add_argument("--max-linker2-edit", type=int, default=6)
    parser.add_argument("--max-linker1-candidates", type=int, default=8)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    competitors = union_bitmaps([canonical_bitmap(path) for path in args.competitor])
    reference = raw_bitmap(args.reference)
    cohorts = {}
    for label, path in args.seqproc:
        if label in cohorts:
            parser.error(f"duplicate seqproc label: {label}")
        emitted = canonical_bitmap(path)
        unique = difference(emitted, competitors)
        cohorts[label] = {
            "path": str(path),
            "emitted": popcount(emitted),
            "bits": unique,
            "unique": popcount(unique),
            "unique_reference": intersection_count(unique, reference),
            "unique_nonreference": popcount(unique) - intersection_count(unique, reference),
            "outcomes": Counter(),
            "reference_outcomes": Counter(),
            "nonreference_outcomes": Counter(),
            "orientations": Counter(),
            "linker1_edits": Counter(),
            "linker2_edits": Counter(),
            "read_lengths": Counter(),
        }

    targets = union_bitmaps([cohort["bits"] for cohort in cohorts.values()])
    bc23, bc23_length = ham1_set(args.bc23_whitelist)
    bc1, bc1_length = ham1_set(args.bc1_whitelist)
    linker1, linker2 = LINKERS["lr"]

    total = 0
    inspected = 0
    with args.fastq.open() as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().strip()
            plus = handle.readline()
            quality = handle.readline()
            if not plus or not quality:
                raise ValueError(f"truncated FASTQ after record {total}: {args.fastq}")
            total += 1
            token = header[1:].split(None, 1)[0]
            try:
                numeric_id = int(token.rsplit(".", 1)[1])
            except (IndexError, ValueError) as error:
                raise ValueError(f"non-numeric accession ID: {token}") from error
            if not contains(targets, numeric_id):
                continue
            inspected += 1
            result = validate_orientation(
                sequence,
                linker1,
                linker2,
                bc23,
                bc1,
                bc23_length,
                bc1_length,
                args.max_linker1_edit,
                args.max_linker2_edit,
                max_candidates=args.max_linker1_candidates,
            )
            orientation = "forward"
            if not result["accepted"]:
                result = validate_orientation(
                    revcomp(sequence),
                    linker1,
                    linker2,
                    bc23,
                    bc1,
                    bc23_length,
                    bc1_length,
                    args.max_linker1_edit,
                    args.max_linker2_edit,
                    max_candidates=args.max_linker1_candidates,
                )
                orientation = "reverse"
            is_reference = contains(reference, numeric_id)
            for cohort in cohorts.values():
                if not contains(cohort["bits"], numeric_id):
                    continue
                reason = result["reason"]
                cohort["outcomes"][reason] += 1
                destination = (
                    cohort["reference_outcomes"]
                    if is_reference
                    else cohort["nonreference_outcomes"]
                )
                destination[reason] += 1
                length_floor = (len(sequence) // 100) * 100
                cohort["read_lengths"][f"{length_floor}-{length_floor + 99}"] += 1
                if result["accepted"]:
                    cohort["orientations"][orientation] += 1
                    cohort["linker1_edits"][result["linker1_edit"]] += 1
                    cohort["linker2_edits"][result["linker2_edit"]] += 1

    output_cohorts = {}
    for label, cohort in cohorts.items():
        if sum(cohort["outcomes"].values()) != cohort["unique"]:
            raise RuntimeError(
                f"did not encounter every {label} target: "
                f"{sum(cohort['outcomes'].values())} of {cohort['unique']}"
            )
        output_cohorts[label] = {
            key: value
            for key, value in cohort.items()
            if key not in {"bits", "outcomes", "reference_outcomes",
                           "nonreference_outcomes", "orientations",
                           "linker1_edits", "linker2_edits", "read_lengths"}
        }
        output_cohorts[label].update(
            {
                "outcomes": counter_dict(cohort["outcomes"]),
                "reference_outcomes": counter_dict(cohort["reference_outcomes"]),
                "nonreference_outcomes": counter_dict(cohort["nonreference_outcomes"]),
                "accepted_orientations": counter_dict(cohort["orientations"]),
                "accepted_linker1_edit_histogram": counter_dict(cohort["linker1_edits"]),
                "accepted_linker2_edit_histogram": counter_dict(cohort["linker2_edits"]),
                "read_length_histogram": counter_dict(cohort["read_lengths"]),
            }
        )

    output = {
        "schema_version": "1.0.0",
        "fastq": str(args.fastq),
        "records_scanned": total,
        "distinct_target_reads_inspected": inspected,
        "competitors": [str(path) for path in args.competitor],
        "reference": str(args.reference),
        "criteria": {
            "max_linker1_edit": args.max_linker1_edit,
            "max_linker2_edit": args.max_linker2_edit,
            "max_linker1_candidates": args.max_linker1_candidates,
        },
        "cohorts": output_cohorts,
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    concise = {
        "schema_version": output["schema_version"],
        "records_scanned": total,
        "distinct_target_reads_inspected": inspected,
        "cohorts": {
            label: {
                "emitted": cohort["emitted"],
                "unique": cohort["unique"],
                "unique_reference": cohort["unique_reference"],
                "unique_nonreference": cohort["unique_nonreference"],
                "outcomes": cohort["outcomes"],
            }
            for label, cohort in output_cohorts.items()
        },
        "out": str(args.out),
    }
    print(json.dumps(concise, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
