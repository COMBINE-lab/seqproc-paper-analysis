#!/usr/bin/env python3
"""Compare FASTQ read sets whose IDs end in a positive numeric accession index."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
from typing import Sequence


POPCOUNT = tuple(bin(value).count("1") for value in range(256))


def read_bitset(path: Path, expected_records: int) -> tuple[bytearray, int, str]:
    bits = bytearray((expected_records + 7) // 8)
    digest = hashlib.sha256()
    count = 0
    duplicates = 0
    with path.open("rb") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not sequence or not plus or not quality:
                raise ValueError(f"truncated FASTQ: {path}")
            digest.update(header)
            digest.update(sequence)
            digest.update(plus)
            digest.update(quality)
            token = header[1:].split(None, 1)[0]
            try:
                index = int(token.rsplit(b".", 1)[1]) - 1
            except (IndexError, ValueError) as error:
                raise ValueError(f"non-numeric accession read ID {token!r} in {path}") from error
            if not 0 <= index < expected_records:
                raise ValueError(f"read index {index + 1} outside 1..{expected_records}: {path}")
            offset, mask = divmod(index, 8)
            flag = 1 << mask
            if bits[offset] & flag:
                duplicates += 1
            else:
                bits[offset] |= flag
                count += 1
    if duplicates:
        raise ValueError(f"{path} contains {duplicates} duplicate read IDs")
    return bits, count, digest.hexdigest()


def intersection_size(left: bytearray, right: bytearray) -> int:
    return sum(POPCOUNT[a & b] for a, b in zip(left, right))


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records", type=int, required=True)
    parser.add_argument("--input", action="append", required=True, metavar="LABEL=FASTQ")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    if args.records <= 0:
        parser.error("--records must be positive")
    entries = []
    seen = set()
    for value in args.input:
        if "=" not in value:
            parser.error(f"--input must be LABEL=FASTQ: {value}")
        label, path_text = value.split("=", 1)
        if not label or label in seen:
            parser.error(f"input label is empty or duplicated: {label!r}")
        seen.add(label)
        path = Path(path_text).resolve()
        bits, count, digest = read_bitset(path, args.records)
        entries.append((label, path, bits, count, digest))
    pairwise = []
    for left, right in itertools.combinations(entries, 2):
        shared = intersection_size(left[2], right[2])
        union = left[3] + right[3] - shared
        pairwise.append(
            {
                "left": left[0],
                "right": right[0],
                "intersection": shared,
                "union": union,
                "jaccard": shared / union if union else 1.0,
                "left_only": left[3] - shared,
                "right_only": right[3] - shared,
            }
        )
    result = {
        "schema_version": "1.0.0",
        "universe_records": args.records,
        "inputs": [
            {
                "label": label,
                "path": str(path),
                "bytes": path.stat().st_size,
                "sha256": digest,
                "unique_read_ids": count,
            }
            for label, path, _, count, digest in entries
        ],
        "pairwise": pairwise,
    }
    content = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(content)
    else:
        print(content, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
