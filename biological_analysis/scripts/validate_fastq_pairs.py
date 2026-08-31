#!/usr/bin/env python3
"""Validate paired FASTQ products and emit a provenance-rich JSON report.

Validation is intentionally streaming: it checks complete four-line records,
sequence/quality lengths, paired identifiers, constant expected read lengths,
record counts, byte sizes, and SHA-256 digests without loading read IDs into
memory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def normalized_id(header: bytes) -> bytes:
    token = header.split(None, 1)[0]
    if token.endswith((b"/1", b"/2")):
        token = token[:-2]
    return token


def read_record(handle, digest: hashlib._Hash, path: Path, number: int):
    lines = [handle.readline() for _ in range(4)]
    if not lines[0]:
        if any(lines[1:]):
            raise ValueError(f"partial FASTQ record at EOF in {path}")
        return None
    if any(not line for line in lines[1:]):
        raise ValueError(f"truncated FASTQ record {number:,} in {path}")
    for line in lines:
        digest.update(line)
    header, sequence, plus, quality = (line.rstrip(b"\r\n") for line in lines)
    if not header.startswith(b"@"):
        raise ValueError(f"record {number:,} in {path} lacks an @ header")
    if not plus.startswith(b"+"):
        raise ValueError(f"record {number:,} in {path} lacks a + separator")
    if len(sequence) != len(quality):
        raise ValueError(
            f"sequence/quality length mismatch in record {number:,} of {path}"
        )
    return header[1:], len(sequence)


def validate_pair(r1: Path, r2: Path, expected_r1: int | None, expected_r2: int | None):
    d1, d2 = hashlib.sha256(), hashlib.sha256()
    count = 0
    lengths1: set[int] = set()
    lengths2: set[int] = set()
    first_id = last_id = None
    with r1.open("rb", buffering=8 * 1024 * 1024) as h1, r2.open(
        "rb", buffering=8 * 1024 * 1024
    ) as h2:
        while True:
            number = count + 1
            a = read_record(h1, d1, r1, number)
            b = read_record(h2, d2, r2, number)
            if a is None or b is None:
                if a is not None or b is not None:
                    raise ValueError(f"paired FASTQs differ in record count: {r1}, {r2}")
                break
            id1, len1 = a
            id2, len2 = b
            nid1, nid2 = normalized_id(id1), normalized_id(id2)
            if nid1 != nid2:
                raise ValueError(
                    f"paired identifier mismatch at record {number:,}: "
                    f"{id1.decode(errors='replace')} != {id2.decode(errors='replace')}"
                )
            if first_id is None:
                first_id = nid1.decode(errors="replace")
            last_id = nid1.decode(errors="replace")
            lengths1.add(len1)
            lengths2.add(len2)
            count += 1
    if expected_r1 is not None and lengths1 != {expected_r1}:
        raise ValueError(f"R1 lengths {sorted(lengths1)} differ from expected {expected_r1}")
    if expected_r2 is not None and lengths2 != {expected_r2}:
        raise ValueError(f"R2 lengths {sorted(lengths2)} differ from expected {expected_r2}")
    return {
        "records": count,
        "first_id": first_id,
        "last_id": last_id,
        "r1": {
            "path": str(r1.resolve()),
            "bytes": r1.stat().st_size,
            "sha256": d1.hexdigest(),
            "lengths": sorted(lengths1),
        },
        "r2": {
            "path": str(r2.resolve()),
            "bytes": r2.stat().st_size,
            "sha256": d2.hexdigest(),
            "lengths": sorted(lengths2),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--r1", type=Path, required=True)
    parser.add_argument("--r2", type=Path, required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-r1-length", type=int)
    parser.add_argument("--expected-r2-length", type=int)
    args = parser.parse_args()
    result = {
        "schema_version": 1,
        "tool": args.name,
        **validate_pair(
            args.r1,
            args.r2,
            args.expected_r1_length,
            args.expected_r2_length,
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
