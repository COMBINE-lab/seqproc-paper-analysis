#!/usr/bin/env python3
"""Compute read-set intersections from compact campaign accession bitmaps."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

MAGIC = b"fastq_numeric_accession_set_v1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load(path: Path):
    fields = path.read_bytes().split(b"\0", 4)
    if len(fields) != 5 or fields[0] != MAGIC:
        raise ValueError(f"unexpected bitmap format: {path}")
    metadata = {
        "mate": int(fields[1]),
        "numeric_id_max": int(fields[2]),
        "accession_prefix": fields[3].decode(),
        "path": str(path.resolve()),
        "sha256": sha256(path),
    }
    return metadata, int.from_bytes(fields[4], "little")


def named(value: str):
    try:
        name, path = value.split(":", 1)
    except ValueError as error:
        raise argparse.ArgumentTypeError("bitmap must be NAME:PATH") from error
    return name, Path(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("bitmap", nargs="+", type=named)
    args = parser.parse_args()
    loaded = {name: load(path) for name, path in args.bitmap}
    domains = {(meta["numeric_id_max"], meta["accession_prefix"]) for meta, _ in loaded.values()}
    if len(domains) != 1:
        raise ValueError("bitmap accession domains differ")
    names = [name for name, _ in args.bitmap]
    result = {"schema_version": 1, "sources": {}, "tools": {}, "pairs": {}}
    for name in names:
        meta, bits = loaded[name]
        result["sources"][name] = meta
        result["tools"][name] = {"records": bits.bit_count()}
    for a, b in itertools.combinations(names, 2):
        ba, bb = loaded[a][1], loaded[b][1]
        intersection = (ba & bb).bit_count()
        union = (ba | bb).bit_count()
        result["pairs"][f"{a}|{b}"] = {
            "intersection": intersection,
            "union": union,
            "jaccard": intersection / union if union else None,
        }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result["pairs"], indent=2))


if __name__ == "__main__":
    main()
