#!/usr/bin/env python3
"""Expand a DNA whitelist to a unique Hamming-radius neighborhood CSV.

The command refuses to choose among multiply owned variants.  This makes an
externally expanded Matchbox sensitivity configuration auditable rather than
silently imposing a first-hit ambiguity policy.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--column", default="barcode")
    parser.add_argument("--distance", type=int, default=1, choices=[0, 1])
    args = parser.parse_args()
    whitelist = [line.strip().upper() for line in args.input.open() if line.strip()]
    if len(whitelist) != len(set(whitelist)):
        raise ValueError("input whitelist contains duplicate sequences")
    lengths = {len(sequence) for sequence in whitelist}
    if len(lengths) != 1:
        raise ValueError("input whitelist contains mixed sequence lengths")

    owners: dict[str, set[str]] = defaultdict(set)
    for sequence in whitelist:
        owners[sequence].add(sequence)
        if args.distance == 1:
            for index, original in enumerate(sequence):
                for base in "ACGT":
                    if base != original:
                        owners[sequence[:index] + base + sequence[index + 1 :]].add(sequence)
    ambiguous = {variant: source for variant, source in owners.items() if len(source) > 1}
    if ambiguous:
        examples = list(sorted(ambiguous.items()))[:5]
        raise ValueError(
            f"{len(ambiguous)} expanded variants have multiple owners; examples: {examples}"
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow([args.column])
        writer.writerows([variant] for variant in sorted(owners))
    print(
        f"wrote {len(owners):,} unique variants from {len(whitelist):,} barcodes "
        f"at Hamming distance <= {args.distance}: {args.output}"
    )


if __name__ == "__main__":
    main()
