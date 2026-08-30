#!/usr/bin/env python3
"""Generate the deterministic SPLiT-seq PE Hamming-one Matchbox list.

Canonical sequences are written first in lexical order, followed by the
remaining unique radius-one variants in lexical order.  The ordering does not
change membership semantics, but makes Matchbox's linear ``contains`` check
return early for the common canonical-barcode case.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT / "configs/matchbox/r2_r3.txt"
DEFAULT_OUTPUT = ROOT / "configs/diagnostics/splitseq_pe_bc8_ham1.csv"


def read_column(path: Path, expected_header: str) -> list[str]:
    lines = path.read_text().splitlines()
    if not lines or lines[0] != expected_header:
        raise ValueError(f"{path} must have the header {expected_header!r}")
    sequences = [line.strip().upper() for line in lines[1:] if line.strip()]
    if len(sequences) != len(set(sequences)):
        raise ValueError(f"{path} contains duplicate sequences")
    if not sequences or any(len(sequence) != 8 for sequence in sequences):
        raise ValueError(f"{path} must contain nonempty 8-base sequences")
    if any(set(sequence) - set("ACGT") for sequence in sequences):
        raise ValueError(f"{path} contains a non-ACGT sequence")
    return sequences


def expand(canonical: list[str]) -> dict[str, set[str]]:
    owners: dict[str, set[str]] = defaultdict(set)
    for sequence in canonical:
        owners[sequence].add(sequence)
        for position, observed in enumerate(sequence):
            for base in "ACGT":
                if base != observed:
                    variant = sequence[:position] + base + sequence[position + 1 :]
                    owners[variant].add(sequence)
    return owners


def render(canonical: list[str], owners: dict[str, set[str]]) -> str:
    canonical_set = set(canonical)
    canonical_order = sorted(canonical)
    variants = sorted(sequence for sequence in owners if sequence not in canonical_set)
    return "barcode\n" + "\n".join([*canonical_order, *variants]) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail if the existing output differs instead of rewriting it",
    )
    args = parser.parse_args()

    canonical = read_column(args.source.resolve(), "round_23")
    owners = expand(canonical)
    content = render(canonical, owners)
    output = args.output.resolve()
    if args.check:
        if not output.is_file() or output.read_text() != content:
            raise SystemExit(f"generated whitelist is stale: {output}")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(content)

    collisions = sum(len(values) > 1 for values in owners.values())
    print(
        json.dumps(
            {
                "source": str(args.source.resolve()),
                "output": str(output),
                "canonical_sequences": len(canonical),
                "expanded_sequences": len(owners),
                "ambiguous_sequences": collisions,
                "ordering": "canonical-lexical-then-noncanonical-lexical",
                "checked": args.check,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
