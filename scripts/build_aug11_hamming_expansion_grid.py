#!/usr/bin/env python3
"""Build deterministic Hamming-one artifacts for the Aug. 11 LR diagnostic grid.

These generated files are sensitivity-analysis inputs, not required seqproc
inputs or proposed user-facing preprocessing.
"""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BC23_SOURCE = ROOT / "configs/seqproc/splitseq_bc23_whitelist.txt"
BC1_SOURCE = ROOT / "configs/seqproc/splitseq_bc1_whitelist_6bp.txt"
BC23_OUTPUT = ROOT / "configs/diagnostics/aug11_lr_splitseq_bc23_ham1.csv"
BC1_OUTPUT = ROOT / "configs/diagnostics/aug11_lr_splitseq_bc1_ham1.csv"
SPLITCODE_OUTPUT = (
    ROOT / "configs/splitcode/diagnostic_aug11_lr_splitseq_expanded_exact.txt"
)
METADATA_OUTPUT = (
    ROOT
    / "publication_results/lr_splitseq_hamming_expansion_grid_2026-08-11"
    / "expanded_whitelist_metadata.json"
)

LINKER1 = "GTGGCCGATGTTTCGCATCGGCGTACGACT"
LINKER2 = "ATCCACGTGCTTGAGACTGTGG"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_sequences(path: Path) -> list[str]:
    sequences = []
    seen = set()
    for line in path.read_text().splitlines():
        sequence = line.strip().upper()
        if not sequence or sequence in seen:
            continue
        if set(sequence) - set("ACGT"):
            raise ValueError(f"invalid sequence in {path}: {sequence!r}")
        seen.add(sequence)
        sequences.append(sequence)
    if not sequences:
        raise ValueError(f"empty whitelist: {path}")
    lengths = {len(sequence) for sequence in sequences}
    if len(lengths) != 1:
        raise ValueError(f"mixed sequence lengths in {path}: {sorted(lengths)}")
    return sequences


def expand_hamming_one(canonical: list[str]) -> dict[str, list[str]]:
    owners: dict[str, set[str]] = defaultdict(set)
    for sequence in canonical:
        owners[sequence].add(sequence)
        for index, current in enumerate(sequence):
            for base in "ACGT":
                if base != current:
                    variant = sequence[:index] + base + sequence[index + 1 :]
                    owners[variant].add(sequence)
    return {
        variant: sorted(canonical_owners)
        for variant, canonical_owners in sorted(owners.items())
    }


def expansion_metadata(
    source: Path, canonical: list[str], owners: dict[str, list[str]]
) -> dict[str, object]:
    collisions = {
        variant: values for variant, values in owners.items() if len(values) > 1
    }
    length = len(canonical[0])
    return {
        "source": str(source.relative_to(ROOT)),
        "source_sha256": sha256(source),
        "sequence_length": length,
        "canonical_records_in_source": len(source.read_text().splitlines()),
        "unique_canonical_sequences": len(canonical),
        "raw_hamming_entries_before_deduplication": len(canonical) * (1 + 3 * length),
        "unique_expanded_sequences": len(owners),
        "ambiguous_expanded_sequences": len(collisions),
        "maximum_canonical_multiplicity": max(map(len, owners.values())),
        "collision_owners": collisions,
    }


def write_csv(path: Path, header: str, owners: dict[str, list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(header + "\n" + "\n".join(owners) + "\n")


def write_splitcode(
    path: Path, bc23: dict[str, list[str]], bc1: dict[str, list[str]]
) -> None:
    lines = [
        "# Aug. 11 LR-SPLiT-seq Hamming-expansion diagnostic.",
        "# Expanded sequences are deduplicated; all collisions share one logical ID.",
        "# Exact expanded tags reproduce accept-on-ambiguity filtering without",
        "# requiring splitcode's approximate tag generation.",
        "ID\tTAG\tGROUP\tDISTANCE\tNEXT\tPREVIOUS\tMINFINDSG\tLOCATION",
        (
            f"linker1\t{LINKER1}\tLINKER1\t3\t"
            "{{BC2}}0-0\t{{BC3}}0-0\t1\t0:18"
        ),
        (
            f"linker2\t{LINKER2}\tLINKER2\t3\t"
            "{{BC1}}0-0\t{{BC2}}0-0\t1\t-"
        ),
    ]
    lines.extend(
        f"bc3\t*{sequence}\tBC3\t0:0:0\t{{linker1}}0-0\t-\t1\t-"
        for sequence in bc23
    )
    lines.extend(
        (
            f"bc2\t{sequence}\tBC2\t0:0:0\t"
            "{linker2}0-0\t{linker1}0-0\t1\t-"
        )
        for sequence in bc23
    )
    lines.extend(
        f"bc1\t{sequence}\tBC1\t0:0:0\t-\t{{linker2}}0-0\t1\t-"
        for sequence in bc1
    )
    lines.extend(
        [
            "@no-chain",
            "@extract\t<prefix[18]>{linker1}",
            "@extract\t{linker1}<bc2[8]>",
            "@extract\t{linker2}<bc1[6]>",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    bc23_canonical = canonical_sequences(BC23_SOURCE)
    bc1_canonical = canonical_sequences(BC1_SOURCE)
    bc23 = expand_hamming_one(bc23_canonical)
    bc1 = expand_hamming_one(bc1_canonical)

    write_csv(BC23_OUTPUT, "round_23", bc23)
    write_csv(BC1_OUTPUT, "barcode", bc1)
    write_splitcode(SPLITCODE_OUTPUT, bc23, bc1)

    metadata = {
        "schema_version": "1.0.0",
        "purpose": "diagnostic sensitivity analysis, not a user-facing requirement",
        "ambiguity_policy": (
            "deduplicate expanded sequences and accept variants with multiple "
            "canonical owners"
        ),
        "bc23": expansion_metadata(BC23_SOURCE, bc23_canonical, bc23),
        "bc1": expansion_metadata(BC1_SOURCE, bc1_canonical, bc1),
        "outputs": {
            "matchbox_bc23": str(BC23_OUTPUT.relative_to(ROOT)),
            "matchbox_bc1": str(BC1_OUTPUT.relative_to(ROOT)),
            "splitcode_config": str(SPLITCODE_OUTPUT.relative_to(ROOT)),
        },
    }
    METADATA_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    METADATA_OUTPUT.write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
