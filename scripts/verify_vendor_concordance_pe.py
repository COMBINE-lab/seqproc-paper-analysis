#!/usr/bin/env python3
"""Compare the final SPLiT-seq PE tool outputs with the archived split-pipe set.

The archived vendor set contains the read IDs accepted by split-pipe on the
first 10,000,000 pairs of SRR6750041. It is an orthogonal vendor comparator,
not biological ground truth. This program restricts the compact accepted-ID
bitmaps from the final full-data publication campaign to that same numeric-ID
range, computes agreement for seqproc, matchbox, and splitcode, and writes
provenance-rich JSON and CSV artifacts.

The program does not require materialized FASTQ outputs from the three tools.
When subset and full FASTQs are supplied, it additionally verifies that each
subset is a byte-identical prefix of the corresponding full campaign input.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import BinaryIO


ROOT = Path(__file__).resolve().parent.parent
DEFAULT_ACCURACY_JSON = (
    ROOT
    / "publication_results"
    / "journal_rerun_2026-08-17"
    / "splitseq_pe_accuracy_metrics.json"
)
DEFAULT_VENDOR_IDS = ROOT / "results_final" / "splitpipe_valid_ids_10M.txt.gz"
DEFAULT_OUTPUT_PREFIX = (
    ROOT
    / "publication_results"
    / "journal_rerun_2026-08-17"
    / "splitseq_pe_splitpipe_vendor_10m"
)
MAGIC = b"fastq_numeric_accession_set_v1"
TOOLS = ("seqproc", "matchbox", "splitcode")
POPCOUNT = tuple(bin(value).count("1") for value in range(256))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def open_maybe_gzip(path: Path) -> BinaryIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rb")
    return path.open("rb")


def load_vendor_bitmap(path: Path, subset_records: int) -> tuple[bytes, dict[str, object]]:
    """Load unique ACCESSION.NUMBER IDs into a dense, one-based bitmap."""
    bitmap = bytearray((subset_records + 7) // 8)
    accession_prefix: bytes | None = None
    records = 0
    minimum = subset_records + 1
    maximum = 0
    with open_maybe_gzip(path) as handle:
        for line_number, line in enumerate(handle, 1):
            token = line.strip().split(None, 1)[0]
            if not token:
                continue
            try:
                prefix, digits = token.rsplit(b".", 1)
                numeric_id = int(digits)
            except (ValueError, IndexError) as error:
                raise ValueError(
                    f"invalid vendor read ID at {path}:{line_number}: {token!r}"
                ) from error
            if not prefix or not 1 <= numeric_id <= subset_records:
                raise ValueError(
                    f"vendor read ID outside 1..{subset_records} at "
                    f"{path}:{line_number}: {token!r}"
                )
            if accession_prefix is None:
                accession_prefix = prefix
            elif prefix != accession_prefix:
                raise ValueError(
                    f"mixed vendor accession prefixes at {path}:{line_number}"
                )
            byte, bit = divmod(numeric_id - 1, 8)
            mask = 1 << bit
            if bitmap[byte] & mask:
                raise ValueError(f"duplicate vendor read ID at {path}:{line_number}: {token!r}")
            bitmap[byte] |= mask
            records += 1
            minimum = min(minimum, numeric_id)
            maximum = max(maximum, numeric_id)
    if records == 0 or accession_prefix is None:
        raise ValueError(f"empty vendor ID set: {path}")
    return bytes(bitmap), {
        "path": str(path.resolve()),
        "sha256": sha256(path),
        "records": records,
        "accession_prefix": accession_prefix.decode(),
        "minimum_numeric_id": minimum,
        "maximum_numeric_id": maximum,
    }


def load_canonical_bitmap(path: Path) -> tuple[dict[str, object], bytes]:
    fields = path.read_bytes().split(b"\0", 4)
    if len(fields) != 5 or fields[0] != MAGIC:
        raise ValueError(f"unexpected accession bitmap format: {path}")
    metadata = {
        "mate": int(fields[1]),
        "numeric_id_max": int(fields[2]),
        "accession_prefix": fields[3].decode(),
    }
    expected_bytes = (int(metadata["numeric_id_max"]) + 7) // 8
    if len(fields[4]) != expected_bytes:
        raise ValueError(
            f"bitmap payload length differs from numeric ID range for {path}"
        )
    return metadata, fields[4]


def restrict_bitmap(bitmap: bytes, subset_records: int) -> bytes:
    required = (subset_records + 7) // 8
    if len(bitmap) < required:
        raise ValueError("tool bitmap is smaller than the requested subset")
    restricted = bytearray(bitmap[:required])
    residual = subset_records % 8
    if residual:
        restricted[-1] &= (1 << residual) - 1
    return bytes(restricted)


def popcount(bitmap: bytes) -> int:
    return sum(POPCOUNT[value] for value in bitmap)


def metrics(emitted: bytes, vendor: bytes, subset_records: int) -> dict[str, object]:
    if len(emitted) != len(vendor):
        raise ValueError("emitted and vendor bitmaps differ in size")
    emitted_records = popcount(emitted)
    vendor_records = popcount(vendor)
    intersection = sum(POPCOUNT[left & right] for left, right in zip(emitted, vendor))
    union = emitted_records + vendor_records - intersection
    precision = intersection / emitted_records if emitted_records else 0.0
    recall = intersection / vendor_records if vendor_records else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {
        "emitted_records": emitted_records,
        "intersection_records": intersection,
        "union_records": union,
        "emitted_fraction": emitted_records / subset_records,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "jaccard": intersection / union if union else 0.0,
    }


def files_are_prefix_equal(subset: Path, full: Path) -> bool:
    subset_size = subset.stat().st_size
    if full.stat().st_size < subset_size:
        return False
    remaining = subset_size
    with subset.open("rb") as left, full.open("rb") as right:
        while remaining:
            size = min(8 * 1024 * 1024, remaining)
            if left.read(size) != right.read(size):
                return False
            remaining -= size
    return True


def input_provenance(subset: Path, full: Path | None) -> dict[str, object]:
    result: dict[str, object] = {
        "subset_path": str(subset.resolve()),
        "subset_bytes": subset.stat().st_size,
        "subset_sha256": sha256(subset),
    }
    if full is not None:
        result.update(
            {
                "full_path": str(full.resolve()),
                "full_bytes": full.stat().st_size,
                "subset_is_byte_prefix_of_full": files_are_prefix_equal(subset, full),
            }
        )
        if not result["subset_is_byte_prefix_of_full"]:
            raise ValueError(f"{subset} is not a byte-identical prefix of {full}")
    return result


def tool_bitmap_path(
    accuracy: dict[str, object], tool: str, bitmap_dir: Path | None
) -> tuple[Path, dict[str, object]]:
    provenance = accuracy["provenance"][f"splitseq_pe/{tool}"]
    if bitmap_dir is not None:
        path = bitmap_dir / f"splitseq_pe.{tool}.accepted.1.bitmap"
    else:
        sources = provenance["sources"]
        if not sources:
            raise ValueError(f"no bitmap source recorded for {tool}")
        path = Path(sources[0]["bitmap"])
    return path, provenance


def write_csv(path: Path, result: dict[str, object]) -> None:
    fieldnames = [
        "dataset",
        "comparator",
        "tool",
        "subset_input_records",
        "vendor_records",
        "emitted_records",
        "intersection_records",
        "union_records",
        "emitted_fraction",
        "precision",
        "recall",
        "f1",
        "jaccard",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for tool in TOOLS:
            row = result["tools"][tool]["metrics"]
            writer.writerow(
                {
                    "dataset": result["dataset"],
                    "comparator": "archived_splitpipe_vendor_accepted_set",
                    "tool": tool,
                    "subset_input_records": result["subset_input_records"],
                    "vendor_records": result["vendor"]["records"],
                    **row,
                }
            )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accuracy-json", type=Path, default=DEFAULT_ACCURACY_JSON)
    parser.add_argument("--vendor-ids", type=Path, default=DEFAULT_VENDOR_IDS)
    parser.add_argument(
        "--bitmap-dir",
        type=Path,
        help="override bitmap paths recorded in --accuracy-json",
    )
    parser.add_argument("--subset-records", type=int, default=10_000_000)
    parser.add_argument("--subset-r1", type=Path)
    parser.add_argument("--subset-r2", type=Path)
    parser.add_argument("--full-r1", type=Path)
    parser.add_argument("--full-r2", type=Path)
    parser.add_argument("--output-prefix", type=Path, default=DEFAULT_OUTPUT_PREFIX)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.subset_records <= 0:
        raise ValueError("--subset-records must be positive")
    if (args.subset_r1 is None) != (args.subset_r2 is None):
        raise ValueError("provide both --subset-r1 and --subset-r2, or neither")
    if (args.full_r1 is None) != (args.full_r2 is None):
        raise ValueError("provide both --full-r1 and --full-r2, or neither")
    if args.full_r1 is not None and args.subset_r1 is None:
        raise ValueError("full inputs require corresponding subset inputs")

    accuracy_path = args.accuracy_json.resolve()
    accuracy = json.loads(accuracy_path.read_text())
    vendor, vendor_provenance = load_vendor_bitmap(
        args.vendor_ids.resolve(), args.subset_records
    )
    result: dict[str, object] = {
        "schema_version": "2.0.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "comparison": "agreement with archived split-pipe vendor accepted-ID set",
        "interpretation": (
            "The split-pipe set is an orthogonal vendor comparator, not biological "
            "ground truth. Reads outside either set remain biologically unadjudicated."
        ),
        "dataset": "splitseq_pe",
        "accession": "SRR6750041",
        "subset_input_records": args.subset_records,
        "vendor": vendor_provenance,
        "accuracy_campaign": {
            "summary_path": str(accuracy_path),
            "summary_sha256": sha256(accuracy_path),
            "campaign": accuracy.get("campaign"),
            "manifest": accuracy.get("manifest"),
        },
        "program": {
            "path": str(Path(__file__).resolve()),
            "sha256": sha256(Path(__file__).resolve()),
            "python": sys.version.split()[0],
        },
        "inputs": {},
        "tools": {},
    }
    if args.subset_r1 is not None:
        full_r1 = args.full_r1.resolve() if args.full_r1 is not None else None
        full_r2 = args.full_r2.resolve() if args.full_r2 is not None else None
        result["inputs"] = {
            "r1": input_provenance(args.subset_r1.resolve(), full_r1),
            "r2": input_provenance(args.subset_r2.resolve(), full_r2),
        }

    expected_prefix = vendor_provenance["accession_prefix"]
    metric_rows = {row["tool"]: row for row in accuracy["metrics"]}
    for tool in TOOLS:
        path, provenance = tool_bitmap_path(accuracy, tool, args.bitmap_dir)
        path = path.resolve()
        bitmap_digest = sha256(path)
        recorded_sources = provenance.get("sources", [])
        if not recorded_sources:
            raise ValueError(f"no recorded bitmap provenance for {tool}")
        recorded_digest = recorded_sources[0].get("bitmap_sha256")
        if bitmap_digest != recorded_digest:
            raise ValueError(
                f"{tool} bitmap SHA-256 {bitmap_digest} differs from the accuracy "
                f"artifact {recorded_digest}"
            )
        metadata, full_bitmap = load_canonical_bitmap(path)
        if metadata["numeric_id_max"] < args.subset_records:
            raise ValueError(f"{tool} bitmap does not cover the requested subset")
        if metadata["accession_prefix"] != expected_prefix:
            raise ValueError(f"{tool} and vendor accession prefixes differ")
        full_count = popcount(full_bitmap)
        expected_full_count = int(metric_rows[tool]["emitted_records"])
        if full_count != expected_full_count:
            raise ValueError(
                f"{tool} full bitmap count {full_count} differs from accuracy summary "
                f"{expected_full_count}"
            )
        restricted = restrict_bitmap(full_bitmap, args.subset_records)
        result["tools"][tool] = {
            "bitmap": {
                "path": str(path),
                "sha256": bitmap_digest,
                "verified_against_accuracy_artifact": True,
                **metadata,
                "full_emitted_records": full_count,
            },
            "condition_id": provenance.get("condition_id"),
            "configs": provenance.get("configs", []),
            "executables": provenance.get("executables", []),
            "repositories": provenance.get("repositories", []),
            "metrics": metrics(restricted, vendor, args.subset_records),
        }

    output_prefix = args.output_prefix.resolve()
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    json_path = output_prefix.with_suffix(".json")
    csv_path = output_prefix.with_suffix(".csv")
    json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    write_csv(csv_path, result)

    print(f"split-pipe vendor accepted IDs: {vendor_provenance['records']:,}")
    print(
        f"{'tool':<10} {'emitted':>11} {'intersection':>13} "
        f"{'precision':>10} {'recall':>10} {'f1':>9} {'jaccard':>9}"
    )
    for tool in TOOLS:
        row = result["tools"][tool]["metrics"]
        print(
            f"{tool:<10} {row['emitted_records']:>11,} "
            f"{row['intersection_records']:>13,} {row['precision']:>10.4%} "
            f"{row['recall']:>10.4%} {row['f1']:>9.6f} {row['jaccard']:>9.6f}"
        )
    print(f"JSON: {json_path}")
    print(f"CSV:  {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
