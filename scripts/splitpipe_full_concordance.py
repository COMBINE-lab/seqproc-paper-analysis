#!/usr/bin/env python3
"""Evaluate final SPLiT-seq PE outputs against a fresh split-pipe run.

This is an agreement analysis against an independent vendor implementation,
not a biological ground truth analysis. It converts split-pipe's
``process/barcode_head.fastq`` to a compact numeric-accession bitmap, verifies
the final publication bitmaps for seqproc, Matchbox, and splitcode, and emits
provenance-rich JSON and CSV artifacts.

For calibration, ``--archived-vendor-ids`` compares the first N records of the
fresh split-pipe set to the archived 10-million-pair accepted-ID set. A valid
calibration should have a zero symmetric difference.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import verify_vendor_concordance_pe as common


ROOT = Path(__file__).resolve().parent.parent
DEFAULT_ACCURACY_JSON = (
    ROOT
    / "publication_results"
    / "journal_rerun_2026-08-17"
    / "splitseq_pe_accuracy_metrics.json"
)
DEFAULT_OUTPUT_PREFIX = (
    ROOT
    / "publication_results"
    / "journal_rerun_2026-08-17"
    / "splitseq_pe_splitpipe_vendor_full"
)
TOOLS = ("seqproc", "matchbox", "splitcode")
ORIGINAL_HEADER_MARKER = b"__OH_@"


def splitpipe_fastq_bitmap(
    path: Path, input_records: int, expected_accession: str | None
) -> tuple[bytes, dict[str, object]]:
    """Parse and validate split-pipe FASTQ into a one-based accession bitmap."""
    bitmap = bytearray((input_records + 7) // 8)
    digest = hashlib.sha256()
    accession_prefix: bytes | None = None
    emitted_records = 0
    minimum = input_records + 1
    maximum = 0

    with path.open("rb") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline()
            separator = handle.readline()
            quality = handle.readline()
            if not sequence or not separator or not quality:
                raise ValueError(f"truncated split-pipe FASTQ after record {emitted_records}")
            for line in (header, sequence, separator, quality):
                digest.update(line)
            if not header.startswith(b"@") or not separator.startswith(b"+"):
                raise ValueError(
                    f"malformed split-pipe FASTQ record {emitted_records + 1}"
                )
            if len(sequence.rstrip(b"\r\n")) != len(quality.rstrip(b"\r\n")):
                raise ValueError(
                    f"sequence/quality length mismatch at record {emitted_records + 1}"
                )
            marker = header.rfind(ORIGINAL_HEADER_MARKER)
            if marker < 0:
                raise ValueError(
                    f"missing original-header marker at record {emitted_records + 1}"
                )
            token = header[marker + len(ORIGINAL_HEADER_MARKER) :].split(None, 1)[0]
            try:
                prefix, digits = token.rsplit(b".", 1)
                numeric_id = int(digits)
            except (ValueError, IndexError) as error:
                raise ValueError(
                    f"invalid original read ID at record {emitted_records + 1}: {token!r}"
                ) from error
            if not prefix or not 1 <= numeric_id <= input_records:
                raise ValueError(
                    f"read ID outside 1..{input_records} at record "
                    f"{emitted_records + 1}: {token!r}"
                )
            if accession_prefix is None:
                accession_prefix = prefix
                if expected_accession is not None and prefix.decode() != expected_accession:
                    raise ValueError(
                        f"split-pipe accession {prefix.decode()} differs from "
                        f"expected {expected_accession}"
                    )
            elif prefix != accession_prefix:
                raise ValueError(
                    f"mixed accessions at split-pipe record {emitted_records + 1}"
                )
            byte, bit = divmod(numeric_id - 1, 8)
            mask = 1 << bit
            if bitmap[byte] & mask:
                raise ValueError(f"duplicate split-pipe read ID: {token.decode()}")
            bitmap[byte] |= mask
            emitted_records += 1
            minimum = min(minimum, numeric_id)
            maximum = max(maximum, numeric_id)

    if emitted_records == 0 or accession_prefix is None:
        raise ValueError(f"empty split-pipe FASTQ: {path}")
    return bytes(bitmap), {
        "path": str(path.resolve()),
        "sha256": digest.hexdigest(),
        "bytes": path.stat().st_size,
        "accepted_records": emitted_records,
        "accepted_fraction": emitted_records / input_records,
        "accession_prefix": accession_prefix.decode(),
        "minimum_numeric_id": minimum,
        "maximum_numeric_id": maximum,
        "fastq_validated": True,
    }


def optional_file(path: Path | None, include_json: bool = False) -> dict[str, object] | None:
    if path is None:
        return None
    resolved = path.resolve()
    result: dict[str, object] = {
        "path": str(resolved),
        "sha256": common.sha256(resolved),
        "bytes": resolved.stat().st_size,
    }
    if include_json:
        result["contents"] = json.loads(resolved.read_text())
    return result


def write_csv(path: Path, result: dict[str, object]) -> None:
    fields = [
        "dataset",
        "comparator",
        "tool",
        "input_records",
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
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for tool in TOOLS:
            writer.writerow(
                {
                    "dataset": result["dataset"],
                    "comparator": "fresh_splitpipe_vendor_accepted_set",
                    "tool": tool,
                    "input_records": result["input_records"],
                    "vendor_records": result["splitpipe"]["accepted_records"],
                    **result["tools"][tool]["metrics"],
                }
            )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--splitpipe-fastq", type=Path, required=True)
    parser.add_argument("--input-records", type=int, required=True)
    parser.add_argument("--accession", default="SRR6750041")
    parser.add_argument("--accuracy-json", type=Path, default=DEFAULT_ACCURACY_JSON)
    parser.add_argument("--bitmap-dir", type=Path)
    parser.add_argument("--archived-vendor-ids", type=Path)
    parser.add_argument("--archived-records", type=int, default=10_000_000)
    parser.add_argument("--input-r1", type=Path)
    parser.add_argument("--input-r2", type=Path)
    parser.add_argument("--splitpipe-run-def", type=Path)
    parser.add_argument("--splitpipe-log", type=Path)
    parser.add_argument("--splitpipe-config", type=Path)
    parser.add_argument("--container-image-id")
    parser.add_argument("--container-image-digest")
    parser.add_argument("--output-prefix", type=Path, default=DEFAULT_OUTPUT_PREFIX)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.input_records <= 0:
        raise ValueError("--input-records must be positive")
    if (args.input_r1 is None) != (args.input_r2 is None):
        raise ValueError("provide both --input-r1 and --input-r2, or neither")
    if args.archived_records <= 0 or args.archived_records > args.input_records:
        raise ValueError("--archived-records must be in 1..--input-records")

    accuracy_path = args.accuracy_json.resolve()
    accuracy = json.loads(accuracy_path.read_text())
    vendor, splitpipe = splitpipe_fastq_bitmap(
        args.splitpipe_fastq.resolve(), args.input_records, args.accession
    )

    output_prefix = args.output_prefix.resolve()
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    vendor_bitmap_path = output_prefix.with_suffix(".vendor.raw")
    vendor_bitmap_path.write_bytes(vendor)
    splitpipe["bitmap"] = {
        "path": str(vendor_bitmap_path),
        "sha256": common.sha256(vendor_bitmap_path),
        "bytes": vendor_bitmap_path.stat().st_size,
        "format": "dense one-based little-bit-order bitmap",
        "numeric_id_max": args.input_records,
    }

    result: dict[str, object] = {
        "schema_version": "1.0.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "comparison": "agreement with a fresh split-pipe vendor accepted-read set",
        "interpretation": (
            "The split-pipe set is an independent vendor comparator, not biological "
            "ground truth. Reads outside either set remain biologically unadjudicated."
        ),
        "dataset": "splitseq_pe",
        "accession": args.accession,
        "input_records": args.input_records,
        "splitpipe": splitpipe,
        "splitpipe_environment": {
            "container_image_id": args.container_image_id,
            "container_image_digest": args.container_image_digest,
            "run_definition": optional_file(args.splitpipe_run_def, include_json=True),
            "log": optional_file(args.splitpipe_log),
            "config": optional_file(args.splitpipe_config),
        },
        "inputs": {},
        "accuracy_campaign": {
            "summary_path": str(accuracy_path),
            "summary_sha256": common.sha256(accuracy_path),
            "campaign": accuracy.get("campaign"),
            "manifest": accuracy.get("manifest"),
        },
        "program": {
            "path": str(Path(__file__).resolve()),
            "sha256": common.sha256(Path(__file__).resolve()),
            "python": sys.version.split()[0],
        },
        "tools": {},
    }
    if args.input_r1 is not None:
        result["inputs"] = {
            "r1": optional_file(args.input_r1),
            "r2": optional_file(args.input_r2),
        }

    metric_rows = {row["tool"]: row for row in accuracy["metrics"]}
    for tool in TOOLS:
        path, provenance = common.tool_bitmap_path(accuracy, tool, args.bitmap_dir)
        path = path.resolve()
        digest = common.sha256(path)
        sources = provenance.get("sources", [])
        if not sources or digest != sources[0].get("bitmap_sha256"):
            raise ValueError(f"{tool} bitmap does not match its accuracy provenance")
        metadata, emitted = common.load_canonical_bitmap(path)
        if metadata["numeric_id_max"] < args.input_records:
            raise ValueError(
                f"{tool} bitmap range {metadata['numeric_id_max']} does not cover "
                f"input range {args.input_records}"
            )
        if metadata["accession_prefix"] != args.accession:
            raise ValueError(f"{tool} bitmap accession differs from {args.accession}")
        expected_count = int(metric_rows[tool]["emitted_records"])
        if common.popcount(emitted) != expected_count:
            raise ValueError(f"{tool} bitmap count differs from accuracy summary")
        restricted = common.restrict_bitmap(emitted, args.input_records)
        result["tools"][tool] = {
            "bitmap": {"path": str(path), "sha256": digest, **metadata},
            "condition_id": provenance.get("condition_id"),
            "configs": provenance.get("configs", []),
            "executables": provenance.get("executables", []),
            "repositories": provenance.get("repositories", []),
            "metrics": common.metrics(restricted, vendor, args.input_records),
        }

    if args.archived_vendor_ids is not None:
        archived, archived_provenance = common.load_vendor_bitmap(
            args.archived_vendor_ids.resolve(), args.archived_records
        )
        fresh_subset = common.restrict_bitmap(vendor, args.archived_records)
        comparison = common.metrics(fresh_subset, archived, args.archived_records)
        comparison.update(
            {
                "fresh_only_records": (
                    comparison["emitted_records"]
                    - comparison["intersection_records"]
                ),
                "archived_only_records": (
                    archived_provenance["records"]
                    - comparison["intersection_records"]
                ),
                "identical": fresh_subset == archived,
            }
        )
        result["archived_calibration"] = {
            "records_considered": args.archived_records,
            "archived": archived_provenance,
            "comparison": comparison,
        }

    json_path = output_prefix.with_suffix(".json")
    csv_path = output_prefix.with_suffix(".csv")
    json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    write_csv(csv_path, result)

    print(
        f"split-pipe accepted: {splitpipe['accepted_records']:,} / "
        f"{args.input_records:,} ({splitpipe['accepted_fraction']:.4%})"
    )
    if "archived_calibration" in result:
        cal = result["archived_calibration"]["comparison"]
        print(
            "archived calibration: "
            f"identical={cal['identical']}, fresh_only={cal['fresh_only_records']:,}, "
            f"archived_only={cal['archived_only_records']:,}"
        )
    print(
        f"{'tool':<10} {'emitted':>12} {'intersection':>13} "
        f"{'precision':>10} {'recall':>10} {'f1':>9} {'jaccard':>9}"
    )
    for tool in TOOLS:
        row = result["tools"][tool]["metrics"]
        print(
            f"{tool:<10} {row['emitted_records']:>12,} "
            f"{row['intersection_records']:>13,} {row['precision']:>10.4%} "
            f"{row['recall']:>10.4%} {row['f1']:>9.6f} {row['jaccard']:>9.6f}"
        )
    print(f"JSON:   {json_path}")
    print(f"CSV:    {csv_path}")
    print(f"bitmap: {vendor_bitmap_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
