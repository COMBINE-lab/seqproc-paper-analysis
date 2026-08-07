#!/usr/bin/env python3
"""Create and document a reverse-complement FASTQ outside benchmark timing."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from collections.abc import Sequence
from pathlib import Path
from typing import BinaryIO

IUPAC = b"ACGTRYKMSWBDHVNacgtrykmswbdhvn"
IUPAC_COMPLEMENT = b"TGCAYRMKSWVHDBNtgcayrmkswvhdbn"
COMPLEMENT = bytes.maketrans(IUPAC, IUPAC_COMPLEMENT)


def write_hashed(handle: BinaryIO, digest: object, value: bytes) -> int:
    handle.write(value)
    digest.update(value)
    return len(value)


def stage(
    input_path: Path, output_path: Path, provenance_path: Path
) -> dict[str, object]:
    if output_path.exists() or provenance_path.exists():
        raise FileExistsError(
            "refusing to overwrite reverse-complement output or provenance"
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    provenance_path.parent.mkdir(parents=True, exist_ok=True)
    partial_output = output_path.with_name(output_path.name + ".partial")
    partial_provenance = provenance_path.with_name(provenance_path.name + ".partial")
    if partial_output.exists() or partial_provenance.exists():
        raise FileExistsError(
            "partial reverse-complement staging output already exists"
        )

    source_sha256 = hashlib.sha256()
    output_sha256 = hashlib.sha256()
    source_bytes = 0
    output_bytes = 0
    records = 0
    started = time.monotonic()
    try:
        with input_path.open("rb") as source, partial_output.open("wb") as sink:
            while True:
                header = source.readline()
                if not header:
                    break
                sequence = source.readline()
                plus = source.readline()
                quality = source.readline()
                if not sequence or not plus or not quality:
                    raise ValueError(f"truncated FASTQ at record {records + 1}")
                for value in (header, sequence, plus, quality):
                    source_sha256.update(value)
                    source_bytes += len(value)
                if not header.startswith(b"@") or not plus.startswith(b"+"):
                    raise ValueError(f"malformed FASTQ at record {records + 1}")
                seq = sequence.rstrip(b"\r\n")
                qual = quality.rstrip(b"\r\n")
                if len(seq) != len(qual):
                    raise ValueError(
                        f"sequence/quality length mismatch at record {records + 1}"
                    )
                output_bytes += write_hashed(
                    sink, output_sha256, header.rstrip(b"\r\n") + b"\n"
                )
                output_bytes += write_hashed(
                    sink, output_sha256, seq.translate(COMPLEMENT)[::-1] + b"\n"
                )
                output_bytes += write_hashed(
                    sink, output_sha256, plus.rstrip(b"\r\n") + b"\n"
                )
                output_bytes += write_hashed(sink, output_sha256, qual[::-1] + b"\n")
                records += 1
            sink.flush()
            os.fsync(sink.fileno())

        payload: dict[str, object] = {
            "schema_version": "1.0.0",
            "transformation": "reverse-complement FASTQ sequence and reverse quality",
            "header_preserved": True,
            "plus_line_preserved": True,
            "reconciliation_performed": False,
            "records": records,
            "source": {
                "path": str(input_path.resolve()),
                "bytes": source_bytes,
                "sha256": source_sha256.hexdigest(),
            },
            "output": {
                "path": str(output_path.resolve()),
                "bytes": output_bytes,
                "sha256": output_sha256.hexdigest(),
            },
            "elapsed_seconds": time.monotonic() - started,
        }
        partial_provenance.write_text(json.dumps(payload, indent=2) + "\n")
        os.replace(partial_output, output_path)
        os.replace(partial_provenance, provenance_path)
        return payload
    except Exception:
        partial_output.unlink(missing_ok=True)
        partial_provenance.unlink(missing_ok=True)
        raise


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    args = parser.parse_args(argv)
    if not args.input.is_file():
        parser.error(f"input FASTQ does not exist: {args.input}")
    payload = stage(args.input, args.output, args.provenance)
    print(json.dumps(payload, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
