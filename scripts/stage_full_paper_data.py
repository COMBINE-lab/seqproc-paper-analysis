#!/usr/bin/env python3
"""Download, verify, decompress, and characterize the full paper FASTQs."""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import hashlib
import json
import os
import shutil
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence
from urllib.parse import urlparse
from urllib.request import Request, urlopen

from stage_paper_subsets import verify_archive_metadata


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "benchmark_specs" / "data" / "paper_accessions.json"
DEFAULT_DEST = Path("/scratch1/seqproc-benchmark-data/full")
OUTPUT_NAMES = {
    ("splitseq_pe", 1): "SRR6750041_R1.fastq",
    ("splitseq_pe", 2): "SRR6750041_R2.fastq",
    ("lr_splitseq", 1): "SRR13948564_full.fastq",
    ("10x_short", 1): "SRR8315379_R1.fastq",
    ("10x_short", 2): "SRR8315379_R2.fastq",
    ("sciseq", 1): "SRR7827254_1.fastq",
    ("sciseq", 2): "SRR7827254_2.fastq",
}


class StagingError(RuntimeError):
    pass


USER_AGENT = "seqproc-paper-full-data/1.0"
MIN_RANGE_BYTES = 16 * 1024 * 1024


def file_digest(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def normalized_name(header: bytes) -> bytes:
    value = header[1:].split(None, 1)[0]
    if value.endswith((b"/1", b"/2")):
        value = value[:-2]
    return value


def median_from_counts(lengths: Counter[int], records: int) -> float | int:
    if records == 0:
        return 0
    targets = ((records - 1) // 2, records // 2)
    values: list[int] = []
    cumulative = 0
    target_index = 0
    for length, count in sorted(lengths.items()):
        cumulative += count
        while target_index < 2 and targets[target_index] < cumulative:
            values.append(length)
            target_index += 1
    value = (values[0] + values[1]) / 2
    return int(value) if value.is_integer() else value


def byte_ranges(total_bytes: int, connections: int) -> list[tuple[int, int]]:
    if total_bytes <= 0:
        raise ValueError("total_bytes must be positive")
    if connections <= 0:
        raise ValueError("connections must be positive")
    count = min(connections, max(1, (total_bytes + MIN_RANGE_BYTES - 1) // MIN_RANGE_BYTES))
    chunk = (total_bytes + count - 1) // count
    return [
        (start, min(total_bytes - 1, start + chunk - 1))
        for start in range(0, total_bytes, chunk)
    ]


def download_range(url: str, path: Path, start: int, end: int) -> None:
    expected = end - start + 1
    observed = path.stat().st_size if path.is_file() else 0
    if observed > expected:
        raise StagingError(f"range part is larger than expected: {path}")
    if observed == expected:
        return

    for attempt in range(1, 9):
        observed = path.stat().st_size if path.is_file() else 0
        if observed == expected:
            return
        request_start = start + observed
        request = Request(
            url,
            headers={
                "Range": f"bytes={request_start}-{end}",
                "User-Agent": USER_AGENT,
                "Accept-Encoding": "identity",
            },
        )
        try:
            with urlopen(request, timeout=120) as response:
                status = getattr(response, "status", None)
                content_range = response.headers.get("Content-Range", "")
                expected_prefix = f"bytes {request_start}-{end}/"
                if status != 206 or not content_range.startswith(expected_prefix):
                    raise StagingError(
                        f"server did not honor range {request_start}-{end}: "
                        f"status={status}, Content-Range={content_range!r}"
                    )
                with path.open("ab") as sink:
                    shutil.copyfileobj(response, sink, length=8 * 1024 * 1024)
            observed = path.stat().st_size
            if observed > expected:
                raise StagingError(f"range part overflow: {path}")
            if observed == expected:
                return
        except Exception:
            if attempt == 8:
                raise
            time.sleep(min(2 ** (attempt - 1), 30))
    raise AssertionError("unreachable")


def download_archive(
    source: dict[str, Any], destination: Path, connections: int = 8
) -> dict[str, Any]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    expected_bytes = int(source["source_bytes"])
    expected_md5 = str(source["source_md5"])
    if destination.is_file():
        if destination.stat().st_size != expected_bytes:
            raise StagingError(f"existing archive has wrong size: {destination}")
        observed_md5 = file_digest(destination, "md5")
        if observed_md5 != expected_md5:
            raise StagingError(f"existing archive has wrong MD5: {destination}")
        return {
            "path": str(destination),
            "bytes": expected_bytes,
            "md5": observed_md5,
            "sha256": file_digest(destination, "sha256"),
            "retrieval": "verified_existing",
        }

    partial = destination.with_suffix(destination.suffix + ".partial")
    parts_dir = destination.with_name(destination.name + ".parts")
    parts_dir.mkdir(parents=True, exist_ok=True)
    ranges = byte_ranges(expected_bytes, connections)
    part_paths = [parts_dir / f"part-{index:03d}" for index in range(len(ranges))]
    with concurrent.futures.ThreadPoolExecutor(max_workers=len(ranges)) as executor:
        futures = [
            executor.submit(
                download_range,
                str(source["url"]),
                part_path,
                start,
                end,
            )
            for part_path, (start, end) in zip(part_paths, ranges)
        ]
        for future in concurrent.futures.as_completed(futures):
            future.result()

    md5_digest = hashlib.md5()
    sha256_digest = hashlib.sha256()
    with partial.open("wb") as sink:
        for part_path in part_paths:
            with part_path.open("rb") as source_handle:
                while chunk := source_handle.read(8 * 1024 * 1024):
                    sink.write(chunk)
                    md5_digest.update(chunk)
                    sha256_digest.update(chunk)
    if partial.stat().st_size != expected_bytes:
        raise StagingError(
            f"download size mismatch for {partial}: {partial.stat().st_size} != {expected_bytes}"
        )
    observed_md5 = md5_digest.hexdigest()
    if observed_md5 != expected_md5:
        raise StagingError(
            f"download MD5 mismatch for {partial}: {observed_md5} != {expected_md5}"
        )
    compressed_sha256 = sha256_digest.hexdigest()
    os.replace(partial, destination)
    for part_path in part_paths:
        part_path.unlink()
    parts_dir.rmdir()
    return {
        "path": str(destination),
        "bytes": expected_bytes,
        "md5": observed_md5,
        "sha256": compressed_sha256,
        "retrieval": "downloaded_multipart",
        "downloader": "seqproc byte-range downloader v1",
        "range_connections": len(ranges),
        "ranges": [[start, end] for start, end in ranges],
    }


def decompress_and_characterize(
    archive: Path, output: Path, expected: dict[str, Any] | None = None
) -> dict[str, Any]:
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        if expected is None:
            raise StagingError(
                f"untracked decompressed output already exists: {output}; provide its provenance or move it"
            )
        if output.stat().st_size != int(expected["bytes"]):
            raise StagingError(f"existing FASTQ has wrong size: {output}")
        observed = file_digest(output, "sha256")
        if observed != expected["sha256"]:
            raise StagingError(f"existing FASTQ has wrong SHA-256: {output}")
        return dict(expected)

    partial = output.with_suffix(output.suffix + ".partial")
    digest = hashlib.sha256()
    names_digest = hashlib.sha256()
    lengths: Counter[int] = Counter()
    records = 0
    sequence_bases = 0
    first_id: bytes | None = None
    last_id: bytes | None = None
    started = time.monotonic()
    with gzip.open(archive, "rb") as source, partial.open("wb") as sink:
        while True:
            header = source.readline()
            if not header:
                break
            sequence = source.readline()
            plus = source.readline()
            quality = source.readline()
            if not sequence or not plus or not quality:
                raise StagingError(f"truncated FASTQ record {records + 1} in {archive}")
            if not header.startswith(b"@") or not plus.startswith(b"+"):
                raise StagingError(f"malformed FASTQ record {records + 1} in {archive}")
            sequence_value = sequence.rstrip(b"\r\n")
            quality_value = quality.rstrip(b"\r\n")
            if len(sequence_value) != len(quality_value):
                raise StagingError(
                    f"sequence/quality mismatch at record {records + 1} in {archive}"
                )
            name = normalized_name(header)
            if first_id is None:
                first_id = name
            last_id = name
            names_digest.update(name)
            names_digest.update(b"\0")
            length = len(sequence_value)
            lengths[length] += 1
            sequence_bases += length
            for line in (header, sequence, plus, quality):
                sink.write(line)
                digest.update(line)
            records += 1
            if records % 5_000_000 == 0:
                print(f"  {output.name}: {records:,} records", flush=True)
    os.replace(partial, output)
    return {
        "path": str(output),
        "bytes": output.stat().st_size,
        "sha256": digest.hexdigest(),
        "records": records,
        "sequence_bases": sequence_bases,
        "length_summary": {
            "minimum": min(lengths) if lengths else 0,
            "median": median_from_counts(lengths, records),
            "maximum": max(lengths) if lengths else 0,
            "distinct_lengths": len(lengths),
        },
        "first_normalized_id": (first_id or b"").decode("utf-8", errors="replace"),
        "last_normalized_id": (last_id or b"").decode("utf-8", errors="replace"),
        "normalized_names_sha256": names_digest.hexdigest(),
        "decompression_wall_seconds": time.monotonic() - started,
    }


def write_provenance(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--dest", type=Path, default=DEFAULT_DEST)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument(
        "--connections-per-file",
        type=int,
        default=8,
        help="simultaneous HTTP byte ranges per archive (default: 8)",
    )
    parser.add_argument("--only", action="append", help="dataset name; repeat to select")
    parser.add_argument(
        "--skip-live-metadata-check",
        action="store_true",
        help="use the checked-in archive manifest without querying ENA first",
    )
    args = parser.parse_args(argv)
    if args.jobs <= 0:
        parser.error("--jobs must be positive")
    if args.connections_per_file <= 0:
        parser.error("--connections-per-file must be positive")

    manifest_path = args.manifest.resolve()
    manifest = json.loads(manifest_path.read_text())
    selected = set(args.only or [])
    datasets = [
        item for item in manifest["datasets"] if not selected or item["name"] in selected
    ]
    unknown = selected - {item["name"] for item in manifest["datasets"]}
    if unknown:
        parser.error(f"unknown datasets: {', '.join(sorted(unknown))}")

    verified_metadata: dict[str, Any] = {}
    if not args.skip_live_metadata_check:
        print("verifying checked archive metadata against ENA", flush=True)
        for dataset in datasets:
            observed = verify_archive_metadata(dataset)
            verified_metadata[str(dataset["name"])] = observed
            print(f"  verified {dataset['accession']}", flush=True)
    dest = args.dest.resolve()
    compressed_root = dest / "compressed"
    fastq_root = dest / "fastq"
    provenance_path = dest / "full_data_provenance.json"
    previous = json.loads(provenance_path.read_text()) if provenance_path.is_file() else {}
    prior_files = {
        (item["dataset"], int(item["lane"])): item
        for item in previous.get("files", [])
    }

    work: list[dict[str, Any]] = []
    for dataset in datasets:
        for source in dataset["files"]:
            key = (str(dataset["name"]), int(source["lane"]))
            work.append(
                {
                    "dataset": dataset,
                    "source": source,
                    "key": key,
                    "archive": compressed_root / Path(urlparse(source["url"]).path).name,
                    "output": fastq_root / OUTPUT_NAMES[key],
                }
            )

    archives: dict[tuple[str, int], dict[str, Any]] = {}
    print(f"downloading/verifying {len(work)} archives with {args.jobs} workers", flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(
                download_archive,
                item["source"],
                item["archive"],
                args.connections_per_file,
            ): item
            for item in work
        }
        for future in concurrent.futures.as_completed(futures):
            item = futures[future]
            archives[item["key"]] = future.result()
            print(f"  verified {item['archive'].name}", flush=True)

    processed_names = {str(item["name"]) for item in datasets}
    records: list[dict[str, Any]] = [
        item for item in previous.get("files", []) if item.get("dataset") not in processed_names
    ]
    generated_at = datetime.now(timezone.utc).isoformat()

    def checkpoint(*, complete: bool, paired_validation: list[dict[str, Any]]) -> None:
        ordered = sorted(records, key=lambda item: (item["dataset"], int(item["lane"])))
        write_provenance(
            provenance_path,
            {
                "schema_version": "1.0.0",
                "complete": complete,
                "generated_at_utc": generated_at,
                "manifest": {
                    "path": str(manifest_path),
                    "sha256": file_digest(manifest_path, "sha256"),
                },
                "live_ena_metadata": verified_metadata,
                "files": ordered,
                "paired_validation": paired_validation,
            },
        )

    print(f"decompressing/characterizing with {args.jobs} workers", flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {}
        for item in work:
            previous_file = prior_files.get(item["key"])
            expected = previous_file.get("fastq") if previous_file else None
            future = executor.submit(
                decompress_and_characterize, item["archive"], item["output"], expected
            )
            futures[future] = item
        for future in concurrent.futures.as_completed(futures):
            item = futures[future]
            dataset = item["dataset"]
            source = item["source"]
            fastq = future.result()
            record = {
                "dataset": dataset["name"],
                "accession": dataset["accession"],
                "lane": source["lane"],
                "source_url": source["url"],
                "source_md5": source["source_md5"],
                "source_bytes": source["source_bytes"],
                "archive": archives[item["key"]],
                "fastq": fastq,
            }
            records.append(record)
            checkpoint(complete=False, paired_validation=[])
            print(f"  characterized {item['output'].name}: {fastq['records']:,} records", flush=True)

    records.sort(key=lambda item: (item["dataset"], int(item["lane"])))
    paired_validation = []
    for dataset in datasets:
        lanes = [item for item in records if item["dataset"] == dataset["name"]]
        if len(lanes) != 2:
            continue
        left, right = lanes
        equal_count = left["fastq"]["records"] == right["fastq"]["records"]
        equal_names = (
            left["fastq"]["normalized_names_sha256"]
            == right["fastq"]["normalized_names_sha256"]
        )
        paired_validation.append(
            {
                "dataset": dataset["name"],
                "record_counts_equal": equal_count,
                "normalized_id_streams_equal": equal_names,
            }
        )
        if not equal_count or not equal_names:
            raise StagingError(f"paired FASTQ mismatch for {dataset['name']}")

    checkpoint(complete=True, paired_validation=paired_validation)
    print(f"wrote {provenance_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
