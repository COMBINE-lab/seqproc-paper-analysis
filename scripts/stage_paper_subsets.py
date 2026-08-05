#!/usr/bin/env python3
"""Stream reproducible FASTQ prefixes for the four manuscript accessions.

The script reads gzip-compressed archive files directly from ENA and stops
after the declared number of complete FASTQ records. It therefore avoids the
much larger transient full-run downloads while recording the immutable archive
MD5/size metadata and SHA-256 digests of every generated subset.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "benchmark_specs" / "data" / "paper_accessions.json"
ENA_REPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"
USER_AGENT = "seqproc-paper-reproduction/1.0"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def normalized_name(header: bytes) -> bytes:
    value = header[1:].split(None, 1)[0]
    if value.endswith((b"/1", b"/2")):
        value = value[:-2]
    return value


def ena_metadata(accession: str) -> dict[str, object]:
    query = urllib.parse.urlencode(
        {
            "accession": accession,
            "result": "read_run",
            "fields": "run_accession,fastq_ftp,fastq_md5,fastq_bytes,library_layout",
            "format": "tsv",
        }
    )
    request = urllib.request.Request(f"{ENA_REPORT}?{query}", headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=60) as response:
        rows = list(csv.DictReader(response.read().decode().splitlines(), delimiter="\t"))
    if len(rows) != 1:
        raise RuntimeError(f"ENA returned {len(rows)} rows for {accession}")
    row = rows[0]
    return {
        "run_accession": row["run_accession"],
        "urls": [f"https://{value}" for value in row["fastq_ftp"].split(";")],
        "md5": row["fastq_md5"].split(";"),
        "bytes": [int(value) for value in row["fastq_bytes"].split(";")],
        "library_layout": row["library_layout"],
    }


def verify_archive_metadata(dataset: dict[str, object]) -> dict[str, object]:
    accession = str(dataset["accession"])
    observed = ena_metadata(accession)
    files = list(dataset["files"])
    expected = {
        "run_accession": accession,
        "urls": [str(item["url"]) for item in files],
        "md5": [str(item["source_md5"]) for item in files],
        "bytes": [int(item["source_bytes"]) for item in files],
        "library_layout": str(dataset["library_layout"]),
    }
    if observed != expected:
        raise RuntimeError(
            f"archive metadata changed for {accession}:\n"
            f"expected={json.dumps(expected, sort_keys=True)}\n"
            f"observed={json.dumps(observed, sort_keys=True)}"
        )
    return observed


def stream_subset(url: str, output: Path, reads: int) -> dict[str, object]:
    output.parent.mkdir(parents=True, exist_ok=True)
    partial = output.with_name(output.name + ".partial")
    for attempt in range(1, 4):
        sequence_bases = 0
        names_digest = hashlib.sha256()
        try:
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(request, timeout=60) as response:
                with gzip.GzipFile(fileobj=response) as source, partial.open("wb") as sink:
                    for index in range(reads):
                        header = source.readline()
                        sequence = source.readline()
                        plus = source.readline()
                        quality = source.readline()
                        if not header:
                            raise RuntimeError(
                                f"source ended after {index:,} records; expected {reads:,}"
                            )
                        if not sequence or not plus or not quality:
                            raise RuntimeError(f"truncated FASTQ record {index} from {url}")
                        if not header.startswith(b"@") or not plus.startswith(b"+"):
                            raise RuntimeError(f"malformed FASTQ record {index} from {url}")
                        if len(sequence.rstrip(b"\r\n")) != len(quality.rstrip(b"\r\n")):
                            raise RuntimeError(f"sequence/quality mismatch at record {index}")
                        sink.writelines((header, sequence, plus, quality))
                        sequence_bases += len(sequence.rstrip(b"\r\n"))
                        names_digest.update(normalized_name(header))
                        names_digest.update(b"\0")
                        if (index + 1) % 1_000_000 == 0:
                            print(f"  {output.name}: {index + 1:,}/{reads:,} records", flush=True)
            partial.replace(output)
            return {
                "path": str(output.resolve()),
                "records": reads,
                "sequence_bases": sequence_bases,
                "bytes": output.stat().st_size,
                "sha256": sha256_file(output),
                "normalized_names_sha256": names_digest.hexdigest(),
            }
        except Exception:
            if partial.exists():
                partial.unlink()
            if attempt == 3:
                raise
            delay = 2 ** (attempt - 1)
            print(f"  retrying {output.name} after {delay}s (attempt {attempt + 1}/3)", flush=True)
            time.sleep(delay)
    raise AssertionError("unreachable")


def existing_record(output: Path, expected: dict[str, object] | None) -> dict[str, object] | None:
    if not output.is_file() or expected is None:
        return None
    if int(expected.get("bytes", -1)) != output.stat().st_size:
        return None
    if str(expected.get("sha256")) != sha256_file(output):
        return None
    return expected


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--dest", type=Path, default=ROOT / "data")
    parser.add_argument("--only", action="append", help="dataset name; repeat to select several")
    parser.add_argument(
        "--limit-reads",
        type=int,
        help="development/smoke-test cap; must not exceed the manifest subset size",
    )
    parser.add_argument(
        "--skip-live-metadata-check",
        action="store_true",
        help="use the checked-in archive manifest without querying ENA first",
    )
    args = parser.parse_args()
    manifest_path = args.manifest.resolve()
    manifest = json.loads(manifest_path.read_text())
    selected = set(args.only or [])
    known = {str(item["name"]) for item in manifest["datasets"]}
    unknown = selected - known
    if unknown:
        parser.error(f"unknown datasets: {', '.join(sorted(unknown))}")
    if args.limit_reads is not None and args.limit_reads <= 0:
        parser.error("--limit-reads must be positive")

    dest = args.dest.resolve()
    dest.mkdir(parents=True, exist_ok=True)
    provenance_path = dest / "paper_data_provenance.json"
    previous = json.loads(provenance_path.read_text()) if provenance_path.is_file() else {}
    prior_outputs = {
        item["output"]: item for item in previous.get("outputs", []) if "output" in item
    }
    datasets = [
        item for item in manifest["datasets"] if not selected or str(item["name"]) in selected
    ]
    processed_names = {str(item["name"]) for item in datasets}
    outputs = [
        dict(item)
        for item in previous.get("outputs", [])
        if str(item.get("dataset")) not in processed_names
    ]
    for dataset in datasets:
        name = str(dataset["name"])
        print(f"{name} ({dataset['accession']}):", flush=True)
        if not args.skip_live_metadata_check:
            verify_archive_metadata(dataset)
            print("  archive metadata matches manifest", flush=True)
        lane_records = []
        for source in dataset["files"]:
            subset_reads = int(dataset["subset_reads"])
            if args.limit_reads is not None:
                if args.limit_reads > subset_reads:
                    parser.error(
                        f"--limit-reads exceeds manifest subset size for {name}: {subset_reads}"
                    )
                subset_reads = args.limit_reads
            relative = str(source["output"])
            output = dest / relative
            prior = existing_record(output, prior_outputs.get(relative))
            if prior is not None:
                record = dict(prior)
                print(f"  verified existing {relative}", flush=True)
            elif output.exists():
                raise RuntimeError(
                    f"existing output does not match prior provenance: {output}; "
                    "move it aside before staging"
                )
            else:
                print(f"  streaming {subset_reads:,} reads to {relative}", flush=True)
                record = stream_subset(str(source["url"]), output, subset_reads)
                record["output"] = relative
            record.update(
                {
                    "dataset": name,
                    "accession": dataset["accession"],
                    "lane": source["lane"],
                    "source_url": source["url"],
                    "source_md5": source["source_md5"],
                    "source_bytes": source["source_bytes"],
                }
            )
            lane_records.append(record)
            outputs.append(record)
        if len(lane_records) == 2 and (
            lane_records[0]["normalized_names_sha256"]
            != lane_records[1]["normalized_names_sha256"]
        ):
            raise RuntimeError(f"paired read names differ for {name}")

    provenance = {
        "schema_version": "1.0.0",
        "manifest": {
            "path": str(manifest_path),
            "sha256": sha256_file(manifest_path),
            "schema_version": manifest["schema_version"],
        },
        "subset_definition": "first N complete FASTQ records in each ENA archive file",
        "outputs": outputs,
    }
    provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(f"wrote {provenance_path}")


if __name__ == "__main__":
    main()
