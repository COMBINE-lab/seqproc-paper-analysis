#!/usr/bin/env python3
"""Benchmark seqproc execution modes on deterministic manuscript geometries.

The generated reads are protocol-shaped synthetic inputs, not substitutes for
the biological datasets used in the paper.  They exercise exact and fuzzy
anchors, whitelist correction, paired/single-end input, and reverse orientation.
Every run is checked against an order-independent digest of the emitted FASTQ.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import itertools
import json
import random
import shutil
import statistics
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONFIGS = ROOT / "configs" / "seqproc"
SCHEMA_VERSION = "1.0.0"
DNA = "ACGT"


@dataclass(frozen=True)
class Protocol:
    name: str
    geometry: str
    paired: bool
    additional: tuple[str, ...] = ()


PROTOCOLS = (
    Protocol("10x_short", "10x_v2.geom", True),
    Protocol("sciseq", "sciseq3_edit.geom", True),
    Protocol(
        "splitseq_pe",
        "splitseq_filter_edit.geom",
        True,
        (
            "splitseq_bc3_seq2seq.tsv",
            "splitseq_bc2_seq2seq.tsv",
            "splitseq_bc1_seq2seq.tsv",
        ),
    ),
    Protocol(
        "splitseq_se",
        "splitseq_singleend_edit_ann.geom",
        False,
        (
            "splitseq_bc3_seq2seq.tsv",
            "splitseq_bc2_seq2seq.tsv",
            "splitseq_bc1_seq2seq.tsv",
        ),
    ),
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def sequence(rng: random.Random, length: int) -> str:
    return "".join(rng.choice(DNA) for _ in range(length))


def substitute(value: str, position: int = 0) -> str:
    current = value[position]
    replacement = DNA[(DNA.index(current) + 1) % len(DNA)]
    return value[:position] + replacement + value[position + 1 :]


def fuzz_anchor(value: str, case: int) -> str:
    """Cycle through exact, substitution, insertion, and deletion cases."""
    mode = case % 4
    if mode == 0:
        return value
    if mode == 1:
        return substitute(value, len(value) // 2)
    if mode == 2:
        return value[: len(value) // 2] + "A" + value[len(value) // 2 :]
    return value[: len(value) // 2] + value[len(value) // 2 + 1 :]


def reverse_complement(value: str) -> str:
    return value.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def load_lines(path: Path) -> list[str]:
    return [line.strip().split()[0] for line in path.read_text().splitlines() if line.strip()]


def write_record(handle, name: str, value: str) -> None:
    handle.write(f"@{name}\n{value}\n+\n{'I' * len(value)}\n")


def generate_protocol(protocol: Protocol, reads: int, data_dir: Path) -> tuple[Path, Path | None]:
    rng = random.Random(f"seqproc-integration-v1:{protocol.name}:{reads}")
    r1_path = data_dir / f"{protocol.name}_R1.fastq"
    r2_path = data_dir / f"{protocol.name}_R2.fastq" if protocol.paired else None
    bc23 = load_lines(CONFIGS / "splitseq_bc23_whitelist.txt")
    bc1 = load_lines(CONFIGS / "splitseq_bc1_whitelist_6bp.txt")
    linker1_pe = "GTGGCCGCTGTTTCGCATCGGCGTACGACT"
    linker1_se = "GTGGCCGATGTTTCGCATCGGCGTACGACT"
    linker2_pe = "ATCCACGTGCTTGAGAGGCCAGAGCATTCG"
    linker2_se = "ATCCACGTGCTTGAGACTGTGG"

    with r1_path.open("w") as r1_handle:
        r2_context = r2_path.open("w") if r2_path else None
        try:
            for index in range(reads):
                name = f"{protocol.name}:{index}"
                if protocol.name == "10x_short":
                    r1 = sequence(rng, 16) + sequence(rng, 10)
                    r2 = sequence(rng, 100)
                elif protocol.name == "sciseq":
                    r1 = (
                        sequence(rng, 9)
                        + (fuzz_anchor("CAGAGC", index) if index % 3 else "CAGAGC")
                        + sequence(rng, 8)
                        + sequence(rng, 10)
                    )
                    r2 = sequence(rng, 100)
                elif protocol.name == "splitseq_pe":
                    chosen3 = bc23[index % len(bc23)]
                    chosen2 = bc23[(index * 7 + 3) % len(bc23)]
                    chosen1 = bc1[(index * 11 + 5) % len(bc1)]
                    if index % 5 == 1:
                        chosen3 = substitute(chosen3, 2)
                    if index % 5 == 2:
                        chosen2 = substitute(chosen2, 5)
                    if index % 5 == 3:
                        chosen1 = substitute(chosen1, 1)
                    r1 = sequence(rng, 100)
                    r2 = (
                        sequence(rng, 2)
                        + sequence(rng, 8)
                        + chosen3
                        + fuzz_anchor(linker1_pe, index)
                        + chosen2
                        + fuzz_anchor(linker2_pe, index + 1)
                        + chosen1
                        + sequence(rng, 20)
                    )
                else:
                    chosen3 = bc23[index % len(bc23)]
                    chosen2 = bc23[(index * 7 + 3) % len(bc23)]
                    chosen1 = bc1[(index * 11 + 5) % len(bc1)] + sequence(rng, 2)
                    r1 = (
                        sequence(rng, 10)
                        + chosen3
                        + fuzz_anchor(linker1_se, index)
                        + chosen2
                        + fuzz_anchor(linker2_se, index + 1)
                        + chosen1
                        + sequence(rng, 50)
                    )
                    if index % 5 == 4:
                        r1 = reverse_complement(r1)
                    r2 = ""
                write_record(r1_handle, name, r1)
                if r2_context is not None:
                    write_record(r2_context, name, r2)
        finally:
            if r2_context is not None:
                r2_context.close()
    return r1_path, r2_path


def open_fastq(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def records(path: Path):
    with open_fastq(path) as handle:
        while True:
            name = handle.readline()
            if not name:
                return
            sequence_line = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not sequence_line or not plus or not quality:
                raise RuntimeError(f"truncated FASTQ: {path}")
            yield name.rstrip(), sequence_line.rstrip(), plus.rstrip(), quality.rstrip()


def output_digest(r1_path: Path, r2_path: Path | None) -> tuple[str, int]:
    hashes: list[bytes] = []
    if r2_path is None:
        iterator = ((left, None) for left in records(r1_path))
    else:
        iterator = itertools.zip_longest(records(r1_path), records(r2_path))
    for left, right in iterator:
        if left is None or (r2_path is not None and right is None):
            raise RuntimeError("paired output FASTQs contain different record counts")
        if right is not None and left[0].split()[0] != right[0].split()[0]:
            raise RuntimeError(f"paired output names diverged: {left[0]} != {right[0]}")
        payload = "\n".join(left).encode()
        if right is not None:
            payload += b"\0" + "\n".join(right).encode()
        hashes.append(hashlib.sha256(payload).digest())
    hashes.sort()
    digest = hashlib.sha256(b"".join(hashes)).hexdigest()
    return digest, len(hashes)


def parse_time_report(path: Path) -> dict[str, float]:
    values: dict[str, float] = {}
    mapping = {
        "User time (seconds)": "user_seconds",
        "System time (seconds)": "system_seconds",
        "Maximum resident set size (kbytes)": "peak_rss_kib",
    }
    for line in path.read_text().splitlines():
        stripped = line.strip()
        for prefix, key in mapping.items():
            if stripped.startswith(prefix + ":"):
                values[key] = float(stripped.rsplit(":", 1)[1].strip())
    return values


def modes(thread_counts: list[int]) -> list[dict[str, object]]:
    result: list[dict[str, object]] = [
        {"name": "legacy_plain", "threads": 1, "flags": []}
    ]
    result.extend(
        {
            "name": "staged_plain",
            "threads": threads,
            "flags": ["--staged-pipeline"],
        }
        for threads in thread_counts
    )
    top = max(thread_counts)
    result.extend(
        (
            {
                "name": "staged_plain_stats",
                "threads": top,
                "flags": ["--staged-pipeline"],
                "summary": True,
            },
            {"name": "ordered_plain", "threads": top, "flags": ["--preserve-order"]},
            {
                "name": "ordered_gzip_serial",
                "threads": top,
                "flags": ["--preserve-order", "--gzip-level", "3"],
                "gzip": True,
            },
            {
                "name": "ordered_gzip_parallel",
                "threads": top,
                "flags": ["--preserve-order", "--gzip-level", "3", "--parallel-gzip"],
                "gzip": True,
            },
        )
    )
    return result


def run_once(
    binary: Path,
    protocol: Protocol,
    r1: Path,
    r2: Path | None,
    mode: dict[str, object],
    replicate: int,
    work_dir: Path,
) -> dict[str, object]:
    label = f"{protocol.name}-{mode['name']}-t{mode['threads']}-r{replicate}"
    gzip_output = bool(mode.get("gzip"))
    suffix = ".fastq.gz" if gzip_output else ".fastq"
    out1 = work_dir / f"{label}_R1{suffix}"
    out2 = work_dir / f"{label}_R2{suffix}" if protocol.paired else None
    summary = work_dir / f"{label}.summary.json"
    time_report = work_dir / f"{label}.time.txt"
    command = [
        "/usr/bin/time",
        "-v",
        "-o",
        str(time_report),
        str(binary),
        "run",
        "--geom",
        str(CONFIGS / protocol.geometry),
        "--file1",
        str(r1),
        "--out1",
        str(out1),
        "--threads",
        str(mode["threads"]),
        *[str(value) for value in mode["flags"]],
    ]
    if bool(mode.get("summary")):
        command.extend(("--summary", str(summary)))
    if r2 is not None and out2 is not None:
        command.extend(("--file2", str(r2), "--out2", str(out2)))
    for additional in protocol.additional:
        command.extend(("--additional", str(CONFIGS / additional)))

    started = time.perf_counter()
    completed = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    elapsed = time.perf_counter() - started
    if completed.returncode != 0:
        raise RuntimeError(
            f"{label} failed ({completed.returncode})\nstdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )
    digest, output_reads = output_digest(out1, out2)
    report = json.loads(summary.read_text()) if summary.is_file() else {}
    result: dict[str, object] = {
        "protocol": protocol.name,
        "geometry": protocol.geometry,
        "mode": mode["name"],
        "threads": mode["threads"],
        "collect_statistics": bool(mode.get("summary")),
        "replicate": replicate,
        "wall_seconds": elapsed,
        "reads_per_second": output_reads / elapsed,
        "output_reads": output_reads,
        "output_bytes": out1.stat().st_size + (out2.stat().st_size if out2 else 0),
        "output_multiset_sha256": digest,
        "summary_schema": report.get("schema_version"),
        "n_processed": report.get("n_processed"),
        "effective_threads": report.get("effective_threads", mode["threads"]),
        **parse_time_report(time_report),
    }
    for path in (out1, out2, summary, time_report):
        if path is not None and path.exists():
            path.unlink()
    return result


def summarize(runs: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, int], list[dict[str, object]]] = {}
    for run in runs:
        key = (str(run["protocol"]), str(run["mode"]), int(run["threads"]))
        grouped.setdefault(key, []).append(run)
    result = []
    for (protocol, mode, threads), values in sorted(grouped.items()):
        throughputs = [float(value["reads_per_second"]) for value in values]
        result.append(
            {
                "protocol": protocol,
                "mode": mode,
                "threads": threads,
                "replicates": len(values),
                "median_reads_per_second": statistics.median(throughputs),
                "min_reads_per_second": min(throughputs),
                "max_reads_per_second": max(throughputs),
                "median_peak_rss_kib": statistics.median(
                    float(value["peak_rss_kib"]) for value in values
                ),
                "median_output_bytes": statistics.median(
                    int(value["output_bytes"]) for value in values
                ),
            }
        )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--reads", type=int, default=100_000)
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 2, 4, 8])
    args = parser.parse_args()
    if args.reads <= 0 or args.replicates <= 0 or any(value <= 0 for value in args.threads):
        parser.error("reads, replicates, and thread counts must be positive")
    binary = args.binary.resolve()
    if not binary.is_file():
        parser.error(f"binary does not exist: {binary}")

    outdir = args.outdir.resolve()
    data_dir = outdir / "inputs"
    work_dir = outdir / "work"
    data_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    inputs = {}
    generated = {}
    for protocol in PROTOCOLS:
        r1, r2 = generate_protocol(protocol, args.reads, data_dir)
        generated[protocol.name] = (r1, r2)
        inputs[protocol.name] = {
            "r1": {"path": str(r1), "bytes": r1.stat().st_size, "sha256": sha256_file(r1)},
            "r2": (
                {"path": str(r2), "bytes": r2.stat().st_size, "sha256": sha256_file(r2)}
                if r2
                else None
            ),
            "geometry_sha256": sha256_file(CONFIGS / protocol.geometry),
        }

    schedule = [
        (protocol, mode, replicate)
        for replicate in range(1, args.replicates + 1)
        for protocol in PROTOCOLS
        for mode in modes(sorted(set(args.threads)))
    ]
    random.Random("seqproc-integration-schedule-v1").shuffle(schedule)
    runs = []
    reference_digests: dict[str, str] = {}
    for index, (protocol, mode, replicate) in enumerate(schedule, 1):
        print(
            f"[{index}/{len(schedule)}] {protocol.name} {mode['name']} "
            f"t={mode['threads']} replicate={replicate}",
            flush=True,
        )
        r1, r2 = generated[protocol.name]
        run = run_once(binary, protocol, r1, r2, mode, replicate, work_dir)
        previous = reference_digests.setdefault(
            protocol.name, str(run["output_multiset_sha256"])
        )
        if run["output_multiset_sha256"] != previous:
            raise RuntimeError(
                f"output divergence for {protocol.name}: {run['output_multiset_sha256']} != {previous}"
            )
        if run["output_reads"] != args.reads or (
            run["n_processed"] is not None and run["n_processed"] != args.reads
        ):
            raise RuntimeError(
                f"unexpected accepted count for {protocol.name}: output={run['output_reads']}, "
                f"summary={run['n_processed']}, expected={args.reads}"
            )
        runs.append(run)

    artifact = {
        "schema_version": SCHEMA_VERSION,
        "scope": "protocol-shaped synthetic integration benchmark; not biological data",
        "binary": {"path": str(binary), "sha256": sha256_file(binary)},
        "parameters": {
            "reads_per_protocol": args.reads,
            "replicates": args.replicates,
            "thread_counts": sorted(set(args.threads)),
            "randomized_schedule": True,
        },
        "inputs": inputs,
        "reference_output_multiset_sha256": reference_digests,
        "runs": runs,
        "summary": summarize(runs),
    }
    artifact_path = outdir / "protocol-integration-results.json"
    artifact_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    shutil.rmtree(work_dir)
    print(f"wrote {artifact_path}")


if __name__ == "__main__":
    main()
