#!/usr/bin/env python3
"""Benchmark seqproc compressed I/O and persistent output on paper FASTQs.

The harness creates deterministic FASTQ prefixes, stages matched plain and
gzip inputs, randomizes external replicates across output filesystems, and
requires decompressed output identity across every mode. Parallel gzip output
is additionally decoded and integrity-checked by Python, gzip, and pigz. The
single-stream mode sweeps compression-worker and deflate-block settings
independently of seqproc's read-batch size.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import random
import statistics
import subprocess
import time
import zlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONFIGS = ROOT / "configs" / "seqproc"

DATASETS = (
    {
        "name": "splitseq_pe",
        "geometry": "publication_splitseq_pe.geom",
        "r1": "SRR6750041_10M_R1.fastq",
        "r2": "SRR6750041_10M_R2.fastq",
        "additional": (),
    },
    {
        "name": "lr_splitseq",
        "geometry": "publication_lr_splitseq_dual_core.geom",
        "r1": "SRR13948564_1M.fastq",
        "r2": None,
        "additional": (),
    },
    {
        "name": "10x_short",
        "geometry": "10x_v2.geom",
        "r1": "10x_short/SRR8315379_1M_R1.fastq",
        "r2": "10x_short/SRR8315379_1M_R2.fastq",
        "additional": (),
    },
    {
        "name": "sciseq",
        "geometry": "sciseq3_edit.geom",
        "r1": "SRR7827254_1M_1.fastq",
        "r2": "SRR7827254_1M_2.fastq",
        "additional": (),
    },
)

MODES = (
    {
        "name": "plain_input_plain_output",
        "gzip_input": False,
        "gzip_output": False,
        "parallel_gzip": False,
    },
    {
        "name": "gzip_input_plain_output",
        "gzip_input": True,
        "gzip_output": False,
        "parallel_gzip": False,
    },
    {
        "name": "gzip_input_plain_output_accelerated",
        "gzip_input": True,
        "gzip_output": False,
        "parallel_gzip": False,
        "accelerated_gzip_input": True,
    },
    {
        "name": "plain_input_gzip_serial",
        "gzip_input": False,
        "gzip_output": True,
        "parallel_gzip": False,
    },
    {
        "name": "plain_input_gzip_parallel",
        "gzip_input": False,
        "gzip_output": True,
        "parallel_gzip": True,
    },
    {
        "name": "gzip_input_gzip_parallel",
        "gzip_input": True,
        "gzip_output": True,
        "parallel_gzip": True,
    },
    {
        "name": "plain_input_gzip_stream",
        "gzip_input": False,
        "gzip_output": True,
        "parallel_gzip": False,
        "parallel_gzip_stream": True,
    },
    {
        "name": "gzip_input_gzip_stream",
        "gzip_input": True,
        "gzip_output": True,
        "parallel_gzip": False,
        "parallel_gzip_stream": True,
    },
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def copy_fastq_prefix(source: Path, destination: Path, records: int) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".partial")
    with source.open("rb") as reader, temporary.open("wb") as writer:
        for index in range(records):
            record = [reader.readline() for _ in range(4)]
            if any(not line for line in record):
                raise RuntimeError(
                    f"{source} ended after {index} records; requested {records}"
                )
            if not record[0].startswith(b"@") or not record[2].startswith(b"+"):
                raise RuntimeError(f"malformed FASTQ record {index + 1} in {source}")
            if len(record[1].rstrip(b"\r\n")) != len(record[3].rstrip(b"\r\n")):
                raise RuntimeError(
                    f"sequence/quality length mismatch at record {index + 1} in {source}"
                )
            writer.writelines(record)
    temporary.replace(destination)


def stage_gzip(source: Path, destination: Path) -> None:
    temporary = destination.with_suffix(destination.suffix + ".partial")
    with temporary.open("wb") as writer:
        completed = subprocess.run(
            ("pigz", "-3", "-c", str(source)), stdout=writer, stderr=subprocess.PIPE
        )
    if completed.returncode != 0:
        raise RuntimeError(
            f"pigz failed for {source}: {completed.stderr.decode(errors='replace')}"
        )
    temporary.replace(destination)


def count_fastq_records(path: Path, compressed: bool) -> int:
    opener = gzip.open if compressed else Path.open
    lines = 0
    if compressed:
        handle_context = opener(path, "rb")
    else:
        handle_context = opener(path, "rb")
    with handle_context as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            lines += chunk.count(b"\n")
    if lines % 4:
        raise RuntimeError(f"output is not complete four-line FASTQ: {path}")
    return lines // 4


def gzip_member_count(path: Path) -> int:
    remaining = path.read_bytes()
    members = 0
    while remaining:
        decoder = zlib.decompressobj(16 + zlib.MAX_WBITS)
        decoder.decompress(remaining)
        decoder.flush()
        if not decoder.eof:
            raise RuntimeError(f"incomplete gzip member in {path}")
        members += 1
        remaining = decoder.unused_data
    return members


def decoded_sha256(path: Path, compressed: bool) -> str:
    digest = hashlib.sha256()
    opener = gzip.open if compressed else Path.open
    with opener(path, "rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def parse_time_report(path: Path) -> dict[str, float]:
    values: dict[str, float] = {}
    mapping = {
        "User time (seconds)": "user_seconds",
        "System time (seconds)": "system_seconds",
        "Maximum resident set size (kbytes)": "peak_rss_kib",
        "File system inputs": "filesystem_inputs",
        "File system outputs": "filesystem_outputs",
    }
    for line in path.read_text().splitlines():
        stripped = line.strip()
        for prefix, key in mapping.items():
            if stripped.startswith(prefix + ":"):
                values[key] = float(stripped.rsplit(":", 1)[1].strip())
    return values


def filesystem_metadata(path: Path) -> dict[str, object]:
    completed = subprocess.run(
        ("findmnt", "-T", str(path), "-n", "-o", "SOURCE,FSTYPE,OPTIONS"),
        check=True,
        text=True,
        capture_output=True,
    )
    mounts = []
    for line in completed.stdout.splitlines():
        fields = line.strip().split(maxsplit=2)
        if fields:
            mounts.append(
                {
                    "source": fields[0],
                    "filesystem_type": fields[1] if len(fields) > 1 else None,
                    "mount_options": fields[2] if len(fields) > 2 else None,
                }
            )
    stats = os.statvfs(path)
    return {
        "path": str(path),
        "mounts": mounts,
        "block_size": stats.f_frsize,
    }


def parse_output_root(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("output roots must use LABEL=PATH")
    label, raw_path = value.split("=", 1)
    if not label or not raw_path:
        raise argparse.ArgumentTypeError("output roots must use nonempty LABEL=PATH")
    return label, Path(raw_path)


def command_for(
    binary: Path,
    dataset: dict[str, object],
    staged: dict[str, Path | None],
    mode: dict[str, object],
    output1: Path,
    output2: Path | None,
    threads: int,
    batch_size: int | None,
    gzip_threads: int | None,
    gzip_block_size: int | None,
    accelerated_gzip_input_threads: int | None,
    accelerated_gzip_input_chunk_size: int | None,
    affinity: str | None,
    time_report: Path,
) -> list[str]:
    input_suffix = "_gzip" if bool(mode["gzip_input"]) else "_plain"
    input1 = staged[f"r1{input_suffix}"]
    input2 = staged[f"r2{input_suffix}"]
    assert input1 is not None
    command = ["/usr/bin/time", "-v", "-o", str(time_report)]
    if affinity:
        command.extend(("taskset", "-c", affinity))
    command.extend(
        (
            str(binary),
            "run",
            "--geom",
            str(CONFIGS / str(dataset["geometry"])),
            "--file1",
            str(input1),
            "--out1",
            str(output1),
            "--threads",
            str(threads),
            "--preserve-order",
        )
    )
    if input2 is not None and output2 is not None:
        command.extend(("--file2", str(input2), "--out2", str(output2)))
    if bool(mode["gzip_output"]):
        command.extend(("--gzip-level", "3"))
    if bool(mode["parallel_gzip"]):
        command.append("--parallel-gzip")
    if bool(mode.get("parallel_gzip_stream", False)):
        assert gzip_threads is not None and gzip_block_size is not None
        command.extend(
            (
                "--parallel-gzip-stream",
                "--gzip-threads",
                str(gzip_threads),
                "--gzip-block-size",
                str(gzip_block_size),
            )
        )
    if bool(mode.get("accelerated_gzip_input", False)):
        assert accelerated_gzip_input_threads is not None
        assert accelerated_gzip_input_chunk_size is not None
        command.extend(
            (
                "--accelerated-gzip-input",
                "--gzip-input-threads",
                str(accelerated_gzip_input_threads),
                "--gzip-input-chunk-size",
                str(accelerated_gzip_input_chunk_size),
            )
        )
    if batch_size is not None:
        command.extend(("--batch-size", str(batch_size)))
    for additional in dataset["additional"]:
        command.extend(("--additional", str(CONFIGS / str(additional))))
    return command


def summarize(runs: list[dict[str, object]]) -> list[dict[str, object]]:
    groups: dict[
        tuple[
            str,
            str,
            str,
            int | None,
            int | None,
            int | None,
            int | None,
            int | None,
        ],
        list[dict[str, object]],
    ] = {}
    for run in runs:
        key = (
            str(run["dataset"]),
            str(run["mode"]),
            str(run["output_root"]),
            int(run["batch_size"]) if run["batch_size"] is not None else None,
            int(run["gzip_threads"]) if run["gzip_threads"] is not None else None,
            int(run["gzip_block_size"])
            if run["gzip_block_size"] is not None
            else None,
            int(run["accelerated_gzip_input_threads"])
            if run["accelerated_gzip_input_threads"] is not None
            else None,
            int(run["accelerated_gzip_input_chunk_size"])
            if run["accelerated_gzip_input_chunk_size"] is not None
            else None,
        )
        groups.setdefault(key, []).append(run)
    summary = []
    for (
        dataset,
        mode,
        root,
        batch_size,
        gzip_threads,
        gzip_block_size,
        accelerated_gzip_input_threads,
        accelerated_gzip_input_chunk_size,
    ), values in sorted(groups.items(), key=lambda item: tuple(str(value) for value in item[0])):
        summary.append(
            {
                "dataset": dataset,
                "mode": mode,
                "output_root": root,
                "batch_size": batch_size,
                "gzip_threads": gzip_threads,
                "gzip_block_size": gzip_block_size,
                "accelerated_gzip_input_threads": accelerated_gzip_input_threads,
                "accelerated_gzip_input_chunk_size": accelerated_gzip_input_chunk_size,
                "replicates": len(values),
                "median_input_fragments_per_second": statistics.median(
                    float(value["input_fragments_per_second"]) for value in values
                ),
                "median_wall_seconds": statistics.median(
                    float(value["wall_seconds"]) for value in values
                ),
                "median_peak_rss_kib": statistics.median(
                    float(value["peak_rss_kib"]) for value in values
                ),
                "median_output_bytes": statistics.median(
                    int(value["output_bytes"]) for value in values
                ),
                "gzip_members": sorted(
                    {
                        int(member)
                        for value in values
                        for member in value["gzip_members"]
                    }
                ),
            }
        )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--input-workdir", type=Path, required=True)
    parser.add_argument(
        "--output-root",
        action="append",
        type=parse_output_root,
        required=True,
        help="repeatable LABEL=PATH persistent output target",
    )
    parser.add_argument("--reads", type=int, default=250_000)
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument(
        "--batch-size",
        type=int,
        nargs="+",
        help="one or more explicit staged batch sizes; omit to use the runtime default",
    )
    parser.add_argument(
        "--gzip-threads",
        type=int,
        nargs="+",
        help="single-stream compression-worker counts (default: --threads)",
    )
    parser.add_argument(
        "--gzip-block-size",
        type=int,
        nargs="+",
        help="single-stream uncompressed deflate block sizes (default: 131072)",
    )
    parser.add_argument(
        "--accelerated-gzip-input-threads",
        type=int,
        nargs="+",
        help="adaptive rapidgzip-core worker ceilings per input (default: 1)",
    )
    parser.add_argument(
        "--accelerated-gzip-input-chunk-size",
        type=int,
        nargs="+",
        help="accelerated input decoded-byte chunks (default: 262144)",
    )
    parser.add_argument(
        "--mode",
        action="append",
        help="mode name; repeat to select several (default: all modes)",
    )
    parser.add_argument("--affinity", default="0-7")
    parser.add_argument("--only", action="append")
    args = parser.parse_args()
    if (
        args.reads <= 0
        or args.replicates <= 0
        or args.threads <= 0
        or any(value <= 0 for value in (args.batch_size or []))
        or any(value <= 0 for value in (args.gzip_threads or []))
        or any(value < 32768 for value in (args.gzip_block_size or []))
        or any(value <= 0 for value in (args.accelerated_gzip_input_threads or []))
        or any(
            value <= 0 for value in (args.accelerated_gzip_input_chunk_size or [])
        )
    ):
        parser.error(
            "reads, replicates, threads, gzip threads, and batch sizes must be positive; "
            "gzip block sizes must be at least 32768"
        )

    binary = args.binary.resolve()
    data_dir = args.data_dir.resolve()
    outdir = args.outdir.resolve()
    input_workdir = args.input_workdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    input_workdir.mkdir(parents=True, exist_ok=True)
    output_roots = []
    labels = set()
    for label, raw_path in args.output_root:
        if label in labels:
            parser.error(f"duplicate output-root label: {label}")
        labels.add(label)
        path = raw_path.resolve()
        path.mkdir(parents=True, exist_ok=True)
        output_roots.append((label, path))

    selected = set(args.only or [])
    known = {str(dataset["name"]) for dataset in DATASETS}
    if selected - known:
        parser.error(f"unknown datasets: {', '.join(sorted(selected - known))}")
    datasets = [
        dataset for dataset in DATASETS if not selected or str(dataset["name"]) in selected
    ]
    selected_modes = set(args.mode or [])
    known_modes = {str(mode["name"]) for mode in MODES}
    if selected_modes - known_modes:
        parser.error(f"unknown modes: {', '.join(sorted(selected_modes - known_modes))}")
    modes = [
        mode for mode in MODES if not selected_modes or str(mode["name"]) in selected_modes
    ]
    batch_sizes: list[int | None] = sorted(set(args.batch_size or [])) or [None]
    gzip_threads = sorted(set(args.gzip_threads or [args.threads]))
    gzip_block_sizes = sorted(set(args.gzip_block_size or [128 * 1024]))
    accelerated_input_chunk_sizes = sorted(
        set(args.accelerated_gzip_input_chunk_size or [256 * 1024])
    )
    accelerated_input_threads = sorted(
        set(args.accelerated_gzip_input_threads or [1])
    )

    staged_inputs: dict[str, dict[str, Path | None]] = {}
    input_manifest: dict[str, object] = {}
    for dataset in datasets:
        name = str(dataset["name"])
        print(f"[stage] {name}: {args.reads} fragments", flush=True)
        staged: dict[str, Path | None] = {}
        records = []
        for lane in ("r1", "r2"):
            relative = dataset[lane]
            if relative is None:
                staged[f"{lane}_plain"] = None
                staged[f"{lane}_gzip"] = None
                continue
            source = data_dir / str(relative)
            if not source.is_file():
                parser.error(f"missing input: {source}")
            plain = input_workdir / f"{name}-{args.reads}-{lane}.fastq"
            compressed = plain.with_suffix(".fastq.gz")
            copy_fastq_prefix(source, plain, args.reads)
            stage_gzip(plain, compressed)
            if count_fastq_records(compressed, True) != args.reads:
                raise RuntimeError(f"compressed staging count mismatch: {compressed}")
            staged[f"{lane}_plain"] = plain
            staged[f"{lane}_gzip"] = compressed
            records.append(
                {
                    "lane": lane,
                    "source": str(source),
                    "source_sha256": sha256_file(source),
                    "plain": {
                        "path": str(plain),
                        "bytes": plain.stat().st_size,
                        "sha256": sha256_file(plain),
                    },
                    "gzip": {
                        "path": str(compressed),
                        "bytes": compressed.stat().st_size,
                        "sha256": sha256_file(compressed),
                    },
                }
            )
        staged_inputs[name] = staged
        input_manifest[name] = {
            "records": records,
            "geometry": str((CONFIGS / str(dataset["geometry"])).resolve()),
            "geometry_sha256": sha256_file(CONFIGS / str(dataset["geometry"])),
        }

    schedule = [
        (
            dataset,
            mode,
            root_label,
            root_path,
            batch_size,
            stream_threads,
            stream_block_size,
            accelerated_threads,
            accelerated_input_chunk_size,
            replicate,
        )
        for replicate in range(1, args.replicates + 1)
        for dataset in datasets
        for mode in modes
        for root_label, root_path in output_roots
        for batch_size in batch_sizes
        for stream_threads, stream_block_size in (
            [
                (candidate_threads, candidate_block_size)
                for candidate_threads in gzip_threads
                for candidate_block_size in gzip_block_sizes
            ]
            if bool(mode.get("parallel_gzip_stream", False))
            else [(None, None)]
        )
        for accelerated_input_chunk_size in (
            accelerated_input_chunk_sizes
            if bool(mode.get("accelerated_gzip_input", False))
            else [None]
        )
        for accelerated_threads in (
            accelerated_input_threads
            if bool(mode.get("accelerated_gzip_input", False))
            else [None]
        )
    ]
    random.Random("seqproc-real-io-v1").shuffle(schedule)
    reference_digests: dict[str, str] = {}
    runs = []
    for index, (
        dataset,
        mode,
        root_label,
        root_path,
        batch_size,
        stream_threads,
        stream_block_size,
        accelerated_threads,
        accelerated_input_chunk_size,
        replicate,
    ) in enumerate(schedule, 1):
        name = str(dataset["name"])
        compressed = bool(mode["gzip_output"])
        suffix = ".fastq.gz" if compressed else ".fastq"
        batch_label = f"b{batch_size}" if batch_size is not None else "bdefault"
        stream_label = (
            f"-g{stream_threads}-z{stream_block_size}"
            if stream_threads is not None
            else ""
        )
        input_label = (
            f"-d{accelerated_threads}-i{accelerated_input_chunk_size}"
            if accelerated_input_chunk_size is not None
            else ""
        )
        label = (
            f"{name}-{mode['name']}-{root_label}-{batch_label}{stream_label}"
            f"{input_label}-r{replicate}"
        )
        output1 = root_path / f"{label}-R1{suffix}"
        output2 = root_path / f"{label}-R2{suffix}" if dataset["r2"] else None
        time_report = outdir / f"{label}.time"
        command = command_for(
            binary,
            dataset,
            staged_inputs[name],
            mode,
            output1,
            output2,
            args.threads,
            batch_size,
            stream_threads,
            stream_block_size,
            accelerated_threads,
            accelerated_input_chunk_size,
            args.affinity or None,
            time_report,
        )
        print(
            f"[{index}/{len(schedule)}] {name} {mode['name']} {root_label} "
            f"batch={batch_size or 'default'} gzip_threads={stream_threads or '-'} "
            f"gzip_block={stream_block_size or '-'} "
            f"gzip_input_threads={accelerated_threads or '-'} "
            f"gzip_input_chunk={accelerated_input_chunk_size or '-'} replicate={replicate}",
            flush=True,
        )
        started = time.perf_counter()
        completed = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
        elapsed = time.perf_counter() - started
        if completed.returncode != 0:
            raise RuntimeError(
                f"{label} failed ({completed.returncode})\nstdout:\n{completed.stdout}\n"
                f"stderr:\n{completed.stderr}"
            )

        output_paths = [output1] + ([output2] if output2 is not None else [])
        lane_digests = []
        output_records = []
        members = []
        for output in output_paths:
            if compressed:
                subprocess.run(("gzip", "-t", str(output)), check=True)
                subprocess.run(("pigz", "-t", str(output)), check=True)
                members.append(gzip_member_count(output))
            lane_digests.append(decoded_sha256(output, compressed))
            output_records.append(count_fastq_records(output, compressed))
        composite = hashlib.sha256("\0".join(lane_digests).encode()).hexdigest()
        previous = reference_digests.setdefault(name, composite)
        if composite != previous:
            raise RuntimeError(f"decompressed output divergence for {name}: {composite} != {previous}")
        if len(set(output_records)) != 1:
            raise RuntimeError(f"paired output record counts diverged for {label}: {output_records}")

        runs.append(
            {
                "dataset": name,
                "mode": mode["name"],
                "output_root": root_label,
                "replicate": replicate,
                "threads": args.threads,
                "batch_size": batch_size,
                "gzip_threads": stream_threads,
                "gzip_block_size": stream_block_size,
                "accelerated_gzip_input_threads": accelerated_threads,
                "accelerated_gzip_input_chunk_size": accelerated_input_chunk_size,
                "accelerated_gzip_input": bool(
                    mode.get("accelerated_gzip_input", False)
                ),
                "command": command,
                "input_fragments": args.reads,
                "output_fragments": output_records[0],
                "wall_seconds": elapsed,
                "input_fragments_per_second": args.reads / elapsed,
                "output_bytes": sum(path.stat().st_size for path in output_paths),
                "decoded_lane_sha256": lane_digests,
                "decoded_composite_sha256": composite,
                "gzip_members": members,
                **parse_time_report(time_report),
            }
        )
        time_report.unlink()
        for output in output_paths:
            output.unlink()

    artifact = {
        "schema_version": "1.0.0",
        "scope": "real-paper FASTQ persistent and compressed I/O development matrix",
        "binary": {"path": str(binary), "sha256": sha256_file(binary)},
        "host": {
            "hostname": os.uname().nodename,
            "machine": os.uname().machine,
            "cpu_affinity": args.affinity or None,
        },
        "parameters": {
            "input_fragments_per_dataset": args.reads,
            "replicates": args.replicates,
            "threads": args.threads,
            "batch_sizes": batch_sizes,
            "gzip_threads": gzip_threads,
            "gzip_block_sizes": gzip_block_sizes,
            "accelerated_gzip_input_threads": accelerated_input_threads,
            "accelerated_gzip_input_chunk_sizes": accelerated_input_chunk_sizes,
            "gzip_level": 3,
            "preserve_order": True,
            "randomized_schedule": True,
            "statistics": False,
        },
        "filesystems": {
            label: filesystem_metadata(path) for label, path in output_roots
        },
        "inputs": input_manifest,
        "reference_decoded_composite_sha256": reference_digests,
        "runs": runs,
        "summary": summarize(runs),
    }
    artifact_path = outdir / "real-seqproc-io-results.json"
    artifact_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {artifact_path}")


if __name__ == "__main__":
    main()
