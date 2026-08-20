#!/usr/bin/env python3
"""Record inputs, versions, parameters, and reference metadata for a downstream run."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def described(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def version(command: list[str]) -> str:
    result = subprocess.run(command, check=True, text=True, capture_output=True)
    return (result.stdout or result.stderr).strip()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--genome", type=Path, required=True)
    parser.add_argument("--star", type=Path, required=True)
    parser.add_argument("--python", type=Path, required=True)
    parser.add_argument("--whitelist", type=Path, required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--min-umi", type=int, required=True)
    parser.add_argument("--cb-match-whitelist-type", default="1MM")
    parser.add_argument("--results-root", type=Path, required=True)
    parser.add_argument("--workflow-file", type=Path, action="append", default=[])
    parser.add_argument("--fastq-provenance", type=Path, action="append", required=True)
    args = parser.parse_args()

    genome_manifest = args.genome / "seqproc_reference_manifest.json"
    output_patterns = [
        "resources.csv",
        "*_Solo.out/Gene/raw/barcodes.tsv*",
        "*_Solo.out/Gene/raw/features.tsv*",
        "*_Solo.out/Gene/raw/matrix.mtx*",
        "*_Solo.out/Gene/Summary.csv",
        "*_Solo.out/Barcodes.stats",
        "analysis/*.json",
        "analysis/*.md",
        "analysis/*.pdf",
        "analysis/*.png",
    ]
    output_paths = sorted(
        {path for pattern in output_patterns for path in args.results_root.glob(pattern)}
    )
    result = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "host": {
            "hostname": platform.node(),
            "platform": platform.platform(),
            "cpu_count": os.cpu_count(),
        },
        "parameters": {
            "threads": args.threads,
            "min_umi": args.min_umi,
            "random_seed": 0,
            "solo_type": "CB_UMI_Complex",
            "cb_positions": ["0_10_0_17", "0_18_0_25", "0_26_0_33"],
            "umi_position": "0_0_0_9",
            "cb_match_whitelist_type": args.cb_match_whitelist_type,
            "features": "Gene",
            "cell_filter": "None",
            "out_sam_type": "None",
        },
        "software": {
            "STAR": version([str(args.star), "--version"]),
            "python": version([str(args.python), "--version"]),
            "python_packages": "\n".join(
                sorted(
                    f"{dist.metadata.get('Name', 'unknown')}=={dist.version}"
                    for dist in importlib.metadata.distributions()
                )
            ),
        },
        "whitelist": {
            "path": str(args.whitelist.resolve()),
            "sha256": sha256(args.whitelist),
            "entries": sum(1 for line in args.whitelist.open() if line.strip()),
        },
        "genome_index": {
            "path": str(args.genome.resolve()),
            "manifest": (
                json.loads(genome_manifest.read_text()) if genome_manifest.exists() else None
            ),
        },
        "fastqs": [json.loads(path.read_text()) for path in args.fastq_provenance],
        "workflow_files": [described(path) for path in args.workflow_file],
        "outputs": [described(path) for path in output_paths],
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
