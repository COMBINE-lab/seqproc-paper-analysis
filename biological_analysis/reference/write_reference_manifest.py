#!/usr/bin/env python3
"""Write checksummed provenance for the Ensembl/STAR reference index."""

from __future__ import annotations

import argparse
import hashlib
import json
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
    return {"path": str(path.resolve()), "bytes": path.stat().st_size, "sha256": sha256(path)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--star", type=Path, required=True)
    parser.add_argument("--index", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--fasta-url", required=True)
    parser.add_argument("--fasta-gzip-sha256", required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--gtf-url", required=True)
    parser.add_argument("--gtf-gzip-sha256", required=True)
    parser.add_argument("--star-source-url", required=True)
    parser.add_argument("--star-source-sha256", required=True)
    parser.add_argument("--sjdb-overhang", type=int, required=True)
    args = parser.parse_args()

    star_version = subprocess.run(
        [str(args.star), "--version"], check=True, text=True, capture_output=True
    ).stdout.strip()
    index_files = sorted(
        path for path in args.index.iterdir() if path.is_file() and path != args.output
    )
    manifest = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "reference": "Mus musculus GRCm38, Ensembl release 102",
        "fasta": {
            **described(args.fasta),
            "url": args.fasta_url,
            "download_gzip_sha256": args.fasta_gzip_sha256,
        },
        "gtf": {
            **described(args.gtf),
            "url": args.gtf_url,
            "download_gzip_sha256": args.gtf_gzip_sha256,
        },
        "STAR": {
            "version": star_version,
            "binary": described(args.star),
            "source_url": args.star_source_url,
            "source_sha256": args.star_source_sha256,
        },
        "genome_generate": {
            "sjdbOverhang": args.sjdb_overhang,
            "command_template": [
                "STAR", "--runMode", "genomeGenerate", "--runThreadN", "32",
                "--genomeDir", "INDEX", "--genomeFastaFiles", "FASTA",
                "--sjdbGTFfile", "GTF", "--sjdbOverhang", str(args.sjdb_overhang),
            ],
        },
        "index_files": [described(path) for path in index_files],
    }
    args.output.write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
