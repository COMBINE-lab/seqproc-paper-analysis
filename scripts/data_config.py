#!/usr/bin/env python3
"""
Centralized data path resolution for all benchmark scripts.

Provides a single source of truth for FASTQ file paths across different
dataset sizes (1m, full). Scripts import resolve_datasets() to get
the correct file paths for the requested --reads level.

Usage from other scripts:
    from data_config import resolve_datasets, add_reads_arg

    # In argparse setup:
    add_reads_arg(parser)

    # After parsing:
    datasets = resolve_datasets(args.reads)
"""

import os
import shutil
from pathlib import Path
from typing import Dict

PROJECT_ROOT = Path(os.environ.get("SEQPROC_PROJECT_ROOT", Path(__file__).parent.parent))
DATA_DIR = Path(os.environ.get("SEQPROC_DATA_DIR", PROJECT_ROOT / "data"))
CONFIGS = PROJECT_ROOT / "configs"

# ============================================================================
# Tool binary auto-discovery (shared by all benchmark scripts)
# ============================================================================

def _find_binary(env_var: str, name: str, relative_paths: list) -> str:
    """Resolve a tool binary from env var, known relative paths, or PATH."""
    # 1. Explicit env var (highest priority)
    val = os.environ.get(env_var)
    if val and os.path.isfile(val) and os.access(val, os.X_OK):
        return val
    # 2. Search candidate directories
    _data_dir = os.environ.get("SEQPROC_DATA_DIR", "")
    search_roots = [PROJECT_ROOT.parent]
    if _data_dir:
        search_roots.append(Path(_data_dir).parent)
    for root in search_roots:
        for rp in relative_paths:
            candidate = root / rp
            if candidate.is_file() and os.access(str(candidate), os.X_OK):
                return str(candidate)
    # 3. Fall back to PATH
    found = shutil.which(name)
    if found:
        return found
    # 4. Return best-guess default for a clear error later.
    return val or str(PROJECT_ROOT.parent / relative_paths[0])


def resolve_binaries() -> Dict[str, str]:
    """Return dict with resolved paths for seqproc, matchbox, splitcode."""
    return {
        "seqproc": _find_binary("SEQPROC_BIN", "seqproc",
            ["combine-lab/seqproc/target/release/seqproc"]),
        "matchbox": _find_binary("MATCHBOX_BIN", "matchbox",
            ["matchbox/target/release/matchbox"]),
        "splitcode": _find_binary("SPLITCODE_BIN", "splitcode",
            ["splitcode/build/src/splitcode"]),
    }

# ============================================================================
# SRA accessions and full read counts (from SRA metadata)
# ============================================================================

SRA_INFO = {
    "SRR6750041": {
        "description": "SPLiT-seq PE (short-read paired-end)",
        "full_reads": 86_820_578,
        "download_cmd": "fasterq-dump --split-files SRR6750041",
    },
    "SRR13948564": {
        "description": "LR-SPLiT-seq (PacBio Sequel II long-read)",
        "full_reads": 4_229_250,
        "download_cmd": "fasterq-dump SRR13948564",
    },
    "SRR8315379": {
        "description": "10x Chromium v2 (short-read paired-end)",
        "full_reads": 56_514_800,
        "download_cmd": "fasterq-dump --split-files SRR8315379",
    },
    "SRR7827254": {
        "description": "sci-RNA-seq3 (short-read paired-end)",
        "full_reads": 10_177_866,
        "download_cmd": "fasterq-dump --split-files SRR7827254",
    },
}


def _p(*parts):
    """Join parts into a Path relative to DATA_DIR."""
    return DATA_DIR.joinpath(*parts)


# Each entry: { "1m": (r1, r2_or_None, reads), "full": (r1, r2_or_None, reads) }
_FILE_MAP = {
    "splitseq_pe": {
        "1m": (
            _p("SRR6750041_1M_R1.fastq"),
            _p("SRR6750041_1M_R2.fastq"),
            1_000_000,
        ),
        "full": (
            _p("SRR6750041_R1.fastq"),
            _p("SRR6750041_R2.fastq"),
            SRA_INFO["SRR6750041"]["full_reads"],
        ),
    },
    "lr_splitseq": {
        "1m": (
            _p("SRR13948564_1M.fastq"),
            None,
            1_000_000,
        ),
        "full": (
            _p("SRR13948564_full.fastq"),
            None,
            SRA_INFO["SRR13948564"]["full_reads"],
        ),
    },
    "10x_short": {
        "1m": (
            _p("10x_short/SRR8315379_1M_R1.fastq"),
            _p("10x_short/SRR8315379_1M_R2.fastq"),
            1_000_000,
        ),
        "full": (
            _p("10x_short/SRR8315379_R1.fastq"),
            _p("10x_short/SRR8315379_R2.fastq"),
            SRA_INFO["SRR8315379"]["full_reads"],
        ),
    },
    "sciseq": {
        "1m": (
            _p("SRR7827254_1M_1.fastq"),
            _p("SRR7827254_1M_2.fastq"),
            1_000_000,
        ),
        "full": (
            _p("SRR7827254_1.fastq"),
            _p("SRR7827254_2.fastq"),
            SRA_INFO["SRR7827254"]["full_reads"],
        ),
    },
}

# ============================================================================
# Tool config definitions (shared between benchmark + concordance scripts)
# ============================================================================

TOOL_CONFIGS = {
    "splitseq_pe": {
        "name": "SPLiT-seq PE",
        "short_name": "SPLiT-seq PE",
        "mode": "paired",
        "seqproc_edit_geom": CONFIGS / "seqproc/splitseq_filter_edit.geom",
        "seqproc_hamming_geom": CONFIGS / "seqproc/splitseq_filter_hamming6.geom",
        "seqproc_geom": CONFIGS / "seqproc/splitseq_filter_edit.geom",
        "seqproc_maps": [
            CONFIGS / "seqproc/splitseq_bc3_seq2seq.tsv",
            CONFIGS / "seqproc/splitseq_bc2_seq2seq.tsv",
            CONFIGS / "seqproc/splitseq_bc1_seq2seq.tsv",
        ],
        "matchbox_config": CONFIGS / "matchbox/splitseq_replacement.mb",
        "matchbox_paired": True,
        "splitcode_config": CONFIGS / "splitcode/splitseq_paper.config",
        "tools": ["seqproc", "matchbox", "splitcode"],
    },
    "lr_splitseq": {
        "name": "LR-SPLiT-seq",
        "short_name": "LR-SPLiT-seq",
        "mode": "single",
        "seqproc_edit_geom": CONFIGS / "seqproc/splitseq_singleend_edit_ann.geom",
        "seqproc_hamming_geom": CONFIGS / "seqproc/splitseq_singleend_ann.geom",
        "seqproc_fw_edit_geom": CONFIGS / "seqproc/splitseq_singleend_edit.geom",
        "seqproc_fw_hamming_geom": CONFIGS / "seqproc/splitseq_singleend.geom",
        "seqproc_geom": CONFIGS / "seqproc/splitseq_singleend_primer_edit.geom",
        "seqproc_maps": [
            CONFIGS / "seqproc/splitseq_bc3_seq2seq.tsv",
            CONFIGS / "seqproc/splitseq_bc2_seq2seq.tsv",
            CONFIGS / "seqproc/splitseq_bc1_seq2seq.tsv",
        ],
        "matchbox_config": CONFIGS / "matchbox/splitseq_singleend_dual.mb",
        "matchbox_paired": False,
        "splitcode_config": CONFIGS / "splitcode/splitseq_singleend.config",
        "tools": ["seqproc", "matchbox", "splitcode"],
    },
    "10x_short": {
        "name": "10x Chromium v2 Short Read",
        "short_name": "10x Short",
        "mode": "paired",
        "seqproc_edit_geom": CONFIGS / "seqproc/10x_v2.geom",
        "seqproc_hamming_geom": None,
        "seqproc_geom": CONFIGS / "seqproc/10x_v2.geom",
        "matchbox_config": CONFIGS / "matchbox/10x_v2.mb",
        "matchbox_paired": False,
        "splitcode_config": CONFIGS / "splitcode/10x_v2.config",
        "tools": ["seqproc", "matchbox", "splitcode"],
    },
    "sciseq": {
        "name": "sci-RNA-seq3",
        "short_name": "sci-RNA-seq3",
        "mode": "paired",
        "seqproc_edit_geom": CONFIGS / "seqproc/sciseq3_edit.geom",
        "seqproc_hamming_geom": CONFIGS / "seqproc/sciseq3.geom",
        "seqproc_geom": CONFIGS / "seqproc/sciseq3_edit.geom",
        "matchbox_config": CONFIGS / "matchbox/sciseq3.mb",
        "matchbox_paired": True,
        "splitcode_config": CONFIGS / "splitcode/sciseq3.config",
        "tools": ["seqproc", "matchbox", "splitcode"],
    },
}

# The four Table 2 datasets (excludes supplementary like gridion/promethion)
TABLE2_DATASETS = ["splitseq_pe", "lr_splitseq", "10x_short", "sciseq"]


def resolve_datasets(reads_level: str = "1m") -> Dict:
    """Resolve dataset configurations for the given reads level.

    Args:
        reads_level: One of "1m" or "full".

    Returns:
        Dict mapping dataset key -> merged config dict with r1, r2, reads fields.
    """
    if reads_level not in ("1m", "full"):
        raise ValueError(
            f"Unknown reads level: {reads_level!r}. Use '1m' or 'full'."
        )

    resolved = {}
    for ds_key in TABLE2_DATASETS:
        cfg = dict(TOOL_CONFIGS[ds_key])
        r1, r2, reads = _FILE_MAP[ds_key][reads_level]
        cfg["r1"] = r1
        cfg["r2"] = r2
        cfg["reads"] = reads
        resolved[ds_key] = cfg

    return resolved


def check_data_availability(datasets: Dict) -> Dict[str, bool]:
    """Check which datasets have their FASTQ files available.

    Returns dict mapping dataset key -> True if all required files exist.
    """
    status = {}
    for ds_key, cfg in datasets.items():
        r1_ok = cfg["r1"].exists()
        r2_ok = cfg["r2"] is None or cfg["r2"].exists()
        status[ds_key] = r1_ok and r2_ok
    return status


def print_data_status(datasets: Dict):
    """Print a table showing which datasets are available."""
    status = check_data_availability(datasets)
    print(f"{'Dataset':<20} {'R1':<50} {'Status'}")
    print("-" * 80)
    for ds_key, cfg in datasets.items():
        ok = "[OK]" if status[ds_key] else "[MISSING]"
        r1_name = str(cfg["r1"].relative_to(DATA_DIR))
        reads_str = f"{cfg['reads']:,}"
        print(
            f"{cfg['short_name']:<20} {r1_name:<50} {ok} ({reads_str} reads)"
        )


def add_reads_arg(parser):
    """Add the --reads argument to an argparse parser."""
    parser.add_argument(
        "--reads",
        type=str,
        choices=["1m", "full"],
        default="1m",
        help=(
            "Dataset size: '1m' for 1M-read subsets (default), "
            "'full' for complete SRA datasets"
        ),
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Check data availability for benchmarks"
    )
    add_reads_arg(parser)
    args = parser.parse_args()

    datasets = resolve_datasets(args.reads)
    print(f"Data status for --reads={args.reads}:\n")
    print_data_status(datasets)

    status = check_data_availability(datasets)
    missing = [k for k, v in status.items() if not v]
    if missing:
        print(f"\nMissing {len(missing)} dataset(s). Download with:")
        print("  cd data/")
        for ds_key in missing:
            accession = {
                "splitseq_pe": "SRR6750041",
                "lr_splitseq": "SRR13948564",
                "10x_short": "SRR8315379",
                "sciseq": "SRR7827254",
            }.get(ds_key, "???")
            cmd = SRA_INFO.get(accession, {}).get(
                "download_cmd", f"fasterq-dump {accession}"
            )
            print(f"  {cmd}")
    else:
        print("\nAll datasets available.")
