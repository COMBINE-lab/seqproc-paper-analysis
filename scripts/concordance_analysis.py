#!/usr/bin/env python3
"""
Concordance Analysis & Hamming vs Edit Comparison

Runs all 3 tools on all 4 Table 2 datasets, extracts read ID sets, computes
pairwise Jaccard indices, characterizes discordant reads, and compares
seqproc hamming vs edit distance recovery.

Usage:
    python scripts/concordance_analysis.py --threads 4
    python scripts/concordance_analysis.py --threads 4 --datasets splitseq_pe lr_splitseq
    python scripts/concordance_analysis.py --threads 4 --skip-runs  # use existing outputs
"""

import subprocess
import shutil
import time
import os
import json
import argparse
from pathlib import Path
from typing import Dict, Set, Tuple

# ============================================================================
# Configuration
# ============================================================================

PROJECT_ROOT = Path(os.environ.get("SEQPROC_PROJECT_ROOT", Path(__file__).parent.parent))
RESULTS_DIR = PROJECT_ROOT / "results" / "concordance"

CONFIGS = PROJECT_ROOT / "configs"

# --------------------------------------------------------------------------
# Dataset definitions -- loaded from centralized data_config module.
# The DATASETS dict is populated at startup by _init_datasets() so that
# existing code (including tests that import DATASETS) continues to work.
# When run as __main__, --reads controls the dataset size.
# --------------------------------------------------------------------------

from data_config import resolve_datasets, add_reads_arg, TOOL_CONFIGS, resolve_binaries, DATA_DIR

_bins = resolve_binaries()
SEQPROC_BIN = _bins["seqproc"]
MATCHBOX_BIN = _bins["matchbox"]
SPLITCODE_BIN = _bins["splitcode"]

# Temp directory base: use scratch space instead of /tmp which fills up
_data_dir = os.environ.get("SEQPROC_DATA_DIR", "")
TMPDIR_BASE = str(Path(_data_dir).parent / "tmp") if _data_dir else None

# Default to 1M subsets; overridden in main() when --reads is specified.
DATASETS = resolve_datasets("1m")


# ============================================================================
# Helper Functions
# ============================================================================

def run_cmd(cmd: str, cwd=None) -> Tuple[float, float, int, str]:
    """Run command, return (runtime_s, peak_mem_mb, returncode, stderr)."""
    env = None
    if TMPDIR_BASE:
        os.makedirs(TMPDIR_BASE, exist_ok=True)
        env = {**os.environ, "TMPDIR": TMPDIR_BASE}
    time_cmd = f"/usr/bin/time -v {cmd}"
    start = time.time()
    result = subprocess.run(time_cmd, shell=True, capture_output=True, text=True, cwd=cwd, env=env)
    runtime = time.time() - start

    peak_mem_kb = 0
    tool_exit_code = result.returncode
    for line in result.stderr.split('\n'):
        if 'Maximum resident set size' in line:
            peak_mem_kb = int(line.split(':')[1].strip())
        elif 'Exit status' in line:
            tool_exit_code = int(line.split(':')[1].strip())

    return runtime, peak_mem_kb / 1024, tool_exit_code, result.stderr


def extract_fastq_ids(filepath: str) -> Set[str]:
    """Extract read IDs from a FASTQ file."""
    ids = set()
    if not os.path.exists(filepath):
        return ids
    with open(filepath, 'r') as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                rid = line.strip().split()[0].lstrip('@')
                ids.add(rid)
    return ids


def count_fastq_reads(filepath: str) -> int:
    """Count reads in a FASTQ file (lines / 4)."""
    if not os.path.exists(filepath):
        return 0
    with open(filepath, 'r') as f:
        lines = sum(1 for _ in f)
    return lines // 4


def extract_tsv_ids(filepath: str) -> Set[str]:
    """Extract read IDs from first column of a TSV file."""
    ids = set()
    if not os.path.exists(filepath):
        return ids
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if parts:
                ids.add(parts[0])
    return ids


def jaccard(a: Set[str], b: Set[str]) -> float:
    """Compute Jaccard index between two sets."""
    if not a and not b:
        return 1.0
    union = len(a | b)
    if union == 0:
        return 1.0
    return len(a & b) / union


def save_id_set(ids: Set[str], filepath: str):
    """Save sorted set of IDs to file."""
    with open(filepath, 'w') as f:
        for rid in sorted(ids):
            f.write(rid + '\n')


# ============================================================================
# Tool Runners (return read ID sets)
# ============================================================================

def run_seqproc(dataset: dict, geom_path: Path, outdir: Path,
                threads: int, label: str) -> Tuple[Set[str], float, float]:
    """Run seqproc and return (read_id_set, runtime, memory_mb)."""
    out1 = outdir / f"seqproc_{label}_R1.fq"
    out2 = outdir / f"seqproc_{label}_R2.fq"

    if dataset['mode'] == 'single':
        cmd = (f"{SEQPROC_BIN} --geom {geom_path} "
               f"--file1 {dataset['r1']} --out1 {out1} --threads {threads}")
    else:
        cmd = (f"{SEQPROC_BIN} --geom {geom_path} "
               f"--file1 {dataset['r1']} --file2 {dataset['r2']} "
               f"--out1 {out1} --out2 {out2} --threads {threads}")

    # Add map files if present
    if 'seqproc_maps' in dataset:
        for m in dataset['seqproc_maps']:
            cmd += f" -a {m}"

    runtime, mem, rc, stderr = run_cmd(cmd, PROJECT_ROOT)
    if rc != 0:
        print(f"    [ERROR] seqproc {label} failed (rc={rc})")
        print(f"    stderr: {stderr[-500:]}")
        return set(), runtime, mem

    # Extract IDs from the appropriate output
    if dataset['mode'] == 'paired':
        ids = extract_fastq_ids(str(out2))
    else:
        ids = extract_fastq_ids(str(out1))

    return ids, runtime, mem


def run_matchbox(dataset: dict, outdir: Path,
                 threads: int) -> Tuple[Set[str], float, float]:
    """Run matchbox and return (read_id_set, runtime, memory_mb)."""
    out_tsv = outdir / "matchbox_out.tsv"

    # Clean stale matchbox output files in outdir
    for fq in ['mb_r1.fq', 'mb_r2.fq']:
        stale = outdir / fq
        if stale.exists():
            stale.unlink()

    # Symlink configs into outdir so matchbox csv() relative paths resolve
    # when we run from outdir (on /fs/) instead of PROJECT_ROOT (NFS).
    configs_link = outdir / "configs"
    if not configs_link.exists():
        configs_link.symlink_to(CONFIGS.parent / "configs")

    # Build args using explicit matchbox_paired flag (not name-based dispatch)
    if dataset.get('matchbox_paired', False):
        args = f"{dataset['r1']} -p {dataset['r2']}"
    else:
        args = str(dataset['r1'])

    cmd = (f"{MATCHBOX_BIN} -e 0.2 -t {threads} "
           f"-s {dataset['matchbox_config']} {args} > {out_tsv}")

    # Run from outdir so .out!() FASTQ files land on scratch, not NFS
    runtime, mem, rc, stderr = run_cmd(cmd, str(outdir))
    if rc != 0:
        print(f"    [ERROR] matchbox failed (rc={rc})")
        print(f"    stderr: {stderr[-500:]}")
        return set(), runtime, mem

    # Collect FASTQ files (.out!() now writes directly to outdir)
    found = {}
    for fq_name in ['mb_r1.fq', 'mb_r2.fq']:
        fq_path = outdir / fq_name
        if fq_path.exists():
            found[fq_name] = fq_path

    # Extract IDs: prefer R2 for paired-end (consistent with seqproc/splitcode)
    ids = set()
    if found:
        prefer = 'mb_r2.fq' if dataset['mode'] == 'paired' else 'mb_r1.fq'
        id_src = found.get(prefer) or next(iter(found.values()))
        ids = extract_fastq_ids(str(id_src))

    # Fallback to TSV
    if not ids:
        ids = extract_tsv_ids(str(out_tsv))

    return ids, runtime, mem


def run_splitcode(dataset: dict, outdir: Path,
                  threads: int) -> Tuple[Set[str], float, float]:
    """Run splitcode and return (read_id_set, runtime, memory_mb)."""
    # Detect assign flag
    with open(dataset['splitcode_config']) as f:
        config_text = f.read()
    has_tags = False
    for line in config_text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith('#') or stripped.startswith('@'):
            continue
        fields = stripped.split('\t')
        if len(fields) >= 3 and fields[0] not in ('ID', 'group'):
            has_tags = True
            break

    assign_flag = "--assign" if has_tags else ""
    mapping = outdir / "splitcode_mapping.txt"
    mapping_flag = f"-m {mapping}" if has_tags else ""

    if dataset['mode'] == 'single':
        out_fq = outdir / "splitcode_out.fq"
        cmd = (f"{SPLITCODE_BIN} -c {dataset['splitcode_config']} "
               f"{assign_flag} -N 1 -t {threads} {mapping_flag} "
               f"-o {out_fq} {dataset['r1']}")
    else:
        out1 = outdir / "splitcode_R1.fq"
        out2 = outdir / "splitcode_R2.fq"
        cmd = (f"{SPLITCODE_BIN} -c {dataset['splitcode_config']} "
               f"{assign_flag} -N 2 -t {threads} {mapping_flag} "
               f"-o {out1},{out2} {dataset['r1']} {dataset['r2']}")

    runtime, mem, rc, stderr = run_cmd(cmd, PROJECT_ROOT)
    if rc != 0:
        print(f"    [ERROR] splitcode failed (rc={rc})")
        print(f"    stderr: {stderr[-500:]}")
        return set(), runtime, mem

    if dataset['mode'] == 'single':
        ids = extract_fastq_ids(str(out_fq))
    else:
        ids = extract_fastq_ids(str(outdir / "splitcode_R2.fq"))

    return ids, runtime, mem


# ============================================================================
# Concordance Analysis
# ============================================================================

def compute_concordance(id_sets: Dict[str, Set[str]]) -> Dict:
    """Compute pairwise concordance for all tool pairs."""
    tools = sorted(id_sets.keys())
    result = {
        "tools": {},
        "pairwise": []
    }

    for t in tools:
        result["tools"][t] = len(id_sets[t])

    for i in range(len(tools)):
        for j in range(i + 1, len(tools)):
            t1, t2 = tools[i], tools[j]
            s1, s2 = id_sets[t1], id_sets[t2]
            both = len(s1 & s2)
            only1 = len(s1 - s2)
            only2 = len(s2 - s1)
            jac = jaccard(s1, s2)

            result["pairwise"].append({
                "tool_a": t1,
                "tool_b": t2,
                "both": both,
                f"{t1}_only": only1,
                f"{t2}_only": only2,
                "jaccard": round(jac, 4),
                "union": both + only1 + only2
            })

    return result


def characterize_discordant(
    id_sets: Dict[str, Set[str]],
    dataset: dict,
    outdir: Path
) -> Dict:
    """Characterize reads found by one tool but not others.
    Returns summary stats about discordant reads."""

    tools = sorted(id_sets.keys())
    disc_info = {}

    # For each tool, find reads unique to that tool (not in any other)
    for t in tools:
        others = set()
        for t2 in tools:
            if t2 != t:
                others |= id_sets[t2]
        unique_to_t = id_sets[t] - others
        disc_info[f"{t}_unique"] = len(unique_to_t)

        # Save unique IDs for later investigation
        if unique_to_t:
            save_id_set(unique_to_t, str(outdir / f"unique_to_{t}.txt"))

    # All-tool consensus (reads found by ALL tools)
    if tools:
        consensus = id_sets[tools[0]]
        for t in tools[1:]:
            consensus = consensus & id_sets[t]
        disc_info["all_tools_consensus"] = len(consensus)

    # Any-tool union
    all_union = set()
    for t in tools:
        all_union |= id_sets[t]
    disc_info["any_tool_union"] = len(all_union)

    return disc_info


# ============================================================================
# Main Analysis
# ============================================================================

def run_dataset_analysis(
    ds_key: str,
    dataset: dict,
    outdir: Path,
    threads: int,
    skip_runs: bool = False
) -> Dict:
    """Run complete analysis for one dataset."""
    ds_outdir = outdir / ds_key
    ds_outdir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"Dataset: {dataset['name']} ({ds_key})")
    print(f"{'='*60}")

    # Check data exists
    if not os.path.exists(dataset['r1']):
        print(f"  [SKIP] Data not found: {dataset['r1']}")
        return {}

    # ----------------------------------------------------------------
    # Step 1: Run all tools with EDIT distance configs (primary)
    # ----------------------------------------------------------------
    id_sets = {}
    perf = {}

    # --- seqproc (edit distance) ---
    ids_file = ds_outdir / "seqproc_edit_ids.txt"
    if skip_runs and ids_file.exists():
        print(f"  seqproc (edit): loading cached IDs...")
        with open(ids_file) as f:
            id_sets["seqproc"] = {line.strip() for line in f}
        perf["seqproc_edit"] = {"runtime": 0, "memory_mb": 0}
    else:
        print(f"  seqproc (edit)...", end=" ", flush=True)
        ids, rt, mem = run_seqproc(dataset, dataset['seqproc_edit_geom'],
                                    ds_outdir, threads, "edit")
        id_sets["seqproc"] = ids
        perf["seqproc_edit"] = {"runtime": round(rt, 2), "memory_mb": round(mem, 1)}
        save_id_set(ids, str(ids_file))
        print(f"{len(ids):,} reads in {rt:.1f}s, {mem:.0f}MB")

    # --- matchbox ---
    ids_file = ds_outdir / "matchbox_ids.txt"
    if skip_runs and ids_file.exists():
        print(f"  matchbox: loading cached IDs...")
        with open(ids_file) as f:
            id_sets["matchbox"] = {line.strip() for line in f}
        perf["matchbox"] = {"runtime": 0, "memory_mb": 0}
    else:
        print(f"  matchbox...", end=" ", flush=True)
        ids, rt, mem = run_matchbox(dataset, ds_outdir, threads)
        id_sets["matchbox"] = ids
        perf["matchbox"] = {"runtime": round(rt, 2), "memory_mb": round(mem, 1)}
        save_id_set(ids, str(ids_file))
        print(f"{len(ids):,} reads in {rt:.1f}s, {mem:.0f}MB")

    # --- splitcode ---
    ids_file = ds_outdir / "splitcode_ids.txt"
    if skip_runs and ids_file.exists():
        print(f"  splitcode: loading cached IDs...")
        with open(ids_file) as f:
            id_sets["splitcode"] = {line.strip() for line in f}
        perf["splitcode"] = {"runtime": 0, "memory_mb": 0}
    else:
        print(f"  splitcode...", end=" ", flush=True)
        ids, rt, mem = run_splitcode(dataset, ds_outdir, threads)
        id_sets["splitcode"] = ids
        perf["splitcode"] = {"runtime": round(rt, 2), "memory_mb": round(mem, 1)}
        save_id_set(ids, str(ids_file))
        print(f"{len(ids):,} reads in {rt:.1f}s, {mem:.0f}MB")

    # ----------------------------------------------------------------
    # Step 2: Concordance
    # ----------------------------------------------------------------
    print(f"\n  --- Concordance ---")
    concordance = compute_concordance(id_sets)

    for pair in concordance["pairwise"]:
        t1, t2 = pair["tool_a"], pair["tool_b"]
        print(f"  {t1} vs {t2}: "
              f"Jaccard={pair['jaccard']:.4f}, "
              f"both={pair['both']:,}, "
              f"{t1}-only={pair[f'{t1}_only']:,}, "
              f"{t2}-only={pair[f'{t2}_only']:,}")

    # ----------------------------------------------------------------
    # Step 3: Discordant read characterization
    # ----------------------------------------------------------------
    print(f"\n  --- Discordant Reads ---")
    discordant = characterize_discordant(id_sets, dataset, ds_outdir)

    print(f"  All-tool consensus: {discordant.get('all_tools_consensus', 0):,}")
    print(f"  Any-tool union:     {discordant.get('any_tool_union', 0):,}")
    for tool in sorted(id_sets.keys()):
        unique_count = discordant.get(f"{tool}_unique", 0)
        if unique_count > 0:
            print(f"  Unique to {tool}: {unique_count:,}")

    # ----------------------------------------------------------------
    # Hamming vs Edit comparison (seqproc only)
    # ----------------------------------------------------------------
    hamming_edit = {}
    hamming_geom = dataset.get('seqproc_hamming_geom')
    if hamming_geom and os.path.exists(hamming_geom):
        print(f"\n  --- seqproc Hamming vs Edit ---")
        ids_file = ds_outdir / "seqproc_hamming_ids.txt"
        if skip_runs and ids_file.exists():
            print(f"  seqproc (hamming): loading cached IDs...")
            with open(ids_file) as f:
                hamming_ids = {line.strip() for line in f}
            ham_rt, ham_mem = 0, 0
        else:
            print(f"  seqproc (hamming)...", end=" ", flush=True)
            hamming_ids, ham_rt, ham_mem = run_seqproc(
                dataset, hamming_geom, ds_outdir, threads, "hamming")
            save_id_set(hamming_ids, str(ids_file))
            print(f"{len(hamming_ids):,} reads in {ham_rt:.1f}s, {ham_mem:.0f}MB")

        edit_ids = id_sets["seqproc"]

        both = len(hamming_ids & edit_ids)
        ham_only = len(hamming_ids - edit_ids)
        edit_only = len(edit_ids - hamming_ids)
        jac = jaccard(hamming_ids, edit_ids)

        hamming_edit = {
            "hamming_reads": len(hamming_ids),
            "edit_reads": len(edit_ids),
            "both": both,
            "hamming_only": ham_only,
            "edit_only": edit_only,
            "edit_gain": edit_only - ham_only,
            "edit_gain_pct": round((edit_only - ham_only) / max(len(hamming_ids), 1) * 100, 2),
            "jaccard": round(jac, 4),
            "hamming_runtime": round(ham_rt, 2),
            "hamming_memory_mb": round(ham_mem, 1),
        }

        print(f"  Hamming: {len(hamming_ids):,} reads")
        print(f"  Edit:    {len(edit_ids):,} reads")
        print(f"  Edit gains: +{edit_only:,} reads ({hamming_edit['edit_gain_pct']}% improvement)")
        print(f"  Hamming-only: {ham_only:,}, Edit-only: {edit_only:,}, Both: {both:,}")
        print(f"  Jaccard: {jac:.4f}")

        # Save edit-only IDs for characterization
        edit_only_ids = edit_ids - hamming_ids
        if edit_only_ids:
            save_id_set(edit_only_ids, str(ds_outdir / "edit_only_ids.txt"))

    else:
        if hamming_geom is None:
            print(f"\n  --- seqproc Hamming vs Edit: N/A (no anchor matching) ---")
        else:
            print(f"\n  --- seqproc Hamming vs Edit: SKIPPED (config not found) ---")

    # ----------------------------------------------------------------
    # Assemble results
    # ----------------------------------------------------------------
    # Count actual input reads from FASTQ (SRA metadata can be wrong)
    actual_input = count_fastq_reads(str(dataset['r1']))
    if actual_input == 0:
        actual_input = dataset["reads"]
    print(f"\n  Input reads (from FASTQ): {actual_input:,}")

    result = {
        "name": dataset["name"],
        "total_reads": actual_input,
        "performance": perf,
        "recovery": {t: len(ids) for t, ids in id_sets.items()},
        "recovery_pct": {t: round(len(ids) / actual_input * 100, 2)
                         for t, ids in id_sets.items()},
        "concordance": concordance,
        "discordant": discordant,
        "hamming_vs_edit": hamming_edit,
    }

    # Save per-dataset results
    with open(ds_outdir / "results.json", 'w') as f:
        json.dump(result, f, indent=2)

    return result


def print_summary_tables(all_results: Dict):
    """Print formatted summary tables."""
    print(f"\n{'='*70}")
    print("SUMMARY: CONCORDANCE MATRIX")
    print(f"{'='*70}")

    for ds_key, res in all_results.items():
        if not res:
            continue
        print(f"\n  {res['name']} ({res['total_reads']:,} reads)")
        print(f"  {'─'*55}")

        # Recovery
        for tool, count in res["recovery"].items():
            pct = res["recovery_pct"][tool]
            print(f"    {tool:12s}: {count:>10,} reads ({pct:5.1f}%)")

        # Pairwise concordance
        print(f"\n  Pairwise Jaccard:")
        for pair in res["concordance"]["pairwise"]:
            t1, t2 = pair["tool_a"], pair["tool_b"]
            print(f"    {t1:12s} vs {t2:12s}: "
                  f"J={pair['jaccard']:.4f}  "
                  f"(both={pair['both']:,}, "
                  f"{t1}-only={pair[f'{t1}_only']:,}, "
                  f"{t2}-only={pair[f'{t2}_only']:,})")

    # Hamming vs Edit summary
    print(f"\n{'='*70}")
    print("SUMMARY: HAMMING vs EDIT DISTANCE (seqproc)")
    print(f"{'='*70}")

    for ds_key, res in all_results.items():
        if not res or not res.get("hamming_vs_edit"):
            continue
        hve = res["hamming_vs_edit"]
        print(f"\n  {res['name']}:")
        print(f"    Hamming: {hve['hamming_reads']:>10,} reads")
        print(f"    Edit:    {hve['edit_reads']:>10,} reads")
        print(f"    Gain:    +{hve['edit_gain']:>9,} reads ({hve['edit_gain_pct']:+.1f}%)")
        print(f"    Jaccard: {hve['jaccard']:.4f}")

    # Discordant summary
    print(f"\n{'='*70}")
    print("SUMMARY: DISCORDANT READS")
    print(f"{'='*70}")

    for ds_key, res in all_results.items():
        if not res:
            continue
        disc = res["discordant"]
        print(f"\n  {res['name']}:")
        print(f"    Consensus (all tools): {disc.get('all_tools_consensus', 0):,}")
        print(f"    Union (any tool):      {disc.get('any_tool_union', 0):,}")
        for key, val in disc.items():
            if key.endswith("_unique") and val > 0:
                tool = key.replace("_unique", "")
                print(f"    Unique to {tool:12s}: {val:,}")


def main():
    global DATASETS
    parser = argparse.ArgumentParser(description='Concordance Analysis')
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--datasets', type=str, nargs='+', default=None,
                        help='Dataset keys to run (default: all)')
    parser.add_argument('--skip-runs', action='store_true',
                        help='Skip tool runs, use cached ID files')
    parser.add_argument('--output', type=str, default=None)
    add_reads_arg(parser)
    args = parser.parse_args()

    # Re-resolve datasets for the requested reads level
    DATASETS = resolve_datasets(args.reads)

    outdir = Path(args.output) if args.output else RESULTS_DIR
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("CONCORDANCE ANALYSIS & HAMMING vs EDIT COMPARISON")
    print("=" * 70)
    print(f"Threads: {args.threads}")
    print(f"Output:  {outdir}")
    if args.skip_runs:
        print("Mode:    Using cached results (--skip-runs)")

    # Determine which datasets to run
    ds_keys = args.datasets if args.datasets else list(DATASETS.keys())
    print(f"Datasets: {', '.join(ds_keys)}")

    all_results = {}
    for ds_key in ds_keys:
        if ds_key not in DATASETS:
            print(f"\n[WARNING] Unknown dataset: {ds_key}")
            continue
        result = run_dataset_analysis(
            ds_key, DATASETS[ds_key], outdir, args.threads, args.skip_runs)
        all_results[ds_key] = result

    # Save combined results
    combined_path = outdir / "concordance_results.json"
    with open(combined_path, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\nSaved combined results: {combined_path}")

    # Print summary tables
    print_summary_tables(all_results)

    print(f"\n{'='*70}")
    print("CONCORDANCE ANALYSIS COMPLETE")
    print(f"{'='*70}")
    print(f"All results in: {outdir}")


if __name__ == "__main__":
    main()
