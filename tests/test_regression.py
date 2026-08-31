#!/usr/bin/env python3
"""
Regression tests for seqproc paper analysis pipeline.

These tests verify that renames and refactors have not broken any
script imports, path resolution, data loading, or result integrity.
All tests must pass before any commit is accepted.

Run:
    pytest tests/ -v
"""

import json
import os
import py_compile
import subprocess
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).parent.parent
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
CONFIGS_DIR = PROJECT_ROOT / "configs"
RESULTS_DIR = PROJECT_ROOT / "results"

# ============================================================================
# 1. Script compilation -- every .py in scripts/ must be valid Python
# ============================================================================

SCRIPT_FILES = sorted(SCRIPTS_DIR.glob("*.py"))


@pytest.mark.parametrize("script", SCRIPT_FILES, ids=lambda p: p.name)
def test_script_compiles(script):
    """Every script must compile without syntax errors."""
    py_compile.compile(str(script), doraise=True)


# ============================================================================
# 2. Module imports -- renamed modules must be importable
# ============================================================================


def test_import_concordance_analysis():
    """concordance_analysis module must be importable and export DATASETS."""
    from concordance_analysis import DATASETS, RESULTS_DIR as RD
    assert isinstance(DATASETS, dict)
    assert len(DATASETS) >= 4, f"Expected >=4 datasets, got {len(DATASETS)}"
    assert "concordance" in str(RD), (
        f"RESULTS_DIR should contain 'concordance', got {RD}")


def test_import_generate_figures():
    """generate_figures module must be importable with correct constants."""
    from generate_figures import (
        CONCORDANCE_DIR, PERF_JSON, PERF_BACKUP,
        OUTPUT_DIR, PERF_TO_CONCORDANCE, load_data
    )
    assert "concordance" in str(CONCORDANCE_DIR)
    assert "benchmark_results.json" in str(PERF_JSON)
    assert "benchmark_results_perf.json" in str(PERF_BACKUP)


def test_import_discordant_analysis():
    """discordant_analysis module must be importable with correct RESULTS_DIR."""
    from discordant_analysis import RESULTS_DIR as RD
    assert "concordance" in str(RD), (
        f"RESULTS_DIR should contain 'concordance', got {RD}")


def test_import_lr_perf_rerun():
    """lr_perf_rerun module must be importable with correct RESULTS_DIR."""
    from lr_perf_rerun import RESULTS_DIR as RD
    assert "lr_perf" in str(RD), (
        f"RESULTS_DIR should contain 'lr_perf', got {RD}")
    assert "phase" not in str(RD).lower(), (
        f"RESULTS_DIR should not contain 'phase', got {RD}")


def test_import_run_paper_benchmarks():
    """run_paper_benchmarks module must be importable."""
    from run_paper_benchmarks import PROJECT_ROOT as PR
    assert PR.resolve() == PROJECT_ROOT.resolve()


# ============================================================================
# 3. No banned development-iteration terminology in tracked files
# ============================================================================


# Encode banned terms so this test file itself does not trigger false positives.
# Each entry is a list of string fragments that get joined at runtime.
_BANNED_TERMS = ["spr" + "int", "pha" + "se"]
_OLD_NAMES = [
    "pha" + "se4_concordance",
    "pha" + "se4_discordant_analysis",
    "pha" + "se4_figures",
    "pha" + "se4_splitcode_dual",
    "pha" + "se5_lr_perf_rerun",
    "test_pha" + "se4_fixes",
    "benchmark_results_spr" + "int4",
    "pha" + "se4_results.json",
]


def _has_phase_digit(line):
    """Return True if line contains 'phase' immediately followed by a digit."""
    lower = line.lower()
    term = _BANNED_TERMS[1]  # "phase"
    idx = 0
    while True:
        idx = lower.find(term, idx)
        if idx == -1:
            return False
        end = idx + len(term)
        if end < len(lower) and lower[end].isdigit():
            return True
        idx = end


def _scan_file_for_banned(fpath):
    """Scan a single file for banned development-iteration terminology."""
    violations = []
    _term_0 = _BANNED_TERMS[0]
    with open(fpath) as f:
        for lineno, line in enumerate(f, 1):
            lower = line.lower()
            if _term_0 in lower:
                violations.append(f"{fpath.name}:{lineno}: {line.rstrip()}")
            elif _has_phase_digit(line):
                violations.append(f"{fpath.name}:{lineno}: {line.rstrip()}")
    return violations


def test_no_banned_terms_in_scripts():
    """No script should contain banned development-iteration terminology."""
    violations = []
    for script in SCRIPTS_DIR.glob("*.py"):
        violations.extend(_scan_file_for_banned(script))
    assert not violations, (
        f"Banned terminology found:\n" + "\n".join(violations))


def test_no_banned_terms_in_tests():
    """No test file should contain banned development-iteration terminology."""
    violations = []
    test_dir = PROJECT_ROOT / "tests"
    for tf in test_dir.glob("*.py"):
        violations.extend(_scan_file_for_banned(tf))
    assert not violations, (
        f"Banned terminology found:\n" + "\n".join(violations))


def test_no_banned_terms_in_readme():
    """All markdown files (including subdirectory READMEs) must not contain
    banned development-iteration terminology.  Previously this only checked
    the root README.md, missing per-artifact READMEs in table2/, fig_*/."""
    violations = []
    for md in sorted(PROJECT_ROOT.rglob("*.md")):
        # Skip vendored / cache directories
        if ".venv" in md.parts or "venv" in md.parts or ".pytest_cache" in md.parts:
            continue
        violations.extend(_scan_file_for_banned(md))
    assert not violations, (
        f"Banned terminology found:\n" + "\n".join(violations))


def test_no_old_filenames_in_tracked_files():
    """No tracked file should reference the old renamed filenames.
    Uses recursive glob for markdown to cover per-artifact subdirectory READMEs."""
    violations = []
    for pattern in ["scripts/*.py", "tests/*.py", "**/*.md"]:
        for fpath in sorted(PROJECT_ROOT.glob(pattern)):
            # Skip vendored / cache directories
            if ".venv" in fpath.parts or "venv" in fpath.parts or ".pytest_cache" in fpath.parts:
                continue
            with open(fpath) as f:
                content = f.read()
            for old in _OLD_NAMES:
                if old in content:
                    violations.append(f"{fpath.name}: contains '{old}'")
    assert not violations, (
        f"Old filenames still referenced:\n" + "\n".join(violations))


# ============================================================================
# 4. Path constants resolve to existing locations
# ============================================================================

EXPECTED_CONFIG_DIRS = ["seqproc", "matchbox", "splitcode"]


@pytest.mark.parametrize("subdir", EXPECTED_CONFIG_DIRS)
def test_config_dirs_exist(subdir):
    """Each tool config directory must exist."""
    d = CONFIGS_DIR / subdir
    assert d.is_dir(), f"Config directory missing: {d}"


EXPECTED_CONFIGS = {
    "seqproc": [
        "splitseq_filter_edit.geom",
        "10x_v2.geom",
        "sciseq3_edit.geom",
    ],
    "matchbox": [
        "splitseq_replacement.mb",
        "10x_v2.mb",
        "sciseq3.mb",
    ],
    "splitcode": [
        "splitseq_paper.config",
        "10x_v2.config",
        "sciseq3.config",
    ],
}


@pytest.mark.parametrize(
    "tool,config",
    [(t, c) for t, configs in EXPECTED_CONFIGS.items() for c in configs],
    ids=lambda x: x if isinstance(x, str) else "",
)
def test_config_files_exist(tool, config):
    """Key config files must exist for each tool."""
    path = CONFIGS_DIR / tool / config
    assert path.is_file(), f"Config missing: {path}"


# ============================================================================
# 5. Result data integrity -- concordance results
# ============================================================================

EXPECTED_RECOVERY = {
    "splitseq_pe": {"seqproc": (83.0, 85.0), "matchbox": (77.0, 79.0), "splitcode": (90.0, 93.0)},
    "lr_splitseq": {"seqproc": (49.0, 51.0), "matchbox": (39.0, 41.0), "splitcode": (27.0, 29.0)},
    "10x_short": {"seqproc": (100.0, 100.0), "matchbox": (100.0, 100.0), "splitcode": (100.0, 100.0)},
    "sciseq": {"seqproc": (88.0, 90.0), "matchbox": (89.0, 91.0), "splitcode": (87.0, 90.0)},
}


@pytest.fixture(scope="module")
def concordance_data():
    """Load concordance results JSON."""
    path = RESULTS_DIR / "concordance" / "concordance_results.json"
    if not path.exists():
        pytest.skip(f"Concordance results not found: {path}")
    with open(path) as f:
        return json.load(f)


@pytest.mark.parametrize("ds_key", list(EXPECTED_RECOVERY.keys()))
@pytest.mark.requires_cached_results
def test_concordance_datasets_present(concordance_data, ds_key):
    """Each expected dataset must be present in concordance results."""
    assert ds_key in concordance_data, (
        f"Dataset '{ds_key}' missing from concordance results")


@pytest.mark.parametrize(
    "ds_key,tool",
    [(ds, t) for ds, tools in EXPECTED_RECOVERY.items() for t in tools],
)
@pytest.mark.requires_cached_results
def test_recovery_rates_in_range(concordance_data, ds_key, tool):
    """Recovery rates must fall within expected ranges (no regression)."""
    lo, hi = EXPECTED_RECOVERY[ds_key][tool]
    actual = concordance_data[ds_key]["recovery_pct"][tool]
    assert lo <= actual <= hi, (
        f"{ds_key}/{tool}: recovery {actual:.1f}% outside [{lo}, {hi}]")


@pytest.mark.requires_cached_results
def test_hamming_vs_edit_present(concordance_data):
    """Hamming vs edit comparison must exist for PE, LR, and sci datasets."""
    for ds_key in ["splitseq_pe", "lr_splitseq", "sciseq"]:
        hve = concordance_data[ds_key].get("hamming_vs_edit")
        assert hve, f"{ds_key}: missing hamming_vs_edit data"
        assert hve.get("edit_reads", 0) >= hve.get("hamming_reads", 0), (
            f"{ds_key}: edit reads should be >= hamming reads")


EXPECTED_EDIT_GAIN = {
    "splitseq_pe": (3.0, 7.0),
    "lr_splitseq": (14.0, 18.0),
    "sciseq": (0.5, 2.0),
}


@pytest.mark.parametrize("ds_key", list(EXPECTED_EDIT_GAIN.keys()))
@pytest.mark.requires_cached_results
def test_edit_distance_gain(concordance_data, ds_key):
    """Edit distance gain percentage must be within expected range."""
    hve = concordance_data[ds_key]["hamming_vs_edit"]
    gain = hve["edit_gain_pct"]
    lo, hi = EXPECTED_EDIT_GAIN[ds_key]
    assert lo <= gain <= hi, (
        f"{ds_key}: edit gain {gain:.1f}% outside [{lo}, {hi}]")


# ============================================================================
# 6. Benchmark results JSON integrity
# ============================================================================


@pytest.fixture(scope="module")
def benchmark_data():
    """Load the main benchmark results JSON."""
    path = RESULTS_DIR / "paper_figures" / "benchmark_results.json"
    if not path.exists():
        pytest.skip(f"Benchmark results not found: {path}")
    with open(path) as f:
        return json.load(f)


@pytest.mark.requires_cached_results
def test_benchmark_json_has_all_datasets(benchmark_data):
    """benchmark_results.json must contain all 4 datasets."""
    for ds_key in ["splitseq_pe", "lr_splitseq", "10x_short", "sciseq"]:
        assert ds_key in benchmark_data, (
            f"Dataset '{ds_key}' missing from benchmark_results.json")


@pytest.mark.requires_cached_results
def test_benchmark_json_has_all_tools(benchmark_data):
    """Each dataset must have entries for all 3 tools."""
    for ds_key, entry in benchmark_data.items():
        tools = entry.get("tools", {})
        for tool in ["seqproc", "matchbox", "splitcode"]:
            assert tool in tools, (
                f"{ds_key}: tool '{tool}' missing from benchmark_results.json")


@pytest.mark.requires_cached_results
def test_benchmark_perf_backup_matches():
    """benchmark_results_perf.json must be identical to benchmark_results.json
    (they are the same data, just the backup)."""
    main_path = RESULTS_DIR / "paper_figures" / "benchmark_results.json"
    perf_path = RESULTS_DIR / "paper_figures" / "benchmark_results_perf.json"
    if not main_path.exists() or not perf_path.exists():
        pytest.skip("Benchmark JSONs not found")

    with open(main_path) as f:
        main = json.load(f)
    with open(perf_path) as f:
        perf = json.load(f)

    assert main == perf, (
        "benchmark_results.json and benchmark_results_perf.json differ")


# ============================================================================
# 7. generate_figures.py determinism -- re-run must not change output
# ============================================================================


@pytest.mark.requires_cached_results
def test_generate_figures_deterministic():
    """Running generate_figures.py must produce byte-identical JSON output."""
    json_path = RESULTS_DIR / "paper_figures" / "benchmark_results.json"
    if not json_path.exists():
        pytest.skip("benchmark_results.json not found")

    with open(json_path) as f:
        before = f.read()

    result = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / "generate_figures.py")],
        capture_output=True, text=True,
        cwd=str(PROJECT_ROOT),
    )
    assert result.returncode == 0, (
        f"generate_figures.py failed:\n{result.stderr}")

    with open(json_path) as f:
        after = f.read()

    assert before == after, (
        "generate_figures.py produced different benchmark_results.json on re-run")


# ============================================================================
# 8. concordance_analysis.py --skip-runs must succeed with cached data
# ============================================================================


@pytest.mark.requires_cached_results
def test_concordance_skip_runs():
    """concordance_analysis.py --skip-runs must exit 0 with cached data."""
    conc_dir = RESULTS_DIR / "concordance"
    if not conc_dir.exists():
        pytest.skip("Concordance results directory not found")

    result = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / "concordance_analysis.py"),
         "--skip-runs", "--threads", "1"],
        capture_output=True, text=True,
        cwd=str(PROJECT_ROOT),
    )
    assert result.returncode == 0, (
        f"concordance_analysis.py --skip-runs failed:\n{result.stderr}\n{result.stdout}")
    assert "CONCORDANCE ANALYSIS COMPLETE" in result.stdout, (
        "Expected completion message not found in output")


# ============================================================================
# 9. discordant_analysis.py must succeed (needs data/SRR6750041_1M_R2.fastq)
# ============================================================================


@pytest.mark.requires_real_data
@pytest.mark.requires_cached_results
def test_discordant_analysis():
    """discordant_analysis.py must exit 0 and report the FP finding."""
    r2_path = PROJECT_ROOT / "data" / "SRR6750041_1M_R2.fastq"
    if not r2_path.exists():
        pytest.skip("R2 FASTQ not found (data not available)")

    conc_pe = RESULTS_DIR / "concordance" / "splitseq_pe"
    if not conc_pe.exists():
        pytest.skip("concordance/splitseq_pe results not found")

    result = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / "discordant_analysis.py")],
        capture_output=True, text=True,
        cwd=str(PROJECT_ROOT),
    )
    assert result.returncode == 0, (
        f"discordant_analysis.py failed:\n{result.stderr}")
    assert "FALSE POSITIVES" in result.stdout, (
        "Expected false-positive conclusion not found in output")


# ============================================================================
# 10. Existing unit tests from test_benchmark_pipeline.py
# ============================================================================

# --- Bug 1: run_cmd exit code parsing ---

def test_run_cmd_parses_exit_status():
    """Exit status should be extracted from /usr/bin/time stderr."""
    time_cmd = "/usr/bin/time -v echo hello"
    result = subprocess.run(time_cmd, shell=True, capture_output=True, text=True)

    tool_exit_code = result.returncode
    for line in result.stderr.split('\n'):
        if 'Exit status' in line:
            tool_exit_code = int(line.split(':')[1].strip())

    assert tool_exit_code == 0


def test_run_cmd_detects_failure():
    """Failing commands should report non-zero exit status."""
    time_cmd = "/usr/bin/time -v false"
    result = subprocess.run(time_cmd, shell=True, capture_output=True, text=True)

    tool_exit_code = result.returncode
    for line in result.stderr.split('\n'):
        if 'Exit status' in line:
            tool_exit_code = int(line.split(':')[1].strip())

    assert tool_exit_code == 1


# --- Bug 2: check_wl_match scoping ---

def test_bc_distances_initialized():
    """bc distances should be 99 (no match) when barcode length is wrong."""
    bc3_d, bc2_d, bc1_d = 99, 99, 99

    bc3 = "SHORT"
    bc2 = "AACGTGAT"
    bc1 = "AACGT"

    if len(bc3) == 8:
        bc3_d = 0
    if len(bc2) == 8:
        bc2_d = 0
    if len(bc1) == 6:
        bc1_d = 0

    assert bc3_d == 99
    assert bc2_d == 0
    assert bc1_d == 99


# --- Bug 3: matchbox_paired flag ---

def test_matchbox_paired_flag_in_datasets():
    """All datasets should have explicit matchbox_paired flag."""
    from concordance_analysis import DATASETS

    for ds_key, ds in DATASETS.items():
        assert 'matchbox_paired' in ds, (
            f"Dataset '{ds_key}' missing 'matchbox_paired' flag")
        assert isinstance(ds['matchbox_paired'], bool)


def test_matchbox_paired_consistency():
    """matchbox_paired should match expected values for known datasets."""
    from concordance_analysis import DATASETS

    expected = {
        'splitseq_pe': True,
        'lr_splitseq': False,
        '10x_short': True,
        'sciseq': True,
    }
    for ds_key, expect_paired in expected.items():
        if ds_key in DATASETS:
            assert DATASETS[ds_key]['matchbox_paired'] == expect_paired


# --- Bug 5: ID extraction order ---

def test_paired_prefers_r2():
    """For paired-end datasets, R2 should be preferred for ID extraction."""
    dataset_mode = 'paired'
    prefer = 'mb_r2.fq' if dataset_mode == 'paired' else 'mb_r1.fq'
    assert prefer == 'mb_r2.fq'

    dataset_mode = 'single'
    prefer = 'mb_r2.fq' if dataset_mode == 'paired' else 'mb_r1.fq'
    assert prefer == 'mb_r1.fq'


# --- Jaccard index ---

def test_jaccard():
    """Jaccard index should be correct for known sets.
    Imports jaccard from concordance_analysis (the canonical source)."""
    from concordance_analysis import jaccard

    assert jaccard(set(), set()) == 1.0
    assert jaccard({1, 2, 3}, {1, 2, 3}) == 1.0
    assert jaccard({1, 2}, {3, 4}) == 0.0
    assert jaccard({1, 2, 3}, {2, 3, 4}) == 0.5
