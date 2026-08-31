import json
import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from summarize_hot_path import summarize


def write_run(root: Path, run_id: str, statistics: bool, seconds, level=None):
    root.mkdir()
    (root / "run.json").write_text(
        json.dumps({"run_id": run_id, "attempt": 1, "success": True})
    )
    payload = {
        "mode": "seeded",
        "statistics": statistics,
        "reads": 100,
        "threads": 4,
        "seconds": seconds,
    }
    if level is not None:
        payload["statistics_level"] = level
    (root / "stdout.txt").write_text(json.dumps(payload))
    return root


def test_summarize_computes_machine_derived_effect(tmp_path):
    off = write_run(tmp_path / "off", "off-id", False, [1.0, 1.0, 1.0])
    on = write_run(tmp_path / "on", "on-id", True, [1.1, 1.1, 1.1])

    result = summarize([off, on])

    assert result["conditions"][0]["n"] == 3
    effect = result["effects"][0]
    assert round(effect["mean_time_reduction_pct"], 6) == round(100 * (0.1 / 1.1), 6)
    assert round(effect["throughput_gain_from_mean_pct"], 6) == 10.0
    assert effect["statistics_off_run_id"] == "off-id"
    assert effect["statistics_on_run_id"] == "on-id"


def test_summarize_compares_basic_and_detailed_to_the_same_off_run(tmp_path):
    off = write_run(tmp_path / "off", "off-id", False, [1.0], "off")
    basic = write_run(tmp_path / "basic", "basic-id", True, [1.01], "basic")
    detailed = write_run(
        tmp_path / "detailed", "detailed-id", True, [1.02], "detailed"
    )

    result = summarize([off, basic, detailed])

    assert [effect["statistics_level"] for effect in result["effects"]] == [
        "basic",
        "detailed",
    ]
    assert round(result["effects"][0]["mean_time_overhead_pct"], 6) == 1.0
    assert round(result["effects"][1]["mean_time_overhead_pct"], 6) == 2.0
