#!/usr/bin/env python3
"""Replace one tool in a provenance-rich publication accuracy artifact.

This is used when a configuration is rerun without repeating unaffected tools.
It replaces the metric and run provenance for one dataset/tool and recomputes
all pairwise read-set intersections from the recorded accession bitmaps.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from generate_emitted_set_upset import canonical_bitmap, tool_bitmap


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def popcount(value: int) -> int:
    """Count set bits on both current and older cluster Python versions."""
    return value.bit_count() if hasattr(value, "bit_count") else bin(value).count("1")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-artifact", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--tool", required=True)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--metrics-key", required=True)
    parser.add_argument("--run-json", type=Path, required=True)
    parser.add_argument("--bitmap", type=Path, required=True)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    return parser.parse_args()


def publication_metric(dataset: str, tool: str, metric: dict) -> dict:
    return {
        "dataset": dataset,
        "emitted_fraction": metric["emitted_fraction"],
        "emitted_records": metric["emitted_records"],
        "f1": metric["f1"],
        "input_records": metric["input_records"],
        "intersection_records": metric["intersection_records"],
        "precision": metric["precision"],
        "recall": metric["recall"],
        "reference_records": metric["reference"]["records"],
        "tool": tool,
    }


def run_provenance(
    dataset: str,
    tool: str,
    attempt: Path,
    run: dict,
    bitmap: Path,
    metric: dict,
) -> dict:
    identity = run["identity"]
    counts_by_path = {Path(item["path"]): item for item in run["output_counts"]}
    bitmap_digest = sha256(bitmap)
    sources = []
    for output in run["outputs"]:
        output_path = Path(output["path"])
        audit = counts_by_path[output_path]
        sources.append(
            {
                "audit": {
                    "max_sequence_length": audit["max_sequence_length"],
                    "min_sequence_length": audit["min_sequence_length"],
                    "records": audit["records"],
                    "sequence_length_counts": audit["sequence_length_counts"],
                },
                "bitmap": str(bitmap.resolve()),
                "bitmap_sha256": bitmap_digest,
                "fastq": str(output_path),
                "fastq_sha256": output["sha256"],
            }
        )

    run_fields = {
        key: run[key]
        for key in (
            "peak_rss_kib",
            "system_cpu_seconds",
            "user_cpu_seconds",
            "wall_seconds",
        )
    }
    metadata = identity["metadata"]
    if metadata["dataset"] != dataset or metadata["tool"] != tool:
        raise ValueError("run metadata does not match requested dataset/tool")
    return {
        "attempt": str(attempt.resolve()),
        "condition_id": metadata["condition_id"],
        "configs": identity["configs"],
        "executables": identity["executables"],
        "group_records": {"accepted": metric["emitted_records"]},
        "repositories": identity["repositories"],
        "run": run_fields,
        "sources": sources,
    }


def main() -> None:
    args = parse_args()
    artifact = json.loads(args.base_artifact.read_text())
    metrics_artifact = json.loads(args.metrics.read_text())
    metric = metrics_artifact["tools"][args.metrics_key]
    metric = {
        **metric,
        "input_records": metrics_artifact["input_records"],
        "reference": metrics_artifact["reference"],
    }
    run = json.loads(args.run_json.read_text())
    attempt = args.run_json.resolve().parent
    _, replacement_bitmap = canonical_bitmap(args.bitmap)

    replacements = [
        index
        for index, item in enumerate(artifact["metrics"])
        if item["dataset"] == args.dataset and item["tool"] == args.tool
    ]
    if len(replacements) != 1:
        raise ValueError(f"expected one metric row to replace, found {len(replacements)}")
    artifact["metrics"][replacements[0]] = publication_metric(
        args.dataset, args.tool, metric
    )
    artifact["provenance"][f"{args.dataset}/{args.tool}"] = run_provenance(
        args.dataset, args.tool, attempt, run, args.bitmap, metric
    )

    # Recompute every pair using the new bitmap for the replaced tool and the
    # checksummed source bitmaps already recorded for unaffected tools.
    bitmap_cache: dict[str, int] = {args.tool: replacement_bitmap}
    for pair in artifact["pairwise"]:
        left = pair["left"]
        right = pair["right"]
        for name in (left, right):
            if name not in bitmap_cache:
                bitmap_cache[name] = tool_bitmap(
                    artifact, args.dataset, name, []
                )[0]
        left_bitmap = bitmap_cache[left]
        right_bitmap = bitmap_cache[right]
        intersection = popcount(left_bitmap & right_bitmap)
        left_count = popcount(left_bitmap)
        right_count = popcount(right_bitmap)
        union = left_count + right_count - intersection
        pair.update(
            {
                "intersection_records": intersection,
                "jaccard": intersection / union,
                "left_only_records": left_count - intersection,
                "right_only_records": right_count - intersection,
                "union_records": union,
            }
        )

    artifact["campaigns"] = [artifact.get("campaign"), str(args.campaign.resolve())]
    artifact["replacement"] = {
        "dataset": args.dataset,
        "tool": args.tool,
        "campaign": str(args.campaign.resolve()),
        "manifest": {
            "path": str(args.manifest.resolve()),
            "sha256": sha256(args.manifest),
        },
        "metrics": {
            "path": str(args.metrics.resolve()),
            "sha256": sha256(args.metrics),
            "key": args.metrics_key,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
