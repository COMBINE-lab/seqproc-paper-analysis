#!/usr/bin/env python3
"""Derive a frozen Matchbox SPLiT-seq PE sensitivity manifest.

The input is a standard publication correctness or timing manifest. This tool
selects one 32-thread Matchbox SPLiT-seq PE condition and substitutes an
explicit sensitivity configuration and its support files. Keeping this
transformation separate prevents externally expanded barcode resources from
being mistaken for the primary, canonical-list workload.
"""

from __future__ import annotations

import argparse
import copy
import json
import subprocess
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import yaml

from benchmark_harness import sha256_file


def git_commit(repository: Path) -> str:
    return subprocess.run(
        ["git", "-C", str(repository), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def frozen(path: Path) -> dict[str, Any]:
    resolved = path.resolve()
    if not resolved.is_file():
        raise FileNotFoundError(resolved)
    return {"path": str(resolved), "sha256": sha256_file(resolved), "verify": False}


def build_manifest(args: argparse.Namespace) -> dict[str, Any]:
    base_path = args.base_manifest.resolve()
    base = yaml.safe_load(base_path.read_text())
    if not isinstance(base, dict) or base.get("mode") != "publication":
        raise ValueError("base manifest must be a publication manifest")

    candidates = [
        run
        for run in base.get("runs", [])
        if run.get("spec", {}).get("metadata", {}).get("dataset") == "splitseq_pe"
        and run.get("spec", {}).get("metadata", {}).get("tool") == "matchbox"
        and run.get("spec", {}).get("metadata", {}).get("measurement_track")
        == args.measurement_track
        and int(run.get("spec", {}).get("metadata", {}).get("threads", -1)) == 32
        and (
            args.measurement_track != "timing"
            or int(run.get("spec", {}).get("metadata", {}).get("replicate", -1))
            == args.replicate
        )
    ]
    if len(candidates) != 1:
        raise ValueError(
            "expected one 32-thread Matchbox SPLiT-seq PE "
            f"{args.measurement_track} run; found {len(candidates)}"
        )

    run = copy.deepcopy(candidates[0])
    spec = run["spec"]
    metadata = spec["metadata"]
    selected_replicate = int(metadata.get("replicate", 1))
    sensitivity_config = args.config.resolve()
    support = [path.resolve() for path in args.support]

    command = list(spec["command"])
    try:
        script_index = command.index("--script-file") + 1
    except ValueError as error:
        raise ValueError("base Matchbox command has no --script-file") from error
    command[script_index] = str(sensitivity_config)
    spec["command"] = command
    spec["configs"] = [frozen(sensitivity_config), *(frozen(path) for path in support)]

    condition_id = (
        f"splitseq_pe-t32-r{selected_replicate}-matchbox-"
        f"ham1-expanded-fuzzy-linkers-{args.measurement_track}"
    )
    block_id = (
        f"splitseq_pe-ham1-expanded-{args.measurement_track}-"
        f"t32-r{selected_replicate}"
    )
    run["id"] = condition_id
    run["block"] = block_id
    metadata.update(
        {
            "condition_id": condition_id,
            "block_id": block_id,
            "analysis_role": "sensitivity",
            "execution_mode": "external-hamming1-expanded-fuzzy-linkers",
            "equivalence_scope": (
                "sensitivity-only; exact matching against an externally generated "
                "Hamming-distance-one-expanded barcode list, edit-distance-three "
                "linkers, fixed component lengths, a terminal exact-barcode anchor, "
                "and post-placement membership checks for the two upstream windows"
            ),
            "tier": f"publication-sensitivity-full-data-{args.measurement_track}",
            "external_resource_preprocessing": "Hamming-distance-one barcode expansion",
        }
    )

    analysis_repo = args.analysis_repo.resolve()
    current_analysis_commit = git_commit(analysis_repo)
    for repository in spec.get("repositories", []):
        if repository.get("name") == "analysis":
            repository["commit"] = current_analysis_commit

    used_artifacts: dict[str, dict[str, Any]] = {}
    for collection in ("executables", "inputs", "configs"):
        for entry in spec.get(collection, []):
            record = dict(entry)
            used_artifacts[str(Path(record["path"]).resolve())] = {
                "path": str(Path(record["path"]).resolve()),
                "sha256": str(record["sha256"]),
            }
    # The harness records and validates argv[0] independently of the explicitly
    # declared executable collection.  Here that is the taskset affinity wrapper.
    command_executable = frozen(Path(spec["command"][0]))
    used_artifacts[command_executable["path"]] = {
        "path": command_executable["path"],
        "sha256": command_executable["sha256"],
    }
    for path in (Path(__file__).resolve(), base_path):
        record = frozen(path)
        used_artifacts[record["path"]] = {
            "path": record["path"],
            "sha256": record["sha256"],
        }

    datasets = [
        copy.deepcopy(dataset)
        for dataset in base.get("datasets", [])
        if dataset.get("name") == "splitseq_pe"
    ]
    if len(datasets) != 1:
        raise ValueError("base manifest must describe SPLiT-seq PE exactly once")

    return {
        "schema_version": "1.0.0",
        "mode": "publication",
        "study": {
            "name": (
                "matchbox-splitseq-pe-hamming1-expanded-sensitivity-"
                f"{args.measurement_track}"
            ),
            "random_seed": args.seed,
            "purpose": (
                f"single 32-thread {args.measurement_track} sensitivity run using externally "
                "Hamming-distance-one-expanded barcodes, fuzzy linkers, a terminal "
                "exact-barcode anchor, and fixed-window membership checks"
            ),
            "require_cross_tool_identity": False,
        },
        "artifacts": sorted(used_artifacts.values(), key=lambda item: item["path"]),
        "source_manifest": frozen(base_path),
        "datasets": datasets,
        "execution": copy.deepcopy(base.get("execution", {})),
        "runs": [run],
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-manifest", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--support", type=Path, nargs="+", required=True)
    parser.add_argument(
        "--measurement-track",
        choices=("correctness", "timing"),
        default="correctness",
    )
    parser.add_argument(
        "--replicate",
        type=int,
        default=1,
        help="source timing replicate to select; ignored for correctness",
    )
    parser.add_argument("--analysis-repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--seed", type=int, default=8282028)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.output.exists():
        parser.error(f"refusing to overwrite {args.output}")
    manifest = build_manifest(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(yaml.safe_dump(manifest, sort_keys=False))
    print(json.dumps({"conditions": len(manifest["runs"]), "output": str(args.output)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
