#!/usr/bin/env python3
"""Derive a frozen Matchbox SPLiT-seq PE sensitivity manifest.

The input is a standard publication correctness manifest.  This tool selects
its single 32-thread Matchbox SPLiT-seq PE condition and substitutes an
explicit sensitivity configuration and its support files.  Keeping this
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
        == "correctness"
        and int(run.get("spec", {}).get("metadata", {}).get("threads", -1)) == 32
    ]
    if len(candidates) != 1:
        raise ValueError(
            f"expected one 32-thread Matchbox SPLiT-seq PE correctness run; found {len(candidates)}"
        )

    run = copy.deepcopy(candidates[0])
    spec = run["spec"]
    metadata = spec["metadata"]
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

    condition_id = "splitseq_pe-t32-r1-matchbox-ham1-expanded-fuzzy-linkers"
    block_id = "splitseq_pe-ham1-expanded-t32-r1"
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
                "linkers, and explicit component-length guards"
            ),
            "tier": "publication-sensitivity-full-data-correctness",
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
            "name": "matchbox-splitseq-pe-hamming1-expanded-sensitivity",
            "random_seed": args.seed,
            "purpose": (
                "single materialized 32-thread sensitivity run using externally "
                "Hamming-distance-one-expanded barcodes, fuzzy linkers, and "
                "component-length guards"
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
