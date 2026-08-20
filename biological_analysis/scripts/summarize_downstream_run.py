#!/usr/bin/env python3
"""Create one machine-readable and one Markdown summary from a completed run."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


TOOLS = ("seqproc", "splitcode", "matchbox")
PREFIX = {"seqproc": "sp", "splitcode": "sc", "matchbox": "mb"}


def read_json(path: Path):
    return json.loads(path.read_text())


def star_summary(path: Path) -> dict[str, object]:
    result = {}
    for key, value in csv.reader(path.open()):
        try:
            result[key] = int(value)
        except ValueError:
            try:
                result[key] = float(value)
            except ValueError:
                result[key] = value
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run", type=Path)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    run = args.run.resolve()
    biological = read_json(run / "analysis/biological_metrics.json")
    counts = read_json(run / "analysis/count_concordance.json")
    read_sets = read_json(run / "analysis/read_set_jaccard.json")
    provenance = read_json(run / "downstream_provenance.json")
    fastq = {item["tool"]: item for item in provenance["fastqs"]}
    resources = list(csv.DictReader((run / "resources.csv").open()))
    resource_by_tool = {
        row["tool"]: {
            "seconds": float(row["seconds"]),
            "peak_ram_mb": float(row["peak_ram_mb"]),
        }
        for row in resources
        if row["step"] == "starsolo"
    }
    summary = {
        "schema_version": 1,
        "run": str(run),
        "parameters": provenance["parameters"],
        "reference": provenance["genome_index"]["manifest"]["reference"],
        "software": {"STAR": provenance["software"]["STAR"]},
        "tools": {},
        "pairs": {},
        "shared_cells": biological["shared_cells"],
        "all_tool_celltype_agreement": biological["celltype_agreement_shared"],
        "celltype_jaccard_mean": biological["celltype_jaccard_mean"],
        "celltype_jaccard_pairwise_mean": biological["celltype_jaccard_pairwise_mean"],
        "celltype_jaccard_pairwise_per_type": biological["celltype_jaccard_pairwise_per_type"],
        "cluster_ari_mean": biological["cluster_ari_mean"],
        "joint_coclustering_agreement": biological["joint_coclustering_agreement"],
        "tool_mixing_entropy": biological["tool_mixing_entropy"],
        "celltype_jaccard_per_type": biological["celltype_jaccard_per_type"],
    }
    for tool in TOOLS:
        summary["tools"][tool] = {
            "input_pairs": fastq[tool]["records"],
            "called_cells_min_umi": biological["cells"][tool],
            "barcode_rank": counts["barcode_rank_inflection"][tool],
            "STARsolo": star_summary(run / f"{PREFIX[tool]}_Solo.out/Gene/Summary.csv"),
            "STARsolo_resources": resource_by_tool[tool],
        }
    for pair, item in read_sets["pairs"].items():
        summary["pairs"][pair] = {
            **item,
            "per_gene_pearson_log1p": counts["per_gene_total_pearson_logspace"][pair],
            "per_gene_spearman": counts["per_gene_total_spearman"][pair],
            "per_barcode_pearson_log1p": counts["per_barcode_umi_pearson_logspace"][pair],
            "per_barcode_spearman": counts["per_barcode_umi_spearman"][pair],
            "celltype_agreement": biological["celltype_agreement_pairwise"][pair],
            "celltype_jaccard_mean": biological["celltype_jaccard_pairwise_mean"][pair],
            "cluster_ari": biological["cluster_ari_pairwise"][pair],
        }
    args.json.write_text(json.dumps(summary, indent=2) + "\n")

    lines = [
        "# Final SPLiT-seq PE downstream summary",
        "",
        "| Tool | Input pairs | Valid barcodes | Called cells | Inflection rank | STARsolo time (s) | Peak RSS (MiB) |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for tool in TOOLS:
        item = summary["tools"][tool]
        lines.append(
            f"| {tool} | {item['input_pairs']:,} | "
            f"{item['STARsolo']['Reads With Valid Barcodes']:.4%} | "
            f"{item['called_cells_min_umi']:,} | {item['barcode_rank']['infl_rank']:,} | "
            f"{item['STARsolo_resources']['seconds']:.2f} | "
            f"{item['STARsolo_resources']['peak_ram_mb']:.1f} |"
        )
    lines.extend(
        [
            "",
            "| Tool pair | Read-set Jaccard | Per-gene Pearson | Per-barcode Pearson | Cell-type agreement | Mean type Jaccard | Cluster ARI |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for pair, item in summary["pairs"].items():
        lines.append(
            f"| {pair.replace('|', ' / ')} | {item['jaccard']:.4f} | "
            f"{item['per_gene_pearson_log1p']:.4f} | "
            f"{item['per_barcode_pearson_log1p']:.4f} | "
            f"{item['celltype_agreement']:.4f} | {item['celltype_jaccard_mean']:.4f} | "
            f"{item['cluster_ari']:.4f} |"
        )
    lines.extend(
        [
            "",
            f"Shared called cells: **{summary['shared_cells']:,}**.  "
            f"All-tool cell-type agreement: **{summary['all_tool_celltype_agreement']:.4f}**.  "
            f"Mean per-type Jaccard: **{summary['celltype_jaccard_mean']:.4f}**.",
        ]
    )
    args.markdown.write_text("\n".join(lines) + "\n")
    print(args.markdown.read_text())


if __name__ == "__main__":
    main()
