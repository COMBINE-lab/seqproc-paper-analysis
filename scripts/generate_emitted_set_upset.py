#!/usr/bin/env python3
"""Generate the manuscript emitted-read UpSet plot from canonical bitmaps.

The figure replaces the percentage table formerly labelled ``tab:concordance``.
It reads the provenance-rich publication accuracy artifacts, verifies the
recorded bitmap checksums, unions multi-product/multi-orientation sources for
each tool, computes the seven mutually exclusive three-tool intersections, and
writes a dependency-free SVG plus machine-readable CSV/JSON.  When
``rsvg-convert`` is available, publication-ready PDF and PNG companions are
also produced.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
from pathlib import Path
from xml.sax.saxutils import escape


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RESULTS = ROOT / "publication_results/journal_rerun_2026-08-17"
DEFAULT_OUTPUT = DEFAULT_RESULTS / "fig_emitted_set_upset"
MAGIC = b"fastq_numeric_accession_set_v1"

TOOLS = ("seqproc", "splitcode", "matchbox")
DATASETS = (
    ("splitseq_pe", "SPLiT-seq PE", "#0072B2"),
    ("lr_splitseq_dual", "LR-SPLiT-seq (dual)", "#D55E00"),
    ("tenx_v2", "10x Chromium v2", "#009E73"),
    ("scirnaseq3", "sci-RNA-seq3", "#CC79A7"),
)

# Membership bits follow TOOLS: seqproc=1, splitcode=2, matchbox=4.
INTERSECTIONS = (
    (7, "All three"),
    (3, "seqproc + splitcode"),
    (5, "seqproc + matchbox"),
    (6, "splitcode + matchbox"),
    (1, "seqproc only"),
    (2, "splitcode only"),
    (4, "matchbox only"),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_bitmap(path: Path) -> tuple[dict[str, object], int]:
    fields = path.read_bytes().split(b"\0", 4)
    if len(fields) != 5 or fields[0] != MAGIC:
        raise RuntimeError(f"unexpected accession bitmap format: {path}")
    metadata = {
        "mate": int(fields[1]),
        "numeric_id_max": int(fields[2]),
        "accession_prefix": fields[3].decode(),
    }
    return metadata, int.from_bytes(fields[4], "little")


def resolve_bitmap(path_text: str, search_roots: list[Path]) -> Path:
    path = Path(path_text)
    if path.is_file():
        return path
    candidates = [root / path.name for root in search_roots if (root / path.name).is_file()]
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise FileNotFoundError(f"bitmap not found: {path}")
    raise RuntimeError(f"ambiguous bitmap basename {path.name}: {candidates}")


def tool_bitmap(
    artifact: dict,
    dataset: str,
    tool: str,
    search_roots: list[Path],
) -> tuple[int, dict[str, object], list[dict[str, object]]]:
    key = f"{dataset}/{tool}"
    provenance = artifact.get("provenance", {}).get(key)
    if provenance is None:
        raise RuntimeError(f"missing provenance for {key}")
    sources = provenance.get("sources", [])
    if not sources:
        raise RuntimeError(f"no bitmap sources for {key}")

    combined = 0
    common_metadata = None
    checked_sources = []
    for source in sources:
        path = resolve_bitmap(str(source["bitmap"]), search_roots)
        observed_sha256 = sha256(path)
        expected_sha256 = source.get("bitmap_sha256")
        if expected_sha256 and observed_sha256 != expected_sha256:
            raise RuntimeError(
                f"bitmap checksum mismatch for {path}: "
                f"expected {expected_sha256}, observed {observed_sha256}"
            )
        metadata, bitmap = canonical_bitmap(path)
        if common_metadata is None:
            common_metadata = metadata
        elif metadata["numeric_id_max"] != common_metadata["numeric_id_max"]:
            raise RuntimeError(f"numeric accession domain differs within {key}")
        elif metadata["accession_prefix"] != common_metadata["accession_prefix"]:
            raise RuntimeError(f"accession prefix differs within {key}")
        combined |= bitmap
        checked_sources.append(
            {"path": str(path), "sha256": observed_sha256, "mate": metadata["mate"]}
        )

    assert common_metadata is not None
    return combined, common_metadata, checked_sources


def exclusive_counts(bitmaps: tuple[int, int, int], universe_size: int) -> dict[int, int]:
    if universe_size <= 0:
        raise ValueError("universe_size must be positive")
    universe = (1 << universe_size) - 1
    counts = {}
    for membership in range(1, 8):
        selected = universe
        for tool_index, bitmap in enumerate(bitmaps):
            if membership & (1 << tool_index):
                selected &= bitmap
            else:
                selected &= ~bitmap
        # Some benchmark hosts still default to Python versions without
        # int.bit_count(); retain a slower compatibility fallback.
        counts[membership] = (
            selected.bit_count()
            if hasattr(selected, "bit_count")
            else bin(selected).count("1")
        )
    return counts


def collect(results_dir: Path, search_roots: list[Path]) -> list[dict[str, object]]:
    datasets = []
    for dataset, label, color in DATASETS:
        artifact_path = results_dir / f"{dataset}_accuracy_metrics.json"
        artifact = json.loads(artifact_path.read_text())
        input_records = int(artifact["metrics"][0]["input_records"])
        tool_values = []
        source_provenance = {}
        domain = None
        for tool in TOOLS:
            bitmap, metadata, sources = tool_bitmap(
                artifact, dataset, tool, search_roots
            )
            if int(metadata["numeric_id_max"]) != input_records:
                raise RuntimeError(
                    f"bitmap domain for {dataset}/{tool} does not equal input_records"
                )
            if domain is None:
                domain = metadata["accession_prefix"]
            elif metadata["accession_prefix"] != domain:
                raise RuntimeError(f"tool accession prefixes differ for {dataset}")
            tool_values.append(bitmap)
            source_provenance[tool] = sources

        counts = exclusive_counts(tuple(tool_values), input_records)
        union = sum(counts.values())
        expected_emitted = {
            str(item["tool"]): int(item["emitted_records"])
            for item in artifact["metrics"]
        }
        for tool_index, tool in enumerate(TOOLS):
            observed = sum(
                count
                for membership, count in counts.items()
                if membership & (1 << tool_index)
            )
            if observed != expected_emitted[tool]:
                raise RuntimeError(
                    f"emitted count mismatch for {dataset}/{tool}: "
                    f"artifact={expected_emitted[tool]}, intersections={observed}"
                )

        datasets.append(
            {
                "dataset": dataset,
                "label": label,
                "color": color,
                "input_records": input_records,
                "accession_prefix": domain,
                "union_records": union,
                "intersections": {
                    str(membership): {
                        "label": intersection_label,
                        "records": counts[membership],
                        "percent_of_union": 100.0 * counts[membership] / union,
                    }
                    for membership, intersection_label in INTERSECTIONS
                },
                "sources": source_provenance,
                "accuracy_artifact": str(artifact_path),
                "accuracy_artifact_sha256": sha256(artifact_path),
            }
        )
    return datasets


def percent_label(records: int, percent: float) -> str:
    if records == 0:
        return ""
    if percent < 0.01:
        return "<0.01%"
    return f"{percent:.2f}%"


def text_element(
    x: float,
    y: float,
    value: str,
    *,
    size: int = 14,
    weight: int = 400,
    anchor: str = "start",
    fill: str = "#202124",
    transform: str | None = None,
) -> str:
    transform_attr = f' transform="{escape(transform)}"' if transform else ""
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" '
        f'font-weight="{weight}" text-anchor="{anchor}" fill="{fill}"'
        f'{transform_attr}>{escape(value)}</text>'
    )


def render_svg(datasets: list[dict[str, object]]) -> str:
    width = 1220
    height = 840
    left = 215
    right = 25
    top = 44
    panel_height = 140
    plot_top_offset = 24
    plot_height = 82
    matrix_rule_y = 620
    matrix_heading_y = 651
    matrix_top = 682
    matrix_row_gap = 38
    centers = [
        left + (index + 0.5) * (width - left - right) / len(INTERSECTIONS)
        for index in range(len(INTERSECTIONS))
    ]
    column_width = (width - left - right) / len(INTERSECTIONS)
    bar_width = min(84.0, column_width * 0.62)

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" role="img">',
        "<title>Exclusive emitted-read intersections across preprocessing tools</title>",
        (
            "<desc>Faceted UpSet plot. Bars show each mutually exclusive tool "
            "intersection as a percentage of the reads emitted by any tool. "
            "The matrix below identifies tool membership.</desc>"
        ),
        '<rect width="100%" height="100%" fill="#FFFFFF"/>',
        (
            '<g font-family="Arial, Helvetica, sans-serif" '
            'shape-rendering="geometricPrecision">'
        ),
        text_element(
            20,
            30,
            "Exclusive intersections (% of any-tool union)",
            size=20,
            weight=600,
        ),
    ]

    for panel_index, dataset in enumerate(datasets):
        panel_top = top + panel_index * panel_height
        plot_top = panel_top + plot_top_offset
        baseline = plot_top + plot_height
        color = str(dataset["color"])
        svg.append(
            text_element(
                20,
                panel_top + 22,
                str(dataset["label"]),
                size=17,
                weight=700,
            )
        )
        svg.append(
            text_element(
                20,
                panel_top + 47,
                f'union = {int(dataset["union_records"]):,}',
                size=13,
                fill="#5F6368",
            )
        )
        for tick in (0, 50, 100):
            y = baseline - plot_height * tick / 100.0
            svg.append(
                f'<line x1="{left}" y1="{y:.1f}" x2="{width-right}" y2="{y:.1f}" '
                f'stroke="{("#9AA0A6" if tick == 0 else "#E3E6E8")}" '
                f'stroke-width="{(1.2 if tick == 0 else 0.8)}"/>'
            )
            if tick in (50, 100):
                svg.append(
                    text_element(
                        left - 10,
                        y + 5,
                        str(tick),
                        size=12,
                        anchor="end",
                        fill="#5F6368",
                    )
                )
        for column_index, (membership, _) in enumerate(INTERSECTIONS):
            entry = dataset["intersections"][str(membership)]
            records = int(entry["records"])
            percent = float(entry["percent_of_union"])
            bar_height = plot_height * percent / 100.0
            x = centers[column_index] - bar_width / 2
            y = baseline - bar_height
            if records:
                svg.append(
                    f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_width:.1f}" '
                    f'height="{max(bar_height, 0.8):.1f}" fill="{color}"/>'
                )
                # The label band above the 100% grid line is reserved inside
                # each panel so a full-height bar cannot collide with the
                # preceding panel or title.
                label_y = max(panel_top + 15, y - 9)
                svg.append(
                    text_element(
                        centers[column_index],
                        label_y,
                        percent_label(records, percent),
                        size=14,
                        weight=700,
                        anchor="middle",
                    )
                )

    svg.append(
        f'<line x1="{left}" y1="{matrix_rule_y}" x2="{width-right}" '
        f'y2="{matrix_rule_y}" stroke="#9AA0A6" stroke-width="1"/>'
    )
    svg.append(
        text_element(
            20,
            matrix_heading_y,
            "Tool membership",
            size=17,
            weight=700,
        )
    )
    row_positions = [matrix_top + index * matrix_row_gap for index in range(len(TOOLS))]
    for tool, y in zip(TOOLS, row_positions):
        svg.append(text_element(left - 22, y + 5, tool, size=14, anchor="end"))
        svg.append(
            f'<line x1="{left}" y1="{y:.1f}" x2="{width-right}" y2="{y:.1f}" '
            'stroke="#EEF0F2" stroke-width="1"/>'
        )

    for column_index, (membership, intersection_label) in enumerate(INTERSECTIONS):
        x = centers[column_index]
        included_rows = [
            row_positions[index]
            for index in range(len(TOOLS))
            if membership & (1 << index)
        ]
        if len(included_rows) > 1:
            svg.append(
                f'<line x1="{x:.1f}" y1="{min(included_rows):.1f}" '
                f'x2="{x:.1f}" y2="{max(included_rows):.1f}" '
                'stroke="#202124" stroke-width="3"/>'
            )
        for tool_index, y in enumerate(row_positions):
            included = bool(membership & (1 << tool_index))
            svg.append(
                f'<circle cx="{x:.1f}" cy="{y:.1f}" r="8" '
                f'fill="{("#202124" if included else "#DADCE0")}"/>'
            )
        words = intersection_label.replace(" + ", "+").split()
        svg.append(
            text_element(
                x,
                matrix_top + 3 * matrix_row_gap + 17,
                " ".join(words),
                size=12,
                anchor="middle",
                fill="#3C4043",
            )
        )

    svg.extend(["</g>", "</svg>"])
    return "\n".join(svg) + "\n"


def write_csv(path: Path, datasets: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "dataset",
                "intersection_mask",
                "intersection",
                *TOOLS,
                "records",
                "union_records",
                "percent_of_union",
            ),
            lineterminator="\n",
        )
        writer.writeheader()
        for dataset in datasets:
            for membership, intersection_label in INTERSECTIONS:
                entry = dataset["intersections"][str(membership)]
                writer.writerow(
                    {
                        "dataset": dataset["dataset"],
                        "intersection_mask": membership,
                        "intersection": intersection_label,
                        **{
                            tool: int(bool(membership & (1 << tool_index)))
                            for tool_index, tool in enumerate(TOOLS)
                        },
                        "records": entry["records"],
                        "union_records": dataset["union_records"],
                        "percent_of_union": f'{entry["percent_of_union"]:.10f}',
                    }
                )


def convert_svg(svg: Path, pdf: Path, png: Path) -> None:
    converter = shutil.which("rsvg-convert")
    if converter is None:
        print("rsvg-convert not found; wrote SVG only")
        return
    subprocess.run(
        [converter, "--format", "pdf", "--output", str(pdf), str(svg)], check=True
    )
    subprocess.run(
        [
            converter,
            "--format",
            "png",
            "--width",
            "2440",
            "--output",
            str(png),
            str(svg),
        ],
        check=True,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--output-prefix", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--bitmap-root",
        action="append",
        type=Path,
        default=[],
        help="fallback directory for relocated bitmap basenames (repeatable)",
    )
    parser.add_argument(
        "--svg-only", action="store_true", help="do not invoke rsvg-convert"
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    prefix = args.output_prefix.resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    datasets = collect(args.results_dir.resolve(), args.bitmap_root)

    json_path = prefix.with_suffix(".json")
    csv_path = prefix.with_suffix(".csv")
    svg_path = prefix.with_suffix(".svg")
    pdf_path = prefix.with_suffix(".pdf")
    png_path = prefix.with_suffix(".png")
    payload = {"schema_version": "1.0.0", "datasets": datasets}
    json_path.write_text(json.dumps(payload, indent=2) + "\n")
    write_csv(csv_path, datasets)
    svg_path.write_text(render_svg(datasets))
    if not args.svg_only:
        convert_svg(svg_path, pdf_path, png_path)

    for dataset in datasets:
        rounded = [
            f'{dataset["intersections"][str(mask)]["percent_of_union"]:.2f}'
            for mask, _ in INTERSECTIONS
        ]
        print(f'{dataset["dataset"]}: union={dataset["union_records"]:,}; ' + ", ".join(rounded))
    print(f"wrote {json_path}")
    print(f"wrote {csv_path}")
    print(f"wrote {svg_path}")
    if pdf_path.exists():
        print(f"wrote {pdf_path}")
    if png_path.exists():
        print(f"wrote {png_path}")


if __name__ == "__main__":
    main()
