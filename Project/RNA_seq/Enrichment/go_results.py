from __future__ import annotations

import csv
import math
import textwrap
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ONTOLOGY_COLORS = {"BP": "#4C78A8", "MF": "#F58518", "CC": "#54A24B"}
GO_RESULT_COLUMNS = (
    "ID",
    "Description",
    "GeneRatio",
    "BgRatio",
    "RichFactor",
    "FoldEnrichment",
    "zScore",
    "pvalue",
    "p.adjust",
    "qvalue",
    "geneID",
    "Count",
    "ontology",
)


def _read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"Missing CSV header in {path}")
        return reader.fieldnames, list(reader)


def _number(value: str | None) -> float:
    if value is None or value.strip() in {"", "NA", "NaN", "nan"}:
        return math.nan
    return float(value)


def _strip_version(identifier: str) -> str:
    return identifier.split(".", 1)[0]


def _unique(values: list[str]) -> list[str]:
    return list(dict.fromkeys(value for value in values if value))


def _write_gene_ids(path: Path, identifiers: list[str]) -> None:
    path.write_text("".join(f"{identifier}\n" for identifier in identifiers), encoding="utf-8")


def prepare_go_inputs(
    de_results_path: Path,
    selected_genes_path: Path,
    universe_genes_path: Path,
    padj_cutoff: float,
    lfc_cutoff: float,
) -> tuple[int, int]:
    if not de_results_path.exists():
        raise FileNotFoundError(f"Missing DESeq2 result table: {de_results_path}")
    _, rows = _read_csv(de_results_path)
    if not rows:
        raise ValueError(f"No DESeq2 results found in {de_results_path}")
    required = {"gene_id", "padj", "log2FoldChange"}
    missing = required.difference(rows[0])
    if missing:
        raise ValueError(f"{de_results_path} is missing columns: {', '.join(sorted(missing))}")

    universe: list[str] = []
    selected: list[str] = []
    for row in rows:
        gene_id = _strip_version(row["gene_id"].strip())
        padj = _number(row["padj"])
        log2fc = _number(row["log2FoldChange"])
        if math.isfinite(padj):
            universe.append(gene_id)
        if math.isfinite(padj) and padj < padj_cutoff and math.isfinite(log2fc) and abs(log2fc) >= lfc_cutoff:
            selected.append(gene_id)

    universe = _unique(universe)
    selected = _unique(selected)
    if not universe:
        raise ValueError("No genes with a finite adjusted p-value are available for the GO universe")
    if not selected:
        raise ValueError("No genes pass the configured adjusted p-value and log2 fold-change cutoffs")

    selected_genes_path.parent.mkdir(parents=True, exist_ok=True)
    _write_gene_ids(selected_genes_path, selected)
    _write_gene_ids(universe_genes_path, universe)
    return len(universe), len(selected)


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _plot_dotplot(rows: list[dict[str, Any]], output_dir: Path) -> None:
    if not rows:
        return
    top_terms = rows[:20]
    descriptions = [textwrap.fill(str(row["Description"]), width=42) for row in top_terms]
    adjusted = [-math.log10(max(row["p.adjust"], float.fromhex("0x0.0000000000001p-1022"))) for row in top_terms]
    counts = [row["Count"] for row in top_terms]
    min_count = min(counts)
    max_count = max(counts)
    span = max_count - min_count
    sizes = [45 + (count - min_count) / span * 180 if span else 100 for count in counts]
    colors = [ONTOLOGY_COLORS.get(str(row["ontology"]), "#777777") for row in top_terms]

    figure, ax = plt.subplots(figsize=(9.5, 7))
    positions = list(range(len(top_terms)))
    ax.scatter(adjusted, positions, s=sizes, c=colors, alpha=0.82, edgecolors="black", linewidths=0.3)
    ax.set_yticks(positions, descriptions, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("-log10 adjusted p-value")
    ax.set_title("GO over-representation analysis: top significant terms")
    ax.grid(axis="x", alpha=0.25)

    ontology_handles = [
        plt.Line2D([0], [0], marker="o", linestyle="", color=color, label=ontology, markersize=7)
        for ontology, color in ONTOLOGY_COLORS.items()
        if ontology in {str(row["ontology"]) for row in top_terms}
    ]
    if ontology_handles:
        ontology_legend = ax.legend(
            handles=ontology_handles,
            title="Ontology",
            frameon=False,
            loc="lower left",
            bbox_to_anchor=(1.01, 0.05),
        )
        ax.add_artist(ontology_legend)
    representative_counts = sorted({min_count, round((min_count + max_count) / 2), max_count})
    size_handles = []
    for count in representative_counts:
        size = 45 + (count - min_count) / span * 180 if span else 100
        size_handles.append(ax.scatter([], [], s=size, color="#999999", alpha=0.7, label=str(count)))
    ax.legend(
        handles=size_handles,
        title="Genes",
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(1.01, 0.95),
    )
    figure.savefig(output_dir / "go_overrepresentation_dotplot.png", dpi=300, bbox_inches="tight")
    plt.close(figure)


def generate_go_outputs(
    all_results_path: Path,
    selected_genes_path: Path,
    universe_genes_path: Path,
    de_results_path: Path,
    padj_cutoff: float,
    lfc_cutoff: float,
    go_padj_cutoff: float = 0.05,
) -> None:
    for path in (all_results_path, selected_genes_path, universe_genes_path):
        if not path.exists():
            raise FileNotFoundError(f"Missing GO result input: {path}")

    fieldnames, rows = _read_csv(all_results_path)
    if rows:
        missing = {"p.adjust", "pvalue", "Count", "Description", "ontology"}.difference(rows[0])
        if missing:
            raise ValueError(f"{all_results_path} is missing columns: {', '.join(sorted(missing))}")
    elif fieldnames == [""]:
        fieldnames = list(GO_RESULT_COLUMNS)

    for row in rows:
        row["p.adjust"] = _number(row["p.adjust"])
        row["pvalue"] = _number(row["pvalue"])
        row["Count"] = int(_number(row["Count"]))
    rows.sort(
        key=lambda row: (
            row["p.adjust"] if math.isfinite(row["p.adjust"]) else math.inf,
            row["pvalue"] if math.isfinite(row["pvalue"]) else math.inf,
        )
    )
    significant = [row for row in rows if math.isfinite(row["p.adjust"]) and row["p.adjust"] < go_padj_cutoff]
    _write_csv(all_results_path, fieldnames, rows)
    _write_csv(all_results_path.parent / "go_overrepresentation_significant.csv", fieldnames, significant)

    universe_count = len(universe_genes_path.read_text(encoding="utf-8").splitlines())
    selected_count = len(selected_genes_path.read_text(encoding="utf-8").splitlines())
    with (all_results_path.parent / "enrichment_summary.txt").open("w", encoding="utf-8") as handle:
        handle.write("GO over-representation analysis\n")
        handle.write(f"Input DE table: {de_results_path}\n")
        handle.write(f"padj cutoff: {padj_cutoff}\n")
        handle.write(f"absolute log2FC cutoff: {lfc_cutoff}\n")
        handle.write(f"Tested genes: {universe_count}\n")
        handle.write(f"Significant genes: {selected_count}\n")
        handle.write(f"GO terms tested: {len(rows)}\n")
        handle.write(f"GO terms significant at padj < {go_padj_cutoff}: {len(significant)}\n")
    _plot_dotplot(significant, all_results_path.parent)
