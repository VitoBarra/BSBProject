from __future__ import annotations

import csv
import math
from pathlib import Path

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
    content = "".join(f"{identifier}\n" for identifier in identifiers)
    if path.is_file() and path.read_text(encoding="utf-8") == content:
        return
    path.write_text(content, encoding="utf-8")


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
        if gene_id:
            universe.append(gene_id)
        if math.isfinite(padj) and padj < padj_cutoff and math.isfinite(log2fc) and abs(log2fc) >= lfc_cutoff:
            selected.append(gene_id)

    universe = _unique(universe)
    selected = _unique(selected)
    if not universe:
        raise ValueError("No gene identifiers are available for the GO universe")
    if not selected:
        raise ValueError("No genes pass the configured adjusted p-value and log2 fold-change cutoffs")

    selected_genes_path.parent.mkdir(parents=True, exist_ok=True)
    _write_gene_ids(selected_genes_path, selected)
    _write_gene_ids(universe_genes_path, universe)
    return len(universe), len(selected)
