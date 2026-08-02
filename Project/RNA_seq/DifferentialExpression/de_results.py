from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.colors as mcolors  # noqa: E402
import matplotlib.patches as mpatches  # noqa: E402

RESULT_COLUMNS = (
    "gene_id",
    "gene_symbol",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
)
NUMERIC_RESULT_COLUMNS = RESULT_COLUMNS[2:]
UP_COLOR = "#B7352D"
DOWN_COLOR = "#5B8DD9"
NS_COLOR = "#BDBDBD"
CONDITION_COLORS = {"normal": "#222222", "tumor": "#F5A000"}


def _read_table(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"Missing header in {path}")
        return list(reader)


def _number(value: str | None) -> float:
    if value is None or value.strip() in {"", "NA", "NaN", "nan"}:
        return math.nan
    return float(value)


def _read_results(path: Path) -> list[dict[str, Any]]:
    rows = _read_table(path)
    if not rows:
        raise ValueError(f"No DESeq2 results found in {path}")
    missing = set(RESULT_COLUMNS).difference(rows[0])
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
    for row in rows:
        for column in NUMERIC_RESULT_COLUMNS:
            row[column] = _number(row[column])
    return rows


def _read_expression_matrix(path: Path, sample_names: list[str]) -> tuple[list[dict[str, str]], np.ndarray]:
    rows = _read_table(path)
    if not rows:
        raise ValueError(f"No expression values found in {path}")
    missing = {"gene_id", "gene_symbol", *sample_names}.difference(rows[0])
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
    matrix = np.asarray([[_number(row[name]) for name in sample_names] for row in rows], dtype=float)
    if not np.isfinite(matrix).all():
        raise ValueError(f"Non-finite expression value found in {path}")
    return rows, matrix


def _label(row: dict[str, Any]) -> str:
    symbol = str(row.get("gene_symbol", "")).strip()
    return symbol if symbol and symbol != "NA" else str(row["gene_id"])


def _direction(row: dict[str, Any], padj_cutoff: float = 0.05, lfc_cutoff: float = 1.0) -> str:
    padj = row["padj"]
    lfc = row["log2FoldChange"]
    if not (math.isfinite(padj) and math.isfinite(lfc) and padj < padj_cutoff):
        return "Not significant"
    if lfc >= lfc_cutoff:
        return "Up"
    if lfc <= -lfc_cutoff:
        return "Down"
    return "Not significant"


def _write_results(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    column: "" if isinstance(row[column], float) and math.isnan(row[column]) else row[column]
                    for column in RESULT_COLUMNS
                }
            )


def _save_figure(figure: plt.Figure, output_dir: Path, stem: str) -> None:
    figure.savefig(output_dir / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(figure)


def _annotate(ax: plt.Axes, rows: list[dict[str, Any]], x_key: str, y_values: dict[str, float]) -> None:
    for index, row in enumerate(rows):
        ax.annotate(
            _label(row),
            (row[x_key], y_values[str(row["gene_id"])]),
            xytext=(4, 5 + (index % 3) * 3),
            textcoords="offset points",
            fontsize=7,
        )


def _plot_volcano(rows: list[dict[str, Any]], output_dir: Path) -> None:
    finite_rows = [row for row in rows if math.isfinite(row["log2FoldChange"]) and math.isfinite(row["padj"])]
    tiny = np.finfo(float).tiny
    y_values = {str(row["gene_id"]): -math.log10(max(row["padj"], tiny)) for row in finite_rows}
    figure, ax = plt.subplots(figsize=(7.2, 6.2))
    for direction, color in (("Down", DOWN_COLOR), ("Not significant", NS_COLOR), ("Up", UP_COLOR)):
        selected = [row for row in finite_rows if _direction(row) == direction]
        ax.scatter(
            [row["log2FoldChange"] for row in selected],
            [y_values[str(row["gene_id"])] for row in selected],
            s=12,
            alpha=0.75,
            color=color,
            label=direction,
            edgecolors="none",
        )
    ax.axvline(-1, color="black", linewidth=0.6, linestyle="--")
    ax.axvline(1, color="black", linewidth=0.6, linestyle="--")
    ax.axhline(-math.log10(0.05), color="black", linewidth=0.6, linestyle="--")
    labels: list[dict[str, Any]] = []
    for direction in ("Down", "Up"):
        labels.extend(sorted((row for row in finite_rows if _direction(row) == direction), key=lambda row: row["padj"])[:10])
    _annotate(ax, labels, "log2FoldChange", y_values)
    ax.set(xlabel="log2 fold change", ylabel="-log10 adjusted p-value")
    ax.set_title("Paired DESeq2 model; padj < 0.05 and |log2FC| >= 1", fontsize=9, pad=9)
    figure.suptitle("Differential expression: tumor vs adjacent normal", fontsize=16, y=0.98)
    figure.subplots_adjust(top=0.88)
    ax.legend(frameon=False)
    _save_figure(figure, output_dir, "volcano_padj")


def _plot_ma(rows: list[dict[str, Any]], output_dir: Path) -> None:
    finite_rows = [
        row
        for row in rows
        if math.isfinite(row["baseMean"]) and row["baseMean"] > 0 and math.isfinite(row["log2FoldChange"])
    ]
    x_values = {str(row["gene_id"]): math.log2(row["baseMean"]) for row in finite_rows}
    figure, ax = plt.subplots(figsize=(8.5, 5.5))
    for direction, color in (("Down", DOWN_COLOR), ("Not significant", NS_COLOR), ("Up", UP_COLOR)):
        selected = [row for row in finite_rows if _direction(row) == direction]
        ax.scatter(
            [x_values[str(row["gene_id"])] for row in selected],
            [row["log2FoldChange"] for row in selected],
            s=9,
            alpha=0.75,
            color=color,
            label=f"{direction}: {len(selected)}" if direction != "Not significant" else "NS",
            edgecolors="none",
        )
    ax.axhline(0, color="black", linewidth=0.7)
    ax.axhline(-1, color="black", linewidth=0.6, linestyle="--")
    ax.axhline(1, color="black", linewidth=0.6, linestyle="--")
    labels: list[dict[str, Any]] = []
    for direction in ("Down", "Up"):
        labels.extend(
            sorted(
                (row for row in finite_rows if _direction(row) == direction),
                key=lambda row: row["padj"] if math.isfinite(row["padj"]) else math.inf,
            )[:9]
        )
    for row in labels:
        ax.annotate(
            _label(row),
            (x_values[str(row["gene_id"])], row["log2FoldChange"]),
            xytext=(4, 5),
            textcoords="offset points",
            fontsize=7,
        )
    ax.set(xlabel="log2 mean normalized expression", ylabel="log2 fold change")
    ax.set_title("Paired DESeq2 model; padj < 0.05 and |log2FC| >= 1", fontsize=9, pad=9)
    figure.suptitle("MA plot: tumor vs adjacent normal", fontsize=16, y=0.98)
    figure.subplots_adjust(top=0.88)
    ax.legend(frameon=False, ncol=3, loc="upper center")
    _save_figure(figure, output_dir, "ma_plot")


def _plot_pca(vst_matrix: np.ndarray, samples: list[dict[str, str]], output_dir: Path) -> None:
    sample_matrix = vst_matrix.T
    centered = sample_matrix - sample_matrix.mean(axis=0, keepdims=True)
    u, singular_values, _ = np.linalg.svd(centered, full_matrices=False)
    coordinates = u[:, :2] * singular_values[:2]
    variances = singular_values**2
    explained = variances / variances.sum() if variances.sum() else np.zeros_like(variances)
    if coordinates.shape[1] < 2:
        coordinates = np.pad(coordinates, ((0, 0), (0, 2 - coordinates.shape[1])))
        explained = np.pad(explained, (0, 2 - len(explained)))

    figure, ax = plt.subplots(figsize=(6, 5))
    conditions = sorted({sample["condition"] for sample in samples})
    for condition in conditions:
        indices = [index for index, sample in enumerate(samples) if sample["condition"] == condition]
        ax.scatter(
            coordinates[indices, 0],
            coordinates[indices, 1],
            color=CONDITION_COLORS.get(condition),
            s=45,
            label=condition,
        )
    for index, sample in enumerate(samples):
        ax.annotate(sample["sample_name"], coordinates[index, :2], xytext=(4, 4), textcoords="offset points", fontsize=7)
    ax.set(
        title="PCA of VST expression",
        xlabel=f"PC1 ({explained[0] * 100:.1f}%)",
        ylabel=f"PC2 ({explained[1] * 100:.1f}%)",
    )
    ax.legend(title="Condition", frameon=False)
    _save_figure(figure, output_dir, "pca_vst")


def _plot_sample_distances(vst_matrix: np.ndarray, sample_names: list[str], output_dir: Path) -> None:
    sample_matrix = vst_matrix.T
    differences = sample_matrix[:, np.newaxis, :] - sample_matrix[np.newaxis, :, :]
    distances = np.sqrt(np.sum(differences**2, axis=2))
    figure, ax = plt.subplots(figsize=(6.5, 5.8))
    image = ax.imshow(distances, cmap="Blues", aspect="equal")
    positions = np.arange(len(sample_names))
    ax.set_xticks(positions, sample_names, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(positions, sample_names, fontsize=8)
    ax.set_title("Sample distance heatmap")
    figure.colorbar(image, ax=ax, label="Distance")
    _save_figure(figure, output_dir, "sample_distance_heatmap")


def _unique_labels(rows: list[dict[str, Any]]) -> list[str]:
    seen: dict[str, int] = {}
    labels: list[str] = []
    for row in rows:
        label = _label(row)
        seen[label] = seen.get(label, 0) + 1
        labels.append(label if seen[label] == 1 else f"{label}.{seen[label]}")
    return labels


@dataclass(slots=True)
class _ClusterNode:
    leaves: list[int]
    height: float = 0.0
    left: _ClusterNode | None = None
    right: _ClusterNode | None = None


def _cluster_tree(matrix: np.ndarray) -> _ClusterNode:
    """Build a complete-linkage clustering tree without extra dependencies."""
    count = matrix.shape[0]
    if count == 0:
        return _ClusterNode([])
    differences = matrix[:, np.newaxis, :] - matrix[np.newaxis, :, :]
    distances = np.sqrt(np.sum(differences**2, axis=2))
    clusters = [_ClusterNode([index]) for index in range(count)]
    while len(clusters) > 1:
        best_pair = (0, 1)
        best_distance = math.inf
        for left in range(len(clusters) - 1):
            for right in range(left + 1, len(clusters)):
                complete = float(np.max(distances[np.ix_(clusters[left].leaves, clusters[right].leaves)]))
                if complete < best_distance:
                    best_distance = complete
                    best_pair = (left, right)
        left, right = best_pair
        left_node = clusters[left]
        right_node = clusters[right]
        merged = _ClusterNode(
            leaves=left_node.leaves + right_node.leaves,
            height=best_distance,
            left=left_node,
            right=right_node,
        )
        clusters = [cluster for index, cluster in enumerate(clusters) if index not in best_pair]
        clusters.append(merged)
    return clusters[0]


def _cluster_order(matrix: np.ndarray) -> list[int]:
    return _cluster_tree(matrix).leaves


def _draw_dendrogram(ax: plt.Axes, root: _ClusterNode, *, orientation: str) -> None:
    positions = {leaf: position for position, leaf in enumerate(root.leaves)}

    def draw(node: _ClusterNode) -> float:
        if node.left is None or node.right is None:
            return float(positions[node.leaves[0]])
        left_position = draw(node.left)
        right_position = draw(node.right)
        if orientation == "row":
            ax.plot(
                [node.left.height, node.height, node.height, node.right.height],
                [left_position, left_position, right_position, right_position],
                color="black",
                linewidth=0.7,
            )
        else:
            ax.plot(
                [left_position, left_position, right_position, right_position],
                [node.left.height, node.height, node.height, node.right.height],
                color="black",
                linewidth=0.7,
            )
        return (left_position + right_position) / 2

    if len(root.leaves) > 1:
        draw(root)
    ax.set_axis_off()
    if orientation == "row":
        ax.set_xlim(max(root.height * 1.05, 1.0), 0)
        ax.set_ylim(len(root.leaves) - 0.5, -0.5)
    else:
        ax.set_xlim(-0.5, len(root.leaves) - 0.5)
        ax.set_ylim(0, max(root.height * 1.05, 1.0))


def _plot_top_gene_heatmap(
    results: list[dict[str, Any]],
    vst_rows: list[dict[str, str]],
    vst_matrix: np.ndarray,
    samples: list[dict[str, str]],
    output_dir: Path,
) -> None:
    sample_names = [sample["sample_name"] for sample in samples]
    ranked = sorted(
        (row for row in results if math.isfinite(row["padj"])),
        key=lambda row: (row["padj"], -abs(row["log2FoldChange"]) if math.isfinite(row["log2FoldChange"]) else 0),
    )[:20]
    if len(ranked) < 2:
        return
    row_indices = {row["gene_id"]: index for index, row in enumerate(vst_rows)}
    selected = [row for row in ranked if str(row["gene_id"]) in row_indices]
    if len(selected) < 2:
        return
    matrix = np.asarray([vst_matrix[row_indices[str(row["gene_id"])]] for row in selected])
    means = matrix.mean(axis=1, keepdims=True)
    standard_deviations = matrix.std(axis=1, keepdims=True)
    row_z = np.divide(matrix - means, standard_deviations, out=np.zeros_like(matrix), where=standard_deviations != 0)
    row_tree = _cluster_tree(row_z)
    column_tree = _cluster_tree(row_z.T)
    row_order = row_tree.leaves
    column_order = column_tree.leaves
    row_z = row_z[np.ix_(row_order, column_order)]
    selected = [selected[index] for index in row_order]
    ordered_sample_names = [sample_names[index] for index in column_order]
    conditions_by_sample = {sample["sample_name"]: sample["condition"] for sample in samples}
    ordered_conditions = [conditions_by_sample[name] for name in ordered_sample_names]

    figure = plt.figure(figsize=(8.2, 8.8))
    grid = figure.add_gridspec(
        3,
        4,
        width_ratios=(1.3, 6.5, 1.1, 0.25),
        height_ratios=(1.2, 0.14, 7),
        hspace=0.015,
        wspace=0.02,
    )
    column_tree_ax = figure.add_subplot(grid[0, 1])
    annotation_ax = figure.add_subplot(grid[1, 1])
    legend_ax = figure.add_subplot(grid[0:2, 2:4])
    row_tree_ax = figure.add_subplot(grid[2, 0])
    ax = figure.add_subplot(grid[2, 1])
    colorbar_ax = figure.add_subplot(grid[2, 3])
    _draw_dendrogram(column_tree_ax, column_tree, orientation="column")
    _draw_dendrogram(row_tree_ax, row_tree, orientation="row")

    annotation_colors = np.asarray(
        [[mcolors.to_rgb(CONDITION_COLORS.get(condition, "#777777")) for condition in ordered_conditions]]
    )
    annotation_ax.imshow(annotation_colors, aspect="auto")
    annotation_ax.set_xticks([])
    annotation_ax.set_yticks([])
    for spine in annotation_ax.spines.values():
        spine.set_linewidth(0.6)

    condition_handles = [
        mpatches.Patch(color=CONDITION_COLORS.get(condition, "#777777"), label=condition.capitalize())
        for condition in dict.fromkeys(ordered_conditions)
    ]
    legend_ax.legend(handles=condition_handles, title="Condition", frameon=False, loc="lower left")
    legend_ax.set_axis_off()

    image = ax.imshow(row_z, cmap="coolwarm", vmin=-2.5, vmax=2.5, aspect="auto")
    ax.set_xticks(np.arange(len(ordered_sample_names)), ordered_sample_names, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(selected)), _unique_labels(selected), fontsize=8)
    ax.yaxis.tick_right()
    ax.tick_params(axis="y", labelright=True, labelleft=False, length=0, pad=4)
    figure.suptitle("Top 20 differentially expressed genes (row z-scores)", fontsize=14, y=0.995)
    figure.colorbar(image, cax=colorbar_ax, label="Row z-score")
    _save_figure(figure, output_dir, "top_de_gene_heatmap")


def generate_de_outputs(results_path: Path, normalized_counts_path: Path, vst_counts_path: Path, samples_path: Path) -> None:
    output_dir = results_path.parent
    for path in (results_path, normalized_counts_path, vst_counts_path, samples_path):
        if not path.exists():
            raise FileNotFoundError(f"Missing DE plotting input: {path}")

    results = _read_results(results_path)
    samples = _read_table(samples_path, delimiter="\t")
    required_sample_columns = {"sample_name", "patient", "condition"}
    if not samples or not required_sample_columns.issubset(samples[0]):
        raise ValueError(f"{samples_path} must contain columns: {', '.join(sorted(required_sample_columns))}")
    sample_names = [sample["sample_name"] for sample in samples]
    _, _ = _read_expression_matrix(normalized_counts_path, sample_names)
    vst_rows, vst_matrix = _read_expression_matrix(vst_counts_path, sample_names)

    significant = [row for row in results if math.isfinite(row["padj"]) and row["padj"] < 0.05]
    ranked = sorted(
        (row for row in results if math.isfinite(row["stat"])),
        key=lambda row: row["stat"],
        reverse=True,
    )
    _write_results(output_dir / "deseq2_significant_genes_padj_0.05.csv", significant)
    _write_results(output_dir / "deseq2_ranked_genes.csv", ranked)
    _plot_volcano(results, output_dir)
    _plot_ma(results, output_dir)
    _plot_pca(vst_matrix, samples, output_dir)
    _plot_sample_distances(vst_matrix, sample_names, output_dir)
    _plot_top_gene_heatmap(results, vst_rows, vst_matrix, samples, output_dir)
