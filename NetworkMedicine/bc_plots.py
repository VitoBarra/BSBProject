"""Visualizations for the breast-cancer network analysis."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx


def plot_network(graph: nx.Graph, rows_by_gene: dict[str, dict[str, object]], output: Path, seed: int) -> None:
    plt.figure(figsize=(14, 11))
    position = nx.spring_layout(graph, seed=seed, weight="score", k=1.2 / math.sqrt(max(len(graph), 1)))
    node_sizes = [250 + 2500 * float(rows_by_gene[node]["pagerank"]) for node in graph]
    node_colors = [int(rows_by_gene[node]["community"]) for node in graph]
    edge_widths = [0.4 + 2.2 * float(data["score"]) for _, _, data in graph.edges(data=True)]
    nx.draw_networkx_edges(graph, position, width=edge_widths, alpha=0.45, edge_color="#64748b")
    nx.draw_networkx_nodes(
        graph, position, node_size=node_sizes, node_color=node_colors,
        cmap=plt.colormaps["tab20"], edgecolors="#172033", linewidths=0.7,
    )
    labels = {
        node: node for node in graph
        if bool(rows_by_gene[node]["selected"]) or graph.degree(node) == 0
    }
    nx.draw_networkx_labels(graph, position, labels=labels, font_size=8)
    plt.title("BC Table S11 induced STRING network")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(output, dpi=220, bbox_inches="tight")
    plt.close()


def plot_largest_component(
    graph: nx.Graph, rows_by_gene: dict[str, dict[str, object]], output: Path, seed: int,
) -> None:
    largest_nodes = max(nx.connected_components(graph), key=len)
    core = graph.subgraph(largest_nodes)
    plt.figure(figsize=(15, 12))
    position = nx.spring_layout(core, seed=seed, weight="score", k=0.75, iterations=300)
    node_sizes = [500 + 5500 * float(rows_by_gene[node]["pagerank"]) for node in core]
    node_colors = [int(rows_by_gene[node]["community"]) for node in core]
    edge_widths = [0.5 + 2.5 * float(data["score"]) for _, _, data in core.edges(data=True)]
    nx.draw_networkx_edges(core, position, width=edge_widths, alpha=0.5, edge_color="#64748b")
    nx.draw_networkx_nodes(
        core, position, node_size=node_sizes, node_color=node_colors, cmap=plt.colormaps["tab10"],
        edgecolors=["#d62728" if bool(rows_by_gene[node]["selected"]) else "#172033" for node in core],
        linewidths=[2.5 if bool(rows_by_gene[node]["selected"]) else 0.8 for node in core],
    )
    nx.draw_networkx_labels(core, position, font_size=9)
    plt.title("Largest connected component of the BC STRING network\n(red border = final module selection)")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(output, dpi=220, bbox_inches="tight")
    plt.close()


def plot_centralities(rows: list[dict[str, object]], output: Path) -> None:
    ranked = sorted(rows, key=lambda row: float(row["composite_percentile"]), reverse=True)[:20]
    genes = [str(row["gene"]) for row in ranked][::-1]
    scores = [float(row["composite_percentile"]) for row in ranked][::-1]
    colors = ["#d95f02" if bool(row["selected"]) else "#4c78a8" for row in ranked][::-1]
    plt.figure(figsize=(10, 8))
    plt.barh(genes, scores, color=colors)
    plt.xlabel("Composite centrality percentile")
    plt.ylabel("Gene")
    plt.title("Top 20 BC genes by composite network rank")
    plt.xlim(0, 1.03)
    plt.tight_layout()
    plt.savefig(output, dpi=220, bbox_inches="tight")
    plt.close()


def plot_modularity_null(null_values: list[float], modularity: float, output: Path) -> None:
    plt.figure(figsize=(8, 5))
    plt.hist(null_values, bins=25, color="#4c78a8", alpha=0.85)
    plt.axvline(modularity, color="#d62728", linewidth=2, label=f"Observed = {modularity:.3f}")
    plt.xlabel("Maximum Louvain modularity")
    plt.ylabel("Randomized networks")
    plt.title("Degree-preserving modularity null model")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output, dpi=220, bbox_inches="tight")
    plt.close()
