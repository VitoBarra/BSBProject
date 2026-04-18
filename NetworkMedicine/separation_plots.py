"""Visualize disease separation and representative network paths."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import networkx as nx


def plot_separation_results(
    rows: list[dict[str, object]], graph: nx.Graph,
    projected: dict[str, set[str]], results_dir: Path,
) -> None:
    """Plot separation, path coverage, and a representative subnetwork."""
    comparisons = [str(row["comparison"]) for row in rows]
    separations = [float(row["separation_s_ab"]) for row in rows]
    plt.figure(figsize=(7, 5))
    bars = plt.bar(comparisons, separations, color=["#7b3294", "#008837"])
    plt.axhline(0, color="#222222", linewidth=0.8)
    plt.ylabel("Network separation $S_{AB}$")
    plt.title("Breast-cancer disease-module separation")
    for bar, value in zip(bars, separations, strict=True):
        plt.text(bar.get_x() + bar.get_width() / 2, value, f"{value:.3f}", ha="center", va="bottom")
    plt.tight_layout()
    plt.savefig(results_dir / "disease_separation.png", dpi=220, bbox_inches="tight")
    plt.close()

    coverage = [100 * float(row["cross_finite_fraction"]) for row in rows]
    plt.figure(figsize=(7, 5))
    bars = plt.bar(comparisons, coverage, color=["#7b3294", "#008837"])
    plt.ylim(0, 105)
    plt.ylabel("Finite cross-module pairs (%)")
    plt.title("Shortest-path coverage")
    for bar, value in zip(bars, coverage, strict=True):
        plt.text(bar.get_x() + bar.get_width() / 2, value, f"{value:.1f}%", ha="center", va="bottom")
    plt.tight_layout()
    plt.savefig(results_dir / "finite_path_coverage.png", dpi=220, bbox_inches="tight")
    plt.close()

    role = {node: "Background" for node in graph}
    for disease, genes in projected.items():
        for gene in genes:
            role[gene] = disease if role[gene] == "Background" else f"{role[gene]}+{disease}"
    core_nodes = max(nx.connected_components(graph), key=len, default=set())
    core = graph.subgraph(core_nodes)
    # Display a compact, interpretable view: for every BC gene in the largest
    # component, retain one deterministic shortest path to its nearest RA gene
    # and one to its nearest DM gene. Calculations and GraphML still use the
    # complete graph above; this reduction affects visualization only.
    display_nodes: set[str] = set()
    for bc_gene in sorted(projected["BC"] & set(core)):
        lengths = nx.single_source_shortest_path_length(core, bc_gene)
        display_nodes.add(bc_gene)
        for disease in ("RA", "DM"):
            candidates = projected[disease] & set(lengths)
            if not candidates:
                continue
            target = min(candidates, key=lambda node: (lengths[node], node))
            display_nodes.update(nx.shortest_path(core, bc_gene, target))
    display = core.subgraph(display_nodes)
    position = nx.spring_layout(display, seed=42, weight="score", k=1.0, iterations=500)
    colors = {
        "Background": "#d9d9d9", "BC": "#e41a1c", "RA": "#984ea3", "DM": "#4daf4a",
        "RA+DM": "#377eb8", "BC+RA": "#ff7f00", "BC+DM": "#a65628", "BC+RA+DM": "#111111",
    }
    plt.figure(figsize=(12, 10))
    nx.draw_networkx_edges(display, position, alpha=0.45, width=1.0, edge_color="#777777")
    nx.draw_networkx_nodes(display, position, node_size=[90 if role[n] == "Background" else 240 for n in display],
                           node_color=[colors.get(role[n], "#377eb8") for n in display], linewidths=0.5,
                           edgecolors="#333333")
    labels = {node: node for node in display if role[node] != "Background"}
    nx.draw_networkx_labels(
        display, position, labels=labels, font_size=7,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.55, "pad": 0.15},
    )
    plt.legend(
        handles=[
            Patch(facecolor=colors["BC"], label="BC"),
            Patch(facecolor=colors["RA"], label="RA"),
            Patch(facecolor=colors["DM"], label="DM"),
            Patch(facecolor=colors["BC+RA"], label="Module overlap"),
            Patch(facecolor=colors["Background"], edgecolor="#333333", label="Background connector"),
        ],
        loc="upper right",
        frameon=False,
    )
    plt.title("Representative shortest paths from BC to nearest RA and DM genes")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(results_dir / "bc_ra_dm_joint_network.png", dpi=220, bbox_inches="tight")
    plt.close()
