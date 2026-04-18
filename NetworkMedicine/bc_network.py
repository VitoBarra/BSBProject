"""Run the BC workflow: download, rank genes, and write analysis outputs."""

from __future__ import annotations

import argparse
import json
import math
import statistics
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
from utils.tsv import write_tsv

from NetworkMedicine.bc_download import download_string_data
from NetworkMedicine.bc_metrics import build_graph, rank_and_select_genes
from NetworkMedicine.bc_null_model import modularity_null_model
from NetworkMedicine.bc_plots import plot_centralities, plot_largest_component, plot_modularity_null, plot_network


DEFAULT_GENES = Path("data/network_medicine/processed/bc_genes_s11.tsv")
DEFAULT_RAW_DIR = Path("data/network_medicine/raw/string_bc")
DEFAULT_RESULTS_DIR = Path("data/network_medicine/results/bc")
DEFAULT_API_URL = "https://version-12-0.string-db.org"


def analyze_network(
    mapping_path: Path,
    edges_path: Path,
    results_dir: Path,
    target_size: int,
    seed: int,
    null_iterations: int,
) -> None:
    graph, mapping_by_gene = build_graph(mapping_path, edges_path)
    results_dir.mkdir(parents=True, exist_ok=True)
    rows, membership, components, selected = rank_and_select_genes(
        graph, mapping_by_gene, target_size, seed
    )

    fields = [
        "gene",
        "ensembl_gene_id",
        "mapping_status",
        "string_id",
        "degree",
        "degree_centrality",
        "betweenness_centrality",
        "harmonic_centrality",
        "pagerank",
        "degree_percentile",
        "betweenness_percentile",
        "harmonic_percentile",
        "pagerank_percentile",
        "composite_percentile",
        "community",
        "component",
        "component_size",
        "selected",
        "selection_reason",
        "biological_rationale",
    ]
    write_tsv(results_dir / "bc_gene_centralities.tsv", rows, fields)
    write_tsv(
        results_dir / "bc_module_selected.tsv",
        [row for row in rows if bool(row["selected"])],
        fields,
    )
    write_tsv(
        results_dir / "bc_communities.tsv",
        sorted(
            [
                {
                    "gene": row["gene"],
                    "community": row["community"],
                    "component": row["component"],
                    "component_size": row["component_size"],
                    "degree": row["degree"],
                    "composite_percentile": row["composite_percentile"],
                }
                for row in rows
            ],
            key=lambda row: (int(row["community"]), -float(row["composite_percentile"])),
        ),
        ["gene", "community", "component", "component_size", "degree", "composite_percentile"],
    )

    modularity, community_sets, null_values = modularity_null_model(
        graph, membership, seed, null_iterations
    )
    null_mean = statistics.mean(null_values) if null_values else math.nan
    null_sd = statistics.stdev(null_values) if len(null_values) > 1 else math.nan
    modularity_z = (modularity - null_mean) / null_sd if null_sd and not math.isnan(null_sd) else math.nan
    empirical_p = (
        (1 + sum(value >= modularity for value in null_values)) / (1 + len(null_values))
        if null_values else math.nan
    )
    write_tsv(
        results_dir / "bc_modularity_null.tsv",
        [{"iteration": index + 1, "modularity": value} for index, value in enumerate(null_values)],
        ["iteration", "modularity"],
    )
    if null_values:
        plot_modularity_null(null_values, modularity, results_dir / "bc_modularity_null.png")
    stats = {
        "node_count": graph.number_of_nodes(),
        "mapped_node_count": sum(
            1 for row in rows if row["mapping_status"] == "mapped"
        ),
        "unmapped_node_count": sum(
            1 for row in rows if row["mapping_status"] != "mapped"
        ),
        "edge_count": graph.number_of_edges(),
        "density": nx.density(graph),
        "average_clustering_weighted": nx.average_clustering(graph, weight="score"),
        "connected_component_count": len(components),
        "component_sizes": [len(component) for component in components],
        "largest_component_size": len(components[0]) if components else 0,
        "isolated_nodes": sorted(nx.isolates(graph)),
        "community_count_including_isolates": len(set(membership.values())),
        "non_isolate_community_count": len(community_sets),
        "louvain_modularity_non_isolates": modularity,
        "modularity_null_model": {
            "method": "degree-preserving double-edge swaps with Louvain redetection and observed edge weights reassigned in edge iteration order",
            "requested_iterations": null_iterations,
            "successful_iterations": len(null_values),
            "mean": null_mean,
            "standard_deviation": null_sd,
            "z_score": modularity_z,
            "empirical_upper_tail_p": empirical_p,
        },
        "selected_module_size": len(selected),
        "selected_genes": [str(row["gene"]) for row in rows if bool(row["selected"])],
        "random_seed": seed,
        "software_versions": {
            "python": sys.version.split()[0],
            "networkx": nx.__version__,
            "matplotlib": plt.matplotlib.__version__,
        },
        "centrality_policy": {
            "degree": "unweighted",
            "betweenness": "unweighted shortest paths",
            "harmonic": "unweighted and divided by N-1",
            "pagerank": "weighted by the class-filtered experiments/databases STRING score",
            "composite": "mean of the four within-network percentile ranks",
        },
        "selection_policy": (
            "Top composite-ranked representative of every Louvain community "
            "containing an edge, then highest composite ranks until target size."
        ),
        "selection_status": (
            "Final reproducible network-based module with documented biological rationale."
        ),
    }
    (results_dir / "bc_network_summary.json").write_text(
        json.dumps(stats, indent=2) + "\n",
        encoding="utf-8",
    )

    rows_by_gene = {str(row["gene"]): row for row in rows}
    for node in graph:
        graph.nodes[node].update(
            community=int(rows_by_gene[node]["community"]),
            composite_percentile=float(rows_by_gene[node]["composite_percentile"]),
            selected=bool(rows_by_gene[node]["selected"]),
        )
    nx.write_graphml(graph, results_dir / "bc_network.graphml")
    plot_network(graph, rows_by_gene, results_dir / "bc_network.png", seed)
    plot_largest_component(
        graph,
        rows_by_gene,
        results_dir / "bc_largest_component.png",
        seed,
    )
    plot_centralities(rows, results_dir / "bc_top20_composite_rank.png")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Download and analyze the 72-gene BC STRING network.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    download = subparsers.add_parser("download", help="Map genes and download the induced STRING network.")
    download.add_argument("--genes", type=Path, default=DEFAULT_GENES)
    download.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR)
    download.add_argument(
        "--processed-dir",
        type=Path,
        default=Path("data/network_medicine/processed"),
    )
    download.add_argument("--score-threshold", type=float, default=0.70)
    download.add_argument("--network-type", choices=("functional", "physical"), default="physical")
    download.add_argument("--api-url", default=DEFAULT_API_URL)

    analyze = subparsers.add_parser("analyze", help="Calculate centralities, community significance, and the final module.")
    analyze.add_argument(
        "--mapping",
        type=Path,
        default=Path("data/network_medicine/processed/bc_string_mapping.tsv"),
    )
    analyze.add_argument(
        "--edges",
        type=Path,
        default=Path("data/network_medicine/processed/bc_ppi_edges.tsv"),
    )
    analyze.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    analyze.add_argument("--target-size", type=int, default=10)
    analyze.add_argument("--seed", type=int, default=42)
    analyze.add_argument("--null-iterations", type=int, default=500)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.command == "download":
        if not 0 <= args.score_threshold <= 1:
            raise ValueError("--score-threshold must be between 0 and 1")
        download_string_data(
            args.genes,
            args.raw_dir,
            args.processed_dir,
            args.score_threshold,
            args.network_type,
            args.api_url,
        )
    elif args.command == "analyze":
        if not 10 <= args.target_size <= 20:
            raise ValueError("--target-size must be between 10 and 20")
        analyze_network(
            args.mapping,
            args.edges,
            args.results_dir,
            args.target_size,
            args.seed,
            args.null_iterations,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
