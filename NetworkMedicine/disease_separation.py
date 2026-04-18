"""Build a joint BC/RA/DM STRING network and calculate disease separation."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import networkx as nx
import matplotlib.pyplot as plt
from utils.tsv import read_tsv, write_tsv

from NetworkMedicine.separation_download import download_network
from NetworkMedicine.separation_metrics import build_joint_graph, calculate_separations
from NetworkMedicine.separation_plots import plot_separation_results


def read_symbols(path: Path, column: str, selected_only: bool = False) -> list[str]:
    rows = read_tsv(path)
    if not rows or column not in rows[0]:
        raise ValueError(f"{path} does not contain {column!r}")
    if selected_only and "selected" in rows[0]:
        rows = [row for row in rows if row["selected"].strip().lower() in {"true", "1", "yes"}]
    symbols = [row[column].strip() for row in rows if row[column].strip()]
    return list(dict.fromkeys(symbols))


def load_modules(bc_path: Path, ra_path: Path, dm_path: Path) -> dict[str, list[str]]:
    modules = {
        "BC": read_symbols(bc_path, "gene", selected_only=True),
        "RA": read_symbols(ra_path, "gene_symbol"),
        "DM": read_symbols(dm_path, "gene_symbol"),
    }
    if not 10 <= len(modules["BC"]) <= 20:
        raise ValueError("The BC module must contain 10–20 selected genes")
    return modules


def analyze(modules: dict[str, list[str]], mapping_path: Path, edges_path: Path, results_dir: Path) -> None:
    graph, mapped = build_joint_graph(mapping_path, edges_path)
    projected = {name: set(genes) & mapped for name, genes in modules.items()}

    rows = calculate_separations(graph, modules, projected)
    results_dir.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    write_tsv(results_dir / "disease_separation.tsv", rows, fields)
    nx.write_graphml(graph, results_dir / "bc_ra_dm_joint_network.graphml")
    plot_separation_results(rows, graph, projected, results_dir)
    write_network_summary(graph, projected, results_dir)


def write_network_summary(graph: nx.Graph, projected: dict[str, set[str]], results_dir: Path) -> None:
    summary = {
        "network_nodes": graph.number_of_nodes(), "network_edges": graph.number_of_edges(),
        "connected_components": nx.number_connected_components(graph),
        "largest_component_size": max(map(len, nx.connected_components(graph)), default=0),
        "isolated_nodes": len(list(nx.isolates(graph))),
        "background_nodes": len(set(graph) - set().union(*projected.values())),
        "software_versions": {
            "python": sys.version.split()[0],
            "networkx": nx.__version__,
            "matplotlib": plt.matplotlib.__version__,
        },
        "distance_policy": "Unweighted shortest paths; all-pairs average distance as taught in class; self-pairs included for within-module distances; unreachable pairs excluded and reported.",
        "separation_formula": "s_AB = d_AB - (d_AA + d_BB) / 2",
        "interpretation": "Negative values indicate topological overlap; positive values indicate separation.",
    }
    (results_dir / "disease_separation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


def parser() -> argparse.ArgumentParser:
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--bc", type=Path, default=Path("data/network_medicine/results/bc/bc_module_selected.tsv"))
    common.add_argument("--ra", type=Path, default=Path("data/network_medicine/processed/ra_genes_disgenet_curated.tsv"))
    common.add_argument("--dm", type=Path, default=Path("data/network_medicine/processed/dm_genes_disgenet_curated.tsv"))
    root = argparse.ArgumentParser(description=__doc__)
    sub = root.add_subparsers(dest="command", required=True)
    download = sub.add_parser("download", parents=[common])
    download.add_argument("--raw-dir", type=Path, default=Path("data/network_medicine/raw/string_joint"))
    download.add_argument("--processed-dir", type=Path, default=Path("data/network_medicine/processed"))
    download.add_argument("--score-threshold", type=float, default=0.70)
    download.add_argument("--network-type", choices=("functional", "physical"), default="physical")
    download.add_argument("--additional-interactors", type=int, default=650)
    download.add_argument("--api-url", default="https://version-12-0.string-db.org")
    analyze_p = sub.add_parser("analyze", parents=[common])
    analyze_p.add_argument("--mapping", type=Path, default=Path("data/network_medicine/processed/joint_string_mapping.tsv"))
    analyze_p.add_argument("--edges", type=Path, default=Path("data/network_medicine/processed/joint_ppi_edges.tsv"))
    analyze_p.add_argument("--results-dir", type=Path, default=Path("data/network_medicine/results/separation"))
    return root


def main() -> int:
    args = parser().parse_args()
    modules = load_modules(args.bc, args.ra, args.dm)
    if args.command == "download":
        if not 0 <= args.score_threshold <= 1:
            raise ValueError("--score-threshold must be between 0 and 1")
        if args.additional_interactors < 0:
            raise ValueError("--additional-interactors must be non-negative")
        download_network(modules, args.raw_dir, args.processed_dir, args.score_threshold,
                         args.network_type, args.additional_interactors, args.api_url)
    else:
        analyze(modules, args.mapping, args.edges, args.results_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
