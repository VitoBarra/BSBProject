"""Build a joint BC/RA/DM STRING network and calculate disease separation."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from NetworkMedicine.bc_network import STRING_VERSION, evidence_score, parse_tsv_text, post_tsv, write_tsv


SPECIES = 9606
CALLER_IDENTITY = "BSBProject_VitoBarra"


def read_symbols(path: Path, column: str, selected_only: bool = False) -> list[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
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


def download_network(
    modules: dict[str, list[str]], raw_dir: Path, processed_dir: Path,
    score_threshold: float, network_type: str, additional_interactors: int,
    api_url: str,
) -> None:
    all_symbols = sorted(set().union(*map(set, modules.values())))
    raw_dir.mkdir(parents=True, exist_ok=True)
    processed_dir.mkdir(parents=True, exist_ok=True)
    common = {
        "identifiers": "\r".join(all_symbols),
        "species": SPECIES,
        "caller_identity": CALLER_IDENTITY,
    }
    mapping_text = post_tsv(api_url, "get_string_ids", {**common, "limit": 1, "echo_query": 1})
    (raw_dir / "string_joint_mapping_raw.tsv").write_text(mapping_text, encoding="utf-8")
    mapping_raw = parse_tsv_text(mapping_text)
    by_query = {row.get("queryItem", ""): row for row in mapping_raw if row.get("queryItem")}

    membership = {symbol: ";".join(name for name, genes in modules.items() if symbol in genes)
                  for symbol in all_symbols}
    mapping = []
    for symbol in all_symbols:
        row = by_query.get(symbol)
        mapping.append({
            "gene_symbol": symbol,
            "disease_modules": membership[symbol],
            "mapping_status": "mapped" if row else "unmapped",
            "string_id": row.get("stringId", "") if row else "",
            "string_preferred_name": row.get("preferredName", "") if row else "",
        })
    write_tsv(processed_dir / "joint_string_mapping.tsv", mapping,
              ["gene_symbol", "disease_modules", "mapping_status", "string_id", "string_preferred_name"])

    mapped_ids = [str(row["string_id"]) for row in mapping if row["string_id"]]
    network_text = post_tsv(api_url, "network", {
        "identifiers": "\r".join(mapped_ids), "species": SPECIES,
        "required_score": 0, "network_type": network_type,
        "add_nodes": additional_interactors,
        "caller_identity": CALLER_IDENTITY,
    })
    (raw_dir / "string_joint_network_raw.tsv").write_text(network_text, encoding="utf-8")
    id_to_symbol = {str(row["string_id"]): str(row["gene_symbol"]) for row in mapping if row["string_id"]}
    edges, seen = [], set()
    for row in parse_tsv_text(network_text):
        a = id_to_symbol.get(row.get("stringId_A", ""), row.get("preferredName_A", ""))
        b = id_to_symbol.get(row.get("stringId_B", ""), row.get("preferredName_B", ""))
        if not a or not b or a == b:
            continue
        class_score = evidence_score(row)
        if class_score < score_threshold:
            continue
        key = tuple(sorted((a, b)))
        if key in seen:
            continue
        seen.add(key)
        edges.append({
            "gene_a": key[0], "gene_b": key[1],
            "combined_score": row.get("score", ""),
            "class_evidence_score": class_score,
            "escore": row.get("escore", ""), "dscore": row.get("dscore", ""),
        })
    edges.sort(key=lambda row: (row["gene_a"], row["gene_b"]))
    write_tsv(processed_dir / "joint_ppi_edges.tsv", edges,
              ["gene_a", "gene_b", "combined_score", "class_evidence_score", "escore", "dscore"])
    metadata = {
        "accessed_at_utc": datetime.now(timezone.utc).isoformat(), "species": SPECIES,
        "api_url": api_url, "string_version": STRING_VERSION,
        "network_type": network_type, "score_threshold": score_threshold,
        "api_required_score": 0,
        "evidence_channels": ["experiments", "databases"],
        "evidence_combination": "STRING probabilistic combination of escore and dscore with prior p=0.041 removed from each channel and restored once",
        "excluded_channels": ["textmining", "coexpression", "neighborhood", "fusion", "cooccurrence"],
        "additional_interactors": additional_interactors, "submitted_unique_genes": len(all_symbols),
        "mapped_genes": len(mapped_ids), "returned_edges": len(edges),
        "module_input_sizes": {name: len(genes) for name, genes in modules.items()},
    }
    (processed_dir / "joint_string_query_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


def all_pairs_distances(graph: nx.Graph, sources: set[str], targets: set[str]) -> dict[str, object]:
    """Class definition: mean shortest path over the Cartesian product."""
    distances: list[float] = []
    unreachable = 0
    for source in sorted(sources):
        lengths = nx.single_source_shortest_path_length(graph, source)
        for target in sorted(targets):
            if target in lengths:
                distances.append(float(lengths[target]))
            else:
                unreachable += 1
    return {
        "distance": sum(distances) / len(distances) if distances else math.nan,
        "finite_pair_distances": len(distances),
        "unreachable_pair_distances": unreachable,
    }


def analyze(modules: dict[str, list[str]], mapping_path: Path, edges_path: Path, results_dir: Path) -> None:
    with mapping_path.open(encoding="utf-8", newline="") as handle:
        mapping = list(csv.DictReader(handle, delimiter="\t"))
    with edges_path.open(encoding="utf-8", newline="") as handle:
        edges = list(csv.DictReader(handle, delimiter="\t"))
    mapped = {row["gene_symbol"] for row in mapping if row["mapping_status"] == "mapped"}
    graph = nx.Graph()
    graph.add_nodes_from(mapped)
    for row in edges:
        score = float(row.get("class_evidence_score") or row["combined_score"])
        graph.add_edge(row["gene_a"], row["gene_b"], score=score)
    projected = {name: set(genes) & mapped for name, genes in modules.items()}

    rows = []
    for disease in ("RA", "DM"):
        bc, other = projected["BC"], projected[disease]
        cross = all_pairs_distances(graph, bc, other)
        within_bc = all_pairs_distances(graph, bc, bc)
        within_other = all_pairs_distances(graph, other, other)
        separation = cross["distance"] - (within_bc["distance"] + within_other["distance"]) / 2
        rows.append({
            "comparison": f"BC-{disease}", "bc_input_genes": len(modules["BC"]),
            "other_input_genes": len(modules[disease]), "bc_mapped_genes": len(bc),
            "other_mapped_genes": len(other), "module_overlap_genes": len(bc & other),
            "jaccard_similarity": len(bc & other) / len(bc | other),
            "d_ab": cross["distance"], "d_aa": within_bc["distance"],
            "d_bb": within_other["distance"], "separation_s_ab": separation,
            "cross_finite_pair_distances": cross["finite_pair_distances"],
            "cross_unreachable_pair_distances": cross["unreachable_pair_distances"],
            "cross_finite_fraction": cross["finite_pair_distances"] / (cross["finite_pair_distances"] + cross["unreachable_pair_distances"]),
            "within_bc_finite_pair_distances": within_bc["finite_pair_distances"],
            "within_bc_unreachable_pair_distances": within_bc["unreachable_pair_distances"],
            "within_other_finite_pair_distances": within_other["finite_pair_distances"],
            "within_other_unreachable_pair_distances": within_other["unreachable_pair_distances"],
            "within_other_finite_fraction": within_other["finite_pair_distances"] / (within_other["finite_pair_distances"] + within_other["unreachable_pair_distances"]),
        })
    results_dir.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    write_tsv(results_dir / "disease_separation.tsv", rows, fields)
    nx.write_graphml(graph, results_dir / "bc_ra_dm_joint_network.graphml")
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
