from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import urllib.parse
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx


DEFAULT_GENES = Path("data/network_medicine/processed/bc_genes_s11.tsv")
DEFAULT_RAW_DIR = Path("data/network_medicine/raw/string_bc")
DEFAULT_RESULTS_DIR = Path("data/network_medicine/results/bc")
DEFAULT_API_URL = "https://string-db.org"
SPECIES = 9606
CALLER_IDENTITY = "BSBProject_VitoBarra"


def read_gene_table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"gene_symbol", "ensembl_gene_id"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"{path} must contain columns: {sorted(required)}")
    for row in rows:
        row["gene_symbol"] = row["gene_symbol"].strip()
        row["ensembl_gene_id"] = row["ensembl_gene_id"].strip()
    symbols = [row["gene_symbol"] for row in rows]
    if any(not symbol for symbol in symbols):
        raise ValueError("Gene symbols must not be empty")
    if len(symbols) != len(set(symbols)):
        raise ValueError("Gene symbols must be unique")
    return rows


def post_tsv(api_url: str, method: str, parameters: dict[str, str | int]) -> str:
    url = f"{api_url.rstrip('/')}/api/tsv/{method}"
    payload = urllib.parse.urlencode(parameters).encode("utf-8")
    request = urllib.request.Request(
        url,
        data=payload,
        headers={"User-Agent": f"{CALLER_IDENTITY}/1.0"},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        return response.read().decode("utf-8")


def parse_tsv_text(text: str) -> list[dict[str, str]]:
    if not text.strip():
        return []
    return list(csv.DictReader(text.splitlines(), delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def download_string_data(
    genes_path: Path,
    raw_dir: Path,
    processed_dir: Path,
    score_threshold: float,
    network_type: str,
    api_url: str,
) -> None:
    genes = read_gene_table(genes_path)
    symbols = [row["gene_symbol"] for row in genes]
    required_score = round(score_threshold * 1000)
    raw_dir.mkdir(parents=True, exist_ok=True)
    processed_dir.mkdir(parents=True, exist_ok=True)

    common_parameters: dict[str, str | int] = {
        "identifiers": "\r".join(symbols),
        "species": SPECIES,
        "caller_identity": CALLER_IDENTITY,
    }
    mapping_text = post_tsv(
        api_url,
        "get_string_ids",
        {**common_parameters, "limit": 1, "echo_query": 1},
    )
    (raw_dir / "string_mapping_raw.tsv").write_text(mapping_text, encoding="utf-8")
    mapping_rows = parse_tsv_text(mapping_text)

    mapping_by_query: dict[str, dict[str, str]] = {}
    for row in mapping_rows:
        query = row.get("queryItem", "")
        if query and query not in mapping_by_query:
            mapping_by_query[query] = row

    clean_mapping: list[dict[str, object]] = []
    for gene in genes:
        symbol = gene["gene_symbol"]
        mapped = mapping_by_query.get(symbol)
        clean_mapping.append(
            {
                "gene_symbol": symbol,
                "ensembl_gene_id": gene["ensembl_gene_id"],
                "mapping_status": "mapped" if mapped else "unmapped",
                "string_id": mapped.get("stringId", "") if mapped else "",
                "string_preferred_name": mapped.get("preferredName", "") if mapped else "",
                "string_annotation": mapped.get("annotation", "") if mapped else "",
            }
        )
    write_tsv(
        processed_dir / "bc_string_mapping.tsv",
        clean_mapping,
        [
            "gene_symbol",
            "ensembl_gene_id",
            "mapping_status",
            "string_id",
            "string_preferred_name",
            "string_annotation",
        ],
    )

    mapped_ids = [str(row["string_id"]) for row in clean_mapping if row["string_id"]]
    if len(mapped_ids) < 2:
        raise RuntimeError("Fewer than two genes mapped to STRING")

    network_text = post_tsv(
        api_url,
        "network",
        {
            "identifiers": "\r".join(mapped_ids),
            "species": SPECIES,
            "required_score": required_score,
            "network_type": network_type,
            "caller_identity": CALLER_IDENTITY,
        },
    )
    (raw_dir / "string_network_raw.tsv").write_text(network_text, encoding="utf-8")
    network_rows = parse_tsv_text(network_text)

    id_to_input = {
        str(row["string_id"]): str(row["gene_symbol"])
        for row in clean_mapping
        if row["string_id"]
    }
    clean_edges: list[dict[str, object]] = []
    seen: set[tuple[str, str]] = set()
    for row in network_rows:
        id_a = row.get("stringId_A", "")
        id_b = row.get("stringId_B", "")
        gene_a = id_to_input.get(id_a, row.get("preferredName_A", ""))
        gene_b = id_to_input.get(id_b, row.get("preferredName_B", ""))
        if gene_a not in symbols or gene_b not in symbols or gene_a == gene_b:
            continue
        edge_key = tuple(sorted((gene_a, gene_b)))
        if edge_key in seen:
            continue
        seen.add(edge_key)
        clean_edges.append(
            {
                "gene_a": edge_key[0],
                "gene_b": edge_key[1],
                "string_id_a": id_a if gene_a == edge_key[0] else id_b,
                "string_id_b": id_b if gene_b == edge_key[1] else id_a,
                "combined_score": row.get("score", ""),
                "nscore": row.get("nscore", ""),
                "fscore": row.get("fscore", ""),
                "pscore": row.get("pscore", ""),
                "ascore": row.get("ascore", ""),
                "escore": row.get("escore", ""),
                "dscore": row.get("dscore", ""),
                "tscore": row.get("tscore", ""),
            }
        )
    clean_edges.sort(key=lambda row: (str(row["gene_a"]), str(row["gene_b"])))
    write_tsv(
        processed_dir / "bc_ppi_edges.tsv",
        clean_edges,
        [
            "gene_a",
            "gene_b",
            "string_id_a",
            "string_id_b",
            "combined_score",
            "nscore",
            "fscore",
            "pscore",
            "ascore",
            "escore",
            "dscore",
            "tscore",
        ],
    )

    metadata = {
        "accessed_at_utc": datetime.now(timezone.utc).isoformat(),
        "api_url": api_url,
        "api_method_mapping": "get_string_ids",
        "api_method_network": "network",
        "species": SPECIES,
        "network_type": network_type,
        "score_threshold": score_threshold,
        "required_score": required_score,
        "additional_interactors": 0,
        "input_gene_count": len(genes),
        "mapped_gene_count": len(mapped_ids),
        "unmapped_gene_count": len(genes) - len(mapped_ids),
        "returned_edge_count": len(clean_edges),
        "caller_identity": CALLER_IDENTITY,
        "interpretation": (
            "Functional STRING associations among the submitted genes only. "
            "No additional interactors were requested."
        ),
    }
    (processed_dir / "bc_string_query_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n",
        encoding="utf-8",
    )


def percentile_ranks(values: dict[str, float]) -> dict[str, float]:
    ordered = sorted(values.items(), key=lambda item: (item[1], item[0]))
    n = len(ordered)
    if n <= 1:
        return {node: 1.0 for node, _ in ordered}
    ranks: dict[str, float] = {}
    start = 0
    while start < n:
        end = start
        while end + 1 < n and ordered[end + 1][1] == ordered[start][1]:
            end += 1
        average_rank = (start + end) / 2
        percentile = average_rank / (n - 1)
        for index in range(start, end + 1):
            ranks[ordered[index][0]] = percentile
        start = end + 1
    return ranks


def build_graph(mapping_path: Path, edges_path: Path) -> tuple[nx.Graph, dict[str, dict[str, str]]]:
    with mapping_path.open(encoding="utf-8", newline="") as handle:
        mapping = list(csv.DictReader(handle, delimiter="\t"))
    with edges_path.open(encoding="utf-8", newline="") as handle:
        edges = list(csv.DictReader(handle, delimiter="\t"))

    mapping_by_gene = {row["gene_symbol"]: row for row in mapping}
    graph = nx.Graph()
    for row in mapping:
        graph.add_node(
            row["gene_symbol"],
            mapping_status=row["mapping_status"],
            string_id=row["string_id"],
            string_preferred_name=row["string_preferred_name"],
        )
    for row in edges:
        score = float(row["combined_score"])
        graph.add_edge(
            row["gene_a"],
            row["gene_b"],
            score=score,
            distance=1.0 / score if score > 0 else math.inf,
        )
    return graph, mapping_by_gene


def detect_communities(graph: nx.Graph, seed: int) -> dict[str, int]:
    membership: dict[str, int] = {}
    connected_nodes = [node for node in graph if graph.degree(node) > 0]
    if connected_nodes:
        connected_graph = graph.subgraph(connected_nodes)
        communities = nx.community.louvain_communities(
            connected_graph,
            weight="score",
            seed=seed,
        )
        ordered = sorted(communities, key=lambda members: (-len(members), sorted(members)[0]))
        for community_id, members in enumerate(ordered, start=1):
            for node in members:
                membership[node] = community_id
    next_id = max(membership.values(), default=0) + 1
    for node in sorted(graph):
        if node not in membership:
            membership[node] = next_id
            next_id += 1
    return membership


def weighted_pagerank(
    graph: nx.Graph,
    alpha: float = 0.85,
    tolerance: float = 1.0e-12,
    max_iterations: int = 1000,
) -> dict[str, float]:
    """Compute weighted PageRank without requiring NetworkX's SciPy backend."""
    nodes = list(graph)
    node_count = len(nodes)
    if node_count == 0:
        return {}
    rank = {node: 1.0 / node_count for node in nodes}
    strengths = {
        node: sum(float(data.get("score", 1.0)) for _, _, data in graph.edges(node, data=True))
        for node in nodes
    }
    teleport = (1.0 - alpha) / node_count
    for _ in range(max_iterations):
        dangling_mass = alpha * sum(
            rank[node] for node in nodes if strengths[node] == 0
        ) / node_count
        updated = {node: teleport + dangling_mass for node in nodes}
        for source in nodes:
            if strengths[source] == 0:
                continue
            contribution = alpha * rank[source] / strengths[source]
            for target, edge_data in graph[source].items():
                updated[target] += contribution * float(edge_data.get("score", 1.0))
        difference = sum(abs(updated[node] - rank[node]) for node in nodes)
        rank = updated
        if difference < tolerance:
            total = sum(rank.values())
            return {node: value / total for node, value in rank.items()}
    raise RuntimeError(f"Weighted PageRank did not converge after {max_iterations} iterations")


def select_module(
    rows: list[dict[str, object]],
    target_size: int,
) -> tuple[set[str], dict[str, str]]:
    non_isolates = [row for row in rows if int(row["degree"]) > 0]
    target_size = min(target_size, len(non_isolates))
    by_community: dict[int, list[dict[str, object]]] = {}
    for row in non_isolates:
        by_community.setdefault(int(row["community"]), []).append(row)
    for community_rows in by_community.values():
        community_rows.sort(
            key=lambda row: (-float(row["composite_percentile"]), str(row["gene"]))
        )

    selected: set[str] = set()
    reasons: dict[str, str] = {}
    for community_id, community_rows in sorted(
        by_community.items(),
        key=lambda item: (-len(item[1]), item[0]),
    ):
        if len(selected) >= target_size:
            break
        if len(community_rows) < 3:
            continue
        gene = str(community_rows[0]["gene"])
        selected.add(gene)
        reasons[gene] = f"top composite-ranked representative of community {community_id}"

    ranked = sorted(
        non_isolates,
        key=lambda row: (
            -float(row["composite_percentile"]),
            -float(row["betweenness_centrality"]),
            str(row["gene"]),
        ),
    )
    for row in ranked:
        if len(selected) >= target_size:
            break
        gene = str(row["gene"])
        if gene not in selected:
            selected.add(gene)
            reasons[gene] = "high composite centrality rank"
    return selected, reasons


def plot_network(graph: nx.Graph, rows_by_gene: dict[str, dict[str, object]], output: Path, seed: int) -> None:
    plt.figure(figsize=(14, 11))
    position = nx.spring_layout(graph, seed=seed, weight="score", k=1.2 / math.sqrt(max(len(graph), 1)))
    node_sizes = [
        250 + 2500 * float(rows_by_gene[node]["pagerank"])
        for node in graph
    ]
    node_colors = [int(rows_by_gene[node]["community"]) for node in graph]
    edge_widths = [0.4 + 2.2 * float(data["score"]) for _, _, data in graph.edges(data=True)]
    nx.draw_networkx_edges(graph, position, width=edge_widths, alpha=0.45, edge_color="#64748b")
    nx.draw_networkx_nodes(
        graph,
        position,
        node_size=node_sizes,
        node_color=node_colors,
        cmap=plt.colormaps["tab20"],
        edgecolors="#172033",
        linewidths=0.7,
    )
    labels = {
        node: node
        for node in graph
        if bool(rows_by_gene[node]["selected"]) or graph.degree(node) == 0
    }
    nx.draw_networkx_labels(graph, position, labels=labels, font_size=8)
    plt.title("BC Table S11 induced STRING network")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(output, dpi=220, bbox_inches="tight")
    plt.close()


def plot_largest_component(
    graph: nx.Graph,
    rows_by_gene: dict[str, dict[str, object]],
    output: Path,
    seed: int,
) -> None:
    largest_nodes = max(nx.connected_components(graph), key=len)
    core = graph.subgraph(largest_nodes)
    plt.figure(figsize=(15, 12))
    position = nx.spring_layout(core, seed=seed, weight="score", k=0.75, iterations=300)
    node_sizes = [
        500 + 5500 * float(rows_by_gene[node]["pagerank"])
        for node in core
    ]
    node_colors = [int(rows_by_gene[node]["community"]) for node in core]
    edge_widths = [0.5 + 2.5 * float(data["score"]) for _, _, data in core.edges(data=True)]
    nx.draw_networkx_edges(core, position, width=edge_widths, alpha=0.5, edge_color="#64748b")
    nx.draw_networkx_nodes(
        core,
        position,
        node_size=node_sizes,
        node_color=node_colors,
        cmap=plt.colormaps["tab10"],
        edgecolors=[
            "#d62728" if bool(rows_by_gene[node]["selected"]) else "#172033"
            for node in core
        ],
        linewidths=[
            2.5 if bool(rows_by_gene[node]["selected"]) else 0.8
            for node in core
        ],
    )
    nx.draw_networkx_labels(core, position, font_size=9)
    plt.title("Largest connected component of the BC STRING network\n(red border = preliminary selection)")
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


def analyze_network(
    mapping_path: Path,
    edges_path: Path,
    results_dir: Path,
    target_size: int,
    seed: int,
) -> None:
    graph, mapping_by_gene = build_graph(mapping_path, edges_path)
    results_dir.mkdir(parents=True, exist_ok=True)

    degree = dict(graph.degree())
    degree_centrality = nx.degree_centrality(graph)
    betweenness = nx.betweenness_centrality(graph, normalized=True, weight=None)
    harmonic_raw = nx.harmonic_centrality(graph)
    harmonic = {
        node: value / (len(graph) - 1) if len(graph) > 1 else 0.0
        for node, value in harmonic_raw.items()
    }
    pagerank = weighted_pagerank(graph, alpha=0.85)
    membership = detect_communities(graph, seed)
    components = sorted(nx.connected_components(graph), key=lambda nodes: (-len(nodes), sorted(nodes)[0]))
    component_id: dict[str, int] = {}
    component_size: dict[str, int] = {}
    for identifier, component in enumerate(components, start=1):
        for node in component:
            component_id[node] = identifier
            component_size[node] = len(component)

    metric_percentiles = {
        "degree": percentile_ranks({node: float(value) for node, value in degree.items()}),
        "betweenness": percentile_ranks(betweenness),
        "harmonic": percentile_ranks(harmonic),
        "pagerank": percentile_ranks(pagerank),
    }
    composite = {
        node: statistics.mean(metric[node] for metric in metric_percentiles.values())
        for node in graph
    }
    rows: list[dict[str, object]] = []
    for node in graph:
        rows.append(
            {
                "gene": node,
                "ensembl_gene_id": mapping_by_gene[node]["ensembl_gene_id"],
                "mapping_status": mapping_by_gene[node]["mapping_status"],
                "string_id": mapping_by_gene[node]["string_id"],
                "degree": degree[node],
                "degree_centrality": degree_centrality[node],
                "betweenness_centrality": betweenness[node],
                "harmonic_centrality": harmonic[node],
                "pagerank": pagerank[node],
                "degree_percentile": metric_percentiles["degree"][node],
                "betweenness_percentile": metric_percentiles["betweenness"][node],
                "harmonic_percentile": metric_percentiles["harmonic"][node],
                "pagerank_percentile": metric_percentiles["pagerank"][node],
                "composite_percentile": composite[node],
                "community": membership[node],
                "component": component_id[node],
                "component_size": component_size[node],
                "selected": False,
                "selection_reason": "",
            }
        )

    selected, reasons = select_module(rows, target_size)
    for row in rows:
        gene = str(row["gene"])
        row["selected"] = gene in selected
        row["selection_reason"] = reasons.get(gene, "")
    rows.sort(key=lambda row: (-float(row["composite_percentile"]), str(row["gene"])))

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
    ]
    write_tsv(results_dir / "bc_gene_centralities.tsv", rows, fields)
    write_tsv(
        results_dir / "bc_module_selected_preliminary.tsv",
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

    non_isolate_graph = graph.subgraph([node for node in graph if graph.degree(node) > 0])
    non_isolate_membership = {
        node: membership[node] for node in non_isolate_graph
    }
    community_sets: list[set[str]] = []
    for community_id in sorted(set(non_isolate_membership.values())):
        community_sets.append(
            {node for node, value in non_isolate_membership.items() if value == community_id}
        )
    modularity = (
        nx.community.modularity(non_isolate_graph, community_sets, weight="score")
        if non_isolate_graph.number_of_edges() and len(community_sets) > 1
        else 0.0
    )
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
        "selected_module_size": len(selected),
        "selected_genes": [str(row["gene"]) for row in rows if bool(row["selected"])],
        "random_seed": seed,
        "centrality_policy": {
            "degree": "unweighted",
            "betweenness": "unweighted shortest paths",
            "harmonic": "unweighted and divided by N-1",
            "pagerank": "weighted by STRING combined score",
            "composite": "mean of the four within-network percentile ranks",
        },
        "selection_policy": (
            "Top composite-ranked representative of every Louvain community "
            "with at least three genes, then highest composite ranks until target size."
        ),
        "selection_status": (
            "Preliminary network-based selection; biological review is required "
            "before treating it as the final BC module."
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
    download.add_argument("--score-threshold", type=float, default=0.80)
    download.add_argument("--network-type", choices=("functional", "physical"), default="functional")
    download.add_argument("--api-url", default=DEFAULT_API_URL)

    analyze = subparsers.add_parser("analyze", help="Calculate centralities, communities, and a preliminary module.")
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
    analyze.add_argument("--target-size", type=int, default=20)
    analyze.add_argument("--seed", type=int, default=42)
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
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
