"""Build the BC graph, compute network metrics, and select the disease module."""

from __future__ import annotations

import math
import statistics
from pathlib import Path

import networkx as nx
from utils.tsv import read_tsv


BIOLOGICAL_RATIONALE = {
    "MSH2": "DNA mismatch-repair gene; selected as a hub of the four-gene repair community.",
    "MSH6": "MSH2 partner in mismatch repair; tied for the highest composite centrality.",
    "RUNX1": "Transcriptional regulator with cancer-gene annotation; hub of its three-gene community.",
    "TAL1": "Transcription factor connected to RUNX1/FLI1 and highly ranked within that community.",
    "CDKN2A": "Canonical cell-cycle tumour suppressor and representative of the CDKN2A-MDM2 component.",
    "CHD4": "Chromatin-remodelling factor and representative of the CHD4-H3F3A component.",
    "EGFR": "Receptor tyrosine kinase with established relevance to breast-cancer signalling.",
    "GNAS": "Signal-transduction gene representing the GNAS-TSHR component.",
    "IL6ST": "Shared cytokine-receptor subunit representing the IL6ST-LIFR component.",
    "MTCP1": "Cancer-gene representative of the MTCP1-TCL1A component.",
}

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
    mapping = read_tsv(mapping_path)
    edges = read_tsv(edges_path)

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
        score = float(row.get("class_evidence_score") or row["combined_score"])
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
    # Represent larger connected communities first, then fill remaining slots
    # with the highest-ranked genes. Isolated genes are not eligible.
    for community_id, community_rows in sorted(
        by_community.items(),
        key=lambda item: (-len(item[1]), item[0]),
    ):
        if len(selected) >= target_size:
            break
        gene = str(community_rows[0]["gene"])
        selected.add(gene)
        reasons[gene] = (
            f"top composite-ranked representative of nontrivial community {community_id}"
        )

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


def rank_and_select_genes(
    graph: nx.Graph, mapping_by_gene: dict[str, dict[str, str]], target_size: int, seed: int,
) -> tuple[list[dict[str, object]], dict[str, int], list[set[str]], set[str]]:
    """Calculate centralities and communities, then select the BC module."""
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

    # Percentile ranks put the four metrics on a common scale; their mean gives
    # each metric equal weight in the final ranking.
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
                "biological_rationale": "",
            }
        )

    selected, reasons = select_module(rows, target_size)
    for row in rows:
        gene = str(row["gene"])
        row["selected"] = gene in selected
        row["selection_reason"] = reasons.get(gene, "")
        row["biological_rationale"] = BIOLOGICAL_RATIONALE.get(gene, "") if gene in selected else ""
    rows.sort(key=lambda row: (-float(row["composite_percentile"]), str(row["gene"])))
    return rows, membership, components, selected
