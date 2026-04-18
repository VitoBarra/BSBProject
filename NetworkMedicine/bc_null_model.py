"""Degree-preserving null model for BC network community structure."""

from __future__ import annotations

import networkx as nx


def modularity_null_model(
    graph: nx.Graph, membership: dict[str, int], seed: int, iterations: int,
) -> tuple[float, list[set[str]], list[float]]:
    """Compare the observed communities with degree-preserving random graphs."""
    non_isolate_graph = graph.subgraph([node for node in graph if graph.degree(node) > 0])
    non_isolate_membership = {node: membership[node] for node in non_isolate_graph}
    community_sets = [
        {node for node, value in non_isolate_membership.items() if value == community_id}
        for community_id in sorted(set(non_isolate_membership.values()))
    ]
    modularity = (
        nx.community.modularity(non_isolate_graph, community_sets, weight="score")
        if non_isolate_graph.number_of_edges() and len(community_sets) > 1
        else 0.0
    )
    null_values: list[float] = []
    if non_isolate_graph.number_of_edges() >= 2 and iterations:
        base = nx.Graph(non_isolate_graph)
        original_weights = [float(data["score"]) for _, _, data in base.edges(data=True)]
        for iteration in range(iterations):
            randomized = base.copy()
            try:
                nx.double_edge_swap(
                    randomized,
                    nswap=max(1, 5 * randomized.number_of_edges()),
                    max_tries=max(100, 100 * randomized.number_of_edges()),
                    seed=seed + iteration,
                )
            except nx.NetworkXAlgorithmError:
                continue
            for edge, weight in zip(randomized.edges(), original_weights, strict=True):
                randomized.edges[edge]["score"] = weight
            null_communities = nx.community.louvain_communities(
                randomized, weight="score", seed=seed + iteration
            )
            null_values.append(nx.community.modularity(randomized, null_communities, weight="score"))
    return modularity, community_sets, null_values
