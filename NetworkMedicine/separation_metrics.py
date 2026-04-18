"""Construct the joint graph and calculate disease-module separation."""

from __future__ import annotations

import math
from pathlib import Path

import networkx as nx
from utils.tsv import read_tsv


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
    # Disconnected pairs cannot contribute a finite distance; report their
    # count alongside the mean so its path coverage remains visible.
    return {
        "distance": sum(distances) / len(distances) if distances else math.nan,
        "finite_pair_distances": len(distances),
        "unreachable_pair_distances": unreachable,
    }


def build_joint_graph(mapping_path: Path, edges_path: Path) -> tuple[nx.Graph, set[str]]:
    """Build the full network and retain mapped disease genes, including isolates."""
    mapping = read_tsv(mapping_path)
    edges = read_tsv(edges_path)
    mapped = {row["gene_symbol"] for row in mapping if row["mapping_status"] == "mapped"}
    graph = nx.Graph()
    graph.add_nodes_from(mapped)
    for row in edges:
        score = float(row.get("class_evidence_score") or row["combined_score"])
        graph.add_edge(row["gene_a"], row["gene_b"], score=score)
    return graph, mapped


def calculate_separations(
    graph: nx.Graph, modules: dict[str, list[str]], projected: dict[str, set[str]]
) -> list[dict[str, object]]:
    """Calculate module distances on the full graph, retaining path coverage."""
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
    return rows
