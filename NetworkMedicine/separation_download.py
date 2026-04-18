"""Download the joint STRING network for the three disease modules."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from utils.tsv import write_tsv
from NetworkMedicine.string_api import (CALLER_IDENTITY, SPECIES, STRING_VERSION, filtered_edges, post_tsv, request_mappings)


def download_network(
    modules: dict[str, list[str]], raw_dir: Path, processed_dir: Path,
    score_threshold: float, network_type: str, additional_interactors: int,
    api_url: str,
) -> None:
    all_symbols = sorted(set().union(*map(set, modules.values())))
    raw_dir.mkdir(parents=True, exist_ok=True)
    processed_dir.mkdir(parents=True, exist_ok=True)
    mapping_text, mapping_raw = request_mappings(api_url, all_symbols)
    (raw_dir / "string_joint_mapping_raw.tsv").write_text(mapping_text, encoding="utf-8")
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
    edges = []
    for key, class_score, row, _, _ in filtered_edges(network_text, id_to_symbol, score_threshold):
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
