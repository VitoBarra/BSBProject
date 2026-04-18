"""Download and prepare the breast-cancer STRING network."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from utils.tsv import read_tsv, write_tsv
from NetworkMedicine.string_api import (CALLER_IDENTITY, SPECIES, STRING_VERSION, filtered_edges, post_tsv, request_mappings)


def read_gene_table(path: Path) -> list[dict[str, str]]:
    rows = read_tsv(path)
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

    mapping_text, mapping_rows = request_mappings(api_url, symbols)
    (raw_dir / "string_mapping_raw.tsv").write_text(mapping_text, encoding="utf-8")

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
            # Apply the threshold locally to experiments + databases only.
            "required_score": 0,
            "network_type": network_type,
            "caller_identity": CALLER_IDENTITY,
        },
    )
    (raw_dir / "string_network_raw.tsv").write_text(network_text, encoding="utf-8")
    id_to_input = {
        str(row["string_id"]): str(row["gene_symbol"])
        for row in clean_mapping
        if row["string_id"]
    }
    clean_edges: list[dict[str, object]] = []
    for edge_key, class_score, row, gene_a, gene_b in filtered_edges(network_text, id_to_input, score_threshold):
        id_a = row.get("stringId_A", "")
        id_b = row.get("stringId_B", "")
        if gene_a not in symbols or gene_b not in symbols:
            continue
        clean_edges.append(
            {
                "gene_a": edge_key[0],
                "gene_b": edge_key[1],
                "string_id_a": id_a if gene_a == edge_key[0] else id_b,
                "string_id_b": id_b if gene_b == edge_key[1] else id_a,
                "combined_score": row.get("score", ""),
                "class_evidence_score": class_score,
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
            "class_evidence_score",
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
        "string_version": STRING_VERSION,
        "api_method_mapping": "get_string_ids",
        "api_method_network": "network",
        "species": SPECIES,
        "network_type": network_type,
        "score_threshold": score_threshold,
        "api_required_score": 0,
        "local_required_score": required_score,
        "evidence_channels": ["experiments", "databases"],
        "evidence_combination": "STRING probabilistic combination of escore and dscore with prior p=0.041 removed from each channel and restored once",
        "excluded_channels": ["textmining", "coexpression", "neighborhood", "fusion", "cooccurrence"],
        "additional_interactors": 0,
        "input_gene_count": len(genes),
        "mapped_gene_count": len(mapped_ids),
        "unmapped_gene_count": len(genes) - len(mapped_ids),
        "returned_edge_count": len(clean_edges),
        "caller_identity": CALLER_IDENTITY,
        "interpretation": (
            "Physical STRING associations among the submitted genes only, "
            "filtered locally to experiments/database evidence score >= threshold. "
            "No additional interactors were requested."
        ),
    }
    (processed_dir / "bc_string_query_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n",
        encoding="utf-8",
    )
