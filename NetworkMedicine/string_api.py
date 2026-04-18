"""Shared STRING API access and evidence-filtering utilities."""

from __future__ import annotations

import csv
import math
import urllib.parse
import urllib.request


SPECIES = 9606
CALLER_IDENTITY = "BSBProject_VitoBarra"
STRING_PRIOR = 0.041
STRING_VERSION = "12.0"


def evidence_score(row: dict[str, str]) -> float:
    """Combine STRING experiments and databases channel scores.

    The API does not expose the web UI's evidence-source checkboxes. Combining
    only these channels makes the class-filtering policy explicit and excludes
    text-mining-only evidence.
    """
    channels = [
        float(row.get("escore", "") or 0.0),
        float(row.get("dscore", "") or 0.0),
    ]
    adjusted = [
        max(0.0, (score - STRING_PRIOR) / (1.0 - STRING_PRIOR))
        for score in channels
        if score > 0.0
    ]
    if not adjusted:
        return 0.0
    combined_without_prior = 1.0 - math.prod(1.0 - score for score in adjusted)
    return combined_without_prior + STRING_PRIOR * (1.0 - combined_without_prior)


def post_tsv(api_url: str, method: str, parameters: dict[str, str | int]) -> str:
    """POST to a STRING TSV endpoint and return its UTF-8 response."""
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
    """Parse a STRING TSV response, accepting an empty response."""
    if not text.strip():
        return []
    return list(csv.DictReader(text.splitlines(), delimiter="\t"))


def request_mappings(api_url: str, symbols: list[str]) -> tuple[str, list[dict[str, str]]]:
    """Request STRING mappings and return both the raw and parsed response."""
    text = post_tsv(api_url, "get_string_ids", {
        "identifiers": "\r".join(symbols),
        "species": SPECIES,
        "caller_identity": CALLER_IDENTITY,
        "limit": 1,
        "echo_query": 1,
    })
    return text, parse_tsv_text(text)


def filtered_edges(text: str, id_to_symbol: dict[str, str], threshold: float):
    """Yield unique STRING edges passing the shared evidence policy."""
    seen: set[tuple[str, str]] = set()
    for row in parse_tsv_text(text):
        a = id_to_symbol.get(row.get("stringId_A", ""), row.get("preferredName_A", ""))
        b = id_to_symbol.get(row.get("stringId_B", ""), row.get("preferredName_B", ""))
        if not a or not b or a == b:
            continue
        score = evidence_score(row)
        key = tuple(sorted((a, b)))
        if score < threshold or key in seen:
            continue
        seen.add(key)
        yield key, score, row, a, b
