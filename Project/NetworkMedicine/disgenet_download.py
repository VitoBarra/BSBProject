"""Download curated RA and DM gene–disease associations from DISGENET.

The Academic DISGENET plan permits API access to CURATED associations. The API
key is read exclusively from the DISGENET_API_KEY environment variable.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


API_URL = "http://api.disgenet.com/api/v1/gda/summary"
PAGE_SIZE = 100
DISEASES = {
    "ra": ("UMLS_C0003873", "Rheumatoid Arthritis"),
    "dm": ("UMLS_C0011849", "Diabetes Mellitus"),
}


def _request_json(url: str, api_key: str, retries: int = 4) -> dict[str, Any]:
    request = urllib.request.Request(
        url,
        headers={
            "Authorization": api_key,
            "Accept": "application/json",
            "User-Agent": "BSBProject-DISGENET-downloader/1.0",
        },
    )
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                return json.load(response)
        except urllib.error.HTTPError as exc:
            if exc.code == 429 and attempt + 1 < retries:
                retry_after = exc.headers.get("x-rate-limit-retry-after-seconds", "5")
                time.sleep(max(float(retry_after), 1.0))
                continue
            body = exc.read().decode("utf-8", errors="replace")
            raise RuntimeError(
                f"DISGENET returned HTTP {exc.code}: {body[:500]}"
            ) from exc
        except urllib.error.URLError as exc:
            if attempt + 1 == retries:
                raise RuntimeError(f"Could not contact DISGENET: {exc.reason}") from exc
            time.sleep(2**attempt)
    raise RuntimeError("DISGENET request failed after retries")


def _scalar(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    return str(value)


def _first(row: dict[str, Any], *names: str) -> Any:
    for name in names:
        if name in row and row[name] not in (None, ""):
            return row[name]
    return ""


def _number(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("-inf")


def fetch_disease(
    disease_id: str, api_key: str, max_pages: int
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    base_params = {
        "disease": disease_id,
        "source": "CURATED",
        "min_score": "0",
        # DISGENET scores can exceed 1 in current releases.
        "max_score": "2",
        "order_by": "score",
    }
    rows_by_key: dict[str, dict[str, Any]] = {}
    first_response: dict[str, Any] | None = None
    total_elements = 0
    pages_requested = 0

    for page in range(max_pages):
        params = dict(base_params, page_number=str(page))
        url = f"{API_URL}?{urllib.parse.urlencode(params)}"
        response = _request_json(url, api_key)
        if first_response is None:
            first_response = response
            total_elements = int(response.get("paging", {}).get("totalElements", 0))
            profile = str(response.get("userinfo", {}).get("profile", "unknown"))
            first_payload = response.get("payload") or []
            if (
                profile.upper() == "TRIAL"
                and isinstance(first_payload, list)
                and total_elements > len(first_payload)
            ):
                raise RuntimeError(
                    "DISGENET authenticated the API key as a TRIAL profile. "
                    f"The query has {total_elements} results, but TRIAL access "
                    f"returns only {len(first_payload)} and forbids pagination. "
                    "Wait for/verify Academic-plan activation in the DISGENET "
                    "profile, then generate or copy the active profile's API key."
                )

        payload = response.get("payload") or []
        if not isinstance(payload, list):
            raise RuntimeError("Unexpected DISGENET response: payload is not a list")

        for row in payload:
            if not isinstance(row, dict):
                continue
            key = str(
                _first(row, "assocID", "associationId", "associationID")
                or json.dumps(row, sort_keys=True, default=str)
            )
            rows_by_key[key] = row

        pages_requested += 1
        if not payload or len(rows_by_key) >= total_elements:
            break
        time.sleep(1.5)
    else:
        raise RuntimeError(
            f"Query requires more than {max_pages} pages; refine the query."
        )

    profile = (first_response or {}).get("userinfo", {}).get("profile", "unknown")
    metadata = {
        "api_endpoint": API_URL,
        "database": "CURATED",
        "disease_query": disease_id,
        "downloaded_at_utc": datetime.now(timezone.utc).isoformat(),
        "max_score": 2,
        "min_score": 0,
        "pages_requested": pages_requested,
        "records_downloaded": len(rows_by_key),
        "total_elements_reported": total_elements,
        "user_profile": profile,
    }
    return list(rows_by_key.values()), metadata


def write_associations(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = sorted({field for row in rows for field in row})
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: _scalar(row.get(field)) for field in fields})


def write_genes(
    path: Path,
    rows: list[dict[str, Any]],
    disease_id: str,
    disease_name: str,
) -> int:
    genes: dict[str, dict[str, Any]] = {}
    for row in rows:
        symbol = str(
            _first(
                row,
                "geneSymbol_keywrod",
                "geneSymbol",
                "symbolOfGene",
                "gene_symbol",
            )
        ).strip()
        if not symbol:
            continue
        candidate = {
            "disease_id": disease_id.removeprefix("UMLS_"),
            "disease_name": disease_name,
            "gene_symbol": symbol,
            "gene_id": _first(row, "geneNcbiID", "geneID", "geneid"),
            "normalized_score": _first(
                row, "normalizedScore", "normalized_score", "scoreNormalized"
            ),
            "score": _first(row, "score"),
        }
        previous = genes.get(symbol)
        if previous is None or (
            _number(candidate["normalized_score"]),
            _number(candidate["score"]),
        ) > (
            _number(previous["normalized_score"]),
            _number(previous["score"]),
        ):
            genes[symbol] = candidate

    ordered = sorted(
        genes.values(),
        key=lambda row: (
            -_number(row["normalized_score"]),
            -_number(row["score"]),
            row["gene_symbol"],
        ),
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        fields = [
            "disease_id",
            "disease_name",
            "gene_symbol",
            "gene_id",
            "normalized_score",
            "score",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(ordered)
    return len(ordered)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("data/network_medicine"),
    )
    parser.add_argument("--max-pages", type=int, default=100)
    args = parser.parse_args()

    api_key = os.environ.get("DISGENET_API_KEY", "").strip()
    if not api_key or api_key == "replace_with_your_disgenet_api_key":
        raise SystemExit(
            "DISGENET_API_KEY is missing. Add it to the repository .env.local file."
        )

    for short_name, (disease_id, disease_name) in DISEASES.items():
        print(f"Downloading {disease_name} ({disease_id}) from CURATED...")
        rows, metadata = fetch_disease(disease_id, api_key, args.max_pages)
        raw_path = (
            args.output_root
            / "raw"
            / f"disgenet_{short_name}_associations_curated.tsv"
        )
        metadata_path = (
            args.output_root / "raw" / f"disgenet_{short_name}_metadata.json"
        )
        gene_path = (
            args.output_root
            / "processed"
            / f"{short_name}_genes_disgenet_curated.tsv"
        )
        write_associations(raw_path, rows)
        metadata_path.write_text(
            json.dumps(metadata, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        gene_count = write_genes(
            gene_path, rows, disease_id=disease_id, disease_name=disease_name
        )
        print(
            f"  {len(rows)} associations; {gene_count} unique genes -> {gene_path}"
        )


if __name__ == "__main__":
    main()
