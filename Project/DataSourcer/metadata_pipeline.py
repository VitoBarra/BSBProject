from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import quote_plus
from urllib.request import urlopen

from Log.log_util import log

EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id={accession}&rettype=runinfo&retmode=text"
ENA_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport?accession={accession}&result=read_run&fields=run_accession,fastq_ftp,fastq_bytes"
NORMALIZED_FIELDS = [
    "sample_id",
    "title",
    "condition",
    "gsm",
    "source_name",
    "platform_id",
    "srx",
    "srr",
    "fastq_url",
    "sample_total_bytes",
    "sample_total_gb",
]


@dataclass
class SampleRecord:
    sample_id: str
    title: str
    condition: str
    gsm: str
    source_name: str
    platform_id: str
    srx: str = ""


def iter_soft_sample_blocks(soft_path: Path):
    text = soft_path.open("rb").read()
    if not text:
        return
    import gzip

    content = gzip.decompress(text).decode("utf-8", "ignore")
    for block in content.split("^SAMPLE = ")[1:]:
        lines = block.splitlines()
        if lines:
            yield lines[0].strip(), lines


def sample_first_value(lines: list[str], prefix: str) -> str:
    for line in lines:
        if line.startswith(prefix):
            return line.split(" = ", 1)[1]
    return ""


def sample_all_values(lines: list[str], prefix: str) -> list[str]:
    values: list[str] = []
    for line in lines:
        if line.startswith(prefix):
            values.append(line.split(" = ", 1)[1])
    return values


def parse_characteristics(lines: list[str]) -> dict[str, str]:
    result: dict[str, str] = {}
    for value in sample_all_values(lines, "!Sample_characteristics_ch1 = "):
        if ": " in value:
            key, raw = value.split(": ", 1)
            result[key.strip()] = raw.strip()
    return result


def download_text(url: str) -> str:
    with urlopen(url) as response:
        return response.read().decode("utf-8", "ignore")


def ensure_soft_file(path: Path, soft_url: str, log_prefix: str) -> Path:
    if path.exists() and path.stat().st_size > 0:
        log(f"Using existing SOFT file: {path}", log_prefix)
        return path
    if path.exists() and path.stat().st_size == 0:
        log(f"Existing SOFT file is empty, re-downloading: {path}", log_prefix)

    log(f"Downloading GEO SOFT file to: {path}", log_prefix)
    path.parent.mkdir(parents=True, exist_ok=True)
    with urlopen(soft_url) as response, path.open("wb") as handle:
        handle.write(response.read())
    log("SOFT download completed", log_prefix)
    return path


def fetch_runinfo(srx: str) -> dict[str, str]:
    text = download_text(EUTILS_URL.format(accession=quote_plus(srx))).strip().splitlines()
    if len(text) < 2:
        raise RuntimeError(f"No runinfo returned for {srx}")
    return next(csv.DictReader(text))


def fetch_ena(srr: str) -> dict[str, str]:
    text = download_text(ENA_URL.format(accession=quote_plus(srr))).strip().splitlines()
    if len(text) < 2:
        raise RuntimeError(f"No ENA fastq info returned for {srr}")
    return next(csv.DictReader(text, delimiter="	"))


def build_rows(samples: list[SampleRecord], selected_sample_ids: list[str], log_prefix: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    selected_samples = [sample for sample in samples if sample.sample_id in selected_sample_ids]
    selected_samples.sort(key=lambda sample: (sample.sample_id, sample.condition, sample.gsm))
    total = len(selected_samples)

    for index, sample in enumerate(selected_samples, start=1):
        base_row = {
            "sample_id": sample.sample_id,
            "title": sample.title,
            "condition": sample.condition,
            "gsm": sample.gsm,
            "source_name": sample.source_name,
            "platform_id": sample.platform_id,
            "srx": sample.srx,
        }
        if not sample.srx:
            rows.append(
                {
                    **base_row,
                    "srr": "",
                    "fastq_url": "",
                    "sample_total_bytes": "0",
                    "sample_total_gb": "0.000",
                }
            )
            continue

        log(
            f"[{index}/{total}] Resolving run metadata for {sample.sample_id} {sample.condition} ({sample.gsm}, {sample.srx})",
            log_prefix,
        )
        runinfo = fetch_runinfo(sample.srx)
        srr = runinfo["Run"]
        log(f"[{index}/{total}] Found run accession: {srr}", log_prefix)
        ena = fetch_ena(srr)
        urls = ["https://" + value for value in ena["fastq_ftp"].split(";") if value]
        sizes = [int(value) for value in ena["fastq_bytes"].split(";") if value]
        total_bytes = sum(sizes)
        log(f"[{index}/{total}] ENA returned {len(urls)} FASTQ file(s), total {total_bytes / 1e9:.3f} GB", log_prefix)
        rows.append(
            {
                **base_row,
                "srr": srr,
                "fastq_url": ";".join(urls),
                "sample_total_bytes": str(total_bytes),
                "sample_total_gb": f"{total_bytes / 1e9:.3f}",
            }
        )

    return rows


def write_tsv(
    rows: list[dict[str, str]],
    output_path: Path,
    log_prefix: str,
    fieldnames: list[str] | None = None,
) -> None:
    log(f"Writing TSV: {output_path}", log_prefix)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    effective_fieldnames = fieldnames or NORMALIZED_FIELDS
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=effective_fieldnames, delimiter="	")
        writer.writeheader()
        writer.writerows(rows)
    log("TSV writing completed", log_prefix)
