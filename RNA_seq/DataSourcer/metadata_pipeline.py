from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import quote_plus
from urllib.request import urlopen

LOGGER = logging.getLogger("build_metadata_table")

EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id={accession}&rettype=runinfo&retmode=text"
ENA_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport?accession={accession}&result=read_run&fields=run_accession,fastq_ftp,fastq_md5"

NORMALIZED_FIELDS = [
    "patient_id",
    "condition",
    "gsm",
    "srx",
    "srr",
    "fastq_url",
    "fastq_filename",
    "fastq_md5",
]


@dataclass
class SampleRecord:
    patient_id: str
    condition: str
    gsm: str
    srx: str = ""


def download_text(url: str) -> str:
    with urlopen(url) as response:
        return response.read().decode("utf-8", "ignore")


def ensure_soft_file(path: Path, soft_url: str) -> Path:
    if path.exists():
        LOGGER.info("Using existing SOFT file: %s", path)
        return path

    LOGGER.info("Downloading GEO SOFT file to: %s", path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with urlopen(soft_url) as response, path.open("wb") as handle:
        handle.write(response.read())
    LOGGER.info("SOFT download completed")
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


def build_rows(samples: list[SampleRecord]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    samples = sorted(samples, key=lambda s: (s.patient_id, s.condition, s.gsm))
    total = len(samples)

    for index, sample in enumerate(samples, start=1):
        base_row = {
            "patient_id": sample.patient_id,
            "condition": sample.condition,
            "gsm": sample.gsm,
            "srx": sample.srx,
        }
        if not sample.srx:
            rows.append(
                {
                    **base_row,
                    "srr": "",
                    "fastq_url": "",
                    "fastq_filename": "",
                    "fastq_md5": "",
                }
            )
            continue

        LOGGER.info(
            "[%d/%d] Resolving run metadata for %s %s (%s, %s)",
            index, total, sample.patient_id, sample.condition, sample.gsm, sample.srx,
        )
        runinfo = fetch_runinfo(sample.srx)
        srr = runinfo["Run"]
        LOGGER.info("[%d/%d] Found run accession: %s", index, total, srr)
        ena = fetch_ena(srr)
        urls = ["https://" + value for value in ena["fastq_ftp"].split(";") if value]
        filenames = [Path(value).name for value in ena["fastq_ftp"].split(";") if value]
        checksums = [value.lower() for value in ena["fastq_md5"].split(";") if value]
        if len(urls) != len(checksums):
            raise RuntimeError(
                f"ENA returned inconsistent FASTQ metadata for {srr}: "
                f"{len(urls)} URLs and {len(checksums)} MD5 checksums"
            )
        LOGGER.info("[%d/%d] ENA returned %d FASTQ file(s)", index, total, len(urls))
        rows.append(
            {
                **base_row,
                "srr": srr,
                "fastq_url": ";".join(urls),
                "fastq_filename": ";".join(filenames),
                "fastq_md5": ";".join(checksums),
            }
        )

    return rows


def write_tsv(
    rows: list[dict[str, str]],
    output_path: Path,
    fieldnames: list[str] | None = None,
) -> None:
    LOGGER.info("Writing TSV: %s", output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    effective_fieldnames = fieldnames or NORMALIZED_FIELDS
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=effective_fieldnames, delimiter="	")
        writer.writeheader()
        writer.writerows(rows)
    LOGGER.info("TSV writing completed")
