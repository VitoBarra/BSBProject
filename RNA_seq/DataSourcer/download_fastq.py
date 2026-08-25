from __future__ import annotations

import csv
import hashlib
import logging
from collections.abc import Iterator
from concurrent.futures import FIRST_EXCEPTION, ThreadPoolExecutor, wait
from pathlib import Path

from . import DataSourceConfig, ensure_metadata_table
from .http_download import download_resumable, format_size

LOGGER = logging.getLogger("download_fastq")
MIN_SPEED_BPS = 256 * 1024
MIN_SPEED_GRACE_SECONDS = 120
USER_AGENT = "BSBProject-GSE103001/1.0"


def _calculate_md5(file_path: Path) -> str:
    """Return the MD5 checksum of a file."""
    digest = hashlib.md5(usedforsecurity=False)
    with file_path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def iter_fastq_files(metadata_path: Path) -> Iterator[tuple[str, str, str, str]]:
    """Yield each FASTQ file's patient ID, URL, filename, and MD5 checksum."""
    with metadata_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            urls = [value.strip() for value in row["fastq_url"].split(";") if value.strip()]
            filenames = [value.strip() for value in row["fastq_filename"].split(";") if value.strip()]
            checksums = [value.strip().lower() for value in row["fastq_md5"].split(";") if value.strip()]
            for url, filename, checksum in zip(urls, filenames, checksums):
                yield row["patient_id"], url, filename, checksum


def _download_fastq_file(
    url: str,
    destination: Path,
    expected_md5: str,
    index: int,
    total_files: int,
    patient_id: str,
) -> None:
    LOGGER.info("[%d/%d] Downloading %s: %s", index, total_files, patient_id, destination.name)
    already_present, downloaded = download_resumable(
        url,
        destination,
        user_agent=USER_AGENT,
        timeout=120,
        minimum_speed_bps=MIN_SPEED_BPS,
        minimum_speed_grace_seconds=MIN_SPEED_GRACE_SECONDS,
    )
    actual_md5 = _calculate_md5(destination)
    if actual_md5 != expected_md5:
        raise RuntimeError(
            f"MD5 mismatch for {destination.name}: expected {expected_md5}, got {actual_md5}"
        )

    if already_present:
        LOGGER.info("[%d/%d] %s %s already present", index, total_files, patient_id, destination.name)
        return
    LOGGER.info(
        "[%d/%d] %s %s completed (%s)",
        index, total_files, patient_id, destination.name, format_size(downloaded),
    )


def download_fastq(config: DataSourceConfig) -> None:
    """Download FASTQ files listed in the configured metadata TSV.

    Downloads files concurrently, resumes incomplete downloads when possible,
    and validates each file against its metadata MD5 checksum.

    Args:
        config: Data-source settings that provide the metadata path, output
            directory, and number of download workers.

    Raises:
        FileNotFoundError: If the configured metadata TSV does not exist.
    """
    ensure_metadata_table(config)
    metadata_path = config.resolved_metadata_path()
    destination_dir = config.resolved_fastq_dest()

    if not metadata_path.exists():
        raise FileNotFoundError(f"Missing metadata file: {metadata_path}")

    download_jobs = [
        (index, patient_id, url, destination_dir / filename, checksum)
        for index, (patient_id, url, filename, checksum) in enumerate(
            iter_fastq_files(metadata_path), start=1
        )
    ]

    total_files = len(download_jobs)
    LOGGER.info("Using metadata TSV: %s", metadata_path)
    LOGGER.info("Destination directory: %s", destination_dir)
    LOGGER.info("Total FASTQ files listed: %d", total_files)
    LOGGER.info("Parallel workers: %d", config.download_workers)

    if total_files == 0:
        LOGGER.info("Done")
        return


    max_workers = min(config.download_workers, total_files)
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(
                _download_fastq_file,
                url,
                destination,
                checksum,
                index,
                total_files,
                patient_id,
            )
            for index, patient_id, url, destination, checksum in download_jobs
        ]
        done, not_done = wait(futures, return_when=FIRST_EXCEPTION)
        for future in done:
            exc = future.exception()
            if exc is not None:
                for pending in not_done:
                    pending.cancel()
                raise exc
        for future in not_done:
            future.result()

    LOGGER.info("Done")
