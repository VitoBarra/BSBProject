from __future__ import annotations

import csv
from concurrent.futures import FIRST_EXCEPTION, ThreadPoolExecutor, wait
from time import perf_counter
from pathlib import Path
from urllib.request import Request, urlopen

from . import DataSourceConfig
from Log import CHUNK_SIZE
from Log.log_util import WorkerProgressDisplay, format_size, format_worker_progress_line, log

LOG_PREFIX = "download_fastq"
REQUIRED_FASTQ_COLUMNS = {"sample_id", "fastq_url"}
MIN_SPEED_BPS = 256 * 1024
MIN_SPEED_GRACE_SECONDS = 30
MAX_RESTARTS_PER_FILE = 3


class SlowDownloadError(RuntimeError):
    pass


PROGRESS = WorkerProgressDisplay()


def iter_fastq_urls(meta_path: Path):
    with meta_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = set(reader.fieldnames or [])
        missing = REQUIRED_FASTQ_COLUMNS - fieldnames
        if missing:
            required = ", ".join(sorted(REQUIRED_FASTQ_COLUMNS))
            found = ", ".join(reader.fieldnames or [])
            raise ValueError(
                f"Metadata file {meta_path} does not support FASTQ download. "
                f"Required columns: {required}. Found: {found or '(none)'}"
            )
        for row in reader:
            for url in row["fastq_url"].split(";"):
                url = url.strip()
                if url:
                    yield row["sample_id"], url


def count_fastq_files(meta_path: Path) -> int:
    return sum(1 for _sample_id, _url in iter_fastq_urls(meta_path))


def _log(message: str) -> None:
    PROGRESS.log(f"[{LOG_PREFIX}] {message}")


def _download_file_once(url: str, destination: Path, index: int, total_files: int, sample_id: str, slot: int) -> None:
    if destination.exists():
        PROGRESS.finish(slot, f"[{index}/{total_files}] {sample_id} {destination.name} | already present")
        return

    _log(f"[{index}/{total_files}] Downloading {sample_id}: {destination.name}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial_destination = destination.with_name(f"{destination.name}.part")
    existing_bytes = partial_destination.stat().st_size if partial_destination.exists() else 0
    request = Request(url)
    if existing_bytes > 0:
        request.add_header("Range", f"bytes={existing_bytes}-")
        _log(
            f"[{index}/{total_files}] Found partial file, trying resume from {format_size(existing_bytes)}: {partial_destination.name}"
        )

    with urlopen(request) as response:
        status_code = getattr(response, "status", None)
        resumed = existing_bytes > 0 and status_code == 206
        if existing_bytes > 0 and not resumed:
            _log(f"[{index}/{total_files}] Server does not support resume, restarting from 0")
            existing_bytes = 0

        content_range = response.headers.get("Content-Range")
        if content_range and "/" in content_range:
            total_part = content_range.rsplit("/", 1)[-1]
            total_bytes = int(total_part) if total_part.isdigit() else None
        else:
            content_length = response.headers.get("Content-Length")
            if content_length is not None and content_length.isdigit():
                total_length = int(content_length)
                total_bytes = (existing_bytes + total_length) if resumed else total_length
            else:
                total_bytes = None

        mode = "ab" if resumed else "wb"
        downloaded = existing_bytes
        session_downloaded = 0
        started_at = perf_counter()
        low_speed_started_at: float | None = None
        with partial_destination.open(mode) as out_handle:
            while True:
                chunk = response.read(CHUNK_SIZE)
                if not chunk:
                    break
                out_handle.write(chunk)
                chunk_size = len(chunk)
                downloaded += chunk_size
                session_downloaded += chunk_size
                elapsed = max(perf_counter() - started_at, 1e-6)
                speed_bps = session_downloaded / elapsed
                PROGRESS.update(slot, format_worker_progress_line(index, total_files, sample_id, destination.name, downloaded, total_bytes, speed_bps))

                now = perf_counter()
                if speed_bps < MIN_SPEED_BPS:
                    if low_speed_started_at is None:
                        low_speed_started_at = now
                    elif now - low_speed_started_at >= MIN_SPEED_GRACE_SECONDS:
                        raise SlowDownloadError(
                            f"Download speed stayed below {format_size(MIN_SPEED_BPS)}/s for "
                            f"{MIN_SPEED_GRACE_SECONDS}s"
                        )
                else:
                    low_speed_started_at = None

    partial_destination.replace(destination)
    PROGRESS.finish(
        slot,
        f"[{index}/{total_files}] {sample_id} {destination.name} | completed {format_size(total_bytes or downloaded)}",
    )


def download_file(url: str, destination: Path, index: int, total_files: int, sample_id: str) -> None:
    with PROGRESS.slot() as slot:
        for attempt in range(1, MAX_RESTARTS_PER_FILE + 1):
            try:
                _download_file_once(url, destination, index, total_files, sample_id, slot)
                return
            except SlowDownloadError as exc:
                _log(
                    f"[{index}/{total_files}] Slow download detected for {destination.name} "
                    f"(attempt {attempt}/{MAX_RESTARTS_PER_FILE}): {exc}"
                )
                if attempt == MAX_RESTARTS_PER_FILE:
                    raise RuntimeError(
                        f"Download kept stalling for {destination.name} after {MAX_RESTARTS_PER_FILE} attempts"
                    ) from exc
                _log(f"[{index}/{total_files}] Restarting download for {destination.name}")


def download_fastq_from_tsv(config: DataSourceConfig) -> Path:
    meta_path = config.resolved_metadata_path()
    dest_dir = config.resolved_fastq_dest()
    if not meta_path.exists():
        raise FileNotFoundError(f"Missing metadata file: {meta_path}")

    jobs = []
    for index, (sample_id, url) in enumerate(iter_fastq_urls(meta_path), start=1):
        filename = url.rsplit("/", 1)[-1]
        jobs.append((index, sample_id, url, dest_dir / filename))

    total_files = len(jobs)
    log(f"Using metadata TSV: {meta_path}", LOG_PREFIX)
    log(f"Destination directory: {dest_dir}", LOG_PREFIX)
    log(f"Total FASTQ files listed: {total_files}", LOG_PREFIX)
    log(f"Parallel workers: {config.download_workers}", LOG_PREFIX)

    if total_files == 0:
        log("Done", LOG_PREFIX)
        return dest_dir

    max_workers = min(config.download_workers, total_files)
    PROGRESS.initialize(max_workers)
    try:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [
                executor.submit(download_file, url, destination, index, total_files, sample_id)
                for index, sample_id, url, destination in jobs
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
    finally:
        PROGRESS.clear()

    log("Done", LOG_PREFIX)
    return dest_dir
