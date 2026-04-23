from __future__ import annotations

from pathlib import Path
from time import perf_counter
from urllib.request import Request, urlopen

from Log import CHUNK_SIZE
from Log.log_util import format_size, log, render_progress

from . import DataSourceConfig

LOG_PREFIX = "download_reference"


def _log(message: str) -> None:
    log(message, LOG_PREFIX)


def _download(url: str, destination: Path) -> None:
    if destination.exists():
        _log(f"Reference transcriptome already present: {destination}")
        return

    destination.parent.mkdir(parents=True, exist_ok=True)
    partial_destination = destination.with_name(f"{destination.name}.part")
    existing_bytes = partial_destination.stat().st_size if partial_destination.exists() else 0

    request = Request(url)
    if existing_bytes > 0:
        request.add_header("Range", f"bytes={existing_bytes}-")
        _log(f"Found partial transcriptome file, trying resume from {format_size(existing_bytes)}")

    with urlopen(request) as response:
        status_code = getattr(response, "status", None)
        resumed = existing_bytes > 0 and status_code == 206
        if existing_bytes > 0 and not resumed:
            _log("Server does not support resume, restarting from 0")
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

        with partial_destination.open(mode) as out_handle:
            while True:
                chunk = response.read(CHUNK_SIZE)
                if not chunk:
                    break
                out_handle.write(chunk)
                downloaded += len(chunk)
                session_downloaded += len(chunk)
                elapsed = max(perf_counter() - started_at, 1e-6)
                render_progress(downloaded, total_bytes, session_downloaded / elapsed)

    print()
    partial_destination.replace(destination)
    _log(f"Saved reference transcriptome to {destination}")


def download_reference_transcriptome(config: DataSourceConfig) -> Path:
    url = config.resolved_transcriptome_url()
    destination = config.resolved_transcriptome_fasta_path()

    _log(f"Dataset: {config.profile.accession}")
    _log(f"Transcriptome URL: {url}")
    _log(f"Destination file: {destination}")
    _download(url, destination)
    _log("Done")
    return destination
