from __future__ import annotations

import logging
import shutil
from math import ceil
from pathlib import Path
from time import monotonic, sleep
from urllib.error import HTTPError

import pycurl

LOGGER = logging.getLogger("http_download")
DEFAULT_MAX_ATTEMPTS = 20
DEFAULT_RETRY_DELAY_SECONDS = 30
DEFAULT_MINIMUM_SPEED_BPS = 1024
DEFAULT_MINIMUM_SPEED_GRACE_SECONDS = 120
PROGRESS_LOG_INTERVAL_SECONDS = 30
RETRYABLE_HTTP_STATUS_CODES = frozenset({403, 408, 429, 500, 502, 503, 504})


class DownloadError(RuntimeError):
    """A transport error reported by libcurl."""


def format_size(num_bytes: int) -> str:
    units = ("B", "KB", "MB", "GB", "TB")
    size = float(num_bytes)
    for unit in units:
        if size < 1024 or unit == units[-1]:
            return f"{size:.1f}{unit}"
        size /= 1024
    return f"{num_bytes}B"


def _configure(
    curl: pycurl.Curl,
    url: str,
    user_agent: str,
    timeout: float | None,
    minimum_speed_bps: int | None,
    minimum_speed_grace_seconds: float,
) -> None:
    curl.setopt(curl.URL, url)
    curl.setopt(curl.USERAGENT, user_agent)
    curl.setopt(curl.FOLLOWLOCATION, True)
    curl.setopt(curl.FAILONERROR, True)
    curl.setopt(curl.NOSIGNAL, True)
    if timeout is not None:
        curl.setopt(curl.CONNECTTIMEOUT, ceil(timeout))
    if minimum_speed_bps is not None:
        curl.setopt(curl.LOW_SPEED_LIMIT, minimum_speed_bps)
        curl.setopt(curl.LOW_SPEED_TIME, ceil(minimum_speed_grace_seconds))


def _perform(curl: pycurl.Curl, url: str) -> None:
    try:
        curl.perform()
    except pycurl.error as exc:
        status = curl.getinfo(curl.RESPONSE_CODE)
        if status >= 400:
            raise HTTPError(url, status, str(exc), hdrs=None, fp=None) from exc
        raise DownloadError(str(exc)) from exc


def remote_size(url: str, user_agent: str, timeout: float = 60) -> int | None:
    curl = pycurl.Curl()
    try:
        _configure(curl, url, user_agent, timeout, None, 0)
        curl.setopt(curl.NOBODY, True)
        _perform(curl, url)
        size = curl.getinfo(curl.CONTENT_LENGTH_DOWNLOAD_T)
        return size if size >= 0 else None
    finally:
        curl.close()


def _partial_path(destination: Path) -> Path:
    return destination.with_name(f"{destination.name}.part")


def _prepare_destination(
    url: str,
    destination: Path,
    user_agent: str,
) -> tuple[bool, int, int | None]:
    """Inspect final and partial files before starting a transfer."""
    expected_bytes = remote_size(url, user_agent)

    # both complete file and partial file are missing, nothing to do
    partial_destination = _partial_path(destination)
    if not destination.exists() and not partial_destination.exists():
        return False, 0, expected_bytes

    # if the final file exists, check its size against the expected size
    if destination.exists():
        actual_bytes = destination.stat().st_size
        if expected_bytes is None or actual_bytes == expected_bytes:
            return True, actual_bytes, expected_bytes

        # if the final file is larger or smaller than expected, discard it and log a warning
        destination.unlink()
        LOGGER.warning(
            "Discarded final file with unexpected size (%s/%s): %s",
            format_size(actual_bytes),
            format_size(expected_bytes),
            destination.name,
        )

    # destination file is now missing after the above check, so check the partial file

    # if the partial file is missing, nothing to do
    if not partial_destination.exists():
        return False, 0, expected_bytes

    # if the partial file is the same size as the expected size, move it to the final destination
    partial_bytes = partial_destination.stat().st_size
    if expected_bytes is not None and partial_bytes == expected_bytes:
        partial_destination.replace(destination)
        return True, partial_bytes, expected_bytes

    # if the partial file is larger than the expected size (corrupted or wrong), discard it and log a warning
    if expected_bytes is not None and partial_bytes > expected_bytes:
        partial_destination.unlink()
        LOGGER.warning(
            "Discarded oversized partial file (%s/%s): %s",
            format_size(partial_bytes),
            format_size(expected_bytes),
            partial_destination.name,
        )
        return False, 0, expected_bytes

    # if the partial file is smaller than the expected size, we can resume from it
    return False, partial_bytes, expected_bytes


def _download_once(
    url: str,
    partial_destination: Path,
    *,
    user_agent: str,
    timeout: float | None,
    minimum_speed_bps: int | None,
    minimum_speed_grace_seconds: float,
    existing_bytes: int,
) -> tuple[int, int | None, int]:
    """Download one request, without corrupting an existing partial file.

    A failed server/proxy can ignore a Range request and start sending the
    resource from byte zero.  Download resumed responses to a separate file
    first, then append only after libcurl confirms the expected 206 response
    and Content-Range.  This preserves the valid prefix for the next retry.
    """
    curl = pycurl.Curl()
    resume_destination = partial_destination.with_name(
        f"{partial_destination.name}.resume"
    )
    response_content_range: str | None = None

    def capture_headers(header: bytes) -> int:
        nonlocal response_content_range
        line = header.decode("iso-8859-1").strip()
        if line.lower().startswith("content-range:"):
            response_content_range = line.split(":", 1)[1].strip()
        return len(header)

    try:
        _configure(curl, url, user_agent, timeout, minimum_speed_bps, minimum_speed_grace_seconds)
        download_destination = resume_destination if existing_bytes else partial_destination
        with download_destination.open("wb") as handle:
            curl.setopt(curl.WRITEDATA, handle)
            curl.setopt(curl.HEADERFUNCTION, capture_headers)
            _configure_progress_logging(curl, partial_destination, existing_bytes)
            if existing_bytes:
                curl.setopt(curl.RESUME_FROM_LARGE, existing_bytes)
            _perform(curl, url)
        response_bytes = curl.getinfo(curl.CONTENT_LENGTH_DOWNLOAD_T)
        status = curl.getinfo(curl.RESPONSE_CODE)
        if existing_bytes:
            expected_range_prefix = f"bytes {existing_bytes}-"
            received_bytes = resume_destination.stat().st_size
            range_length_matches = False
            if response_content_range and response_content_range.lower().startswith(expected_range_prefix):
                try:
                    byte_range = response_content_range.split(" ", 1)[1].split("/", 1)[0]
                    range_start, range_end = (int(value) for value in byte_range.split("-", 1))
                    range_length_matches = range_end - range_start + 1 == received_bytes
                except (IndexError, ValueError):
                    pass
            if status not in (200, 206) or (
                status == 206
                and (
                    not response_content_range
                    or not range_length_matches
                )
            ):
                raise DownloadError(
                    "Server returned an invalid range response while resuming "
                    f"at {format_size(existing_bytes)} (HTTP {status}, "
                    f"Content-Range: {response_content_range or 'missing'})"
                )
            with resume_destination.open("rb") as source, partial_destination.open("ab") as destination_handle:
                shutil.copyfileobj(source, destination_handle)
        total_bytes = existing_bytes + response_bytes if status == 206 else response_bytes
        return partial_destination.stat().st_size, total_bytes if total_bytes >= 0 else None, status
    finally:
        curl.close()
        resume_destination.unlink(missing_ok=True)


def _configure_progress_logging(curl: pycurl.Curl, destination: Path, existing_bytes: int) -> None:
    """Log an active-download heartbeat without flooding the console."""
    last_log_time = monotonic()

    def report_progress(download_total: float, downloaded_now: float, _upload_total: float, _uploaded_now: float) -> int:
        nonlocal last_log_time
        now = monotonic()
        if now - last_log_time < PROGRESS_LOG_INTERVAL_SECONDS:
            return 0

        downloaded = existing_bytes + int(downloaded_now)
        total = existing_bytes + int(download_total) if download_total > 0 else None
        progress = f" of {format_size(total)} ({downloaded / total:.1%})" if total else ""
        LOGGER.info("Download still in progress: %s downloaded%s: %s", format_size(downloaded), progress, destination.name)
        last_log_time = now
        return 0

    curl.setopt(curl.NOPROGRESS, False)
    curl.setopt(curl.XFERINFOFUNCTION, report_progress)



def _download_resumable_once(
    url: str,
    destination: Path,
    *,
    user_agent: str,
    timeout: float | None,
    minimum_speed_bps: int | None,
    minimum_speed_grace_seconds: float,
) -> tuple[bool, int]:
    download_completed, downloaded_bytes, expected_bytes = _prepare_destination(
        url, destination, user_agent
    )
    if download_completed:
        return True, downloaded_bytes

    destination.parent.mkdir(parents=True, exist_ok=True)

    partial_destination = _partial_path(destination)
    existing_bytes = downloaded_bytes
    if existing_bytes:
        LOGGER.info(
            "Found partial file, trying resume from %s: %s",
            format_size(existing_bytes),
            partial_destination.name,
        )
    expected_size = f" of {format_size(expected_bytes)}" if expected_bytes is not None else ""
    LOGGER.info(
        "Starting download%s; progress will be reported every %ds: %s",
        expected_size,
        PROGRESS_LOG_INTERVAL_SECONDS,
        destination.name,
    )
    try:
        downloaded, response_total_bytes, status = _download_once(
            url, partial_destination, user_agent=user_agent, timeout=timeout,
            minimum_speed_bps=minimum_speed_bps, minimum_speed_grace_seconds=minimum_speed_grace_seconds,
            existing_bytes=existing_bytes,
        )
    except HTTPError as exc:
        if not existing_bytes or exc.code != 416:
            raise

        LOGGER.info("Server rejected the partial range, restarting from 0: %s", partial_destination.name)
        existing_bytes = 0
        downloaded, response_total_bytes, status = _download_once(
            url, partial_destination, user_agent=user_agent, timeout=timeout,
            minimum_speed_bps=minimum_speed_bps, minimum_speed_grace_seconds=minimum_speed_grace_seconds,
            existing_bytes=0,
        )

    if existing_bytes and status == 200:
        LOGGER.info("Server ignored the partial range, restarting from 0: %s", partial_destination.name)
        downloaded, response_total_bytes, status = _download_once(
            url, partial_destination, user_agent=user_agent, timeout=timeout,
            minimum_speed_bps=minimum_speed_bps, minimum_speed_grace_seconds=minimum_speed_grace_seconds,
            existing_bytes=0,
        )
    elif existing_bytes and status != 206:
        raise DownloadError(f"Expected 206 for resumed download, got HTTP {status}")

    total_bytes = expected_bytes if expected_bytes is not None else response_total_bytes
    if total_bytes is not None and downloaded != total_bytes:
        raise DownloadError(f"Connection ended after {format_size(downloaded)}; expected {format_size(total_bytes)}")

    partial_destination.replace(destination)
    return False, downloaded


def download_resumable(
    url: str,
    destination: Path,
    *,
    user_agent: str,
    timeout: float | None = None,
    minimum_speed_bps: int | None = DEFAULT_MINIMUM_SPEED_BPS,
    minimum_speed_grace_seconds: float = DEFAULT_MINIMUM_SPEED_GRACE_SECONDS,
    max_attempts: int = DEFAULT_MAX_ATTEMPTS,
    retry_delay_seconds: float = DEFAULT_RETRY_DELAY_SECONDS,
) -> tuple[bool, int]:
    """Download atomically via libcurl, retrying transient failures.

    Incomplete transfers are stored in a ``.part`` file and resumed on a
    subsequent attempt.

    Returns:
        A pair containing whether the completed file was already present and
        the file size in bytes.
    """

    if max_attempts < 1:
        raise ValueError("max_attempts must be at least 1")

    if retry_delay_seconds < 0:
        raise ValueError("retry_delay_seconds cannot be negative")

    for attempt in range(1, max_attempts + 1):
        try:
            return _download_resumable_once(
                url, destination, user_agent=user_agent, timeout=timeout,
                minimum_speed_bps=minimum_speed_bps,
                minimum_speed_grace_seconds=minimum_speed_grace_seconds,
            )

        except (DownloadError, HTTPError) as exc:
            retryable = not isinstance(exc, HTTPError) or exc.code in RETRYABLE_HTTP_STATUS_CODES
            if not retryable or attempt == max_attempts:
                raise
            LOGGER.warning(
                "Transient download failure for %s (attempt %d/%d): %s",
                destination.name, attempt, max_attempts, exc,
            )
            LOGGER.info("Retrying in %gs: %s", retry_delay_seconds, destination.name)
            sleep(retry_delay_seconds)

    raise AssertionError("unreachable")
