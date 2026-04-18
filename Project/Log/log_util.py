from __future__ import annotations

import os
import sys
import threading
from contextlib import contextmanager
from queue import SimpleQueue
from time import perf_counter

from Log import BAR_WIDTH


def log(message: str, prefix: str | None = None) -> None:
    if prefix:
        print(f"[{prefix}] {message}", flush=True)
    else:
        print(message, flush=True)


def format_size(num_bytes: int) -> str:
    units = ["B", "KB", "MB", "GB", "TB"]
    size = float(num_bytes)
    for unit in units:
        if size < 1024 or unit == units[-1]:
            return f"{size:.1f}{unit}"
        size /= 1024
    return f"{num_bytes}B"


def format_speed(bytes_per_second: float | None) -> str:
    if bytes_per_second is None or bytes_per_second <= 0:
        return "--B/s"
    return f"{format_size(int(bytes_per_second))}/s"


def format_eta(seconds: float | None) -> str:
    if seconds is None or seconds < 0:
        return "N/A"
    total_seconds = int(round(seconds))
    hours, remainder = divmod(total_seconds, 3600)
    minutes, secs = divmod(remainder, 60)
    if hours > 0:
        return f"{hours:02d}:{minutes:02d}:{secs:02d}"
    return f"{minutes:02d}:{secs:02d}"


def render_progress(downloaded: int, total: int | None, speed_bps: float | None = None) -> None:
    speed_text = format_speed(speed_bps)
    if total and total > 0:
        ratio = min(downloaded / total, 1.0)
        filled = int(BAR_WIDTH * ratio)
        bar = "#" * filled + "-" * (BAR_WIDTH - filled)
        percent = ratio * 100
        remaining_bytes = max(total - downloaded, 0)
        eta_seconds = (remaining_bytes / speed_bps) if speed_bps and speed_bps > 0 else None
        eta_text = format_eta(eta_seconds)
        message = (
            f"\rProgress: [{bar}] {percent:6.2f}% | "
            f"Downloaded: {format_size(downloaded)} | "
            f"Total: {format_size(total)} | "
            f"Speed: {speed_text} | "
            f"ETA: {eta_text}"
        )
    else:
        message = (
            f"\rDownloaded: {format_size(downloaded)} | "
            f"Speed: {speed_text} | "
            f"ETA: N/A"
        )

    try:
        print(message, end="", flush=True)
    except OSError:
        return


def format_worker_progress_line(
    index: int,
    total_files: int,
    sample_id: str,
    filename: str,
    downloaded: int,
    total_bytes: int | None,
    speed_bps: float | None,
    bar_width: int = 20,
) -> str:
    prefix = f"[{index}/{total_files}] {sample_id} {filename}"
    speed_text = format_speed(speed_bps)
    if total_bytes and total_bytes > 0:
        ratio = min(downloaded / total_bytes, 1.0)
        filled = int(bar_width * ratio)
        bar = "#" * filled + "-" * (bar_width - filled)
        percent = ratio * 100
        remaining_bytes = max(total_bytes - downloaded, 0)
        eta_seconds = (remaining_bytes / speed_bps) if speed_bps and speed_bps > 0 else None
        return (
            f"{prefix} [{bar}] {percent:6.2f}% | "
            f"{format_size(downloaded)}/{format_size(total_bytes)} | "
            f"{speed_text} | ETA {format_eta(eta_seconds)}"
        )
    return f"{prefix} | {format_size(downloaded)} | {speed_text} | ETA N/A"


class WorkerProgressDisplay:
    FALLBACK_UPDATE_INTERVAL_SECONDS = 1.0

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._lines: list[str] = []
        self._active = False
        self._slot_pool: SimpleQueue[int] = SimpleQueue()
        self._enabled = self._detect_console_support()
        self._last_fallback_lines: dict[int, str] = {}
        self._last_fallback_times: dict[int, float] = {}

    def _detect_console_support(self) -> bool:
        if os.name != "nt":
            return sys.stdout.isatty()
        try:
            import ctypes

            kernel32 = ctypes.windll.kernel32
            handle = kernel32.GetStdHandle(-11)
            if handle in (0, -1):
                return False
            mode = ctypes.c_uint()
            if not kernel32.GetConsoleMode(handle, ctypes.byref(mode)):
                return False
            kernel32.SetConsoleMode(handle, mode.value | 0x0004)
            return True
        except Exception:
            return False

    def initialize(self, workers: int) -> None:
        with self._lock:
            self._lines = [""] * workers
            self._active = workers > 0 and self._enabled
            self._last_fallback_lines.clear()
            self._last_fallback_times.clear()
            while not self._slot_pool.empty():
                self._slot_pool.get()
            for slot in range(workers):
                self._slot_pool.put(slot)
            if self._active:
                sys.stdout.write("\n" * workers)
                sys.stdout.flush()

    def _move_to_block_start(self) -> None:
        if self._active and self._lines:
            sys.stdout.write(f"\x1b[{len(self._lines)}F")

    def _redraw(self) -> None:
        if not self._active:
            return
        self._move_to_block_start()
        for line in self._lines:
            sys.stdout.write("\x1b[2K")
            sys.stdout.write((line or "") + "\n")
        sys.stdout.flush()

    def _fallback_print(self, slot: int, line: str, force: bool = False) -> None:
        now = perf_counter()
        last_line = self._last_fallback_lines.get(slot)
        last_time = self._last_fallback_times.get(slot, 0.0)
        if force or line != last_line or now - last_time >= self.FALLBACK_UPDATE_INTERVAL_SECONDS:
            print(line, flush=True)
            self._last_fallback_lines[slot] = line
            self._last_fallback_times[slot] = now

    @contextmanager
    def slot(self):
        slot = self._slot_pool.get()
        try:
            yield slot
        finally:
            with self._lock:
                if self._active and 0 <= slot < len(self._lines):
                    self._lines[slot] = ""
                    self._redraw()
                else:
                    self._last_fallback_lines.pop(slot, None)
                    self._last_fallback_times.pop(slot, None)
            self._slot_pool.put(slot)

    def update(self, slot: int, line: str) -> None:
        with self._lock:
            if self._active:
                self._lines[slot] = line
                self._redraw()
            else:
                self._fallback_print(slot, line)

    def finish(self, slot: int, line: str) -> None:
        with self._lock:
            if self._active:
                self._lines[slot] = line
                self._redraw()
            else:
                self._fallback_print(slot, line, force=True)

    def log(self, message: str) -> None:
        with self._lock:
            if self._active and self._lines:
                self._move_to_block_start()
                sys.stdout.write("\x1b[J")
                print(message, flush=True)
                sys.stdout.write("\n" * len(self._lines))
                self._redraw()
            else:
                print(message, flush=True)

    def clear(self) -> None:
        with self._lock:
            if self._active and self._lines:
                self._move_to_block_start()
                for _ in self._lines:
                    sys.stdout.write("\x1b[2K\n")
                self._move_to_block_start()
                sys.stdout.flush()
            self._active = False
            self._lines = []
            self._last_fallback_lines.clear()
            self._last_fallback_times.clear()
            while not self._slot_pool.empty():
                self._slot_pool.get()
