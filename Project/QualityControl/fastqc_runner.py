from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

from Log.log_util import log

from . import QualityControlConfig

LOG_PREFIX = "fastqc"
FASTQ_SUFFIX = ".fastq.gz"
PAIRED_FASTQ_RE = re.compile(r"^(?P<stem>.+)_(?P<mate>[12])\.fastq\.gz$")


def _log(message: str) -> None:
    log(message, LOG_PREFIX)


def _windows_to_wsl_path(path: Path) -> str:
    drive = path.drive.rstrip(":").lower()
    parts = [part.replace("'", "'\"'\"'") for part in path.parts[1:]]
    return f"/mnt/{drive}/{'/'.join(parts)}"


def collect_fastqc_inputs(fastq_dir: Path) -> list[Path]:
    if not fastq_dir.exists():
        raise FileNotFoundError(f"Missing FASTQ directory: {fastq_dir}")

    paired: dict[str, dict[str, Path]] = {}
    standalone: dict[str, Path] = {}

    for path in sorted(fastq_dir.glob(f"*{FASTQ_SUFFIX}")):
        match = PAIRED_FASTQ_RE.match(path.name)
        if match:
            paired.setdefault(match.group("stem"), {})[match.group("mate")] = path
        else:
            standalone[path.name[: -len(FASTQ_SUFFIX)]] = path

    selected: list[Path] = []
    for stem in sorted(paired):
        mates = paired[stem]
        if "1" in mates and "2" in mates:
            selected.extend([mates["1"], mates["2"]])
            if stem in standalone:
                _log(f"Skipping standalone FASTQ because paired files are available: {standalone[stem].name}")
            continue
        selected.extend(mates[mate] for mate in sorted(mates))

    for stem in sorted(standalone):
        if stem not in paired:
            selected.append(standalone[stem])

    if not selected:
        raise FileNotFoundError(f"No {FASTQ_SUFFIX} files found in {fastq_dir}")
    return selected


def run_fastqc(
    config: QualityControlConfig,
    executable: str = "fastqc",
    threads: int = 2,
    use_trimmed_reads: bool = False,
) -> Path:
    if threads <= 0:
        raise ValueError("threads must be positive")

    if use_trimmed_reads:
        fastq_dir = config.resolved_trimmed_fastq_dir()
        out_dir = config.resolved_fastqc_trimmed_report_out()
    else:
        fastq_dir = config.resolved_fastq_dir()
        out_dir = config.resolved_fastqc_report_out()
    inputs = collect_fastqc_inputs(fastq_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    _log(f"FASTQ directory: {fastq_dir}")
    _log(f"FastQC output directory: {out_dir}")
    _log(f"Input FASTQ files selected: {len(inputs)}")

    wsl_executable = shutil.which("wsl") or r"C:\Windows\System32\wsl.exe"
    quoted_inputs = " ".join(f"'{_windows_to_wsl_path(path)}'" for path in inputs)
    command = [
        wsl_executable,
        "bash",
        "-lc",
        f"command -v {executable} >/dev/null 2>&1 || {{ echo 'FastQC not found in WSL PATH' >&2; exit 127; }}; "
        f"{executable} --threads {threads} --outdir '{_windows_to_wsl_path(out_dir)}' {quoted_inputs}",
    ]

    _log(f"Running command: {' '.join(command)}")
    try:
        subprocess.run(command, check=True)
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            "Unable to start the WSL FastQC command. Verify that WSL is installed and reachable from Python."
        ) from exc
    _log("Done")
    return out_dir
