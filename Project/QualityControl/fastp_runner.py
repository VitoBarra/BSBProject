from __future__ import annotations

import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from Log.log_util import log

from . import QualityControlConfig
from .fastqc_runner import _windows_to_wsl_path

LOG_PREFIX = "fastp"
FASTQ_SUFFIX = ".fastq.gz"
PAIRED_FASTQ_RE = re.compile(r"^(?P<stem>.+)_(?P<mate>[12])\.fastq\.gz$")


@dataclass(slots=True, frozen=True)
class FastpJob:
    name: str
    input_fastq_1: Path
    output_fastq_1: Path
    html_report: Path
    json_report: Path
    input_fastq_2: Path | None = None
    output_fastq_2: Path | None = None


def _log(message: str) -> None:
    log(message, LOG_PREFIX)


def _shell_quote(value: str) -> str:
    return "'" + value.replace("'", "'\"'\"'") + "'"


def _job_name(path: Path) -> str:
    if path.name.endswith(FASTQ_SUFFIX):
        return path.name[: -len(FASTQ_SUFFIX)]
    return path.stem


def collect_fastp_inputs(fastq_dir: Path) -> list[Path]:
    if not fastq_dir.exists():
        raise FileNotFoundError(f"Missing FASTQ directory: {fastq_dir}")

    fastq_files = sorted(fastq_dir.glob(f"*{FASTQ_SUFFIX}"))
    if not fastq_files:
        raise FileNotFoundError(f"No {FASTQ_SUFFIX} files found in {fastq_dir}")
    return fastq_files

def _build_fastp_jobs(fastq_dir: Path, trimmed_dir: Path, report_dir: Path) -> list[FastpJob]:
    input_fastq_files = collect_fastp_inputs(fastq_dir)

    paired: dict[str, dict[str, Path]] = {}
    standalone: dict[str, Path] = {}
    for path in input_fastq_files:
        match = PAIRED_FASTQ_RE.match(path.name)
        if match:
            paired.setdefault(match.group("stem"), {})[match.group("mate")] = path
        else:
            standalone[path.name[: -len(FASTQ_SUFFIX)]] = path

    jobs: list[FastpJob] = []
    for stem in sorted(paired):
        mates = paired[stem]
        if "1" in mates and "2" in mates:
            jobs.append(
                FastpJob(
                    name=stem,
                    input_fastq_1=mates["1"],
                    output_fastq_1=trimmed_dir / mates["1"].name,
                    input_fastq_2=mates["2"],
                    output_fastq_2=trimmed_dir / mates["2"].name,
                    html_report=report_dir / f"{stem}.html",
                    json_report=report_dir / f"{stem}.json",
                )
            )
            if stem in standalone:
                _log(f"Skipping standalone FASTQ because paired files are available: {standalone[stem].name}")
            continue

        for mate in sorted(mates):
            path = mates[mate]
            name = _job_name(path)
            _log(f"Only one mate found for {path.name}; trimming it as single-end input")
            jobs.append(
                FastpJob(
                    name=name,
                    input_fastq_1=path,
                    output_fastq_1=trimmed_dir / path.name,
                    html_report=report_dir / f"{name}.html",
                    json_report=report_dir / f"{name}.json",
                )
            )

    for stem in sorted(standalone):
        if stem in paired:
            continue
        path = standalone[stem]
        jobs.append(
            FastpJob(
                name=stem,
                input_fastq_1=path,
                output_fastq_1=trimmed_dir / path.name,
                html_report=report_dir / f"{stem}.html",
                json_report=report_dir / f"{stem}.json",
            )
        )

    return jobs


def _build_command(job: FastpJob, executable: str, threads: int) -> str:
    command_parts = [
        executable,
        "--thread",
        str(threads),
        "--html",
        _shell_quote(_windows_to_wsl_path(job.html_report)),
        "--json",
        _shell_quote(_windows_to_wsl_path(job.json_report)),
    ]

    command_parts.extend([
        "--in1",
        _shell_quote(_windows_to_wsl_path(job.input_fastq_1)),
        "--out1",
        _shell_quote(_windows_to_wsl_path(job.output_fastq_1)),
    ])

    if job.input_fastq_2 is not None and job.output_fastq_2 is not None:
        command_parts.extend([
            "--in2",
            _shell_quote(_windows_to_wsl_path(job.input_fastq_2)),
            "--out2",
            _shell_quote(_windows_to_wsl_path(job.output_fastq_2)),
        ])

    return " ".join(command_parts)


def run_fastp(
    config: QualityControlConfig,
    executable: str = "fastp",
    threads: int = 2,
) -> Path:
    if threads <= 0:
        raise ValueError("threads must be positive")

    fastq_dir = config.resolved_fastq_dir()
    trimmed_dir = config.resolved_trimmed_fastq_dir()
    report_dir = config.resolved_fastp_report_out()
    if not fastq_dir.exists():
        raise FileNotFoundError(f"Missing FASTQ directory: {fastq_dir}")

    jobs = _build_fastp_jobs(fastq_dir, trimmed_dir, report_dir)
    trimmed_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)

    _log(f"FASTQ input directory: {fastq_dir}")
    _log(f"Trimmed FASTQ output directory: {trimmed_dir}")
    _log(f"fastp report directory: {report_dir}")
    _log(f"Input datasets selected: {len(jobs)}")

    wsl_executable = shutil.which("wsl") or r"C:\Windows\System32\wsl.exe"
    for index, job in enumerate(jobs, start=1):
        command = [
            wsl_executable,
            "bash",
            "-lc",
            f"command -v {executable} >/dev/null 2>&1 || {{ echo 'fastp not found in WSL PATH' >&2; exit 127; }}; "
            f"{_build_command(job, executable=executable, threads=threads)}",
        ]
        _log(f"[{index}/{len(jobs)}] Trimming {job.name}")
        _log(f"Running command: {' '.join(command)}")
        try:
            subprocess.run(command, check=True)
        except FileNotFoundError as exc:
            raise FileNotFoundError(
                "Unable to start the WSL fastp command. Verify that WSL is installed and reachable from Python."
            ) from exc

    _log("Done")
    return trimmed_dir
