from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from ExternalTools import ExternalToolRunner
from Log.log_util import log

from . import QuantificationConfig

LOG_PREFIX = "salmon"
FASTQ_SUFFIX = ".fastq.gz"
PAIRED_FASTQ_RE = re.compile(r"^(?P<stem>.+)_(?P<mate>[12])\.fastq\.gz$")


@dataclass(slots=True, frozen=True)
class SalmonJob:
    name: str
    read_1: Path
    read_2: Path
    output_dir: Path


def _log(message: str) -> None:
    log(message, LOG_PREFIX)


def collect_salmon_inputs(trimmed_fastq_dir: Path, quant_dir: Path) -> list[SalmonJob]:
    if not trimmed_fastq_dir.exists():
        raise FileNotFoundError(f"Missing trimmed FASTQ directory: {trimmed_fastq_dir}")

    paired: dict[str, dict[str, Path]] = {}
    standalone: dict[str, Path] = {}
    for path in sorted(trimmed_fastq_dir.glob(f"*{FASTQ_SUFFIX}")):
        match = PAIRED_FASTQ_RE.match(path.name)
        if match:
            paired.setdefault(match.group("stem"), {})[match.group("mate")] = path
        else:
            standalone[path.name[: -len(FASTQ_SUFFIX)]] = path

    jobs: list[SalmonJob] = []
    for stem in sorted(paired):
        mates = paired[stem]
        if "1" not in mates or "2" not in mates:
            missing = "1" if "1" not in mates else "2"
            raise FileNotFoundError(f"Missing mate _{missing} for paired-end sample '{stem}' in {trimmed_fastq_dir}")
        jobs.append(
            SalmonJob(
                name=stem,
                read_1=mates["1"],
                read_2=mates["2"],
                output_dir=quant_dir / stem,
            )
        )
        if stem in standalone:
            _log(f"Skipping standalone FASTQ because paired files are available: {standalone[stem].name}")

    unexpected_standalone = sorted(stem for stem in standalone if stem not in paired)
    if unexpected_standalone:
        names = ", ".join(f"{stem}{FASTQ_SUFFIX}" for stem in unexpected_standalone)
        raise FileNotFoundError(
            f"Found standalone FASTQ files without paired mates in {trimmed_fastq_dir}: {names}"
        )

    if not jobs:
        raise FileNotFoundError(f"No paired trimmed FASTQ files found in {trimmed_fastq_dir}")

    return jobs


def _build_args(
    job: SalmonJob,
    runner: ExternalToolRunner,
    index_dir: Path,
    libtype: str,
    threads: int,
) -> list[str]:
    return [
        "quant",
        "--index",
        runner.path_arg(index_dir),
        "--libType",
        libtype,
        "--threads",
        str(threads),
        "--validateMappings",
        "--mates1",
        runner.path_arg(job.read_1),
        "--mates2",
        runner.path_arg(job.read_2),
        "--output",
        runner.path_arg(job.output_dir),
    ]


def _build_index_args(
    runner: ExternalToolRunner,
    transcriptome_fasta: Path,
    index_dir: Path,
    threads: int,
) -> list[str]:
    return [
        "index",
        "--transcripts",
        runner.path_arg(transcriptome_fasta),
        "--index",
        runner.path_arg(index_dir),
        "--threads",
        str(threads),
    ]


def _ensure_salmon_index(
    config: QuantificationConfig,
    executable: str,
    threads: int,
) -> Path:
    index_dir = config.resolved_salmon_index_dir()
    if index_dir.exists():
        return index_dir

    transcriptome_fasta = config.resolved_salmon_transcriptome_fasta()
    if transcriptome_fasta is None:
        default_reference_dir = config.reference_dir()
        raise FileNotFoundError(
            "Missing Salmon index directory and no transcriptome FASTA was found to build it. "
            f"Expected an index at {index_dir} or a transcriptome FASTA under {default_reference_dir} "
            "named one of: transcriptome.fa(.gz), transcriptome.fasta(.gz), transcripts.fa(.gz), transcripts.fasta(.gz). "
            "Alternatively, pass --salmon-index-dir to a prebuilt index or --salmon-transcriptome-fasta plus --build-salmon-index."
        )

    if not transcriptome_fasta.exists():
        raise FileNotFoundError(f"Missing Salmon transcriptome FASTA: {transcriptome_fasta}")

    index_dir.parent.mkdir(parents=True, exist_ok=True)
    _log(f"Building Salmon index from transcriptome FASTA: {transcriptome_fasta}")
    _log(f"Salmon index output directory: {index_dir}")

    runner = ExternalToolRunner(executable=executable, display_name="Salmon", log=_log)
    runner.run(
        _build_index_args(
            runner=runner,
            transcriptome_fasta=transcriptome_fasta,
            index_dir=index_dir,
            threads=threads,
        ),
        missing_message="salmon not found in PATH",
    )

    return index_dir


def run_salmon(
    config: QuantificationConfig,
    executable: str = "salmon",
    threads: int = 6,
    libtype: str = "A",
    build_index_if_missing: bool = False,
) -> Path:
    if threads <= 0:
        raise ValueError("threads must be positive")

    trimmed_fastq_dir = config.resolved_trimmed_fastq_dir()
    quant_dir = config.resolved_salmon_quant_dir()

    if build_index_if_missing:
        index_dir = _ensure_salmon_index(config, executable=executable, threads=threads)
    else:
        index_dir = config.resolved_salmon_index_dir()
        if not index_dir.exists():
            transcriptome_fasta = config.resolved_salmon_transcriptome_fasta()
            transcriptome_hint = (
                f" A transcriptome FASTA is available at {transcriptome_fasta}; rerun with --build-salmon-index to create the index automatically."
                if transcriptome_fasta is not None
                else " Provide --salmon-index-dir for a prebuilt index, or pass --salmon-transcriptome-fasta together with --build-salmon-index."
            )
            raise FileNotFoundError(f"Missing Salmon index directory: {index_dir}.{transcriptome_hint}")

    jobs = collect_salmon_inputs(trimmed_fastq_dir, quant_dir)
    quant_dir.mkdir(parents=True, exist_ok=True)

    _log(f"Trimmed FASTQ input directory: {trimmed_fastq_dir}")
    _log(f"Salmon index directory: {index_dir}")
    _log(f"Salmon output directory: {quant_dir}")
    _log(f"Input datasets selected: {len(jobs)}")

    runner = ExternalToolRunner(executable=executable, display_name="Salmon", log=_log)
    for index, job in enumerate(jobs, start=1):
        job.output_dir.mkdir(parents=True, exist_ok=True)
        _log(f"[{index}/{len(jobs)}] Quantifying {job.name}")
        runner.run(
            _build_args(job, runner=runner, index_dir=index_dir, libtype=libtype, threads=threads),
            missing_message="salmon not found in PATH",
        )

    _log("Done")
    return quant_dir
