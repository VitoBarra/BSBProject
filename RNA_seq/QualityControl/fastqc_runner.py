from __future__ import annotations

import logging
import shutil
import zipfile
from dataclasses import dataclass
from pathlib import Path

from ExternalTools import ExternalToolRunner
from ..fastq_manifest import FASTQ_SUFFIX, load_paired_fastqs

from . import QualityControlConfig

LOGGER = logging.getLogger("fastqc")


@dataclass(frozen=True, slots=True)
class FastQCJob:
    input_fastq: Path
    output_dir: Path

    @property
    def output_stem(self) -> str:
        if not self.input_fastq.name.endswith(FASTQ_SUFFIX):
            raise ValueError(f"Expected a {FASTQ_SUFFIX} input: {self.input_fastq}")
        return f"{self.input_fastq.name.removesuffix(FASTQ_SUFFIX)}_fastqc"

    @property
    def html_report(self) -> Path:
        return self.output_dir / f"{self.output_stem}.html"

    @property
    def zip_report(self) -> Path:
        return self.output_dir / f"{self.output_stem}.zip"

    def validate_outputs(self) -> None:
        for path in (self.html_report, self.zip_report):
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError(f"Missing or empty FastQC output: {path}")

        try:
            with zipfile.ZipFile(self.zip_report) as archive:
                names = archive.namelist()
                if archive.testzip() is not None:
                    raise ValueError(f"Corrupt FastQC ZIP output: {self.zip_report}")
        except (OSError, zipfile.BadZipFile) as exc:
            raise ValueError(f"Invalid FastQC ZIP output: {self.zip_report}: {exc}") from exc
        if not any(name.endswith("/fastqc_data.txt") for name in names):
            raise ValueError(f"FastQC ZIP has no fastqc_data.txt: {self.zip_report}")
        if not any(name.endswith("/summary.txt") for name in names):
            raise ValueError(f"FastQC ZIP has no summary.txt: {self.zip_report}")



def collect_fastqc_inputs(metadata_path: Path, fastq_dir: Path) -> list[Path]:
    return [read for pair in load_paired_fastqs(metadata_path, fastq_dir) for read in (pair.read_1, pair.read_2)]


def run_fastqc(
    config: QualityControlConfig,
    fastq_dir: Path,
    out_dir: Path,
) -> None:
    """Run FastQC for each FASTQ file in ``fastq_dir``.

    Reports are written to ``out_dir``. Existing valid reports are reused;
    missing, empty, or malformed reports are rebuilt in a temporary directory
    and then replaced.
    """

    inputs = collect_fastqc_inputs(config.resolved_metadata_path(), fastq_dir)
    jobs = [FastQCJob(path, out_dir) for path in inputs]
    out_dir.mkdir(parents=True, exist_ok=True)

    LOGGER.info("FASTQ directory: %s", fastq_dir)
    LOGGER.info("FastQC output directory: %s", out_dir)
    LOGGER.info("Input FASTQ files selected: %d", len(inputs))

    jobs_to_rerun: list[FastQCJob] = []
    for job in jobs:
        try:
            job.validate_outputs()  # Check that the output files are valid and complete.
            LOGGER.info("Skipping completed FastQC: %s", job.input_fastq.name)
        except ValueError as exc:
            LOGGER.warning("FastQC output for %s requires rebuilding: %s", job.input_fastq.name, exc)
            jobs_to_rerun.append(job)

    if not jobs_to_rerun:
        LOGGER.info("All FastQC outputs are complete")
        return

    rebuild_dir = out_dir / ".fastqc-rebuild"
    if rebuild_dir.exists():
        shutil.rmtree(rebuild_dir)
    rebuild_dir.mkdir(parents=True)
    runner = ExternalToolRunner(executable="fastqc", display_name="FastQC", logger=LOGGER)
    runner.run(
        [
            # Run this many FastQC worker threads.
            "--threads",
            str(config.fastqc_threads),
            # Write temporary HTML and ZIP reports to the rebuild directory.
            "--outdir",
            runner.path_arg(rebuild_dir),
            # Analyze every FASTQ whose report must be rebuilt.
            *(runner.path_arg(job.input_fastq) for job in jobs_to_rerun),
        ],
        missing_message="FastQC not found in PATH",
    )

    for job in jobs_to_rerun:
        rebuilt = FastQCJob(job.input_fastq, rebuild_dir)
        # Validate and replace the rebuilt outputs.
        rebuilt.validate_outputs()
        rebuilt.html_report.replace(job.html_report)
        rebuilt.zip_report.replace(job.zip_report)

    shutil.rmtree(rebuild_dir)

    LOGGER.info("Done")
