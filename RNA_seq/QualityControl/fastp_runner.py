from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ExternalTools import ExternalToolRunner
from ..fastq_manifest import PairedFastq, load_paired_fastqs

from . import QualityControlConfig

LOGGER = logging.getLogger("fastp")


@dataclass(slots=True, frozen=True)
class FastpJob:
    name: str
    input_fastq_1: Path
    input_fastq_2: Path
    trimmed_dir: Path
    report_dir: Path

    @classmethod
    def from_pair(cls, pair: PairedFastq, trimmed_dir: Path, report_dir: Path) -> FastpJob:
        return cls(
            name=pair.srr,
            input_fastq_1=pair.read_1,
            input_fastq_2=pair.read_2,
            trimmed_dir=trimmed_dir,
            report_dir=report_dir,
        )

    @property
    def output_fastq_1(self) -> Path:
        return self.trimmed_dir / self.input_fastq_1.name

    @property
    def output_fastq_2(self) -> Path:
        return self.trimmed_dir / self.input_fastq_2.name

    @property
    def html_report(self) -> Path:
        return self.report_dir / f"{self.name}.html"

    @property
    def json_report(self) -> Path:
        return self.report_dir / f"{self.name}.json"

    @property
    def output_paths(self) -> tuple[Path, Path, Path, Path]:
        return self.output_fastq_1, self.output_fastq_2, self.html_report, self.json_report

    def validate_outputs(self) -> dict[str, Any]:
        for path in self.output_paths:
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError(f"Missing or empty fastp output: {path}")
        try:
            report = json.loads(self.json_report.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ValueError(f"Invalid fastp JSON report: {self.json_report}: {exc}") from exc
        try:
            before_reads = int(report["summary"]["before_filtering"]["total_reads"])
            after_reads = int(report["summary"]["after_filtering"]["total_reads"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"Incomplete fastp summary: {self.json_report}") from exc
        if before_reads <= 0 or after_reads < 0 or after_reads > before_reads:
            raise ValueError(f"Invalid fastp read counts in {self.json_report}")
        return report


def run_fastp(
    config: QualityControlConfig,
) -> None:
    """Trim each paired FASTQ sample.

    Valid existing trimmed reads and reports are reused. Missing or invalid
    outputs are regenerated with the configured number of worker threads.
    """
    fastq_dir = config.resolved_fastq_dir()
    trimmed_dir = config.resolved_trimmed_fastq_dir()
    report_dir = config.resolved_fastp_report_out()

    if not fastq_dir.exists():
        raise FileNotFoundError(f"Missing FASTQ directory: {fastq_dir}")

    jobs = [
        FastpJob.from_pair(pair, trimmed_dir, report_dir)
        for pair in load_paired_fastqs(config.resolved_metadata_path(), fastq_dir)
    ]
    trimmed_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)

    LOGGER.info("FASTQ input directory: %s", fastq_dir)
    LOGGER.info("Trimmed FASTQ output directory: %s", trimmed_dir)
    LOGGER.info("fastp report directory: %s", report_dir)
    LOGGER.info("Input datasets selected: %d", len(jobs))

    jobs_to_rerun: list[FastpJob] = []
    for index, job in enumerate(jobs, start=1):
        try:
            report = job.validate_outputs()
        except ValueError as exc:
            LOGGER.warning("fastp outputs for %s require rebuilding: %s", job.name, exc)
            jobs_to_rerun.append(job)
            continue

        retained = report["summary"]["after_filtering"]["total_reads"]
        LOGGER.info("[%d/%d] Skipping completed cleaning for %s (%s reads retained)",
                    index,len(jobs),job.name,retained,)

    if not jobs_to_rerun:
        LOGGER.info("All fastp outputs are complete")
        return

    runner = ExternalToolRunner(executable="fastp", display_name="fastp", logger=LOGGER)
    for index, job in enumerate(jobs_to_rerun, start=1):
        LOGGER.info("[%d/%d] Trimming %s", index, len(jobs_to_rerun), job.name)
        runner.run(
            [
                # Run this many fastp worker threads.
                "--thread",
                str(config.fastp_threads),
                # Write the per-sample HTML report here.
                "--html",
                runner.path_arg(job.html_report),
                # Write the machine-readable per-sample report here.
                "--json",
                runner.path_arg(job.json_report),
                # Read the first mate of the paired FASTQ sample.
                "--in1",
                runner.path_arg(job.input_fastq_1),
                # Write the trimmed first mate here.
                "--out1",
                runner.path_arg(job.output_fastq_1),
                # Read the second mate of the paired FASTQ sample.
                "--in2",
                runner.path_arg(job.input_fastq_2),
                # Write the trimmed second mate here.
                "--out2",
                runner.path_arg(job.output_fastq_2),
            ],
            missing_message="fastp not found in PATH",
        )
        job.validate_outputs()

    LOGGER.info("Done")
