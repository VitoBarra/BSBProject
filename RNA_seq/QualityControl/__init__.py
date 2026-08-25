from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .. import DEFAULT_CONFIG, RNASeqConfig


@dataclass(slots=True)
class QualityControlConfig:
    workflow: RNASeqConfig = DEFAULT_CONFIG
    fastqc_threads: int = 2
    fastp_threads: int = 2

    def __post_init__(self) -> None:
        if self.fastqc_threads <= 0:
            raise ValueError("fastqc_threads must be positive")
        if self.fastp_threads <= 0:
            raise ValueError("fastp_threads must be positive")

    def resolved_metadata_path(self) -> Path:
        return self.workflow.resolved_metadata_path()

    def resolved_fastq_dir(self) -> Path:
        return self.workflow.paths.raw_fastq

    def resolved_trimmed_fastq_dir(self) -> Path:
        return self.workflow.paths.trimmed_fastq

    def resolved_fastp_report_out(self) -> Path:
        return self.workflow.paths.fastp_reports

    def resolved_fastqc_report_out(self) -> Path:
        return self.workflow.paths.raw_fastqc_reports

    def resolved_fastqc_trimmed_report_out(self) -> Path:
        return self.workflow.paths.trimmed_fastqc_reports

    def resolved_multiqc_raw_report_out(self) -> Path:
        return self.workflow.paths.raw_multiqc_reports

    def resolved_multiqc_trimmed_report_out(self) -> Path:
        return self.workflow.paths.trimmed_multiqc_reports


from .fastp_runner import run_fastp
from .fastqc_runner import run_fastqc
from .multiqc_runner import run_multiqc
