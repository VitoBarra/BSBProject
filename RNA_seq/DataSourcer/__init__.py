from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .. import DEFAULT_CONFIG, RNASeqConfig
@dataclass(slots=True)
class DataSourceConfig:
    workflow: RNASeqConfig = DEFAULT_CONFIG
    download_workers: int = 3

    def __post_init__(self) -> None:
        if self.download_workers <= 0:
            raise ValueError("download_workers must be positive")

    @property
    def num_pairs(self) -> int:
        return self.workflow.num_pairs

    @property
    def paths(self):
        return self.workflow.paths

    def resolved_soft_path(self) -> Path:
        return self.paths.dataset_root / self.workflow.datasetInfo.soft_filename

    def resolved_metadata_path(self) -> Path:
        return self.workflow.resolved_metadata_path()

    def resolved_fastq_dest(self) -> Path:
        return self.paths.raw_fastq

from .GSE103001_metadata_extractor import ensure_metadata_table
from .download_fastq import download_fastq
from .download_reference_transcriptome import download_reference_transcriptome
