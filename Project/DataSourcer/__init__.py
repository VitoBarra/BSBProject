from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

from .datasets import (
    DatasetProfile,
    GSE103001_PROFILE,
    GSE210787_PROFILE,
    GSE31210_PROFILE,
    PROFILES,
)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATASOURCER_ROOT = Path(__file__).resolve().parent

DATA_DIR_NAME = "data"
DATA_ROOT = PROJECT_ROOT / DATA_DIR_NAME

MetadataBuilder = Callable[["DataSourceConfig"], Path]


@dataclass(slots=True)
class DataSourceConfig:
    dataset_key: str = GSE103001_PROFILE.key
    num_pairs: int = 4
    soft_path: Path | None = None
    metadata_path: Path | None = None
    fastq_dest: Path | None = None
    download_workers: int = 3
    profile: DatasetProfile = field(init=False)

    def __post_init__(self) -> None:
        if self.download_workers <= 0:
            raise ValueError("download_workers must be positive")
        try:
            self.profile = PROFILES[self.dataset_key.lower()]
        except KeyError as exc:
            supported = ", ".join(sorted(PROFILES))
            raise ValueError(f"Unsupported dataset '{self.dataset_key}'. Supported datasets: {supported}") from exc

    def dataset_root(self) -> Path:
        return DATA_ROOT / self.profile.accession

    def resolved_soft_path(self) -> Path:
        return self.soft_path or (self.dataset_root() / self.profile.soft_filename())

    def resolved_metadata_path(self) -> Path:
        return self.metadata_path or (self.dataset_root() / self.profile.metadata_filename(self.num_pairs))

    def resolved_fastq_dest(self) -> Path:
        return self.fastq_dest or (self.dataset_root() / "raw_fastq")


def build_metadata_table(config: DataSourceConfig) -> Path:
    builder = METADATA_BUILDERS[config.profile.key]
    return builder(config)


from .GSE103001_metadata_extractor import build_metadata_table as build_gse103001_metadata_table
from .GSE210787_metadata_extractor import build_metadata_table as build_gse210787_metadata_table
from .GSE31210_metadata_extractor import build_metadata_table as build_gse31210_metadata_table

METADATA_BUILDERS: dict[str, MetadataBuilder] = {
    GSE103001_PROFILE.key: build_gse103001_metadata_table,
    GSE210787_PROFILE.key: build_gse210787_metadata_table,
    GSE31210_PROFILE.key: build_gse31210_metadata_table,
}
