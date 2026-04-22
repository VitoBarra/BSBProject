from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from DataSourcer.datasets import GSE103001_PROFILE, PROFILES, DatasetProfile

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_ROOT = PROJECT_ROOT / "data"


@dataclass(slots=True)
class QualityControlConfig:
    dataset_key: str = GSE103001_PROFILE.key
    fastq_dir: Path | None = None
    fastqc_report_out: Path | None = None
    multiqc_report_out: Path | None = None
    profile: DatasetProfile = field(init=False)

    def __post_init__(self) -> None:
        try:
            self.profile = PROFILES[self.dataset_key.lower()]
        except KeyError as exc:
            supported = ", ".join(sorted(PROFILES))
            raise ValueError(f"Unsupported dataset '{self.dataset_key}'. Supported datasets: {supported}") from exc

    def dataset_root(self) -> Path:
        return DATA_ROOT / self.profile.accession

    def resolved_fastq_dir(self) -> Path:
        return self.fastq_dir or (self.dataset_root() / "raw_fastq")

    def resolved_fastqc_report_out(self) -> Path:
        return self.fastqc_report_out or (self.dataset_root() / "qc" / "fastqc")

    def resolved_multiqc_report_out(self) -> Path:
        return self.multiqc_report_out or (self.dataset_root() / "qc" / "multiqc")

    def resolved_qc_root(self) -> Path:
        return self.dataset_root() / "qc"


from .fastqc_runner import run_fastqc
from .multiqc_runner import run_multiqc
