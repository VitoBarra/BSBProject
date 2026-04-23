from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from DataSourcer.datasets import GSE103001_PROFILE, PROFILES, DatasetProfile

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_ROOT = PROJECT_ROOT / "data"


@dataclass(slots=True)
class QuantificationConfig:
    dataset_key: str = GSE103001_PROFILE.key
    trimmed_fastq_dir: Path | None = None
    salmon_index_dir: Path | None = None
    salmon_quant_dir: Path | None = None
    salmon_transcriptome_fasta: Path | None = None
    profile: DatasetProfile = field(init=False)

    def __post_init__(self) -> None:
        try:
            self.profile = PROFILES[self.dataset_key.lower()]
        except KeyError as exc:
            supported = ", ".join(sorted(PROFILES))
            raise ValueError(f"Unsupported dataset '{self.dataset_key}'. Supported datasets: {supported}") from exc

    def dataset_root(self) -> Path:
        return DATA_ROOT / self.profile.accession

    def reference_dir(self) -> Path:
        return self.dataset_root() / "reference"

    def resolved_trimmed_fastq_dir(self) -> Path:
        return self.trimmed_fastq_dir or (self.dataset_root() / "trimmed_fastq")

    def resolved_salmon_index_dir(self) -> Path:
        return self.salmon_index_dir or (self.reference_dir() / "salmon_index")

    def resolved_salmon_quant_dir(self) -> Path:
        return self.salmon_quant_dir or (self.dataset_root() / "quant" / "salmon")

    def resolved_salmon_transcriptome_fasta(self) -> Path | None:
        if self.salmon_transcriptome_fasta is not None:
            return self.salmon_transcriptome_fasta

        for filename in (
            "transcriptome.fa.gz",
            "transcriptome.fasta.gz",
            "transcriptome.fa",
            "transcriptome.fasta",
            "transcripts.fa.gz",
            "transcripts.fasta.gz",
            "transcripts.fa",
            "transcripts.fasta",
        ):
            candidate = self.reference_dir() / filename
            if candidate.exists():
                return candidate
        return None


from .salmon_runner import run_salmon
