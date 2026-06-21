from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from DataSourcer.datasets import GSE103001_PROFILE, PROFILES, DatasetProfile

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_ROOT = PROJECT_ROOT / "data"


@dataclass(slots=True)
class DifferentialExpressionConfig:
    dataset_key: str = GSE103001_PROFILE.key
    gene_counts_path: Path | None = None
    sample_table_path: Path | None = None
    de_results_dir: Path | None = None
    min_count: int = 10
    min_samples: int = 2
    profile: DatasetProfile = field(init=False)

    def __post_init__(self) -> None:
        try:
            self.profile = PROFILES[self.dataset_key.lower()]
        except KeyError as exc:
            supported = ", ".join(sorted(PROFILES))
            raise ValueError(f"Unsupported dataset '{self.dataset_key}'. Supported datasets: {supported}") from exc
        if self.min_count < 0:
            raise ValueError("min_count must be non-negative")
        if self.min_samples <= 0:
            raise ValueError("min_samples must be positive")

    def dataset_root(self) -> Path:
        return DATA_ROOT / self.profile.accession

    def de_root(self) -> Path:
        return self.dataset_root() / "de"

    def resolved_gene_counts_path(self) -> Path:
        return self.gene_counts_path or (self.de_root() / "salmon_gene_counts.tsv")

    def resolved_sample_table_path(self) -> Path:
        return self.sample_table_path or (self.de_root() / "sample_table.tsv")

    def resolved_de_results_dir(self) -> Path:
        return self.de_results_dir or (self.de_root() / "results")


@dataclass(slots=True)
class SalmonGeneAggregationConfig:
    dataset_key: str = GSE103001_PROFILE.key
    metadata_path: Path | None = None
    salmon_quant_dir: Path | None = None
    transcriptome_fasta_path: Path | None = None
    output_dir: Path | None = None
    tx2gene_path: Path | None = None
    profile: DatasetProfile = field(init=False)

    def __post_init__(self) -> None:
        try:
            self.profile = PROFILES[self.dataset_key.lower()]
        except KeyError as exc:
            supported = ", ".join(sorted(PROFILES))
            raise ValueError(f"Unsupported dataset '{self.dataset_key}'. Supported datasets: {supported}") from exc

    def dataset_root(self) -> Path:
        return DATA_ROOT / self.profile.accession

    def resolved_metadata_path(self) -> Path:
        if self.metadata_path is not None:
            return self.metadata_path
        candidates = sorted(self.dataset_root().glob(f"{self.profile.accession}_selected_*pairs.tsv"))
        if len(candidates) != 1:
            raise FileNotFoundError(
                f"Expected exactly one selected-pairs metadata TSV under {self.dataset_root()}, found {len(candidates)}. "
                "Pass --metadata-path explicitly."
            )
        return candidates[0]

    def resolved_salmon_quant_dir(self) -> Path:
        return self.salmon_quant_dir or (self.dataset_root() / "quant" / "salmon")

    def resolved_transcriptome_fasta_path(self) -> Path:
        if self.transcriptome_fasta_path is not None:
            return self.transcriptome_fasta_path
        reference_dir = self.dataset_root() / "reference"
        for filename in ("transcriptome.fa.gz", "transcriptome.fasta.gz", "transcriptome.fa", "transcriptome.fasta"):
            candidate = reference_dir / filename
            if candidate.exists():
                return candidate
        raise FileNotFoundError(f"No transcriptome FASTA found under {reference_dir}")

    def resolved_output_dir(self) -> Path:
        return self.output_dir or (self.dataset_root() / "de")

    def resolved_tx2gene_path(self) -> Path:
        return self.tx2gene_path or (self.resolved_output_dir() / "tx2gene.tsv")


from .salmon_to_gene import build_salmon_gene_matrices
from .deseq2_runner import run_deseq2
