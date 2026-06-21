from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from DataSourcer.datasets import GSE103001_PROFILE, PROFILES, DatasetProfile

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_ROOT = PROJECT_ROOT / "data"


@dataclass(slots=True)
class EnrichmentConfig:
    dataset_key: str = GSE103001_PROFILE.key
    de_results_path: Path | None = None
    enrichment_dir: Path | None = None
    padj_cutoff: float = 0.05
    lfc_cutoff: float = 0.0
    profile: DatasetProfile = field(init=False)

    def __post_init__(self) -> None:
        try:
            self.profile = PROFILES[self.dataset_key.lower()]
        except KeyError as exc:
            supported = ", ".join(sorted(PROFILES))
            raise ValueError(f"Unsupported dataset '{self.dataset_key}'. Supported datasets: {supported}") from exc
        if not 0 < self.padj_cutoff <= 1:
            raise ValueError("padj_cutoff must be in the interval (0, 1]")
        if self.lfc_cutoff < 0:
            raise ValueError("lfc_cutoff must be non-negative")

    def dataset_root(self) -> Path:
        return DATA_ROOT / self.profile.accession

    def resolved_de_results_path(self) -> Path:
        return self.de_results_path or (self.dataset_root() / "de" / "results" / "deseq2_all_genes.csv")

    def resolved_enrichment_dir(self) -> Path:
        return self.enrichment_dir or (self.dataset_root() / "enrichment")


from .go_runner import run_go_enrichment
