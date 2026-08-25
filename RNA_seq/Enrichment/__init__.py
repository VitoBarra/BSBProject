from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .. import DEFAULT_CONFIG, RNASeqConfig


@dataclass(slots=True)
class EnrichmentConfig:
    workflow: RNASeqConfig = DEFAULT_CONFIG
    padj_cutoff: float = 0.05
    lfc_cutoff: float = 0.0

    def __post_init__(self) -> None:
        if not 0 < self.padj_cutoff <= 1:
            raise ValueError("padj_cutoff must be in the interval (0, 1]")
        if self.lfc_cutoff < 0:
            raise ValueError("lfc_cutoff must be non-negative")

    def resolved_de_results_path(self) -> Path:
        from ..DifferentialExpression import DifferentialExpressionConfig

        return DifferentialExpressionConfig(workflow=self.workflow).resolved_deseq2_results_path()

    def resolved_enrichment_dir(self) -> Path:
        return self.workflow.paths.enrichment

    def resolved_go_enrichment_script_path(self) -> Path:
        script_path = Path(__file__).resolve().parent / "scripts" / "run_go_enrichment.R"
        if not script_path.is_file():
            raise FileNotFoundError(f"Missing GO enrichment script: {script_path}")
        return script_path

    def resolved_selected_genes_path(self) -> Path:
        return self.resolved_enrichment_dir() / "selected_genes_ensembl.txt"

    def resolved_universe_genes_path(self) -> Path:
        return self.resolved_enrichment_dir() / "universe_genes_ensembl.txt"

    def resolved_go_all_results_path(self) -> Path:
        return self.resolved_enrichment_dir() / "go_overrepresentation_all.csv"


from .go_runner import generate_go_enrichment_plots, run_go_enrichment_analysis
