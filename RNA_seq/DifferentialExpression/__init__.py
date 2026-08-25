from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .. import DEFAULT_CONFIG, RNASeqConfig

DESEQ2_RESULT_COLUMNS = (
    "gene_id",
    "gene_symbol",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
)
DESEQ2_NUMERIC_RESULT_COLUMNS = DESEQ2_RESULT_COLUMNS[2:]


@dataclass(slots=True)
class DifferentialExpressionConfig:
    workflow: RNASeqConfig = DEFAULT_CONFIG
    min_count: int = 10
    min_samples: int = 2

    def __post_init__(self) -> None:
        if self.min_count < 0:
            raise ValueError("min_count must be non-negative")
        if self.min_samples <= 0:
            raise ValueError("min_samples must be positive")

    def de_root(self) -> Path:
        return self.workflow.paths.de_root

    def resolved_metadata_path(self) -> Path:
        return self.workflow.resolved_metadata_path()

    def resolved_salmon_quant_dir(self) -> Path:
        return self.workflow.paths.salmon_quant

    def resolved_transcriptome_fasta_path(self) -> Path:
        for filename in ("transcriptome.fa.gz", "transcriptome.fasta.gz", "transcriptome.fa", "transcriptome.fasta"):
            candidate = self.workflow.paths.reference_root / filename
            if candidate.exists():
                return candidate
        raise FileNotFoundError(f"No transcriptome FASTA found under {self.workflow.paths.reference_root}")

    def resolved_output_dir(self) -> Path:
        return self.de_root()

    def resolved_deseq2_sample_sheet_path(self) -> Path:
        return self.de_root() / "deseq2_sample_sheet.tsv"

    def resolved_transcript_to_gene_map_path(self) -> Path:
        return self.de_root() / "transcript_to_gene_map.tsv"

    def resolved_deseq2_script_path(self) -> Path:
        script_path = Path(__file__).resolve().parent / "scripts" / "run_deseq2.R"
        if not script_path.is_file():
            raise FileNotFoundError(f"Missing DESeq2 script: {script_path}")
        return script_path

    def resolved_de_results_dir(self) -> Path:
        return self.de_root() / "results"

    def resolved_deseq2_results_path(self) -> Path:
        return self.resolved_de_results_dir() / "deseq2_all_genes.csv"

    def resolved_vst_counts_path(self) -> Path:
        return self.resolved_de_results_dir() / "vst_counts.csv"

    def resolved_deseq2_summary_path(self) -> Path:
        return self.resolved_de_results_dir() / "deseq2_summary.txt"


from .prepare_deseq2_inputs import prepare_deseq2_inputs
from .deseq2_runner import generate_deseq2_plots, run_deseq2_analysis
