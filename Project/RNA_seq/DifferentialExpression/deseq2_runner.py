from __future__ import annotations

from pathlib import Path

from ExternalTools import ExternalToolRunner
from Log.log_util import log

from . import DifferentialExpressionConfig
from .de_results import generate_de_outputs

LOG_PREFIX = "deseq2"
SCRIPT_PATH = Path(__file__).resolve().parent / "scripts" / "run_deseq2.R"


def _log(message: str) -> None:
    log(message, LOG_PREFIX)


def run_deseq2_analysis(
    config: DifferentialExpressionConfig,
    executable: str = "Rscript",
) -> Path:
    sample_table_path = config.resolved_sample_table_path()
    tx2gene_path = config.resolved_tx2gene_path()
    output_dir = config.resolved_de_results_dir()

    if not sample_table_path.exists():
        raise FileNotFoundError(f"Missing sample table: {sample_table_path}")
    if not tx2gene_path.exists():
        raise FileNotFoundError(f"Missing transcript-to-gene table: {tx2gene_path}")
    if not SCRIPT_PATH.exists():
        raise FileNotFoundError(f"Missing DESeq2 script: {SCRIPT_PATH}")

    output_dir.mkdir(parents=True, exist_ok=True)
    _log(f"Sample table: {sample_table_path}")
    _log(f"Transcript-to-gene table: {tx2gene_path}")
    _log(f"Output directory: {output_dir}")
    _log("Design formula: ~ patient + condition")

    runner = ExternalToolRunner(executable=executable, display_name="Rscript", log=_log)
    runner.run(
        [
            runner.path_arg(SCRIPT_PATH),
            "--samples",
            runner.path_arg(sample_table_path),
            "--tx2gene",
            runner.path_arg(tx2gene_path),
            "--outdir",
            runner.path_arg(output_dir),
            "--min-count",
            str(config.min_count),
            "--min-samples",
            str(config.min_samples),
        ],
        missing_message="Rscript not found in PATH",
    )

    _log("Done")
    return output_dir


def generate_deseq2_plots(config: DifferentialExpressionConfig) -> Path:
    output_dir = config.resolved_de_results_dir()
    sample_table_path = config.resolved_sample_table_path()

    _log("Generating DE result tables and plots with Python")
    generate_de_outputs(
        results_path=output_dir / "deseq2_all_genes.csv",
        normalized_counts_path=output_dir / "normalized_counts.csv",
        vst_counts_path=output_dir / "vst_counts.csv",
        samples_path=sample_table_path,
    )

    _log("Done")
    return output_dir


def run_deseq2(
    config: DifferentialExpressionConfig,
    executable: str = "Rscript",
) -> Path:
    run_deseq2_analysis(config, executable=executable)
    return generate_deseq2_plots(config)
