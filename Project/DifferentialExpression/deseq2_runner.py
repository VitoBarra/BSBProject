from __future__ import annotations

from pathlib import Path

from ExternalTools import ExternalToolRunner
from Log.log_util import log

from . import DifferentialExpressionConfig

LOG_PREFIX = "deseq2"
SCRIPT_PATH = Path(__file__).resolve().parent / "scripts" / "run_deseq2.R"


def _log(message: str) -> None:
    log(message, LOG_PREFIX)


def run_deseq2(
    config: DifferentialExpressionConfig,
    executable: str = "Rscript",
) -> Path:
    counts_path = config.resolved_gene_counts_path()
    sample_table_path = config.resolved_sample_table_path()
    output_dir = config.resolved_de_results_dir()

    if not counts_path.exists():
        raise FileNotFoundError(f"Missing gene count matrix: {counts_path}")
    if not sample_table_path.exists():
        raise FileNotFoundError(f"Missing sample table: {sample_table_path}")
    if not SCRIPT_PATH.exists():
        raise FileNotFoundError(f"Missing DESeq2 script: {SCRIPT_PATH}")

    output_dir.mkdir(parents=True, exist_ok=True)
    _log(f"Gene count matrix: {counts_path}")
    _log(f"Sample table: {sample_table_path}")
    _log(f"Output directory: {output_dir}")
    _log("Design formula: ~ patient + condition")

    runner = ExternalToolRunner(executable=executable, display_name="Rscript", log=_log)
    runner.run(
        [
            runner.path_arg(SCRIPT_PATH),
            "--counts",
            runner.path_arg(counts_path),
            "--samples",
            runner.path_arg(sample_table_path),
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
