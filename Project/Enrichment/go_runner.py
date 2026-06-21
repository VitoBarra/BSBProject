from __future__ import annotations

from pathlib import Path

from ExternalTools import ExternalToolRunner
from Log.log_util import log

from . import EnrichmentConfig

LOG_PREFIX = "enrichment"
SCRIPT_PATH = Path(__file__).resolve().parent / "scripts" / "run_go_enrichment.R"


def _log(message: str) -> None:
    log(message, LOG_PREFIX)


def run_go_enrichment(
    config: EnrichmentConfig,
    executable: str = "Rscript",
) -> Path:
    de_results_path = config.resolved_de_results_path()
    output_dir = config.resolved_enrichment_dir()

    if not de_results_path.exists():
        raise FileNotFoundError(f"Missing DESeq2 result table: {de_results_path}")
    if not SCRIPT_PATH.exists():
        raise FileNotFoundError(f"Missing enrichment script: {SCRIPT_PATH}")

    output_dir.mkdir(parents=True, exist_ok=True)
    _log(f"DESeq2 result table: {de_results_path}")
    _log(f"Output directory: {output_dir}")

    runner = ExternalToolRunner(executable=executable, display_name="Rscript", log=_log)
    runner.run(
        [
            runner.path_arg(SCRIPT_PATH),
            "--de-results",
            runner.path_arg(de_results_path),
            "--outdir",
            runner.path_arg(output_dir),
            "--padj-cutoff",
            str(config.padj_cutoff),
            "--lfc-cutoff",
            str(config.lfc_cutoff),
        ],
        missing_message="Rscript not found in PATH",
    )

    _log("Done")
    return output_dir
