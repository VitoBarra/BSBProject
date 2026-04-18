from __future__ import annotations

import logging

from utils.runner import ExternalToolRunner

from . import EnrichmentConfig


LOGGER = logging.getLogger("enrichment")


def _run_r_script(config: EnrichmentConfig, plots_only: bool = False) -> None:
    runner = ExternalToolRunner(executable="Rscript", display_name="Rscript", logger=LOGGER)
    arguments = [
        runner.path_arg(config.resolved_go_enrichment_script_path()),
        "--results-path",
        runner.path_arg(config.resolved_go_all_results_path()),
        "--de-results-path",
        runner.path_arg(config.resolved_de_results_path()),
        "--padj-cutoff",
        str(config.padj_cutoff),
        "--lfc-cutoff",
        str(config.lfc_cutoff),
    ]
    if plots_only:
        arguments.append("--plots-only")
    runner.run(arguments, missing_message="Rscript not found in PATH")


def run_go_enrichment_analysis(
    config: EnrichmentConfig,
) -> None:
    de_results_path = config.resolved_de_results_path()
    output_dir = config.resolved_enrichment_dir()
    all_results_path = config.resolved_go_all_results_path()

    if not de_results_path.exists():
        raise FileNotFoundError(f"Missing DESeq2 result table: {de_results_path}")

    output_dir.mkdir(parents=True, exist_ok=True)
    LOGGER.info("DESeq2 result table: %s", de_results_path)
    LOGGER.info("Output directory: %s", output_dir)

    if all_results_path.is_file() and all_results_path.stat().st_size > 0:
        LOGGER.info("GO enrichment results already exist; refreshing derived outputs")
        _run_r_script(config, plots_only=True)
        return

    _run_r_script(config)
    if not all_results_path.is_file() or all_results_path.stat().st_size == 0:
        raise RuntimeError(f"GO enrichment completed without producing results in {output_dir}")

    LOGGER.info("Done")
