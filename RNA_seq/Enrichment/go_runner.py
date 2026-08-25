from __future__ import annotations

import csv
import logging
import math
from pathlib import Path

from ExternalTools import ExternalToolRunner

from . import EnrichmentConfig
from .go_results import GO_RESULT_COLUMNS, prepare_go_inputs

LOGGER = logging.getLogger("enrichment")

def _completed_results(result_path: Path) -> bool:
    if not result_path.is_file() or result_path.stat().st_size == 0:
        return False
    try:
        with result_path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            if not set(GO_RESULT_COLUMNS).issubset(reader.fieldnames or ()):
                return False
            for row in reader:
                if not row["ID"].strip() or not row["ontology"].strip():
                    return False
                if not math.isfinite(float(row["p.adjust"])) or int(row["Count"]) < 0:
                    return False
    except (OSError, csv.Error, TypeError, ValueError):
        return False
    return True


def run_go_enrichment_analysis(
    config: EnrichmentConfig,
) -> None:
    de_results_path = config.resolved_de_results_path()
    output_dir = config.resolved_enrichment_dir()
    go_enrichment_script_path = config.resolved_go_enrichment_script_path()
    selected_genes_path = config.resolved_selected_genes_path()
    universe_genes_path = config.resolved_universe_genes_path()
    all_results_path = config.resolved_go_all_results_path()

    if not de_results_path.exists():
        raise FileNotFoundError(f"Missing DESeq2 result table: {de_results_path}")

    output_dir.mkdir(parents=True, exist_ok=True)
    LOGGER.info("DESeq2 result table: %s", de_results_path)
    LOGGER.info("Output directory: %s", output_dir)

    universe_count, selected_count = prepare_go_inputs(
        de_results_path=de_results_path,
        selected_genes_path=selected_genes_path,
        universe_genes_path=universe_genes_path,
        padj_cutoff=config.padj_cutoff,
        lfc_cutoff=config.lfc_cutoff,
    )
    LOGGER.info("GO universe: %d genes", universe_count)
    LOGGER.info("Selected DE genes: %d", selected_count)

    if _completed_results(all_results_path):
        LOGGER.info("Complete GO enrichment results already exist; skipping the R analysis")
        return

    runner = ExternalToolRunner(executable="Rscript", display_name="Rscript", logger=LOGGER)
    runner.run(
        [
            runner.path_arg(go_enrichment_script_path),
            "--selected-genes",
            runner.path_arg(selected_genes_path),
            "--universe-genes",
            runner.path_arg(universe_genes_path),
            "--results-path",
            runner.path_arg(all_results_path),
            "--de-results-path",
            runner.path_arg(de_results_path),
            "--padj-cutoff",
            str(config.padj_cutoff),
            "--lfc-cutoff",
            str(config.lfc_cutoff),
        ],
        missing_message="Rscript not found in PATH",
    )
    if not _completed_results(all_results_path):
        raise RuntimeError(f"GO enrichment completed without producing valid outputs in {output_dir}")

    LOGGER.info("Done")

def generate_go_enrichment_plots(config: EnrichmentConfig) -> None:
    LOGGER.info("Generating GO result tables, summary, and plot with R")
    runner = ExternalToolRunner(executable="Rscript", display_name="Rscript", logger=LOGGER)
    runner.run(
        [
            runner.path_arg(config.resolved_go_enrichment_script_path()),
            "--selected-genes",
            runner.path_arg(config.resolved_selected_genes_path()),
            "--universe-genes",
            runner.path_arg(config.resolved_universe_genes_path()),
            "--results-path",
            runner.path_arg(config.resolved_go_all_results_path()),
            "--de-results-path",
            runner.path_arg(config.resolved_de_results_path()),
            "--padj-cutoff",
            str(config.padj_cutoff),
            "--lfc-cutoff",
            str(config.lfc_cutoff),
            "--plots-only",
        ],
        missing_message="Rscript not found in PATH",
    )
    LOGGER.info("Done")
