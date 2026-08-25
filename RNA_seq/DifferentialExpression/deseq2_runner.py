from __future__ import annotations

import csv
import logging
from pathlib import Path

from ExternalTools import ExternalToolRunner

from . import DESEQ2_RESULT_COLUMNS, DifferentialExpressionConfig

LOGGER = logging.getLogger("deseq2")
def _csv_has_rows(path: Path, required_columns: set[str]) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            return required_columns.issubset(reader.fieldnames or ()) and next(reader, None) is not None
    except (OSError, csv.Error):
        return False


def _completed_results(
    results_path: Path,
    vst_counts_path: Path,
    summary_path: Path,
    deseq2_sample_sheet_path: Path | None = None,
) -> bool:
    if not _csv_has_rows(results_path, set(DESEQ2_RESULT_COLUMNS)):
        return False
    if not _csv_has_rows(vst_counts_path, {"gene_id", "gene_symbol"}):
        return False
    if not summary_path.is_file() or summary_path.stat().st_size == 0:
        return False
    if deseq2_sample_sheet_path is not None:
        try:
            with deseq2_sample_sheet_path.open("r", encoding="utf-8-sig", newline="") as handle:
                sample_names = [row["sample_name"] for row in csv.DictReader(handle, delimiter="\t")]
            with vst_counts_path.open("r", encoding="utf-8-sig", newline="") as handle:
                vst_columns = list(csv.DictReader(handle).fieldnames or ())
        except (OSError, csv.Error, KeyError):
            return False
        if vst_columns[2:] != sample_names:
            return False
    return True


def run_deseq2_analysis(
    config: DifferentialExpressionConfig,
) -> None:
    deseq2_sample_sheet_path = config.resolved_deseq2_sample_sheet_path()
    transcript_to_gene_map_path = config.resolved_transcript_to_gene_map_path()
    deseq2_script_path = config.resolved_deseq2_script_path()
    output_dir = config.resolved_de_results_dir()
    results_path = config.resolved_deseq2_results_path()
    vst_counts_path = config.resolved_vst_counts_path()
    summary_path = config.resolved_deseq2_summary_path()


    # check that the required input files exist
    if not deseq2_sample_sheet_path.exists():
        raise FileNotFoundError(f"Missing DESeq2 sample sheet: {deseq2_sample_sheet_path}")
    if not transcript_to_gene_map_path.exists():
        raise FileNotFoundError(f"Missing transcript-to-gene map: {transcript_to_gene_map_path}")
    output_dir.mkdir(parents=True, exist_ok=True)

    # check if the DESeq2 results already exist and are complete
    if _completed_results(results_path, vst_counts_path, summary_path, deseq2_sample_sheet_path):
        LOGGER.info("Complete DESeq2 results already exist; skipping the R analysis")
        return

    # run the DESeq2 R script with the provided sample sheet and transcript-to-gene map
    LOGGER.info("DESeq2 sample sheet: %s", deseq2_sample_sheet_path)
    LOGGER.info("Transcript-to-gene map: %s", transcript_to_gene_map_path)
    LOGGER.info("Output directory: %s", output_dir)
    LOGGER.info("Design formula: ~ patient + condition")

    runner = ExternalToolRunner(executable="Rscript", display_name="Rscript", logger=LOGGER)
    runner.run(
        [
            runner.path_arg(deseq2_script_path),
            "--samples",
            runner.path_arg(deseq2_sample_sheet_path),
            "--transcript-to-gene-map",
            runner.path_arg(transcript_to_gene_map_path),
            "--results-path",
            runner.path_arg(results_path),
            "--vst-counts-path",
            runner.path_arg(vst_counts_path),
            "--summary-path",
            runner.path_arg(summary_path),
            "--min-count",
            str(config.min_count),
            "--min-samples",
            str(config.min_samples),
        ],
        missing_message="Rscript not found in PATH",
    )
    # check if the DESeq2 results are complete after running the R script
    if not _completed_results(results_path, vst_counts_path, summary_path, deseq2_sample_sheet_path):
        raise RuntimeError(f"DESeq2 completed without producing valid outputs in {output_dir}")
    LOGGER.info("Done")


def generate_deseq2_plots(config: DifferentialExpressionConfig) -> None:
    LOGGER.info("Generating DE result tables and plots with R")
    runner = ExternalToolRunner(executable="Rscript", display_name="Rscript", logger=LOGGER)
    runner.run(
        [
            runner.path_arg(config.resolved_deseq2_script_path()),
            "--samples",
            runner.path_arg(config.resolved_deseq2_sample_sheet_path()),
            "--transcript-to-gene-map",
            runner.path_arg(config.resolved_transcript_to_gene_map_path()),
            "--results-path",
            runner.path_arg(config.resolved_deseq2_results_path()),
            "--vst-counts-path",
            runner.path_arg(config.resolved_vst_counts_path()),
            "--summary-path",
            runner.path_arg(config.resolved_deseq2_summary_path()),
            "--plots-only",
        ],
        missing_message="Rscript not found in PATH",
    )
    LOGGER.info("Done")
