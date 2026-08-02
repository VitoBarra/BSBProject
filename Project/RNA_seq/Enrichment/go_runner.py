from __future__ import annotations

from pathlib import Path

from ExternalTools import ExternalToolRunner
from Log.log_util import log

from . import EnrichmentConfig
from .go_results import generate_go_outputs, prepare_go_inputs

LOG_PREFIX = "enrichment"
SCRIPT_PATH = Path(__file__).resolve().parent / "scripts" / "run_go_enrichment.R"


def _log(message: str) -> None:
    log(message, LOG_PREFIX)


def run_go_enrichment_analysis(
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

    selected_genes_path = output_dir / "selected_genes_ensembl.txt"
    universe_genes_path = output_dir / "universe_genes_ensembl.txt"
    universe_count, selected_count = prepare_go_inputs(
        de_results_path=de_results_path,
        selected_genes_path=selected_genes_path,
        universe_genes_path=universe_genes_path,
        padj_cutoff=config.padj_cutoff,
        lfc_cutoff=config.lfc_cutoff,
    )
    _log(f"GO universe: {universe_count} genes")
    _log(f"Selected DE genes: {selected_count}")

    runner = ExternalToolRunner(executable=executable, display_name="Rscript", log=_log)
    runner.run(
        [
            runner.path_arg(SCRIPT_PATH),
            "--selected-genes",
            runner.path_arg(selected_genes_path),
            "--universe-genes",
            runner.path_arg(universe_genes_path),
            "--outdir",
            runner.path_arg(output_dir),
        ],
        missing_message="Rscript not found in PATH",
    )

    _log("Done")
    return output_dir


def generate_go_enrichment_plots(config: EnrichmentConfig) -> Path:
    de_results_path = config.resolved_de_results_path()
    output_dir = config.resolved_enrichment_dir()
    selected_genes_path = output_dir / "selected_genes_ensembl.txt"
    universe_genes_path = output_dir / "universe_genes_ensembl.txt"

    _log("Generating GO result tables, summary, and plot with Python")
    generate_go_outputs(
        all_results_path=output_dir / "go_overrepresentation_all.csv",
        selected_genes_path=selected_genes_path,
        universe_genes_path=universe_genes_path,
        de_results_path=de_results_path,
        padj_cutoff=config.padj_cutoff,
        lfc_cutoff=config.lfc_cutoff,
    )

    _log("Done")
    return output_dir


def run_go_enrichment(
    config: EnrichmentConfig,
    executable: str = "Rscript",
) -> Path:
    run_go_enrichment_analysis(config, executable=executable)
    return generate_go_enrichment_plots(config)
