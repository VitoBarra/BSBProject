from __future__ import annotations

import argparse
import logging

from RNA_seq import RNASeqConfig
from RNA_seq.DataSourcer import (
    DataSourceConfig, download_fastq, download_reference_transcriptome,
)
from RNA_seq.DifferentialExpression import (
    DifferentialExpressionConfig,
    generate_deseq2_plots,
    prepare_deseq2_inputs,
    run_deseq2_analysis,
)
from RNA_seq.Enrichment import EnrichmentConfig, generate_go_enrichment_plots, run_go_enrichment_analysis
from RNA_seq.QualityControl import QualityControlConfig, run_fastp, run_fastqc, run_multiqc
from RNA_seq.Quantification import QuantificationConfig, build_salmon_index, run_salmon


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the GSE103001 RNA-seq workflow.")
    parser.add_argument(
        "--download-RnaSeq-data",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Prepare metadata, paired FASTQ files, and the reference transcriptome; reuse completed files.",
    )

    parser.add_argument("--workers", type=int, default=6, help="Number of FASTQ downloads to run in parallel.")
    parser.add_argument("--num-pairs", type=int, default=4, help="Number of matched patient pairs to include.")

    parser.add_argument(
        "--run-raw-fastqc",
        action="store_true",
        help="Run FastQC on raw reads.",
    )
    parser.add_argument(
        "--run-trimmed-fastqc",
        action="store_true",
        help="Run FastQC on trimmed reads.",
    )

    parser.add_argument("--fastqc-threads", type=int, default=6, help="Worker threads to pass to FastQC.")

    parser.add_argument("--run-fastp", action="store_true", help="Run fastp on the FASTQ files after raw FastQC.")
    parser.add_argument("--fastp-threads", type=int, default=6, help="Worker threads to pass to fastp.")
    parser.add_argument("--run-raw-multiqc", action="store_true", help="Aggregate raw-read FastQC results with MultiQC.")
    parser.add_argument(
        "--run-trimmed-multiqc",
        action="store_true",
        help="Aggregate trimmed-read FastQC and fastp results with MultiQC.",
    )

    parser.add_argument("--run-salmon", action="store_true", help="Run Salmon quantification on trimmed paired-end reads.")
    parser.add_argument(
        "--build-salmon-index",
        action="store_true",
        help="Build or validate the decoy-aware Salmon index without quantifying reads.",
    )

    parser.add_argument(
        "--run-dea",
        action="store_true",
        help="Prepare inputs, run paired DESeq2 when results are missing, and regenerate DEA plots.",
    )
    parser.add_argument("--de-min-count", type=int, default=10, help="Minimum count used for DESeq2 prefiltering.")
    parser.add_argument("--de-min-samples", type=int, default=2, help="Minimum number of samples passing --de-min-count.")

    parser.add_argument(
        "--run-enrichment",
        action="store_true",
        help="Run GO enrichment when results are missing and regenerate its tables, summary, and plot.",
    )
    parser.add_argument("--enrichment-padj-cutoff", type=float, default=0.05, help="Adjusted p-value cutoff for enrichment input genes.")
    parser.add_argument("--enrichment-lfc-cutoff", type=float, default=0.0, help="Absolute log2 fold-change cutoff for enrichment input genes.")
    return parser.parse_args()


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="[%(name)s] %(message)s")
    args = parse_args()
    workflow = RNASeqConfig(num_pairs=args.num_pairs)
    data_config = DataSourceConfig(
        workflow=workflow,
        download_workers=args.workers,
    )
    qc_config = QualityControlConfig(
        workflow=workflow,
        fastqc_threads=args.fastqc_threads,
        fastp_threads=args.fastp_threads,
    )
    quant_config = QuantificationConfig(workflow=workflow)

    de_config = DifferentialExpressionConfig(
        workflow=workflow,
        min_count=args.de_min_count,
        min_samples=args.de_min_samples,
    )

    enrichment_config = EnrichmentConfig(
        workflow=workflow,
        padj_cutoff=args.enrichment_padj_cutoff,
        lfc_cutoff=args.enrichment_lfc_cutoff,
    )


    if args.download_RnaSeq_data:
        download_fastq(data_config)
        download_reference_transcriptome(data_config.workflow)

    if args.run_raw_fastqc:
        run_fastqc(
            qc_config,
            fastq_dir=qc_config.resolved_fastq_dir(),
            out_dir=qc_config.resolved_fastqc_report_out(),
        )

    if args.run_fastp:
        run_fastp(qc_config)

    if args.run_trimmed_fastqc:
        run_fastqc(
            qc_config,
            fastq_dir=qc_config.resolved_trimmed_fastq_dir(),
            out_dir=qc_config.resolved_fastqc_trimmed_report_out(),
        )

    if args.run_raw_multiqc:
        run_multiqc(
            input_dirs=[qc_config.resolved_fastqc_report_out()],
            output_dir=qc_config.resolved_multiqc_raw_report_out(),
            report_name="multiqc_raw.html",
        )

    if args.run_trimmed_multiqc:
        run_multiqc(
            input_dirs=[qc_config.resolved_fastqc_trimmed_report_out(), qc_config.resolved_fastp_report_out()],
            output_dir=qc_config.resolved_multiqc_trimmed_report_out(),
            report_name="multiqc_trimmed.html",
        )

    if args.run_salmon:
        build_salmon_index(quant_config)
        run_salmon(quant_config)

    if args.run_dea:
        prepare_deseq2_inputs(de_config)
        run_deseq2_analysis(de_config)
        generate_deseq2_plots(de_config)

    if args.run_enrichment:
        run_go_enrichment_analysis(enrichment_config)
        generate_go_enrichment_plots(enrichment_config)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
