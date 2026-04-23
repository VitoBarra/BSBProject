from __future__ import annotations

import argparse
from pathlib import Path

from DataSourcer import DataSourceConfig, PROFILES, build_metadata_table, download_reference_transcriptome
from Quantification import QuantificationConfig, run_salmon
from QualityControl import QualityControlConfig, run_fastp, run_fastqc, run_multiqc
from DataSourcer.download_fastq import download_fastq_from_tsv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build dataset metadata tables and optionally download FASTQ files.")
    parser.add_argument(
        "--dataset",
        default="gse103001",
        choices=sorted(PROFILES),
        help="Dataset key to process.",
    )
    parser.add_argument(
        "--build-metadata-table",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Build the dataset metadata TSV before downstream steps.",
    )
    parser.add_argument(
        "--download-fastq",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Download FASTQ files listed in the metadata TSV.",
    )
    parser.add_argument(
        "--download-reference-transcriptome",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Download the reference transcriptome FASTA used to build a Salmon index.",
    )
    parser.add_argument("--workers", type=int, default=6, help="Number of FASTQ downloads to run in parallel.")
    parser.add_argument("--num-pairs", type=int, default=4, help="Number of matched patient pairs to include.")
    parser.add_argument("--soft-path", type=Path, default=None, help="Optional custom path for the GEO SOFT file.")
    parser.add_argument("--metadata-path", type=Path, default=None, help="Optional custom output/input TSV path.")

    parser.add_argument(
        "--run-raw-fastqc",
        default=False,
        action="store_true",
        help="Run FastQC on raw reads.",
    )
    parser.add_argument(
        "--run-trimmed-fastqc",
        default=False,
        action="store_true",
        help="Run FastQC on trimmed reads.",
    )

    parser.add_argument("--fastqc-threads", type=int, default=6, help="Worker threads to pass to FastQC.")
    parser.add_argument("--fastq-dest", type=Path, default=None, help="Optional custom destination directory for FASTQ files.")
    parser.add_argument("--reference-dir", type=Path, default=None, help="Optional destination directory for reference assets.")
    parser.add_argument(
        "--transcriptome-fasta-path",
        type=Path,
        default=None,
        help="Optional path for the reference transcriptome FASTA(.gz).",
    )
    parser.add_argument(
        "--transcriptome-url",
        default=None,
        help="Optional transcriptome FASTA URL. Overrides the dataset default.",
    )
    parser.add_argument("--fastqc-report-out", type=Path, default=None, help="Optional FastQC output directory for raw reads.")
    parser.add_argument("--fastqc-trimmed-report-out", type=Path, default=None, help="Optional FastQC output directory for trimmed reads.")

    parser.add_argument("--run-fastp", default= False, action="store_true", help="Run fastp on the FASTQ files after raw FastQC.")
    parser.add_argument("--trimmed-fastq-dest", type=Path, default=None, help="Optional custom destination directory for trimmed FASTQ files.")
    parser.add_argument("--fastp-report-out", type=Path, default=None, help="Optional fastp report directory.")
    parser.add_argument("--fastp-threads", type=int, default=6, help="Worker threads to pass to fastp.")

    parser.add_argument("--run-salmon", default=False, action="store_true", help="Run Salmon quantification on trimmed paired-end reads.")
    parser.add_argument("--salmon-threads", type=int, default=6, help="Worker threads to pass to Salmon.")
    parser.add_argument(
        "--build-salmon-index",
        default=False,
        action="store_true",
        help="Build the Salmon index automatically when it is missing.",
    )
    parser.add_argument("--salmon-index-dir", type=Path, default=None, help="Optional path to a prebuilt Salmon index directory.")
    parser.add_argument(
        "--salmon-transcriptome-fasta",
        type=Path,
        default=None,
        help="Optional transcriptome FASTA(.gz) used to build a Salmon index when requested.",
    )
    parser.add_argument("--salmon-quant-dir", type=Path, default=None, help="Optional output directory for Salmon quantification results.")
    parser.add_argument("--salmon-libtype", default="A", help="Salmon library type, for example A or ISR.")

    parser.add_argument("--run-multiqc", default=False, action="store_true", help="Run MultiQC on all QC results, including FastQC and fastp reports.")
    parser.add_argument("--multiqc-report-out", type=Path, default=None, help="Optional MultiQC output directory.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    config = DataSourceConfig(
        dataset_key=args.dataset,
        num_pairs=args.num_pairs,
        soft_path=args.soft_path,
        metadata_path=args.metadata_path,
        fastq_dest=args.fastq_dest,
        reference_dir=args.reference_dir,
        transcriptome_fasta_path=args.transcriptome_fasta_path,
        transcriptome_url=args.transcriptome_url,
        download_workers=args.workers,
    )

    if args.build_metadata_table:
        build_metadata_table(config)

    if args.download_fastq:
        download_fastq_from_tsv(config)

    if args.download_reference_transcriptome:
        download_reference_transcriptome(config)

    qc_config = QualityControlConfig(
        dataset_key=args.dataset,
        fastq_dir=args.fastq_dest,
        trimmed_fastq_dir=args.trimmed_fastq_dest,
        fastp_report_out=args.fastp_report_out,
        fastqc_report_out=args.fastqc_report_out,
        fastqc_trimmed_report_out=args.fastqc_trimmed_report_out,
        multiqc_report_out=args.multiqc_report_out,
    )

    quant_config = QuantificationConfig(
        dataset_key=args.dataset,
        trimmed_fastq_dir=args.trimmed_fastq_dest,
        salmon_index_dir=args.salmon_index_dir,
        salmon_quant_dir=args.salmon_quant_dir,
        salmon_transcriptome_fasta=args.salmon_transcriptome_fasta or args.transcriptome_fasta_path,
    )

    if args.run_raw_fastqc:
        run_fastqc(
            qc_config,
            threads=args.fastqc_threads,
            use_trimmed_reads=False
        )

    if args.run_fastp:
        run_fastp(
            qc_config,
            threads=args.fastp_threads,
        )

    if args.run_trimmed_fastqc:
        run_fastqc(
            qc_config,
            threads=args.fastqc_threads,
            use_trimmed_reads=True,
        )

    if args.run_multiqc:
        run_multiqc(qc_config)

    if args.run_salmon:
        run_salmon(
            quant_config,
            threads=args.salmon_threads,
            libtype=args.salmon_libtype,
            build_index_if_missing=args.build_salmon_index,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
