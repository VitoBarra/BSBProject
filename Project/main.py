from __future__ import annotations

import argparse
from pathlib import Path

from DataSourcer import DataSourceConfig, PROFILES, build_metadata_table
from DataSourcer.download_fastq import download_fastq_from_tsv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build dataset metadata tables and optionally download FASTQ files.")
    parser.add_argument(
        "--dataset",
        default="gse103001",
        choices=sorted(PROFILES),
        help="Dataset key to process.",
    )
    parser.add_argument("--num-pairs", type=int, default=4, help="Number of matched patient pairs to include.")
    parser.add_argument("--soft-path", type=Path, default=None, help="Optional custom path for the GEO SOFT file.")
    parser.add_argument("--metadata-path", type=Path, default=None, help="Optional custom output/input TSV path.")
    parser.add_argument("--fastq-dest", type=Path, default=None, help="Optional custom destination directory for FASTQ files.")
    parser.add_argument("--workers", type=int, default=6, help="Number of FASTQ downloads to run in parallel.")
    parser.add_argument("--skip-build", action="store_true", help="Skip metadata table generation and use the existing TSV.")
    parser.add_argument("--skip-download", action="store_true", help="Skip FASTQ download after metadata table generation.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    config = DataSourceConfig(
        dataset_key=args.dataset,
        num_pairs=args.num_pairs,
        soft_path=args.soft_path,
        metadata_path=args.metadata_path,
        fastq_dest=args.fastq_dest,
        download_workers=args.workers,
    )

    if not args.skip_build:
        build_metadata_table(config)

    if not args.skip_download:
        download_fastq_from_tsv(config)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
