from __future__ import annotations

import csv
import gzip
import logging
import re
from pathlib import Path
from typing import TextIO

from ..path_utils import portable_path
from . import DifferentialExpressionConfig

LOGGER = logging.getLogger("prepare-deseq2-inputs")
FASTA_GENE_RE = re.compile(r"(?:^|\s)gene:([^\s]+)")
FASTA_SYMBOL_RE = re.compile(r"(?:^|\s)gene_symbol:([^\s]+)")


def _strip_version(identifier: str) -> str:
    return identifier.split(".", 1)[0]


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def build_transcript_to_gene_map(transcriptome_path: Path, output_path: Path) -> dict[str, tuple[str, str]]:
    mapping: dict[str, tuple[str, str]] = {}
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with _open_text(transcriptome_path) as fasta, output_path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.writer(output, delimiter="\t")
        writer.writerow(["transcript_id", "gene_id", "gene_symbol"])

        for line in fasta:
            # skip any non-header lines in the FASTA file
            if not line.startswith(">"):
                continue

            # Extract transcript ID, gene ID, and gene symbol from the FASTA header

            # Remove the leading '>' and any leading/trailing whitespace from the header line
            header = line[1:].strip()

            # extract the transcript ID from the header (the first word after the '>')
            transcript_id = _strip_version(header.split(maxsplit=1)[0])

            # extract the gene ID from the header using a regex search
            gene_match = FASTA_GENE_RE.search(header)
            if gene_match is None: # If no gene identifier is found, raise an error
                raise ValueError(f"Transcript FASTA header has no gene identifier: {header}")
            gene_id = _strip_version(gene_match.group(1))

            # Extract the gene symbol from the header using a regex search (if present)
            symbol_match = FASTA_SYMBOL_RE.search(header)
            gene_symbol = symbol_match.group(1) if symbol_match else ""

            # Check for conflicting mappings: if the transcript ID has already been mapped to a different gene ID or symbol, raise an error
            previous = mapping.get(transcript_id)
            current = (gene_id, gene_symbol)
            if previous is not None and previous != current:
                raise ValueError(f"Conflicting gene mappings for transcript {transcript_id}: {previous} and {current}")

            if previous is None:
                mapping[transcript_id] = current
                writer.writerow([transcript_id, gene_id, gene_symbol])

    if not mapping:
        raise ValueError(f"No transcript-to-gene mappings found in {transcriptome_path}")

    return mapping


def build_deseq2_sample_records(metadata_path: Path, quant_dir: Path) -> list[dict[str, str | Path]]:
    with metadata_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"patient_id", "condition", "srr"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{metadata_path} must contain columns: {', '.join(sorted(required))}")

        samples: list[dict[str, str | Path]] = []
        seen_names: set[str] = set()
        for row in reader:
            sample_name = f"{row['patient_id']}_{row['condition']}"
            if sample_name in seen_names:
                raise ValueError(f"Duplicate sample name in metadata: {sample_name}")
            seen_names.add(sample_name)
            samples.append(
                {
                    "sample_name": sample_name,
                    "patient": row["patient_id"],
                    "condition": row["condition"],
                    "srr": row["srr"],
                    "quant_sf": quant_dir / row["srr"] / "quant.sf",
                }
            )
    if not samples:
        raise ValueError(f"No samples found in {metadata_path}")
    return samples


def write_deseq2_sample_sheet(path: Path, sample_records: list[dict[str, str | Path]]) -> None:
    fieldnames = ["sample_name", "patient", "condition", "srr", "quant_sf"]

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(
            {**sample, "quant_sf": portable_path(sample["quant_sf"])}
            for sample in sample_records
        )


def prepare_transcript_to_gene_map(config: DifferentialExpressionConfig) -> None:
    """Create the transcript-to-gene map required by tximport, if needed."""
    transcript_to_gene_map_path = config.resolved_transcript_to_gene_map_path()

    if transcript_to_gene_map_path.exists():
        LOGGER.info("Using transcript-to-gene map: %s", transcript_to_gene_map_path)
        return

    transcriptome_path = config.resolved_transcriptome_fasta_path()
    LOGGER.info("Building transcript-to-gene map from: %s", transcriptome_path)
    build_transcript_to_gene_map(transcriptome_path, transcript_to_gene_map_path)


def prepare_deseq2_sample_sheet(config: DifferentialExpressionConfig) -> None:
    """Create the DESeq2 sample sheet and verify its Salmon quantification files."""
    metadata_path = config.resolved_metadata_path()
    quant_dir = config.resolved_salmon_quant_dir()
    deseq2_sample_sheet_path = config.resolved_deseq2_sample_sheet_path()

    # Verify that the metadata file and Salmon quantification directory exist
    for path, label in ((metadata_path, "metadata"), (quant_dir, "Salmon quantification directory")):
        if not path.exists():
            raise FileNotFoundError(f"Missing {label}: {path}")

    # Create the parent directory for the DESeq2 sample sheet if it doesn't exist
    deseq2_sample_sheet_path.parent.mkdir(parents=True, exist_ok=True)

    # Build the DESeq2 sample records from the metadata and check for missing Salmon quantification files
    sample_records = build_deseq2_sample_records(metadata_path, quant_dir)
    missing_quant_files = [
        Path(sample["quant_sf"])
        for sample in sample_records
        if not Path(sample["quant_sf"]).exists()
    ]
    if missing_quant_files:
        missing = ", ".join(str(path) for path in missing_quant_files)
        raise FileNotFoundError(f"Missing Salmon quantification files: {missing}")

    # Write the DESeq2 sample sheet to the specified path
    write_deseq2_sample_sheet(deseq2_sample_sheet_path, sample_records)


def prepare_deseq2_inputs(config: DifferentialExpressionConfig) -> None:
    """Prepare the transcript-to-gene map and DESeq2 sample sheet."""
    output_dir = config.resolved_output_dir()
    output_dir.mkdir(parents=True, exist_ok=True)

    prepare_transcript_to_gene_map(config)
    prepare_deseq2_sample_sheet(config)

    LOGGER.info("Prepared DESeq2 input files in: %s", output_dir)
