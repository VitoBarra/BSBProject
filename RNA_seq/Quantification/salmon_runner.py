from __future__ import annotations

import logging
import gzip
import json
import shutil
from dataclasses import dataclass
from pathlib import Path

from ExternalTools import ExternalToolRunner
from ..fastq_manifest import load_paired_fastqs

from . import QuantificationConfig

LOGGER = logging.getLogger("salmon")
@dataclass(slots=True, frozen=True)
class SalmonJob:
    name: str
    read_1: Path
    read_2: Path
    output_dir: Path


def collect_salmon_inputs(metadata_path: Path, trimmed_fastq_dir: Path, quant_dir: Path) -> list[SalmonJob]:
    return [
        SalmonJob(name=pair.srr, read_1=pair.read_1, read_2=pair.read_2, output_dir=quant_dir / pair.srr)
        for pair in load_paired_fastqs(metadata_path, trimmed_fastq_dir)
    ]

def build_salmon_index(config: QuantificationConfig) -> None:
    """
    Build a decoy-aware Salmon index from the transcriptome and GRCh38 genome.

    Reuses an existing decoy-aware index without rebuilding it.

    Args:
        config: Quantification settings used to resolve reference paths.
    """
    index_dir = config.resolved_salmon_index_dir()
    if index_dir.exists():
        try:
            num_decoys = json.loads((index_dir / "info.json").read_text(encoding="utf-8"))["num_decoys"]
        except (FileNotFoundError, KeyError, json.JSONDecodeError) as exc:
            raise RuntimeError(f"Existing Salmon index is not verifiably decoy-aware: {index_dir}") from exc
        if num_decoys > 0:
            return
        raise RuntimeError(
            "The existing Salmon index is transcriptome-only."
        )

    transcriptome_fasta = config.resolved_salmon_transcriptome_fasta()
    genome_fasta = config.resolved_salmon_genome_fasta()
    if transcriptome_fasta is None or genome_fasta is None:
        raise FileNotFoundError(
            "Missing reference FASTA required for a decoy-aware Salmon index. "
            f"Run `make download`; expected transcriptome and genome FASTAs under {config.reference_dir()}."
        )

    index_dir.parent.mkdir(parents=True, exist_ok=True)
    gentrome_fasta, decoys_path = _prepare_decoy_aware_reference(config, transcriptome_fasta, genome_fasta)
    LOGGER.info("Building decoy-aware Salmon index from gentrome FASTA: %s", gentrome_fasta)
    LOGGER.info("Salmon index output directory: %s", index_dir)
    sparse = config.resolved_salmon_sparse()
    threads = config.resolved_salmon_threads()
    LOGGER.info("Salmon index mode: %s; worker threads: %d", "sparse" if sparse else "dense", threads)

    runner = ExternalToolRunner(executable="salmon", display_name="Salmon", logger=LOGGER)
    runner.run(
        [
            # Build a Salmon index
            "index",
            # Index cDNA plus genomic decoy sequences
            "--transcripts",
            runner.path_arg(gentrome_fasta),
            # Mark all genomic contigs as decoys
            "--decoys",
            runner.path_arg(decoys_path),
            # Write the generated index to the configured directory
            "--index",
            runner.path_arg(index_dir),
            # Sparse pufferfish uses substantially less RAM for whole-genome
            # decoy-aware indexes, at the cost of somewhat slower mapping.
            *(["--sparse"] if sparse else []),
            # Run this many Salmon worker threads
            "--threads",
            str(threads),
        ],
        missing_message="salmon not found in PATH",
    )


def _prepare_decoy_aware_reference(
    config: QuantificationConfig, transcriptome_fasta: Path, genome_fasta: Path
) -> tuple[Path, Path]:
    """Create the concatenated reference and one-decoy-contig-per-line list."""
    gentrome_fasta = config.resolved_salmon_gentrome_fasta()
    decoys_path = config.resolved_salmon_decoys_path()

    if not decoys_path.exists():
        decoys_path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(genome_fasta, "rt", encoding="utf-8") as genome, decoys_path.open("w", encoding="utf-8") as decoys:
            for line in genome:
                if line.startswith(">"):
                    decoys.write(line[1:].split(maxsplit=1)[0] + "\n")

    if not gentrome_fasta.exists():
        temporary_path = gentrome_fasta.with_suffix(gentrome_fasta.suffix + ".tmp")
        # This is an intermediate consumed locally by Salmon. Fast compression
        # keeps its FASTA content identical while avoiding a lengthy CPU-bound
        # preprocessing step for the primary human assembly.
        with gzip.open(temporary_path, "wb", compresslevel=1) as destination:
            for source_path in (transcriptome_fasta, genome_fasta):
                with gzip.open(source_path, "rb") as source:
                    shutil.copyfileobj(source, destination, length=1024 * 1024)
        temporary_path.replace(gentrome_fasta)
    return gentrome_fasta, decoys_path


def run_salmon(
    config: QuantificationConfig,
) -> None:
    """Quantify trimmed paired-end reads with Salmon.

    Existing sample output directories are reused. The Salmon index must
    already exist; use ``build_salmon_index`` to create it when needed.

    Args:
        config: Quantification settings containing the input, index, output,
            and thread configuration.

    """
    trimmed_fastq_dir = config.resolved_trimmed_fastq_dir()
    quant_dir = config.resolved_salmon_quant_dir()

    index_dir = config.resolved_salmon_index_dir()
    if not index_dir.exists():
        raise FileNotFoundError(
            f"Missing Salmon index directory: {index_dir}. "
            "Run build_salmon_index before quantification."
        )
    jobs = collect_salmon_inputs(config.resolved_metadata_path(), trimmed_fastq_dir, quant_dir)
    quant_dir.mkdir(parents=True, exist_ok=True)

    LOGGER.info("Trimmed FASTQ input directory: %s", trimmed_fastq_dir)
    LOGGER.info("Salmon index directory: %s", index_dir)
    LOGGER.info("Salmon output directory: %s", quant_dir)
    LOGGER.info("Input datasets selected: %d", len(jobs))
    threads = config.resolved_salmon_threads()
    LOGGER.info("Salmon worker threads: %d", threads)

    runner = ExternalToolRunner(executable="salmon", display_name="Salmon", logger=LOGGER)
    for index, job in enumerate(jobs, start=1):
        if job.output_dir.exists():
            LOGGER.info("[%d/%d] Skipping existing quantification for %s", index, len(jobs), job.name)
            continue
        LOGGER.info("[%d/%d] Quantifying %s", index, len(jobs), job.name)
        runner.run(
            [
                # Run Salmon quantification
                "quant",
                # Use the configured transcriptome index
                "--index",
                runner.path_arg(index_dir),
                # Set the expected library orientation (Allows Salmon to autodetect)
                "--libType",
                "A",
                # Run this many Salmon worker threads.
                "--threads",
                str(threads),
                # Use selective alignment to validate read mappings.
                "--validateMappings",
                # Provide the paired-end trimmed reads
                "--mates1",
                runner.path_arg(job.read_1),
                "--mates2",
                runner.path_arg(job.read_2),
                # Write this sample's quantification results here
                "--output",
                runner.path_arg(job.output_dir),
            ],
            missing_message="salmon not found in PATH",
        )
    LOGGER.info("Done")
