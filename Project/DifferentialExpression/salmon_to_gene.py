from __future__ import annotations

import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path
from typing import TextIO

from Log.log_util import log

from . import SalmonGeneAggregationConfig

LOG_PREFIX = "salmon-to-gene"
FASTA_GENE_RE = re.compile(r"(?:^|\s)gene:([^\s]+)")
FASTA_SYMBOL_RE = re.compile(r"(?:^|\s)gene_symbol:([^\s]+)")


def _log(message: str) -> None:
    log(message, LOG_PREFIX)


def _strip_version(identifier: str) -> str:
    return identifier.split(".", 1)[0]


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _build_tx2gene(transcriptome_path: Path, output_path: Path) -> dict[str, tuple[str, str]]:
    mapping: dict[str, tuple[str, str]] = {}
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with _open_text(transcriptome_path) as fasta, output_path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.writer(output, delimiter="\t")
        writer.writerow(["transcript_id", "gene_id", "gene_symbol"])

        for line in fasta:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            transcript_id = _strip_version(header.split(maxsplit=1)[0])
            gene_match = FASTA_GENE_RE.search(header)
            if gene_match is None:
                raise ValueError(f"Transcript FASTA header has no gene identifier: {header}")
            gene_id = _strip_version(gene_match.group(1))
            symbol_match = FASTA_SYMBOL_RE.search(header)
            gene_symbol = symbol_match.group(1) if symbol_match else ""

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


def _read_tx2gene(path: Path) -> dict[str, tuple[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"transcript_id", "gene_id", "gene_symbol"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain columns: {', '.join(sorted(required))}")
        mapping: dict[str, tuple[str, str]] = {}
        for row in reader:
            symbol = row["gene_symbol"]
            if symbol.upper() == "NA":
                symbol = ""
            mapping[_strip_version(row["transcript_id"])] = (_strip_version(row["gene_id"]), symbol)
        return mapping


def _read_samples(metadata_path: Path, quant_dir: Path) -> list[dict[str, str | Path]]:
    with metadata_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_id", "condition", "srr"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{metadata_path} must contain columns: {', '.join(sorted(required))}")

        samples: list[dict[str, str | Path]] = []
        seen_names: set[str] = set()
        for row in reader:
            sample_name = f"{row['sample_id']}_{row['condition']}"
            if sample_name in seen_names:
                raise ValueError(f"Duplicate sample name in metadata: {sample_name}")
            seen_names.add(sample_name)
            samples.append(
                {
                    "sample_name": sample_name,
                    "patient": row["sample_id"],
                    "condition": row["condition"],
                    "srr": row["srr"],
                    "quant_sf": quant_dir / row["srr"] / "quant.sf",
                }
            )
    if not samples:
        raise ValueError(f"No samples found in {metadata_path}")
    return samples


def _aggregate_quant(
    quant_path: Path,
    tx2gene: dict[str, tuple[str, str]],
) -> tuple[dict[str, float], dict[str, float], dict[str, str]]:
    if not quant_path.exists():
        raise FileNotFoundError(f"Missing Salmon quantification file: {quant_path}")

    counts: defaultdict[str, float] = defaultdict(float)
    tpms: defaultdict[str, float] = defaultdict(float)
    symbols: dict[str, str] = {}
    unmapped = 0

    with quant_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Name", "TPM", "NumReads"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{quant_path} must contain columns: {', '.join(sorted(required))}")

        for row in reader:
            transcript_id = _strip_version(row["Name"])
            gene = tx2gene.get(transcript_id)
            if gene is None:
                unmapped += 1
                continue
            gene_id, gene_symbol = gene
            counts[gene_id] += float(row["NumReads"])
            tpms[gene_id] += float(row["TPM"])
            symbols.setdefault(gene_id, gene_symbol)

    if unmapped:
        _log(f"Warning: skipped {unmapped} transcripts without a gene mapping in {quant_path}")
    return dict(counts), dict(tpms), symbols


def _write_matrix(
    path: Path,
    sample_names: list[str],
    gene_ids: list[str],
    symbols: dict[str, str],
    values: list[dict[str, float]],
    decimals: int,
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene_id", "gene_symbol", *sample_names])
        for gene_id in gene_ids:
            writer.writerow(
                [
                    gene_id,
                    symbols.get(gene_id, ""),
                    *(f"{sample.get(gene_id, 0.0):.{decimals}f}" for sample in values),
                ]
            )


def build_salmon_gene_matrices(config: SalmonGeneAggregationConfig) -> Path:
    metadata_path = config.resolved_metadata_path()
    quant_dir = config.resolved_salmon_quant_dir()
    transcriptome_path = config.resolved_transcriptome_fasta_path()
    output_dir = config.resolved_output_dir()
    tx2gene_path = config.resolved_tx2gene_path()

    for path, label in ((metadata_path, "metadata"), (quant_dir, "Salmon quantification directory")):
        if not path.exists():
            raise FileNotFoundError(f"Missing {label}: {path}")
    output_dir.mkdir(parents=True, exist_ok=True)

    if tx2gene_path.exists():
        _log(f"Using transcript-to-gene map: {tx2gene_path}")
        tx2gene = _read_tx2gene(tx2gene_path)
    else:
        if not transcriptome_path.exists():
            raise FileNotFoundError(f"Missing transcriptome FASTA: {transcriptome_path}")
        _log(f"Building transcript-to-gene map from: {transcriptome_path}")
        tx2gene = _build_tx2gene(transcriptome_path, tx2gene_path)

    samples = _read_samples(metadata_path, quant_dir)
    sample_names = [str(sample["sample_name"]) for sample in samples]
    count_values: list[dict[str, float]] = []
    tpm_values: list[dict[str, float]] = []
    symbols: dict[str, str] = {}

    for index, sample in enumerate(samples, start=1):
        _log(f"[{index}/{len(samples)}] Aggregating {sample['sample_name']}")
        counts, tpms, sample_symbols = _aggregate_quant(Path(sample["quant_sf"]), tx2gene)
        count_values.append(counts)
        tpm_values.append(tpms)
        symbols.update({gene_id: symbol for gene_id, symbol in sample_symbols.items() if symbol})

    gene_ids = sorted(set().union(*(values.keys() for values in count_values)))
    _write_matrix(output_dir / "salmon_gene_counts.tsv", sample_names, gene_ids, symbols, count_values, decimals=6)
    _write_matrix(output_dir / "salmon_gene_tpm.tsv", sample_names, gene_ids, symbols, tpm_values, decimals=6)
    _write_matrix(output_dir / "salmon_gene_counts_rounded.tsv", sample_names, gene_ids, symbols, count_values, decimals=0)

    with (output_dir / "sample_table.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_name", "patient", "condition", "srr", "quant_sf"])
        for sample in samples:
            writer.writerow(
                [
                    sample["sample_name"],
                    sample["patient"],
                    sample["condition"],
                    sample["srr"],
                    Path(sample["quant_sf"]).resolve(),
                ]
            )

    _log(f"Wrote DEA inputs to: {output_dir}")
    return output_dir
