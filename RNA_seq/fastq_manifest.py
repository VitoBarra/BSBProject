from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

FASTQ_SUFFIX = ".fastq.gz"


@dataclass(frozen=True, slots=True)
class PairedFastq:
    patient_id: str
    srr: str
    read_1: Path
    read_2: Path


def load_paired_fastqs(metadata_path: Path, fastq_dir: Path) -> list[PairedFastq]:
    if not metadata_path.exists():
        raise FileNotFoundError(f"Missing sample metadata: {metadata_path}")
    if not fastq_dir.exists():
        raise FileNotFoundError(f"Missing FASTQ directory: {fastq_dir}")

    with metadata_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"patient_id", "srr", "fastq_filename"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{metadata_path} must contain columns: {', '.join(sorted(required))}")

        pairs: list[PairedFastq] = []

        for row in reader:
            patient_id = row["patient_id"].strip()
            srr = row["srr"].strip()

            filenames = [value.strip() for value in row["fastq_filename"].split(";") if value.strip()]
            read_1, read_2 = (fastq_dir / filename for filename in filenames)

            # Check file existence.
            for fastq_path in (read_1, read_2):
                if not fastq_path.is_file():
                    raise FileNotFoundError(f"Missing expected FASTQ file for {srr}: {fastq_path}")

            pairs.append(PairedFastq(patient_id=patient_id, srr=srr, read_1=read_1, read_2=read_2))

    if not pairs:
        raise ValueError(f"No paired FASTQ samples found in {metadata_path}")

    return pairs
