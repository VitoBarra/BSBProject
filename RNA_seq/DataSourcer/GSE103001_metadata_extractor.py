from __future__ import annotations

import csv
import gzip
import logging
import re
from pathlib import Path

from Bio import Geo

from .. import GSE103001
from . import DataSourceConfig
from .metadata_pipeline import (
    NORMALIZED_FIELDS,
    SampleRecord,
    build_rows,
    ensure_soft_file,
    write_tsv,
)
LOGGER = logging.getLogger("build_metadata_table")



def parse_samples(soft_path: Path) -> list[SampleRecord]:
    LOGGER.info("Parsing samples from SOFT: %s", soft_path)
    samples: list[SampleRecord] = []

    with gzip.open(soft_path, "rt", encoding="utf-8", errors="ignore") as handle:
        for record in Geo.parse(handle):
            if record.entity_type != "SAMPLE":
                continue

            attributes = record.entity_attributes
            title = str(attributes.get("Sample_title", ""))
            relations = attributes.get("Sample_relation", [])
            if isinstance(relations, str):
                relations = [relations]
            srx_match = next(
                (match for relation in relations if (match := re.search(r"SRX\d+", relation))),
                None,
            )
            srx = srx_match.group(0) if srx_match else ""

            match = re.match(r"Pat_(\d+-\d+)_(normal|tumor)", title)
            if match:
                samples.append(
                    SampleRecord(
                        patient_id=match.group(1),  # Patient ID, e.g., "1-1".
                        condition=match.group(2),  # Either "normal" or "tumor".
                        gsm=record.entity_id,
                        srx=srx,
                    )
                )
    LOGGER.info("Parsed %d sample entries", len(samples))
    return samples


def choose_samples(
        samples: list[SampleRecord],
        num_pairs: int,
        required_conditions: set[str],
        excluded_sample_ids: set[str] | None = None,
        ) -> list[SampleRecord]:

    excluded = excluded_sample_ids or set()
    patients = sorted({sample.patient_id for sample in samples})
    selected: list[str] = []

    for patient in patients:

        # Skip excluded patients known to have incomplete data.
        if patient in excluded:
            continue

        # Exclude patients without every required condition (normal and tumor).
        conditions = {s.condition for s in samples if s.patient_id == patient}
        if required_conditions.issubset(conditions):
            selected.append(patient)

        # Stop early once enough pairs have been selected.
        if len(selected) == num_pairs:
            break

    LOGGER.info("Selected %d matched pairs: %s", len(selected), ", ".join(selected))
    if len(selected) < num_pairs:
        raise RuntimeError(f"Only found {len(selected)} complete matched pairs")

    return [sample for sample in samples if sample.patient_id in selected]


def build_metadata_table(config: DataSourceConfig) -> None:
    LOGGER.info("Starting metadata build for %s", GSE103001.accession)
    LOGGER.info("Requested matched pairs: %d", config.num_pairs)
    soft_path = ensure_soft_file(config.resolved_soft_path(), GSE103001.soft_url)
    samples = parse_samples(soft_path)
    selected_samples = choose_samples(
        samples,
        config.num_pairs,
        set(GSE103001.required_conditions),
        set(GSE103001.excluded_patients),
    )
    rows = build_rows(selected_samples)
    output_path = config.resolved_metadata_path()
    write_tsv(rows, output_path)
    LOGGER.info("Series: %s", GSE103001.accession)
    LOGGER.info("Selected samples: %s", ", ".join(sorted({sample.patient_id for sample in selected_samples})))
    LOGGER.info("Rows written: %d", len(rows))
    LOGGER.info("Output: %s", output_path)
    LOGGER.info("Done")


def ensure_metadata_table(config: DataSourceConfig) -> None:
    """Return the current sample metadata table, rebuilding it when needed."""
    metadata_path = config.resolved_metadata_path()
    if metadata_path.exists():
        with metadata_path.open(encoding="utf-8-sig", newline="") as handle:
            current_fields = next(csv.reader(handle, delimiter="\t"), [])
        if current_fields == NORMALIZED_FIELDS:
            LOGGER.info("Metadata already present: %s", metadata_path)
            return
        LOGGER.info("Metadata schema is outdated, rebuilding: %s", metadata_path)

    build_metadata_table(config)
