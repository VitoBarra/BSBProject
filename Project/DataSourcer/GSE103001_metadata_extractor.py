from __future__ import annotations

import re
from pathlib import Path

from . import DataSourceConfig
from .metadata_pipeline import (
    SampleRecord,
    build_rows,
    ensure_soft_file,
    iter_soft_sample_blocks,
    sample_first_value,
    write_tsv,
)
from Log.log_util import log

LOG_PREFIX = "build_metadata_table"


def parse_samples(soft_path: Path) -> list[SampleRecord]:
    log(f"Parsing samples from SOFT: {soft_path}", LOG_PREFIX)
    samples: list[SampleRecord] = []
    for gsm, lines in iter_soft_sample_blocks(soft_path):
        title = sample_first_value(lines, "!Sample_title = ")
        source_name = sample_first_value(lines, "!Sample_source_name_ch1 = ")
        platform_id = sample_first_value(lines, "!Sample_platform_id = ")
        srx = ""
        for line in lines:
            if line.startswith("!Sample_relation = SRA: "):
                match = re.search(r"SRX\d+", line)
                if match:
                    srx = match.group(0)
        match = re.match(r"Pat_(\d+-\d+)_(normal|tumor)", title)
        if match:
            samples.append(
                SampleRecord(
                    sample_id=match.group(1),
                    title=title,
                    condition=match.group(2),
                    gsm=gsm,
                    source_name=source_name,
                    platform_id=platform_id,
                    srx=srx,
                )
            )
    log(f"Parsed {len(samples)} sample entries", LOG_PREFIX)
    return samples


def choose_patients(samples: list[SampleRecord], num_pairs: int, required_conditions: set[str]) -> list[str]:
    patients = sorted({sample.sample_id for sample in samples})
    selected: list[str] = []
    for patient in patients:
        conditions = {s.condition for s in samples if s.sample_id == patient}
        if required_conditions.issubset(conditions):
            selected.append(patient)
        if len(selected) == num_pairs:
            break
    log(f"Selected {len(selected)} matched pairs: {', '.join(selected)}", LOG_PREFIX)
    return selected


def build_metadata_table(config: DataSourceConfig) -> Path:
    if config.num_pairs <= 0:
        raise ValueError("num_pairs must be positive")

    profile = config.profile
    log(f"Starting metadata build for {profile.accession}", LOG_PREFIX)
    log(f"Requested matched pairs: {config.num_pairs}", LOG_PREFIX)
    soft_path = ensure_soft_file(config.resolved_soft_path(), profile.soft_url, LOG_PREFIX)
    samples = parse_samples(soft_path)
    selected_patients = choose_patients(samples, config.num_pairs, set(profile.required_conditions))
    if len(selected_patients) < config.num_pairs:
        raise RuntimeError(f"Only found {len(selected_patients)} complete matched pairs in {profile.accession}")

    rows = build_rows(samples, selected_patients, LOG_PREFIX)
    output_path = config.resolved_metadata_path()
    write_tsv(rows, output_path, LOG_PREFIX)

    log(f"Series: {profile.accession}", LOG_PREFIX)
    log(f"Selected samples: {', '.join(selected_patients)}", LOG_PREFIX)
    log(f"Rows written: {len(rows)}", LOG_PREFIX)
    log(f"Output: {output_path}", LOG_PREFIX)
    log(f"Total GB: {sum(float(row['sample_total_gb']) for row in rows):.3f}", LOG_PREFIX)
    log("Done", LOG_PREFIX)
    return output_path
