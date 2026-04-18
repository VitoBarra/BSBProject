from __future__ import annotations

from pathlib import Path

from . import DataSourceConfig
from .metadata_pipeline import (
    SampleRecord,
    ensure_soft_file,
    iter_soft_sample_blocks,
    sample_first_value,
    write_tsv,
)
from Log.log_util import log

LOG_PREFIX = "build_metadata_table"


def infer_condition(title: str, source_name: str) -> str:
    normalized = f"{title} {source_name}".strip().lower()
    if "normal" in normalized:
        return "normal"
    if "tumor" in normalized or "adenocarcinoma" in normalized:
        return "tumor"
    return ""


def parse_rows(soft_path: Path) -> list[SampleRecord]:
    log(f"Parsing samples from SOFT: {soft_path}", LOG_PREFIX)
    rows: list[SampleRecord] = []

    for gsm, lines in iter_soft_sample_blocks(soft_path):
        title = sample_first_value(lines, "!Sample_title = ")
        source_name = sample_first_value(lines, "!Sample_source_name_ch1 = ")
        platform_id = sample_first_value(lines, "!Sample_platform_id = ")
        rows.append(
            SampleRecord(
                sample_id=title,
                title=title,
                condition=infer_condition(title, source_name),
                gsm=gsm,
                source_name=source_name,
                platform_id=platform_id,
            )
        )

    log(f"Parsed {len(rows)} samples", LOG_PREFIX)
    return rows


def build_metadata_table(config: DataSourceConfig) -> Path:
    if config.num_pairs <= 0:
        raise ValueError("num_pairs must be positive")

    profile = config.profile
    log(f"Starting metadata build for {profile.accession}", LOG_PREFIX)
    log(f"Requested samples: {config.num_pairs}", LOG_PREFIX)

    soft_path = ensure_soft_file(config.resolved_soft_path(), profile.soft_url, LOG_PREFIX)
    samples = parse_rows(soft_path)
    selected_rows = samples[: config.num_pairs]
    if len(selected_rows) < config.num_pairs:
        raise RuntimeError(f"Only found {len(selected_rows)} samples in {profile.accession}")

    output_path = config.resolved_metadata_path()
    write_tsv([
        {
            "sample_id": sample.sample_id,
            "title": sample.title,
            "condition": sample.condition,
            "gsm": sample.gsm,
            "source_name": sample.source_name,
            "platform_id": sample.platform_id,
            "srx": "",
            "srr": "",
            "fastq_url": "",
            "sample_total_bytes": "0",
            "sample_total_gb": "0.000",
        }
        for sample in selected_rows
    ], output_path, LOG_PREFIX)

    log(f"Series: {profile.accession}", LOG_PREFIX)
    log(f"Rows written: {len(selected_rows)}", LOG_PREFIX)
    log(f"Output: {output_path}", LOG_PREFIX)
    log("Done", LOG_PREFIX)
    return output_path
