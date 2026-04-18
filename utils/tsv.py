"""Helpers for reading and writing UTF-8 tab-separated tables."""

from __future__ import annotations

import csv
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path


def read_tsv(
    path: Path, *, encoding: str = "utf-8", required_columns: set[str] | None = None
) -> list[dict[str, str]]:
    with path.open(encoding=encoding, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if required_columns and not required_columns.issubset(reader.fieldnames or ()):
            raise ValueError(f"{path} must contain columns: {', '.join(sorted(required_columns))}")
        return list(reader)


def write_tsv(
    path: Path,
    rows: Iterable[Mapping[str, object]],
    fieldnames: Sequence[str],
    *,
    lineterminator: str = "\n",
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator=lineterminator)
        writer.writeheader()
        writer.writerows(rows)
