from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path

from .. import DEFAULT_CONFIG, RNASeqConfig

_GIB = 1024**3
_DENSE_INDEX_MIN_MEMORY_BYTES = 48 * _GIB


def available_memory_bytes() -> int | None:
    """Return memory available to this process, respecting Linux cgroup limits."""
    cgroup_limit = Path("/sys/fs/cgroup/memory.max")
    cgroup_current = Path("/sys/fs/cgroup/memory.current")
    try:
        limit = cgroup_limit.read_text().strip()
        if limit != "max":
            return max(0, int(limit) - int(cgroup_current.read_text().strip()))
    except (OSError, ValueError):
        pass

    meminfo = Path("/proc/meminfo")
    try:
        values = dict(
            line.replace(":", "").split(maxsplit=1)
            for line in meminfo.read_text().splitlines()
        )
        return int(values["MemAvailable"].split()[0]) * 1024
    except (OSError, KeyError, ValueError):
        pass

    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_AVPHYS_PAGES")
    except (AttributeError, ValueError, OSError):
        return None


@dataclass(slots=True)
class QuantificationConfig:
    workflow: RNASeqConfig = DEFAULT_CONFIG

    def resolved_salmon_sparse(self) -> bool:
        """Use sparse indexing unless enough memory is available for dense."""
        memory = available_memory_bytes()
        return memory is None or memory < _DENSE_INDEX_MIN_MEMORY_BYTES

    def resolved_salmon_threads(self) -> int:
        """Choose workers from the memory available to this process."""
        memory = available_memory_bytes()
        if memory is None:
            return 1
        memory_limited_workers = max(1, memory // (8 * _GIB))
        return min(os.cpu_count() or 1, memory_limited_workers)

    def resolved_metadata_path(self) -> Path:
        return self.workflow.resolved_metadata_path()

    def reference_dir(self) -> Path:
        return self.workflow.paths.reference_root

    def resolved_trimmed_fastq_dir(self) -> Path:
        return self.workflow.paths.trimmed_fastq

    def resolved_salmon_index_dir(self) -> Path:
        return self.workflow.paths.salmon_index

    def resolved_salmon_quant_dir(self) -> Path:
        return self.workflow.paths.salmon_quant

    def resolved_salmon_transcriptome_fasta(self) -> Path | None:
        for filename in (
            "transcriptome.fa.gz",
            "transcriptome.fasta.gz",
            "transcriptome.fa",
            "transcriptome.fasta",
            "transcripts.fa.gz",
            "transcripts.fasta.gz",
            "transcripts.fa",
            "transcripts.fasta",
        ):
            candidate = self.reference_dir() / filename
            if candidate.exists():
                return candidate
        return None

    def resolved_salmon_genome_fasta(self) -> Path | None:
        candidate = self.reference_dir() / "genome.primary_assembly.fa.gz"
        return candidate if candidate.exists() else None

    def resolved_salmon_gentrome_fasta(self) -> Path:
        return self.reference_dir() / "gentrome.fa.gz"

    def resolved_salmon_decoys_path(self) -> Path:
        return self.reference_dir() / "decoys.txt"

from .salmon_runner import build_salmon_index, run_salmon
