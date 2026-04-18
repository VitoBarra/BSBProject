from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class DatasetInfo:
    accession: str
    soft_url: str
    soft_filename: str
    required_conditions: frozenset[str]
    excluded_patients: frozenset[str]


@dataclass(frozen=True, slots=True)
class ReferenceInfo:
    name: str
    transcriptome_url: str
    transcriptome_filename: str
    genome_url: str
    genome_filename: str


GSE103001 = DatasetInfo(
    accession="GSE103001",
    soft_url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103001/soft/GSE103001_family.soft.gz",
    soft_filename="GSE103001_family.soft.gz",
    required_conditions=frozenset({"normal", "tumor"}),
    # 12-02 tumor has an anomalous mixed ENA layout (standalone plus _1/_2).
    # Use the next complete patient so all eight selected runs are paired-only.
    excluded_patients=frozenset({"12-02"}),
)

ENSEMBL_GRCH38_CDNA = ReferenceInfo(
    name="Ensembl GRCh38 cDNA release 115",
    transcriptome_url=(
        "https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cdna/"
        "Homo_sapiens.GRCh38.cdna.all.fa.gz"
    ),
    transcriptome_filename="transcriptome.fa.gz",
    genome_url=(
        "https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/"
        "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    ),
    genome_filename="genome.primary_assembly.fa.gz",
)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_ROOT = PROJECT_ROOT / "data"
RNA_SEQ_DATA_ROOT = DATA_ROOT / "RNA-seq"
DATASET_ROOT = RNA_SEQ_DATA_ROOT / GSE103001.accession
REFERENCE_ROOT = RNA_SEQ_DATA_ROOT / "reference"


@dataclass(frozen=True, slots=True)
class RNASeqPaths:
    dataset_root: Path = DATASET_ROOT
    reference_root: Path = REFERENCE_ROOT

    @property
    def raw_fastq(self) -> Path:
        return self.dataset_root / "raw_fastq"

    @property
    def trimmed_fastq(self) -> Path:
        return self.dataset_root / "trimmed_fastq"

    @property
    def fastp_reports(self) -> Path:
        return self.dataset_root / "qc" / "fastp"

    @property
    def raw_fastqc_reports(self) -> Path:
        return self.dataset_root / "qc" / "fastqc"

    @property
    def trimmed_fastqc_reports(self) -> Path:
        return self.dataset_root / "qc" / "fastqc_trimmed"

    @property
    def raw_multiqc_reports(self) -> Path:
        return self.dataset_root / "qc" / "multiqc_raw"

    @property
    def trimmed_multiqc_reports(self) -> Path:
        return self.dataset_root / "qc" / "multiqc_trimmed"

    @property
    def salmon_index(self) -> Path:
        return self.dataset_root / "quant" / "salmon_index"

    @property
    def salmon_quant(self) -> Path:
        return self.dataset_root / "quant" / "salmon"

    @property
    def de_root(self) -> Path:
        return self.dataset_root / "de"

    @property
    def de_results(self) -> Path:
        return self.de_root / "results"

    @property
    def enrichment(self) -> Path:
        return self.dataset_root / "enrichment"


DEFAULT_PATHS = RNASeqPaths()


@dataclass(frozen=True, slots=True)
class RNASeqConfig:
    num_pairs: int = 4
    datasetInfo: DatasetInfo = GSE103001
    paths: RNASeqPaths = DEFAULT_PATHS

    def __post_init__(self) -> None:
        if self.num_pairs <= 0:
            raise ValueError("num_pairs must be positive")

    def resolved_metadata_path(self) -> Path:
        return self.paths.dataset_root / f"{self.datasetInfo.accession}_selected_{self.num_pairs}pairs.tsv"


DEFAULT_CONFIG = RNASeqConfig()
