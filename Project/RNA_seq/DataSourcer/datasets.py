from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True, frozen=True)
class DatasetProfile:
    key: str
    accession: str
    soft_url: str
    required_conditions: tuple[str, ...] = ()
    excluded_sample_ids: tuple[str, ...] = ()
    metadata_filename_template: str = "{accession}_selected_{num_pairs}samples.tsv"
    soft_filename_template: str = "{accession}_family.soft.gz"
    transcriptome_url: str | None = None
    transcriptome_filename: str = "transcriptome.fa.gz"

    def soft_filename(self) -> str:
        return self.soft_filename_template.format(accession=self.accession)

    def metadata_filename(self, num_pairs: int) -> str:
        return self.metadata_filename_template.format(accession=self.accession, num_pairs=num_pairs)


GSE103001_PROFILE = DatasetProfile(
    key="gse103001",
    accession="GSE103001",
    soft_url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103001/soft/GSE103001_family.soft.gz",
    metadata_filename_template="{accession}_selected_{num_pairs}pairs.tsv",
    required_conditions=("normal", "tumor"),
    # 12-02 tumor has an anomalous mixed ENA layout (standalone plus _1/_2).
    # Use the next complete patient so all eight selected runs are paired-only.
    excluded_sample_ids=("12-02",),
    transcriptome_url="https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",
)


PROFILES = {
    GSE103001_PROFILE.key: GSE103001_PROFILE,
}
