from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True, frozen=True)
class DatasetProfile:
    key: str
    accession: str
    soft_url: str
    required_conditions: tuple[str, ...] = ()
    metadata_filename_template: str = "{accession}_selected_{num_pairs}samples.tsv"
    soft_filename_template: str = "{accession}_family.soft.gz"

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
)

GSE210787_PROFILE = DatasetProfile(
    key="gse210787",
    accession="GSE210787",
    soft_url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE210nnn/GSE210787/soft/GSE210787_family.soft.gz",
)

GSE31210_PROFILE = DatasetProfile(
    key="gse31210",
    accession="GSE31210",
    soft_url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE31nnn/GSE31210/soft/GSE31210_family.soft.gz",
)

PROFILES = {
    GSE103001_PROFILE.key: GSE103001_PROFILE,
    GSE210787_PROFILE.key: GSE210787_PROFILE,
    GSE31210_PROFILE.key: GSE31210_PROFILE,
}
