from __future__ import annotations

import logging
from pathlib import Path

from .. import ENSEMBL_GRCH38_CDNA, RNASeqConfig
from .http_download import download_resumable

LOGGER = logging.getLogger("download_reference")
USER_AGENT = "BSBProject-GSE103001/1.0"


def _download(url: str, destination: Path) -> None:
    already_present, _downloaded = download_resumable(
        url,
        destination,
        user_agent=USER_AGENT,
    )
    if already_present:
        LOGGER.info("Reference transcriptome already present: %s", destination)
        return
    LOGGER.info("Saved reference transcriptome to %s", destination)


def download_reference_transcriptome(config: RNASeqConfig) -> None:
    """Download the transcriptome and genome needed for a decoy-aware index."""
    LOGGER.info("Reference: %s", ENSEMBL_GRCH38_CDNA.name)
    for label, url, filename in (
        ("Transcriptome", ENSEMBL_GRCH38_CDNA.transcriptome_url, ENSEMBL_GRCH38_CDNA.transcriptome_filename),
        ("Primary-assembly genome", ENSEMBL_GRCH38_CDNA.genome_url, ENSEMBL_GRCH38_CDNA.genome_filename),
    ):
        destination = config.paths.reference_root / filename
        LOGGER.info("%s URL: %s", label, url)
        LOGGER.info("Destination file: %s", destination)
        _download(url, destination)
    LOGGER.info("Done")
