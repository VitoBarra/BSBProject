# Bioinformatics and Systems Biology — Assignment

This repository reproduces a gene-level RNA-seq comparison of ER-positive
breast tumours and matched adjacent tissue, then applies a network-medicine
analysis to breast cancer, rheumatoid arthritis, and diabetes mellitus.

The assignment brief is available in
[docs/Assignment_Barra.md](docs/Assignment_Barra.md). The finished material is
organised here:

- [Part 1 report: paired RNA-seq analysis](report/part1/report.md)
- [Part 2 report: network medicine analysis](report/part2/report.md)
- [Report index and supporting notes](report/README.md)

## Requirements

The reproducible environment is defined in `Dockerfile`. The fully pinned
`requirements.txt` is the authoritative list of Python dependencies and is
installed by the Docker build. Install Docker and run commands from the
repository root. Build the image once:

```bash
make docker-build
```

Use `make help` to see all available targets.

## Part 1 — RNA-seq workflow

The workflow uses four matched ER-positive tumour/adjacent-normal pairs from
GSE103001. It downloads the reads, runs QC and fastp trimming, quantifies with
Salmon against a decoy-aware GRCh38 reference, performs paired DESeq2 analysis,
and runs GO enrichment.

```bash
make analysis
```

Individual stages can be run when required:

```bash
make download-data
make qc-raw
make trim
make qc-trimmed
make quantify
make dea
make enrichment
```

The full download, trimmed reads, and Salmon index require substantial storage
(approximately 74 GB in the completed local analysis), network access, and
time. Raw and intermediate data are intentionally excluded from Git; compact
result tables and figures required to inspect the reports are tracked.

## Part 2 — Network medicine workflow

The first workflow creates a STRING network from the 72 breast-cancer genes in
Supplemental Table S11, calculates centralities and communities, and selects a
10-gene breast-cancer module:

```bash
make network-bc
```

The disease-separation workflow maps the breast-cancer module, DisGeNET
rheumatoid-arthritis genes, and diabetes-mellitus control genes to a common
physical STRING background:

```bash
make network-separation
```

To refresh the curated DisGeNET associations, create the ignored local
credential file and set its API key:

```bash
make init-env
# Add DISGENET_API_KEY=... to .env.local
make network-disgenet-download
```

## Reproducibility and outputs

The Makefile runs analysis tools inside Docker and is the supported execution
interface. The reports link directly to the machine-readable output tables,
figures, and GraphML networks under `data/`. Analysis choices, limitations,
software versions, database versions, and access dates are documented in the
relevant report and workflow note.
