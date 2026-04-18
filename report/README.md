# Project reports

This directory contains the final reports and supporting documentation for the
two project parts. Final deliverables are separated from workflow and source
notes so the report material is easy to navigate.

## Part 1 — RNA-seq analysis

- Final report: [Paired RNA-seq analysis of ER-positive breast cancer](part1/report.md)
- The report includes the workflow, data sources and reproducibility instructions.

## Part 2 — Network medicine

- Final report: [Network medicine analysis](part2/report.md)
- The report describes the network workflow, inputs, results and limitations.

## Shared assets

`assets/` contains the MultiQC overview images used by the Part 1 report.
Analysis data, figures, and machine-readable results remain in `../data/` so
the reports link directly to their reproducible outputs.

## Delivery checklist

The assignment explicitly requests a short report, analysis scripts or notebook,
DEG table, and enrichment results and discussion for Part 1. Part 2 requires
the network analysis described in its report. For a repository-based delivery,
include:

- Both reports: `report/part1/report.md` and `report/part2/report.md`.
- The analysis code, `Makefile`, `Dockerfile`, and `requirements.txt`.
- The complete DEG table:
  `data/RNA-seq/GSE103001/de/results/deseq2_all_genes.csv`.
- The GO enrichment tables under `data/RNA-seq/GSE103001/enrichment/`.
- The Part 2 inputs, rankings, selected module, separation results, figures,
  and networks under `data/network_medicine/processed/` and `results/`.
- The report figures and other linked tracked files, keeping their relative
  paths intact. A clean Git clone includes these compact deliverables.

Raw FASTQs, intermediate QC reports, the Salmon index, raw API responses, and
local credentials are not tracked. Do not include credentials in a submission.
The assignment does not specify a submission format or location; confirm these
with the course before submitting. If PDFs are required, export both reports
and check their figures and links, while still providing the requested tables
and code separately.
