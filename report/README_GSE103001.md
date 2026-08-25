# GSE103001 workflow

Run the reproducible RNA-seq workflow from the repository root. Analysis tools
run inside the project Docker image.

```bash
make docker-build
make analysis
```

The workflow selects four complete tumor/normal patient pairs, downloads their
paired FASTQ files, verifies their ENA MD5 checksums, performs raw and
cleaned-read QC, quantifies the Ensembl
release 115 GRCh38 transcriptome with Salmon, imports `quant.sf` files through
`tximport`, runs paired DESeq2 analysis, and performs GO enrichment.

The Salmon index is checked for its required files and tied to the reference
transcriptome by SHA-256 checksum. Existing quantifications are reused only
when both `quant.sf` and `aux_info/meta_info.json` pass structural and numeric
validation; incomplete outputs stop the workflow with an explicit error.
Completed fastp outputs are likewise reused only when both cleaned FASTQ files
and the HTML and JSON reports are present and the reported read counts are
consistent.

Useful individual targets are:

```bash
make download-data
make qc-raw
make trim
make qc-trimmed
make quantify
make dea
make enrichment
```

The selected samples and resolved download URLs are recorded in
`data/RNA-seq/GSE103001/GSE103001_selected_4pairs.tsv`. Patient 12-02 is excluded
because its ENA record mixes singleton and paired FASTQ files.

See `report/Part1_Report_Barra.md` for methods, results, limitations, and the
principal output files.
