# Analysis Pipeline

The DEA and enrichment steps are run inside Docker. The Python code does not need to know that it is inside Docker; the container provides `python3`, `Rscript`, DESeq2, and the enrichment packages.

Build the R/Bioconductor image from the repository root:

```bash
make docker-build
```

Run paired differential expression analysis:

```bash
make prepare-dea
make dea
```

`prepare-dea` reads every Salmon `quant.sf`, derives the transcript-to-gene
mapping from the Ensembl transcriptome FASTA when needed, and writes:

- `Project/data/GSE103001/de/tx2gene.tsv`
- `Project/data/GSE103001/de/salmon_gene_counts.tsv`
- `Project/data/GSE103001/de/salmon_gene_counts_rounded.tsv`
- `Project/data/GSE103001/de/salmon_gene_tpm.tsv`
- `Project/data/GSE103001/de/sample_table.tsv`

To revise or regenerate only the volcano plot, MA plot, and clustered DEG
heatmap from existing DESeq2 outputs:

```bash
make plots
```

Run GO enrichment from the DESeq2 table:

```bash
make enrichment
```

Run both steps:

```bash
make analysis
```

The Makefile mounts `Project/` into the container as `/project`, matching the paths already used in the sample table.

Main DEA outputs:

- `Project/data/GSE103001/de/results/deseq2_all_genes.csv`
- `Project/data/GSE103001/de/results/deseq2_significant_genes_padj_0.05.csv`
- `Project/data/GSE103001/de/results/deseq2_ranked_genes.csv`
- `Project/data/GSE103001/de/results/volcano_padj.pdf`
- `Project/data/GSE103001/de/results/ma_plot.pdf`
- `Project/data/GSE103001/de/results/sample_distance_heatmap.pdf`
- `Project/data/GSE103001/de/results/top_de_gene_heatmap.pdf`
- `Project/data/GSE103001/de/results/pca_vst.pdf`

Main enrichment outputs:

- `Project/data/GSE103001/enrichment/go_overrepresentation_all.csv`
- `Project/data/GSE103001/enrichment/go_overrepresentation_significant.csv`
- `Project/data/GSE103001/enrichment/go_overrepresentation_dotplot.pdf`
