#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x)) y else x

parse_args <- function(args) {
  values <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key)
    }
    if (i == length(args) || startsWith(args[[i + 1]], "--")) {
      values[[substring(key, 3)]] <- TRUE
      i <- i + 1
    } else {
      values[[substring(key, 3)]] <- args[[i + 1]]
      i <- i + 2
    }
  }
  values
}

for (package in c("DESeq2", "tximport")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "Missing required R package '", package, "'. Install it with BiocManager before running this step.",
      call. = FALSE
    )
  }
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
samples_path <- args[["samples"]]
tx2gene_path <- args[["tx2gene"]]
outdir <- args[["outdir"]]
min_count <- as.integer(args[["min-count"]] %||% "10")
min_samples <- as.integer(args[["min-samples"]] %||% "2")

if (is.null(samples_path) || is.null(tx2gene_path) || is.null(outdir)) {
  stop("Usage: run_deseq2.R --samples sample_table.tsv --tx2gene tx2gene.tsv --outdir results_dir")
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

samples <- read.delim(samples_path, check.names = FALSE, stringsAsFactors = FALSE)
tx2gene <- read.delim(tx2gene_path, check.names = FALSE, stringsAsFactors = FALSE)

required_sample_cols <- c("sample_name", "patient", "condition")
missing_sample_cols <- setdiff(required_sample_cols, colnames(samples))
if (length(missing_sample_cols) > 0) {
  stop("Sample table is missing columns: ", paste(missing_sample_cols, collapse = ", "))
}

required_tx2gene_cols <- c("transcript_id", "gene_id", "gene_symbol")
missing_tx2gene_cols <- setdiff(required_tx2gene_cols, colnames(tx2gene))
if (length(missing_tx2gene_cols) > 0) {
  stop("Transcript-to-gene table is missing columns: ", paste(missing_tx2gene_cols, collapse = ", "))
}

sample_names <- samples$sample_name
quant_files <- setNames(samples$quant_sf, sample_names)
missing_quant_files <- quant_files[!file.exists(quant_files)]
if (length(missing_quant_files) > 0) {
  stop("Missing Salmon quant.sf files: ", paste(missing_quant_files, collapse = ", "))
}

tx2gene_import <- unique(tx2gene[, c("transcript_id", "gene_id")])
txi <- tximport::tximport(
  files = quant_files,
  type = "salmon",
  tx2gene = tx2gene_import,
  ignoreTxVersion = TRUE
)

gene_annot <- unique(tx2gene[, c("gene_id", "gene_symbol")])
gene_annot <- gene_annot[!duplicated(gene_annot$gene_id), , drop = FALSE]
rownames(gene_annot) <- gene_annot$gene_id

rownames(samples) <- samples$sample_name
samples <- samples[sample_names, , drop = FALSE]
samples$patient <- factor(samples$patient)
samples$condition <- relevel(factor(samples$condition), ref = "normal")

keep <- rowSums(txi$counts >= min_count) >= min_samples
txi$counts <- txi$counts[keep, , drop = FALSE]
txi$abundance <- txi$abundance[keep, , drop = FALSE]
txi$length <- txi$length[keep, , drop = FALSE]
gene_annot <- gene_annot[rownames(txi$counts), , drop = FALSE]
if (nrow(txi$counts) == 0) {
  stop("No genes remain after count filtering.")
}

dds <- DESeq2::DESeqDataSetFromTximport(
  txi = txi,
  colData = samples,
  design = ~ patient + condition
)
dds <- DESeq2::DESeq(dds)

res <- DESeq2::results(dds, contrast = c("condition", "tumor", "normal"))
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- merge(gene_annot, res_df, by = "gene_id", all.y = TRUE, sort = FALSE)
res_df <- res_df[, c("gene_id", "gene_symbol", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
res_df <- res_df[order(res_df$padj, res_df$pvalue, na.last = TRUE), ]
write.csv(res_df, file.path(outdir, "deseq2_all_genes.csv"), row.names = FALSE)

summary_path <- file.path(outdir, "deseq2_summary.txt")
sink(summary_path)
cat("DESeq2 paired analysis: tumor vs normal\n")
cat("Import: tximport gene-level Salmon estimates with average transcript-length offset\n")
cat("Design: ~ patient + condition\n\n")
cat("Samples:\n")
print(samples[, intersect(c("sample_name", "patient", "condition", "srr", "quant_sf"), colnames(samples)), drop = FALSE])
cat("\nGenes retained after filtering:", nrow(txi$counts), "\n")
cat("Significant genes padj < 0.05:", sum(!is.na(res_df$padj) & res_df$padj < 0.05), "\n\n")
print(summary(res))
sink()

normalized_counts <- as.data.frame(DESeq2::counts(dds, normalized = TRUE))
normalized_counts$gene_id <- rownames(normalized_counts)
normalized_counts <- merge(gene_annot, normalized_counts, by = "gene_id", all.y = TRUE, sort = FALSE)
write.csv(normalized_counts, file.path(outdir, "normalized_counts.csv"), row.names = FALSE)

vst <- DESeq2::vst(dds, blind = FALSE)
vst_counts <- as.data.frame(SummarizedExperiment::assay(vst))
vst_counts$gene_id <- rownames(vst_counts)
vst_counts <- merge(gene_annot, vst_counts, by = "gene_id", all.y = TRUE, sort = FALSE)
write.csv(vst_counts, file.path(outdir, "vst_counts.csv"), row.names = FALSE)
