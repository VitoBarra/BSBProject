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

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop(
    "Missing required R package 'DESeq2'. Install it with BiocManager before running this step.",
    call. = FALSE
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
counts_path <- args[["counts"]]
samples_path <- args[["samples"]]
outdir <- args[["outdir"]]
min_count <- as.integer(args[["min-count"]] %||% "10")
min_samples <- as.integer(args[["min-samples"]] %||% "2")

if (is.null(counts_path) || is.null(samples_path) || is.null(outdir)) {
  stop("Usage: run_deseq2.R --counts counts.tsv --samples sample_table.tsv --outdir results_dir")
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

counts_raw <- read.delim(counts_path, check.names = FALSE, stringsAsFactors = FALSE)
samples <- read.delim(samples_path, check.names = FALSE, stringsAsFactors = FALSE)

required_sample_cols <- c("sample_name", "patient", "condition")
missing_sample_cols <- setdiff(required_sample_cols, colnames(samples))
if (length(missing_sample_cols) > 0) {
  stop("Sample table is missing columns: ", paste(missing_sample_cols, collapse = ", "))
}

required_count_cols <- c("gene_id", "gene_symbol")
missing_count_cols <- setdiff(required_count_cols, colnames(counts_raw))
if (length(missing_count_cols) > 0) {
  stop("Count table is missing columns: ", paste(missing_count_cols, collapse = ", "))
}

sample_names <- samples$sample_name
missing_counts <- setdiff(sample_names, colnames(counts_raw))
if (length(missing_counts) > 0) {
  stop("Count matrix is missing sample columns: ", paste(missing_counts, collapse = ", "))
}

gene_annot <- counts_raw[, c("gene_id", "gene_symbol")]
counts <- as.matrix(counts_raw[, sample_names, drop = FALSE])
storage.mode(counts) <- "numeric"
counts <- round(counts)
rownames(counts) <- gene_annot$gene_id

rownames(samples) <- samples$sample_name
samples <- samples[sample_names, , drop = FALSE]
samples$patient <- factor(samples$patient)
samples$condition <- relevel(factor(samples$condition), ref = "normal")

keep <- rowSums(counts >= min_count) >= min_samples
counts <- counts[keep, , drop = FALSE]
gene_annot <- gene_annot[keep, , drop = FALSE]
if (nrow(counts) == 0) {
  stop("No genes remain after count filtering.")
}

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts,
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
cat("Design: ~ patient + condition\n\n")
cat("Samples:\n")
print(samples[, intersect(c("sample_name", "patient", "condition", "srr", "quant_sf"), colnames(samples)), drop = FALSE])
cat("\nGenes retained after filtering:", nrow(counts), "\n")
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
