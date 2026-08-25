#!/usr/bin/env Rscript

for (package in c("DESeq2", "tximport", "optparse")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "Missing required R package '", package, "'. Install it with BiocManager before running this step.",
      call. = FALSE
    )
  }
}

parser <- optparse::OptionParser(option_list = list(
  optparse::make_option("--samples", type = "character", help = "DESeq2 sample sheet TSV"),
  optparse::make_option("--transcript-to-gene-map", dest = "transcript_to_gene_map", type = "character", help = "Transcript-to-gene map TSV"),
  optparse::make_option("--results-path", dest = "results_path", type = "character", help = "DESeq2 results CSV"),
  optparse::make_option("--vst-counts-path", dest = "vst_counts_path", type = "character", help = "VST counts CSV"),
  optparse::make_option("--summary-path", dest = "summary_path", type = "character", help = "DESeq2 summary TXT"),
  optparse::make_option("--min-count", dest = "min_count", type = "integer", default = 10L),
  optparse::make_option("--min-samples", dest = "min_samples", type = "integer", default = 2L),
  optparse::make_option("--plots-only", dest = "plots_only", action = "store_true", default = FALSE,
                         help = "Regenerate figures from existing DESeq2 outputs")
))
args <- optparse::parse_args(parser)
samples_path <- args$samples
transcript_to_gene_map_path <- args$transcript_to_gene_map
results_path <- args$results_path
vst_counts_path <- args$vst_counts_path
summary_path <- args$summary_path
min_count <- args$min_count
min_samples <- args$min_samples

if (is.null(samples_path) || is.null(transcript_to_gene_map_path) || is.null(results_path) || is.null(vst_counts_path) || is.null(summary_path)) {
  optparse::print_help(parser)
  stop("--samples, --transcript-to-gene-map, --results-path, --vst-counts-path, and --summary-path are required", call. = FALSE)
}

for (output_path in c(results_path, vst_counts_path, summary_path)) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
}

plot_de_outputs <- function(results_path, vst_counts_path, samples_path) {
  output_dir <- dirname(results_path)
  results <- read.csv(results_path, check.names = FALSE, stringsAsFactors = FALSE)
  vst_counts <- read.csv(vst_counts_path, check.names = FALSE, stringsAsFactors = FALSE)
  samples <- read.delim(samples_path, check.names = FALSE, stringsAsFactors = FALSE)
  required_results <- c("gene_id", "gene_symbol", "baseMean", "log2FoldChange", "padj")
  if (!all(required_results %in% colnames(results))) stop("DESeq2 result table is missing required plotting columns.")
  if (!all(c("gene_id", "gene_symbol", samples$sample_name) %in% colnames(vst_counts))) stop("VST matrix is missing required plotting columns.")
  write.csv(results[!is.na(results$padj) & results$padj < 0.05, , drop = FALSE],
            file.path(output_dir, "deseq2_significant_genes_padj_0.05.csv"), row.names = FALSE)

  label <- function(rows) ifelse(is.na(rows$gene_symbol) | rows$gene_symbol == "", rows$gene_id, rows$gene_symbol)
  direction <- function(rows) {
    ifelse(!is.na(rows$padj) & !is.na(rows$log2FoldChange) & rows$padj < 0.05 & rows$log2FoldChange >= 1, "Up",
      ifelse(!is.na(rows$padj) & !is.na(rows$log2FoldChange) & rows$padj < 0.05 & rows$log2FoldChange <= -1, "Down", "Not significant"))
  }
  colors <- c("Down" = "#5B8DD9", "Not significant" = "#BDBDBD", "Up" = "#B7352D")
  groups <- direction(results)
  finite <- is.finite(results$log2FoldChange) & is.finite(results$padj)
  png(file.path(output_dir, "volcano_padj.png"), width = 2160, height = 1860, res = 300)
  plot(results$log2FoldChange[finite], -log10(pmax(results$padj[finite], .Machine$double.xmin)), pch = 16, cex = .45,
       col = colors[groups[finite]], xlab = "log2 fold change", ylab = "-log10 adjusted p-value",
       main = "Differential expression: tumor vs adjacent normal")
  abline(v = c(-1, 1), h = -log10(.05), lty = 2, lwd = .7)
  legend("topright", legend = names(colors), col = colors, pch = 16, bty = "n")
  for (group in c("Down", "Up")) {
    selected <- which(finite & groups == group)
    selected <- selected[order(results$padj[selected])][seq_len(min(10, length(selected)))]
    if (length(selected)) text(results$log2FoldChange[selected], -log10(pmax(results$padj[selected], .Machine$double.xmin)), label(results)[selected], cex = .55, pos = 3)
  }
  dev.off()

  finite <- is.finite(results$baseMean) & results$baseMean > 0 & is.finite(results$log2FoldChange)
  png(file.path(output_dir, "ma_plot.png"), width = 2550, height = 1650, res = 300)
  plot(log2(results$baseMean[finite]), results$log2FoldChange[finite], pch = 16, cex = .4, col = colors[groups[finite]],
       xlab = "log2 mean normalized expression", ylab = "log2 fold change", main = "MA plot: tumor vs adjacent normal")
  abline(h = c(-1, 0, 1), lty = c(2, 1, 2), lwd = c(.6, .7, .6))
  legend("topright", legend = names(colors), col = colors, pch = 16, bty = "n")
  dev.off()

  sample_names <- samples$sample_name
  matrix <- as.matrix(vst_counts[, sample_names, drop = FALSE])
  storage.mode(matrix) <- "double"
  pca <- prcomp(t(matrix), scale. = FALSE)
  condition_colors <- c(normal = "#222222", tumor = "#F5A000")
  png(file.path(output_dir, "pca_vst.png"), width = 1800, height = 1500, res = 300)
  plot(pca$x[, 1], pca$x[, 2], pch = 16, cex = 1.3, col = condition_colors[samples$condition],
       xlab = sprintf("PC1 (%.1f%%)", 100 * summary(pca)$importance[2, 1]),
       ylab = sprintf("PC2 (%.1f%%)", 100 * summary(pca)$importance[2, 2]), main = "PCA of VST expression")
  text(pca$x[, 1], pca$x[, 2], labels = sample_names, pos = 3, cex = .65)
  legend("topright", legend = names(condition_colors), col = condition_colors, pch = 16, bty = "n")
  dev.off()

  png(file.path(output_dir, "sample_distance_heatmap.png"), width = 1950, height = 1740, res = 300)
  heatmap(as.matrix(dist(t(matrix))), symm = TRUE, col = hcl.colors(100, "YlOrRd", rev = TRUE), main = "Sample distances (VST)")
  dev.off()

  ranked <- results[is.finite(results$padj), , drop = FALSE]
  ranked <- ranked[order(ranked$padj, -abs(ranked$log2FoldChange)), , drop = FALSE]
  ranked <- ranked[seq_len(min(20, nrow(ranked))), , drop = FALSE]
  indices <- match(ranked$gene_id, vst_counts$gene_id)
  indices <- indices[!is.na(indices)]
  if (length(indices) >= 2) {
    top_matrix <- matrix[indices, , drop = FALSE]
    rownames(top_matrix) <- make.unique(label(ranked)[!is.na(match(ranked$gene_id, vst_counts$gene_id))])
    row_z <- t(scale(t(top_matrix)))
    row_z[!is.finite(row_z)] <- 0
    png(file.path(output_dir, "top_de_gene_heatmap.png"), width = 2460, height = 2640, res = 300)
    heatmap(row_z, col = hcl.colors(101, "Blue-Red 3"), scale = "none", margins = c(8, 11),
            main = "Top 20 differentially expressed genes (row z-scores)")
    dev.off()
  }
}

if (isTRUE(args$plots_only)) {
  plot_de_outputs(results_path, vst_counts_path, samples_path)
  quit(save = "no")
}

samples <- read.delim(samples_path, check.names = FALSE, stringsAsFactors = FALSE)
transcript_to_gene_map <- read.delim(transcript_to_gene_map_path, check.names = FALSE, stringsAsFactors = FALSE)

required_sample_cols <- c("sample_name", "patient", "condition")
missing_sample_cols <- setdiff(required_sample_cols, colnames(samples))
if (length(missing_sample_cols) > 0) {
  stop("Sample table is missing columns: ", paste(missing_sample_cols, collapse = ", "))
}

required_transcript_to_gene_map_columns <- c("transcript_id", "gene_id", "gene_symbol")
missing_transcript_to_gene_map_columns <- setdiff(required_transcript_to_gene_map_columns, colnames(transcript_to_gene_map))
if (length(missing_transcript_to_gene_map_columns) > 0) {
  stop("Transcript-to-gene map is missing columns: ", paste(missing_transcript_to_gene_map_columns, collapse = ", "))
}

sample_names <- samples$sample_name
if (anyNA(samples[, required_sample_cols]) || any(!nzchar(trimws(unlist(samples[, required_sample_cols]))))) {
  stop("Sample names, patients, and conditions must all be non-empty.")
}
if (anyDuplicated(sample_names)) {
  stop("Sample names must be unique: ", paste(unique(sample_names[duplicated(sample_names)]), collapse = ", "))
}
expected_conditions <- c("normal", "tumor")
observed_conditions <- sort(unique(samples$condition))
if (!identical(observed_conditions, expected_conditions)) {
  stop(
    "Expected exactly the conditions 'normal' and 'tumor'; found: ",
    paste(observed_conditions, collapse = ", ")
  )
}
patient_condition_counts <- table(samples$patient, samples$condition)
if (any(patient_condition_counts != 1L)) {
  invalid_patients <- rownames(patient_condition_counts)[apply(patient_condition_counts != 1L, 1, any)]
  stop(
    "Paired design requires exactly one normal and one tumor sample per patient. Invalid patients: ",
    paste(invalid_patients, collapse = ", ")
  )
}
quant_files <- setNames(samples$quant_sf, sample_names)
missing_quant_files <- quant_files[!file.exists(quant_files)]
if (length(missing_quant_files) > 0) {
  stop("Missing Salmon quant.sf files: ", paste(missing_quant_files, collapse = ", "))
}

transcript_to_gene_map_import <- unique(transcript_to_gene_map[, c("transcript_id", "gene_id")])
txi <- tximport::tximport(
  files = quant_files,
  type = "salmon",
  tx2gene = transcript_to_gene_map_import,
  ignoreTxVersion = TRUE
)

gene_annot <- unique(transcript_to_gene_map[, c("gene_id", "gene_symbol")])
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
write.csv(res_df, results_path, row.names = FALSE)

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

# Library-size-normalized VST values for visualization, not differential-expression testing.
vst <- DESeq2::vst(dds, blind = FALSE)
vst_counts <- as.data.frame(SummarizedExperiment::assay(vst))
vst_counts$gene_id <- rownames(vst_counts)
vst_counts <- merge(gene_annot, vst_counts, by = "gene_id", all.y = TRUE, sort = FALSE)
write.csv(vst_counts, vst_counts_path, row.names = FALSE)

plot_de_outputs(results_path, vst_counts_path, samples_path)
