# Plotting helpers for the paired DESeq2 analysis.
read_plot_inputs <- function(results_path, vst_counts_path, samples_path) {
  results <- read.csv(results_path, check.names = FALSE, stringsAsFactors = FALSE)
  vst_counts <- read.csv(vst_counts_path, check.names = FALSE, stringsAsFactors = FALSE)
  samples <- read.delim(samples_path, check.names = FALSE, stringsAsFactors = FALSE)
  required_results <- c("gene_id", "gene_symbol", "baseMean", "log2FoldChange", "padj")
  if (!all(required_results %in% colnames(results))) stop("DESeq2 result table is missing required plotting columns.")
  if (!all(c("gene_id", "gene_symbol", samples$sample_name) %in% colnames(vst_counts))) stop("VST matrix is missing required plotting columns.")
  list(results = results, vst_counts = vst_counts, samples = samples)
}

gene_labels <- function(rows) {
  ifelse(is.na(rows$gene_symbol) | rows$gene_symbol == "", rows$gene_id, rows$gene_symbol)
}

de_direction <- function(rows) {
  ifelse(!is.na(rows$padj) & !is.na(rows$log2FoldChange) & rows$padj < 0.05 & rows$log2FoldChange >= 1, "Up",
    ifelse(!is.na(rows$padj) & !is.na(rows$log2FoldChange) & rows$padj < 0.05 & rows$log2FoldChange <= -1, "Down", "Not significant"))
}

write_significant_genes <- function(results, output_dir) {
  significant <- results[!is.na(results$padj) & results$padj < 0.05, , drop = FALSE]
  write.csv(significant, file.path(output_dir, "deseq2_significant_genes_padj_0.05.csv"), row.names = FALSE)
}

plot_volcano <- function(results, groups, colors, output_dir) {
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
    if (length(selected)) text(results$log2FoldChange[selected], -log10(pmax(results$padj[selected], .Machine$double.xmin)), gene_labels(results)[selected], cex = .55, pos = 3)
  }
  dev.off()
}

plot_ma <- function(results, groups, colors, output_dir) {
  finite <- is.finite(results$baseMean) & results$baseMean > 0 & is.finite(results$log2FoldChange)
  png(file.path(output_dir, "ma_plot.png"), width = 2550, height = 1650, res = 300)
  plot(log2(results$baseMean[finite]), results$log2FoldChange[finite], pch = 16, cex = .4, col = colors[groups[finite]],
       xlab = "log2 mean normalized expression", ylab = "log2 fold change", main = "MA plot: tumor vs adjacent normal")
  abline(h = c(-1, 0, 1), lty = c(2, 1, 2), lwd = c(.6, .7, .6))
  legend("topright", legend = names(colors), col = colors, pch = 16, bty = "n")
  dev.off()
}

prepare_vst_matrix <- function(vst_counts, samples) {
  sample_names <- samples$sample_name
  matrix <- as.matrix(vst_counts[, sample_names, drop = FALSE])
  storage.mode(matrix) <- "double"
  matrix
}

plot_pca <- function(matrix, samples, output_dir) {
  sample_names <- samples$sample_name
  pca <- prcomp(t(matrix), scale. = FALSE)
  condition_colors <- c(normal = "#222222", tumor = "#F5A000")
  png(file.path(output_dir, "pca_vst.png"), width = 1800, height = 1500, res = 300)
  plot(pca$x[, 1], pca$x[, 2], pch = 16, cex = 1.3, col = condition_colors[samples$condition],
       xlab = sprintf("PC1 (%.1f%%)", 100 * summary(pca)$importance[2, 1]),
       ylab = sprintf("PC2 (%.1f%%)", 100 * summary(pca)$importance[2, 2]), main = "PCA of VST expression")
  text(pca$x[, 1], pca$x[, 2], labels = sample_names, pos = 3, cex = .65)
  legend("topright", legend = names(condition_colors), col = condition_colors, pch = 16, bty = "n")
  dev.off()
}

plot_sample_distance_heatmap <- function(matrix, output_dir) {
  png(file.path(output_dir, "sample_distance_heatmap.png"), width = 1950, height = 1740, res = 300)
  heatmap(as.matrix(dist(t(matrix))), symm = TRUE, col = hcl.colors(100, "YlOrRd", rev = TRUE), main = "Sample distances (VST)")
  dev.off()
}

plot_top_de_gene_heatmap <- function(results, vst_counts, matrix, output_dir) {
  ranked <- results[is.finite(results$padj), , drop = FALSE]
  ranked <- ranked[order(ranked$padj, -abs(ranked$log2FoldChange)), , drop = FALSE]
  ranked <- ranked[seq_len(min(20, nrow(ranked))), , drop = FALSE]
  indices <- match(ranked$gene_id, vst_counts$gene_id)
  found <- !is.na(indices)
  indices <- indices[found]
  if (length(indices) >= 2) {
    top_matrix <- matrix[indices, , drop = FALSE]
    rownames(top_matrix) <- make.unique(gene_labels(ranked)[found])
    row_z <- t(scale(t(top_matrix)))
    row_z[!is.finite(row_z)] <- 0
    png(file.path(output_dir, "top_de_gene_heatmap.png"), width = 2460, height = 2640, res = 300)
    heatmap(row_z, col = hcl.colors(101, "Blue-Red 3"), scale = "none", margins = c(8, 11),
            main = "Top 20 differentially expressed genes (row z-scores)")
    dev.off()
  }
}

plot_de_outputs <- function(results_path, vst_counts_path, samples_path) {
  output_dir <- dirname(results_path)
  inputs <- read_plot_inputs(results_path, vst_counts_path, samples_path)
  colors <- c("Down" = "#5B8DD9", "Not significant" = "#BDBDBD", "Up" = "#B7352D")
  groups <- de_direction(inputs$results)
  matrix <- prepare_vst_matrix(inputs$vst_counts, inputs$samples)
  write_significant_genes(inputs$results, output_dir)
  plot_volcano(inputs$results, groups, colors, output_dir)
  plot_ma(inputs$results, groups, colors, output_dir)
  plot_pca(matrix, inputs$samples, output_dir)
  plot_sample_distance_heatmap(matrix, output_dir)
  plot_top_de_gene_heatmap(inputs$results, inputs$vst_counts, matrix, output_dir)
}
