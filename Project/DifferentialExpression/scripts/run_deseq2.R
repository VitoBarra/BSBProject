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

require_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Missing required R package '", pkg, "'. Install it with BiocManager before running this step.",
      call. = FALSE
    )
  }
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

require_package("DESeq2")
require_package("ggplot2")
require_package("SummarizedExperiment")

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
write.csv(
  res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ],
  file.path(outdir, "deseq2_significant_genes_padj_0.05.csv"),
  row.names = FALSE
)
ranked <- res_df[!is.na(res_df$stat), ]
ranked <- ranked[order(ranked$stat, decreasing = TRUE), ]
write.csv(ranked, file.path(outdir, "deseq2_ranked_genes.csv"), row.names = FALSE)

normalized_counts <- as.data.frame(DESeq2::counts(dds, normalized = TRUE))
normalized_counts$gene_id <- rownames(normalized_counts)
normalized_counts <- merge(gene_annot, normalized_counts, by = "gene_id", all.y = TRUE, sort = FALSE)
write.csv(normalized_counts, file.path(outdir, "normalized_counts.csv"), row.names = FALSE)

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

volcano_df <- res_df
volcano_df$neg_log10_padj <- -log10(pmax(volcano_df$padj, .Machine$double.xmin))
volcano_df$direction <- "Not significant"
volcano_df$direction[
  !is.na(volcano_df$padj) &
    volcano_df$padj < 0.05 &
    !is.na(volcano_df$log2FoldChange) &
    volcano_df$log2FoldChange >= 1
] <- "Up"
volcano_df$direction[
  !is.na(volcano_df$padj) &
    volcano_df$padj < 0.05 &
    !is.na(volcano_df$log2FoldChange) &
    volcano_df$log2FoldChange <= -1
] <- "Down"
volcano_df$direction <- factor(volcano_df$direction, levels = c("Down", "Not significant", "Up"))
volcano_labels <- rbind(
  head(volcano_df[volcano_df$direction == "Down", ][order(volcano_df$padj[volcano_df$direction == "Down"]), ], 10),
  head(volcano_df[volcano_df$direction == "Up", ][order(volcano_df$padj[volcano_df$direction == "Up"]), ], 10)
)
volcano_labels$plot_label <- ifelse(
  is.na(volcano_labels$gene_symbol) | volcano_labels$gene_symbol == "",
  volcano_labels$gene_id,
  volcano_labels$gene_symbol
)

volcano <- ggplot2::ggplot(volcano_df, ggplot2::aes(x = log2FoldChange, y = neg_log10_padj, color = direction)) +
  ggplot2::geom_point(alpha = 0.75, size = 1.25, na.rm = TRUE) +
  ggplot2::scale_color_manual(values = c("Down" = "#5B8DD9", "Not significant" = "grey75", "Up" = "#B7352D")) +
  ggplot2::geom_vline(xintercept = c(-1, 1), linewidth = 0.35, linetype = "dashed") +
  ggplot2::geom_hline(yintercept = -log10(0.05), linewidth = 0.3, linetype = "dashed") +
  ggplot2::geom_label(
    data = volcano_labels,
    ggplot2::aes(label = plot_label),
    size = 3,
    label.padding = grid::unit(0.12, "lines"),
    show.legend = FALSE
  ) +
  ggplot2::labs(
    title = "Differential expression: tumor vs adjacent normal",
    subtitle = "Paired DESeq2 model; padj < 0.05 and |log2FC| >= 1",
    x = "log2 fold change",
    y = "-log10 adjusted p-value",
    color = NULL
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(legend.position = "right")
ggplot2::ggsave(file.path(outdir, "volcano_padj.pdf"), volcano, width = 7.2, height = 6.2)
ggplot2::ggsave(file.path(outdir, "volcano_padj.png"), volcano, width = 7.2, height = 6.2, dpi = 300)

ma_df <- res_df[
  !is.na(res_df$baseMean) &
    res_df$baseMean > 0 &
    !is.na(res_df$log2FoldChange),
]
ma_df$log2_mean_expression <- log2(ma_df$baseMean)
ma_df$direction <- "NS"
ma_df$direction[!is.na(ma_df$padj) & ma_df$padj < 0.05 & ma_df$log2FoldChange >= 1] <- "Up"
ma_df$direction[!is.na(ma_df$padj) & ma_df$padj < 0.05 & ma_df$log2FoldChange <= -1] <- "Down"
ma_df$direction <- factor(ma_df$direction, levels = c("Up", "Down", "NS"))
ma_labels <- rbind(
  head(ma_df[ma_df$direction == "Down", ][order(ma_df$padj[ma_df$direction == "Down"]), ], 9),
  head(ma_df[ma_df$direction == "Up", ][order(ma_df$padj[ma_df$direction == "Up"]), ], 9)
)
ma_labels$plot_label <- ifelse(
  is.na(ma_labels$gene_symbol) | ma_labels$gene_symbol == "",
  ma_labels$gene_id,
  ma_labels$gene_symbol
)
direction_counts <- table(ma_df$direction)
ma_legend_labels <- c(
  paste0("Up: ", direction_counts[["Up"]]),
  paste0("Down: ", direction_counts[["Down"]]),
  "NS"
)

ma_plot <- ggplot2::ggplot(
  ma_df,
  ggplot2::aes(x = log2_mean_expression, y = log2FoldChange, color = direction)
) +
  ggplot2::geom_point(alpha = 0.75, size = 0.8) +
  ggplot2::scale_color_manual(
    values = c("Up" = "#CC1F27", "Down" = "#2166AC", "NS" = "grey70"),
    labels = ma_legend_labels
  ) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.45) +
  ggplot2::geom_hline(yintercept = c(-1, 1), linewidth = 0.35, linetype = "dashed") +
  ggplot2::geom_text(
    data = ma_labels,
    ggplot2::aes(label = plot_label),
    size = 3,
    fontface = "bold",
    check_overlap = TRUE,
    vjust = -0.6,
    show.legend = FALSE
  ) +
  ggplot2::labs(
    title = "MA plot: tumor vs adjacent normal",
    subtitle = "Paired DESeq2 model; padj < 0.05 and |log2FC| >= 1",
    x = "log2 mean normalized expression",
    y = "log2 fold change",
    color = NULL
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(legend.position = "top")
ggplot2::ggsave(file.path(outdir, "ma_plot.pdf"), ma_plot, width = 8.5, height = 5.5)
ggplot2::ggsave(file.path(outdir, "ma_plot.png"), ma_plot, width = 8.5, height = 5.5, dpi = 300)

vst <- DESeq2::vst(dds, blind = FALSE)

pdf(file.path(outdir, "pca_vst.pdf"), width = 6, height = 5)
print(DESeq2::plotPCA(vst, intgroup = c("patient", "condition")) + ggplot2::theme_minimal(base_size = 11))
dev.off()

sample_dist <- as.matrix(dist(t(SummarizedExperiment::assay(vst))))
sample_dist_long <- as.data.frame(as.table(sample_dist))
colnames(sample_dist_long) <- c("sample_1", "sample_2", "distance")
heatmap_plot <- ggplot2::ggplot(sample_dist_long, ggplot2::aes(sample_1, sample_2, fill = distance)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient(low = "white", high = "#2166AC") +
  ggplot2::coord_fixed() +
  ggplot2::labs(title = "Sample distance heatmap", x = NULL, y = NULL, fill = "Distance") +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
ggplot2::ggsave(file.path(outdir, "sample_distance_heatmap.pdf"), heatmap_plot, width = 6.5, height = 5.8)

top_de <- res_df[!is.na(res_df$padj), ]
top_de <- top_de[order(top_de$padj, -abs(top_de$log2FoldChange), na.last = TRUE), ]
top_de <- head(top_de$gene_id, 50)
if (length(top_de) > 1) {
  vst_matrix <- SummarizedExperiment::assay(vst)
  top_de <- head(top_de, 20)
  top_matrix <- vst_matrix[top_de[top_de %in% rownames(vst_matrix)], , drop = FALSE]
  row_z <- t(scale(t(top_matrix)))
  row_z[is.na(row_z)] <- 0
  gene_labels <- gene_annot$gene_symbol[match(rownames(row_z), gene_annot$gene_id)]
  gene_labels[is.na(gene_labels) | gene_labels == ""] <- rownames(row_z)[is.na(gene_labels) | gene_labels == ""]
  rownames(row_z) <- make.unique(gene_labels)
  condition_colors <- c(normal = "#222222", tumor = "#F5A000")
  column_side_colors <- condition_colors[as.character(samples[colnames(row_z), "condition"])]
  draw_heatmap <- function() {
    heatmap(
      row_z,
      Rowv = TRUE,
      Colv = TRUE,
      scale = "none",
      col = grDevices::colorRampPalette(c("#173B9A", "white", "#F44336"))(101),
      breaks = seq(-2.5, 2.5, length.out = 102),
      ColSideColors = column_side_colors,
      margins = c(9, 8),
      cexRow = 0.8,
      cexCol = 0.75,
      main = "Top 20 differentially expressed genes (row z-scores)"
    )
    legend(
      "topright",
      legend = c("Normal", "Tumor"),
      fill = condition_colors[c("normal", "tumor")],
      title = "Condition",
      border = NA,
      bty = "n",
      xpd = TRUE
    )
  }
  pdf(file.path(outdir, "top_de_gene_heatmap.pdf"), width = 7.5, height = 8)
  draw_heatmap()
  dev.off()
  png(file.path(outdir, "top_de_gene_heatmap.png"), width = 2250, height = 2400, res = 300)
  draw_heatmap()
  dev.off()
}
