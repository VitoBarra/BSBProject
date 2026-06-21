#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    "Usage: plot_de_results.R deseq2_all_genes.csv normalized_counts.csv sample_table.tsv output_dir",
    call. = FALSE
  )
}

results_path <- args[[1]]
normalized_counts_path <- args[[2]]
samples_path <- args[[3]]
outdir <- args[[4]]

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Missing required R package 'ggplot2'.", call. = FALSE)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
de <- read.csv(results_path, stringsAsFactors = FALSE)
normalized <- read.csv(normalized_counts_path, check.names = FALSE, stringsAsFactors = FALSE)
samples <- read.delim(samples_path, check.names = FALSE, stringsAsFactors = FALSE)
sample_names <- samples$sample_name

gene_label <- function(symbol, gene_id) {
  ifelse(is.na(symbol) | symbol == "", gene_id, symbol)
}

classify_direction <- function(padj, log2fc) {
  direction <- rep("Not significant", length(log2fc))
  direction[!is.na(padj) & padj < 0.05 & !is.na(log2fc) & log2fc >= 1] <- "Up"
  direction[!is.na(padj) & padj < 0.05 & !is.na(log2fc) & log2fc <= -1] <- "Down"
  factor(direction, levels = c("Down", "Not significant", "Up"))
}

volcano_df <- de
volcano_df$neg_log10_padj <- -log10(pmax(volcano_df$padj, .Machine$double.xmin))
volcano_df$direction <- classify_direction(volcano_df$padj, volcano_df$log2FoldChange)
volcano_labels <- rbind(
  head(volcano_df[volcano_df$direction == "Down", ][order(volcano_df$padj[volcano_df$direction == "Down"]), ], 10),
  head(volcano_df[volcano_df$direction == "Up", ][order(volcano_df$padj[volcano_df$direction == "Up"]), ], 10)
)
volcano_labels$plot_label <- gene_label(volcano_labels$gene_symbol, volcano_labels$gene_id)

volcano <- ggplot2::ggplot(
  volcano_df,
  ggplot2::aes(x = log2FoldChange, y = neg_log10_padj, color = direction)
) +
  ggplot2::geom_point(alpha = 0.75, size = 1.25, na.rm = TRUE) +
  ggplot2::scale_color_manual(
    values = c("Down" = "#5B8DD9", "Not significant" = "grey75", "Up" = "#B7352D")
  ) +
  ggplot2::geom_vline(xintercept = c(-1, 1), linewidth = 0.35, linetype = "dashed") +
  ggplot2::geom_hline(yintercept = -log10(0.05), linewidth = 0.35, linetype = "dashed") +
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
  ggplot2::theme_classic(base_size = 11)
ggplot2::ggsave(file.path(outdir, "volcano_padj.pdf"), volcano, width = 7.2, height = 6.2)
ggplot2::ggsave(file.path(outdir, "volcano_padj.png"), volcano, width = 7.2, height = 6.2, dpi = 300)

ma_df <- de[!is.na(de$baseMean) & de$baseMean > 0 & !is.na(de$log2FoldChange), ]
ma_df$log2_mean_expression <- log2(ma_df$baseMean)
ma_df$direction <- classify_direction(ma_df$padj, ma_df$log2FoldChange)
levels(ma_df$direction) <- c("Down", "NS", "Up")
ma_df$direction <- factor(ma_df$direction, levels = c("Up", "Down", "NS"))
ma_labels <- rbind(
  head(ma_df[ma_df$direction == "Down", ][order(ma_df$padj[ma_df$direction == "Down"]), ], 9),
  head(ma_df[ma_df$direction == "Up", ][order(ma_df$padj[ma_df$direction == "Up"]), ], 9)
)
ma_labels$plot_label <- gene_label(ma_labels$gene_symbol, ma_labels$gene_id)
direction_counts <- table(ma_df$direction)

ma_plot <- ggplot2::ggplot(
  ma_df,
  ggplot2::aes(x = log2_mean_expression, y = log2FoldChange, color = direction)
) +
  ggplot2::geom_point(alpha = 0.75, size = 0.8) +
  ggplot2::scale_color_manual(
    values = c("Up" = "#CC1F27", "Down" = "#2166AC", "NS" = "grey70"),
    labels = c(
      paste0("Up: ", direction_counts[["Up"]]),
      paste0("Down: ", direction_counts[["Down"]]),
      "NS"
    )
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

top_de <- de[!is.na(de$padj), ]
top_de <- head(top_de[order(top_de$padj, -abs(top_de$log2FoldChange)), ], 20)
matrix_rows <- match(top_de$gene_id, normalized$gene_id)
top_matrix <- as.matrix(normalized[matrix_rows, sample_names, drop = FALSE])
storage.mode(top_matrix) <- "numeric"
top_matrix <- log2(top_matrix + 1)
row_z <- t(scale(t(top_matrix)))
row_z[is.na(row_z)] <- 0
rownames(row_z) <- make.unique(gene_label(top_de$gene_symbol, top_de$gene_id))

condition_colors <- c(normal = "#222222", tumor = "#F5A000")
column_side_colors <- condition_colors[samples$condition[match(colnames(row_z), samples$sample_name)]]
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
