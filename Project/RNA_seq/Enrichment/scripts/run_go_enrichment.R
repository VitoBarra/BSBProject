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
de_results_path <- args[["de-results"]]
outdir <- args[["outdir"]]
padj_cutoff <- as.numeric(args[["padj-cutoff"]] %||% "0.05")
lfc_cutoff <- as.numeric(args[["lfc-cutoff"]] %||% "0")

if (is.null(de_results_path) || is.null(outdir)) {
  stop("Usage: run_go_enrichment.R --de-results deseq2_all_genes.csv --outdir enrichment_dir")
}

require_package("clusterProfiler")
require_package("org.Hs.eg.db")
require_package("AnnotationDbi")
require_package("ggplot2")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

de <- read.csv(de_results_path, stringsAsFactors = FALSE)
required_cols <- c("gene_id", "padj", "log2FoldChange")
missing_cols <- setdiff(required_cols, colnames(de))
if (length(missing_cols) > 0) {
  stop("DE table is missing columns: ", paste(missing_cols, collapse = ", "))
}

clean_ensembl <- sub("\\..*$", "", de$gene_id)
sig <- de[
  !is.na(de$padj) &
    de$padj < padj_cutoff &
    !is.na(de$log2FoldChange) &
    abs(de$log2FoldChange) >= lfc_cutoff,
]
sig_ensembl <- unique(sub("\\..*$", "", sig$gene_id))
universe_ensembl <- unique(clean_ensembl[!is.na(de$padj)])

run_one <- function(ontology) {
  clusterProfiler::enrichGO(
    gene = sig_ensembl,
    universe = universe_ensembl,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )
}

ego_list <- list(BP = run_one("BP"), MF = run_one("MF"), CC = run_one("CC"))
all_results <- do.call(
  rbind,
  lapply(names(ego_list), function(ont) {
    df <- as.data.frame(ego_list[[ont]])
    if (nrow(df) == 0) {
      return(NULL)
    }
    df$ontology <- ont
    df
  })
)

if (is.null(all_results)) {
  all_results <- data.frame()
}

if (nrow(all_results) > 0) {
  all_results <- all_results[order(all_results$p.adjust, all_results$pvalue, na.last = TRUE), ]
}

write.csv(all_results, file.path(outdir, "go_overrepresentation_all.csv"), row.names = FALSE)
significant <- all_results[!is.na(all_results$p.adjust) & all_results$p.adjust < 0.05, , drop = FALSE]
write.csv(significant, file.path(outdir, "go_overrepresentation_significant.csv"), row.names = FALSE)

summary_path <- file.path(outdir, "enrichment_summary.txt")
sink(summary_path)
cat("GO over-representation analysis\n")
cat("Input DE table:", de_results_path, "\n")
cat("padj cutoff:", padj_cutoff, "\n")
cat("absolute log2FC cutoff:", lfc_cutoff, "\n")
cat("Tested genes:", length(universe_ensembl), "\n")
cat("Significant genes:", length(sig_ensembl), "\n")
cat("GO terms tested:", nrow(all_results), "\n")
cat("GO terms significant at padj < 0.05:", nrow(significant), "\n")
sink()

if (nrow(significant) > 0) {
  top_terms <- head(significant, 20)
  top_terms$Description <- factor(top_terms$Description, levels = rev(top_terms$Description))
  plot <- ggplot2::ggplot(
    top_terms,
    ggplot2::aes(x = Description, y = -log10(p.adjust), size = Count, color = ontology)
  ) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "-log10 adjusted p-value", size = "Genes", color = "Ontology") +
    ggplot2::theme_minimal(base_size = 10)
  ggplot2::ggsave(file.path(outdir, "go_overrepresentation_dotplot.pdf"), plot, width = 8, height = 6)
}
