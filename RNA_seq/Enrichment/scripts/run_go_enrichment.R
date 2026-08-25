#!/usr/bin/env Rscript

require_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Missing required R package '", pkg, "'. Install it with BiocManager before running this step.",
      call. = FALSE
    )
  }
}

require_package("optparse")
parser <- optparse::OptionParser(option_list = list(
  optparse::make_option("--selected-genes", dest = "selected_genes", type = "character"),
  optparse::make_option("--universe-genes", dest = "universe_genes", type = "character"),
  optparse::make_option("--results-path", dest = "results_path", type = "character"),
  optparse::make_option("--de-results-path", dest = "de_results_path", type = "character", default = ""),
  optparse::make_option("--padj-cutoff", dest = "padj_cutoff", type = "double", default = 0.05),
  optparse::make_option("--lfc-cutoff", dest = "lfc_cutoff", type = "double", default = 0),
  optparse::make_option("--plots-only", dest = "plots_only", action = "store_true", default = FALSE)
))
args <- optparse::parse_args(parser)
selected_genes_path <- args$selected_genes
universe_genes_path <- args$universe_genes
results_path <- args$results_path

if (is.null(selected_genes_path) || is.null(universe_genes_path) || is.null(results_path)) {
  optparse::print_help(parser)
  stop("--selected-genes, --universe-genes, and --results-path are required", call. = FALSE)
}

require_package("clusterProfiler")
require_package("org.Hs.eg.db")

selected_genes <- unique(readLines(selected_genes_path, warn = FALSE))
universe_genes <- unique(readLines(universe_genes_path, warn = FALSE))
selected_genes <- selected_genes[nzchar(selected_genes)]
universe_genes <- universe_genes[nzchar(universe_genes)]
if (length(selected_genes) == 0) {
  stop("Selected gene list is empty")
}
if (length(universe_genes) == 0) {
  stop("Gene universe is empty")
}

dir.create(dirname(results_path), recursive = TRUE, showWarnings = FALSE)

write_go_outputs <- function(all_results) {
  output_dir <- dirname(results_path)
  all_results <- all_results[order(all_results$p.adjust, all_results$pvalue, na.last = TRUE), , drop = FALSE]
  write.csv(all_results, results_path, row.names = FALSE)
  significant <- all_results[is.finite(all_results$p.adjust) & all_results$p.adjust < 0.05, , drop = FALSE]
  write.csv(significant, file.path(output_dir, "go_overrepresentation_significant.csv"), row.names = FALSE)
  summary_path <- file.path(output_dir, "enrichment_summary.txt")
  cat(
    "GO over-representation analysis\n",
    "Input DE table: ", args$de_results_path, "\n",
    "padj cutoff: ", args$padj_cutoff, "\n",
    "absolute log2FC cutoff: ", args$lfc_cutoff, "\n",
    "Tested genes: ", length(universe_genes), "\n",
    "Significant genes: ", length(selected_genes), "\n",
    "GO terms tested: ", nrow(all_results), "\n",
    "GO terms significant at padj < 0.05: ", nrow(significant), "\n",
    file = summary_path, sep = ""
  )
  plot_path <- file.path(output_dir, "go_overrepresentation_dotplot.png")
  if (nrow(significant) == 0) {
    if (file.exists(plot_path)) unlink(plot_path)
    return(invisible(NULL))
  }
  top_terms <- head(significant, 20)
  descriptions <- vapply(top_terms$Description, function(x) paste(strwrap(x, width = 42), collapse = "\n"), character(1))
  adjusted <- -log10(pmax(top_terms$p.adjust, .Machine$double.xmin))
  counts <- top_terms$Count
  sizes <- 1.2 + 2.8 * (counts - min(counts)) / max(1, max(counts) - min(counts))
  colors <- c(BP = "#4C78A8", MF = "#F58518", CC = "#54A24B")
  png(plot_path, width = 2850, height = 2100, res = 300)
  par(mar = c(5, 15, 3, 2))
  plot(adjusted, rev(seq_along(adjusted)), pch = 16, cex = sizes, col = colors[top_terms$ontology],
       xlab = "-log10 adjusted p-value", ylab = "", yaxt = "n", main = "GO over-representation analysis: top significant terms")
  axis(2, at = rev(seq_along(adjusted)), labels = descriptions, las = 1, cex.axis = .72)
  legend("bottomright", legend = names(colors)[names(colors) %in% top_terms$ontology], col = colors[names(colors) %in% top_terms$ontology], pch = 16, bty = "n", title = "Ontology")
  dev.off()
}

if (isTRUE(args$plots_only)) {
  write_go_outputs(read.csv(results_path, check.names = FALSE, stringsAsFactors = FALSE))
  quit(save = "no")
}

run_one <- function(ontology) {
  clusterProfiler::enrichGO(
    gene = selected_genes,
    universe = universe_genes,
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
  lapply(names(ego_list), function(ontology) {
    result <- as.data.frame(ego_list[[ontology]])
    if (nrow(result) == 0) {
      return(NULL)
    }
    result$ontology <- ontology
    result
  })
)

if (is.null(all_results)) {
  all_results <- data.frame(
    ID = character(), Description = character(), GeneRatio = character(), BgRatio = character(),
    RichFactor = numeric(), FoldEnrichment = numeric(), zScore = numeric(), pvalue = numeric(),
    p.adjust = numeric(), qvalue = numeric(), geneID = character(), Count = integer(),
    ontology = character()
  )
}

write_go_outputs(all_results)
