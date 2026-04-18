#!/usr/bin/env Rscript

require_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
  {
    stop(
      "Missing required R package '", pkg, "'. Install it with BiocManager before running this step.",
      call. = FALSE
    )
  }
}

require_package("optparse")
parser <- optparse::OptionParser(option_list = list(
  optparse::make_option("--results-path",    dest = "results_path", type = "character"),
  optparse::make_option("--de-results-path", dest = "de_results_path", type = "character", default = ""),
  optparse::make_option("--padj-cutoff",     dest = "padj_cutoff", type = "double", default = 0.05),
  optparse::make_option("--lfc-cutoff",      dest = "lfc_cutoff", type = "double", default = 0),
  optparse::make_option("--plots-only",      dest = "plots_only", action = "store_true", default = FALSE)
))
args <- optparse::parse_args(parser)
results_path        <- args$results_path

validate_required_paths <- function(results_path, de_results_path, parser) {
  if (is.null(results_path) || !nzchar(de_results_path)) {
    optparse::print_help(parser)
    stop("--results-path and --de-results-path are required", call. = FALSE)
  }
}

select_go_genes <- function(de_results_path, padj_cutoff, lfc_cutoff) {
  # Check if the DESeq2 results file exists and read it
  if (!file.exists(de_results_path)) stop("Missing DESeq2 result table: ", de_results_path)
  # Read the DESeq2 results and validate required columns
  de_results <- read.csv(de_results_path, check.names = FALSE, stringsAsFactors = FALSE)
  required_columns <- c("gene_id", "padj", "log2FoldChange")
  missing_columns <- setdiff(required_columns, names(de_results))
  # If any required columns are missing, stop execution and report the missing columns
  if (length(missing_columns) > 0)
  {
    stop("DESeq2 results are missing columns: ", paste(missing_columns, collapse = ", "))
  }
  # Keep the selected genes and enrichment background in memory.
  genes <- sub("\\..*$", "", trimws(as.character(de_results$gene_id)))
  valid_genes <- !is.na(genes) & nzchar(genes)
  # Create a unique list of gene identifiers for the universe and selected genes based on the specified cutoffs
  universe_genes <- unique(genes[valid_genes])
  padj <- suppressWarnings(as.numeric(de_results$padj))
  log2fc <- suppressWarnings(as.numeric(de_results$log2FoldChange))
  selected_genes <- unique(genes[valid_genes & is.finite(padj) & padj < padj_cutoff &
                                 is.finite(log2fc) & abs(log2fc) >= lfc_cutoff])

  if (length(universe_genes) == 0) stop("No gene identifiers are available for the GO universe")
  if (length(selected_genes) == 0) stop("No genes pass the configured adjusted p-value and log2 fold-change cutoffs")

  list(selected = selected_genes, universe = universe_genes)
}

write_enrichment_summary <- function(summary_path, de_results_path, padj_cutoff, lfc_cutoff,
                                     universe_genes, selected_genes, all_results, significant) {
  cat(
    "GO over-representation analysis\n",
    "Input DE table: ", de_results_path, "\n",
    "padj cutoff: ", padj_cutoff, "\n",
    "absolute log2FC cutoff: ", lfc_cutoff, "\n",
    "Tested genes: ", length(universe_genes), "\n",
    "Significant genes: ", length(selected_genes), "\n",
    "GO terms tested: ", nrow(all_results), "\n",
    "GO terms significant at padj < 0.05: ", nrow(significant), "\n",
    file = summary_path, sep = ""
  )
}

plot_go_dotplot <- function(significant, plot_path)
{
  if (nrow(significant) == 0)
  {
    if (file.exists(plot_path)) unlink(plot_path)
    # If there are no significant GO terms skipp plotting and return NULL to indicate that no plot was created
    return(invisible(NULL))
  }

  # Select the top 20 significant GO terms based on adjusted p-value and prepare data for plotting
  top_terms <- head(significant, 20)
  descriptions <- vapply(top_terms$Description, function(x) paste(strwrap(x, width = 42), collapse = "\n"), character(1))
  adjusted <- -log10(pmax(top_terms$p.adjust, .Machine$double.xmin))
  counts <- top_terms$Count
  sizes <- 1.2 + 2.8 * (counts - min(counts)) / max(1, max(counts) - min(counts))
  colors <- c(BP = "#4C78A8", MF = "#F58518", CC = "#54A24B")

  # Create a dot plot of the top significant GO terms, with point sizes representing gene counts and colors representing ontology categories
  png(plot_path, width = 2850, height = 2100, res = 300)
  par(mar = c(5, 15, 3, 2))
  plot(adjusted, rev(seq_along(adjusted)), pch = 16, cex = sizes, col = colors[top_terms$ontology],
       xlab = "-log10 adjusted p-value", ylab = "", yaxt = "n", main = "GO over-representation analysis: top significant terms")
  axis(2, at = rev(seq_along(adjusted)), labels = descriptions, las = 1, cex.axis = .72)
  legend("bottomright", legend = names(colors)[names(colors) %in% top_terms$ontology], col = colors[names(colors) %in% top_terms$ontology], pch = 16, bty = "n", title = "Ontology")
  dev.off()
}

write_go_outputs <- function(all_results, output_dir, results_path, de_results_path, padj_cutoff,
                             lfc_cutoff, universe_genes, selected_genes) {
  all_results <- all_results[order(all_results$p.adjust, all_results$pvalue, na.last = TRUE), , drop = FALSE]
  write.csv(all_results, results_path, row.names = FALSE)
  significant <- all_results[is.finite(all_results$p.adjust) & all_results$p.adjust < 0.05, , drop = FALSE]
  write.csv(significant, file.path(output_dir, "go_overrepresentation_significant.csv"), row.names = FALSE)
  write_enrichment_summary(file.path(output_dir, "enrichment_summary.txt"), de_results_path, padj_cutoff,
                           lfc_cutoff, universe_genes, selected_genes, all_results, significant)
  plot_go_dotplot(significant, file.path(output_dir, "go_overrepresentation_dotplot.png"))
}

run_go_for_ontology <- function(ontology, selected_genes, universe_genes) {
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

combine_go_results <- function(enrichment_by_ontology)
{
  result_tables <- lapply(names(enrichment_by_ontology), function(ontology)
  {
    result <- as.data.frame(enrichment_by_ontology[[ontology]])
    if (nrow(result) == 0) {
      return(NULL)
    }
    #  Add the ontology column to the result data frame
    result$ontology <- ontology
    result
  })
  # Remove NULL dataframe from the list of result tables and combine them into a single data frame
  result_tables <- Filter(Negate(is.null), result_tables)
  all_results <- if (length(result_tables) > 0) do.call(rbind, result_tables) else NULL

  # If no results were found for any ontology, create an empty data frame with the expected columns
  if (is.null(all_results)) {
    all_results <- data.frame(
      ID = character(), Description = character(), GeneRatio = character(), BgRatio = character(),
      RichFactor = numeric(), FoldEnrichment = numeric(), zScore = numeric(), pvalue = numeric(),
      p.adjust = numeric(), qvalue = numeric(), geneID = character(), Count = integer(),
      ontology = character()
    )
  }
  all_results
}

validate_required_paths(results_path, args$de_results_path, parser)
output_dir <- dirname(results_path)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
go_inputs <- select_go_genes(args$de_results_path, args$padj_cutoff, args$lfc_cutoff)
selected_genes <- go_inputs$selected
universe_genes <- go_inputs$universe

if (isTRUE(args$plots_only)) {
  all_results <- read.csv(results_path, check.names = FALSE, stringsAsFactors = FALSE)
} else {
  require_package("clusterProfiler")
  require_package("org.Hs.eg.db")

  ontology_names <- c("BP", "MF", "CC")
  enrichment_by_ontology <- lapply(ontology_names, run_go_for_ontology,
                                   selected_genes = selected_genes, universe_genes = universe_genes)
  names(enrichment_by_ontology) <- ontology_names
  all_results <- combine_go_results(enrichment_by_ontology)
}

write_go_outputs(all_results, output_dir, results_path, args$de_results_path, args$padj_cutoff,
                 args$lfc_cutoff, universe_genes, selected_genes)
