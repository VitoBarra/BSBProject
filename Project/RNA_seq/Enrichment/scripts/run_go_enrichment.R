#!/usr/bin/env Rscript

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
selected_genes_path <- args[["selected-genes"]]
universe_genes_path <- args[["universe-genes"]]
outdir <- args[["outdir"]]

if (is.null(selected_genes_path) || is.null(universe_genes_path) || is.null(outdir)) {
  stop(
    "Usage: run_go_enrichment.R --selected-genes selected.txt ",
    "--universe-genes universe.txt --outdir enrichment_dir"
  )
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

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

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

write.csv(all_results, file.path(outdir, "go_overrepresentation_all.csv"), row.names = FALSE)
