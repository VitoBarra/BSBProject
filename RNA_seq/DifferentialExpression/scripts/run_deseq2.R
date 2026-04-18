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

# Check that required arguments are provided
if (is.null(samples_path) || is.null(transcript_to_gene_map_path) || is.null(results_path) || is.null(vst_counts_path) || is.null(summary_path)) {
  optparse::print_help(parser)
  stop("--samples, --transcript-to-gene-map, --results-path, --vst-counts-path, and --summary-path are required", call. = FALSE)
}

# check that output directories exist, and create them if they don't
for (output_path in c(results_path, vst_counts_path, summary_path)) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
}

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_argument)))
source(file.path(script_dir, "deseq2_plots.R"))

if (isTRUE(args$plots_only)) {
  plot_de_outputs(results_path, vst_counts_path, samples_path)
  quit(save = "no")
}

# read the dataframe form the sample sheet and transcript-to-gene map
# - The check.names  = FALSE argument prevents R from modifying column names
# - stringsAsFactors = FALSE ensures that string columns are not automatically converted to factors (hard categorical).
samples <- read.delim(samples_path, check.names = FALSE, stringsAsFactors = FALSE)
transcript_to_gene_map <- read.delim(transcript_to_gene_map_path, check.names = FALSE, stringsAsFactors = FALSE)


# Validate that the sample sheet and transcript-to-gene map contain the required columns
required_sample_cols <- c("sample_name", "patient", "condition")
missing_sample_cols <- setdiff(required_sample_cols, colnames(samples))
if (length(missing_sample_cols) > 0)
{
  stop("Sample table is missing columns: ", paste(missing_sample_cols, collapse = ", "))
}

required_transcript_to_gene_map_columns <- c("transcript_id", "gene_id", "gene_symbol")
missing_transcript_to_gene_map_columns <- setdiff(required_transcript_to_gene_map_columns, colnames(transcript_to_gene_map))
if (length(missing_transcript_to_gene_map_columns) > 0)
{
  stop("Transcript-to-gene map is missing columns: ", paste(missing_transcript_to_gene_map_columns, collapse = ", "))
}

# Validate that the sample sheet contains non-empty, unique sample names, and that each patient has exactly one normal and one tumor sample
sample_names <- samples$sample_name
if (anyNA(samples[, required_sample_cols]) || any(!nzchar(trimws(unlist(samples[, required_sample_cols])))))
{
  stop("Sample names, patients, and conditions must all be non-empty.")
}

# validate if sample names are unique
if (anyDuplicated(sample_names)) {
  stop("Sample names must be unique: ", paste(unique(sample_names[duplicated(sample_names)]), collapse = ", "))
}
# Validate that the sample sheet contains exactly two conditions: "normal" and "tumor"
expected_conditions <- c("normal", "tumor")
observed_conditions <- sort(unique(samples$condition))
if (!identical(observed_conditions, expected_conditions)) {
  stop(
    "Expected exactly the conditions 'normal' and 'tumor'; found: ",
    paste(observed_conditions, collapse = ", ")
  )
}

# Validate that each patient has exactly one normal and one tumor sample
patient_condition_counts <- table(samples$patient, samples$condition)
if (any(patient_condition_counts != 1L)) {
  invalid_patients <- rownames(patient_condition_counts)[apply(patient_condition_counts != 1L, 1, any)]
  stop(
    "Paired design requires exactly one normal and one tumor sample per patient. Invalid patients: ",
    paste(invalid_patients, collapse = ", ")
  )
}

# Validate that the Salmon quant.sf files exist for each sample
quant_files <- setNames(samples$quant_sf, sample_names)
missing_quant_files <- quant_files[!file.exists(quant_files)]
if (length(missing_quant_files) > 0)
{
  stop("Missing Salmon quant.sf files: ", paste(missing_quant_files, collapse = ", "))
}

transcript_to_gene_map_import <- unique(transcript_to_gene_map[, c("transcript_id", "gene_id")])
txi <- tximport::tximport(
  files = quant_files,
  type = "salmon",
  tx2gene = transcript_to_gene_map_import,
  ignoreTxVersion = TRUE
)

# Create a gene annotation dataframe with unique gene_id and gene_symbol pairs, ensuring no duplicates in gene_id
gene_annot <- unique(transcript_to_gene_map[, c("gene_id", "gene_symbol")])
gene_annot <- gene_annot[!duplicated(gene_annot$gene_id), , drop = FALSE]
rownames(gene_annot) <- gene_annot$gene_id

# Create a DESeq2 dataset using the tximport results and the sample metadata, specifying the design formula for paired analysis
rownames(samples) <- samples$sample_name
samples <- samples[sample_names, , drop = FALSE]
samples$patient <- factor(samples$patient)
samples$condition <- relevel(factor(samples$condition), ref = "normal")

# Filter out genes with low counts across samples based on the specified minimum count and minimum sample thresholds
keep <- rowSums(txi$counts >= min_count) >= min_samples
txi$counts <- txi$counts[keep, , drop = FALSE]
txi$abundance <- txi$abundance[keep, , drop = FALSE]
txi$length <- txi$length[keep, , drop = FALSE]
gene_annot <- gene_annot[rownames(txi$counts), , drop = FALSE]

# Check if any genes remain after filtering; if not, stop the analysis with an error message
if (nrow(txi$counts) == 0)
{
  stop("No genes remain after count filtering.")
}

# Create a DESeq2 dataset from the tximport results and sample metadata, specifying the design formula for paired analysis
dds <- DESeq2::DESeqDataSetFromTximport(
  txi = txi,
  colData = samples,
  design = ~ patient + condition
)

dds    <- DESeq2::DESeq(dds)
res    <- DESeq2::results(dds, contrast = c("condition", "tumor", "normal"))
# Order the results by adjusted p-value and p-value, and merge with gene annotation for output
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- merge(gene_annot, res_df, by = "gene_id", all.y = TRUE, sort = FALSE)
res_df <- res_df[, c("gene_id", "gene_symbol", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
res_df <- res_df[order(res_df$padj, res_df$pvalue, na.last = TRUE), ]
write.csv(res_df, results_path, row.names = FALSE)

# Write a summary of the DESeq2 analysis to a text file, including sample information and counts of retained and significant genes
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
