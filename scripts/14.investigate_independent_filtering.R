#!/usr/bin/env Rscript

# Study 1: Investigate DESeq2 independent filtering
# Purpose:
#   Examine the independent-filtering threshold and determine
#   which genes have missing adjusted p-values.
#
# This script reads existing results and does not refit DESeq2.
# It does not modify Script 10 or overwrite existing outputs.

# ------------------------------------------------------------------
# 1. Project paths
# ------------------------------------------------------------------

project_root <- getwd()

results_file <- file.path(
  project_root,
  "analysis",
  "study1",
  "study1_deseq2_results.csv"
)

output_dir <- file.path(
  project_root,
  "analysis",
  "study1",
  "filtering_investigation"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(results_file)) {
  stop("DESeq2 results file not found: ", results_file)
}

# ------------------------------------------------------------------
# 2. Load existing results
# ------------------------------------------------------------------

results_df <- read.csv(
  results_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "GeneID",
  "baseMean",
  "pvalue",
  "padj"
)

if (!all(required_columns %in% names(results_df))) {
  stop("The results file is missing required columns.")
}

# ------------------------------------------------------------------
# 3. Summarise missing values
# ------------------------------------------------------------------

total_genes <- nrow(results_df)

missing_pvalue <- is.na(results_df$pvalue)
missing_padj <- is.na(results_df$padj)

valid_pvalue_missing_padj <- (
  !missing_pvalue &
    missing_padj
)

# ------------------------------------------------------------------
# 4. Examine expression levels by adjusted-p-value availability
# ------------------------------------------------------------------

expression_summary <- data.frame(
  Group = c(
    "Adjusted p-value available",
    "Raw p-value available but adjusted p-value missing",
    "Raw p-value missing"
  ),
  Gene_count = c(
    sum(!missing_padj),
    sum(valid_pvalue_missing_padj),
    sum(missing_pvalue)
  ),
  Median_baseMean = c(
    median(results_df$baseMean[!missing_padj], na.rm = TRUE),
    median(
      results_df$baseMean[valid_pvalue_missing_padj],
      na.rm = TRUE
    ),
    median(results_df$baseMean[missing_pvalue], na.rm = TRUE)
  ),
  Maximum_baseMean = c(
    max(results_df$baseMean[!missing_padj], na.rm = TRUE),
    max(
      results_df$baseMean[valid_pvalue_missing_padj],
      na.rm = TRUE
    ),
    max(results_df$baseMean[missing_pvalue], na.rm = TRUE)
  )
)

write.csv(
  expression_summary,
  file.path(output_dir, "study1_filtering_expression_summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 5. Save genes with missing adjusted p-values
# ------------------------------------------------------------------

missing_padj_genes <- results_df[
  missing_padj,
  ,
  drop = FALSE
]

write.csv(
  missing_padj_genes,
  file.path(output_dir, "study1_missing_adjusted_pvalues.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 6. Summarise the findings
# ------------------------------------------------------------------

summary_lines <- c(
  "STUDY 1: INDEPENDENT-FILTERING INVESTIGATION",
  "============================================",
  "",
  paste("Analysis date:", as.character(Sys.Date())),
  paste("Total genes in results:", total_genes),
  paste("Genes with adjusted p-values:", sum(!missing_padj)),
  paste("Genes with missing adjusted p-values:", sum(missing_padj)),
  paste("Genes with missing raw p-values:", sum(missing_pvalue)),
  paste(
    "Genes with raw p-values but missing adjusted p-values:",
    sum(valid_pvalue_missing_padj)
  ),
  "",
  "Expression summary by adjusted-p-value availability:",
  capture.output(print(expression_summary, row.names = FALSE)),
  "",
  "Important limitation:",
  paste(
    "This script examines the existing results and expression levels.",
    "It does not recover DESeq2's exact independent-filtering",
    "threshold because that threshold is not recorded in the CSV."
  ),
  paste(
    "To recover the exact threshold, the fitted DESeqResults object",
    "or a controlled rerun of the DESeq2 results calculation is needed."
  )
)

writeLines(
  summary_lines,
  file.path(output_dir, "study1_filtering_investigation_report.txt")
)

cat(
  paste(summary_lines, collapse = "\n"),
  "\n\nFiltering investigation completed successfully.\n"
)

cat("Output directory: ", output_dir, "\n", sep = "")