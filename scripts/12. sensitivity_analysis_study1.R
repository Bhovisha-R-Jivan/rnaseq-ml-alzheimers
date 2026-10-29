#!/usr/bin/env Rscript

# Study 1: Alzheimer's disease RNA-seq
# Purpose:
#   Assess the sensitivity of differential-expression results
#   to DESeq2 independent filtering.
#
# Design:
#   ~ AGE + sex + diagnosis
#
# Contrast:
#   Alzheimer's disease versus control
#
# Important:
#   This script does not modify Script 10 or overwrite its outputs.
#   Both analyses use the same samples, counts, and design.
#   Cook's-distance filtering remains enabled in both analyses.

# ------------------------------------------------------------------
# 1. Project paths
# ------------------------------------------------------------------

project_root <- getwd()

count_file <- file.path(
  project_root,
  "raw_data",
  "Study1_gene_counts.gz"
)

metadata_file <- file.path(
  project_root,
  "raw_data",
  "Study1_MetaData.csv"
)

output_dir <- file.path(
  project_root,
  "analysis",
  "study1",
  "sensitivity"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------------
# 2. Check inputs and packages
# ------------------------------------------------------------------

if (!file.exists(count_file)) {
  stop("Count matrix not found: ", count_file)
}

if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file)
}

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("DESeq2 is not installed in the active R environment.")
}

suppressPackageStartupMessages({
  library(DESeq2)
})

# ------------------------------------------------------------------
# 3. Load and validate the count matrix
# ------------------------------------------------------------------

counts_df <- read.delim(
  gzfile(count_file),
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!"GeneID" %in% names(counts_df)) {
  stop("The count matrix must contain a GeneID column.")
}

gene_ids <- counts_df$GeneID

if (anyNA(gene_ids) || anyDuplicated(gene_ids)) {
  stop("Gene IDs are missing or duplicated.")
}

rownames(counts_df) <- gene_ids
counts_df$GeneID <- NULL

counts <- as.matrix(counts_df)

if (!is.numeric(counts)) {
  stop("The count matrix contains non-numeric values.")
}

if (
  anyNA(counts) ||
  any(!is.finite(counts)) ||
  any(counts < 0) ||
  any(counts != floor(counts))
) {
  stop("Counts must be finite, non-negative integers.")
}

storage.mode(counts) <- "integer"

# ------------------------------------------------------------------
# 4. Load and validate metadata
# ------------------------------------------------------------------

metadata <- read.csv(
  metadata_file,
  sep = ";",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

required_columns <- c(
  "GEO_Accession (exp)",
  "diagnosis",
  "AGE",
  "sex"
)

if (!all(required_columns %in% names(metadata))) {
  stop("One or more required metadata columns are missing.")
}

sample_ids <- metadata[["GEO_Accession (exp)"]]

if (anyNA(sample_ids) || anyDuplicated(sample_ids)) {
  stop("Metadata sample IDs are missing or duplicated.")
}

if (
  anyDuplicated(colnames(counts)) ||
  !setequal(colnames(counts), sample_ids)
) {
  stop("Count matrix and metadata sample IDs do not match.")
}

# Align metadata to the count matrix.
metadata <- metadata[
  match(colnames(counts), sample_ids),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(counts),
    metadata[["GEO_Accession (exp)"]]
  )
)

# ------------------------------------------------------------------
# 5. Clean covariates and define reference levels
# ------------------------------------------------------------------

metadata$AGE <- trimws(as.character(metadata$AGE))

age_placeholders <- c(
  "",
  "na",
  "n/a",
  "not applicable",
  "not available",
  "unknown"
)

age_missing <- (
  is.na(metadata$AGE) |
    tolower(metadata$AGE) %in% age_placeholders
)

age_numeric <- suppressWarnings(
  as.numeric(metadata$AGE)
)

invalid_age <- (
  !age_missing &
    (is.na(age_numeric) | !is.finite(age_numeric))
)

if (any(invalid_age)) {
  stop("Unexpected non-numeric age values found.")
}

age_numeric[age_missing] <- NA_real_
metadata$AGE <- age_numeric

metadata$sex <- trimws(as.character(metadata$sex))

sex_placeholders <- c(
  "",
  "na",
  "n/a",
  "not applicable",
  "not available",
  "unknown"
)

sex_missing <- (
  is.na(metadata$sex) |
    tolower(metadata$sex) %in% sex_placeholders
)

metadata$sex[sex_missing] <- NA_character_

expected_diagnoses <- c(
  "Alzheimer's disease",
  "control"
)

if (
  anyNA(metadata$diagnosis) ||
    !all(metadata$diagnosis %in% expected_diagnoses)
) {
  stop("Missing or unexpected diagnosis labels found.")
}

observed_sexes <- unique(
  metadata$sex[!is.na(metadata$sex)]
)

if (!all(observed_sexes %in% c("female", "male"))) {
  stop("Unexpected sex labels found.")
}

metadata$diagnosis <- factor(
  metadata$diagnosis,
  levels = c("control", "Alzheimer's disease")
)

metadata$sex <- factor(
  metadata$sex,
  levels = c("female", "male")
)

# ------------------------------------------------------------------
# 6. Select complete cases
# ------------------------------------------------------------------

complete_case <- complete.cases(
  metadata[, c("AGE", "sex", "diagnosis")]
)

excluded_samples <- metadata[
  !complete_case,
  "GEO_Accession (exp)"
]

metadata <- metadata[
  complete_case,
  ,
  drop = FALSE
]

counts <- counts[
  ,
  complete_case,
  drop = FALSE
]

if (ncol(counts) != nrow(metadata)) {
  stop("Count matrix and metadata dimensions do not match.")
}

if (
  !identical(
    colnames(counts),
    metadata[["GEO_Accession (exp)"]]
  )
) {
  stop("Sample alignment failed after complete-case filtering.")
}

metadata$diagnosis <- droplevels(metadata$diagnosis)
metadata$sex <- droplevels(metadata$sex)

if (
  nlevels(metadata$diagnosis) != 2 ||
    nlevels(metadata$sex) != 2
) {
  stop("Both diagnosis groups and both sex categories are required.")
}

# ------------------------------------------------------------------
# 7. Apply the same zero-count filter and design as Script 10
# ------------------------------------------------------------------

genes_before_filtering <- nrow(counts)

keep_genes <- rowSums(counts) > 0

counts <- counts[
  keep_genes,
  ,
  drop = FALSE
]

genes_after_filtering <- nrow(counts)

design_formula <- ~ AGE + sex + diagnosis

design_matrix <- model.matrix(
  design_formula,
  data = metadata
)

if (qr(design_matrix)$rank != ncol(design_matrix)) {
  stop("The design matrix is not full rank.")
}

# ------------------------------------------------------------------
# 8. Fit the DESeq2 model
# ------------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = design_formula
)

dds <- DESeq(dds)

# ------------------------------------------------------------------
# 9. Extract results under two independent-filtering settings
# ------------------------------------------------------------------

# Analysis A: DESeq2 defaults.
# Independent filtering and Cook's-distance filtering are enabled.

res_default <- results(
  dds,
  contrast = c(
    "diagnosis",
    "Alzheimer's disease",
    "control"
  ),
  alpha = 0.05
)

# Analysis B: Disable independent filtering only.
# Cook's-distance filtering remains enabled.

res_no_independent_filtering <- results(
  dds,
  contrast = c(
    "diagnosis",
    "Alzheimer's disease",
    "control"
  ),
  alpha = 0.05,
  independentFiltering = FALSE
)

# ------------------------------------------------------------------
# 10. Assemble and save separate result tables
# ------------------------------------------------------------------

default_df <- as.data.frame(res_default)
default_df$GeneID <- rownames(default_df)

default_df <- default_df[
  ,
  c("GeneID", setdiff(names(default_df), "GeneID")),
  drop = FALSE
]

no_filter_df <- as.data.frame(
  res_no_independent_filtering
)

no_filter_df$GeneID <- rownames(no_filter_df)

no_filter_df <- no_filter_df[
  ,
  c("GeneID", setdiff(names(no_filter_df), "GeneID")),
  drop = FALSE
]

write.csv(
  default_df,
  file.path(
    output_dir,
    "study1_sensitivity_default_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  no_filter_df,
  file.path(
    output_dir,
    "study1_sensitivity_no_independent_filtering.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 11. Compare results gene by gene
# ------------------------------------------------------------------

comparison_df <- data.frame(
  GeneID = rownames(res_default),
  baseMean = res_default$baseMean,
  log2FoldChange_default = res_default$log2FoldChange,
  pvalue_default = res_default$pvalue,
  padj_default = res_default$padj,
  pvalue_no_independent_filtering =
    res_no_independent_filtering$pvalue,
  padj_no_independent_filtering =
    res_no_independent_filtering$padj
)

comparison_df$significant_default <- (
  !is.na(comparison_df$padj_default) &
    comparison_df$padj_default < 0.05
)

comparison_df$significant_no_independent_filtering <- (
  !is.na(comparison_df$padj_no_independent_filtering) &
    comparison_df$padj_no_independent_filtering < 0.05
)

comparison_df$significance_changed <- (
  comparison_df$significant_default !=
    comparison_df$significant_no_independent_filtering
)

write.csv(
  comparison_df,
  file.path(
    output_dir,
    "study1_sensitivity_comparison.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 12. Summarise the comparison
# ------------------------------------------------------------------

default_valid_padj <- !is.na(res_default$padj)

no_filter_valid_padj <- !is.na(
  res_no_independent_filtering$padj
)

default_significant <- (
  default_valid_padj &
    res_default$padj < 0.05
)

no_filter_significant <- (
  no_filter_valid_padj &
    res_no_independent_filtering$padj < 0.05
)

both_valid <- (
  default_valid_padj &
    no_filter_valid_padj
)

summary_lines <- c(
  "STUDY 1: DESEQ2 SENSITIVITY ANALYSIS",
  "====================================",
  "",
  paste("Analysis date:", as.character(Sys.Date())),
  paste("Design:", deparse(design_formula)),
  "Contrast: Alzheimer's disease versus control",
  "Reference diagnosis: control",
  "Reference sex: female",
  "",
  paste("Original samples:", length(sample_ids)),
  paste("Samples excluded for missing covariates:",
        length(excluded_samples)),
  paste("Samples analysed:", ncol(dds)),
  paste("Control samples:", sum(metadata$diagnosis == "control")),
  paste("Alzheimer's disease samples:",
        sum(metadata$diagnosis == "Alzheimer's disease")),
  "",
  paste("Genes before zero-count filtering:",
        genes_before_filtering),
  paste("Genes after zero-count filtering:",
        genes_after_filtering),
  "",
  "Analysis A: DESeq2 default independent filtering",
  paste("Genes with adjusted p-values:",
        sum(default_valid_padj)),
  paste("Genes with adjusted p-value < 0.05:",
        sum(default_significant)),
  paste("Genes with missing raw p-values:",
        sum(is.na(res_default$pvalue))),
  paste("Genes with missing adjusted p-values:",
        sum(is.na(res_default$padj))),
  "",
  "Analysis B: Independent filtering disabled",
  paste("Genes with adjusted p-values:",
        sum(no_filter_valid_padj)),
  paste("Genes with adjusted p-value < 0.05:",
        sum(no_filter_significant)),
  paste("Genes with missing raw p-values:",
        sum(is.na(res_no_independent_filtering$pvalue))),
  paste("Genes with missing adjusted p-values:",
        sum(is.na(res_no_independent_filtering$padj))),
  "",
  "Comparison",
  paste("Genes with valid adjusted p-values in both analyses:",
        sum(both_valid)),
  paste("Genes significant in both analyses:",
        sum(default_significant & no_filter_significant)),
  paste("Genes significant only with default filtering:",
        sum(default_significant & !no_filter_significant)),
  paste("Genes significant only without independent filtering:",
        sum(!default_significant & no_filter_significant)),
  paste("Genes whose significance classification changed:",
        sum(comparison_df$significance_changed)),
  "",
  "Interpretation note:",
  paste(
    "Independent filtering changes adjusted p-value availability.",
    "Cook's-distance filtering remains enabled in both analyses."
  ),
  paste(
    "Differences in significance counts should be interpreted",
    "alongside adjusted p-values and the independent-filtering",
    "threshold, not as evidence that one analysis is automatically"
  ),
  "correct and the other incorrect."
)

writeLines(
  summary_lines,
  con = file.path(
    output_dir,
    "study1_sensitivity_report.txt"
  )
)

cat(
  paste(summary_lines, collapse = "\n"),
  "\n\nSensitivity analysis completed successfully.\n"
)

cat("Sensitivity results directory: ", output_dir, "\n", sep = "")