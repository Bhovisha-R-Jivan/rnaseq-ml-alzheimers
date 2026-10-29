#!/usr/bin/env Rscript

# Study 1: Alzheimer's disease RNA-seq
# Purpose:
#   Perform differential gene expression analysis using DESeq2.
#   Generate exploratory and differential-expression plots.
#
# Candidate design:
#   ~ AGE + sex + diagnosis
#
# Reference levels:
#   Diagnosis: control
#   Sex: female

# -------------------------------------------------------------------------------------------------------------------------------------------------
#                                                      1. Project paths
# -------------------------------------------------------------------------------------------------------------------------------------------------

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
  "study1"
)

plot_dir <- file.path(
  output_dir,
  "plots"
)

dir.create(
  plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# -------------------------------------------------------------------------------------------------------------------------------------------------
#                                                      2. Check required inputs and packages
# -------------------------------------------------------------------------------------------------------------------------------------------------

if (!file.exists(count_file)) {
  stop("Count matrix not found: ", count_file)
}

if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file)
}

required_packages <- c(
  "DESeq2",
  "ggplot2",
  "pheatmap"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Missing R packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script."
  )
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

# -------------------------------------------------------------------------------------------------------------------------------------------------
#                                                      3. Load count matrix and metadata
# -------------------------------------------------------------------------------------------------------------------------------------------------

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

if (anyNA(counts) || any(!is.finite(counts))) {
  stop("The count matrix contains missing or non-finite values.")
}

if (any(counts < 0) || any(counts != floor(counts))) {
  stop("DESeq2 requires non-negative integer counts.")
}

storage.mode(counts) <- "integer"

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

if (anyDuplicated(colnames(counts))) {
  stop("Count matrix sample IDs are duplicated.")
}

if (!setequal(colnames(counts), sample_ids)) {
  stop(
    "Count matrix sample IDs and metadata sample IDs do not match."
  )
}

# Reorder metadata to match count matrix column order.
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

# -------------------------------------------------------------------------------------------------------------------------------------------------
#                                                      4. Clean and validate covariates
# -------------------------------------------------------------------------------------------------------------------------------------------------

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

# -------------------------------------------------------------------------------------------------------------------------------------------------
#                                                      5. Select complete cases for the adjusted model
# -------------------------------------------------------------------------------------------------------------------------------------------------

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

if (!identical(colnames(counts), metadata[["GEO_Accession (exp)"]])) {
  stop("Sample alignment failed after complete-case filtering.")
}

if (nlevels(droplevels(metadata$diagnosis)) != 2) {
  stop("Both diagnosis groups must be represented.")
}

if (nlevels(droplevels(metadata$sex)) != 2) {
  stop("Both sex categories must be represented.")
}

metadata$diagnosis <- droplevels(metadata$diagnosis)
metadata$sex <- droplevels(metadata$sex)

# -------------------------------------------------------------------------------------------------------------------------------------------------
#                                                      6. Filter genes with zero counts across all samples
# -------------------------------------------------------------------------------------------------------------------------------------------------

genes_before_filtering <- nrow(counts)

keep_genes <- rowSums(counts) > 0

counts <- counts[
  keep_genes,
  ,
  drop = FALSE
]

genes_after_filtering <- nrow(counts)

if (genes_after_filtering == 0) {
  stop("No genes remain after filtering.")
}

# -------------------------------------------------------------------------------------------------------------------------------------------------
#                                                      7. Validate the design matrix
# -------------------------------------------------------------------------------------------------------------------------------------------------

design_formula <- ~ AGE + sex + diagnosis

design_matrix <- model.matrix(
  design_formula,
  data = metadata
)

if (qr(design_matrix)$rank != ncol(design_matrix)) {
  stop("The design matrix is not full rank.")
}

# -------------------------------------------------------------------------------------------------------------------------------------------------
#                                                      8. Construct and fit the DESeq2 model
# -------------------------------------------------------------------------------------------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = design_formula
)

dds <- DESeq(dds)

# ============================================================
# DIAGNOSTIC: Investigate missing DESeq2 p-values
# ============================================================

# Original results using DESeq2 defaults
res_default <- results(
  dds,
  contrast = c("diagnosis", "Alzheimer's disease", "control"),
  alpha = 0.05
)

# Results without Cook's-distance filtering or independent filtering
res_unfiltered <- results(
  dds,
  contrast = c("diagnosis", "Alzheimer's disease", "control"),
  alpha = 0.05,
  cooksCutoff = FALSE,
  independentFiltering = FALSE
)

# Identify genes with missing raw p-values in the original results
missing_pvalue <- is.na(res_default$pvalue)

diagnostic_results <- data.frame(
  GeneID = rownames(res_default),
  baseMean = res_default$baseMean,
  log2FoldChange = res_default$log2FoldChange,
  stat = res_default$stat,
  original_pvalue = res_default$pvalue,
  unfiltered_pvalue = res_unfiltered$pvalue,
  original_padj = res_default$padj,
  maxCooks = mcols(dds)$maxCooks,
  missing_original_pvalue = missing_pvalue
)

# Keep only genes with missing original raw p-values
missing_pvalue_genes <- diagnostic_results[
  diagnostic_results$missing_original_pvalue,
]

# Save the diagnostic separately from the original results
write.csv(
  missing_pvalue_genes,
  "analysis/study1/study1_missing_pvalues_diagnostic.csv",
  row.names = FALSE
)

# Summarise what happened
cat("\n--- Missing p-value diagnostic ---\n")
cat("Genes tested:", nrow(res_default), "\n")
cat("Missing original p-values:", sum(is.na(res_default$pvalue)), "\n")
cat(
  "Missing p-values restored without Cook's filtering:",
  sum(
    is.na(res_default$pvalue) &
    !is.na(res_unfiltered$pvalue)
  ),
  "\n"
)
cat(
  "Missing adjusted p-values:",
  sum(is.na(res_default$padj)),
  "\n"
)

print(missing_pvalue_genes)

# -------------------------------------------------------------------------------------------------------------------------------------------------
#                                                      9. Extract Alzheimer's disease versus control results
# -------------------------------------------------------------------------------------------------------------------------------------------------

results_table <- results(
  dds,
  contrast = c(
    "diagnosis",
    "Alzheimer's disease",
    "control"
  ),
  alpha = 0.05
)

results_df <- as.data.frame(results_table)

results_df$GeneID <- rownames(results_df)

results_df <- results_df[
  ,
  c("GeneID", setdiff(names(results_df), "GeneID")),
  drop = FALSE
]

results_df <- results_df[
  order(results_df$padj, na.last = TRUE),
  ,
  drop = FALSE
]

write.csv(
  results_df,
  file = file.path(
    output_dir,
    "study1_deseq2_results.csv"
  ),
  row.names = FALSE
)

# -------------------------------------------------------------------------------------------------------------------------------------------------
#                                                      10. Save normalised counts
# -------------------------------------------------------------------------------------------------------------------------------------------------

normalised_counts <- counts(dds, normalized = TRUE)

normalised_df <- data.frame(
  GeneID = rownames(normalised_counts),
  normalised_counts,
  check.names = FALSE
)

write.csv(
  normalised_df,
  file = file.path(
    output_dir,
    "study1_normalised_counts.csv"
  ),
  row.names = FALSE
)

# -------------------------------------------------------------------------------------------------------------------------------------------------
#                                                      11. Variance-stabilising transformation
# -------------------------------------------------------------------------------------------------------------------------------------------------

vsd <- vst(
  dds,
  blind = FALSE
)

vst_matrix <- assay(vsd)

# -------------------------------------------------------------------------------------------------------------------------------------------------
#                                                      12. PCA plot
# -------------------------------------------------------------------------------------------------------------------------------------------------

pca_data <- plotPCA(
  vsd,
  intgroup = c("diagnosis"),
  returnData = TRUE
)

percent_variance <- round(
  100 * attr(pca_data, "percentVar")
)

pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    colour = diagnosis
  )
) +
  geom_point(
    size = 2.5,
    alpha = 0.8
  ) +
  labs(
    title = "Study 1: PCA of Alzheimer's disease and control samples",
    x = paste0("PC1: ", percent_variance[1], "% variance"),
    y = paste0("PC2: ", percent_variance[2], "% variance"),
    colour = "Diagnosis"
  ) +
  theme_classic()

ggsave(
  filename = file.path(plot_dir, "study1_pca.png"),
  plot = pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# -------------------------------------------------------------------------------------------------------------------------------------------------
#                                                      13. Sample correlation heatmap
# -------------------------------------------------------------------------------------------------------------------------------------------------
sample_correlations <- cor(
  vst_matrix,
  method = "pearson"
)

annotation_col <- data.frame(
  Diagnosis = metadata$diagnosis,
  Sex = metadata$sex,
  row.names = colnames(vst_matrix)
)

png(
  filename = file.path(
    plot_dir,
    "study1_sample_correlation_heatmap.png"
  ),
  width = 2400,
  height = 2200,
  res = 200
)

pheatmap(
  sample_correlations,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Study 1: Sample correlation heatmap"
)

dev.off()

# -------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                  14. MA plot
# -------------------------------------------------------------------------------------------------------------------------------------------------

png(
  filename = file.path(plot_dir, "study1_MA_plot.png"),
  width = 1800,
  height = 1500,
  res = 200
)

plotMA(
  results_table,
  ylim = c(-5, 5),
  alpha = 0.05,
  main = "Study 1: Alzheimer's disease versus control"
)

dev.off()

# -------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                  15. Volcano plot
# -------------------------------------------------------------------------------------------------------------------------------------------------

volcano_df <- results_df

volcano_df$neg_log10_padj <- -log10(
  pmax(volcano_df$padj, .Machine$double.xmin)
)

volcano_df$significance <- "Not significant"

volcano_df$significance[
  !is.na(volcano_df$padj) &
    volcano_df$padj < 0.05 &
    volcano_df$log2FoldChange > 0
] <- "Higher in Alzheimer's disease"

volcano_df$significance[
  !is.na(volcano_df$padj) &
    volcano_df$padj < 0.05 &
    volcano_df$log2FoldChange < 0
] <- "Lower in Alzheimer's disease"

volcano_df$significance[
  is.na(volcano_df$padj) |
    is.na(volcano_df$log2FoldChange)
] <- NA_character_

volcano_plot <- ggplot(
  volcano_df,
  aes(
    x = log2FoldChange,
    y = neg_log10_padj,
    colour = significance
  )
) +
  geom_point(
    alpha = 0.65,
    size = 1.2,
    na.rm = TRUE
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted"
  ) +
  labs(
    title = "Study 1: Differential expression",
    subtitle = "Alzheimer's disease versus control",
    x = "Log2 fold change",
    y = "-Log10 adjusted p-value",
    colour = "Result"
  ) +
  theme_classic()

ggsave(
  filename = file.path(plot_dir, "study1_volcano_plot.png"),
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 16. Summarise results
# ------------------------------------------------------------

valid_padj <- !is.na(results_df$padj)

significant <- (
  valid_padj &
    results_df$padj < 0.05
)

summary_lines <- c(
  "STUDY 1 DESEQ2 ANALYSIS REPORT",
  "==============================",
  "",
  paste("Analysis date:", as.character(Sys.Date())),
  paste("Candidate design:", deparse(design_formula)),
  "Contrast: Alzheimer's disease versus control",
  "Reference diagnosis: control",
  "Reference sex: female",
  "",
  paste("Samples in original count matrix:", length(sample_ids)),
  paste("Samples excluded for missing covariates:", length(excluded_samples)),
  paste("Samples analysed:", ncol(dds)),
  paste("Control samples analysed:", sum(metadata$diagnosis == "control")),
  paste("Alzheimer's samples analysed:", sum(metadata$diagnosis == "Alzheimer's disease")),
  paste("Genes before zero-count filtering:", genes_before_filtering),
  paste("Genes after zero-count filtering:", genes_after_filtering),
  paste("Genes with adjusted p-values:", sum(valid_padj)),
  paste("Genes with adjusted p-value < 0.05:", sum(significant)),
  "",
  "Excluded sample IDs:",
  if (length(excluded_samples) > 0) {
    excluded_samples
  } else {
    "None"
  }
)

writeLines(
  summary_lines,
  con = file.path(
    output_dir,
    "study1_deseq2_report.txt"
  )
)

cat(
  paste(summary_lines, collapse = "\n"),
  "\n\nDESeq2 analysis completed successfully.\n"
)

cat("Results directory: ", output_dir, "\n", sep = "")
cat("Plot directory: ", plot_dir, "\n", sep = "")