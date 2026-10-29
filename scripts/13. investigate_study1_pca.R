# usr/bin/env Rscript

# Study 1: Investigate PCA sample structure
# Purpose:
#   Assess associations between PCA coordinates and available
#   metadata variables: age, sex, and diagnosis.
#
# This script does not modify Script 10 or its outputs.

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
  "pca_investigation"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------------
# 2. Check packages and inputs
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

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 is not installed in the active R environment.")
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

# ------------------------------------------------------------------
# 3. Load count matrix and metadata
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

if (
  anyNA(sample_ids) ||
  anyDuplicated(sample_ids) ||
  anyDuplicated(colnames(counts)) ||
  !setequal(colnames(counts), sample_ids)
) {
  stop("Count matrix and metadata sample IDs do not match.")
}

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
# 4. Clean covariates using Script 10 conventions
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

if (
  !all(
    unique(metadata$sex[!is.na(metadata$sex)]) %in%
      c("female", "male")
  )
) {
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
# 5. Apply the same complete-case selection as Script 10
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

# ------------------------------------------------------------------
# 6. Apply the same zero-count gene filter
# ------------------------------------------------------------------

keep_genes <- rowSums(counts) > 0

counts <- counts[
  keep_genes,
  ,
  drop = FALSE
]

# ------------------------------------------------------------------
# 7. Transform counts and calculate PCA
# ------------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ AGE + sex + diagnosis
)

# Use the same variance-stabilising approach as Script 10.
vsd <- vst(
  dds,
  blind = FALSE
)

pca <- prcomp(
  t(assay(vsd)),
  center = TRUE,
  scale. = FALSE
)

variance_explained <- (
  pca$sdev^2 / sum(pca$sdev^2)
) * 100

pca_df <- data.frame(
  SampleID = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  AGE = metadata$AGE,
  sex = metadata$sex,
  diagnosis = metadata$diagnosis
)

write.csv(
  pca_df,
  file.path(output_dir, "study1_pca_coordinates.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 8. Test associations with PC1 and PC2
# ------------------------------------------------------------------

# Age: Pearson correlation with each principal component.
age_pc1 <- cor.test(
  pca_df$AGE,
  pca_df$PC1,
  method = "pearson"
)

age_pc2 <- cor.test(
  pca_df$AGE,
  pca_df$PC2,
  method = "pearson"
)

# Sex and diagnosis: one-way ANOVA for each component.
sex_pc1 <- summary(
  aov(PC1 ~ sex, data = pca_df)
)[[1]]

sex_pc2 <- summary(
  aov(PC2 ~ sex, data = pca_df)
)[[1]]

diagnosis_pc1 <- summary(
  aov(PC1 ~ diagnosis, data = pca_df)
)[[1]]

diagnosis_pc2 <- summary(
  aov(PC2 ~ diagnosis, data = pca_df)
)[[1]]

association_results <- data.frame(
  Variable = c(
    "AGE",
    "AGE",
    "sex",
    "sex",
    "diagnosis",
    "diagnosis"
  ),
  Component = c(
    "PC1",
    "PC2",
    "PC1",
    "PC2",
    "PC1",
    "PC2"
  ),
  Statistic = c(
    unname(age_pc1$estimate),
    unname(age_pc2$estimate),
    sex_pc1$`F value`[1],
    sex_pc2$`F value`[1],
    diagnosis_pc1$`F value`[1],
    diagnosis_pc2$`F value`[1]
  ),
  P_value = c(
    age_pc1$p.value,
    age_pc2$p.value,
    sex_pc1$`Pr(>F)`[1],
    sex_pc2$`Pr(>F)`[1],
    diagnosis_pc1$`Pr(>F)`[1],
    diagnosis_pc2$`Pr(>F)`[1]
  )
)

# Correct for the six exploratory association tests.
association_results$Adjusted_P_value <- p.adjust(
  association_results$P_value,
  method = "BH"
)

write.csv(
  association_results,
  file.path(output_dir, "study1_pca_associations.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 9. Generate PCA plots coloured by metadata
# ------------------------------------------------------------------

plot_pca <- function(data, colour_variable, colour_label, filename) {

  p <- ggplot(
    data,
    aes(
      x = PC1,
      y = PC2,
      colour = .data[[colour_variable]]
    )
  ) +
    geom_point(
      size = 2.5,
      alpha = 0.8
    ) +
    labs(
      title = paste(
        "Study 1: PCA coloured by",
        colour_label
      ),
      x = paste0(
        "PC1: ",
        round(variance_explained[1], 1),
        "% variance"
      ),
      y = paste0(
        "PC2: ",
        round(variance_explained[2], 1),
        "% variance"
      ),
      colour = colour_label
    ) +
    theme_classic()

  ggsave(
    filename = file.path(output_dir, filename),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}

plot_pca(
  pca_df,
  "sex",
  "Sex",
  "study1_pca_by_sex.png"
)

plot_pca(
  pca_df,
  "diagnosis",
  "Diagnosis",
  "study1_pca_by_diagnosis.png"
)

# Age is continuous, so use a continuous colour scale.
age_plot <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    colour = AGE
  )
) +
  geom_point(
    size = 2.5,
    alpha = 0.8
  ) +
  labs(
    title = "Study 1: PCA coloured by age",
    x = paste0(
      "PC1: ",
      round(variance_explained[1], 1),
      "% variance"
    ),
    y = paste0(
      "PC2: ",
      round(variance_explained[2], 1),
      "% variance"
    ),
    colour = "Age"
  ) +
  theme_classic()

ggsave(
  filename = file.path(output_dir, "study1_pca_by_age.png"),
  plot = age_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------------
# 10. Write summary report
# ------------------------------------------------------------------

report_lines <- c(
  "STUDY 1: PCA SAMPLE-STRUCTURE INVESTIGATION",
  "==========================================",
  "",
  paste("Analysis date:", as.character(Sys.Date())),
  paste("Samples in original count matrix:", length(sample_ids)),
  paste("Samples excluded for missing covariates:",
        length(excluded_samples)),
  paste("Samples analysed:", nrow(pca_df)),
  paste("PC1 variance explained:",
        round(variance_explained[1], 2), "%"),
  paste("PC2 variance explained:",
        round(variance_explained[2], 2), "%"),
  "",
  "Associations between metadata and PCA coordinates",
  capture.output(print(association_results)),
  "",
  "Interpretation:",
  paste(
    "These exploratory associations do not establish causation",
    "or prove that a variable caused sample clustering."
  ),
  paste(
    "Interpret p-values alongside effect sizes, sample sizes,",
    "the PCA plots, and the study design."
  )
)

writeLines(
  report_lines,
  con = file.path(output_dir, "study1_pca_investigation_report.txt")
)

cat(
  paste(report_lines, collapse = "\n"),
  "\n\nPCA investigation completed successfully.\n"
)

cat("Results directory: ", output_dir, "\n", sep = "")