#!/usr/bin/env Rscript

# Study 1: Alzheimer's disease RNA-seq
# Purpose: Assess metadata and candidate covariates before DESeq2.

project_root <- getwd()

metadata_file <- file.path(
  project_root,
  "raw_data",
  "Study1_MetaData.csv"
)

if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file)
}

metadata <- read.csv(
  metadata_file,
  sep = ";",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

# Convert known placeholders to missing values for analysis.
metadata$AGE[
  tolower(trimws(metadata$AGE)) %in%
    c("not applicable", "not available", "unknown", "n/a")
] <- NA

metadata$AGE <- as.numeric(metadata$AGE)

# Report age and sex distributions by diagnosis.
cat("\n--- Sample counts by diagnosis ---\n")
print(table(metadata$diagnosis, useNA = "ifany"))

cat("\n--- Age summary by diagnosis ---\n")
print(
  aggregate(
    AGE ~ diagnosis,
    data = metadata,
    FUN = function(x) {
      c(
        n = sum(!is.na(x)),
        mean = mean(x, na.rm = TRUE),
        median = median(x, na.rm = TRUE),
        sd = sd(x, na.rm = TRUE)
      )
    }
  )
)

cat("\n--- Sex by diagnosis ---\n")
print(
  table(
    metadata$diagnosis,
    metadata$sex,
    useNA = "ifany"
  )
)

cat("\n--- Missing age and sex ---\n")
print(
  metadata[
    is.na(metadata$AGE) | is.na(metadata$sex),
    c("GEO_Accession (exp)", "diagnosis", "AGE", "sex")
  ]
)
