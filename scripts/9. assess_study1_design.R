#!/usr/bin/env Rscript

# Study 1: Alzheimer's disease RNA-seq
# Purpose: Assess metadata and validate the candidate DESeq2 design.

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                   1. Project paths and input validation
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------

project_root <- getwd()

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

report_file <- file.path(
  output_dir,
  "design_assessment_report.txt"
)

if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file)
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                   2. Load metadata
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

missing_columns <- setdiff(
  required_columns,
  names(metadata)
)

if (length(missing_columns) > 0) {
  stop(
    "Required metadata columns are missing: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (anyDuplicated(metadata[["GEO_Accession (exp)"]])) {
  stop("Duplicate GEO sample identifiers found.")
}

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                   3. Validate and clean age values
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------

age_original <- trimws(as.character(metadata$AGE))

age_placeholders <- c(
  "",
  "na",
  "n/a",
  "not applicable",
  "not available",
  "unknown"
)

age_is_placeholder <- (
  is.na(age_original) |
    tolower(age_original) %in% age_placeholders
)

age_numeric <- suppressWarnings(
  as.numeric(age_original)
)

invalid_age <- (
  !age_is_placeholder &
    (is.na(age_numeric) | !is.finite(age_numeric))
)

if (any(invalid_age)) {
  cat(
    "Invalid age values found:\n",
    file = stderr()
  )

  print(
    data.frame(
      sample_id = metadata[["GEO_Accession (exp)"]][invalid_age],
      original_age = age_original[invalid_age]
    )
  )

  stop(
    "Age validation failed. Resolve invalid values before proceeding."
  )
}

age_numeric[age_is_placeholder] <- NA_real_
metadata$AGE <- age_numeric

# Clean sex placeholders.
metadata$sex <- trimws(as.character(metadata$sex))

sex_placeholders <- c(
  "",
  "na",
  "n/a",
  "not applicable",
  "not available",
  "unknown"
)

sex_is_placeholder <- (
  is.na(metadata$sex) |
    tolower(metadata$sex) %in% sex_placeholders
)

metadata$sex[sex_is_placeholder] <- NA_character_

# Validate diagnosis labels.
expected_diagnoses <- c(
  "Alzheimer's disease",
  "control"
)

unexpected_diagnoses <- setdiff(
  unique(metadata$diagnosis),
  expected_diagnoses
)

if (
  anyNA(metadata$diagnosis) ||
    length(unexpected_diagnoses) > 0
) {
  stop(
    "Missing or unexpected diagnosis labels found. ",
    "Review the metadata before proceeding."
  )
}

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                   4. Prepare the candidate design
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Explicit reference levels:
# - control is the diagnosis reference
# - female is the sex reference

metadata$diagnosis <- factor(
  metadata$diagnosis,
  levels = c("control", "Alzheimer's disease")
)

observed_sexes <- unique(
  metadata$sex[!is.na(metadata$sex)]
)

if (!all(observed_sexes %in% c("female", "male"))) {
  stop(
    "Unexpected sex labels found: ",
    paste(
      setdiff(observed_sexes, c("female", "male")),
      collapse = ", "
    )
  )
}

metadata$sex <- factor(
  metadata$sex,
  levels = c("female", "male")
)

candidate_formula <- ~ AGE + sex + diagnosis

# Complete cases are used only for design assessment here.
# This does not modify the original metadata file.

complete_case <- complete.cases(
  metadata[, c("AGE", "sex", "diagnosis")]
)

complete_metadata <- droplevels(
  metadata[complete_case, , drop = FALSE]
)

if (nrow(complete_metadata) == 0) {
  stop("No complete samples remain for design assessment.")
}

if (nlevels(complete_metadata$diagnosis) < 2) {
  stop(
    "Both diagnosis groups must be represented ",
    "among complete cases."
  )
}

if (nlevels(complete_metadata$sex) < 2) {
  stop(
    "Both sex categories must be represented ",
    "among complete cases."
  )
}

# Build the design matrix and assess estimability.
design_matrix <- model.matrix(
  candidate_formula,
  data = complete_metadata
)

design_rank <- qr(design_matrix)$rank
design_columns <- ncol(design_matrix)

design_full_rank <- design_rank == design_columns

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                   5. Generate the design assessment report
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------

report_lines <- capture.output({

  cat("STUDY 1 DESIGN ASSESSMENT REPORT\n")
  cat("================================\n\n")

  cat("Project directory:", project_root, "\n")
  cat("Metadata file:", metadata_file, "\n")
  cat("Candidate design: ~ AGE + sex + diagnosis\n")
  cat("Assessment date:", as.character(Sys.Date()), "\n\n")

  cat("SAMPLE COUNTS BY DIAGNOSIS\n")
  cat("--------------------------\n")
  print(table(metadata$diagnosis, useNA = "ifany"))

  cat("\nAGE SUMMARY BY DIAGNOSIS\n")
  cat("-----------------------\n")

  age_summary <- aggregate(
    AGE ~ diagnosis,
    data = metadata,
    FUN = function(x) {
      c(
        n = sum(!is.na(x)),
        mean = if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE),
        median = if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE),
        sd = if (sum(!is.na(x)) < 2) NA_real_ else sd(x, na.rm = TRUE),
        min = if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE),
        max = if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
      )
    }
  )

  print(age_summary)

  cat("\nSEX BY DIAGNOSIS\n")
  cat("---------------\n")
  print(
    table(
      metadata$diagnosis,
      metadata$sex,
      useNA = "ifany"
    )
  )

  cat("\nAGE SUMMARY BY DIAGNOSIS AND SEX\n")
  cat("-------------------------------\n")

  age_sex_summary <- aggregate(
    AGE ~ diagnosis + sex,
    data = metadata,
    FUN = function(x) {
      c(
        n = sum(!is.na(x)),
        min = if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE),
        max = if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
      )
    }
  )

  print(age_sex_summary)

  cat("\nMISSING COVARIATES\n")
  cat("------------------\n")

  cat("Samples with missing age:", sum(is.na(metadata$AGE)), "\n")
  cat("Samples with missing sex:", sum(is.na(metadata$sex)), "\n")
  cat("Samples with either missing age or sex:",
      sum(!complete_case), "\n")

  missing_covariates <- metadata[
    !complete_case,
    c("GEO_Accession (exp)", "diagnosis", "AGE", "sex"),
    drop = FALSE
  ]

  if (nrow(missing_covariates) > 0) {
    print(missing_covariates)
  } else {
    cat("No samples with missing candidate covariates.\n")
  }

  cat("\nCOMPLETE-CASE GROUP COUNTS\n")
  cat("--------------------------\n")
  print(table(complete_metadata$diagnosis))

  cat("\nDESIGN MATRIX ASSESSMENT\n")
  cat("-----------------------\n")

  cat("Formula:", deparse(candidate_formula), "\n")
  cat("Reference diagnosis: control\n")
  cat("Reference sex: female\n")
  cat("Complete-case samples:", nrow(complete_metadata), "\n")
  cat("Design matrix rows:", nrow(design_matrix), "\n")
  cat("Design matrix columns:", design_columns, "\n")
  cat("Design matrix rank:", design_rank, "\n")
  cat("Full rank:", design_full_rank, "\n\n")

  cat("Design matrix column names:\n")
  print(colnames(design_matrix))

  cat("\nINTERPRETATION AND LIMITATIONS\n")
  cat("------------------------------\n")

  cat(
    "This assessment evaluates metadata completeness and ",
    "the mathematical rank of the candidate design matrix.\n",
    sep = ""
  )

  cat(
    "A full-rank matrix is necessary but does not, by itself, ",
    "establish that the model is scientifically appropriate.\n",
    sep = ""
  )

  cat(
    "The candidate design and complete-case sample set must ",
    "be reviewed before differential expression analysis.\n",
    sep = ""
  )

  cat(
    "The original metadata has not been modified or saved over.\n"
  )
})

writeLines(
  report_lines,
  con = report_file
)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                     6. Print summary and finish
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------

cat(paste(report_lines, collapse = "\n"), "\n")

cat("\nReport saved to:\n", report_file, "\n", sep = "")

if (!design_full_rank) {
  stop(
    "The candidate design matrix is not full rank. ",
    "Do not proceed to DESeq2 with this design until investigated."
  )
}

cat("\nDesign matrix is full rank.\n")
cat("Script 9 completed successfully.\n")