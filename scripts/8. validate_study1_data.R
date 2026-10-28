#!/usr/bin/env Rscript

# Study 1: Alzheimer's disease RNA-seq
# Purpose: Validate the input count matrix and sample metadata.
#
# Inputs:
#   raw_data/Study1_gene_counts.gz
#   raw_data/Study1_MetaData.csv
#
# Output:
#   analysis/study1/data_validation_report.txt
#
# Run from the repository root:
#   Rscript "scripts/8. validate_study1_data.R"

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                             1. Identify the project directory
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

report_file <- file.path(
  output_dir,
  "data_validation_report.txt"
)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                          2. Check that the input files exist
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------

if (!file.exists(count_file)) {
  stop("Count matrix not found: ", count_file)
}

if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file)
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                          3. Import the count matrix and metadata
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------

message("Importing count matrix...")

counts <- read.delim(
  gzfile(count_file),
  header = TRUE,
  row.names = NULL,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

message("Importing sample metadata...")

metadata <- read.csv(
  metadata_file,
  sep = ";",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------   
#                                                         4. Validate basic structure
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------

if (!"GeneID" %in% colnames(counts)) {
  stop("The count matrix must contain a GeneID column.")
}

gene_ids <- counts$GeneID
sample_ids <- setdiff(colnames(counts), "GeneID")

if (anyDuplicated(gene_ids)) {
  stop("Duplicate gene identifiers detected.")
}

if (anyDuplicated(sample_ids)) {
  stop("Duplicate sample identifiers detected in the count matrix.")
}

if (anyDuplicated(metadata[["GEO_Accession (exp)"]])) {
  stop("Duplicate sample identifiers detected in the metadata.")
}

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                        5. Validate count values
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------

count_values <- as.matrix(
  counts[, sample_ids, drop = FALSE]
)

if (!all(vapply(
  as.data.frame(count_values),
  is.numeric,
  logical(1)
))) {
  stop("One or more sample columns contain non-numeric values.")
}

if (anyNA(count_values)) {
  stop("Missing values detected in the count matrix.")
}

if (any(!is.finite(count_values))) {
  stop("Non-finite values detected in the count matrix.")
}

if (any(count_values < 0)) {
  stop("Negative counts detected.")
}

if (any(count_values != floor(count_values))) {
  stop("Non-integer values detected in the count matrix.")
}

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                        6. Validate sample identifiers
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------

metadata_ids <- metadata[["GEO_Accession (exp)"]]

if (anyNA(metadata_ids) || any(metadata_ids == "")) {
  stop("Missing sample identifiers detected in the metadata.")
}

missing_metadata <- setdiff(sample_ids, metadata_ids)
extra_metadata <- setdiff(metadata_ids, sample_ids)

if (length(missing_metadata) > 0 ||
    length(extra_metadata) > 0) {
  stop(
    "Count matrix and metadata sample IDs do not match.\n",
    "Samples missing from metadata: ",
    paste(missing_metadata, collapse = ", "),
    "\nSamples missing from count matrix: ",
    paste(extra_metadata, collapse = ", ")
  )
}

# Reorder metadata to match the count matrix.
metadata <- metadata[
  match(sample_ids, metadata_ids),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    sample_ids,
    metadata[["GEO_Accession (exp)"]]
  )
)

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                        7. Summarise the dataset
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------

diagnosis_counts <- table(
  metadata$diagnosis,
  useNA = "ifany"
)

report_lines <- c(
  "STUDY 1 DATA VALIDATION REPORT",
  "================================",
  "",
  paste("Project directory:", normalizePath(project_root)),
  paste("Count matrix:", count_file),
  paste("Metadata:", metadata_file),
  paste("Validation date:", Sys.Date()),
  "",
  "COUNT MATRIX",
  "------------",
  paste("Gene rows:", nrow(counts)),
  paste("Sample columns:", length(sample_ids)),
  paste("Unique gene IDs:", length(unique(gene_ids))),
  paste("Unique sample IDs:", length(unique(sample_ids))),
  paste("Minimum count:", min(count_values)),
  paste("Maximum count:", max(count_values)),
  "Count values: non-negative integers",
  "",
  "SAMPLE METADATA",
  "---------------",
  paste("Metadata records:", nrow(metadata)),
  paste("Sample IDs match:", TRUE),
  "",
  "DIAGNOSIS GROUPS",
  "----------------",
  capture.output(print(diagnosis_counts)),
  "",
  "MISSING METADATA VALUES",
  "-----------------------",
  paste(
    capture.output(
      print(colSums(is.na(metadata)))
    ),
    collapse = "\n"
  ),
  "",
  "OVERALL RESULT: PASS",
  "",
  "Note: Passing these checks confirms structural consistency",
  "and count-value formatting. It does not establish that the",
  "experimental design or statistical assumptions are valid."
)

writeLines(report_lines, report_file)

message("Validation completed successfully.")
message("Report saved to: ", report_file)