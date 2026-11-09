# ============================================================
# Script 11: Annotate Study 1 DESeq2 Results#
# Purpose:
#   Add gene symbols and descriptions to the existing DESeq2
#   results using the matching NCBI GRCh38.p13 annotation.
#
# Inputs:
#   analysis/study1/study1_deseq2_results.csv
#   raw_data/Human.GRCh38.p13.annot.tsv.gz
#
# Output:
#   analysis/study1/study1_deseq2_results_annotated.csv
#
# The original DESeq2 results are not modified.

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                  # 1. Define file paths
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_file <- "analysis/study1/study1_deseq2_results.csv"

annotation_file <- "raw_data/Human.GRCh38.p13.annot.tsv.gz"

output_file <- "analysis/study1/study1_deseq2_results_annotated.csv"

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                  # 2. Check input files
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------

if (!file.exists(results_file)) {
  stop("DESeq2 results file not found: ", results_file)
}

if (!file.exists(annotation_file)) {
  stop("Annotation file not found: ", annotation_file)
}

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                  # 3. Read DESeq2 results
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------
results <- read.csv(
  results_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!"GeneID" %in% names(results)) {
  stop("The DESeq2 results must contain a GeneID column.")
}

# Keep gene identifiers as character strings.
results$GeneID <- as.character(results$GeneID)

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                 # 4. Read NCBI annotation
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------

annotation <- read.delim(
  gzfile(annotation_file),
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "GeneID",
  "Symbol",
  "Description",
  "GeneType",
  "EnsemblGeneID"
)

missing_columns <- setdiff(
  required_columns,
  names(annotation)
)

if (length(missing_columns) > 0) {
  stop(
    "Required annotation columns are missing: ",
    paste(missing_columns, collapse = ", "),
    "\nCheck the annotation file's column names."
  )
}

annotation$GeneID <- as.character(annotation$GeneID)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                     # 5. Check annotation identifiers
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

if (anyDuplicated(annotation$GeneID)) {
  warning(
    "Duplicate GeneIDs found in the annotation. ",
    "Keeping the first annotation record for each GeneID."
  )

  annotation <- annotation[
    !duplicated(annotation$GeneID),
  ]
}

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                    # 6. Match annotations to DESeq2 results
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Match by NCBI GeneID, preserving the original result order.

annotation_index <- match(
  results$GeneID,
  annotation$GeneID
)

annotation_columns <- c(
  "Symbol",
  "Description",
  "GeneType",
  "EnsemblGeneID"
)

for (column in annotation_columns) {
  results[[column]] <- annotation[[column]][annotation_index]
}

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                   # 7. Calculate annotation coverage
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

total_genes <- nrow(results)

matched_genes <- sum(
  !is.na(annotation_index)
)

unmatched_genes <- total_genes - matched_genes

annotation_percentage <- if (total_genes > 0) {
  100 * matched_genes / total_genes
} else {
  NA_real_
}

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                           # 8. Save annotated results
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

write.csv(
  results,
  output_file,
  row.names = FALSE,
  na = ""
)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                            # 9. Report summary
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

cat("\n")
cat("============================================\n")
cat("Study 1 Gene Annotation Summary\n")
cat("============================================\n")

cat("DESeq2 results loaded: ", results_file, "\n", sep = "")
cat("Annotation loaded: ", annotation_file, "\n", sep = "")
cat("Total result rows: ", total_genes, "\n", sep = "")
cat("Genes matched: ", matched_genes, "\n", sep = "")
cat("Genes unmatched: ", unmatched_genes, "\n", sep = "")

cat(
  "Annotation coverage: ",
  sprintf("%.2f%%", annotation_percentage),
  "\n",
  sep = ""
)

cat("Output saved: ", output_file, "\n", sep = "")

cat("============================================\n\n")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                             # 10. Display the first annotated results
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

cat("First 10 annotated results:\n\n")

display_columns <- c(
  "GeneID",
  "Symbol",
  "Description",
  "log2FoldChange",
  "pvalue",
  "padj"
)

available_display_columns <- intersect(
  display_columns,
  names(results)
)

print(
  head(results[, available_display_columns, drop = FALSE], 10),
  row.names = FALSE
)

cat("\nAnnotation complete.\n")