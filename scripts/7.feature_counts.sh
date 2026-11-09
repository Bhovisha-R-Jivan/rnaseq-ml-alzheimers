#!/bin/bash

# =====================================================================
# SCRIPT 7: GENE-LEVEL QUANTIFICATION WITH FEATURECOUNTS
# =====================================================================

set -euo pipefail

# Load pipeline configuration.
source scripts/0.config.sh

# Define input and output directories.
BAM_DIR="${ANALYSIS_DIR}/4_bamfiles"
QUANT_DIR="${ANALYSIS_DIR}/6_quantification/featureCounts"

# Define output files.
COUNTS_FILE="${QUANT_DIR}/featurecounts.txt"
COUNTS_TSV="${QUANT_DIR}/featurecounts.tsv"
SUMMARY_FILE="${COUNTS_FILE}.summary"

# Create output directory.
mkdir -p "${QUANT_DIR}"


# ---------------------------------------------------------------------
# 1. ACTIVATE CONDA ENVIRONMENT
# ---------------------------------------------------------------------

# Activate Conda environment.
if [[ ! -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
    echo "ERROR: Conda initialization script not found."
    exit 1
fi

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate rnaseq

# Check that featureCounts and SAMtools are installed.
for TOOL in featureCounts samtools; do
    if ! command -v "${TOOL}" > /dev/null 2>&1; then
        echo "ERROR: ${TOOL} is not installed or unavailable."
        exit 1
    fi
done


# ---------------------------------------------------------------------
# 2. INPUT VALIDATION
# ---------------------------------------------------------------------

# Check that the annotation file exists and is not empty.
if [[ ! -s "${ANNOTATION}" ]]; then
    echo "ERROR: Annotation GTF not found or empty: ${ANNOTATION}"
    exit 1
fi

# Check that the BAM directory exists.
if [[ ! -d "${BAM_DIR}" ]]; then
    echo "ERROR: BAM directory not found: ${BAM_DIR}"
    exit 1
fi

# Validate sequencing layout.
if [[ "${READ_TYPE}" != "paired" && "${READ_TYPE}" != "single" ]]; then
    echo "ERROR: Unsupported READ_TYPE: ${READ_TYPE}"
    exit 1
fi

# Enable nullglob to handle missing BAM patterns.
shopt -s nullglob

# Identify STAR coordinate-sorted BAM files.
BAM_FILES=("${BAM_DIR}"/*_Aligned.sortedByCoord.out.bam)

# Check that BAM files exist.
if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No STAR sorted BAM files found in ${BAM_DIR}"
    exit 1
fi


# ---------------------------------------------------------------------
# 3. BAM VALIDATION
# ---------------------------------------------------------------------

echo "Validating BAM files..."

# Validate each BAM file.
for BAM_FILE in "${BAM_FILES[@]}"; do

    if ! samtools quickcheck -v "${BAM_FILE}"; then
        echo "ERROR: Invalid BAM file: ${BAM_FILE}"
        exit 1
    fi

done

echo "BAM validation completed."


# ---------------------------------------------------------------------
# 4. FEATURECOUNTS CONFIGURATION
# ---------------------------------------------------------------------

# Define featureCounts arguments.
FEATURECOUNTS_CMD=(
    featureCounts
    -T "${THREADS}"
    -t "${FEATURE}"
    -g "${ATTRIBUTE}"
    -s "${STRANDEDNESS}"
    -a "${ANNOTATION}"
    -o "${COUNTS_FILE}"
)

# Configure paired-end fragment counting.
if [[ "${READ_TYPE}" == "paired" ]]; then

    FEATURECOUNTS_CMD+=(
        -p
        --countReadPairs
    )

fi


# ---------------------------------------------------------------------
# 5. GENE-LEVEL QUANTIFICATION
# ---------------------------------------------------------------------

echo "============================================================"
echo "Starting featureCounts gene-level quantification."
echo "Number of BAM files: ${#BAM_FILES[@]}"
echo "Read type: ${READ_TYPE}"
echo "Feature: ${FEATURE}"
echo "Attribute: ${ATTRIBUTE}"
echo "Strandedness: ${STRANDEDNESS}"
echo "============================================================"

# Run featureCounts.
"${FEATURECOUNTS_CMD[@]}" "${BAM_FILES[@]}"


# ---------------------------------------------------------------------
# 6. OUTPUT VALIDATION
# ---------------------------------------------------------------------

# Check that the raw featureCounts output exists and is not empty.
if [[ ! -s "${COUNTS_FILE}" ]]; then
    echo "ERROR: featureCounts output is missing or empty."
    exit 1
fi

# Check that the featureCounts summary file exists.
if [[ ! -s "${SUMMARY_FILE}" ]]; then
    echo "ERROR: featureCounts summary file is missing or empty."
    exit 1
fi

echo "Raw featureCounts output verified."
echo "FeatureCounts summary verified."


# ---------------------------------------------------------------------
# 7. GENERATE PROCESSED COUNT TABLE
# ---------------------------------------------------------------------

# Remove comment lines and simplify BAM sample column names.
# Preserve the original gene annotation columns.

awk '
    BEGIN {
        FS = OFS = "\t"
    }

    /^#/ {
        next
    }

    {
        for (i = 7; i <= NF; i++) {
            sub(/^.*\//, "", $i)
            sub(/_Aligned\.sortedByCoord\.out\.bam$/, "", $i)
        }

        print
    }
' "${COUNTS_FILE}" > "${COUNTS_TSV}"


# Check that the processed count table exists and is not empty.
if [[ ! -s "${COUNTS_TSV}" ]]; then
    echo "ERROR: Processed count table is missing or empty."
    exit 1
fi


# ---------------------------------------------------------------------
# 8. COMPLETION
# ---------------------------------------------------------------------

shopt -u nullglob

echo "============================================================"
echo "FEATURECOUNTS QUANTIFICATION COMPLETED"
echo "Number of BAM files: ${#BAM_FILES[@]}"
echo "Raw output: ${COUNTS_FILE}"
echo "Processed count table: ${COUNTS_TSV}"
echo "FeatureCounts summary: ${SUMMARY_FILE}"
echo "============================================================"