```bash
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

# Create output directory.
mkdir -p "${QUANT_DIR}"


# =====================================================================
# 1. INPUT VALIDATION
# =====================================================================

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

# Activate Conda environment.
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate rnaseq

# Check that featureCounts is installed.
if ! command -v featureCounts > /dev/null 2>&1; then
    echo "ERROR: featureCounts is not installed or unavailable."
    exit 1
fi

# Enable nullglob to handle missing BAM patterns.
shopt -s nullglob

# Identify STAR coordinate-sorted BAM files.
BAM_FILES=( "${BAM_DIR}"/*_Aligned.sortedByCoord.out.bam )

# Check that BAM files exist.
if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No STAR sorted BAM files found in ${BAM_DIR}"
    exit 1
fi


# =====================================================================
# 2. BAM VALIDATION
# =====================================================================

# Check that SAMtools is available for BAM integrity checks.
if ! command -v samtools > /dev/null 2>&1; then
    echo "ERROR: SAMtools is not installed or unavailable."
    exit 1
fi

# Validate each BAM file.
for BAM_FILE in "${BAM_FILES[@]}"; do

    if ! samtools quickcheck -v "${BAM_FILE}"; then
        echo "ERROR: Invalid BAM file: ${BAM_FILE}"
        exit 1
    fi

done


# =====================================================================
# 3. FEATURECOUNTS CONFIGURATION
# =====================================================================

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


# =====================================================================
# 4. GENE-LEVEL QUANTIFICATION
# =====================================================================

echo "============================================================"
echo "Starting featureCounts gene-level quantification."
echo "Number of BAM files: ${#BAM_FILES[@]}"
echo "Read type: ${READ_TYPE}"
echo "Strandedness: ${STRANDEDNESS}"
echo "============================================================"

# Run featureCounts.
"${FEATURECOUNTS_CMD[@]}" "${BAM_FILES[@]}"


# =====================================================================
# 5. OUTPUT VALIDATION
# =====================================================================

# Check that the raw featureCounts output exists and is not empty.
if [[ ! -s "${COUNTS_FILE}" ]]; then
    echo "ERROR: featureCounts output is missing or empty."
    exit 1
fi

# Generate a processed count table.
# Remove comment lines and simplify BAM sample column names.

awk '
    BEGIN {
        FS = OFS = "\t"
    }

    /^#/ {
        next
    }

    {
        for (i = 1; i <= NF; i++) {
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


# =====================================================================
# 6. COMPLETION
# =====================================================================

echo "============================================================"
echo "FEATURECOUNTS QUANTIFICATION COMPLETED"
echo "Raw output: ${COUNTS_FILE}"
echo "Processed count table: ${COUNTS_TSV}"
echo "============================================================"