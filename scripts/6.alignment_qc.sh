#!/bin/bash

# =====================================================================
# SCRIPT 6: ALIGNMENT QUALITY ASSESSMENT WITH QUALIMAP
# =====================================================================

set -euo pipefail

# Load pipeline configuration.
source scripts/0.config.sh

# Define input and output directories.
BAM_DIR="${ANALYSIS_DIR}/4_bamfiles"
QC_DIR="${ANALYSIS_DIR}/5_alignmentQC"

# Create output directory.
mkdir -p "${QC_DIR}"


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

# Check that Qualimap and SAMtools are installed.
for TOOL in qualimap samtools; do
    if ! command -v "${TOOL}" > /dev/null 2>&1; then
        echo "ERROR: ${TOOL} is not installed or unavailable."
        exit 1
    fi
done


# ---------------------------------------------------------------------
# 2. INPUT VALIDATION
# ---------------------------------------------------------------------

# Check that the BAM directory exists.
if [[ ! -d "${BAM_DIR}" ]]; then
    echo "ERROR: BAM directory not found: ${BAM_DIR}"
    exit 1
fi

# Check that the annotation file exists and is not empty.
if [[ ! -s "${ANNOTATION}" ]]; then
    echo "ERROR: Annotation GTF not found or empty: ${ANNOTATION}"
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

echo "BAM files detected: ${#BAM_FILES[@]}"


# ---------------------------------------------------------------------
# 3. ALIGNMENT QUALITY ASSESSMENT
# ---------------------------------------------------------------------

for BAM_FILE in "${BAM_FILES[@]}"; do

    # Extract sample name from BAM filename.
    BAM_NAME="$(basename "${BAM_FILE}" _Aligned.sortedByCoord.out.bam)"

    # Define sample-specific output directory.
    SAMPLE_QC_DIR="${QC_DIR}/${BAM_NAME}"

    echo "============================================================"
    echo "Running Qualimap for sample: ${BAM_NAME}"
    echo "============================================================"

    # Check BAM integrity.
    if ! samtools quickcheck -v "${BAM_FILE}"; then
        echo "ERROR: BAM file failed integrity check: ${BAM_FILE}"
        exit 1
    fi

    # Check that BAM index exists.
    if [[ ! -s "${BAM_FILE}.bai" && ! -s "${BAM_FILE%.bam}.bai" ]]; then
        echo "ERROR: BAM index not found for ${BAM_NAME}"
        echo "Run SAMtools indexing before alignment QC."
        exit 1
    fi

    # Create sample-specific output directory.
    mkdir -p "${SAMPLE_QC_DIR}"

    # Build Qualimap command.
    QUALIMAP_CMD=(
        qualimap rnaseq
        -bam "${BAM_FILE}"
        -gtf "${ANNOTATION}"
        -outdir "${SAMPLE_QC_DIR}"
        -outfile "${BAM_NAME}"
        -outformat HTML
    )

    # Specify paired-end sequencing when applicable.
    if [[ "${READ_TYPE}" == "paired" ]]; then
        QUALIMAP_CMD+=( -pe )
    fi

    # Run Qualimap.
    "${QUALIMAP_CMD[@]}"

    echo "Completed alignment QC for sample: ${BAM_NAME}"

done

shopt -u nullglob


# ---------------------------------------------------------------------
# 4. COMPLETION
# ---------------------------------------------------------------------

echo "============================================================"
echo "ALIGNMENT QUALITY ASSESSMENT COMPLETED"
echo "Samples assessed: ${#BAM_FILES[@]}"
echo "Qualimap reports: ${QC_DIR}"
echo "============================================================"