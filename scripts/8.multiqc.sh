#!/bin/bash

# =====================================================================
# SCRIPT 8: FINAL RNA-SEQ QUALITY CONTROL SUMMARY WITH MULTIQC
# =====================================================================

set -euo pipefail

# Load pipeline configuration.
source scripts/0.config.sh

# Define input and output directories.
FASTQC_DIR="${ANALYSIS_DIR}/1_fastqc"
TRIMMED_DIR="${ANALYSIS_DIR}/2_trimmed"
POSTQC_DIR="${ANALYSIS_DIR}/3_postqc_fastqc"
BAM_DIR="${ANALYSIS_DIR}/4_bamfiles"
ALIGNMENT_QC_DIR="${ANALYSIS_DIR}/5_alignmentQC"
COUNTS_DIR="${ANALYSIS_DIR}/6_quantification/featureCounts"

# Define final MultiQC output directory.
MULTIQC_DIR="${ANALYSIS_DIR}/7_multiqc"

# Create output directory.
mkdir -p "${MULTIQC_DIR}"


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

# Check that MultiQC is installed.
if ! command -v multiqc > /dev/null 2>&1; then
    echo "ERROR: MultiQC is not installed or unavailable."
    exit 1
fi


# ---------------------------------------------------------------------
# 2. INPUT VALIDATION
# ---------------------------------------------------------------------

# Check that required analysis directories exist.
for DIR in \
    "${FASTQC_DIR}" \
    "${TRIMMED_DIR}" \
    "${POSTQC_DIR}" \
    "${BAM_DIR}" \
    "${ALIGNMENT_QC_DIR}" \
    "${COUNTS_DIR}"; do

    if [[ ! -d "${DIR}" ]]; then
        echo "ERROR: Required analysis directory not found: ${DIR}"
        exit 1
    fi

done

# Check that raw FastQC reports exist.
if ! compgen -G "${FASTQC_DIR}/*_fastqc.zip" > /dev/null; then
    echo "ERROR: Raw FastQC reports not found."
    exit 1
fi

# Check that Trim Galore reports exist.
if ! compgen -G "${TRIMMED_DIR}/*_trimming_report.txt" > /dev/null; then
    echo "ERROR: Trim Galore reports not found."
    exit 1
fi

# Check that post-trimming FastQC reports exist.
if ! compgen -G "${POSTQC_DIR}/*_fastqc.zip" > /dev/null; then
    echo "ERROR: Post-trimming FastQC reports not found."
    exit 1
fi

# Check that STAR alignment logs exist.
if ! compgen -G "${BAM_DIR}/*_Log.final.out" > /dev/null; then
    echo "ERROR: STAR alignment logs not found."
    exit 1
fi

# Check that featureCounts summary exists.
if [[ ! -s "${COUNTS_DIR}/featurecounts.txt.summary" ]]; then
    echo "ERROR: featureCounts summary file not found or empty."
    exit 1
fi


# ---------------------------------------------------------------------
# 3. DEFINE MULTIQC INPUT DIRECTORIES
# ---------------------------------------------------------------------

# MultiQC searches these directories recursively for supported reports.
MULTIQC_INPUTS=(
    "${FASTQC_DIR}"
    "${TRIMMED_DIR}"
    "${POSTQC_DIR}"
    "${BAM_DIR}"
    "${ALIGNMENT_QC_DIR}"
    "${COUNTS_DIR}"
)


# ---------------------------------------------------------------------
# 4. GENERATE FINAL MULTIQC REPORT
# ---------------------------------------------------------------------

echo "============================================================"
echo "Generating final RNA-seq MultiQC report."
echo "============================================================"

multiqc \
    "${MULTIQC_INPUTS[@]}" \
    --outdir "${MULTIQC_DIR}" \
    --filename multiqc_report.html \
    --title "RNA-seq Pipeline Quality Control Summary" \
    --force


# ---------------------------------------------------------------------
# 5. OUTPUT VALIDATION
# ---------------------------------------------------------------------

# Check that the final MultiQC HTML report exists.
if [[ ! -s "${MULTIQC_DIR}/multiqc_report.html" ]]; then
    echo "ERROR: Final MultiQC HTML report was not generated."
    exit 1
fi

echo "Final MultiQC report generated successfully."


# ---------------------------------------------------------------------
# 6. COMPLETION
# ---------------------------------------------------------------------

echo "============================================================"
echo "FINAL MULTIQC REPORT COMPLETED"
echo "Report: ${MULTIQC_DIR}/multiqc_report.html"
echo "Output directory: ${MULTIQC_DIR}"
echo "============================================================"