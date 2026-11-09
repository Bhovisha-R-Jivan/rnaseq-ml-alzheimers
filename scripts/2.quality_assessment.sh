#!/bin/bash

# =====================================================================
# SCRIPT 2: FASTQ QUALITY ASSESSMENT
# =====================================================================

set -euo pipefail

# Load pipeline configuration.
source scripts/0.config.sh

# Define FastQC and MultiQC output directories.
QC_DIR="${ANALYSIS_DIR}/1_fastqc"
MULTIQC_DIR="${QC_DIR}/multiqc_report"

# Create output directories.
mkdir -p "${QC_DIR}"
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

# Check that FastQC is installed.
if ! command -v fastqc > /dev/null 2>&1; then
    echo "ERROR: FastQC is not installed or unavailable."
    exit 1
fi

# Check that MultiQC is installed.
if ! command -v multiqc > /dev/null 2>&1; then
    echo "ERROR: MultiQC is not installed or unavailable."
    exit 1
fi


# ---------------------------------------------------------------------
# 2. CHECK INPUT DATA
# ---------------------------------------------------------------------

# Check that raw data directory exists.
if [[ ! -d "${RAW_DATA_DIR}" ]]; then
    echo "ERROR: Raw data directory not found: ${RAW_DATA_DIR}"
    exit 1
fi

# Enable nullglob so unmatched patterns produce empty arrays.
shopt -s nullglob


# =====================================================================
# 3. PAIRED-END DATA
# =====================================================================

if [[ "${READ_TYPE}" == "paired" ]]; then

    echo "Processing paired-end FASTQ files."

    # Identify forward and reverse reads.
    R1_FILES=("${RAW_DATA_DIR}"/*"${R1_SUFFIX}")
    R2_FILES=("${RAW_DATA_DIR}"/*"${R2_SUFFIX}")

    # Check that forward and reverse reads exist.
    if [[ ${#R1_FILES[@]} -eq 0 || ${#R2_FILES[@]} -eq 0 ]]; then
        echo "ERROR: Paired-end FASTQ files are missing."
        exit 1
    fi

    # Check that the number of forward and reverse files matches.
    if [[ ${#R1_FILES[@]} -ne ${#R2_FILES[@]} ]]; then
        echo "ERROR: Number of forward and reverse FASTQ files does not match."
        exit 1
    fi

    # Verify that every forward read has a corresponding reverse read.
    for R1 in "${R1_FILES[@]}"; do

        SAMPLE="${R1%"${R1_SUFFIX}"}"
        R2="${SAMPLE}${R2_SUFFIX}"

        if [[ ! -f "${R2}" ]]; then
            echo "ERROR: Missing reverse read for ${R1}"
            exit 1
        fi

    done

    # Run FastQC on all paired-end FASTQ files.
    echo "Running FastQC..."

    fastqc \
        --outdir "${QC_DIR}" \
        --threads "${THREADS}" \
        "${R1_FILES[@]}" \
        "${R2_FILES[@]}"


# =====================================================================
# 4. SINGLE-END DATA
# =====================================================================

elif [[ "${READ_TYPE}" == "single" ]]; then

    echo "Processing single-end FASTQ files."

    # Identify FASTQ files.
    FASTQ_FILES=("${RAW_DATA_DIR}"/*.fastq.gz)

    # Check that FASTQ files exist.
    if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
        echo "ERROR: No FASTQ files found in ${RAW_DATA_DIR}"
        exit 1
    fi

    # Run FastQC on all single-end FASTQ files.
    echo "Running FastQC..."

    fastqc \
        --outdir "${QC_DIR}" \
        --threads "${THREADS}" \
        "${FASTQ_FILES[@]}"

else

    echo "ERROR: Unsupported READ_TYPE: ${READ_TYPE}"
    exit 1

fi

shopt -u nullglob


# =====================================================================
# 5. MULTIQC REPORT
# =====================================================================

echo "Generating MultiQC report..."

# Summarise FastQC results.
multiqc \
    "${QC_DIR}" \
    --outdir "${MULTIQC_DIR}" \
    --force

echo "MultiQC report generated."


# =====================================================================
# 6. COMPLETION
# =====================================================================

echo "============================================================"
echo "FASTQ QUALITY ASSESSMENT COMPLETED"
echo "FastQC output: ${QC_DIR}"
echo "MultiQC report: ${MULTIQC_DIR}"
echo "============================================================"