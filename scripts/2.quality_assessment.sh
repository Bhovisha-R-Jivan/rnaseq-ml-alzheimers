#!/bin/bash

# =====================================================================
# SCRIPT 2: FASTQ QUALITY ASSESSMENT
# =====================================================================

set -euo pipefail

# Load pipeline configuration.
source scripts/0.config.sh

# Define FastQC output directory.
QC_DIR="${ANALYSIS_DIR}/1_fastqc"

# Create output directory.
mkdir -p "${QC_DIR}"

# Activate Conda environment.
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate rnaseq

# Check that FastQC is installed.
if ! command -v fastqc > /dev/null 2>&1; then
    echo "ERROR: FastQC is not installed or unavailable."
    exit 1
fi

# Check that raw data directory exists.
if [[ ! -d "${RAW_DATA_DIR}" ]]; then
    echo "ERROR: Raw data directory not found: ${RAW_DATA_DIR}"
    exit 1
fi

# Enable nullglob so unmatched patterns produce empty arrays.
shopt -s nullglob

# =====================================================================
# 1. PAIRED-END DATA
# =====================================================================

if [[ "${READ_TYPE}" == "paired" ]]; then

    echo "Processing paired-end FASTQ files."

    # Identify forward and reverse reads.
    R1_FILES=( "${RAW_DATA_DIR}"/*"${R1_SUFFIX}" )
    R2_FILES=( "${RAW_DATA_DIR}"/*"${R2_SUFFIX}" )

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
    fastqc \
        --outdir "${QC_DIR}" \
        --threads "${THREADS}" \
        "${R1_FILES[@]}" \
        "${R2_FILES[@]}"


# =====================================================================
# 2. SINGLE-END DATA
# =====================================================================

elif [[ "${READ_TYPE}" == "single" ]]; then

    echo "Processing single-end FASTQ files."

    # Identify FASTQ files.
    FASTQ_FILES=( "${RAW_DATA_DIR}"/*.fastq.gz )

    # Check that FASTQ files exist.
    if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
        echo "ERROR: No FASTQ files found in ${RAW_DATA_DIR}"
        exit 1
    fi

    # Run FastQC on all single-end FASTQ files.
    fastqc \
        --outdir "${QC_DIR}" \
        --threads "${THREADS}" \
        "${FASTQ_FILES[@]}"

else

    echo "ERROR: Unsupported READ_TYPE: ${READ_TYPE}"
    exit 1

fi


# =====================================================================
# 3. COMPLETION
# =====================================================================

echo "============================================================"
echo "FASTQC QUALITY ASSESSMENT COMPLETED"
echo "Output directory: ${QC_DIR}"
echo "============================================================"