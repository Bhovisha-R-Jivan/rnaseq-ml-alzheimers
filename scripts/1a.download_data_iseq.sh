#!/bin/bash

# =====================================================================
# SCRIPT 1: RNA-SEQ DATA ACQUISITION
# =====================================================================

set -euo pipefail

# Load pipeline configuration.
source scripts/0.config.sh

# Create raw data directory if it does not exist.
mkdir -p "${RAW_DATA_DIR}"


# ---------------------------------------------------------------------
# 1. LOCAL FASTQ INPUT
# ---------------------------------------------------------------------

if [[ "${DATA_SOURCE}" == "local" ]]; then

    echo "Using locally supplied FASTQ files."
    echo "FASTQ directory: ${RAW_DATA_DIR}"

    # Check that at least one gzipped FASTQ file exists.
    if ! compgen -G "${RAW_DATA_DIR}/*.fastq.gz" > /dev/null; then
        echo "ERROR: No .fastq.gz files found in ${RAW_DATA_DIR}"
        exit 1
    fi

    echo "Local FASTQ files detected."
    echo "No download required."


# ---------------------------------------------------------------------
# 2. SRA ACCESSION DOWNLOAD
# ---------------------------------------------------------------------

elif [[ "${DATA_SOURCE}" == "sra" ]]; then

    # Check that the accession file exists and is not empty.
    if [[ ! -f "${ACCESSION_FILE}" ]]; then
        echo "ERROR: Accession file not found: ${ACCESSION_FILE}"
        exit 1
    fi

    if [[ ! -s "${ACCESSION_FILE}" ]]; then
        echo "ERROR: Accession file is empty."
        exit 1
    fi

    # Activate Conda environment.
    if [[ ! -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
        echo "ERROR: Conda initialization script not found."
        exit 1
    fi

    source "$HOME/miniconda3/etc/profile.d/conda.sh"
    conda activate rnaseq

    # Check iSeq availability.
    if ! command -v iseq > /dev/null 2>&1; then
        echo "ERROR: iSeq is not installed or unavailable."
        exit 1
    fi

    # Download each SRA run accession individually.
    while IFS= read -r accession || [[ -n "${accession}" ]]; do

        # Remove Windows carriage returns and surrounding whitespace.
        accession="${accession//$'\r'/}"
        accession="${accession#"${accession%%[![:space:]]*}"}"
        accession="${accession%"${accession##*[![:space:]]}"}"

        # Skip empty lines.
        [[ -z "${accession}" ]] && continue

        echo "Downloading ${accession}..."

        iseq \
            --gzip \
            --fastq \
            --input "${accession}" \
            --database sra \
            --output "${RAW_DATA_DIR}"

    done < "${ACCESSION_FILE}"

    echo "SRA downloads completed."

else

    echo "ERROR: Invalid DATA_SOURCE."
    echo "Allowed values: local or sra."
    exit 1

fi


# ---------------------------------------------------------------------
# 3. FINAL INPUT CHECK
# ---------------------------------------------------------------------

echo "Checking FASTQ input files..."

# Check that gzipped FASTQ files were produced or detected.
if ! compgen -G "${RAW_DATA_DIR}/*.fastq.gz" > /dev/null; then
    echo "ERROR: No .fastq.gz files found after data acquisition."
    exit 1
fi

# Validate paired-end FASTQ files.
if [[ "${READ_TYPE}" == "paired" ]]; then

    shopt -s nullglob

    R1_FILES=("${RAW_DATA_DIR}"/*"${R1_SUFFIX}")
    R2_FILES=("${RAW_DATA_DIR}"/*"${R2_SUFFIX}")

    shopt -u nullglob

    if [[ ${#R1_FILES[@]} -eq 0 || ${#R2_FILES[@]} -eq 0 ]]; then
        echo "ERROR: Paired-end FASTQ files are incomplete."
        echo "R1 files found: ${#R1_FILES[@]}"
        echo "R2 files found: ${#R2_FILES[@]}"
        exit 1
    fi

    # Check that each R1 file has a corresponding R2 file.
    for R1_FILE in "${R1_FILES[@]}"; do

        SAMPLE_NAME="${R1_FILE%"${R1_SUFFIX}"}"
        R2_FILE="${SAMPLE_NAME}${R2_SUFFIX}"

        if [[ ! -f "${R2_FILE}" ]]; then
            echo "ERROR: Missing R2 file for ${R1_FILE}"
            exit 1
        fi

    done

    # Check that R1 and R2 file counts match.
    if [[ ${#R1_FILES[@]} -ne ${#R2_FILES[@]} ]]; then
        echo "ERROR: R1 and R2 file counts do not match."
        exit 1
    fi

    echo "Paired-end FASTQ files validated."
    echo "R1 files: ${#R1_FILES[@]}"
    echo "R2 files: ${#R2_FILES[@]}"

fi


# ---------------------------------------------------------------------
# 4. COMPLETION
# ---------------------------------------------------------------------

echo "============================================"
echo "DATA ACQUISITION COMPLETE"
echo "FASTQ directory: ${RAW_DATA_DIR}"
echo "Read type: ${READ_TYPE}"
echo "============================================"