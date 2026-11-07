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

echo "============================================"
echo "DATA ACQUISITION COMPLETE"
echo "FASTQ directory: ${RAW_DATA_DIR}"
echo "============================================"