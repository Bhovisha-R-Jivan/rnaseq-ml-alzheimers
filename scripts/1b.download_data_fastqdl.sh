#!/bin/bash

# =====================================================================
# SCRIPT 1B: RNA-SEQ DATA ACQUISITION USING FASTQ-DL
# =====================================================================

set -euo pipefail

# Load pipeline configuration.
source scripts/0.config.sh

# Create raw data directory if it does not exist.
mkdir -p "${RAW_DATA_DIR}"


# ---------------------------------------------------------------------
# 1. CHECK ACCESSION FILE
# ---------------------------------------------------------------------

# Check that the accession file exists and is not empty.
if [[ ! -f "${ACCESSION_FILE}" ]]; then
    echo "ERROR: Accession file not found: ${ACCESSION_FILE}"
    exit 1
fi

if [[ ! -s "${ACCESSION_FILE}" ]]; then
    echo "ERROR: Accession file is empty."
    exit 1
fi


# ---------------------------------------------------------------------
# 2. ACTIVATE CONDA ENVIRONMENT
# ---------------------------------------------------------------------

# Activate Conda environment.
if [[ ! -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
    echo "ERROR: Conda initialization script not found."
    exit 1
fi

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate rnaseq

# Check whether fastq-dl is available.
if ! command -v fastq-dl > /dev/null 2>&1; then
    echo "ERROR: fastq-dl is not installed or unavailable."
    exit 1
fi


# ---------------------------------------------------------------------
# 3. DOWNLOAD SRA RUN ACCESSIONS
# ---------------------------------------------------------------------

# Read and download each individual SRA run accession.
while IFS= read -r accession || [[ -n "${accession}" ]]; do

    # Remove Windows carriage returns and surrounding whitespace.
    accession="${accession//$'\r'/}"
    accession="${accession#"${accession%%[![:space:]]*}"}"
    accession="${accession%"${accession##*[![:space:]]}"}"

    # Skip empty lines.
    [[ -z "${accession}" ]] && continue

    echo "Downloading ${accession}..."

    fastq-dl \
        --accession "${accession}" \
        --outdir "${RAW_DATA_DIR}"

done < "${ACCESSION_FILE}"


# ---------------------------------------------------------------------
# 4. FINAL INPUT CHECK
# ---------------------------------------------------------------------

echo "Checking downloaded FASTQ files..."

# Check that gzipped FASTQ files were downloaded.
if ! compgen -G "${RAW_DATA_DIR}/*.fastq.gz" > /dev/null; then
    echo "ERROR: No .fastq.gz files found after downloading."
    exit 1
fi

echo "FASTQ files detected."


# ---------------------------------------------------------------------
# 5. COMPLETION
# ---------------------------------------------------------------------

echo "============================================"
echo "FASTQ-DL DOWNLOAD COMPLETE"
echo "FASTQ directory: ${RAW_DATA_DIR}"
echo "============================================"