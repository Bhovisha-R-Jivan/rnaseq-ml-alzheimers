#!/bin/bash

# =====================================================================
# SCRIPT 1C: RNA-SEQ DATA ACQUISITION USING SRA TOOLKIT
# =====================================================================

set -euo pipefail

# Load pipeline configuration.
source scripts/0.config.sh

# Create raw data and temporary directories.
mkdir -p "${RAW_DATA_DIR}"
TEMP_DIR="${ANALYSIS_DIR}/tmp"
mkdir -p "${TEMP_DIR}"


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

# Check whether required SRA Toolkit programs are available.
if ! command -v prefetch > /dev/null 2>&1; then
    echo "ERROR: prefetch is not installed or unavailable."
    exit 1
fi

if ! command -v fasterq-dump > /dev/null 2>&1; then
    echo "ERROR: fasterq-dump is not installed or unavailable."
    exit 1
fi


# ---------------------------------------------------------------------
# 3. DOWNLOAD SRA FILES
# ---------------------------------------------------------------------

# Download each SRA run accession individually.
while IFS= read -r accession || [[ -n "${accession}" ]]; do

    # Remove Windows carriage returns and surrounding whitespace.
    accession="${accession//$'\r'/}"
    accession="${accession#"${accession%%[![:space:]]*}"}"
    accession="${accession%"${accession##*[![:space:]]}"}"

    # Skip empty lines.
    [[ -z "${accession}" ]] && continue

    echo "Downloading ${accession}..."

    prefetch \
        --output-directory "${RAW_DATA_DIR}" \
        "${accession}"

    # Check whether the downloaded SRA file exists.
    if [[ ! -f "${RAW_DATA_DIR}/${accession}/${accession}.sra" ]]; then
        echo "ERROR: SRA file for ${accession} not found."
        exit 1
    fi

done < "${ACCESSION_FILE}"

echo "SRA downloads completed."


# ---------------------------------------------------------------------
# 4. CONVERT SRA FILES INTO FASTQ
# ---------------------------------------------------------------------

# Convert each downloaded SRA file into FASTQ format.
while IFS= read -r accession || [[ -n "${accession}" ]]; do

    # Remove Windows carriage returns and surrounding whitespace.
    accession="${accession//$'\r'/}"
    accession="${accession#"${accession%%[![:space:]]*}"}"
    accession="${accession%"${accession##*[![:space:]]}"}"

    # Skip empty lines.
    [[ -z "${accession}" ]] && continue

    SRA_FILE="${RAW_DATA_DIR}/${accession}/${accession}.sra"

    echo "Converting ${accession} into FASTQ..."

    # Convert paired-end reads.
    if [[ "${READ_TYPE}" == "paired" ]]; then

        fasterq-dump \
            --split-files \
            --threads "${THREADS}" \
            --temp "${TEMP_DIR}" \
            --outdir "${RAW_DATA_DIR}" \
            "${SRA_FILE}"

    # Convert single-end reads.
    elif [[ "${READ_TYPE}" == "single" ]]; then

        fasterq-dump \
            --threads "${THREADS}" \
            --temp "${TEMP_DIR}" \
            --outdir "${RAW_DATA_DIR}" \
            "${SRA_FILE}"

    fi

    # Identify FASTQ files generated for this accession.
    shopt -s nullglob
    FASTQ_FILES=("${RAW_DATA_DIR}/${accession}"*.fastq)
    shopt -u nullglob

    # Check that FASTQ files were generated.
    if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
        echo "ERROR: No FASTQ files generated for ${accession}."
        exit 1
    fi

    # Compress generated FASTQ files.
    gzip -- "${FASTQ_FILES[@]}"

    echo "Conversion and compression completed for ${accession}."

done < "${ACCESSION_FILE}"


# ---------------------------------------------------------------------
# 5. FINAL INPUT CHECK
# ---------------------------------------------------------------------

echo "Checking generated FASTQ files..."

if ! compgen -G "${RAW_DATA_DIR}/*.fastq.gz" > /dev/null; then
    echo "ERROR: No .fastq.gz files found after conversion."
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

    if [[ ${#R1_FILES[@]} -ne ${#R2_FILES[@]} ]]; then
        echo "ERROR: R1 and R2 file counts do not match."
        exit 1
    fi

    echo "Paired-end FASTQ files validated."
    echo "R1 files: ${#R1_FILES[@]}"
    echo "R2 files: ${#R2_FILES[@]}"

fi


# ---------------------------------------------------------------------
# 6. COMPLETION
# ---------------------------------------------------------------------

echo "============================================"
echo "SRA DOWNLOAD AND CONVERSION COMPLETE"
echo "FASTQ directory: ${RAW_DATA_DIR}"
echo "Read type: ${READ_TYPE}"
echo "============================================"