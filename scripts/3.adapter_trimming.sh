#!/bin/bash

# =====================================================================
# SCRIPT 3: ADAPTER TRIMMING AND POST-TRIMMING QUALITY ASSESSMENT
# =====================================================================

set -euo pipefail

# Load pipeline configuration.
source scripts/0.config.sh

# Define output directories.
TRIMMED_DIR="${ANALYSIS_DIR}/2_trimmed"
POSTQC_DIR="${ANALYSIS_DIR}/3_postqc_fastqc"
POSTMULTIQC_DIR="${ANALYSIS_DIR}/3_postqc_multiqc"

# Create output directories.
mkdir -p "${TRIMMED_DIR}" "${POSTQC_DIR}" "${POSTMULTIQC_DIR}"


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

# Check that required tools are installed.
for TOOL in trim_galore fastqc multiqc; do
    if ! command -v "${TOOL}" > /dev/null 2>&1; then
        echo "ERROR: ${TOOL} is not installed or unavailable."
        exit 1
    fi
done


# ---------------------------------------------------------------------
# 2. CHECK INPUT DATA
# ---------------------------------------------------------------------

# Check that the raw data directory exists.
if [[ ! -d "${RAW_DATA_DIR}" ]]; then
    echo "ERROR: Raw data directory not found: ${RAW_DATA_DIR}"
    exit 1
fi

# Enable nullglob to handle missing file patterns.
shopt -s nullglob


# =====================================================================
# 3. PAIRED-END DATA
# =====================================================================

if [[ "${READ_TYPE}" == "paired" ]]; then

    echo "Starting paired-end adapter trimming."

    # Identify forward and reverse FASTQ files.
    R1_FILES=("${RAW_DATA_DIR}"/*"${R1_SUFFIX}")
    R2_FILES=("${RAW_DATA_DIR}"/*"${R2_SUFFIX}")

    # Check that both read sets exist.
    if [[ ${#R1_FILES[@]} -eq 0 || ${#R2_FILES[@]} -eq 0 ]]; then
        echo "ERROR: Paired-end FASTQ files are missing."
        exit 1
    fi

    # Check that the number of forward and reverse files matches.
    if [[ ${#R1_FILES[@]} -ne ${#R2_FILES[@]} ]]; then
        echo "ERROR: Number of forward and reverse FASTQ files does not match."
        exit 1
    fi

    # Process each paired-end sample.
    for R1 in "${R1_FILES[@]}"; do

        SAMPLE="${R1%"${R1_SUFFIX}"}"
        R2="${SAMPLE}${R2_SUFFIX}"

        # Check that the corresponding reverse read exists.
        if [[ ! -f "${R2}" ]]; then
            echo "ERROR: Missing reverse read: ${R2}"
            exit 1
        fi

        # Extract sample name.
        SAMPLE_NAME="${SAMPLE##*/}"

        echo "Trimming paired-end sample: ${SAMPLE_NAME}"

        trim_galore \
            --paired \
            --cores "${THREADS}" \
            --quality "${TRIM_QUALITY}" \
            --length "${MIN_LENGTH}" \
            --gzip \
            --output_dir "${TRIMMED_DIR}" \
            --basename "${SAMPLE_NAME}" \
            "${R1}" \
            "${R2}"

    done


# =====================================================================
# 4. SINGLE-END DATA
# =====================================================================

elif [[ "${READ_TYPE}" == "single" ]]; then

    echo "Starting single-end adapter trimming."

    # Identify single-end FASTQ files.
    FASTQ_FILES=("${RAW_DATA_DIR}"/*.fastq.gz)

    # Check that FASTQ files exist.
    if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
        echo "ERROR: No FASTQ files found in ${RAW_DATA_DIR}"
        exit 1
    fi

    # Process each single-end sample.
    for INPUT_FASTQ in "${FASTQ_FILES[@]}"; do

        SAMPLE_NAME="$(basename "${INPUT_FASTQ}" .fastq.gz)"

        echo "Trimming single-end sample: ${SAMPLE_NAME}"

        trim_galore \
            --cores "${THREADS}" \
            --quality "${TRIM_QUALITY}" \
            --length "${MIN_LENGTH}" \
            --gzip \
            --output_dir "${TRIMMED_DIR}" \
            --basename "${SAMPLE_NAME}" \
            "${INPUT_FASTQ}"

    done

else

    echo "ERROR: Unsupported READ_TYPE: ${READ_TYPE}"
    exit 1

fi


# =====================================================================
# 5. POST-TRIMMING FASTQ CHECK
# =====================================================================

echo "Checking trimmed FASTQ files..."

# Identify trimmed FASTQ files.
TRIMMED_FILES=(
    "${TRIMMED_DIR}"/*.fq.gz
    "${TRIMMED_DIR}"/*.fastq.gz
)

# Check that trimming produced FASTQ files.
if [[ ${#TRIMMED_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No trimmed FASTQ files found."
    exit 1
fi

echo "Trimmed FASTQ files detected: ${#TRIMMED_FILES[@]}"


# =====================================================================
# 6. POST-TRIMMING QUALITY ASSESSMENT
# =====================================================================

echo "Starting post-trimming FastQC."

# Run FastQC on trimmed reads.
fastqc \
    --outdir "${POSTQC_DIR}" \
    --threads "${THREADS}" \
    "${TRIMMED_FILES[@]}"


# =====================================================================
# 7. POST-TRIMMING MULTIQC REPORT
# =====================================================================

echo "Generating post-trimming MultiQC report..."

# Summarise post-trimming FastQC results.
multiqc \
    "${POSTQC_DIR}" \
    --outdir "${POSTMULTIQC_DIR}" \
    --force

echo "Post-trimming MultiQC report generated."

shopt -u nullglob


# =====================================================================
# 8. COMPLETION
# =====================================================================

echo "============================================================"
echo "ADAPTER TRIMMING COMPLETED"
echo "POST-TRIMMING FASTQC COMPLETED"
echo "POST-TRIMMING MULTIQC COMPLETED"
echo "Trimmed FASTQ files: ${TRIMMED_DIR}"
echo "FastQC reports: ${POSTQC_DIR}"
echo "MultiQC report: ${POSTMULTIQC_DIR}"
echo "============================================================"