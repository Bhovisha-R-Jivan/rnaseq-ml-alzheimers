#!/bin/bash

# =====================================================================
# SCRIPT 5: STAR RNA-SEQ ALIGNMENT
# =====================================================================

set -euo pipefail

# Load pipeline configuration.
source scripts/0.config.sh

# Define input and output directories.
TRIMMED_DIR="${ANALYSIS_DIR}/2_trimmed"
BAM_DIR="${ANALYSIS_DIR}/4_bamfiles"

# Create BAM output directory.
mkdir -p "${BAM_DIR}"

# =====================================================================
# 1. INPUT VALIDATION
# =====================================================================

# Check that the trimmed FASTQ directory exists.
if [[ ! -d "${TRIMMED_DIR}" ]]; then
    echo "ERROR: Trimmed FASTQ directory not found: ${TRIMMED_DIR}"
    exit 1
fi

# Check that the STAR genome index exists.
if [[ ! -s "${STAR_INDEX_DIR}/Genome" ]]; then
    echo "ERROR: STAR genome index not found or incomplete:"
    echo "${STAR_INDEX_DIR}"
    echo "Build the STAR index before running alignment."
    exit 1
fi

# Check that STAR and SAMtools are installed.
for TOOL in STAR samtools; do
    if ! command -v "${TOOL}" > /dev/null 2>&1; then
        echo "ERROR: ${TOOL} is not installed or unavailable."
        exit 1
    fi
done

# Enable nullglob to handle missing FASTQ patterns.
shopt -s nullglob


# =====================================================================
# 2. PAIRED-END ALIGNMENT
# =====================================================================

if [[ "${READ_TYPE}" == "paired" ]]; then

    echo "Starting paired-end STAR alignment."

    # Identify trimmed forward and reverse reads.
    R1_FILES=( "${TRIMMED_DIR}"/*_val_1.fq.gz )
    R2_FILES=( "${TRIMMED_DIR}"/*_val_2.fq.gz )

    # Check that both read sets exist.
    if [[ ${#R1_FILES[@]} -eq 0 || ${#R2_FILES[@]} -eq 0 ]]; then
        echo "ERROR: Trimmed paired-end FASTQ files are missing."
        exit 1
    fi

    # Check that forward and reverse read counts match.
    if [[ ${#R1_FILES[@]} -ne ${#R2_FILES[@]} ]]; then
        echo "ERROR: Number of forward and reverse reads does not match."
        exit 1
    fi

    # Process each paired-end sample.
    for R1 in "${R1_FILES[@]}"; do

        # Identify the corresponding reverse read.
        SAMPLE_NAME="$(basename "${R1}" _val_1.fq.gz)"
        R2="${TRIMMED_DIR}/${SAMPLE_NAME}_val_2.fq.gz"

        # Check that the reverse read exists.
        if [[ ! -s "${R2}" ]]; then
            echo "ERROR: Matching reverse read not found for ${R1}"
            exit 1
        fi

        # Define STAR output prefix and BAM output.
        OUT_PREFIX="${BAM_DIR}/${SAMPLE_NAME}_"
        SORTED_BAM="${BAM_DIR}/${SAMPLE_NAME}_Aligned.sortedByCoord.out.bam"

        echo "Aligning paired-end sample: ${SAMPLE_NAME}"

        # Align paired-end reads and generate coordinate-sorted BAM.
        STAR \
            --runThreadN "${THREADS}" \
            --genomeDir "${STAR_INDEX_DIR}" \
            --readFilesIn "${R1}" "${R2}" \
            --readFilesCommand zcat \
            --outSAMunmapped Within \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix "${OUT_PREFIX}"

        # Verify BAM generation.
        if [[ ! -s "${SORTED_BAM}" ]]; then
            echo "ERROR: STAR did not produce a valid BAM for ${SAMPLE_NAME}"
            exit 1
        fi

        # Index sorted BAM file.
        samtools index \
            --threads "${THREADS}" \
            "${SORTED_BAM}"

        echo "Completed alignment for sample: ${SAMPLE_NAME}"

    done


# =====================================================================
# 3. SINGLE-END ALIGNMENT
# =====================================================================

elif [[ "${READ_TYPE}" == "single" ]]; then

    echo "Starting single-end STAR alignment."

    # Identify trimmed single-end FASTQ files.
    FASTQ_FILES=( "${TRIMMED_DIR}"/*_trimmed.fq.gz )

    # Check that FASTQ files exist.
    if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
        echo "ERROR: No trimmed single-end FASTQ files found."
        exit 1
    fi

    # Process each single-end sample.
    for INPUT_FASTQ in "${FASTQ_FILES[@]}"; do

        # Extract sample name.
        SAMPLE_NAME="$(basename "${INPUT_FASTQ}" _trimmed.fq.gz)"

        # Define STAR output prefix and BAM output.
        OUT_PREFIX="${BAM_DIR}/${SAMPLE_NAME}_"
        SORTED_BAM="${BAM_DIR}/${SAMPLE_NAME}_Aligned.sortedByCoord.out.bam"

        echo "Aligning single-end sample: ${SAMPLE_NAME}"

        # Align single-end reads and generate coordinate-sorted BAM.
        STAR \
            --runThreadN "${THREADS}" \
            --genomeDir "${STAR_INDEX_DIR}" \
            --readFilesIn "${INPUT_FASTQ}" \
            --readFilesCommand zcat \
            --outSAMunmapped Within \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix "${OUT_PREFIX}"

        # Verify BAM generation.
        if [[ ! -s "${SORTED_BAM}" ]]; then
            echo "ERROR: STAR did not produce a valid BAM for ${SAMPLE_NAME}"
            exit 1
        fi

        # Index sorted BAM file.
        samtools index \
            --threads "${THREADS}" \
            "${SORTED_BAM}"

        echo "Completed alignment for sample: ${SAMPLE_NAME}"

    done

else

    echo "ERROR: Unsupported READ_TYPE: ${READ_TYPE}"
    exit 1

fi


# =====================================================================
# 4. COMPLETION
# =====================================================================

echo "============================================================"
echo "STAR ALIGNMENT COMPLETED SUCCESSFULLY"
echo "BAM output directory: ${BAM_DIR}"
echo "============================================================"