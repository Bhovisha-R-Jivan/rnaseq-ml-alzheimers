#!/bin/bash

# =====================================================================
# SCRIPT 4: STAR GENOME INDEX GENERATION
# =====================================================================

set -euo pipefail

# Load pipeline configuration.
source scripts/0.config.sh


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

# Check that STAR is installed.
if ! command -v STAR > /dev/null 2>&1; then
    echo "ERROR: STAR is not installed or unavailable."
    exit 1
fi


# ---------------------------------------------------------------------
# 2. INPUT VALIDATION
# ---------------------------------------------------------------------

# Check that the reference genome exists and is not empty.
if [[ ! -s "${GENOME}" ]]; then
    echo "ERROR: Reference genome FASTA not found or empty: ${GENOME}"
    exit 1
fi

# Check that the annotation exists and is not empty.
if [[ ! -s "${ANNOTATION}" ]]; then
    echo "ERROR: Annotation GTF not found or empty: ${ANNOTATION}"
    exit 1
fi

# Check that the reference genome is in FASTA format.
if ! grep -q '^>' "${GENOME}"; then
    echo "ERROR: Reference genome does not appear to be in FASTA format."
    exit 1
fi

# Check that the annotation contains GTF entries.
if ! grep -qv '^#' "${ANNOTATION}"; then
    echo "ERROR: Annotation GTF does not contain any feature entries."
    exit 1
fi


# ---------------------------------------------------------------------
# 3. INDEX DIRECTORY VALIDATION
# ---------------------------------------------------------------------

# Do not overwrite an existing STAR index.
if [[ -f "${STAR_INDEX_DIR}/Genome" ]]; then
    echo "ERROR: A STAR genome index already exists in:"
    echo "${STAR_INDEX_DIR}"
    echo "Remove the existing index directory if you intend to rebuild it."
    exit 1
fi

# Do not build over an incomplete or partially existing index.
if [[ -d "${STAR_INDEX_DIR}" ]] &&
   [[ -n "$(ls -A "${STAR_INDEX_DIR}")" ]]; then
    echo "ERROR: STAR index directory is not empty."
    echo "Check or remove the existing directory before proceeding."
    exit 1
fi

# Create the STAR index directory.
mkdir -p "${STAR_INDEX_DIR}"


# ---------------------------------------------------------------------
# 4. CALCULATE GENOME SIZE
# ---------------------------------------------------------------------

# Calculate the total number of nucleotide bases in the FASTA file.
GENOME_SIZE=$(awk '
    /^>/ { next }
    {
        gsub(/[[:space:]]/, "")
        total += length($0)
    }
    END { print total }
' "${GENOME}")

# Check that the FASTA contains sequence data.
if [[ ! "${GENOME_SIZE}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: Unable to determine a valid genome size."
    exit 1
fi

# Calculate STAR genomeSAindexNbases.
# Formula: min(14, log2(genome size)/2 - 1)
SA_INDEX_BASES=$(awk -v size="${GENOME_SIZE}" '
    BEGIN {
        value = log(size) / log(2) / 2 - 1
        if (value > 14) value = 14
        if (value < 1) value = 1
        print int(value)
    }
')

echo "Genome size: ${GENOME_SIZE} bases"
echo "STAR SA index bases: ${SA_INDEX_BASES}"


# ---------------------------------------------------------------------
# 5. BUILD STAR GENOME INDEX
# ---------------------------------------------------------------------

echo "Starting STAR genome index generation."

STAR \
    --runMode genomeGenerate \
    --runThreadN "${THREADS}" \
    --genomeDir "${STAR_INDEX_DIR}" \
    --genomeFastaFiles "${GENOME}" \
    --sjdbGTFfile "${ANNOTATION}" \
    --sjdbOverhang "${SJDB_OVERHANG}" \
    --genomeSAindexNbases "${SA_INDEX_BASES}"


# ---------------------------------------------------------------------
# 6. VERIFY INDEX GENERATION
# ---------------------------------------------------------------------

if [[ ! -s "${STAR_INDEX_DIR}/Genome" ]]; then
    echo "ERROR: STAR genome index generation could not be verified."
    exit 1
fi

if [[ ! -s "${STAR_INDEX_DIR}/SA" ]]; then
    echo "ERROR: STAR suffix array file was not generated."
    exit 1
fi

echo "STAR index files verified."


# ---------------------------------------------------------------------
# 7. COMPLETION
# ---------------------------------------------------------------------

echo "============================================================"
echo "STAR GENOME INDEX GENERATION COMPLETED"
echo "Reference genome: ${GENOME}"
echo "Annotation: ${ANNOTATION}"
echo "Genome size: ${GENOME_SIZE} bases"
echo "SA index bases: ${SA_INDEX_BASES}"
echo "STAR index directory: ${STAR_INDEX_DIR}"
echo "============================================================"