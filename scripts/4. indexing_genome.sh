#!/bin/bash
set -euo pipefail

# Project directories
PROJECT_DIR="$HOME/rnaseq-ml-alzheimers"
RESOURCE_DIR="${PROJECT_DIR}/resources"
REFERENCE_DIR="${RESOURCE_DIR}/reference_genome"
STAR_INDEX_DIR="${RESOURCE_DIR}/star_index"

# Reference files: GRCh38 genome assembly and GENCODE v47 annotation
GENOME_FASTA="${REFERENCE_DIR}/GRCh38.primary_assembly.genome.fa"
ANNOTATION_GTF="${REFERENCE_DIR}/gencode.v47.primary_assembly.annotation.gtf"

# STAR settings
N_THREADS=2
SJDB_OVERHANG=49

# Create directories
mkdir -p "$REFERENCE_DIR"
mkdir -p "$STAR_INDEX_DIR"

# Activate the conda environment
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate rnaseq

# Check that the reference files exist
if [[ ! -f "$GENOME_FASTA" ]]; then
    echo "ERROR: Reference genome FASTA not found: $GENOME_FASTA"
    echo "Obtain the GRCh38 primary assembly FASTA before running this script."
    exit 1
fi

if [[ ! -f "$ANNOTATION_GTF" ]]; then
    echo "ERROR: GENCODE annotation GTF not found: $ANNOTATION_GTF"
    echo "Obtain the GENCODE v47 primary assembly GTF before running this script."
    exit 1
fi

# Check that STAR is installed
if ! command -v STAR >/dev/null 2>&1; then
    echo "ERROR: STAR is not installed in the rnaseq environment."
    exit 1
fi

# Build the STAR genome index
STAR \
    --runMode genomeGenerate \
    --runThreadN "$N_THREADS" \
    --genomeDir "$STAR_INDEX_DIR" \
    --genomeFastaFiles "$GENOME_FASTA" \
    --sjdbGTFfile "$ANNOTATION_GTF" \
    --sjdbOverhang "$SJDB_OVERHANG"

echo "STAR genome indexing completed successfully."