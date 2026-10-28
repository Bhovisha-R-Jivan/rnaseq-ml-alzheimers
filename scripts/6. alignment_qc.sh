#!/bin/bash
set -euo pipefail

# Project directories
PROJECT_DIR="$HOME/rnaseq-ml-alzheimers"
RESOURCE_DIR="${PROJECT_DIR}/resources"
REFERENCE_DIR="${RESOURCE_DIR}/reference_genome"
BAM_DIR="${PROJECT_DIR}/analysis/4_bamfiles"
QC_DIR="${PROJECT_DIR}/analysis/5_alignmentQC"

# Reference annotation
ANNOTATION_GTF="${REFERENCE_DIR}/gencode.v47.primary_assembly.annotation.gtf"

# Sample accession list
ACCESSION_FILE="${PROJECT_DIR}/docs/SRR_Acc_List.txt"

# Number of threads
N_THREADS=2

# Create output directory
mkdir -p "$QC_DIR"

# Activate the Conda environment
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate rnaseq

# Check that the annotation file exists
if [[ ! -f "$ANNOTATION_GTF" ]]; then
    echo "ERROR: GENCODE annotation file not found: $ANNOTATION_GTF"
    exit 1
fi

# Check that the accession list exists
if [[ ! -f "$ACCESSION_FILE" ]]; then
    echo "ERROR: Accession list not found: $ACCESSION_FILE"
    exit 1
fi

# Check that Qualimap is installed
if ! command -v qualimap >/dev/null 2>&1; then
    echo "ERROR: Qualimap is not installed in the rnaseq environment."
    exit 1
fi

# Run alignment quality control for each sample
while IFS= read -r accession || [[ -n "$accession" ]]; do

    # Skip empty lines
    [[ -z "$accession" ]] && continue

    # Input sorted BAM file
    BAM_FILE="${BAM_DIR}/${accession}_Aligned.sorted.bam"

    # Output directory for this sample
    SAMPLE_QC_DIR="${QC_DIR}/${accession}"

    # Check that the BAM file exists
    if [[ ! -f "$BAM_FILE" ]]; then
        echo "ERROR: BAM file not found: $BAM_FILE"
        exit 1
    fi

    # Create sample-specific output directory
    mkdir -p "$SAMPLE_QC_DIR"

    echo "Running Qualimap for sample: ${accession}"

    # Assess RNA-seq alignment quality
    qualimap rnaseq \
        -bam "$BAM_FILE" \
        -gtf "$ANNOTATION_GTF" \
        -outdir "$SAMPLE_QC_DIR" \
        -outfile "${accession}" \
        -outformat HTML \
        -sorted

    echo "Completed alignment QC for sample: ${accession}"

done < "$ACCESSION_FILE"

echo "Alignment quality-control workflow completed successfully."