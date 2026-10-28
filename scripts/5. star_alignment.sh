#!/bin/bash
set -euo pipefail

# Project directories
PROJECT_DIR="$HOME/rnaseq-ml-alzheimers"
REFERENCE_DIR="${PROJECT_DIR}/resources/reference_genome"
STAR_INDEX_DIR="${PROJECT_DIR}/resources/star_index"
TRIMMED_DIR="${PROJECT_DIR}/analysis/2_trimmed"
BAM_DIR="${PROJECT_DIR}/analysis/4_bamfiles"

# Reference genome
GENOME_FASTA="${REFERENCE_DIR}/GRCh38.primary_assembly.genome.fa"

# Sample accession list
ACCESSION_FILE="${PROJECT_DIR}/docs/SRR_Acc_List.txt"

# STAR settings
N_THREADS=2

# Create output directory
mkdir -p "$BAM_DIR"

# Activate the Conda environment
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate rnaseq

# Check that the STAR index exists
if [[ ! -f "${STAR_INDEX_DIR}/Genome" ]]; then
    echo "ERROR: STAR genome index not found in ${STAR_INDEX_DIR}"
    echo "Build the STAR index before running alignment."
    exit 1

# Check that the reference genome exists
if [[ ! -f "$GENOME_FASTA" ]]; then
    echo "ERROR: Reference genome FASTA not found: $GENOME_FASTA"
    exit 1

# Check that the accession list exists
if [[ ! -f "$ACCESSION_FILE" ]]; then
    echo "ERROR: Accession list not found: $ACCESSION_FILE"
    exit 1

# Check that STAR and SAMtools are installed
if ! command -v STAR >/dev/null 2>&1; then
    echo "ERROR: STAR is not installed in the rnaseq environment."
    exit 1

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: SAMtools is not installed in the rnaseq environment."
    exit 1

# Align single-end, quality-controlled reads using STAR
while IFS= read -r accession || [[ -n "$accession" ]]; do

    # Skip empty lines
    [[ -z "$accession" ]] && continue

    # Input FASTQ file produced by Trim Galore
    INPUT_FASTQ="${TRIMMED_DIR}/${accession}_trimmed.fq.gz"

    # Output file prefix
    OUT_PREFIX="${BAM_DIR}/${accession}_"

    # Final sorted BAM file
    SORTED_BAM="${BAM_DIR}/${accession}_Aligned.sorted.bam"

    # Check that the trimmed FASTQ file exists
    if [[ ! -f "$INPUT_FASTQ" ]]; then
        echo "ERROR: Trimmed FASTQ file not found: $INPUT_FASTQ"
        exit 1

    echo "Aligning sample: ${accession}"

    # Align reads to the reference genome
    STAR \
        --genomeDir "$STAR_INDEX_DIR" \
        --readFilesIn "$INPUT_FASTQ" \
        --readFilesCommand zcat \
        --outSAMunmapped Within \
        --outFileNamePrefix "$OUT_PREFIX" \
        --runThreadN "$N_THREADS"

    # Sort the STAR alignment and convert SAM to BAM
    samtools sort \
        --reference "$GENOME_FASTA" \
        --threads "$N_THREADS" \
        -O BAM \
        -o "$SORTED_BAM" \
        "${OUT_PREFIX}Aligned.out.sam"

    # Index the sorted BAM file
    samtools index \
        --threads "$N_THREADS" \
        "$SORTED_BAM"

    # Remove the intermediate SAM file to save disk space
    rm "${OUT_PREFIX}Aligned.out.sam"

    echo "Completed alignment for sample: ${accession}"

done < "$ACCESSION_FILE"

echo "STAR alignment workflow completed successfully."
