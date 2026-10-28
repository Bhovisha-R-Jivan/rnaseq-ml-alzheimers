#!/bin/bash
set -euo pipefail

# Project directories
PROJECT_DIR="$HOME/rnaseq-ml-alzheimers"
REFERENCE_DIR="${PROJECT_DIR}/resources/reference_genome"
BAM_DIR="${PROJECT_DIR}/analysis/4_bamfiles"
QUANT_DIR="${PROJECT_DIR}/analysis/6_quantification/featureCounts"

# Reference annotation
ANNOTATION_GTF="${REFERENCE_DIR}/gencode.v47.primary_assembly.annotation.gtf"

# Output files
COUNTS_FILE="${QUANT_DIR}/featurecounts.txt"
COUNTS_TSV="${QUANT_DIR}/featurecounts.tsv"

# featureCounts settings
N_THREADS=2
FEATURE="exon"
ATTRIBUTE="gene_id"
STRAND=0  # Provisional: unstranded counting

# Create output directory
mkdir -p "$QUANT_DIR"

# Activate the Conda environment
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate rnaseq

# Check that the annotation file exists
if [[ ! -f "$ANNOTATION_GTF" ]]; then
    echo "ERROR: GENCODE annotation file not found: $ANNOTATION_GTF"
    exit 1
fi

# Check that the BAM directory exists
if [[ ! -d "$BAM_DIR" ]]; then
    echo "ERROR: BAM directory not found: $BAM_DIR"
    exit 1
fi

# Check that at least one sorted BAM file exists
shopt -s nullglob
BAM_FILES=("${BAM_DIR}"/*_Aligned.sorted.bam)

if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No sorted BAM files found in $BAM_DIR"
    exit 1
fi

# Check that featureCounts is installed
if ! command -v featureCounts >/dev/null 2>&1; then
    echo "ERROR: featureCounts is not installed in the rnaseq environment."
    exit 1
fi

echo "Running featureCounts..."

# Count reads assigned to genes
featureCounts \
    -T "$N_THREADS" \
    -t "$FEATURE" \
    -g "$ATTRIBUTE" \
    -s "$STRAND" \
    -a "$ANNOTATION_GTF" \
    -o "$COUNTS_FILE" \
    "${BAM_FILES[@]}"

# Verify that featureCounts produced its output
if [[ ! -s "$COUNTS_FILE" ]]; then
    echo "ERROR: featureCounts output is missing or empty: $COUNTS_FILE"
    exit 1
fi

# Remove comment lines and simplify sample column names
awk '
    BEGIN { FS = OFS = "\t" }
    /^#/ { next }
    {
        for (i = 1; i <= NF; i++) {
            sub(/^.*\//, "", $i)
            sub(/_Aligned\.sorted\.bam$/, "", $i)
        }
        print
    }
' "$COUNTS_FILE" > "$COUNTS_TSV"

echo "featureCounts completed successfully."
echo "Raw output: $COUNTS_FILE"
echo "Processed count table: $COUNTS_TSV"