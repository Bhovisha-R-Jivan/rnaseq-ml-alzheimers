#!/bin/bash
# Bash is the interpreter for this script.

set -euo pipefail
# Helps catch failed commands, unset variables and pipeline failures.

# Define the project directory in WSL.
project_folder="$HOME/rnaseq-ml-alzheimers"

# Define the accession list containing the SRR run accessions.
ACCESSION_FILE="${project_folder}/docs/SRR_Acc_List.txt"

# Define the directory for downloaded SRA files.
RAW_DATA_DIR="${project_folder}/raw_data"

# Create the raw data directory if it does not already exist.
mkdir -p "$RAW_DATA_DIR"

# Activate the conda environment containing the required tools.
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate rnaseq

# Check whether prefetch is available.
if ! command -v prefetch >/dev/null 2>&1; then
    echo "Error: prefetch is not installed or not in PATH."
    exit 1
fi

# Check whether fasterq-dump is available.
if ! command -v fasterq-dump >/dev/null 2>&1; then
    echo "Error: fasterq-dump is not installed or not in PATH."
    exit 1
fi

# Check whether the accession list exists.
if [[ ! -f "$ACCESSION_FILE" ]]; then
    echo "Error: Accession list not found."
    exit 1
fi

# Download SRA files for each run accession listed in the accession file.

while IFS= read -r accession; do
    # Skip empty lines.
    [[ -z "$accession" ]] && continue

    echo "Downloading ${accession}..."

    prefetch \
        --output-directory "$RAW_DATA_DIR" \
        "$accession"

done < "$ACCESSION_FILE"

echo "SRA downloads completed."

# fasterq-dump can use a susbtantial amount of temporary disk space, and by creating a temporary directory it becomes eaiser to monitor.

TEMP_DIR="${project_folder}/analysis/tmp"
mkdir -p "$TEMP_DIR"

# Convert downloaded SRA files into FASTQ files.
while IFS= read -r accession; do
    [[ -z "$accession" ]] && continue

    echo "Converting ${accession} into FASTQ..."

    fasterq-dump \
        --split-files \
        --threads 4 \
        --temp "$TEMP_DIR" \
        --outdir "$RAW_DATA_DIR" \
        "${RAW_DATA_DIR}/${accession}/${accession}.sra"

done < "$ACCESSION_FILE"

# Compress the generated FASTQ files using gzip.
gzip "$RAW_DATA_DIR"/*.fastq

# After downloading, assess read quality with FastQC.
# Continue with adapter trimming, post-trimming QC, STAR alignment,
# alignment quality assessment, and gene quantification with featureCounts.