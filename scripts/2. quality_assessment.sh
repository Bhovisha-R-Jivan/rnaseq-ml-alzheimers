#!/bin/bash
# Bash is the interpreter for this script

set -euo pipefail
# Helps catch failed commands and unset variables, making the script more robust.

# Define the directories for the project, raw data, and quality control (QC) outputs.
raw_data="$HOME/rnaseq-ml-alzheimers/raw_data"
QC_DIR="$HOME/rnaseq-ml-alzheimers/analysis/1_fastqc"
ACCESSION_FILE="$HOME/rnaseq-ml-alzheimers/docs/SRR_Acc_List.txt"

mkdir -p "$QC_DIR"

# Activate the conda environment containing the necessary tools.
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate rnaseq

# Check that the accession list exists.
if [[ ! -f "$ACCESSION_FILE" ]]; then
    echo "ERROR: Accession list not found: $ACCESSION_FILE"
    exit 1
fi

# Check that FastQC is installed.
if ! command -v fastqc >/dev/null 2>&1; then
    echo "ERROR: FastQC is not installed in the rnaseq environment."
    exit 1
fi

# Set the number of threads FastQC will use.
N_THREADS=2

# Loop through the SRA run accessions in the accession list.
while IFS= read -r accession || [[ -n "$accession" ]]; do

    # Skip blank lines.
    [[ -z "$accession" ]] && continue

    # Define the single-end FASTQ filename produced by iSeq.
    INPUT_FASTQ="${raw_data}/${accession}_1.fastq.gz"

    # Check that the FASTQ file exists before running FastQC.
    if [[ ! -f "$INPUT_FASTQ" ]]; then
        echo "ERROR: FASTQ file not found: $INPUT_FASTQ"
        exit 1
    fi

    # Run FastQC on the single-end FASTQ file for the current accession.
    fastqc \
        --outdir "$QC_DIR" \
        --threads "$N_THREADS" \
        "$INPUT_FASTQ"

done < "$ACCESSION_FILE"

echo "FastQC quality assessment completed successfully."