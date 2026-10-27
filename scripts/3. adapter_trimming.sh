```bash
#!/bin/bash
# Bash is the interpreter for this script.

set -euo pipefail
# Helps catch failed commands and unset variables, making the script more robust.

# Define the directories for the project, raw data, and trimming outputs.
raw_data="$HOME/rnaseq-ml-alzheimers/raw_data"
ACCESSION_FILE="$HOME/rnaseq-ml-alzheimers/docs/SRR_Acc_List.txt"
TRIMMED_DIR="$HOME/rnaseq-ml-alzheimers/analysis/2_trimmed"
POSTQC_DIR="$HOME/rnaseq-ml-alzheimers/analysis/3_postqc_fastqc"

# Create the output directories if they do not already exist.
mkdir -p "$TRIMMED_DIR"
mkdir -p "$POSTQC_DIR"

# Activate the conda environment containing Trim Galore and FastQC.
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate rnaseq

# Set the number of CPU cores for processing.
N_THREADS=2

# Loop through the SRA run accessions in the accession list.
while IFS= read -r accession || [[ -n "$accession" ]]; do

    # Define the input FASTQ file for the current single-end run.
    INPUT_FASTQ="${raw_data}/${accession}.fastq.gz"

    # Trim adapters and low-quality bases, then assess the trimmed reads with FastQC.
    trim_galore \
        --cores "$N_THREADS" \
        --quality 30 \
        --length 20 \
        --gzip \
        --output_dir "$TRIMMED_DIR" \
        --basename "$accession" \
        --fastqc \
        --fastqc_args "--outdir $POSTQC_DIR" \
        "$INPUT_FASTQ"

done < "$ACCESSION_FILE"
