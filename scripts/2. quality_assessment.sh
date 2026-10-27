#!/bin/bash 
# Bash is the interpreter for this script

set -euo pipefail
# Helps catch failed commands and unset variables, making the script more robust.

# Define the directories for the project, raw data, and quality control (QC) outputs.
raw_data="$HOME/rnaseq-ml-alzheimers/raw_data" 
QC_DIR="$HOME/rnaseq-ml-alzheimers/analysis/1_fastqc"
mkdir -p "$QC_DIR"
ACCESSION_FILE="$HOME/rnaseq-ml-alzheimers/docs/SRR_Acc_List.txt"

# Activate your conda environment containing the necessary tools for downloading and processing the data.
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate rnaseq

# Set the number of threads FastQC will use.
N_THREADS=2

# Loop through the SRA run accessions in the accession list. 
while IFS= read -r accession || [[ -n "$accession" ]]; do

    # Run FastQC on the single-end FASTQ file for the current accession.
    fastqc \
        --outdir "$QC_DIR" \
        --threads "$N_THREADS" \
        "${raw_data}/${accession}.fastq.gz" 
        
done < "$ACCESSION_FILE"