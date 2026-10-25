#!/bin/bash
# Bash is the interpreter for this script.

set -euo pipefail
# Stop if a command fails, an unset variable is used, or a pipeline fails.


# Define your desktop path 
desktop_path='/mnt/c/Users/Bhovisha/OneDrive/Desktop/'
project_folder="$HOME/rnaseq-ml-alzheimers"
raw_data_dir="${project_folder}/raw_data"
# Navigate to the project_folder
# This can be easily done using the 'cd' command
cd "$project_folder"

# Create a directory for FastQC reports to be saved
QC_DIR="${project_folder}/analysis/1_fastqc"
mkdir -p "$QC_DIR"

# Activate your conda environment containing the necessary tools for downloading and processing the data.
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate rnaseq

# Set the number of threads FastQC will use.
N_THREADS=2

# Identify the sample to assess.
SAMPLE_NAME="ERR188044_chrX"


# Run FastQC on both paired-end raw FASTQ files.
fastqc \
    --outdir "$QC_DIR" \
    --threads "$N_THREADS" \
    "${raw_data_dir}/${SAMPLE_NAME}_1.fastq.gz" \
    "${raw_data_dir}/${SAMPLE_NAME}_2.fastq.gz"
    