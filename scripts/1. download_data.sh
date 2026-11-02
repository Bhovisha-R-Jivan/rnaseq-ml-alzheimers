#!/bin/bash
# Bash is the interpreter for this script.

set -euo pipefail
# Helps catch failed commands and unset variables, making the script more robust.

# Define the project directory in WSL.
project_folder="$HOME/rnaseq-ml-alzheimers"

# Create the raw data directory if it does not already exist.
mkdir -p "${project_folder}/raw_data"

# Activate the conda environment containing the required tools.
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate rnaseq

# Download all sequencing runs associated with the selected BioProject.
# iSeq retrieves the associated run accessions and downloads their sequencing data.
# FASTQ files are compressed using gzip and saved in the raw_data directory.

iseq \
    --gzip \
    --fastq \
    --input PRJNA516886 \
    --database sra \
    --output "${project_folder}/raw_data"

# After downloading, assess read quality with FastQC.
# Continue with adapter trimming, post-trimming QC, STAR alignment,
# alignment quality assessment, and gene quantification with featureCounts.