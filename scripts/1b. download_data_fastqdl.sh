#!/bin/bash
# Bash is the interpreter for this script.

set -euo pipefail
# Helps catch failed commands, unset variables and pipeline failures.

# Define the project directory in WSL.
project_folder="$HOME/rnaseq-ml-alzheimers"

# Create the raw data directory if it does not already exist.
mkdir -p "${project_folder}/raw_data"

# Activate the conda environment containing the required tools.
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate rnaseq

# Check whether fastq-dl is available in the active environment.
if ! command -v fastq-dl >/dev/null 2>&1; then
    echo "Error: fastq-dl is not installed or not in PATH."
    exit 1
fi

# Download all sequencing runs associated with the selected BioProject.
# fastq-dl retrieves the associated run accessions and downloads their sequencing data.
# FASTQ files are compressed using gzip and saved in the raw_data directory.

fastq-dl \
    --accession PRJNA516886 \
    --outdir "${project_folder}/raw_data"

# After downloading, assess read quality with FastQC.
# Continue with adapter trimming, post-trimming QC, STAR alignment,
# alignment quality assessment, and gene quantification with featureCounts.