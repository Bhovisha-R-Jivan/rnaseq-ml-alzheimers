#!/bin/bash 
# Bash is the interpreter for this script

# Project: RNA-seq and Machine Learning for Alzheimer's Disease
# Purpose: Download sequencing data for RNA-seq analysis.
# Author: Bhovisha Jivan

set -euo pipefail
# Helps catch failed commands and unset variables, making the script more robust.

# Define your desktop path 
# Setting variables like desktop_path allows for easy modification of the script without hardcoding paths.
desktop_path='/mnt/c/Users/Bhovisha/OneDrive/Desktop/'

# Create a directory for the project if it doesn't exist to store all project related files.
project_folder="${desktop_path}rnaseq-ml-alzheimers"
mkdir -p "$project_folder"

# Navigate to the created directory.
# This can be easily done using the 'cd' command
cd "$project_folder"

# Create a subdirectory for the data to keep the project organised and separate from other files.
mkdir -p "${project_folder}/data/raw_data"

# Activate your conda environment containing the necessary tools for downloading and processing the data.
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate rnaseq

# There are 3 different ways in which FASTQ files can be downloaded.
# *Method 1:* Using the Sequence Read Archive (SRA) tool.
# Replace "SRR8482616" with your desired SRA accession number for the dataset to be downloaded.
# By specifying the "--outdir" option, you can choose where to save the downloaded FASTQ files. 

#prefetch SRR32699855  \
    #-O "${project_folder}/data/raw"

#fasterq-dump \
   # --outdir ${project_folder}/data/raw/ \
   #"${project_folder}/data/raw/SRR8482616"  # Replace SRR8482616 with your SRA accession number


# *Method 2:* Using the 'fastq-dl' tool.
# This approach can also be used to download FASTQ files from the SRA database or the European Nucleotide Archive (ENA) using accession numbers.
# It simplifies the process by directly retrieving FASTQ files, skipping the SRA-to-FASTQ conversion step.
# fastq-dl is particularly useful for users who want a straightforward and versatile way to obtain FASTQ files as it supports both SRA, ERR and ENA accession numbers.
# *Commands*:
    # fastq-dl
    # --accession ERS5684710 : This is the accession number as displayed on the ENA website
    # --outdir "${project_folder}/data/raw/" : By specifying the "--outdir" option, you can choose where to save the downloaded FASTQ files. 
    #
    # Take not of the following before running the command:
    # fastq-dl needs to be installed in the rnaseq_env in order to work.
    # It can be installed using the command: conda install fastq-dl
    # Make sure your accession number (e.g., ERS5684710) is correct and corresponds to the dataset you want to download.
    # Check that the output directory ('${project_folder}/data/raw/') exists, to ensure all downloaded files are in the same place for easier access and analysis.
    # Sometimes fastq files can be downloaded as zip files but these can easily be uncompressed if needed or remain zipped for further downstream tools.
#fastq-dl \
    #--accession ERS5684710 \
    #--outdir "${project_folder}/data/raw/"
 
# *Method 3:* Using the 'iSeq' tool.
# Also downloads FASTQ files from the SRA or other databases.
# iSeq is more flexible because it supportss various sequencing data and compression formats
# iSeq allows the specification of output formats and compression options saving space as files can be quite large.
# *Commands*:
    # - 'iSeq
    # --gzip : Compresseses FASTQ files in gzip format (.fastq.gz)
    # --fastq: Ensures the output file is FASTQ which is the standard for RNA-seq data.
    # --input SRR32701816: The relevant SRA accession number for the dataset to be downloaded.
    # --database sra: Indicates that the data should be retrieved from the SRA database.
    # -- output "${project_folder}/data/raw/": Specifies the output directory for the downloaded files.
#
# Take not of the following before running the command:
    # iSeq needs to be installed in the rnaseq_env in order to work.
    # Make sure your accession number (e.g., SRR8482616) is correct and corresponds to the dataset you want to download.
    # Check that the output directory ('${project_folder}/data/raw/') exists, to ensure all downloaded files are in the same place for easier access and analysis. 
    # If you are using the --gzip option ensure it is compatible with downstream tools.

#iseq \
#    --gzip \
#    --fastq \
#    --input SRR32701816 \
#    --database sra \
#    --output "${project_folder}/data/raw/"

# After downloading, assess read quality with FastQC.
# Continue with adapter trimming, post-trimming QC, STAR alignment,
# alignment quality assessment, and gene quantification with featureCounts.