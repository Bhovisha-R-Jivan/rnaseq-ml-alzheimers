#!/bin/bash 
# Bash is the interpreter for this script

set -euo pipefail
# Helps catch failed commands and unset variables, making the script more robust.

# Define the project directory in WSL
# Setting variables like project_folder allows for easy modification of the script without hardcoding paths.
project_folder="$HOME/rnaseq-ml-alzheimers"

# Create a subdirectory for the data to keep the project organised and separate from other files.
mkdir -p "${project_folder}/raw_data"

# The accession list has been sourced from the NCBI SRA run selector and BioProject PRJNA516886
# All 289 SRA run accession numbers associated with the selected BioProject are in the file SRR_Acc_List.txt.
# The list allows the download commands to iterate over the runs without hardcoding every accession into the script.
# Depending on the method used, 3 loops have been provided to download the FASTQ file for each sample.
ACCESSION_FILE="${project_folder}/docs/SRR_Acc_List.txt"

# Activate your conda environment containing the necessary tools for downloading and processing the data.
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate rnaseq

# There are 3 different ways in which FASTQ files can be downloaded.
# *Method 1:* Using the Sequence Read Archive (SRA) tool. The data is first downloaded with prefetch and then converted to FASTQ format with fasterq-dump.
# By specifying the "--outdir" option, you can choose where to save the downloaded FASTQ files. 

# while IFS= read -r accession || [[ -n "$accession" ]]; do
# Download the current SRA run.
    # prefetch "$accession" \
        # --output-directory "${project_folder}/raw_data"

    # Convert the downloaded SRA run into FASTQ files.
    # fasterq-dump \
        # --outdir "${project_folder}/raw_data" \
        # --split-files \
        # "${project_folder}/raw_data/${accession}/${accession}.sra"

# done < "$ACCESSION_FILE"

# *Method 2:* Using the 'fastq-dl' tool.
# This approach can also be used to download FASTQ files from the SRA database or the European Nucleotide Archive (ENA) using accession numbers.
# It simplifies the process by directly retrieving FASTQ files, skipping the SRA-to-FASTQ conversion step.
# fastq-dl is particularly useful for users who want a straightforward and versatile way to obtain FASTQ files as it supports both SRA, ERR and ENA accession numbers.
# *Commands*:
    # fastq-dl
    # --accession SRR8482616 : This is the accession number as displayed on the NCBI SRA run selection accession list.
    # --outdir "${project_folder}/raw_data/" : By specifying the "--outdir" option, you can choose where to save the downloaded FASTQ files. 

# while IFS= read -r accession || [[ -n "$accession" ]]; do

    # Download FASTQ files for the current SRA run.
    # fastq-dl \
        # --accession "$accession" \
        # --outdir "${project_folder}/raw_data"

# done < "$ACCESSION_FILE"
    
    # Take note of the following before running the command:
    # fastq-dl needs to be installed in the rnaseq in order to work.
    # It can be installed using the command: conda install fastq-dl
    # Make sure your accession number (e.g., SRR8482616) is correct and corresponds to the dataset you want to download.
    # Check that the output directory ('${project_folder}/raw_data/') exists, to ensure all downloaded files are in the same place for easier access and analysis.
    # Sometimes fastq files can be downloaded as zip files but these can easily be uncompressed if needed or remain zipped for further downstream tools.
 
# *Method 3:* Using the 'iSeq' tool.
# Also downloads FASTQ files from the SRA or other databases.
# iSeq is more flexible because it supports various sequencing data and compression formats
# iSeq allows the specification of output formats and compression options saving space as files can be quite large.
# *Commands*:
    # - 'iSeq
    # --gzip : Compresses FASTQ files in gzip format (.fastq.gz)
    # --fastq: Ensures the output file is FASTQ which is the standard for RNA-seq data.
    # --input SRR8482616: The relevant SRA accession number for the dataset to be downloaded.
    # --database sra: Indicates that the data should be retrieved from the SRA database.
    # --output "${project_folder}/raw_data/": Specifies the output directory for the downloaded files.
#
# Take note of the following before running the command:
    # iSeq needs to be installed in the rnaseq in order to work.
    # Make sure your accession number (e.g., SRR8482616) is correct and corresponds to the dataset you want to download.
    # Check that the output directory ('${project_folder}/raw_data/') exists, to ensure all downloaded files are in the same place for easier access and analysis. 
    # If you are using the --gzip option ensure it is compatible with downstream tools.

while IFS= read -r accession || [[ -n "$accession" ]]; do

    # Download the current SRA run as a compressed FASTQ file.
    iseq \
        --gzip \
        --fastq \
        --input "$accession" \
        --database sra \
        --output "${project_folder}/raw_data"

done < "$ACCESSION_FILE"

# After downloading, assess read quality with FastQC.
# Continue with adapter trimming, post-trimming QC, STAR alignment,
# alignment quality assessment, and gene quantification with featureCounts.