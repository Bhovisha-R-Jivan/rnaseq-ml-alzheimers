#!/bin/bash

# ===========================================================================================================================================================================
# RNA-seq PIPELINE CONFIGURATION
# ===========================================================================================================================================================================
# Configure the pipeline for a particular bulk RNA-seq dataset.
# Supported input: Gzipped FASTQ (.fastq.gz)
# Supported read types: Paired-end and single-end
# Reference: User-supplied genome FASTA and matching annotation GTF
# ===========================================================================================================================================================================
# 1. PROJECT LOCATION
# ===========================================================================================================================================================================
# Main project directory.

# This is the OneDrive Desktop location accessed through WSL.
PROJECT_DIR="${PROJECT_DIR:-/mnt/c/Users/Bhovisha/OneDrive/Desktop/rnaseq-ml-alzheimers}"

# Main project directories.
RAW_DATA_DIR="${PROJECT_DIR}/raw_data"
REFERENCE_DIR="${PROJECT_DIR}/reference"
SCRIPTS_DIR="${PROJECT_DIR}/scripts"
ANALYSIS_DIR="${PROJECT_DIR}/analysis"

# STAR genome index directory.
STAR_INDEX_DIR="${REFERENCE_DIR}/star_index"

# ===========================================================================================================================================================================
# 2. DATA ACQUISITION SETTINGS
# ===========================================================================================================================================================================
# Text file containing individual SRA run accession IDs.
ACCESSION_FILE="${PROJECT_DIR}/accession_ids.txt"

# ===========================================================================================================================================================================
# 3. INPUT DATA SETTINGS
# ===========================================================================================================================================================================
# Sequencing read type:
# "paired" = Paired-end sequencing
# "single" = Single-end sequencing
READ_TYPE="paired"
# Paired-end FASTQ naming convention.

# Example:
# sample_1.fastq.gz
# sample_2.fastq.gz
R1_SUFFIX="_1.fastq.gz"
R2_SUFFIX="_2.fastq.gz"

# Expected input format: .fastq.gz
# ===========================================================================================================================================================================
# 4. REFERENCE GENOME AND ANNOTATION
# ===========================================================================================================================================================================
# User-supplied reference genome and matching annotation.
# Both files must correspond to the same organism and genome assembly.
GENOME="${REFERENCE_DIR}/genome.fa"
ANNOTATION="${REFERENCE_DIR}/annotation.gtf"

# ============================================================================================================================================================================
# 5. COMPUTATIONAL RESOURCES
# ============================================================================================================================================================================
# Number of CPU threads allocated to analysis tools.
THREADS=2

# ===========================================================================================================================================================================
# 6. READ TRIMMING PARAMETERS
# ===========================================================================================================================================================================
# Minimum Phred quality threshold.
TRIM_QUALITY=30
# Minimum retained read length.
MIN_LENGTH=20

# ===========================================================================================================================================================================
# 7. ALIGNMENT PARAMETERS
# ===========================================================================================================================================================================
# STAR genome index overhang.
# Typically:
# Maximum read length - 1
# This must be selected according to the sequencing read length.
SJDB_OVERHANG=49

# ===========================================================================================================================================================================
# 8. QUANTIFICATION PARAMETERS
# ===========================================================================================================================================================================
# GTF feature type.
FEATURE="exon"
# GTF attribute used to group features by gene.
ATTRIBUTE="gene_id"
# Library strandedness:
# 0 = Unstranded
# 1 = Forward-stranded
# 2 = Reverse-stranded
STRANDEDNESS=0

# ===========================================================================================================================================================================
# 9. CONFIGURATION VALIDATION
# ===========================================================================================================================================================================
# Validate sequencing layout.
if [[ "${READ_TYPE}" != "paired" && "${READ_TYPE}" != "single" ]]; then
echo "ERROR: READ_TYPE must be either 'paired' or 'single'."
return 1 2>/dev/null || exit 1
fi

# Validate computational resources.
if ! [[ "${THREADS}" =~ ^[1-9][0-9]*$ ]]; then
echo "ERROR: THREADS must be a positive integer."
return 1 2>/dev/null || exit 1
fi

# Validate trimming parameters.
if ! [[ "${TRIM_QUALITY}" =~ ^[0-9]+$ ]] ||
! [[ "${MIN_LENGTH}" =~ ^[1-9][0-9]*$ ]]; then
echo "ERROR: Invalid trimming parameters."
return 1 2>/dev/null || exit 1
fi

# Validate STAR overhang.
if ! [[ "${SJDB_OVERHANG}" =~ ^[0-9]+$ ]]; then
echo "ERROR: SJDB_OVERHANG must be a non-negative integer."
return 1 2>/dev/null || exit 1
fi

# Validate featureCounts strandedness.
if [[ "${STRANDEDNESS}" != "0" &&
"${STRANDEDNESS}" != "1" &&
"${STRANDEDNESS}" != "2" ]]; then
echo "ERROR: STRANDEDNESS must be 0, 1 or 2."
return 1 2>/dev/null || exit 1
fi

# ===========================================================================================================================================================================
# END OF CONFIGURATION