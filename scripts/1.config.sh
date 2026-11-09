#!/bin/bash

# ===========================================================================================================================================================================
# RNA-seq PIPELINE CONFIGURATION
# ===========================================================================================================================================================================
# Configure the pipeline for a particular bulk RNA-seq dataset.
# Supported input: Gzipped FASTQ (.fastq.gz)
# Supported read types: Paired-end and single-end
# Reference: User-supplied genome FASTA and matching annotation GTF
# Supports organisms with suitable reference genome and annotation resources.
# ===========================================================================================================================================================================


# ===========================================================================================================================================================================
# 1. PROJECT LOCATION
# ===========================================================================================================================================================================

# Main project directory.
PROJECT_DIR="${PROJECT_DIR:-/mnt/c/Users/Bhovisha/OneDrive/Desktop/rnaseq-ml-alzheimers}"

# Input FASTQ directory containing the sequencing data.
RAW_DATA_DIR="${RAW_DATA_DIR:-/mnt/c/Users/Bhovisha/OneDrive/Desktop/rnaseq-ml-alzheimers/raw_data}"

# Main project directories.
REFERENCE_DIR="${PROJECT_DIR}/resource_folder"
SCRIPTS_DIR="${PROJECT_DIR}/scripts"
ANALYSIS_DIR="${PROJECT_DIR}/analysis"

# STAR genome index directory.
STAR_INDEX_DIR="${REFERENCE_DIR}/star_index"


# ===========================================================================================================================================================================
# 2. DATA ACQUISITION SETTINGS
# ===========================================================================================================================================================================

# Text file containing individual SRA run accession IDs.
ACCESSION_FILE="${PROJECT_DIR}/accession_ids.txt"

# Data acquisition method:
# "sra"   = Download SRA accessions using iSeq.
# "local" = Use FASTQ files already present in raw_data/.

DATA_SOURCE="${DATA_SOURCE:-local}"


# ===========================================================================================================================================================================
# 3. INPUT DATA SETTINGS
# ===========================================================================================================================================================================

# Sequencing read type:
# "paired" = Paired-end sequencing
# "single" = Single-end sequencing

READ_TYPE="${READ_TYPE:-paired}"

# Paired-end FASTQ naming convention.
R1_SUFFIX="${R1_SUFFIX:-_1.fastq.gz}"
R2_SUFFIX="${R2_SUFFIX:-_2.fastq.gz}"


# ===========================================================================================================================================================================
# 4. REFERENCE GENOME AND ANNOTATION
# ===========================================================================================================================================================================

# User-supplied reference genome and matching annotation.
# Both files must correspond to the same organism and genome assembly.
#
# Expected default filenames:
# resource_folder/genome.fa
# resource_folder/annotation.gtf

GENOME="${GENOME:-${REFERENCE_DIR}/genome.fa}"
ANNOTATION="${ANNOTATION:-${REFERENCE_DIR}/annotation.gtf}"


# ===========================================================================================================================================================================
# 5. COMPUTATIONAL RESOURCES
# ===========================================================================================================================================================================

THREADS="${THREADS:-2}"

# ===========================================================================================================================================================================
# 6. READ TRIMMING PARAMETERS
# ===========================================================================================================================================================================

TRIM_QUALITY="${TRIM_QUALITY:-30}"
MIN_LENGTH="${MIN_LENGTH:-20}"

# ===========================================================================================================================================================================
# 7. ALIGNMENT PARAMETERS
# ===========================================================================================================================================================================

# Maximum read length minus one.
# Default of 49 is for the current 50 bp test reads.
# Update this value for other datasets.
# Changing this parameter requires rebuilding the STAR index.

SJDB_OVERHANG="${SJDB_OVERHANG:-49}"


# ===========================================================================================================================================================================
# 8. QUANTIFICATION PARAMETERS
# ===========================================================================================================================================================================

FEATURE="${FEATURE:-exon}"
ATTRIBUTE="${ATTRIBUTE:-gene_id}"

# Library strandedness:
# 0 = Unstranded
# 1 = Forward-stranded
# 2 = Reverse-stranded

STRANDEDNESS="${STRANDEDNESS:-0}"


# ===========================================================================================================================================================================
# 9. CONFIGURATION VALIDATION
# ===========================================================================================================================================================================

# Validate data acquisition method.
if [[ "${DATA_SOURCE}" != "local" && "${DATA_SOURCE}" != "sra" ]]; then
    echo "ERROR: DATA_SOURCE must be either 'local' or 'sra'."
    return 1 2>/dev/null || exit 1
fi

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
# ===========================================================================================================================================================================