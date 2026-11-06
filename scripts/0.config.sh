#!/bin/bash

# ========================================================================================================================================================
# RNA-seq PIPELINE CONFIGURATION
# ========================================================================================================================================================
#
# Edit the settings in this file to configure the pipeline for a particular RNA-seq dataset.
#
# The individual analysis scripts will read these variables.
#
# Supported input: gzipped FASTQ files (.fastq.gz)
# Supported read types: paired-end and single-end
# Reference: User-supplied genome FASTA and annotation GTF
#
# ===========================================================================================================================================================


# ===========================================================================================================================================================
# 1. PROJECT LOCATION
# ===========================================================================================================================================================

# Default project location.
# Change this path if your Desktop is stored in a cloud folder, such as OneDrive.
# Alternatively, PROJECT_DIR can be set before running the pipeline to override this default.

PROJECT_DIR="${PROJECT_DIR:-$HOME/Desktop/RNAseq_Pipeline}"

# Main project directories

RAW_DATA_DIR="${PROJECT_DIR}/raw_data"
REFERENCE_DIR="${PROJECT_DIR}/reference"
SCRIPTS_DIR="${PROJECT_DIR}/scripts"
ANALYSIS_DIR="${PROJECT_DIR}/analysis"

# STAR index directory

STAR_INDEX_DIR="${REFERENCE_DIR}/star_index"

# ===========================================================================================================================================================
# 2. DATA ACQUISITION SETTINGS
# ===========================================================================================================================================================

# Text file containing SRA run accession IDs.
ACCESSION_FILE="${PROJECT_DIR}/accession_ids.txt"

# ===========================================================================================================================================================
# 2. INPUT DATA SETTINGS
# ===========================================================================================================================================================

# Specify the sequencing read type:
# "paired"  = Paired-end sequencing
# "single"  = Single-end sequencing

READ_TYPE="paired"

# FASTQ file naming conventions for paired-end data.
# Example:
# sample_1.fastq.gz
# sample_2.fastq.gz

R1_SUFFIX="_1.fastq.gz"
R2_SUFFIX="_2.fastq.gz"

# Expected FASTQ file format: .fastq.gz

# =============================================================================================================================================================
# 3. REFERENCE GENOME AND ANNOTATION
# =============================================================================================================================================================

# User-supplied reference genome and matching annotation.
# These files must belong to the same organism and genome assembly.

GENOME="${REFERENCE_DIR}/genome.fa"
ANNOTATION="${REFERENCE_DIR}/annotation.gtf"

# ============================================================================================================================================================
# 4. COMPUTATIONAL RESOURCES
# ============================================================================================================================================================

# Number of CPU threads allocated to analysis tools.
THREADS=2

# ============================================================================================================================================================
# 5. READ TRIMMING PARAMETERS
# ============================================================================================================================================================  

# Minimum Phred quality threshold for trimming.
TRIM_QUALITY=30

# Minimum read length retained after trimming.
MIN_LENGTH=20

# ============================================================================================================================================================
# 6. ALIGNMENT PARAMETERS
# ============================================================================================================================================================

# STAR genome index overhang.
#
# Typically calculated as:
# Maximum read length - 1
#
# This value will need to be adjusted according to the sequencing read length.

SJDB_OVERHANG=49

# ============================================================================================================================================================
# 7. QUANTIFICATION PARAMETERS
# ============================================================================================================================================================

# featureCounts feature type.
# "exon" selects exon features from the GTF annotation.

FEATURE="exon"

# Group exon counts by gene identifier.

ATTRIBUTE="gene_id"

# Library strandedness:
# 0 = Unstranded
# 1 = Forward-stranded
# 2 = Reverse-stranded

STRANDEDNESS=0

# =============================================================================================================================================================
# END OF CONFIGURATION