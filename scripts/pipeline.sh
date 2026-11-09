
#!/bin/bash

# =====================================================================
# RNA-SEQ PIPELINE EXECUTION
# =====================================================================

set -euo pipefail

# Identify the project directory from the location of this script.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"

# Move into the project directory.
cd "${PROJECT_DIR}"

# Load pipeline configuration.
source scripts/0.config.sh


# =====================================================================
# 1. DATA ACQUISITION SETTINGS
# =====================================================================

# Select the data acquisition method:
# "local"    = Use FASTQ files already present in raw_data/.
# "iseq"     = Download individual SRA run accessions using iSeq.
# "fastqdl"  = Download individual SRA run accessions using fastq-dl.
# "prefetch" = Download and convert SRA run accessions using SRA Toolkit.

ACQUISITION_METHOD="${ACQUISITION_METHOD:-local}"

# Validate acquisition method.
if [[ "${ACQUISITION_METHOD}" != "local" &&
      "${ACQUISITION_METHOD}" != "iseq" &&
      "${ACQUISITION_METHOD}" != "fastqdl" &&
      "${ACQUISITION_METHOD}" != "prefetch" ]]; then

    echo "ERROR: Invalid acquisition method: ${ACQUISITION_METHOD}"
    echo "Allowed values: local, iseq, fastqdl, prefetch."
    exit 1

fi


# =====================================================================
# 2. DATA ACQUISITION
# =====================================================================

echo ""
echo "============================================================"
echo "STEP 1: DATA ACQUISITION"
echo "Acquisition method: ${ACQUISITION_METHOD}"
echo "============================================================"

case "${ACQUISITION_METHOD}" in

    local)

        echo "Using locally supplied FASTQ files."

        ;;

    iseq)

        echo "Downloading SRA accessions using iSeq."

        DATA_SOURCE=sra bash scripts/1a.download_data_iseq.sh

        ;;

    fastqdl)

        echo "Downloading SRA accessions using fastq-dl."

        bash scripts/1b.download_data_fastqdl.sh

        ;;

    prefetch)

        echo "Downloading and converting SRA files using SRA Toolkit."

        bash scripts/1c.download_data_prefetch.sh

        ;;

esac


# Validate FASTQ inputs after acquisition.
# Use the local input mode to check that the FASTQ files exist
# and that paired-end reads have matching R1 and R2 files.

echo ""
echo "Validating FASTQ input files."

DATA_SOURCE=local bash scripts/1a.download_data_iseq.sh


# =====================================================================
# 3. RAW READ QUALITY ASSESSMENT
# =====================================================================

echo ""
echo "============================================================"
echo "STEP 2: RAW READ QUALITY ASSESSMENT"
echo "Running FastQC and MultiQC."
echo "============================================================"

bash scripts/2.quality_assessment.sh


# =====================================================================
# 4. ADAPTER TRIMMING AND POST-TRIMMING QC
# =====================================================================

echo ""
echo "============================================================"
echo "STEP 3: ADAPTER TRIMMING"
echo "Running Trim Galore, FastQC and MultiQC."
echo "============================================================"

bash scripts/3.adapter_trimming.sh


# =====================================================================
# 5. STAR GENOME INDEX GENERATION
# =====================================================================

echo ""
echo "============================================================"
echo "STEP 4: STAR GENOME INDEX GENERATION"
echo "============================================================"

# Check whether a complete STAR index already exists.
if [[ -s "${STAR_INDEX_DIR}/Genome" &&
      -s "${STAR_INDEX_DIR}/SA" ]]; then

    echo "A STAR genome index already exists."
    echo "Skipping index generation."
    echo "Index directory: ${STAR_INDEX_DIR}"

else

    # Check whether an incomplete index directory exists.
    if [[ -d "${STAR_INDEX_DIR}" ]] &&
       [[ -n "$(ls -A "${STAR_INDEX_DIR}")" ]]; then

        echo "ERROR: An incomplete STAR index directory exists."
        echo "Check or remove it before rebuilding the index."
        exit 1

    fi

    echo "Building STAR genome index."

    bash scripts/4.indexing_genome.sh

fi


# =====================================================================
# 6. STAR ALIGNMENT
# =====================================================================

echo ""
echo "============================================================"
echo "STEP 5: STAR ALIGNMENT"
echo "============================================================"

bash scripts/5.star_alignment.sh


# =====================================================================
# 7. ALIGNMENT QUALITY ASSESSMENT
# =====================================================================

echo ""
echo "============================================================"
echo "STEP 6: ALIGNMENT QUALITY ASSESSMENT"
echo "Running Qualimap."
echo "============================================================"

bash scripts/6.alignment_qc.sh


# =====================================================================
# 8. GENE-LEVEL QUANTIFICATION
# =====================================================================

echo ""
echo "============================================================"
echo "STEP 7: GENE-LEVEL QUANTIFICATION"
echo "Running featureCounts."
echo "============================================================"

bash scripts/7.feature_counts.sh

# =====================================================================
# 9. FINAL PIPELINE QUALITY CONTROL SUMMARY
# =====================================================================

echo ""
echo "============================================================"
echo "STEP 8: FINAL MULTIQC REPORT"
echo "============================================================"

bash scripts/15.multiqc.sh

# =====================================================================
# 10. PIPELINE COMPLETION
# =====================================================================

echo ""
echo "============================================================"
echo "RNA-SEQ PIPELINE COMPLETED SUCCESSFULLY"
echo "============================================================"

echo "Project directory: ${PROJECT_DIR}"
echo "Raw FASTQ files: ${RAW_DATA_DIR}"
echo "STAR genome index: ${STAR_INDEX_DIR}"
echo "Analysis outputs: ${ANALYSIS_DIR}"
echo "Gene counts: ${ANALYSIS_DIR}/6_quantification/featureCounts/featurecounts.tsv"

echo ""
echo "Completed analysis stages:"
echo "1. Data acquisition"
echo "2. Raw read quality assessment"
echo "3. Adapter trimming and post-trimming QC"
echo "4. STAR genome indexing"
echo "5. STAR alignment"
echo "6. Alignment quality assessment"
echo "7. Gene-level quantification"

echo "============================================================"