#!/bin/bash
# #! = shebang, which specifies the interpreter used to execute the script.
# /bin/bash = location of the Bash interpreter.

set -euo pipefail
# -e = terminates the script if a command fails (returns a non-zero exit status).
# -u = treats the use of an undefined variable as an error.
# -o pipefail = causes a pipeline to return a failure if any command within it fails.

# Define project directories.
# Double quotes are used to protect paths containing spaces.
RAW_DATA="$HOME/testing_code/raw_data"
QC_DIR="$HOME/testing_code/analysis/1_fastqc"

# Create the output directory if it does not exist.
# mkdir = creates a directory.
# -p = creates parent directories as needed and does not produce an error if the directory already exists.
mkdir -p "$QC_DIR"

# Activate the Conda environment.
# source = executes a file within the current shell.
# conda.sh = configures the shell so that Conda commands can be used.
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate rnaseq

# Check that FastQC is installed.
# command -v fastqc = checks whether the FastQC command is available in the current environment.
# ! = logical NOT, reverses the command's exit status.
# >/dev/null = redirects standard output to /dev/null, where it is discarded.
# 2>&1 = redirects standard error (2) to standard output (1), so both are discarded.
# then = marks the beginning of the commands executed if the condition is true.
# exit 1 = terminates the script with a non-zero exit status, indicating an error.
if ! command -v fastqc >/dev/null 2>&1; then
    echo "ERROR: FastQC is not installed in the rnaseq environment."
    exit 1
fi

# Set the number of threads.
# Threads allow a programme to perform multiple tasks simultaneously.
N_THREADS=2

# Identify forward-read FASTQ files.
# shopt = Bash command used to configure shell options.
# -s = enables the specified option.
# nullglob = causes unmatched filename patterns to expand to nothing instead of remaining as literal text.
# * = wildcard that matches any sequence of characters.
# _1.fastq.gz = restricts the match to forward-read FASTQ files.
shopt -s nullglob
FASTQ_FILES=("$RAW_DATA"/*_1.fastq.gz)

# Check that forward-read files exist.
# [[ ]] = Bash conditional-testing syntax.
# ${} = expands a variable or parameter.
# # = retrieves the number of characters or elements, depending on the variable.
# array = variable that can hold multiple values.
# FASTQ_FILES = the array containing the forward-read filenames.
# [@] = refers to all elements of the array.
# Together, ${#FASTQ_FILES[@]} gives the number of elements in the array.
# -eq = compares integers to check whether they are equal.
# 0 = checks whether the array contains zero files.
if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No forward-read FASTQ files found in $RAW_DATA"
    exit 1
fi

# Validate that every forward read has a corresponding reverse read.
# for = starts a loop.
# INPUT_FASTQ1 = variable holding the current forward-read filename.
# in = specifies the collection over which the loop iterates.
# "${FASTQ_FILES[@]}" = expands all elements of the array.
# do = marks the beginning of the loop body.
for INPUT_FASTQ1 in "${FASTQ_FILES[@]}"; do

    # Construct the corresponding reverse-read filename.
    # % = removes the shortest matching pattern from the end of the variable's value.
    # _1.fastq.gz = suffix being removed.
    # _2.fastq.gz = suffix being appended.
    INPUT_FASTQ2="${INPUT_FASTQ1%_1.fastq.gz}_2.fastq.gz"

    # Check that both FASTQ files exist.
    # ! = logical NOT.
    # -f = checks whether the specified path is an existing regular file.
    # || = logical OR, meaning the condition is true if either test is true.
    # If either FASTQ file is missing, an error message is printed and the script terminates.
    if [[ ! -f "$INPUT_FASTQ1" || ! -f "$INPUT_FASTQ2" ]]; then
        echo "ERROR: Missing paired FASTQ file for ${INPUT_FASTQ1}"
        exit 1
    fi

    # echo = prints a confirmation message.
    echo "Both files exist: $INPUT_FASTQ1"

done
# done = marks the end of the for loop.

# Run FastQC on all paired-end FASTQ files.
# A second loop is used because the first loop only validates the input files.

for INPUT_FASTQ1 in "${FASTQ_FILES[@]}"; do

    # Construct the corresponding reverse-read filename.
    INPUT_FASTQ2="${INPUT_FASTQ1%_1.fastq.gz}_2.fastq.gz"

    # Run FastQC on both reads.
    # fastqc = executes the FastQC programme.
    # \ = indicates that the command continues onto the next line.
    # --outdir = specifies the directory where FastQC reports will be saved.
    # --threads = specifies the number of processing threads to use.
    # "$INPUT_FASTQ1" and "$INPUT_FASTQ2" = the forward and reverse FASTQ input files.
    # FastQC processes both files for each iteration of the loop.
    fastqc \
        --outdir "$QC_DIR" \
        --threads "$N_THREADS" \
        "$INPUT_FASTQ1" "$INPUT_FASTQ2"

done

# Print a completion message after all samples have been processed.
echo "FastQC analysis completed successfully."

######################################################################################################################################################

# Trimming

#!/bin/bash
# Shebang specifies Bash as the interpreter.

set -euo pipefail
# -e = stops the script if a command fails.
# -u = treats undefined variables as errors.
# -o pipefail = makes a pipeline fail if any command within it fails.

# Define project directories.
RAW_DATA="$HOME/testing_code/raw_data"
QC_DIR="$HOME/testing_code/analysis/1_fastqc"
TRIMMED_DIR="$HOME/testing_code/analysis/2_trimmed"
POSTQC_DIR="$HOME/testing_code/analysis/3_postqc_fastqc"

# Create the output directories if they do not already exist.
# mkdir = creates a directory.
# -p = creates parent directories if needed and does not produce an error if the directory already exists.
mkdir -p "$TRIMMED_DIR"
mkdir -p "$POSTQC_DIR"

# Activate the Conda environment.
# source = executes a file within the current shell.
# conda.sh = configures the shell so Conda commands can be used.
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate rnaseq

# Check that required tools are installed.
# command -v = checks whether a programme is available.
# ! = reverses the command's result.
# >/dev/null 2>&1 = discards normal output and error messages.
# then = begins the commands executed if the condition is true.
# exit 1 = stops the script with an error status.
# fi = closes the conditional statement.

if ! command -v trim_galore >/dev/null 2>&1; then
    echo "ERROR: Trim Galore is not installed in the rnaseq environment."
    exit 1
fi

if ! command -v fastqc >/dev/null 2>&1; then
    echo "ERROR: FastQC is not installed in the rnaseq environment."
    exit 1
fi

# Check that the raw data directory exists.
# [[ ]] = checks a condition.
# ! = NOT.
# -d = checks whether a directory exists.
if [[ ! -d "$RAW_DATA" ]]; then
    echo "ERROR: Raw data directory not found: $RAW_DATA"
    exit 1
fi

# Identify forward-read FASTQ files.
# shopt = changes how Bash behaves.
# -s = enables an option.
# nullglob = makes unmatched wildcards expand to nothing.
# FASTQ_FILES = an array storing the matching filenames.
# () = creates the array.
# * = matches any sequence of characters.
# _1.fastq.gz = restricts the search to forward-read FASTQ files.
shopt -s nullglob
FASTQ_FILES=("$RAW_DATA"/*_1.fastq.gz)

# Check that forward-read files exist.
# ${} = accesses a variable.
# # = counts the elements in the array in this expression.
# [@] = refers to all elements of the array.
# -eq = checks whether two numbers are equal.
# 0 = checks whether the array contains zero files.
if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No forward-read FASTQ files found in $RAW_DATA"
    exit 1
fi

# Validate that every forward read has a corresponding reverse read.
# for = starts a loop.
# INPUT_FASTQ1 = temporarily stores the current forward-read filename.
# in = specifies the collection being processed.
# "${FASTQ_FILES[@]}" = supplies all filenames in the array.
# do = starts the commands inside the loop.
for INPUT_FASTQ1 in "${FASTQ_FILES[@]}"; do

    # Construct the corresponding reverse-read filename.
    # ${INPUT_FASTQ1} = accesses the forward-read filename.
    # % = removes the specified ending from the filename.
    # _1.fastq.gz = the ending being removed.
    # _2.fastq.gz = the ending being added.
    INPUT_FASTQ2="${INPUT_FASTQ1%_1.fastq.gz}_2.fastq.gz"

    # Check that both FASTQ files exist.
    # -f = checks whether a file exists.
    # || = logical OR; the condition is true if either test is true.
    if [[ ! -f "$INPUT_FASTQ1" || ! -f "$INPUT_FASTQ2" ]]; then
        echo "ERROR: Missing paired FASTQ file for ${INPUT_FASTQ1}"
        exit 1
    fi

    echo "Both files exist: $INPUT_FASTQ1"

done
# done = closes the validation loop.

# Set the number of processing threads.
N_THREADS=2

# Trim adapters and low-quality bases, then assess the trimmed reads with FastQC.
# This second loop performs the actual trimming after the first loop validates the files.
for INPUT_FASTQ1 in "${FASTQ_FILES[@]}"; do

    # Construct the corresponding reverse-read filename.
    INPUT_FASTQ2="${INPUT_FASTQ1%_1.fastq.gz}_2.fastq.gz"

    # Run Trim Galore on both paired-end FASTQ files.
    # trim_galore = runs the Trim Galore programme.
    # \ = allows the command to continue onto the next line.
    # --cores = specifies the number of processing cores.
    # --paired = tells Trim Galore that the input consists of paired-end reads.
    # --quality 30 = sets the Phred quality threshold for trimming.
    # --length 20 = discards reads shorter than 20 bases after trimming.
    # --gzip = compresses the resulting FASTQ files.
    # --output_dir = specifies where trimmed FASTQ files are saved.
    # --fastqc = runs FastQC on the trimmed reads.
    # --fastqc_args = passes additional options to FastQC.
    # "$INPUT_FASTQ1" and "$INPUT_FASTQ2" = the forward and reverse input files.

    trim_galore \
        --cores "$N_THREADS" \
        --paired \
        --quality 30 \
        --length 20 \
        --gzip \
        --output_dir "$TRIMMED_DIR" \
        --fastqc \
        --fastqc_args "--outdir $POSTQC_DIR" \
        "$INPUT_FASTQ1" "$INPUT_FASTQ2"

done
# done = closes the trimming loop.

# Print a completion message.
echo "Adapter trimming completed successfully."

##############################################################################################################################################################
# Star index

# Specifies that this script should be run using Bash.
#!/bin/bash

# Stop the script if a command fails, an undefined variable is used,
# or a command in a pipeline fails.
set -euo pipefail

# Define the main project directory.
PROJECT_DIR="$HOME/testing_code"

# Define the directory containing reference files and indexing resources.
RESOURCE_DIR="${PROJECT_DIR}/resource_folder"

# Define the directory containing the reference genome files.
REFERENCE_DIR="${RESOURCE_DIR}/reference_genome"

# Define where STAR will save the generated genome index.
STAR_INDEX_DIR="${RESOURCE_DIR}/star_index"

# Define the reference genome FASTA file.
GENOME_FASTA="${REFERENCE_DIR}/GRCh38.primary_assembly.genome.fa"

# Define the gene annotation GTF file.
ANNOTATION_GTF="${REFERENCE_DIR}/gencode.v47.primary_assembly.annotation.gtf"

# Set the number of CPU threads STAR can use.
N_THREADS=2

# Set the STAR splice-junction overhang parameter.
# Has a usual recommendation of maximum read length - 1
# A value of 49 corresponds to a read length of 50 bases.
SJDB_OVERHANG=49

# Create the reference genome and STAR index directories if they do not exist.
mkdir -p "$REFERENCE_DIR"
mkdir -p "$STAR_INDEX_DIR"

# Load Conda so that it can be used in this Bash session.
source "$HOME/miniconda3/etc/profile.d/conda.sh"

# Activate the environment containing the required RNA-seq tools.
conda activate rnaseq

# Check whether the reference genome FASTA file exists.
# -f = checks if a file exists
if [[ ! -f "$GENOME_FASTA" ]]; then

    # Display an error if the FASTA file cannot be found.
    echo "ERROR: Reference genome FASTA not found: $GENOME_FASTA"

    # Stop the script because the reference genome is missing.
    exit 1
fi

# Check whether the annotation GTF file exists.
if [[ ! -f "$ANNOTATION_GTF" ]]; then

    # Display an error if the GTF file cannot be found.
    echo "ERROR: GENCODE annotation GTF not found: $ANNOTATION_GTF"

    # Stop the script because the annotation file is missing.
    exit 1
fi

# Check whether SAMtools is installed and available.
if ! command -v samtools >/dev/null 2>&1; then

    # Display an error if SAMtools cannot be found.
    echo "ERROR: SAMtools is not installed in the rnaseq environment."

    # Stop the script because SAMtools is required.
    exit 1
fi

# Check whether STAR is installed and available.
if ! command -v STAR >/dev/null 2>&1; then

    # Display an error if STAR cannot be found.
    echo "ERROR: STAR is not installed in the rnaseq environment."

    # Stop the script because STAR is required.
    exit 1
fi

# Generate an index for the reference genome FASTA file.
# This allows tools such as SAMtools to access specific sequences efficiently.
samtools faidx "$GENOME_FASTA"

# Generate the STAR genome index using the reference genome and annotation.
STAR \
    --runMode genomeGenerate \ # Tells star to generate a genome index instead of align
    --runThreadN "$N_THREADS" \
    --genomeDir "$STAR_INDEX_DIR" \
    --genomeFastaFiles "$GENOME_FASTA" \
    --sjdbGTFfile "$ANNOTATION_GTF" \
    --sjdbOverhang "$SJDB_OVERHANG"

# Display a message confirming that indexing has completed.
echo "Reference genome indexing completed successfully."

####################################################################################################################################
#!/bin/bash

# Stop the script if a command fails, an undefined variable is used,
# or a command in a pipeline fails.
set -euo pipefail

# Define the main project directory.
PROJECT_DIR="$HOME/testing_code"

# Define the directory containing the STAR genome index.
STAR_INDEX_DIR="${PROJECT_DIR}/resource_folder/star_index"

# Define the directory containing trimmed FASTQ files.
TRIMMED_DIR="${PROJECT_DIR}/analysis/2_trimmed"

# Define where the resulting BAM files will be saved.
BAM_DIR="${PROJECT_DIR}/analysis/4_bamfiles"

# Set the number of processing threads.
N_THREADS=2

# Create the BAM output directory if it does not already exist.
mkdir -p "$BAM_DIR"

# Load Conda so that it can be used in this Bash session.
source "$HOME/miniconda3/etc/profile.d/conda.sh"

# Activate the environment containing STAR and SAMtools.
conda activate rnaseq

# Check whether the STAR genome index exists.
if [[ ! -f "${STAR_INDEX_DIR}/Genome" ]]; then
    echo "ERROR: STAR genome index not found: $STAR_INDEX_DIR"
    exit 1
fi

# Check whether STAR is installed and available.
if ! command -v STAR >/dev/null 2>&1; then
    echo "ERROR: STAR is not installed in the rnaseq environment."
    exit 1
fi

# Check whether SAMtools is installed and available.
if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: SAMtools is not installed in the rnaseq environment."
    exit 1
fi

# Allow unmatched wildcards to expand to nothing.
shopt -s nullglob

# Store all trimmed forward-read FASTQ filenames in an array.
TRIMMED_READ1_FILES=("$TRIMMED_DIR"/*_1_val_1.fq.gz)

# Check whether any trimmed forward-read files were found.
# ${#TRIMMED_READ1_FILES[@]} = counts the number of elements in the array.
# == 0 = the number must equal to zero for the condition to be true
if (( ${#TRIMMED_READ1_FILES[@]} == 0 )); then
    echo "ERROR: No trimmed FASTQ files found in $TRIMMED_DIR"
    exit 1
fi

# Check that every forward-read file has a matching reverse-read file.
# STAR needs both FASTQ files for paired-end alignment, this check catches missing mates before alignment
for INPUT_FASTQ1 in "${TRIMMED_READ1_FILES[@]}"; do

    # Construct the corresponding reverse-read filename.
    INPUT_FASTQ2="${INPUT_FASTQ1%_1_val_1.fq.gz}_2_val_2.fq.gz"

    # Stop if the matching reverse-read file does not exist.
    if [[ ! -f "$INPUT_FASTQ2" ]]; then
        echo "ERROR: Matching reverse-read file not found: $INPUT_FASTQ2"
        exit 1
    fi

done

# Align paired-end reads using STAR.
for INPUT_FASTQ1 in "${TRIMMED_READ1_FILES[@]}"; do

    # Construct the corresponding reverse-read filename.
    INPUT_FASTQ2="${INPUT_FASTQ1%_1_val_1.fq.gz}_2_val_2.fq.gz"

    # Extract the sample name from the FASTQ filename.
    SAMPLE_NAME="$(basename "${INPUT_FASTQ1%_1_val_1.fq.gz}")"

    # Define the STAR output filename prefix.
    OUT_PREFIX="${BAM_DIR}/${SAMPLE_NAME}_"

    # Define the name of the final sorted BAM file.
    SORTED_BAM="${BAM_DIR}/${SAMPLE_NAME}_Aligned.sorted.bam"

    echo "Aligning sample: ${SAMPLE_NAME}"

    # Align paired-end reads against the reference genome.
    STAR \
        --genomeDir "$STAR_INDEX_DIR" \
        --readFilesIn "$INPUT_FASTQ1" "$INPUT_FASTQ2" \
        --readFilesCommand zcat \ # decompresses gzipped FASTQ files while STAR reads them
        --outSAMunmapped Within \ # includes unmapped reads in the output BAM file
        --outFileNamePrefix "$OUT_PREFIX" \
        --runThreadN "$N_THREADS"

    # Sort the alignments and convert the SAM file into BAM format.
    samtools sort \
        --threads "$N_THREADS" \
        -O BAM \
        -o "$SORTED_BAM" \
        "${OUT_PREFIX}Aligned.out.sam"

    # Generate an index for the sorted BAM file.
    samtools index \
        --threads "$N_THREADS" \
        "$SORTED_BAM"

    # Remove the intermediate SAM file after successful conversion.
    rm "${OUT_PREFIX}Aligned.out.sam"

    echo "Completed alignment for sample: ${SAMPLE_NAME}"

done

# Display a message confirming that the alignment workflow is complete.
echo "STAR alignment workflow completed successfully."

##############################################################################################################################################################################
#!/bin/bash

# Stop the script if a command fails, an undefined variable is used,
# or a command in a pipeline fails.
set -euo pipefail

# Define the main project directory.
PROJECT_DIR="$HOME/testing_code"

# Define the directory containing reference resources.
RESOURCE_DIR="${PROJECT_DIR}/resource_folder"

# Define the directory containing sorted BAM files.
BAM_DIR="${PROJECT_DIR}/analysis/4_bamfiles"

# Define where alignment QC reports will be saved.
QC_DIR="${PROJECT_DIR}/analysis/5_alignmentQC"

# Define the chromosome X gene annotation file.
ANNOTATION_GTF="${RESOURCE_DIR}/GRCh38_chrX_gencode.v47.annotation.gtf"

# Set the number of processing threads.
N_THREADS=2

# Create the output directory if it does not already exist.
mkdir -p "$QC_DIR"

# Load Conda so that it can be used in this Bash session.
source "$HOME/miniconda3/etc/profile.d/conda.sh"

# Activate the environment containing Qualimap.
conda activate rnaseq

# Check whether the annotation GTF file exists.
# -f = checks if a file exists
# -d = checks if a directory exists
# ! = reverses the condition
if [[ ! -f "$ANNOTATION_GTF" ]]; then
    echo "ERROR: GENCODE annotation file not found: $ANNOTATION_GTF"
    exit 1
fi

# Check whether Qualimap is installed and available.
if ! command -v qualimap >/dev/null 2>&1; then
    echo "ERROR: Qualimap is not installed in the rnaseq environment."
    exit 1
fi

# Check whether the BAM directory exists.
if [[ ! -d "$BAM_DIR" ]]; then
    echo "ERROR: BAM directory not found: $BAM_DIR"
    exit 1
fi

# Allow unmatched wildcards to expand to nothing.
shopt -s nullglob

# Store all sorted BAM filenames in an array.
BAM_FILES=("$BAM_DIR"/*_Aligned.sorted.bam)

# Check whether any sorted BAM files were found.
if (( ${#BAM_FILES[@]} == 0 )); then
    echo "ERROR: No sorted BAM files found in $BAM_DIR"
    exit 1
fi

# Run alignment quality control for each sample.
for BAM_FILE in "${BAM_FILES[@]}"; do

    # Extract the sample name from the BAM filename.
    SAMPLE_NAME="$(basename "${BAM_FILE%_Aligned.sorted.bam}")"

    # Define a separate output directory for this sample.
    SAMPLE_QC_DIR="${QC_DIR}/${SAMPLE_NAME}"

    # Check whether the BAM index file exists.
    if [[ ! -f "${BAM_FILE}.bai" ]]; then
        echo "ERROR: BAM index not found for ${SAMPLE_NAME}"
        exit 1
    fi

    # Create the sample-specific output directory.
    mkdir -p "$SAMPLE_QC_DIR"

    echo "Running Qualimap for sample: ${SAMPLE_NAME}"

    # Assess the quality of the RNA-seq alignments.
    qualimap rnaseq \
        -bam "$BAM_FILE" \
        -gtf "$ANNOTATION_GTF" \
        -outdir "$SAMPLE_QC_DIR" \
        -outfile "${SAMPLE_NAME}" \
        -outformat HTML \
        -sorted

    echo "Completed alignment QC for sample: ${SAMPLE_NAME}"

done

# Display a message confirming that the QC workflow is complete.
echo "Alignment quality-control workflow completed successfully."

###########################################################################################################################################
# Feature Counts
#!/bin/bash

# Stop the script if a command fails, an undefined variable is used,
# or a command in a pipeline fails.
set -euo pipefail

# Define the main project directory.
PROJECT_DIR="$HOME/testing_code"

# Define the directory containing reference resources.
RESOURCE_DIR="${PROJECT_DIR}/resource_folder"

# Define the directory containing sorted BAM files.
BAM_DIR="${PROJECT_DIR}/analysis/4_bamfiles"

# Define where gene-counting results will be saved.
QUANT_DIR="${PROJECT_DIR}/analysis/6_quantification/featureCounts"

# Define the reference gene annotation file.
ANNOTATION_GTF="${RESOURCE_DIR}/GRCh38_chrX_gencode.v47.annotation.gtf"

# Define the filename for the original featureCounts output.
COUNTS_FILE="${QUANT_DIR}/featurecounts.txt"

# Define the filename for the cleaned count table.
COUNTS_TSV="${QUANT_DIR}/featurecounts.tsv"

# Set the featureCounts parameters.
N_THREADS=2
FEATURE="exon" # count read overlapping annotated exons
ATTRIBUTE="gene_id" # groups exons belonging to the same gene
STRAND=0  # Provisional: unstranded counting.

# Create the output directory if it does not already exist.
mkdir -p "$QUANT_DIR"

# Load Conda so that it can be used in this Bash session.
source "$HOME/miniconda3/etc/profile.d/conda.sh"

# Activate the environment containing featureCounts.
conda activate rnaseq

# Check whether the annotation GTF file exists.
if [[ ! -f "$ANNOTATION_GTF" ]]; then
    echo "ERROR: GTF annotation file not found: $ANNOTATION_GTF"
    exit 1
fi

# Check whether the BAM directory exists.
if [[ ! -d "$BAM_DIR" ]]; then
    echo "ERROR: BAM directory not found: $BAM_DIR"
    exit 1
fi

# Allow unmatched wildcards to expand to nothing.
shopt -s nullglob

# Store all sorted BAM filenames in an array.
BAM_FILES=("${BAM_DIR}"/*_Aligned.sorted.bam)

# Check whether any sorted BAM files were found.
if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No sorted BAM files found in $BAM_DIR"
    exit 1
fi

# Check whether featureCounts is installed and available.
if ! command -v featureCounts >/dev/null 2>&1; then
    echo "ERROR: featureCounts is not installed in the rnaseq environment."
    exit 1
fi

# Display the number of BAM files found.
echo "Found ${#BAM_FILES[@]} BAM files."

echo "Running featureCounts..."

# Count paired-end reads assigned to genes.
featureCounts \
    -T "$N_THREADS" \
    -p \ # paired-end input
    --countReadPairs \ # counts read pairs rather than treating each mate as an independent read
    -t "$FEATURE" \
    -g "$ATTRIBUTE" \
    -s "$STRAND" \
    -a "$ANNOTATION_GTF" \
    -o "$COUNTS_FILE" \
    "${BAM_FILES[@]}"

# Check whether featureCounts produced a non-empty output file.
# if there is no usable output, script stops
if [[ ! -s "$COUNTS_FILE" ]]; then
    echo "ERROR: featureCounts output is missing or empty: $COUNTS_FILE"
    exit 1
fi

# Remove comment lines and simplify sample column names.
awk ' 
    BEGIN { FS = OFS = "\t" } #AWK field separator and output field separator set to tab
    /^#/ { next } # identify lines beginning with # so it is skipped and moved to the next one
    {
        for (i = 1; i <= NF; i++) { #i=1 starts at the first column and continues until the last column going one at a time
            sub(/^.*\//, "", $i) #removes the directory path from a column
            sub(/_Aligned\.sorted\.bam$/, "", $i) #removes _Aligned.sorted.bam ending
        }
        print
    }
' "$COUNTS_FILE" > "$COUNTS_TSV"

# Display completion messages and output file locations.
echo "featureCounts completed successfully."
echo "Raw output: $COUNTS_FILE"
echo "Processed count table: $COUNTS_TSV"
