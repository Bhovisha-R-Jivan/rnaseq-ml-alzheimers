#!/bin/bash - Bash is the interpreter for this script

# Project: RNA-seq and Machine Learning for Alzheimer's Disease
# Purpose: Download sequencing data for RNA-seq analysis.
# Author: Bhovisha Jivan

set -euo pipefail
# Helps catch failed commands and unset variables, making the script more robust.

# Identify the root directory of this GitHub project
PROJECT_DIR="$HOME/rnaseq-ml-alzheimers"

# Define the directory for raw sequencing data
RAW_DATA_DIR="${PROJECT_DIR}/data/raw"

# Define the directory for downloaded sequencing metadata
METADATA_DIR="${PROJECT_DIR}/data/metadata"

# Create the directories if they do not already exist
mkdir -p "$RAW_DATA_DIR" "$METADATA_DIR"

# Report the locations being used.
echo "Project directory: $PROJECT_DIR"
echo "Raw data directory: $RAW_DATA_DIR"
echo "Metadata directory: $METADATA_DIR"