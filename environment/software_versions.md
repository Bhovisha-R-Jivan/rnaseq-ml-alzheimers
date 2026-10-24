# Software Environment

## Overview

This document records the software and package versions used in the RNA-seq processing, statistical analysis, and machine learning workflows. Recording these versions supports reproducibility and helps identify differences between computational environments.

The project uses a local Ubuntu/WSL Conda environment for RNA-seq processing and R-based statistical analysis, while Python-based data analysis and machine learning are performed using Google Colab.

## RNA-seq Processing

The following software is installed in the local `rnaseq` Conda environment.

| Software      | Version | Purpose                                                |
| ------------- | ------- | ------------------------------------------------------ |
| STAR          | 2.7.11b | Alignment of RNA-seq reads to the reference genome     |
| FastQC        | 0.12.1  | Quality assessment of raw and trimmed sequencing reads |
| MultiQC       | 1.34    | Aggregation and reporting of quality control results   |
| Trim Galore   | 0.6.11  | Adapter removal and quality trimming                   |
| SAMtools      | 1.23.1  | Processing and manipulation of alignment files         |
| featureCounts | 2.1.1   | Gene-level read quantification                         |

## Statistical Analysis

R and DESeq2 are installed in the local `rnaseq` Conda environment.

| Software | Version | Purpose                                 |
| -------- | ------- | --------------------------------------- |
| R        | 4.4.3   | Statistical computing and data analysis |
| DESeq2   | 1.46.0  | Differential gene expression analysis   |

## Python and Machine Learning

Python-based data analysis and machine learning will be conducted using Google Colab.

| Software     | Version        | Purpose                                |
| ------------ | -------------- | -------------------------------------- |
| Python       | To be recorded | Programming language for data analysis |
| pandas       | To be recorded | Data manipulation and analysis         |
| NumPy        | To be recorded | Numerical computing                    |
| scikit-learn | To be recorded | Machine learning models and evaluation |

Versions will be recorded from the Google Colab runtime used for the analysis.

## Computational Environments

* **RNA-seq processing:** Ubuntu/WSL using the local Conda environment named `rnaseq`.
* **Statistical analysis:** R and DESeq2 within the local `rnaseq` environment.
* **Python and machine learning:** Google Colab.

