# Software Environment

## Overview

This document records the software and package versions installed for the RNA-seq processing and statistical analysis workflows. Recording these versions supports reproducibility and helps identify differences between computational environments.

RNA-seq processing and R-based statistical analysis are conducted in a local Ubuntu/WSL Conda environment named `rnaseq`. Python-based data analysis and machine learning are planned for Google Colab.

The file `environment/rnaseq_environment.yml` contains the exported Conda environment configuration, including package dependencies and versions.

## RNA-seq Processing

The following tools are installed in the local `rnaseq` Conda environment.

| Software                | Version | Purpose                                                    |
| ----------------------- | ------- | ---------------------------------------------------------- |
| SRA Toolkit             | 3.1.1   | Downloading sequencing data from the Sequence Read Archive |
| FastQC                  | 0.12.1  | Quality assessment of sequencing reads                     |
| MultiQC                 | 1.34    | Aggregating quality control reports                        |
| Trim Galore             | 0.6.11  | Adapter removal and quality trimming                       |
| HISAT2                  | 2.2.2   | Alignment of RNA-seq reads to a reference genome           |
| STAR                    | 2.7.11b | Alignment of RNA-seq reads to a reference genome           |
| SAMtools                | 1.23.1  | Processing and manipulation of alignment files             |
| Subread / featureCounts | 2.1.1   | Gene-level read quantification                             |
| HTSeq                   | 2.1.2   | Read counting for gene expression analysis                 |
| Qualimap                | 2.3     | Quality assessment of sequencing alignments                |
| RSeQC                   | 5.0.4   | RNA-seq quality control and alignment assessment           |
| Salmon                  | 1.11.4  | Transcript abundance quantification                        |
| Kallisto                | 0.52.0  | Transcript abundance quantification                        |
| fastq-dl                | 3.0.1   | Downloading sequencing data                                |
| iseq                    | 1.9.8   | Supporting sequencing-data workflows                       |

These tools are installed in the environment but are not necessarily all used in the final analysis. The primary processing workflow uses the tools required by the project's RNA-seq scripts.

## Statistical Analysis

R and the Bioconductor package DESeq2 are installed in the local `rnaseq` Conda environment.

| Software        | Version | Purpose                                 |
| --------------- | ------- | --------------------------------------- |
| R               | 4.4.3   | Statistical computing and data analysis |
| DESeq2          | 1.46.0  | Differential gene expression analysis   |
| clusterProfiler | 4.14.0  | Functional enrichment analysis          |
| ComplexHeatmap  | 2.22.0  | Visualisation of complex data matrices  |
| pheatmap        | 1.0.13  | Heatmap generation                      |

The primary differential gene expression analysis uses DESeq2.

## Python and Machine Learning

Python-based data analysis and machine learning are planned using Google Colab.

| Software     | Version        | Purpose                                |
| ------------ | -------------- | -------------------------------------- |
| Python       | To be recorded | Programming language for data analysis |
| pandas       | To be recorded | Data manipulation and analysis         |
| NumPy        | To be recorded | Numerical computing                    |
| scikit-learn | To be recorded | Machine learning models and evaluation |

Versions will be recorded from the Google Colab runtime used for the machine-learning analysis.

## Computational Environments

* **RNA-seq processing:** Ubuntu/WSL using the local Conda environment named `rnaseq`.
* **Statistical analysis:** R and DESeq2 within the local `rnaseq` environment.
* **Python and machine learning:** Google Colab, with runtime versions to be documented when the analysis begins.

## Environment Reproducibility

The exported Conda environment configuration is stored at:

`environment/rnaseq_environment.yml`

On a system with Conda installed, the environment can be recreated using:

```bash
conda env create -f environment/rnaseq_environment.yml
```

The environment name is specified in the YAML file. If an environment named `rnaseq` already exists, use a different environment name or remove the existing environment only if you are certain it is no longer needed.

The exported configuration records the packages in the current environment. Successful recreation may depend on package availability and compatibility with the operating system.