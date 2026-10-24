# RNA-seq and Machine Learning for Alzheimer's disease

## Project Overview

This project investigates gene expression patterns associated with Alzheimer's Disease (AD) using publicly available RNA-sequencing (RNA-seq) datasets

The project will involve RNA-seq data processing, differential gene expression analysis and machine learning (ML) to explore whether gene expression profiles can distinguish AD samples from controls

## Objectives
1. Process RNA-seq data using a reproducible bioinformatics workflow.
2. Perform quality control, read alignment, and gene-level quantification.
3. Investigate publicly available gene expression data for AD
4. Conduct differential gene expression analysis.
5. Explore ML approaches for distinguishing AD samples from controls.
6. Document the analysis, code, and decisions using Git and GitHub.

## Dataset strategy

Due to the computational demands of processing large-scale RNA-seq datasets, a smaller, dataset will be used to establish and demonstrate the RNA-seq processing workflow. This will include quality control, read alignment and gene-level quantification.

The disease-specific analysis will focus the identified AD datasets. Gene-level count data will be obtained or generated, subject to availability and suitability, for differential gene expression analysis and subsequent ML.

This approach separates the development and testing of the RNA-seq processing workflow from the primary research objective, while ensuring that the downstream statistical and machine learning analyses remain focused on AD datasets.

## Datasets

### RNA-seq workflow practice
A smaller RNA-seq dataset will be used to develop and demonstrate the processing workflow, from raw-read quality control to gene-level quantification.

### Alzheimer's Disease Analysis

Two publicly available datasets have been selected for investigation:

- **GSE125583**
- **GSE261050**

The availability, format, and suitability of gene expression count data will be assessed for each dataset before downstream analysis. Sample metadata will also be examined to establish the sample composition, relevant biological variables, and data quality.

## Tools

- Linux / Bash
- STAR
- FastQC and MultiQC
- Trim Galore
- SAMtools
- featureCounts
- R / DESeq2
- Python / scikit-learn
- Git and GitHub

## Project Status

**Current stage:** Repository setup and dataset assessment.

- [x] Initial project documentation
- [ ] RNA-seq processing workflow
- [ ] Alzheimer's dataset count-data assessment
- [ ] Metadata quality assessment
- [ ] Differential gene expression analysis
- [ ] Machine learning analysis