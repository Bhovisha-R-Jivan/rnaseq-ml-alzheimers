# RNA-seq and Machine Learning for Alzheimer's Disease

## Project Overview

This project investigates gene expression patterns associated with Alzheimer's disease (AD) using publicly available RNA-sequencing (RNA-seq) datasets.

The project combines a reproducible RNA-seq processing workflow with downstream differential gene expression analysis. Machine learning (ML) will be explored as a subsequent stage to investigate whether gene expression profiles can distinguish AD samples from controls.

## Objectives

1. Develop and document a reproducible RNA-seq processing workflow.
2. Perform quality control, read alignment, and gene-level quantification.
3. Investigate publicly available gene expression datasets for Alzheimer's disease.
4. Perform differential gene expression analysis using DESeq2.
5. Explore machine learning approaches for distinguishing AD samples from controls.
6. Document the analysis, code, and methodological decisions using Git and GitHub.

## Dataset Strategy

Due to the computational demands of processing large-scale RNA-seq datasets, a smaller dataset is used to demonstrate the RNA-seq processing workflow, including read quality assessment, adapter trimming, genome alignment, alignment quality control, and gene-level quantification.

The primary disease-specific analysis uses publicly available Alzheimer's disease gene-count data. This separates the demonstration of the RNA-seq processing workflow from the downstream statistical analysis of the AD dataset.

## Datasets

### RNA-seq Workflow Practice

A smaller RNA-seq dataset is used to demonstrate the processing workflow, from raw-read quality assessment to gene-level quantification.

The workflow includes FastQC, MultiQC, Trim Galore, STAR, SAMtools, and featureCounts.

### Alzheimer's Disease Analysis

**Primary dataset: GSE125583**

GSE125583 provides gene-level count data and sample metadata for the Alzheimer's disease differential gene expression analysis.

The analysis includes:

* Count-data and metadata validation.
* Assessment of the experimental design and sample composition.
* Differential gene expression analysis using DESeq2.
* Gene annotation.
* Sensitivity analysis of independent filtering.
* Exploratory principal component analysis (PCA) and sample-correlation analysis.
* Investigation of missing p-values and adjusted p-values.

A second dataset, **GSE261050**, was also identified for investigation. Its suitability for any additional analysis should be established before incorporating it into the workflow.

## Differential Gene Expression Analysis

The primary differential expression analysis uses DESeq2 in R.

The experimental design is:

`~ AGE + sex + diagnosis`

The model accounts for age and sex while estimating the association between diagnosis and gene expression.

Samples with missing values in the required model variables are excluded before fitting the model. Differential expression results are generated for the Alzheimer's disease group relative to the control group.

The analysis uses DESeq2's default independent filtering and Cook's-distance handling. Sensitivity analysis and diagnostic investigations are used to examine the effect of these procedures on the results.

## Tools

* Linux / Bash
* SRA Toolkit
* FastQC and MultiQC
* Trim Galore
* STAR
* SAMtools
* Subread / featureCounts
* R
* DESeq2
* Git and GitHub
* Python / scikit-learn (planned for the machine learning stage)

## Project Status

### Completed

* [x] Initial project documentation and repository setup
* [x] RNA-seq processing scripts developed
* [x] Study 1 count-data and metadata validation
* [x] Experimental design assessment
* [x] DESeq2 differential gene expression analysis
* [x] Gene annotation
* [x] Sensitivity analysis
* [x] PCA investigation
* [x] Investigation of missing p-values and adjusted p-values

### Remaining

* [ ] Final review of RNA-seq outputs and interpretation
* [ ] Finalise RNA-seq documentation and reproducibility notes
* [ ] Review Git tracking and repository organisation
* [ ] Investigate the suitability of GSE261050 for any additional analysis
* [ ] Machine learning analysis

## Reproducibility

The project uses numbered scripts in the `scripts/` directory to organise the workflow. Analysis outputs are stored under `analysis/`, while raw input data are stored under `raw_data/`.

Large raw datasets and generated analysis outputs should not be committed to GitHub unless there is a specific reason to include them. The repository should retain the scripts and documentation required to understand and reproduce the analysis, subject to data availability and computational requirements.