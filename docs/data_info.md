# Dataset Documentation

## Study 1: Alzheimer's Disease RNA-seq

### Study Identification

* GEO accession: GSE125583
* BioProject accession: PRJNA516886
* Organism: *Homo sapiens*
* Tissue: Fusiform gyrus
* Intended analysis: Differential gene expression, followed by downstream analyses where appropriate.

### Input Files

| File                                     | Description                                                       |
| ---------------------------------------- | ----------------------------------------------------------------- |
| `raw_data/Study1_gene_counts.gz`         | Combined gene-level count matrix                                  |
| `raw_data/Study1_MetaData.csv`           | Sample metadata; semicolon-delimited despite the `.csv` extension |
| `raw_data/Human.GRCh38.p13.annot.tsv.gz` | Gene annotation reference                                         |

Raw input data are stored locally under `raw_data/` and excluded from Git tracking.

### Count Matrix

* Gene identifiers are stored in the first column, `GeneID`.
* Sample identifiers are stored in the remaining columns.
* Original matrix dimensions: 39,376 genes and 289 samples.
* All inspected count values were non-negative integers.
* Sample identifiers match between the count matrix and metadata.
* The working count matrix was compared with `GSE125583_raw_counts_GRCh38.p13_NCBI.tsv.gz`; the decompressed files were confirmed to be identical.

The source dataset also provides individual sample files containing gene identifiers, normalised expression values (`VALUE`), and uniquely aligning read counts (`count`). Normalised expression values must not be substituted for raw counts in DESeq2.

### Sample Metadata

The original dataset contains 289 samples:

| Diagnosis           | Original samples |
| ------------------- | ---------------: |
| Alzheimer's disease |              219 |
| Controls            |               70 |
| Total               |              289 |

One Alzheimer's disease sample, `GSM3577605`, has missing age and sex information.

For the differential expression analysis, samples missing values in the required model variables were excluded. Consequently, the final analysis included 288 samples:

| Diagnosis           | Samples analysed |
| ------------------- | ---------------: |
| Alzheimer's disease |              218 |
| Controls            |               70 |
| Total               |              288 |

The mean age was approximately 83.38 years for the Alzheimer's disease group and 87.00 years for controls, calculated using samples with recorded ages.

The recorded sex distribution was:

* Alzheimer's disease: 97 female, 121 male, and one sample with missing sex.
* Controls: 33 female and 37 male.

These demographic summaries describe the available metadata and do not, by themselves, establish that all potential confounding factors have been addressed.

### Data Validation

The following checks were completed using `scripts/8. validate_study1_data.R`:

* Count matrix dimensions checked.
* Gene and sample identifiers checked for uniqueness.
* Count values checked for non-negative integer formatting.
* Sample identifiers checked against the metadata.
* Diagnosis group counts examined.
* Missing age and sex information identified.
* Working count matrix compared with the NCBI reference file.

The validation report is available at:

`analysis/study1/data_validation_report.txt`

### Experimental Design Assessment

The experimental design was assessed using `scripts/9. assess_study1_design.R`.

The differential expression model uses:

`~ AGE + sex + diagnosis`

The design assessment confirmed that the complete-case dataset contained 288 samples and that the model matrix was full rank.

The design assessment report is available at:

`analysis/study1/design_assessment_report.txt`

### Differential Gene Expression Analysis

Differential gene expression analysis was performed using DESeq2 in `scripts/10. deseq2_study1.R`.

The analysis:

* Excluded the sample with missing age and sex information.
* Removed genes with zero counts across all analysed samples.
* Retained 38,192 genes for the DESeq2 analysis.
* Modelled age and sex alongside diagnosis.
* Compared Alzheimer's disease samples with controls.
* Used DESeq2's default independent filtering and Cook's-distance handling.

The primary results are available at:

* `analysis/study1/study1_deseq2_results.csv`
* `analysis/study1/study1_deseq2_results_annotated.csv`
* `analysis/study1/study1_normalised_counts.csv`
* `analysis/study1/study1_deseq2_report.txt`

The analysis produced 14,084 genes with adjusted p-values below 0.05 under the default DESeq2 results procedure. This is a statistical result, not proof that every identified gene has a biologically meaningful or independently validated association with Alzheimer's disease.

### Additional Analyses

The following supplementary analyses were completed:

**Gene annotation**

* Script: `scripts/11. annotate_study1_results.R`
* Output: `analysis/study1/study1_deseq2_results_annotated.csv`

**Sensitivity analysis**

* Script: `scripts/12. sensitivity_analysis_study1.R`
* Outputs are stored under `analysis/study1/sensitivity/`.
* The analysis compared the default DESeq2 results procedure with results calculated without independent filtering.

**PCA investigation**

* Script: `scripts/13. investigate_study1_pca.R`
* Outputs are stored under `analysis/study1/pca_investigation/`.
* Associations between metadata variables and principal components were explored. These are exploratory associations and do not establish causation.

**Missing p-value investigation**

* Script: `scripts/14. investigate_independent_filtering.R`
* Outputs are stored under `analysis/study1/filtering_investigation/`.
* The investigation summarised expression levels for genes with missing adjusted p-values and documented the limitations of interpreting missing p-values from the exported results alone.

### Reproducibility Notes

The workflow is organised into numbered scripts under `scripts/`.

* Raw input data are stored under `raw_data/`.
* Study 1 analysis outputs are stored under `analysis/study1/`.
* Generated analysis outputs and raw data are excluded from Git tracking by the repository's ignore rules.
* Scripts and documentation are version-controlled.
* Software versions are documented in `environment/software_versions.md`.

The scripts document the analytical steps, while the reports and result files record the outputs produced during execution. Reproducing the analysis may require access to the original input data and a compatible software environment.