# Dataset Documentation

## Study 1: Alzheimer's Disease RNA-seq

### Study identification

- GEO accession: GSE125583
- BioProject accession: PRJNA516886
- Tissue: Fusiform gyrus
- Organism: Homo sapiens
- Intended analysis: Differential gene expression, followed by downstream
  analyses where appropriate.

### Input files

| File | Description |
|---|---|
| `Study1_gene_counts.gz` | Combined gene-level count matrix |
| `Study1_MetaData.csv` | Sample metadata; semicolon-delimited despite the `.csv` extension |

The files are currently stored in the user's local Downloads directory.
They should be placed in the project's local `data/raw/` directory when
the data-access and storage arrangements have been established.

### Count matrix

- Gene identifiers are stored in the first column, `GeneID`.
- Sample identifiers are stored in the remaining columns.
- The matrix contains 289 sample columns.
- The matrix contains 39,376 gene rows.
- All inspected count values were non-negative integers.
- The working matrix was compared with
  `GSE125583_raw_counts_GRCh38.p13_NCBI.tsv.gz`.
- The two decompressed files were confirmed to be identical.

The source dataset also provides individual sample files with columns
for gene identifier, normalised expression (`VALUE`), and uniquely
aligning read counts (`count`). The normalised expression column must
not be substituted for raw counts in DESeq2.

### Sample metadata

- Total samples: 289
- Alzheimer's disease: 219
- Controls: 70
- Sample identifiers match between the count matrix and metadata.
- One sample has missing sex information.
- The same sample has age recorded as `not applicable`.
- The treatment of this sample in adjusted analyses remains undecided.

### Validation status

Completed:
- Count matrix dimensions checked.
- Sample identifiers checked against metadata.
- Count values checked for integer formatting and non-negativity.
- Working count matrix compared with the NCBI reference file.

Pending:
- Automated validation script.
- Investigation of missing demographic information.
- Assessment of group balance and covariate suitability.
- Finalisation of the statistical design.
- Confirmation of all relevant source-processing details.

### Reproducibility notes

Validation steps should be implemented in version-controlled scripts.
Input data and generated results are excluded from Git by the repository's
current ignore rules. Data sources, software versions, analysis decisions,
and execution instructions should be documented separately.

### Script Paths
scripts/8. validate_study1_data.R

analysis/study1/data_validation_report.txt

- The validation script imports the combined count matrix and semicolon-delimited metadata.
- The matrix contains 39,376 gene rows and 289 samples.
- All gene identifiers and sample identifiers are unique.
- All count values passed the checks for non-negative integer formatting.
- All 289 sample identifiers match between the count matrix and metadata.
- The dataset contains 219 Alzheimer's disease samples and 70 controls.
- One Alzheimer's disease sample, GSM3577605, has missing age and sex information.
- The mean age is 83.38 years for the Alzheimer's disease group and 87.00 years for controls, calculated from samples with recorded ages.
- The sex distribution is 97 female, 121 male and one missing for Alzheimer's disease, compared with 33 female and 37 male controls.
- These checks establish structural consistency, but further assessment of the experimental design is required before differential expression analysis.