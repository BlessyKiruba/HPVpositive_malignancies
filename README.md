# HPV-Positive Malignancies: Integrative Multi-Omics Analysis

## Overview

Human papillomavirus (HPV) is a major global health concern responsible for nearly all cervical cancers and a growing proportion of head and neck cancers. Although these malignancies share HPV-driven molecular mechanisms, conserved oncogenic pathways and shared biomarkers across cancer types remain incompletely understood.

This repository contains the full computational workflow used to identify shared epithelial gene signatures across HPV-positive cervical cancer (CC) and head and neck cancer (HNC) using an integrative multi-omics framework combining bulk transcriptomics, network analysis, machine learning, and single-cell RNA sequencing.

> **Key finding:** A conserved four-gene signature — **ASF1B, DTL, CDKN2A, CLSPN** — is consistently upregulated in HPV-associated tumors and enriched within epithelial cells of HPV-positive malignancies.

---

## Paper

**AI-Driven Discovery and Single-Cell Validation of a Conserved Epithelial Signature in HPV-Positive Malignancies**

> *Manuscript in preparation / under review.*

---

## Study Workflow

The analysis pipeline integrates the following computational approaches in sequence:

| Step | Method | Script(s) |
|------|--------|-----------|
| 1 | Differential Gene Expression (DEG) Analysis | `DEG_Code.R` |
| 2 | Weighted Gene Co-expression Network Analysis (WGCNA) | `WGCNA.R` |
| 3 | Machine Learning-based Feature Selection | `RF_single_feature.py`, `SVM_single_feature.py`, `LR_single_feature.py`, `KNN_single_feature.py`, `LGBM_single_feature.py`, `MLP_single_feature.py` |
| 4 | Single-cell RNA-seq Validation (Pseudobulk DEG) | `pseudo-bulk degs.R` |
| 5 | Copy Number Variation Inference (inferCNV) | `Infercnv.R` |
| 6 | Tumor Microenvironment & Cell–Cell Communication | `cell chat.R` |

---

## Repository Structure

```
HPVpositive_malignancies/
│
├── DEG_Code.R               # Bulk RNA-seq differential gene expression (DESeq2)
├── WGCNA.R                  # Co-expression network construction and module–trait analysis
├── pseudo-bulk degs.R       # Pseudobulk DEG analysis from scRNA-seq data
├── Infercnv.R               # Copy number variation inference from scRNA-seq (inferCNV)
├── cell chat.R              # Cell–cell communication analysis (CellChat)
│
├── RF_single_feature.py     # Random Forest classifier (single-feature evaluation)
├── SVM_single_feature.py    # Support Vector Machine classifier
├── LR_single_feature.py     # Logistic Regression classifier
├── KNN_single_feature.py    # K-Nearest Neighbor classifier
├── LGBM_single_feature.py   # LightGBM classifier
└── MLP_single_feature.py    # Multi-Layer Perceptron (neural network) classifier
```

---

## Prerequisites

### R (≥ 4.0)

```r
install.packages(c("ggplot2", "dplyr", "reshape2", "gridExtra", "writexl", "readxl", "cowplot"))
BiocManager::install(c("DESeq2", "WGCNA", "GEOquery", "Seurat", "infercnv", "rtracklayer"))
# CellChat — install from GitHub (see https://github.com/sqjin/CellChat for full instructions)
devtools::install_github("sqjin/CellChat")
```

### Python (≥ 3.8)

```bash
pip install pandas numpy scikit-learn lightgbm tqdm
```

---

## Data Requirements

The scripts expect the following input files (not included in the repository):

| File | Used by | Description |
|------|---------|-------------|
| `count_matrix.xlsx` | `DEG_Code.R` | Raw count matrix with a `SYMBOL` column |
| `Metadata.xlsx` | `DEG_Code.R` | Sample metadata with `Barcode` and `tissue_type` columns |
| `CC.csv` | ML scripts | Feature matrix for cervical cancer samples |
| `HNC.csv` | ML scripts | Feature matrix for head and neck cancer samples |
| `gencode.v48.basic.annotation.gtf.gz` | `Infercnv.R` | GENCODE v48 gene annotation file — download from [GENCODE](https://www.gencodegenes.org/human/release_48.html) |

---

## Usage

### 1. Differential Gene Expression

Run `DEG_Code.R` to perform DESeq2-based differential expression between tumor and normal samples. Outputs a volcano plot and an Excel file of DEG results.

### 2. Co-expression Network Analysis

Run `WGCNA.R` to construct a signed co-expression network, identify modules correlated with tumor status, and export gene lists for downstream analysis.

### 3. Machine Learning Feature Selection

Each Python script evaluates every candidate gene as a single-feature classifier across CC and HNC datasets. Results are ranked by average accuracy across both datasets:

```bash
python RF_single_feature.py
python SVM_single_feature.py
# ... repeat for other classifiers
```

Output CSVs are saved to a dedicated results folder (e.g., `RF_Results/` for Random Forest, `SVM_Results/` for SVM, and so on for each classifier).

### 4. Single-cell Validation

Run `pseudo-bulk degs.R` for pseudobulk DEG analysis, `Infercnv.R` for copy number inference, and `cell chat.R` for cell–cell communication profiling within the tumor microenvironment.

---

## License

This project is licensed under the **GNU General Public License v3.0 (GPL-3.0)**. See [LICENSE](https://www.gnu.org/licenses/gpl-3.0.html) for details.

---

## Contact

For questions or correspondence regarding the code, please contact:
**Naisarg Patel** — naisargbpatel14 [at] gmail [dot] com
