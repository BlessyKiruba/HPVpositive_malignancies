# HPV-Positive Malignancies: Integrative Multi-Omics Analysis

## Overview

Human papillomavirus (HPV) is a major global health concern responsible for nearly all cervical cancers and a growing proportion of head and neck cancers. Although these malignancies share HPV-driven molecular mechanisms, conserved oncogenic pathways and shared biomarkers across cancer types remain incompletely understood.

This repository contains the full computational workflow used to identify shared epithelial gene signatures across HPV-positive cervical cancer (CC) and head and neck cancer (HNC) using an integrative multi-omics framework combining bulk transcriptomics, network analysis, model-free gene separability, multivariate panel selection, and single-cell RNA sequencing.

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
| 3 | Model-free Per-gene Separability (CC/HNC consistency) | `gene_separability.py` |
| 4 | Multivariate Gene Panel Selection + Panel Size Sweep | `gene_panel_selection.py`, `panel_size_sweep.py`, `plot_panel_sweep.py` |
| 5 | External GEO Validation of Selected Panel | `geo_external_test.py` |
| 6 | Manuscript Figure Generation | `make_figures.py` |
| 7 | Tumor Microenvironment & Cell–Cell Communication | `cell chat.R` |

---

## Repository Structure

```
HPVpositive_malignancies/
│
├── DEG_Code.R               # Bulk RNA-seq differential gene expression (DESeq2)
├── WGCNA.R                  # Co-expression network construction and module–trait analysis
├── cell chat.R              # Cell–cell communication analysis (CellChat)
│
├── standardise.py           # Standardises labels/columns and writes CC_data.csv + HNC_data.csv
├── gene_separability.py     # Model-free per-gene separability, FDR, and cross-cohort consistency
├── gene_panel_selection.py  # Elastic-net stability + mRMR composite ranking for compact panel selection
├── panel_size_sweep.py      # Sweep panel size (k=1..all genes) and evaluate transfer performance
├── plot_panel_sweep.py      # Plots panel sweep metrics (AUC/BalAcc/CV)
├── geo_external_test.py     # External validation on independent GEO cohorts
└── make_figures.py          # Manuscript-ready ROC/PR, boxplots, stability, and confusion figures
```

---

## Prerequisites

### R (≥ 4.0)

```r
install.packages(c("ggplot2", "dplyr", "reshape2", "gridExtra", "writexl", "readxl", "cowplot"))
BiocManager::install(c("DESeq2", "WGCNA", "GEOquery", "Seurat", "rtracklayer"))
# CellChat — install from GitHub (see https://github.com/sqjin/CellChat for full instructions)
devtools::install_github("sqjin/CellChat")
```

### Python (≥ 3.8)

```bash
pip install pandas numpy scipy scikit-learn matplotlib openpyxl
```

---

## Data Requirements

The scripts expect the following input files (not included in the repository):

| File | Used by | Description |
|------|---------|-------------|
| `count_matrix.xlsx` | `DEG_Code.R` | Raw count matrix with a `SYMBOL` column |
| `Metadata.xlsx` | `DEG_Code.R` | Sample metadata with `Barcode` and `tissue_type` columns |
| `CC_data.csv` | Python panel scripts | Standardised cervical cohort matrix with `Sample_ID`, `Label`, and gene columns |
| `HNC_data.csv` | Python panel scripts | Standardised head and neck cohort matrix with `Sample_ID`, `Label`, and gene columns |
| `CC CESC_HPV_VST_52genes (1).csv` | `standardise.py` | Raw CC matrix used to generate `CC_data.csv` |
| `HNC VST_genes_of_interest_HPVpos_vs_HPVneg.xlsx` | `standardise.py` | Raw HNC matrix used to generate `HNC_data.csv` |
| `CC GSE151666_HPV_VST_52genes.csv` | `geo_external_test.py` | External GEO cervical validation cohort |
| `HNC VST_genes_of_interest_GSE74927_HPVpos_vs_HPVneg.xlsx` | `geo_external_test.py` | External GEO head and neck validation cohort |
| `gencode.v48.basic.annotation.gtf.gz` | `Infercnv.R` | GENCODE v48 gene annotation file — download from [GENCODE](https://www.gencodegenes.org/human/release_48.html) |

---

## Usage

### 1. Differential Gene Expression

Run `DEG_Code.R` to perform DESeq2-based differential expression between tumor and normal samples. Outputs a volcano plot and an Excel file of DEG results.

### 2. Co-expression Network Analysis

Run `WGCNA.R` to construct a signed co-expression network, identify modules correlated with tumor status, and export gene lists for downstream analysis.

### 3. Prepare and Standardise Input Matrices

Use `standardise.py` to harmonise sample ID and label columns and generate `CC_data.csv` and `HNC_data.csv`:

```bash
python standardise.py
```

### 4. Model-free Separability Ranking

Run per-gene separability and consistency analysis across CC and HNC:

```bash
python gene_separability.py
```

This writes `Separability_Results/separability_{CC,HNC}.csv` and `Separability_Results/separability_combined.csv`.

### 5. Multivariate Panel Selection and Panel-size Sweep

Build the compact panel using stability selection + mRMR, then evaluate performance across panel sizes:

```bash
python gene_panel_selection.py
python panel_size_sweep.py
python plot_panel_sweep.py
```

Outputs are written to `Panel_Selection_Results/`.

### 6. External GEO Validation and Figure Generation

Validate the selected panel on external GEO cohorts and generate manuscript figures:

```bash
python geo_external_test.py
python make_figures.py
```

Outputs are written to `GEO_test/` and `Manuscript_Results/figures/`.

### 7. Single-cell Validation

 `cell chat.R` for cell–cell communication profiling within the tumor microenvironment.

---

## License

This project is licensed under the **GNU General Public License v3.0 (GPL-3.0)**. See [LICENSE](https://www.gnu.org/licenses/gpl-3.0.html) for details.

---
