**Overview**

Human papillomavirus (HPV) is a major global health concern and is responsible for nearly all cervical cancers and a growing proportion of head and neck cancers. Although these malignancies share HPV-driven molecular mechanisms, conserved oncogenic pathways and shared biomarkers across these cancer types remain incompletely understood.
This repository contains the computational workflow used to identify shared epithelial gene signatures across HPV-positive cervical cancer and head and neck cancer using an integrative multi-omics framework combining bulk transcriptomics, network analysis, machine learning, and single-cell RNA sequencing.
The study identifies a conserved four-gene signature (ASF1B, DTL, CDKN2A, CLSPN) that is consistently upregulated in HPV-associated tumors and enriched within epithelial cells of HPV-positive malignancies.

**Study Workflow**

The analysis pipeline integrates multiple computational approaches:
  - Differential Gene Expression (DEG) Analysis
  - Weighted Gene Co-expression Network Analysis (WGCNA)
  - Machine Learning-based Feature Selection
  - Single-cell RNA-seq Validation
  - Tumor Microenvironment and Cellular Interaction Analysis

These analyses collectively identify conserved oncogenic mechanisms across HPV-driven cancers.

**Repository Structure**

HPVpositive_malignancies/
│
├── DEG_Code.R # Differential gene expression analysis
├── WGCNA.R # Co-expression network analysis
├── pseudobulk_degs.R # Pseudobulk DEG analysis from scRNA-seq
├── Infercnv.R # Copy number variation inference from scRNA-seq
├── cell_chat.R # Cell-cell communication analysis
│
├── RF_single_feature.py # Random Forest classifier
├── SVM_single_feature.py # Support Vector Machine classifier
├── LR_single_feature.py # Logistic Regression classifier
├── KNN_single_feature.py # K-Nearest Neighbor classifier
├── LGBM_single_feature.py # LightGBM classifier
├── MLP_single_feature.py # Neural network classifier
├── MLP_single_feature.py     # Neural network classifier


