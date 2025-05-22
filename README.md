# 💡 The gene expression signature of electrical stimulation in the human brain – Multi-Omics Analysis

This repository contains analysis scripts and relevant code used for the integrative analysis of **single-nucleus multi-omics** (snRNA-seq + snATAC-seq) and **bulk RNA-seq** data to study the molecular effects of **direct electrical stimulation (DES)** in the **human brain**. The project focuses on identifying the **molecular signature** to DES and characterizing transcriptional and epigenomic changes at both bulk and single-cell resolution.

## 🧠 Project Overview

**Direct electrical stimulation (DES)** has long been used as a clinical tool to map brain function in neurosurgical settings. However, the molecular consequences of DES on the human brain remained largely unexplored. In this study, we use **state-of-the-art transcriptomic and epigenomic sequencing techniques** to investigate the **cell type-specific** and **bulk tissue** responses to DES.

Key findings include:
- Robust and reproducible gene expression changes in **microglia**, particularly in genes involved in **cytokine activity**.
- Cross-validation in a mouse model.
- Application of a custom **deep learning model** to identify cell-type-specific molecular responses.


## 🧬 Data Types

- **Bulk RNA-Seq** from two different experimental paradigms (stimulated (human and mouse) & unstimulated (human))
- **Single-nucleus multi-omics** of snRNA and snATAC data (human)

## 📁 Repository Structure

.
├── README.md                                # This file
└── src/                                     # Analysis scripts
    ├── bulk_rnaseq/
    │   ├── 01_differential_expression.R     # DESeq2 analysis for bulk RNA-seq data
    │   └── 02_deconvolution.R               # Cell-type deconvolution using reference profiles
    └── sn_multiomics/
        ├── 01_fastq-processing.sh           # Preprocessing pipeline (Cell Ranger, etc.)
        ├── 02_samples-seurat-processing.R   # Quality control and clustering in Seurat
        ├── 03_prepare-for-DCA.R             # Prepares Seurat object for deep count autoencoder (DCA)
        ├── 04_DCA-imputation.sh             # Runs DCA (deep learning-based imputation)
        ├── 05_01_RNA-pseudobulk-calculations.R   # Pseudobulk expression from snRNA
        ├── 05_02_ATAC-pseudobulk-calculations.R  # Pseudobulk accessibility from snATAC
        ├── 06_GEX-differential-expression.R # DE analysis on imputed gene expression data
        ├── 07_enrichment-analysis.R         # Gene set enrichment (GO, ActivePathways)
        ├── 09_comparison-w-lit.R            # Compares DEGs with previously published datasets
        ├── 10_ATAC-motif-analysis.R         # Transcription factor motif enrichment
        ├── 11_neuroestimator.R              # Deep learning-based tool to estimate neuroactivity
        ├── 12_final-manuscript-plots.R      # Generates plots used in the manuscript
        └── 98_manuscript-color-palettes.R   # Centralized color theme used across plots


## ⚙️ Setup

### Data Access
The sequencing datasets used in this study are available via the NCBI Gene Expression Omnibus (GEO) under accession number:
- 🔗 [GSE224952](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224952)
This repository uses placeholder paths and dummy data where necessary.

## 📊 Main Tools & Libraries

- **R / Bioconductor**: DESeq2, Seurat, Signac, ActivePathways, goSeq
- **Deep learning models**: NEUROeSTIMator and DCA
- **Cell Ranger (10x Genomics)**: Preprocessing of snRNA/snATAC data

## 📈 Reproducibility

Final manuscript plots and color schemes:

- `12_final-manuscript-plots.R`: Generates all figures
- `98_manuscript-color-palettes.R`: Unified color settings

Each script is standalone and modular. Input/output files are annotated inside the scripts. 

## 🧾 Citation

If you use or adapt this code, please cite the corresponding manuscript:

> **[Title Placeholder]**  
> *placeholder*, et al.  
> *Journal Placeholder*, 2025.  

## 📬 Contact

For questions, collaboration, or issues with the code, please contact:

**Muhammad Elsadany**  
📧 muhammad-elsadany@uiowa.edu  
🧪 Lab: [Michaelson lab]  
🏛 Institution: [University of Iowa]
