# TNBC_GABA Project

This repository contains the final project for BS831 titled **"TNBC_GABA"** by Omnia Abdelrahman.

## Files Included
- [Full Report (PDF)](Omnia%20Abdelrahman%20BS831%20Final%20Project.pdf)  
- [Poster Presentation (PDF)](TNBC_Poster.pdf)  
- [Code Samples (PDF)](TNBC_Full_Coding.PDF)

## Summary
Triple-negative breast cancer (TNBC) is a highly aggressive subtype with limited treatment options due to the absence of hormone receptors. This project investigated the potential role of **GABA receptor signaling**, particularly the **GABRA3 subunit**, in driving tumor proliferation and epithelial-to-mesenchymal transition (EMT).

**Data Source:**  
Public RNA-seq data from **213 TCGA-BRCA samples** (97 TNBC tumors, 116 normal breast tissues). TNBC cases were defined by negative expression of ESR1, PGR, and ERBB2.

**Analytic Workflow:**  
- Differential expression analysis of GABA subunits using `DESeq2`  
- Pathway enrichment with `fgsea` and MSigDB hallmark sets  
- Visualization with `ggplot2` and `pheatmap`  
- Subgroup stratification based on EMT markers to explore downstream gene expression programs  

## Tools & Packages
- **R 4.4**  
- `DESeq2`, `TCGAbiolinks`, `fgsea`, `ggplot2`, `pheatmap`  
- Data wrangling and QC with `dplyr` and `tidyr`  

## Author
**Omnia Abdelrahman, BDS, MPH**  
Graduate Student, Epidemiology & Biostatistics  
Boston University School of Public Health

## Contact
For questions or collaborations, feel free to open an issue or email me at [omnia@bu.edu](mailto:omnia@bu.edu).
