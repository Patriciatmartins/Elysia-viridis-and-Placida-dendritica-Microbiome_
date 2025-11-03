# Microbiome Characterization of *Elysia viridis* and *Placida dendritica*

This repository contains the R script and associated data for reproducing the analyses presented in the manuscript:

**"Microbiome characterization of the sea slugs *Elysia viridis* and *Placida dendritica*: insights into potential roles in kleptoplasty"**

---

## 📜 Description

This script performs microbiome analyses and generates the figures included in the study. Specifically, it covers:

- **Rarefaction curves** – Figure 1  
- **Venn diagram of shared zOTUs** – Figure 2  
- **Alpha diversity analysis** – Figure 3  
- **Beta diversity (PCoA & PERMANOVA)** – Figure 4  
- **Comparative taxonomic composition** – Figure 5  
- **Core microbiome analysis**  

All steps use standard R packages for microbiome data processing and visualization.

---

## 📁 Required Files

To run the script, ensure the following input files are placed in your working directory:

- `zotus.OTU.table_new.txt` – zOTU abundance and taxonomy table  
- `Labels.txt` – Sample metadata file  

---

## 🔧 How to Run

1. Open the script `R_script_updated.R` in RStudio.
2. Make sure the input files are in the same folder.
3. The script will automatically set the working directory to the script location using:

```r
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
