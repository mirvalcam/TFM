# 🧬 Microarray and Proteomic Analysis in Aging

## 📄 Description
This repository contains transcriptomic and proteomic data analysis from wild-type and galectin-3 knockout mouse microglia at different ages (6 and 24 months). The project investigates how aging alters microglial function—particularly inflammation, synaptic regulation, and lipid metabolism—and how deletion of galectin-3 modulates these processes. Results highlight galectin-3 as a key regulator of age-related neuroinflammation and a potential therapeutic target for neurodegenerative diseases.

## ⚙️ Requirements
- **R version:** ≥ 4.4.2 (tested on Windows 10)  
- **Required packages in microarrays analysis:**  
```r
install.packages(c(
  "oligo", 
  "clariomsmousetranscriptcluster.db", 
  "AnnotationDbi", 
  "arrayQualityMetrics", 
  "FactoMineR", 
  "factoextra", 
  "limma", 
  "tidyverse", 
  "clusterProfiler", 
  "org.Mm.eg.db", 
  "pathview", 
  "kableExtra"
))
```
- **Required packages in proteomic analysis:**  
```r
install.packages(c(
  "readxl",
  "dplyr",
  "tidyr",
  "ggplot2",
  "VennDiagram",
  "clusterProfiler",
  "org.Mm.eg.db",
  "kableExtra"
))
```
📥 Data access: Due to data policy, datasets are available upon request. Please contact miriamvaldayo99@gmail.com to obtain access.

## 📊 Results
Aging induces a pro-inflammatory microglial profile and downregulates genes involved in synaptic function and lipid metabolism.

Galectin-3 knockout partially restores synaptic remodeling and reduces expression of pro-inflammatory markers.

Proteomic analysis confirms the recovery of key molecular functions and reduced lipid droplet accumulation in the knockout condition.

Galectin-3 emerges as a potential therapeutic target for age-related neurodegenerative diseases.

## 📬 Contact
For questions or issues related to the code or data access, please contact:
📧 miriamvaldayo99@gmail.com
