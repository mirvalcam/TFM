# Transcriptomic Ageing — Microarray Analysis



## Table of Contents
- [Description](#description)
- Requirements
  R version: “4.4.2.”
  install.packages(c(
  "oligo", "clariomsmousetranscriptcluster.db", "AnnotationDbi",
  "arrayQualityMetrics", "FactoMineR", "factoextra", "limma",
  "tidyverse", "clusterProfiler", "org.Mm.eg.db", "pathview",
  "kableExtra"))
- [Installation](#installation)
  git clone <tu-repo-url>
  cd <tu-repo>
  Rscript -e "install.packages(c(...))"
  
  After installing all the packages, create a folder (e.g., data/) and deposit the microarray files provided by contacting Rocío Martínez de Pablos (depablos@us.es).

- [Repository Structure](#repository-structure)
  data/         # Archivos CEL con datos de microarrays
  scripts/      # Código de análisis (preprocesamiento, DE, visualización)
  results/      # Figuras y tablas generadas

- [Quick Start](#quick-start)
- [Results](#results)
- [License](#license)
- [Citation](#citation)
- [Contributing](#contributing)
- [Contact](#contact)

## Description


## Requirements
* R ≥ 4.4.2 (tested on Windows 10)

