# Genetic dissection of GLP-1 signalling reveals distinct regulatory mechanisms at GLP1R and ARRB1 influencing insulin secretion in type 2 diabetes

## Overview

This repository contains a complete bioinformatics pipeline for:

* Multi-trait GWAS analysis (MAGIC consortium)
* Genetic correlation using LDSC
* Regional association analysis (GLP1R & ARRB1 loci)
* Colocalization analysis with GTEx eQTL data
* Functional annotation (VEP, CADD, SIFT, PolyPhen)
* Candidate variant prioritization

---

## Biological Motivation

GLP1R and ARRB1 are key regulators of insulin secretion and GLP-1 signaling pathways.
This project aims to:

* Identify shared genetic signals between GWAS traits and gene expression
* Distinguish regulatory mechanisms at each locus
* Prioritize functional variants influencing type 2 diabetes

---

## ⚙️ Pipeline Overview

1. **GWAS Data Processing**

   * Import VCF (MAGIC)
   * Clean and harmonize summary statistics

2. **Genetic Correlation (LDSC)**

   * Estimate trait–trait correlations
   * Visualize correlation matrix

3. **Regional Analysis**

   * Extract loci:

     * GLP1R (chr6)
     * ARRB1 (chr11)

4. **Colocalization (COLOC)**

   * Integrate GWAS with GTEx eQTL
   * Compute posterior probabilities (PP.H4)

5. **Functional Annotation**

   * Variant Effect Predictor (VEP)
   * CADD, SIFT, PolyPhen scores

6. **Visualization**

   * Regional association plots
   * Heatmaps and candidate tables

---

## Project Structure

```
scripts/        # Analysis scripts
data/           # Input data (excluded from repo)
results/        # Output results
docs/           # Documentation
```

---

## 📦 Requirements

```r
install.packages(c(
  "data.table", "dplyr", "ggplot2", "coloc",
  "qqman", "pheatmap", "corrplot", "LDlinkR", "ieugwasr"
))
```

---

##  Environment Setup

Set your OpenGWAS token:

```r
Sys.setenv(OPENGWAS_JWT = Sys.getenv("OPENGWAS_JWT"))
```

---

## Output

* Genetic correlation matrix
* Colocalization results (PP.H4)
* Candidate variant table
* Publication-ready figures

---

## Author

Mahmoud M. Omran
Nile University – Bioinformatics

---

##  Notes

* GTEx Portal filtered data is **not suitable for colocalization**
* Use **GTEx allpairs files** for full analysis

---
