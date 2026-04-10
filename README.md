# Genetic dissection of GLP-1 signalling reveals distinct regulatory mechanisms at GLP1R and ARRB1 influencing insulin secretion in type 2 diabetes

This repository contains the code used for integrative genetic analysis of GLP1R and ARRB1 loci in insulin-related traits.

## Overview
We analysed GWAS summary statistics for nine insulin-related traits to:
- Estimate SNP-based heritability and genetic correlations (LDSC)
- Perform locus-specific association analysis
- Conduct Bayesian colocalisation (coloc)
- Perform fine-mapping (SuSiE)
- Integrate functional annotation (VEP, HaploReg, RegulomeDB)

## Key Findings
- GLP1R shows a diffuse association pattern
- ARRB1 shows structured genetic architecture
- Fine-mapping identified two candidate causal variants:
  - rs899115
  - rs1789685

## Data Sources
- OpenGWAS: https://gwas.mrcieu.ac.uk/
- MAGIC consortium (Chen et al., 2021)
- GTEx v8

## Pipeline

1. Data harmonisation
2. LDSC heritability & correlation
3. Regional association analysis
4. Colocalisation (coloc)
5. Fine-mapping (susieR)
6. Functional annotation

## Reproducibility

To reproduce results:

```bash
Rscript scripts/run_pipeline.R
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


---
