# ================================
# 06_coloc_analysis.R
# Colocalization analysis (GWAS vs eQTL)
# ================================

library(coloc)
library(dplyr)
library(data.table)

# ----------------------------
# Function: Load region safely
# ----------------------------
load_region <- function(dir, trait){

  files <- list.files(dir, pattern = paste0("^", trait, "_.*\\.tsv$"), full.names = TRUE)

  if(length(files) == 0){
    stop(paste("No files found for", trait))
  }

  # choose wider region
  file <- files[which.max(nchar(files))]

  cat("Loading:", file, "\n")

  df <- fread(file) %>% as.data.frame()

  return(df)
}

# ----------------------------
# Function: Run coloc
# ----------------------------
run_coloc <- function(gwas, eqtl, trait_name, eqtl_name){

  # Ensure numeric
  gwas <- gwas %>%
    mutate(
      beta = as.numeric(beta),
      varbeta = as.numeric(se^2),
      pval = as.numeric(p),
      af = ifelse(is.na(af) | af <= 0 | af >= 1, 0.3, af)
    )

  eqtl <- eqtl %>%
    mutate(
      beta = as.numeric(beta),
      varbeta = as.numeric(varbeta),
      pval = as.numeric(pval),
      af = ifelse(is.na(af) | af <= 0 | af >= 1, 0.3, af)
    )

  # SNP overlap
  common_snps <- intersect(gwas$rsid, eqtl$rsid)

  cat("SNPs used:", length(common_snps), "\n")

  if(length(common_snps) < 20){
    warning("Too few SNPs for reliable coloc")
    return(NULL)
  }

  gwas_f <- gwas %>% filter(rsid %in% common_snps)
  eqtl_f <- eqtl %>% filter(rsid %in% common_snps)

  # Run coloc
  result <- coloc.abf(
    dataset1 = list(
      beta = gwas_f$beta,
      varbeta = gwas_f$varbeta,
      snp = gwas_f$rsid,
      MAF = gwas_f$af,
      N = gwas_f$n,
      type = "quant"
    ),
    dataset2 = list(
      beta = eqtl_f$beta,
      varbeta = eqtl_f$varbeta,
      snp = eqtl_f$rsid,
      MAF = eqtl_f$af,
      N = 500,
      type = "quant"
    )
  )

  return(result)
}

# ----------------------------
# Example usage
# ----------------------------

# Load GWAS region
gwas_dir <- "data/regions/GLP1R"
gwas <- load_region(gwas_dir, "peak_ins")

# Example eQTL (user should prepare this separately)
# eqtl <- read.csv("data/eqtl/glp1r_pancreas.csv")

# result <- run_coloc(gwas, eqtl, "peak_ins", "GLP1R_pancreas")

# print(result$summary)
