# ================================
# 04_gwas_qc_plots.R
# GWAS QC and Manhattan plots
# ================================

library(qqman)
library(ggplot2)
library(dplyr)

# Load cleaned data
load("data/clean/gwas_clean.RData")

# ----------------------------
# Function: Calculate Lambda GC
# ----------------------------
calc_lambda <- function(p){
  chisq <- qchisq(1 - p, 1)
  lambda <- median(chisq, na.rm = TRUE) / 0.456
  return(lambda)
}

# ----------------------------
# Function: Manhattan plot
# ----------------------------
plot_manhattan <- function(df, trait_name){

  df <- df %>%
    mutate(P = p)

  png(paste0("results/", trait_name, "_manhattan.png"), width=1200, height=600)

  manhattan(
    df,
    chr = "chr",
    bp = "pos",
    p = "P",
    snp = "rsid",
    main = trait_name,
    genomewideline = -log10(5e-8),
    suggestiveline = -log10(1e-5)
  )

  dev.off()
}

# ----------------------------
# Run QC for all traits
# ----------------------------
traits <- list(
  clean_ins_sec,
  clean_peak_ins,
  clean_proins,
  clean_isi,
  clean_homa_b,
  clean_homa_ir
)

trait_names <- c(
  "ins_sec",
  "peak_ins",
  "proins",
  "isi",
  "homa_b",
  "homa_ir"
)

# Create results folder
dir.create("results", showWarnings = FALSE)

# Loop
for (i in seq_along(traits)){

  df <- traits[[i]]
  name <- trait_names[i]

  # Lambda GC
  lambda <- calc_lambda(df$p)
  cat(name, "Lambda GC:", lambda, "\n")

  # Manhattan
  plot_manhattan(df, name)
}

cat("✔ QC and Manhattan plots completed\n")
