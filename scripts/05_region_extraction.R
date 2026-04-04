# ================================
# 05_region_extraction.R
# Extract GWAS regions for GLP1R and ARRB1
# ================================

library(dplyr)
library(data.table)

# Load cleaned GWAS
load("data/clean/gwas_clean.RData")

# Create output folders
dir.create("data/regions", showWarnings = FALSE)
dir.create("data/regions/GLP1R", recursive = TRUE, showWarnings = FALSE)
dir.create("data/regions/ARRB1", recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Define gene regions (±250kb)
# ----------------------------
regions <- list(
  GLP1R = list(chr = 6, start = 39048781 - 250000, end = 39091303 + 250000),
  ARRB1 = list(chr = 11, start = 68296181 - 250000, end = 68496181 + 250000)
)

# ----------------------------
# Function: extract region
# ----------------------------
extract_region <- function(df, region, trait_name, gene_name){

  subset <- df %>%
    filter(
      chr == region$chr,
      pos >= region$start,
      pos <= region$end
    )

  file_name <- paste0(
    "data/regions/", gene_name, "/",
    trait_name, "_chr", region$chr, "_",
    region$start, "_", region$end, ".tsv"
  )

  fwrite(subset, file_name, sep = "\t")

  cat("Saved:", file_name, " | SNPs:", nrow(subset), "\n")
}

# ----------------------------
# Traits list
# ----------------------------
traits <- list(
  ins_sec = clean_ins_sec,
  peak_ins = clean_peak_ins,
  proins = clean_proins,
  isi = clean_isi,
  homa_b = clean_homa_b,
  homa_ir = clean_homa_ir
)

# ----------------------------
# Run extraction
# ----------------------------
for (trait_name in names(traits)){

  df <- traits[[trait_name]]

  # GLP1R
  extract_region(df, regions$GLP1R, trait_name, "GLP1R")

  # ARRB1
  extract_region(df, regions$ARRB1, trait_name, "ARRB1")
}

cat("✔ Region extraction completed\n")
