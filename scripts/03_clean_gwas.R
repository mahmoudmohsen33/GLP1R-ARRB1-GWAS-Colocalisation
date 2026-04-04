# ================================
# 03_clean_gwas.R
# Clean and harmonize GWAS datasets
# ================================

library(dplyr)
library(data.table)

# Load raw GWAS
load("data/raw/gwas_raw.RData")

# ----------------------------
# Sample sizes (from metadata)
# ----------------------------
n_ins_sec   <- 527
n_peak_ins  <- 2337
n_proins    <- 10701
n_isi       <- 16753
n_homa_b    <- 36466
n_homa_ir   <- 37037

# ----------------------------
# Cleaning function
# ----------------------------
clean_magic <- function(df, study_n, trait_label){

  df2 <- df

  # Keep autosomes only
  df2$chrom <- as.character(df2$chrom)
  df2 <- df2[df2$chrom %in% as.character(1:22), ]

  # Remove multi-allelic variants
  df2 <- df2[!grepl(",", df2$alt, fixed = TRUE), ]

  # Ensure numeric
  df2$ES <- as.numeric(df2$ES)
  df2$SE <- as.numeric(df2$SE)
  df2$pval <- as.numeric(df2$pval)
  df2$AF <- as.numeric(df2$AF)

  # Remove invalid rows
  df2 <- df2 %>%
    filter(!is.na(ES), !is.na(SE), SE > 0) %>%
    filter(!is.na(pval), pval > 0, pval <= 1)

  # Sample size
  df2$N <- ifelse("SS" %in% names(df2), df2$SS, study_n)

  # Standard format
  df2 <- df2 %>%
    select(rsid, chrom, pos, ref, alt, ES, SE, pval, AF, N) %>%
    rename(
      chr = chrom,
      other_allele = ref,
      effect_allele = alt,
      beta = ES,
      se = SE,
      p = pval,
      af = AF,
      n = N
    )

  # Add trait label
  df2$trait <- trait_label

  # Order
  df2$chr <- as.integer(df2$chr)
  df2 <- df2[order(df2$chr, df2$pos), ]

  return(df2)
}

# ----------------------------
# Apply cleaning
# ----------------------------
clean_ins_sec  <- clean_magic(ins_sec,   n_ins_sec,  "Insulin secretion")
clean_peak_ins <- clean_magic(peak_ins,  n_peak_ins, "Peak insulin")
clean_proins   <- clean_magic(proins,    n_proins,   "Proinsulin")
clean_isi      <- clean_magic(isi,       n_isi,      "ISI")
clean_homa_b   <- clean_magic(homa_b,    n_homa_b,   "HOMA-B")
clean_homa_ir  <- clean_magic(homa_ir,   n_homa_ir,  "HOMA-IR")

# ----------------------------
# Save cleaned data
# ----------------------------
dir.create("data/clean", showWarnings = FALSE)

save(
  clean_ins_sec,
  clean_peak_ins,
  clean_proins,
  clean_isi,
  clean_homa_b,
  clean_homa_ir,
  file = "data/clean/gwas_clean.RData"
)

cat("✔ Cleaning completed successfully\n")
