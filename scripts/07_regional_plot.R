# ================================
# 07_regional_plot.R
# Regional association plot (GWAS + eQTL)
# ================================

library(ggplot2)
library(dplyr)
library(data.table)

# ----------------------------
# Prepare GWAS
# ----------------------------
prepare_gwas <- function(df){

  df <- as.data.frame(df)

  df <- df %>%
    mutate(
      pval = as.numeric(p),
      logp = -log10(pval),
      source = "GWAS"
    ) %>%
    select(rsid, pos, logp, source)

  return(df)
}

# ----------------------------
# Prepare eQTL
# ----------------------------
prepare_eqtl <- function(df){

  df <- as.data.frame(df)

  df <- df %>%
    mutate(
      logp = -log10(pval),
      source = "eQTL"
    ) %>%
    select(rsid, pos, logp, source)

  return(df)
}

# ----------------------------
# Plot function
# ----------------------------
plot_region <- function(gwas_df, eqtl_df, title){

  combined <- bind_rows(gwas_df, eqtl_df)

  p <- ggplot(combined, aes(x = pos, y = logp, color = source)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = c("GWAS" = "blue", "eQTL" = "red")) +
    theme_minimal() +
    labs(
      title = title,
      x = "Genomic Position",
      y = "-log10(p-value)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.title = element_blank()
    )

  return(p)
}

# ----------------------------
# Example usage
# ----------------------------

# Load GWAS region
gwas <- fread("data/regions/GLP1R/peak_ins_chr6_*.tsv") %>% as.data.frame()

# Load eQTL (user prepared)
# eqtl <- read.csv("data/eqtl/glp1r_pancreas.csv")

# Prepare
gwas_p <- prepare_gwas(gwas)
# eqtl_p <- prepare_eqtl(eqtl)

# Plot
# p <- plot_region(gwas_p, eqtl_p, "GLP1R - Peak Insulin Colocalization")

# Save
# ggsave("results/glp1r_peak_ins_coloc.png", p, width = 8, height = 5)
