# ================================
# 02_load_gwas_data.R
# Load all MAGIC GWAS datasets
# ================================

source("scripts/01_read_magic_vcf.R")

gwas_path <- "data/raw/"

ins_sec   <- read_magic_vcf(paste0(gwas_path, "ebi-a-GCST006679.vcf.gz"))
peak_ins  <- read_magic_vcf(paste0(gwas_path, "ebi-a-GCST006676.vcf.gz"))
proins    <- read_magic_vcf(paste0(gwas_path, "ebi-a-GCST001212.vcf.gz"))
isi       <- read_magic_vcf(paste0(gwas_path, "ebi-a-GCST003658.vcf.gz"))

homa_b    <- read_magic_vcf(paste0(gwas_path, "ieu-b-117.vcf.gz"))
homa_ir   <- read_magic_vcf(paste0(gwas_path, "ieu-b-118.vcf.gz"))

# Save objects
save(ins_sec, peak_ins, proins, isi, homa_b, homa_ir,
     file = "data/raw/gwas_raw.RData")
